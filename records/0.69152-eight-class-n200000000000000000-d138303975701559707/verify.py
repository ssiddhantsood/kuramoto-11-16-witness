#!/usr/bin/env python3
"""Standalone replay of the compact eight-class record specification.

This intentionally imports no approach-156 producer or certificate code.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import mpmath as mp
import numpy as np


def parse(spec: dict[str, Any]):
    sizes = [int(row["size"]) for row in spec["classes"]]
    if len(sizes) != 8 or sum(sizes) != int(spec["vertex_count"]):
        raise AssertionError("invalid class sizes")
    blocks = {}
    directed = [[0] * 8 for _ in range(8)]
    edges = sum(size * (size - 1) // 2 for size in sizes)
    for row in spec["cross_blocks"]:
        i, j = map(int, row["classes"])
        if (i, j) in blocks or not 0 <= i < j < 8:
            raise AssertionError("bad block list")
        blocks[i, j] = row
        if row["type"] == "absent":
            dij = dji = count = 0
        elif row["type"] == "complete":
            dij, dji = sizes[j], sizes[i]
            count = sizes[i] * sizes[j]
        else:
            modulus = int(row["modulus"])
            selected = int(row["shift_count"])
            if (
                row["type"] != "residue_orbits"
                or modulus != math.gcd(sizes[i], sizes[j])
                or not 0 < selected < modulus
            ):
                raise AssertionError("invalid residue block")
            dij = selected * sizes[j] // modulus
            dji = selected * sizes[i] // modulus
            count = sizes[i] * dij
            if count != sizes[j] * dji:
                raise AssertionError("non-biregular block")
        directed[i][j], directed[j][i] = dij, dji
        edges += count
    if len(blocks) != 28:
        raise AssertionError("not all class pairs specified")
    return sizes, blocks, directed, edges


def oracle(
    class_u: int,
    local_u: int,
    class_v: int,
    local_v: int,
    blocks: dict[tuple[int, int], dict[str, Any]],
) -> bool:
    if class_u == class_v:
        return local_u != local_v
    if class_u > class_v:
        class_u, class_v = class_v, class_u
        local_u, local_v = local_v, local_u
    block = blocks[class_u, class_v]
    if block["type"] == "absent":
        return False
    if block["type"] == "complete":
        return True
    modulus = int(block["modulus"])
    return (
        local_v - local_u - int(block["shift_start"])
    ) % modulus < int(block["shift_count"])


def solve(spec: dict[str, Any], directed: list[list[int]], dps: int):
    mp.mp.dps = dps
    start = [mp.mpf(value) for value in spec["equilibrium"]["phase_seed"][:4]]

    def equations(a, b, c, d):
        phase = [a, b, c, d, -a, -b, -c, -d]
        return tuple(
            mp.fsum(
                directed[i][j] * mp.sin(phase[i] - phase[j])
                for j in range(8)
            )
            for i in range(4)
        )

    angles = mp.findroot(
        equations,
        tuple(start),
        solver="mnewton",
        tol=mp.mpf(10) ** (-(dps - 15)),
        maxsteps=80,
    )
    phase = list(angles) + [-value for value in angles]
    residual = max(
        abs(
            mp.fsum(
                directed[i][j] * mp.sin(phase[i] - phase[j])
                for j in range(8)
            )
        )
        for i in range(8)
    )
    return phase, residual


def quotient(directed: list[list[int]], phase: list[mp.mpf]):
    matrix = mp.matrix(8)
    for i in range(8):
        for j in range(8):
            matrix[i, j] = 0
    for i in range(8):
        for j in range(i + 1, 8):
            if not directed[i][j]:
                continue
            cosine = mp.cos(phase[i] - phase[j])
            matrix[i, i] += directed[i][j] * cosine
            matrix[j, j] += directed[j][i] * cosine
            value = -mp.sqrt(directed[i][j] * directed[j][i]) * cosine
            matrix[i, j] = matrix[j, i] = value
    values, _ = mp.eigsy(matrix)
    return [values[i] for i in range(8)]


def transverse(
    sizes: list[int],
    blocks: dict[tuple[int, int], dict[str, Any]],
    directed: list[list[int]],
    phase: list[mp.mpf],
):
    matrix = mp.matrix(8)
    for i in range(8):
        for j in range(8):
            matrix[i, j] = 0
        matrix[i, i] = sizes[i] + mp.fsum(
            directed[i][j] * mp.cos(phase[i] - phase[j])
            for j in range(8)
            if j != i
        )
    for (i, j), block in blocks.items():
        if block["type"] != "residue_orbits":
            continue
        present = directed[i][j] * directed[j][i]
        omitted = (sizes[j] - directed[i][j]) * (
            sizes[i] - directed[j][i]
        )
        value = -abs(mp.cos(phase[i] - phase[j])) * mp.sqrt(
            min(present, omitted)
        )
        matrix[i, j] = matrix[j, i] = value
    values, _ = mp.eigsy(matrix)
    return [values[i] for i in range(8)]


def connected(blocks: dict[tuple[int, int], dict[str, Any]]) -> bool:
    adjacency = [set() for _ in range(8)]
    for (i, j), row in blocks.items():
        if row["type"] != "absent":
            adjacency[i].add(j)
            adjacency[j].add(i)
    seen = {0}
    stack = [0]
    while stack:
        i = stack.pop()
        for j in adjacency[i] - seen:
            seen.add(j)
            stack.append(j)
    return len(seen) == 8


def main() -> None:
    directory = Path(__file__).resolve().parent.parent
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--spec",
        type=Path,
        default=directory / "artifacts" / "highest_record_graph_spec.json",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=directory / "artifacts" / "independent_replay.json",
    )
    parser.add_argument("--dps", type=int, default=110)
    args = parser.parse_args()
    raw = args.spec.read_bytes()
    spec = json.loads(raw)
    sizes, blocks, directed, edge_count = parse(spec)
    phase, residual = solve(spec, directed, args.dps)
    qvalues = quotient(directed, phase)
    tvalues = transverse(sizes, blocks, directed, phase)
    total = sum(sizes)
    degrees = [
        sizes[i] - 1 + sum(directed[i]) for i in range(8)
    ]
    minimum = min(degrees)
    order = abs(
        mp.fsum(sizes[i] * mp.exp(1j * phase[i]) for i in range(8))
        / total
    )

    rng = np.random.default_rng(156)
    oracle_checks = 0
    for (i, j), _block in blocks.items():
        probes = [
            (0, 0),
            (sizes[i] - 1, sizes[j] - 1),
            (
                int(rng.integers(0, sizes[i])),
                int(rng.integers(0, sizes[j])),
            ),
        ]
        for u, v in probes:
            if oracle(i, u, j, v, blocks) != oracle(j, v, i, u, blocks):
                raise AssertionError("asymmetric adjacency oracle")
            oracle_checks += 1

    normalized_gap = qvalues[1] / total
    normalized_transverse = tvalues[0] / total
    one_zero = (
        abs(qvalues[0]) < mp.mpf("1e-80")
        and all(value > 0 for value in qvalues[1:])
        and tvalues[0] > 0
    )
    accepted = (
        minimum * 460_799 > 316_802 * (total - 1)
        and residual < mp.mpf("1e-14")
        and normalized_gap > mp.mpf("1e-6")
        and one_zero
        and connected(blocks)
        and order < mp.mpf("0.99")
    )
    report = {
        "schema": "independent-eight-class-record-replay-v1",
        "outcome": "CONFIRMED" if accepted else "REJECTED",
        "independence": {
            "imports_producer_code": False,
            "input_spec": str(args.spec.resolve()),
            "input_sha256": hashlib.sha256(raw).hexdigest(),
        },
        "graph": {
            "simple": True,
            "connected": connected(blocks),
            "vertex_count": total,
            "edge_count": edge_count,
            "class_sizes": sizes,
            "class_degrees": degrees,
            "minimum_degree": minimum,
            "minimum_degree_ratio": f"{minimum}/{total - 1}",
            "oracle_symmetry_checks": oracle_checks,
            "beats_approach147_cross_product": (
                minimum * 460_799 - 316_802 * (total - 1)
            ),
        },
        "equilibrium": {
            "phases": [mp.nstr(value, 90) for value in phase],
            "torque_residual": mp.nstr(residual, 40),
            "order_parameter": mp.nstr(order, 50),
        },
        "spectrum": {
            "quotient_values": [mp.nstr(value, 50) for value in qvalues],
            "normalized_gap": mp.nstr(normalized_gap, 50),
            "transverse_comparison_values": [
                mp.nstr(value, 50) for value in tvalues
            ],
            "normalized_transverse_gap": mp.nstr(
                normalized_transverse, 50
            ),
            "exactly_one_zero": one_zero,
            "dimension_split": f"8+{total - 8}={total}",
        },
        "accepted": accepted,
    }
    args.out.write_text(json.dumps(report, indent=2) + "\n")
    print(
        json.dumps(
            {
                "outcome": report["outcome"],
                "N": total,
                "degree": minimum,
                "mu": minimum / (total - 1),
                "normalized_gap": mp.nstr(normalized_gap, 20),
                "normalized_transverse": mp.nstr(
                    normalized_transverse, 20
                ),
                "torque_residual": mp.nstr(residual, 10),
                "output": str(args.out),
            },
            indent=2,
        )
    )
    if not accepted:
        raise SystemExit(1)


if __name__ == "__main__":
    main()

