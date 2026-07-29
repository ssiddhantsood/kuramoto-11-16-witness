#!/usr/bin/env python3
"""Standalone verification of the finite Kuramoto witness.

The verifier reads only ``graph_spec.json``.  It checks exact graph
combinatorics, resolves the six-class equilibrium at high precision, computes
the complete class-constant Hessian spectrum, and proves positivity on the
class-zero-sum subspace using a block-operator comparison bound.
"""

from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction
from pathlib import Path
from typing import Any

import mpmath as mp


HERE = Path(__file__).resolve().parent


def parse_spec(path: Path) -> tuple[dict[str, Any], list[int], dict]:
    spec = json.loads(path.read_text())
    if spec.get("schema") != "six-class-circulant-witness-v1":
        raise ValueError("unsupported graph specification")
    sizes = [int(c["size"]) for c in spec["classes"]]
    if len(sizes) != 6 or sum(sizes) != int(spec["vertex_count"]):
        raise ValueError("invalid class sizes")

    blocks = {}
    for block in spec["cross_blocks"]:
        i, j = map(int, block["classes"])
        if not 0 <= i < j < 6 or (i, j) in blocks:
            raise ValueError(f"invalid or duplicate pair {(i, j)}")
        if block["type"] == "residue_orbits":
            modulus = int(block["modulus"])
            if modulus != math.gcd(sizes[i], sizes[j]):
                raise ValueError(f"wrong modulus for {(i, j)}")
            count = int(block["shift_count"])
            if not 0 <= count <= modulus:
                raise ValueError(f"wrong orbit count for {(i, j)}")
        elif block["type"] not in {"complete", "absent"}:
            raise ValueError(f"unknown block type {block['type']}")
        blocks[i, j] = block
    if len(blocks) != 15:
        raise ValueError("all 15 cross-class pairs must be specified")
    return spec, sizes, blocks


def block_data(sizes: list[int], blocks: dict) -> tuple[list[list[int]], int]:
    degrees = [[0] * 6 for _ in range(6)]
    edge_count = sum(n * (n - 1) // 2 for n in sizes)
    for (i, j), block in blocks.items():
        if block["type"] == "absent":
            edges = 0
        elif block["type"] == "complete":
            edges = sizes[i] * sizes[j]
        else:
            edges = (
                sizes[i]
                * sizes[j]
                // int(block["modulus"])
                * int(block["shift_count"])
            )
        if edges % sizes[i] or edges % sizes[j]:
            raise ValueError(f"block {(i, j)} is not biregular")
        degrees[i][j] = edges // sizes[i]
        degrees[j][i] = edges // sizes[j]
        edge_count += edges
    return degrees, edge_count


def support_connected(blocks: dict) -> bool:
    adjacency = [[] for _ in range(6)]
    for (i, j), block in blocks.items():
        present = block["type"] == "complete" or (
            block["type"] == "residue_orbits"
            and int(block["shift_count"]) > 0
        )
        if present:
            adjacency[i].append(j)
            adjacency[j].append(i)
    seen, stack = {0}, [0]
    while stack:
        i = stack.pop()
        for j in adjacency[i]:
            if j not in seen:
                seen.add(j)
                stack.append(j)
    return len(seen) == 6


def edge(
    class_u: int,
    local_u: int,
    class_v: int,
    local_v: int,
    blocks: dict,
) -> bool:
    """Exact adjacency oracle for any pair of vertices in local coordinates."""
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
    delta = (
        local_v - local_u - int(block["shift_start"])
    ) % int(block["modulus"])
    return delta < int(block["shift_count"])


def solve_equilibrium(
    spec: dict[str, Any], degrees: list[list[int]], dps: int
) -> tuple[list[mp.mpf], mp.mpf]:
    mp.mp.dps = dps
    seed = [mp.mpf(x) for x in spec["equilibrium"]["phase_seed"]]

    def phases(s: mp.mpf, x: mp.mpf, y: mp.mpf) -> list[mp.mpf]:
        return [mp.mpf(0), s, x, s - x, s - y, y]

    def torque(i: int, s: mp.mpf, x: mp.mpf, y: mp.mpf) -> mp.mpf:
        phi = phases(s, x, y)
        return mp.fsum(
            degrees[i][j] * mp.sin(phi[i] - phi[j]) for j in range(6)
        )

    s, x, y = mp.findroot(
        lambda ss, xx, yy: (
            torque(0, ss, xx, yy),
            torque(2, ss, xx, yy),
            torque(4, ss, xx, yy),
        ),
        (seed[1], seed[2], seed[5]),
        solver="mdnewton",
        tol=mp.mpf(10) ** (-(dps - 15)),
        maxsteps=100,
    )
    phi = phases(s, x, y)
    residual = max(
        abs(
            mp.fsum(
                degrees[i][j] * mp.sin(phi[i] - phi[j])
                for j in range(6)
            )
        )
        for i in range(6)
    )
    phi[2] += 2 * mp.pi
    return phi, residual


def quotient_hessian(
    degrees: list[list[int]], phases: list[mp.mpf]
) -> mp.matrix:
    """Hessian in the orthonormal class-constant basis."""
    matrix = mp.matrix(6)
    for i in range(6):
        for j in range(6):
            matrix[i, j] = mp.mpf(0)
    for i in range(6):
        for j in range(i + 1, 6):
            if not degrees[i][j]:
                continue
            cosine = mp.cos(phases[i] - phases[j])
            matrix[i, i] += degrees[i][j] * cosine
            matrix[j, j] += degrees[j][i] * cosine
            value = -mp.sqrt(degrees[i][j] * degrees[j][i]) * cosine
            matrix[i, j] = matrix[j, i] = value
    return matrix


def transverse_comparison(
    sizes: list[int],
    degrees: list[list[int]],
    phases: list[mp.mpf],
    blocks: dict,
) -> mp.matrix:
    """Comparison matrix proving positivity on class-zero-sum vectors.

    A partial biregular block has adjacency norm at most
    sqrt(d_ij*d_ji).  Complete blocks vanish between zero-sum class
    subspaces.  Each clique contributes n_i on its zero-sum subspace.
    """
    matrix = mp.matrix(6)
    for i in range(6):
        for j in range(6):
            matrix[i, j] = mp.mpf(0)
    for i in range(6):
        matrix[i, i] = sizes[i] + mp.fsum(
            degrees[i][j] * mp.cos(phases[i] - phases[j])
            for j in range(6)
        )
    for (i, j), block in blocks.items():
        if block["type"] != "residue_orbits":
            continue
        value = -abs(mp.cos(phases[i] - phases[j])) * mp.sqrt(
            degrees[i][j] * degrees[j][i]
        )
        matrix[i, j] = matrix[j, i] = value
    return matrix


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--spec", type=Path, default=HERE / "graph_spec.json")
    parser.add_argument("--out", type=Path)
    parser.add_argument("--dps", type=int, default=100)
    args = parser.parse_args()

    spec, sizes, blocks = parse_spec(args.spec)
    degrees, edge_count = block_data(sizes, blocks)
    phases, residual = solve_equilibrium(spec, degrees, args.dps)
    quotient = quotient_hessian(degrees, phases)
    comparison = transverse_comparison(sizes, degrees, phases, blocks)
    quotient_values = list(mp.eigsy(quotient, eigvals_only=True))
    comparison_values = list(mp.eigsy(comparison, eigvals_only=True))

    class_degrees = [
        sizes[i] - 1 + sum(degrees[i]) for i in range(6)
    ]
    vertex_count = sum(sizes)
    minimum_degree = min(class_degrees)
    excess = 16 * minimum_degree - 11 * (vertex_count - 1)
    order_parameter = abs(
        mp.fsum(sizes[i] * mp.exp(1j * phases[i]) for i in range(6))
        / vertex_count
    )
    dimension_count = 6 + sum(n - 1 for n in sizes)

    accepted = (
        excess > 0
        and residual < mp.mpf("1e-60")
        and abs(quotient_values[0]) < mp.mpf("1e-70")
        and all(v > 0 for v in quotient_values[1:])
        and comparison_values[0] > 0
        and support_connected(blocks)
        and dimension_count == vertex_count
    )
    report = {
        "outcome": "WITNESS" if accepted else "REJECTED",
        "graph": {
            "vertex_count": vertex_count,
            "edge_count": edge_count,
            "class_sizes": sizes,
            "class_degrees": class_degrees,
            "minimum_degree": minimum_degree,
            "minimum_degree_ratio_exact": str(
                Fraction(minimum_degree, vertex_count - 1)
            ),
            "11_16_cross_product_excess": excess,
            "simple": True,
            "connected": support_connected(blocks),
        },
        "equilibrium": {
            "phases": [mp.nstr(v, 90) for v in phases],
            "torque_residual": mp.nstr(residual, 30),
            "order_parameter": mp.nstr(order_parameter, 30),
        },
        "hessian": {
            "quotient_eigenvalues": [
                mp.nstr(v, 40) for v in quotient_values
            ],
            "quotient_gap": mp.nstr(quotient_values[1], 40),
            "normalized_gap": mp.nstr(
                quotient_values[1] / vertex_count, 40
            ),
            "transverse_comparison_eigenvalues": [
                mp.nstr(v, 40) for v in comparison_values
            ],
            "transverse_lower_bound": mp.nstr(
                comparison_values[0], 40
            ),
            "dimension_accounting": {
                "quotient": 6,
                "transverse": sum(n - 1 for n in sizes),
                "total": dimension_count,
            },
            "exactly_one_zero": bool(
                abs(quotient_values[0]) < mp.mpf("1e-70")
                and all(v > 0 for v in quotient_values[1:])
                and comparison_values[0] > 0
            ),
        },
        "accepted": accepted,
    }
    text = json.dumps(report, indent=2) + "\n"
    if args.out:
        args.out.write_text(text)
    print(text, end="")
    if not accepted:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
