#!/usr/bin/env python3
"""Independent 1024-bit Arb audit for ultrasmall-gap graph specifications.

This verifier imports no producer or weighted-optimization code.  From each
declarative graph spec it independently reconstructs exact graph arithmetic,
proves the reflection rule for every vertex label, resolves the equilibrium
at 280 decimal digits, proves local root uniqueness by a Krawczyk inclusion,
proves a strictly positive quotient gap with interval LDL (without a fixed
1e-6 shift), proves all transverse modes positive by a Fourier-free block
comparison, and performs a rescaled local nonlinear return test.
"""

from __future__ import annotations

import argparse
from collections import deque
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Iterable, Sequence

import mpmath as mp
import numpy as np
from flint import arb, arb_mat, ctx


MP_DPS = 280
ARB_BITS = 1024
ROOT_RADIUS = "1e-170"
mp.mp.dps = MP_DPS

HERE = Path(__file__).resolve().parent
APPROACH = HERE.parent
DEFAULT_SPEC_DIR = APPROACH / "finite"
DEFAULT_CERT_DIR = APPROACH / "artifacts" / "certificates"
DEFAULT_LEDGER = APPROACH / "artifacts" / "certified_record_ledger.json"
BASELINE = Fraction(
    138304542182060000915026158,
    199999999999999999999999999,
)


def decimal(value: mp.mpf, digits: int = 260) -> str:
    return mp.nstr(value, digits, strip_zeros=False)


def arb_record(value: arb, digits: int = 180) -> dict[str, str]:
    return {
        "ball": value.str(digits, radius=True, more=True),
        "lower": value.lower().str(digits, radius=False, more=True),
        "upper": value.upper().str(digits, radius=False, more=True),
    }


def fraction_record(value: Fraction) -> dict[str, str]:
    return {
        "numerator": str(value.numerator),
        "denominator": str(value.denominator),
        "fraction": f"{value.numerator}/{value.denominator}",
        "decimal": decimal(
            mp.mpf(value.numerator) / value.denominator, 110
        ),
    }


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_integer(value: Any) -> int:
    return int(str(value))


def load_spec(path: Path) -> dict[str, Any]:
    spec = json.loads(path.read_text())
    if spec.get("schema") != "ultrasmall-gap-reflection-witness-v1":
        raise ValueError(f"unsupported graph schema in {path}")
    return spec


def reconstruct_graph(spec: dict[str, Any]) -> dict[str, Any]:
    classes = sorted(spec["classes"], key=lambda row: int(row["id"]))
    if [int(row["id"]) for row in classes] != list(range(6)):
        raise AssertionError("expected exactly class ids 0..5")
    sizes = [parse_integer(row["size"]) for row in classes]
    if any(size <= 1 for size in sizes):
        raise AssertionError("every clique fiber must be nontrivial")
    total = parse_integer(spec["vertex_count"])
    if sum(sizes) != total:
        raise AssertionError("class sizes do not sum to N")
    if any(sizes[2 * i] != sizes[2 * i + 1] for i in range(3)):
        raise AssertionError("reflected classes have unequal sizes")

    expected = {(i, j) for i in range(6) for j in range(i + 1, 6)}
    seen: set[tuple[int, int]] = set()
    by_pair: dict[tuple[int, int], dict[str, Any]] = {}
    r = [[0 for _ in range(6)] for _ in range(6)]
    edge = [[0 for _ in range(6)] for _ in range(6)]
    block_rows = []
    for raw in spec["cross_blocks"]:
        item = dict(raw)
        i, j = map(int, item["classes"])
        if not 0 <= i < j < 6 or (i, j) in seen:
            raise AssertionError(f"invalid or repeated block {(i, j)}")
        seen.add((i, j))
        by_pair[(i, j)] = item
        modulus = math.gcd(sizes[i], sizes[j])
        kind = item["type"]
        if kind == "absent":
            count = 0
            start = None
        elif kind == "complete":
            count = modulus
            start = None
        elif kind == "residue_orbits":
            if parse_integer(item["modulus"]) != modulus:
                raise AssertionError("residue modulus is not the exact gcd")
            count = parse_integer(item["shift_count"])
            start = parse_integer(item["shift_start"])
            if not 0 < count < modulus or not 0 <= start < modulus:
                raise AssertionError("invalid nontrivial residue interval")
        else:
            raise AssertionError(f"unsupported block type {kind}")
        edges = sizes[i] * sizes[j] // modulus * count
        if edges % sizes[i] or edges % sizes[j]:
            raise AssertionError("block is not biregular")
        r[i][j], r[j][i] = edges // sizes[i], edges // sizes[j]
        edge[i][j] = edge[j][i] = edges
        block_rows.append(
            {
                "classes": [i, j],
                "type": kind,
                "gcd_modulus": str(modulus),
                "shift_count": str(count),
                "shift_start": None if start is None else str(start),
                "degree_i_to_j": str(r[i][j]),
                "degree_j_to_i": str(r[j][i]),
                "edge_count": str(edges),
                "biregular_identity": (
                    f"{sizes[i]}*{r[i][j]}="
                    f"{sizes[j]}*{r[j][i]}={edges}"
                ),
            }
        )
    if seen != expected or len(spec["cross_blocks"]) != 15:
        raise AssertionError("cross-block partition is incomplete")

    degrees = [sizes[i] - 1 + sum(r[i]) for i in range(6)]
    total_edges = sum(size * (size - 1) // 2 for size in sizes)
    total_edges += sum(edge[i][j] for i in range(6) for j in range(i + 1, 6))
    handshake = sum(sizes[i] * degrees[i] for i in range(6))
    if handshake != 2 * total_edges:
        raise AssertionError("exact handshake identity failed")

    support = [[] for _ in range(6)]
    for i, j in expected:
        if edge[i][j]:
            support[i].append(j)
            support[j].append(i)
    reached = {0}
    queue = deque([0])
    while queue:
        i = queue.popleft()
        for j in support[i]:
            if j not in reached:
                reached.add(j)
                queue.append(j)
    if reached != set(range(6)):
        raise AssertionError("class support is disconnected")

    minimum_degree = min(degrees)
    mu = Fraction(minimum_degree, total - 1)
    if mu <= BASELINE:
        raise AssertionError("graph does not improve the current record")
    return {
        "sizes": sizes,
        "N": total,
        "r": r,
        "edge": edge,
        "degree": degrees,
        "E": total_edges,
        "blocks": by_pair,
        "report": {
            "vertex_count": str(total),
            "class_sizes": [str(value) for value in sizes],
            "class_degrees": [str(value) for value in degrees],
            "minimum_degree": str(minimum_degree),
            "edge_count": str(total_edges),
            "handshake_sum": str(handshake),
            "handshake_exact": True,
            "simple_unweighted": True,
            "connected": True,
            "block_arithmetic": block_rows,
            "mu": fraction_record(mu),
            "beats_current_exactly": True,
            "improvement_over_current": fraction_record(mu - BASELINE),
        },
    }


def reflection_proof(
    sizes: Sequence[int],
    blocks: dict[tuple[int, int], dict[str, Any]],
    spec: dict[str, Any],
) -> dict[str, Any]:
    involution = list(map(int, spec["reflection"]["class_involution"]))
    if involution != [1, 0, 3, 2, 5, 4]:
        raise AssertionError("unexpected class reflection")
    rows = []
    for (i, j), source in sorted(blocks.items()):
        ri, rj = involution[i], involution[j]
        order_preserved = ri < rj
        target_key = (ri, rj) if order_preserved else (rj, ri)
        target = blocks[target_key]
        if source["type"] != target["type"]:
            raise AssertionError("reflection changes a block type")
        row = {
            "source": [i, j],
            "target": list(target_key),
            "type": source["type"],
            "canonical_endpoint_order_preserved": order_preserved,
        }
        if source["type"] == "residue_orbits":
            modulus = math.gcd(sizes[i], sizes[j])
            count = parse_integer(source["shift_count"])
            start = parse_integer(source["shift_start"])
            expected_start = (
                -(start + count - 1) if order_preserved else start
            ) % modulus
            if (
                parse_integer(target["modulus"]) != modulus
                or parse_integer(target["shift_count"]) != count
                or parse_integer(target["shift_start"]) != expected_start
            ):
                raise AssertionError("reflection interval identity failed")
            row.update(
                {
                    "modulus": str(modulus),
                    "count": str(count),
                    "source_start": str(start),
                    "required_target_start": str(expected_start),
                    "universal_residue_identity": True,
                }
            )
        rows.append(row)
    return {
        "class_involution": involution,
        "local_involution": "u -> -u mod class_size",
        "class_sizes_preserved": all(
            sizes[i] == sizes[involution[i]] for i in range(6)
        ),
        "adjacency_preserved_for_every_vertex_pair": True,
        "proof": (
            "Exact block partition and modular interval identities cover all "
            "local labels; no vertex sampling is used."
        ),
        "block_orbits": rows,
    }


def reflected_phases(theta: Sequence[Any]) -> list[Any]:
    result = []
    for value in theta:
        result.extend((value, -value))
    return result


def phase_coefficients(pair_count: int) -> list[list[int]]:
    rows = []
    for pair in range(pair_count):
        plus = [0] * pair_count
        minus = [0] * pair_count
        plus[pair], minus[pair] = 1, -1
        rows.extend((plus, minus))
    return rows


def torque_mp(
    theta: Sequence[mp.mpf], r: Sequence[Sequence[int]], total: int
) -> list[mp.mpf]:
    phase = reflected_phases(theta)
    result = []
    for pair in range(len(theta)):
        i = 2 * pair
        result.append(
            sum(
                (
                    mp.mpf(r[i][j])
                    / total
                    * mp.sin(phase[j] - phase[i])
                    for j in range(len(phase))
                    if i != j and r[i][j]
                ),
                mp.mpf(0),
            )
        )
    return result


def full_torque_mp(
    theta: Sequence[mp.mpf], r: Sequence[Sequence[int]], total: int
) -> list[mp.mpf]:
    phase = reflected_phases(theta)
    return [
        sum(
            (
                mp.mpf(r[i][j])
                / total
                * mp.sin(phase[j] - phase[i])
                for j in range(len(phase))
                if i != j and r[i][j]
            ),
            mp.mpf(0),
        )
        for i in range(len(phase))
    ]


def jacobian_mp(
    theta: Sequence[mp.mpf], r: Sequence[Sequence[int]], total: int
) -> list[list[mp.mpf]]:
    pair_count = len(theta)
    phase = reflected_phases(theta)
    coefficients = phase_coefficients(pair_count)
    matrix = [
        [mp.mpf(0) for _ in range(pair_count)] for _ in range(pair_count)
    ]
    for pair in range(pair_count):
        i = 2 * pair
        for j in range(2 * pair_count):
            if i == j or not r[i][j]:
                continue
            common = (
                mp.mpf(r[i][j])
                / total
                * mp.cos(phase[j] - phase[i])
            )
            for variable in range(pair_count):
                derivative = (
                    coefficients[j][variable] - coefficients[i][variable]
                )
                matrix[pair][variable] += common * derivative
    return matrix


def solve_root(
    seed: Sequence[str], r: Sequence[Sequence[int]], total: int
) -> dict[str, Any]:
    values = [mp.mpf(value) for value in seed]
    initial = list(values)
    history = []
    for iteration in range(30):
        function = mp.matrix(torque_mp(values, r, total))
        jacobian = mp.matrix(jacobian_mp(values, r, total))
        step = mp.lu_solve(jacobian, -function)
        values = [values[i] + step[i] for i in range(len(values))]
        residual = max(abs(value) for value in function)
        step_norm = max(abs(step[i]) for i in range(step.rows))
        history.append(
            {
                "iteration": iteration,
                "residual_max": decimal(residual, 45),
                "step_max": decimal(step_norm, 45),
            }
        )
        if step_norm < mp.mpf("1e-260"):
            break
    else:
        raise AssertionError("high-precision torque solve did not converge")
    all_torques = full_torque_mp(values, r, total)
    residual = max(abs(value) for value in all_torques)
    if residual > mp.mpf("1e-245"):
        raise AssertionError("high-precision torque residual is too large")
    return {
        "root": values,
        "report": {
            "working_decimal_digits": MP_DPS,
            "seed": list(seed),
            "seed_motion": [
                decimal(values[i] - initial[i], 100)
                for i in range(len(values))
            ],
            "angles": [decimal(value, 270) for value in values],
            "newton_iterations": len(history),
            "iteration_history": history,
            "all_six_normalized_torques": [
                decimal(value, 255) for value in all_torques
            ],
            "maximum_normalized_torque_residual": decimal(residual, 255),
        },
    }


def torque_arb(
    theta: Sequence[arb], r: Sequence[Sequence[int]], total: int
) -> list[arb]:
    phase = reflected_phases(theta)
    result = []
    for pair in range(len(theta)):
        i = 2 * pair
        value = arb(0)
        for j in range(len(phase)):
            if i != j and r[i][j]:
                value += (
                    arb(r[i][j])
                    / arb(total)
                    * (phase[j] - phase[i]).sin()
                )
        result.append(value)
    return result


def full_torque_arb(
    theta: Sequence[arb], r: Sequence[Sequence[int]], total: int
) -> list[arb]:
    phase = reflected_phases(theta)
    result = []
    for i in range(len(phase)):
        value = arb(0)
        for j in range(len(phase)):
            if i != j and r[i][j]:
                value += (
                    arb(r[i][j])
                    / arb(total)
                    * (phase[j] - phase[i]).sin()
                )
        result.append(value)
    return result


def jacobian_arb(
    theta: Sequence[arb], r: Sequence[Sequence[int]], total: int
) -> list[list[arb]]:
    pair_count = len(theta)
    phase = reflected_phases(theta)
    coefficients = phase_coefficients(pair_count)
    matrix = [[arb(0) for _ in range(pair_count)] for _ in range(pair_count)]
    for pair in range(pair_count):
        i = 2 * pair
        for j in range(2 * pair_count):
            if i == j or not r[i][j]:
                continue
            common = (
                arb(r[i][j])
                / arb(total)
                * (phase[j] - phase[i]).cos()
            )
            for variable in range(pair_count):
                derivative = (
                    coefficients[j][variable] - coefficients[i][variable]
                )
                if derivative:
                    matrix[pair][variable] += common * derivative
    return matrix


def krawczyk_certificate(
    root: Sequence[mp.mpf], r: Sequence[Sequence[int]], total: int
) -> dict[str, Any]:
    ctx.prec = ARB_BITS
    pair_count = len(root)
    midpoint = [arb(decimal(value, 275)).mid() for value in root]
    radius = arb(ROOT_RADIUS).upper()
    delta = [arb(0, radius) for _ in root]
    box = [midpoint[i] + delta[i] for i in range(pair_count)]
    inverse_mp = mp.matrix(jacobian_mp(root, r, total)) ** -1
    inverse = [
        [
            arb(decimal(inverse_mp[i, j], 275)).mid()
            for j in range(pair_count)
        ]
        for i in range(pair_count)
    ]
    determinant = arb_mat(inverse).det()
    if determinant.contains(0):
        raise AssertionError("Krawczyk preconditioner may be singular")
    function = torque_arb(midpoint, r, total)
    interval_jacobian = jacobian_arb(box, r, total)
    remainder = [
        [arb(int(i == j)) for j in range(pair_count)]
        for i in range(pair_count)
    ]
    for i in range(pair_count):
        for j in range(pair_count):
            remainder[i][j] -= sum(
                (
                    inverse[i][k] * interval_jacobian[k][j]
                    for k in range(pair_count)
                ),
                arb(0),
            )
    image = []
    for i in range(pair_count):
        center = midpoint[i] - sum(
            (inverse[i][j] * function[j] for j in range(pair_count)),
            arb(0),
        )
        image.append(
            center
            + sum(
                (
                    remainder[i][j] * delta[j]
                    for j in range(pair_count)
                ),
                arb(0),
            )
        )
    inclusion = [
        box[i].contains_interior(image[i]) for i in range(pair_count)
    ]
    if not all(inclusion):
        raise AssertionError(f"Krawczyk inclusion failed: {inclusion}")
    full_box = full_torque_arb(box, r, total)
    return {
        "box": box,
        "report": {
            "arb_precision_bits": ARB_BITS,
            "nominal_box_radius": ROOT_RADIUS,
            "root_box": [arb_record(value) for value in box],
            "preconditioner_determinant": arb_record(determinant),
            "krawczyk_image": [arb_record(value) for value in image],
            "strict_interior_inclusion": inclusion,
            "unique_root_in_box": True,
            "root_box_all_six_torque_images": [
                arb_record(value) for value in full_box
            ],
        },
    }


def quotient_mp(
    theta: Sequence[mp.mpf],
    sizes: Sequence[int],
    r: Sequence[Sequence[int]],
    edge: Sequence[Sequence[int]],
    total: int,
) -> list[list[mp.mpf]]:
    phase = reflected_phases(theta)
    matrix = [[mp.mpf(0) for _ in range(6)] for _ in range(6)]
    for i in range(6):
        matrix[i][i] = sum(
            (
                mp.mpf(r[i][j])
                / total
                * mp.cos(phase[j] - phase[i])
                for j in range(6)
                if i != j and r[i][j]
            ),
            mp.mpf(0),
        )
        for j in range(i + 1, 6):
            if edge[i][j]:
                value = (
                    -mp.mpf(edge[i][j])
                    / (mp.mpf(total) * mp.sqrt(mp.mpf(sizes[i]) * sizes[j]))
                    * mp.cos(phase[j] - phase[i])
                )
                matrix[i][j] = matrix[j][i] = value
    return matrix


def quotient_arb(
    theta: Sequence[arb],
    sizes: Sequence[int],
    r: Sequence[Sequence[int]],
    edge: Sequence[Sequence[int]],
    total: int,
) -> list[list[arb]]:
    phase = reflected_phases(theta)
    matrix = [[arb(0) for _ in range(6)] for _ in range(6)]
    for i in range(6):
        diagonal = arb(0)
        for j in range(6):
            if i != j and r[i][j]:
                diagonal += (
                    arb(r[i][j])
                    / arb(total)
                    * (phase[j] - phase[i]).cos()
                )
        matrix[i][i] = diagonal
        for j in range(i + 1, 6):
            if edge[i][j]:
                scale = (
                    arb(edge[i][j])
                    / arb(sizes[i] * sizes[j]).sqrt()
                    / arb(total)
                )
                matrix[i][j] = matrix[j][i] = (
                    -scale * (phase[j] - phase[i]).cos()
                )
    return matrix


def parity_blocks(matrix: Sequence[Sequence[Any]]):
    even = [[None for _ in range(3)] for _ in range(3)]
    odd = [[None for _ in range(3)] for _ in range(3)]
    for a in range(3):
        p, n = 2 * a, 2 * a + 1
        for b in range(3):
            q, z = 2 * b, 2 * b + 1
            even[a][b] = (
                matrix[p][q]
                + matrix[p][z]
                + matrix[n][q]
                + matrix[n][z]
            ) / 2
            odd[a][b] = (
                matrix[p][q]
                - matrix[p][z]
                - matrix[n][q]
                + matrix[n][z]
            ) / 2
    return even, odd


def ldl_certificate(
    matrix: Sequence[Sequence[arb]], threshold: Fraction
) -> tuple[bool, list[arb]]:
    count = len(matrix)
    work = [[matrix[i][j] for j in range(count)] for i in range(count)]
    shift = arb(f"{threshold.numerator}/{threshold.denominator}")
    for i in range(count):
        work[i][i] -= shift
    lower = [[arb(0) for _ in range(count)] for _ in range(count)]
    diagonal = [arb(0) for _ in range(count)]
    for i in range(count):
        lower[i][i] = arb(1)
        pivot = work[i][i]
        for k in range(i):
            pivot -= lower[i][k] * lower[i][k] * diagonal[k]
        diagonal[i] = pivot
        if not (pivot.lower() > 0):
            return False, diagonal[: i + 1]
        for j in range(i + 1, count):
            numerator = work[j][i]
            for k in range(i):
                numerator -= lower[j][k] * lower[i][k] * diagonal[k]
            lower[j][i] = numerator / diagonal[i]
    return True, diagonal


def downward_decimal(value: mp.mpf, places: int) -> Fraction:
    scale = 10**places
    return Fraction(int(mp.floor(value * scale)), scale)


def matrix_mp_record(
    matrix: Sequence[Sequence[mp.mpf]], digits: int = 90
) -> list[list[str]]:
    return [[decimal(value, digits) for value in row] for row in matrix]


def matrix_arb_record(
    matrix: Sequence[Sequence[arb]], digits: int = 80
) -> list[list[str]]:
    return [
        [value.str(digits, radius=True, more=True) for value in row]
        for row in matrix
    ]


def quotient_certificate(
    root: Sequence[mp.mpf],
    root_box: Sequence[arb],
    sizes: Sequence[int],
    r: Sequence[Sequence[int]],
    edge: Sequence[Sequence[int]],
    total: int,
) -> dict[str, Any]:
    point = quotient_mp(root, sizes, r, edge, total)
    interval = quotient_arb(root_box, sizes, r, edge, total)
    even_point, odd_point = parity_blocks(point)
    even_interval, odd_interval = parity_blocks(interval)
    full_values_matrix, full_vectors = mp.eigsy(mp.matrix(point))
    full_values = [mp.mpf(full_values_matrix[i]) for i in range(6)]
    odd_values_matrix, _ = mp.eigsy(mp.matrix(odd_point))
    odd_values = [mp.mpf(odd_values_matrix[i]) for i in range(3)]
    rotation_point = [
        mp.sqrt(mp.mpf(2 * sizes[2 * i]) / total) for i in range(3)
    ]
    rotation_interval = [
        (arb(2 * sizes[2 * i]) / arb(total)).sqrt() for i in range(3)
    ]
    lifted_point = [
        [
            even_point[i][j]
            + 2 * rotation_point[i] * rotation_point[j]
            for j in range(3)
        ]
        for i in range(3)
    ]
    lifted_interval = [
        [
            even_interval[i][j]
            + 2 * rotation_interval[i] * rotation_interval[j]
            for j in range(3)
        ]
        for i in range(3)
    ]
    lifted_values_matrix, _ = mp.eigsy(mp.matrix(lifted_point))
    lifted_values = [
        mp.mpf(lifted_values_matrix[i]) for i in range(3)
    ]
    numerical_gap = min(min(lifted_values), min(odd_values))
    if numerical_gap <= 0:
        raise AssertionError("numerical quotient gap is not positive")
    candidate = numerical_gap * (1 - mp.mpf("1e-8"))
    threshold = downward_decimal(candidate, 60)
    if threshold <= 0:
        threshold = downward_decimal(numerical_gap / 2, 80)
    ok_even, pivots_even = ldl_certificate(lifted_interval, threshold)
    ok_odd, pivots_odd = ldl_certificate(odd_interval, threshold)
    if not ok_even or not ok_odd or threshold <= 0:
        raise AssertionError("Arb quotient LDL failed at positive threshold")

    sqrt_sizes = [mp.sqrt(value) for value in sizes]
    rotation_residual = [
        sum(
            (point[i][j] * sqrt_sizes[j] for j in range(6)),
            mp.mpf(0),
        )
        for i in range(6)
    ]
    weak_index = next(
        i for i, value in enumerate(full_values) if value > threshold / 2
    )
    weak = [mp.mpf(full_vectors[i, weak_index]) for i in range(6)]
    return {
        "Qmp": point,
        "weak_vector": weak,
        "numerical_gap": numerical_gap,
        "report": {
            "normalization": "mass-orthonormal class quotient divided by N",
            "full_6x6_midpoint": matrix_mp_record(point),
            "full_6x6_arb_enclosure": matrix_arb_record(interval),
            "full_eigenvalues_numerical": [
                decimal(value, 120) for value in full_values
            ],
            "lifted_even_eigenvalues_numerical": [
                decimal(value, 120) for value in lifted_values
            ],
            "odd_eigenvalues_numerical": [
                decimal(value, 120) for value in odd_values
            ],
            "actual_numerical_normalized_gap": decimal(
                numerical_gap, 120
            ),
            "certified_normalized_gap_lower": fraction_record(threshold),
            "certified_gap_strictly_positive": True,
            "fixed_acceptance_shift_used": False,
            "arb_precision_bits": ARB_BITS,
            "lifted_even_ldl_pivots": [
                arb_record(value, 130) for value in pivots_even
            ],
            "odd_ldl_pivots": [
                arb_record(value, 130) for value in pivots_odd
            ],
            "rotation_residual_numerical_max": decimal(
                max(abs(value) for value in rotation_residual), 100
            ),
            "rotation_zero_exact": True,
            "zero_count": 1,
            "zero_count_argument": (
                "The class-constant rotation vector is an exact zero by "
                "pairwise biregular cancellation. Arb LDL proves the "
                "rank-one-lifted even block and odd block exceed a strictly "
                "positive rational threshold."
            ),
            "weak_mode_mass_orthonormal_vector": [
                decimal(value, 100) for value in weak
            ],
        },
    }


def choose_upper(a: arb, b: arb) -> tuple[arb, str]:
    if a.upper() < b.upper():
        return a.upper(), "complement_count"
    return b.upper(), "cosecant"


def transverse_certificate(
    root_box: Sequence[arb],
    sizes: Sequence[int],
    r: Sequence[Sequence[int]],
    blocks: Iterable[dict[str, Any]],
    total: int,
) -> dict[str, Any]:
    phase = reflected_phases(root_box)
    beta = []
    for i in range(6):
        value = arb(sizes[i])
        for j in range(6):
            if i != j and r[i][j]:
                value += arb(r[i][j]) * (phase[j] - phase[i]).cos()
        beta.append(value / arb(total))
    comparison = [[arb(0) for _ in range(6)] for _ in range(6)]
    for i in range(6):
        comparison[i][i] = beta[i]
    norm_rows = []
    for item in blocks:
        if item["type"] != "residue_orbits":
            continue
        i, j = map(int, item["classes"])
        modulus = math.gcd(sizes[i], sizes[j])
        count = parse_integer(item["shift_count"])
        start = parse_integer(item["shift_start"])
        complement = min(count, modulus - count)
        count_bound = arb(complement)
        csc_bound = 1 / (arb.pi() / arb(modulus)).sin()
        fourier, selected = choose_upper(count_bound, csc_bound)
        multiplicity = arb(sizes[i] * sizes[j]).sqrt() / arb(modulus)
        normalized = (multiplicity * fourier / arb(total)).upper()
        weighted = (
            normalized * abs((phase[j] - phase[i]).cos())
        ).upper()
        comparison[i][j] = comparison[j][i] = -weighted
        norm_rows.append(
            {
                "classes": [i, j],
                "modulus": str(modulus),
                "count": str(count),
                "start": str(start),
                "complement_count": str(complement),
                "selected_bound": selected,
                "normalized_unweighted_operator_norm_upper": (
                    normalized.str(100, radius=False, more=True)
                ),
                "normalized_cos_weighted_operator_norm_upper": (
                    weighted.str(100, radius=False, more=True)
                ),
                "covers_all_nonconstant_characters": True,
            }
        )
    midpoint = np.array(
        [[float(value.mid()) for value in row] for row in comparison]
    )
    approximate_values = np.linalg.eigvalsh(midpoint)
    threshold = downward_decimal(
        mp.mpf(str(approximate_values[0])) - mp.mpf("1e-12"), 20
    )
    ok, pivots = ldl_certificate(comparison, threshold)
    if not ok or threshold <= 0:
        raise AssertionError("transverse Arb comparison failed")
    dimensions = [size - 1 for size in sizes]
    if sum(dimensions) != total - 6:
        raise AssertionError("transverse dimensions do not sum to N-6")
    return {
        "report": {
            "decomposition": (
                "R^N is the six-dimensional class quotient direct-sum the "
                "six classwise zero-sum spaces."
            ),
            "classwise_zero_sum_dimensions": [
                str(value) for value in dimensions
            ],
            "transverse_dimension_total": str(sum(dimensions)),
            "expected_N_minus_6": str(total - 6),
            "normalized_diagonal_beta_intervals": [
                arb_record(value, 110) for value in beta
            ],
            "comparison_matrix_arb": matrix_arb_record(comparison),
            "comparison_eigenvalues_numerical": [
                f"{value:.17g}" for value in approximate_values
            ],
            "certified_normalized_lower": fraction_record(threshold),
            "arb_precision_bits": ARB_BITS,
            "comparison_ldl_pivots": [
                arb_record(value, 110) for value in pivots
            ],
            "fourier_free": True,
            "operator_comparison_argument": (
                "Complete blocks annihilate classwise zero-sum vectors. For "
                "each residue block, the nonconstant-character norm is "
                "bounded uniformly by min(k,g-k,csc(pi/g)); Arb LDL proves "
                "the resulting scalar block comparison positive."
            ),
            "fractional_block_norms": norm_rows,
        }
    }


def order_parameter(
    root: Sequence[mp.mpf],
    root_box: Sequence[arb],
    sizes: Sequence[int],
    total: int,
) -> dict[str, Any]:
    point = sum(
        (
            mp.mpf(2 * sizes[2 * i])
            / total
            * mp.cos(root[i])
            for i in range(3)
        ),
        mp.mpf(0),
    )
    interval = sum(
        (
            arb(2 * sizes[2 * i])
            / arb(total)
            * root_box[i].cos()
            for i in range(3)
        ),
        arb(0),
    )
    return {
        "signed_real_point": decimal(point, 150),
        "magnitude_point": decimal(abs(point), 150),
        "signed_real_interval": arb_record(interval),
        "magnitude_interval": arb_record(abs(interval)),
        "imaginary_part_exact": "0",
        "nonsynchronous": not abs(interval).contains(1),
    }


def phase_separation_certificate(
    root: Sequence[mp.mpf], root_box: Sequence[arb]
) -> dict[str, Any]:
    point_phases = [
        mp.fmod(value, 2 * mp.pi) for value in reflected_phases(root)
    ]
    point_phases = [
        value + 2 * mp.pi if value < 0 else value
        for value in point_phases
    ]
    point_phases.sort()
    point_gaps = [
        point_phases[i + 1] - point_phases[i]
        for i in range(len(point_phases) - 1)
    ]
    point_gaps.append(
        point_phases[0] + 2 * mp.pi - point_phases[-1]
    )
    interval_phases = reflected_phases(root_box)
    rows = []
    all_distinct = True
    for i in range(6):
        for j in range(i + 1, 6):
            half_chord = abs(
                ((interval_phases[i] - interval_phases[j]) / 2).sin()
            )
            excludes_zero = half_chord.lower() > 0
            all_distinct = all_distinct and excludes_zero
            rows.append(
                {
                    "classes": [i, j],
                    "absolute_sine_half_difference": arb_record(
                        half_chord, 100
                    ),
                    "distinct_modulo_2pi": excludes_zero,
                }
            )
    if not all_distinct:
        raise AssertionError("phase boxes do not prove all phases distinct")
    return {
        "minimum_circular_separation_point": decimal(
            min(point_gaps), 130
        ),
        "all_six_phases_distinct_modulo_2pi": True,
        "arb_precision_bits": ARB_BITS,
        "proof": (
            "For every class pair Arb proves "
            "|sin((phi_i-phi_j)/2)|>0 over the complete root box."
        ),
        "pairwise_interval_checks": rows,
    }


def nonlinear_return(
    root: Sequence[mp.mpf],
    quotient: Sequence[Sequence[mp.mpf]],
    weak_vector: Sequence[mp.mpf],
    numerical_gap: mp.mpf,
    sizes: Sequence[int],
    r: Sequence[Sequence[int]],
    edge: Sequence[Sequence[int]],
    total: int,
) -> dict[str, Any]:
    phases = reflected_phases(root)
    mass = [mp.mpf(size) / total for size in sizes]
    weak = [mp.mpf(value) for value in weak_vector]
    rotation = [mp.sqrt(value) for value in mass]
    rotation_component = sum(
        (weak[i] * rotation[i] for i in range(6)), mp.mpf(0)
    )
    weak = [
        weak[i] - rotation_component * rotation[i] for i in range(6)
    ]
    weak_norm = mp.sqrt(sum((value * value for value in weak), mp.mpf(0)))
    weak = [value / weak_norm for value in weak]
    physical = [weak[i] / mp.sqrt(mass[i]) for i in range(6)]
    rate = mp.mpf(numerical_gap)
    relative_amplitudes = (
        mp.mpf("1e-4"),
        mp.mpf("3e-4"),
        mp.mpf("1e-3"),
    )
    scaled_horizon = mp.mpf(20)
    rows = []
    largest = mp.mpf(0)

    def aligned_norm(displacement: Sequence[mp.mpf]) -> mp.mpf:
        center = sum(
            (mass[i] * displacement[i] for i in range(6)), mp.mpf(0)
        )
        return mp.sqrt(
            sum(
                (
                    mass[i] * (displacement[i] - center) ** 2
                    for i in range(6)
                ),
                mp.mpf(0),
            )
        )

    def energy_difference(displacement: Sequence[mp.mpf]) -> mp.mpf:
        terms = []
        for i in range(6):
            for j in range(i + 1, 6):
                if edge[i][j]:
                    shift = displacement[i] - displacement[j]
                    delta = phases[i] - phases[j]
                    weight = mp.mpf(edge[i][j]) / (mp.mpf(total) * total)
                    terms.append(
                        2
                        * weight
                        * mp.sin(delta + shift / 2)
                        * mp.sin(shift / 2)
                    )
        return mp.fsum(terms)

    def flow_difference(
        displacement: Sequence[mp.mpf],
    ) -> list[mp.mpf]:
        output = []
        for i in range(6):
            value = mp.mpf(0)
            for j in range(6):
                if i == j or not r[i][j]:
                    continue
                delta = phases[j] - phases[i]
                shift = displacement[j] - displacement[i]
                sine_change = (
                    2
                    * mp.cos(delta + shift / 2)
                    * mp.sin(shift / 2)
                )
                value += mp.mpf(r[i][j]) / total * sine_change
            output.append(value)
        return output

    def flow_jacobian(
        displacement: Sequence[mp.mpf],
    ) -> list[list[mp.mpf]]:
        matrix = [[mp.mpf(0) for _ in range(6)] for _ in range(6)]
        for i in range(6):
            diagonal = mp.mpf(0)
            for j in range(6):
                if i == j or not r[i][j]:
                    continue
                value = (
                    mp.mpf(r[i][j])
                    / total
                    * mp.cos(
                        phases[j]
                        - phases[i]
                        + displacement[j]
                        - displacement[i]
                    )
                )
                matrix[i][j] = value
                diagonal -= value
            matrix[i][i] = diagonal
        return matrix

    def implicit_return(
        amplitude: mp.mpf,
        sign: int,
        step_count: int,
    ) -> dict[str, Any]:
        step = scaled_horizon / step_count
        previous = [mp.mpf(sign) * value for value in physical]
        initial = [amplitude * value for value in previous]
        initial_flow = flow_difference(initial)
        projected = sum(
            (
                weak[i]
                * mp.sqrt(mass[i])
                * initial_flow[i]
                / (amplitude * rate)
                for i in range(6)
            ),
            mp.mpf(0),
        )
        maximum_newton_residual = mp.mpf(0)
        maximum_newton_steps = 0
        for _ in range(step_count):
            current = list(previous)
            for newton_step in range(20):
                displacement = [amplitude * value for value in current]
                flow = flow_difference(displacement)
                residual = [
                    current[i]
                    - previous[i]
                    - step * flow[i] / (amplitude * rate)
                    for i in range(6)
                ]
                residual_norm = max(abs(value) for value in residual)
                maximum_newton_residual = max(
                    maximum_newton_residual, residual_norm
                )
                jacobian = flow_jacobian(displacement)
                newton_matrix = mp.matrix(
                    [
                        [
                            mp.mpf(int(i == j))
                            - step * jacobian[i][j] / rate
                            for j in range(6)
                        ]
                        for i in range(6)
                    ]
                )
                correction = mp.lu_solve(
                    newton_matrix, -mp.matrix(residual)
                )
                current = [
                    current[i] + correction[i] for i in range(6)
                ]
                correction_norm = max(
                    abs(correction[i]) for i in range(correction.rows)
                )
                maximum_newton_steps = max(
                    maximum_newton_steps, newton_step + 1
                )
                if correction_norm < mp.mpf("1e-80"):
                    break
            else:
                raise AssertionError("implicit nonlinear step did not converge")
            center = sum(
                (mass[i] * current[i] for i in range(6)), mp.mpf(0)
            )
            previous = [value - center for value in current]
        final = [amplitude * value for value in previous]
        initial_distance = aligned_norm(initial)
        final_distance = aligned_norm(final)
        return {
            "projected": projected,
            "initial": initial,
            "final": final,
            "return_ratio": final_distance / initial_distance,
            "energy_initial": energy_difference(initial),
            "energy_final": energy_difference(final),
            "step_count": step_count,
            "scaled_step": step,
            "maximum_newton_steps": maximum_newton_steps,
            "maximum_pre_correction_residual": maximum_newton_residual,
        }

    for relative in relative_amplitudes:
        amplitude = rate * relative
        both = True
        for sign in (-1, 1):
            coarse = implicit_return(amplitude, sign, 40)
            fine = implicit_return(amplitude, sign, 80)
            return_ratio = fine["return_ratio"]
            convergence_error = abs(
                fine["return_ratio"] - coarse["return_ratio"]
            )
            restoring = sign * fine["projected"] < 0
            returned = bool(
                restoring
                and return_ratio < mp.mpf("1e-4")
                and convergence_error < mp.mpf("2e-7")
                and fine["energy_initial"] > 0
                and fine["energy_final"] >= 0
                and fine["energy_final"] < fine["energy_initial"]
            )
            both = both and returned
            rows.append(
                {
                    "sign": sign,
                    "relative_amplitude_to_gap": decimal(relative, 20),
                    "weighted_rms_amplitude": decimal(amplitude, 30),
                    "scaled_time_horizon": decimal(scaled_horizon, 20),
                    "physical_time_horizon": decimal(
                        scaled_horizon / rate, 30
                    ),
                    "initial_projected_scaled_flow": decimal(
                        fine["projected"], 40
                    ),
                    "restoring_sign": restoring,
                    "return_ratio_fine": decimal(return_ratio, 40),
                    "return_ratio_coarse": decimal(
                        coarse["return_ratio"], 40
                    ),
                    "step_refinement_difference": decimal(
                        convergence_error, 40
                    ),
                    "energy_above_root_initial": decimal(
                        fine["energy_initial"], 60
                    ),
                    "energy_above_root_final": decimal(
                        fine["energy_final"], 60
                    ),
                    "integrator": (
                        "high-precision L-stable backward Euler, verified "
                        "at 40 and 80 scaled-time steps"
                    ),
                    "maximum_newton_steps": fine[
                        "maximum_newton_steps"
                    ],
                    "returned": returned,
                }
            )
        if both:
            largest = amplitude
    if largest == 0:
        raise AssertionError("no symmetric nonlinear return was observed")
    weak_rate = mp.fsum(
        (
            weak[i] * mp.fsum(
                (
                    quotient[i][j] * weak[j]
                    for j in range(6)
                )
            )
            for i in range(6)
        )
    )
    return {
        "model": (
            "Exact six-class equitable Kuramoto flow integrated in scaled "
            "time s=lambda_min*t and scaled displacement u/amplitude with "
            "high-precision L-stable implicit steps."
        ),
        "weak_linear_rate_numerical": decimal(weak_rate, 50),
        "tests": rows,
        "symmetric_finite_basin_demonstrated_weighted_rms": (
            decimal(largest, 30)
        ),
        "tiny_basin_disclosure": (
            "The demonstrated basin is deliberately scaled below the "
            "ultrasmall linear gap; no macroscopic basin is claimed."
        ),
        "verdict": "LOCAL_RETURN_NUMERICALLY_OBSERVED_BOTH_SIGNS",
    }


def verify_one(spec_path: Path, output_path: Path) -> dict[str, Any]:
    spec = load_spec(spec_path)
    graph = reconstruct_graph(spec)
    reflection = reflection_proof(
        graph["sizes"], graph["blocks"], spec
    )
    root_result = solve_root(
        [str(value) for value in spec["equilibrium"]["angle_seed"]],
        graph["r"],
        graph["N"],
    )
    root = root_result["root"]
    krawczyk = krawczyk_certificate(root, graph["r"], graph["N"])
    quotient = quotient_certificate(
        root,
        krawczyk["box"],
        graph["sizes"],
        graph["r"],
        graph["edge"],
        graph["N"],
    )
    transverse = transverse_certificate(
        krawczyk["box"],
        graph["sizes"],
        graph["r"],
        spec["cross_blocks"],
        graph["N"],
    )
    order = order_parameter(
        root, krawczyk["box"], graph["sizes"], graph["N"]
    )
    phase_geometry = phase_separation_certificate(root, krawczyk["box"])
    nonlinear = nonlinear_return(
        root,
        quotient["Qmp"],
        quotient["weak_vector"],
        quotient["numerical_gap"],
        graph["sizes"],
        graph["r"],
        graph["edge"],
        graph["N"],
    )
    result = {
        "schema": "ultrasmall-gap-independent-arb-certificate-v1",
        "audit_identity": {
            "method": "independent reconstruction from graph spec only",
            "graph_spec_path": str(spec_path.resolve()),
            "graph_spec_sha256": sha256(spec_path),
            "producer_code_imported": False,
            "mpmath_decimal_digits": MP_DPS,
            "arb_precision_bits": ARB_BITS,
        },
        "source_gamma": spec["source_weighted_point"]["gamma"],
        "graph": graph["report"],
        "reflection": reflection,
        "root": root_result["report"],
        "krawczyk": krawczyk["report"],
        "order_parameter": order,
        "phase_geometry": phase_geometry,
        "quotient": quotient["report"],
        "transverse": transverse["report"],
        "nonlinear": nonlinear,
        "acceptance": {
            "exact_mu_improves_current": True,
            "simple_connected_unweighted": True,
            "unique_local_nonsynchronous_root": True,
            "all_required_phases_distinct": True,
            "exactly_one_hessian_zero": True,
            "strictly_positive_quotient_gap": True,
            "strictly_positive_all_mode_transverse_bound": True,
            "local_nonlinear_return_both_signs": True,
            "accepted": True,
        },
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(result, indent=2) + "\n")
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--spec-dir", type=Path, default=DEFAULT_SPEC_DIR)
    parser.add_argument("--certificate-dir", type=Path, default=DEFAULT_CERT_DIR)
    parser.add_argument("--ledger", type=Path, default=DEFAULT_LEDGER)
    parser.add_argument("--spec", type=Path)
    args = parser.parse_args()
    ctx.prec = ARB_BITS

    if args.spec is not None:
        paths = [args.spec]
    else:
        paths = sorted(args.spec_dir.glob("*_graph_spec.json"))
    if not paths:
        raise RuntimeError("no graph specs found")

    rows = []
    for path in paths:
        output = args.certificate_dir / path.name.replace(
            "_graph_spec.json", "_certificate.json"
        )
        result = verify_one(path, output)
        row = {
            "gamma": result["source_gamma"],
            "graph_spec": str(path.resolve()),
            "graph_spec_sha256": result["audit_identity"][
                "graph_spec_sha256"
            ],
            "certificate": str(output.resolve()),
            "certificate_sha256": sha256(output),
            "vertex_count": result["graph"]["vertex_count"],
            "minimum_degree": result["graph"]["minimum_degree"],
            "mu": result["graph"]["mu"],
            "actual_numerical_normalized_gap": result["quotient"][
                "actual_numerical_normalized_gap"
            ],
            "certified_normalized_gap_lower": result["quotient"][
                "certified_normalized_gap_lower"
            ],
            "certified_transverse_lower": result["transverse"][
                "certified_normalized_lower"
            ],
            "demonstrated_basin": result["nonlinear"][
                "symmetric_finite_basin_demonstrated_weighted_rms"
            ],
            "accepted": result["acceptance"]["accepted"],
        }
        rows.append(row)
        print(
            f"gamma={row['gamma']:>6} mu={row['mu']['decimal']} "
            f"gap>{row['certified_normalized_gap_lower']['decimal']} "
            f"accepted={int(row['accepted'])}",
            flush=True,
        )
    rows.sort(
        key=lambda row: Fraction(
            int(row["mu"]["numerator"]), int(row["mu"]["denominator"])
        )
    )
    ledger = {
        "schema": "ultrasmall-gap-certified-record-ledger-v1",
        "independent_verifier": str(Path(__file__).resolve()),
        "mpmath_decimal_digits": MP_DPS,
        "arb_precision_bits": ARB_BITS,
        "fixed_gap_floor_used": False,
        "records": rows,
        "all_accepted": all(row["accepted"] for row in rows),
        "all_strict_successive_records": all(
            Fraction(
                int(rows[i]["mu"]["numerator"]),
                int(rows[i]["mu"]["denominator"]),
            )
            > (
                BASELINE
                if i == 0
                else Fraction(
                    int(rows[i - 1]["mu"]["numerator"]),
                    int(rows[i - 1]["mu"]["denominator"]),
                )
            )
            for i in range(len(rows))
        ),
        "highest": rows[-1],
    }
    args.ledger.parent.mkdir(parents=True, exist_ok=True)
    args.ledger.write_text(json.dumps(ledger, indent=2) + "\n")
    print(str(args.ledger.resolve()))


if __name__ == "__main__":
    main()
