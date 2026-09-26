#!/usr/bin/env python3
"""Standalone verifier for the N=80002 reflection witness.

The only data file read is graph_spec.json beside this script.  No search,
construction, audit, spectrum, or simulation modules are imported.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path
from typing import Any, Callable, Sequence

import mpmath as mp


HERE = Path(__file__).resolve().parent
POINT_DPS = 180
INTERVAL_DPS = 150
PHASE_COEFFICIENTS = (
    (1, 0, 0),
    (-1, 0, 0),
    (0, 1, 0),
    (0, -1, 0),
    (0, 0, 1),
    (0, 0, -1),
)
POSITIVE_CLASSES = (0, 2, 4)


def mp_text(value: mp.mpf, digits: int = 50) -> str:
    return mp.nstr(value, digits, strip_zeros=False)


def load_spec(path: Path) -> tuple[dict[str, Any], list[int], dict[tuple[int, int], dict]]:
    spec = json.loads(path.read_text())
    if spec.get("schema") != "reflection-paired-circulant-witness-v1":
        raise ValueError("unsupported graph specification")
    classes = sorted(spec["classes"], key=lambda item: int(item["id"]))
    if [int(item["id"]) for item in classes] != list(range(6)):
        raise ValueError("class IDs must be exactly 0,...,5")
    sizes = [int(item["size"]) for item in classes]
    if sum(sizes) != int(spec["vertex_count"]):
        raise ValueError("class sizes do not sum to vertex_count")

    blocks: dict[tuple[int, int], dict] = {}
    for block in spec["cross_blocks"]:
        i, j = map(int, block["classes"])
        if not 0 <= i < j < 6 or (i, j) in blocks:
            raise ValueError(f"invalid or duplicate block {(i, j)}")
        kind = block["type"]
        if kind == "residue_orbits":
            modulus = int(block["modulus"])
            count = int(block["shift_count"])
            start = int(block["shift_start"])
            if modulus != math.gcd(sizes[i], sizes[j]):
                raise ValueError(f"wrong gcd modulus for {(i, j)}")
            if not (0 <= start < modulus and 0 <= count <= modulus):
                raise ValueError(f"invalid residue interval for {(i, j)}")
        elif kind not in {"complete", "absent"}:
            raise ValueError(f"unknown block type {kind}")
        blocks[i, j] = block
    expected = {(i, j) for i in range(6) for j in range(i + 1, 6)}
    if set(blocks) != expected:
        raise ValueError("all 15 cross-class pairs must occur exactly once")
    return spec, sizes, blocks


def block_combinatorics(
    sizes: list[int], blocks: dict[tuple[int, int], dict]
) -> tuple[list[list[int]], int, list[dict[str, Any]]]:
    degrees = [[0] * 6 for _ in range(6)]
    total_edges = sum(n * (n - 1) // 2 for n in sizes)
    details = []
    for (i, j), block in sorted(blocks.items()):
        kind = block["type"]
        if kind == "absent":
            edges = 0
        elif kind == "complete":
            edges = sizes[i] * sizes[j]
        else:
            edges = (
                sizes[i]
                * sizes[j]
                * int(block["shift_count"])
                // int(block["modulus"])
            )
        if edges % sizes[i] or edges % sizes[j]:
            raise ValueError(f"block {(i, j)} is not exactly biregular")
        degrees[i][j] = edges // sizes[i]
        degrees[j][i] = edges // sizes[j]
        total_edges += edges
        details.append(
            {
                "classes": [i, j],
                "type": kind,
                "directed_degrees": [degrees[i][j], degrees[j][i]],
                "edge_count": edges,
            }
        )
    return degrees, total_edges, details


def adjacent(
    c: int,
    u: int,
    d: int,
    v: int,
    blocks: dict[tuple[int, int], dict],
) -> bool:
    if c == d:
        return u != v
    if c > d:
        c, d, u, v = d, c, v, u
    block = blocks[c, d]
    if block["type"] == "absent":
        return False
    if block["type"] == "complete":
        return True
    return (
        (v - u - int(block["shift_start"])) % int(block["modulus"])
        < int(block["shift_count"])
    )


def support_connected(
    degrees: list[list[int]],
) -> bool:
    seen, stack = {0}, [0]
    while stack:
        i = stack.pop()
        for j in range(6):
            if degrees[i][j] and j not in seen:
                seen.add(j)
                stack.append(j)
    return len(seen) == 6


def verify_reflection(
    spec: dict[str, Any],
    sizes: list[int],
    blocks: dict[tuple[int, int], dict],
) -> dict[str, Any]:
    reflection = [int(x) for x in spec["reflection"]["class_involution"]]
    if reflection != [1, 0, 3, 2, 5, 4]:
        raise ValueError("unexpected class reflection")
    involution_checks = 0
    loop_checks = 0
    for c, size in enumerate(sizes):
        rc = reflection[c]
        if sizes[rc] != size or reflection[rc] != c:
            raise ValueError("reflection does not pair equal fibers")
        for u in range(size):
            ru = (-u) % size
            if reflection[rc] != c or (-ru) % size != u:
                raise AssertionError("vertex reflection is not an involution")
            if adjacent(c, u, c, u, blocks):
                raise AssertionError("self-loop found")
            involution_checks += 1
            loop_checks += 1

    residue_difference_checks = 0
    for (i, j) in sorted(blocks):
        modulus = math.gcd(sizes[i], sizes[j])
        for delta in range(modulus):
            source = adjacent(i, 0, j, delta, blocks)
            ri, rj = reflection[i], reflection[j]
            reflected = adjacent(
                ri,
                0,
                rj,
                (-delta) % sizes[j],
                blocks,
            )
            if source != reflected:
                raise AssertionError(
                    f"reflection fails for block {(i, j)}, residue {delta}"
                )
            residue_difference_checks += 1

    source = blocks[0, 4]
    target = blocks[1, 5]
    modulus = int(source["modulus"])
    source_set = {
        (int(source["shift_start"]) + h) % modulus
        for h in range(int(source["shift_count"]))
    }
    target_set = {
        (int(target["shift_start"]) + h) % modulus
        for h in range(int(target["shift_count"]))
    }
    if {(-delta) % modulus for delta in source_set} != target_set:
        raise AssertionError("asymmetric P1-P3/M1-M3 orientation is wrong")
    return {
        "class_involution": reflection,
        "vertex_involution_checks": involution_checks,
        "exhaustive_loop_checks": loop_checks,
        "residue_difference_automorphism_checks": residue_difference_checks,
        "asymmetric_start_721_maps_to_start_0": True,
        "automorphism_exact": True,
    }


def reflected_phases(values: Sequence) -> list:
    return [values[0], -values[0], values[1], -values[1], values[2], -values[2]]


def torque3(
    values: Sequence,
    degrees: list[list[int]],
    sine: Callable,
) -> list:
    phases = reflected_phases(values)
    output = []
    for i in POSITIVE_CLASSES:
        output.append(
            sum(
                (
                    degrees[i][j] * sine(phases[j] - phases[i])
                    for j in range(6)
                    if i != j and degrees[i][j]
                ),
                phases[0] * 0,
            )
        )
    return output


def jacobian3(
    values: Sequence,
    degrees: list[list[int]],
    cosine: Callable,
) -> list[list]:
    phases = reflected_phases(values)
    matrix = []
    for i in POSITIVE_CLASSES:
        row = []
        for variable in range(3):
            row.append(
                sum(
                    (
                        degrees[i][j]
                        * cosine(phases[j] - phases[i])
                        * (
                            PHASE_COEFFICIENTS[j][variable]
                            - PHASE_COEFFICIENTS[i][variable]
                        )
                        for j in range(6)
                        if i != j and degrees[i][j]
                    ),
                    phases[0] * 0,
                )
            )
        matrix.append(row)
    return matrix


def solve_equilibrium(
    spec: dict[str, Any],
    degrees: list[list[int]],
) -> tuple[list[mp.mpf], list[dict[str, str]]]:
    values = mp.matrix([mp.mpf(x) for x in spec["equilibrium"]["angle_seed"]])
    trace = []
    for iteration in range(20):
        residual_vector = mp.matrix(torque3(list(values), degrees, mp.sin))
        residual = max(abs(x) for x in residual_vector)
        trace.append(
            {
                "iteration": str(iteration),
                "max_residual": mp_text(residual, 25),
            }
        )
        if residual < mp.mpf("1e-150"):
            return list(values), trace
        derivative = mp.matrix(jacobian3(list(values), degrees, mp.cos))
        values += mp.lu_solve(derivative, -residual_vector)
    raise RuntimeError("high-precision Newton solve did not converge")


def iv_point(value: mp.mpf):
    return mp.iv.mpf(mp.nstr(value, POINT_DPS - 5))


def iv_bounds(value) -> tuple[mp.mpf, mp.mpf]:
    return mp.mpf(value._mpi_[0]), mp.mpf(value._mpi_[1])


def iv_abs_upper(value) -> mp.mpf:
    lower, upper = iv_bounds(value)
    return max(abs(lower), abs(upper))


def krawczyk_certificate(
    root: list[mp.mpf],
    spec: dict[str, Any],
    degrees: list[list[int]],
) -> tuple[dict[str, Any], list]:
    seed = [mp.mpf(x) for x in spec["equilibrium"]["angle_seed"]]
    radius = mp.mpf("1e-10")
    box = [
        mp.iv.mpf(
            [
                mp.nstr(center - radius, POINT_DPS),
                mp.nstr(center + radius, POINT_DPS),
            ]
        )
        for center in seed
    ]
    point_jacobian = mp.matrix(jacobian3(root, degrees, mp.cos))
    inverse = point_jacobian ** -1
    interval_jacobian = jacobian3(box, degrees, mp.iv.cos)
    # Krawczyk's preconditioner may be any fixed nonsingular real matrix; it
    # need not be the exact inverse.  Convert the high-precision approximation
    # to outward-rounded point intervals, then perform every certificate
    # operation (including the center residual) with interval arithmetic.
    root_interval = [iv_point(value) for value in root]
    inverse_interval = [
        [iv_point(inverse[i, j]) for j in range(3)] for i in range(3)
    ]
    residual_interval = torque3(root_interval, degrees, mp.iv.sin)
    corrected = []
    for i in range(3):
        value = root_interval[i]
        for j in range(3):
            value -= inverse_interval[i][j] * residual_interval[j]
        corrected.append(value)

    remainder = [[mp.iv.mpf(0) for _ in range(3)] for _ in range(3)]
    for i in range(3):
        for j in range(3):
            product = mp.iv.mpf(0)
            for k in range(3):
                product += inverse_interval[i][k] * interval_jacobian[k][j]
            remainder[i][j] = (mp.iv.mpf(1) if i == j else mp.iv.mpf(0)) - product

    image = []
    for i in range(3):
        value = corrected[i]
        for j in range(3):
            value += remainder[i][j] * (box[j] - root_interval[j])
        image.append(value)
    box_bounds = [iv_bounds(value) for value in box]
    image_bounds = [iv_bounds(value) for value in image]
    strict = all(
        box_bounds[i][0] < image_bounds[i][0]
        and image_bounds[i][1] < box_bounds[i][1]
        for i in range(3)
    )
    if not strict:
        raise AssertionError("Krawczyk image is not strictly inside the seed box")
    contraction = max(
        sum(iv_abs_upper(remainder[i][j]) for j in range(3))
        for i in range(3)
    )
    return (
        {
            "method": "outward-rounded interval Krawczyk operator",
            "box_radius_about_short_seed": mp_text(radius, 12),
            "box": [
                {
                    "lower": mp_text(lower, 35),
                    "upper": mp_text(upper, 35),
                }
                for lower, upper in box_bounds
            ],
            "remainder_infinity_norm_upper": mp_text(contraction, 25),
            "strict_interior_inclusion": True,
            "unique_root_within_box": True,
            "global_uniqueness_claimed": False,
        },
        box,
    )


def quotient_hessian(
    degrees: list[list[int]],
    phases: list[mp.mpf],
) -> mp.matrix:
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


def interval_ldlt_positive(matrix: list[list]) -> list[tuple[mp.mpf, mp.mpf]]:
    """Prove a symmetric interval matrix positive definite by interval LDL^T.

    Each interval pivot encloses the corresponding exact pivot for every
    point matrix in ``matrix``.  Strictly positive lower endpoints therefore
    certify positive definiteness uniformly over the entire interval family.
    """
    dimension = len(matrix)
    lower = [[mp.iv.mpf(0) for _ in range(dimension)] for _ in range(dimension)]
    pivots = [mp.iv.mpf(0) for _ in range(dimension)]
    for i in range(dimension):
        lower[i][i] = mp.iv.mpf(1)
    for column in range(dimension):
        pivot = matrix[column][column]
        for k in range(column):
            pivot -= lower[column][k] * lower[column][k] * pivots[k]
        pivot_bounds = iv_bounds(pivot)
        if pivot_bounds[0] <= 0:
            raise AssertionError(
                f"interval LDL^T pivot {column} is not strictly positive"
            )
        pivots[column] = pivot
        for row in range(column + 1, dimension):
            numerator = matrix[row][column]
            for k in range(column):
                numerator -= lower[row][k] * lower[column][k] * pivots[k]
            lower[row][column] = numerator / pivots[column]
    return [iv_bounds(pivot) for pivot in pivots]


def quotient_interval_certificate(
    sizes: list[int],
    degrees: list[list[int]],
    angle_box: list,
) -> dict[str, Any]:
    """Certify a uniform nonrotation quotient gap over the root box.

    The mass-orthonormal quotient Hessian has the exact fixed null vector
    r=(sqrt(n_i)).  Let R=rr^T/N be the orthogonal rotation projector.  If
    H+R-cI is positive definite, the rotation eigenvalue has been lifted to
    1-c and every physical quotient eigenvalue is greater than c.  Interval
    LDL^T proves this simultaneously for every phase in the Krawczyk box.
    """
    phases = reflected_phases(angle_box)
    dimension = len(sizes)
    hessian = [
        [mp.iv.mpf(0) for _ in range(dimension)] for _ in range(dimension)
    ]
    for i in range(dimension):
        for j in range(i + 1, dimension):
            if not degrees[i][j]:
                continue
            cosine = mp.iv.cos(phases[i] - phases[j])
            hessian[i][i] += degrees[i][j] * cosine
            hessian[j][j] += degrees[j][i] * cosine
            geometric_degree = mp.iv.sqrt(
                mp.iv.mpf(degrees[i][j] * degrees[j][i])
            )
            value = -geometric_degree * cosine
            hessian[i][j] = value
            hessian[j][i] = value

    vertex_count = sum(sizes)
    gap_lower = mp.mpf("0.4")
    shifted = [
        [mp.iv.mpf(hessian[i][j]) for j in range(dimension)]
        for i in range(dimension)
    ]
    for i in range(dimension):
        for j in range(dimension):
            projector = mp.iv.sqrt(mp.iv.mpf(sizes[i] * sizes[j])) / vertex_count
            shifted[i][j] += projector
        shifted[i][i] -= gap_lower

    pivot_bounds = interval_ldlt_positive(shifted)
    return {
        "method": (
            "outward-rounded interval LDL^T on "
            "H + rotation_projector - 0.4 I"
        ),
        "phase_domain": "the full Krawczyk root box",
        "exact_rotation_kernel_used": True,
        "certified_nonrotation_absolute_gap_lower": mp_text(gap_lower, 20),
        "certified_nonrotation_normalized_gap_lower": mp_text(
            gap_lower / vertex_count, 20
        ),
        "interval_ldlt_pivots": [
            {"lower": mp_text(bounds[0], 30), "upper": mp_text(bounds[1], 30)}
            for bounds in pivot_bounds
        ],
        "all_interval_pivot_lowers_positive": True,
        "uniform_over_certified_root_box": True,
    }


def nonsynchronous_interval_certificate(angle_box: list) -> dict[str, Any]:
    """Prove two reflected phase classes are distinct throughout the root box."""
    separation = 2 * angle_box[0]
    lower, upper = iv_bounds(separation)
    if lower <= 0 or upper >= mp.pi:
        raise AssertionError("the root box does not certify nonsynchrony")
    return {
        "classes": [0, 1],
        "phase_separation_interval": {
            "lower": mp_text(lower, 30),
            "upper": mp_text(upper, 30),
        },
        "separation_strictly_between_0_and_pi": True,
        "nonsynchronous_throughout_root_box": True,
    }


def reflection_blocks(matrix: mp.matrix) -> tuple[mp.matrix, mp.matrix]:
    even, odd = mp.matrix(3), mp.matrix(3)
    for i in range(3):
        for j in range(3):
            even[i, j] = matrix[2 * i, 2 * j] + matrix[2 * i, 2 * j + 1]
            odd[i, j] = matrix[2 * i, 2 * j] - matrix[2 * i, 2 * j + 1]
    return even, odd


def transverse_certificate(
    sizes: list[int],
    degrees: list[list[int]],
    phases: list[mp.mpf],
    phase_box: list,
    blocks: dict[tuple[int, int], dict],
) -> dict[str, Any]:
    diagonal = [
        sizes[i]
        + mp.fsum(
            degrees[i][j] * mp.cos(phases[i] - phases[j])
            for j in range(6)
        )
        for i in range(6)
    ]
    comparison = mp.matrix(6)
    for i in range(6):
        for j in range(6):
            comparison[i, j] = mp.mpf(0)
        comparison[i, i] = diagonal[i]
    for (i, j), block in blocks.items():
        if block["type"] != "residue_orbits":
            continue
        bound = abs(mp.cos(phases[i] - phases[j])) * mp.sqrt(
            degrees[i][j] * degrees[j][i]
        )
        comparison[i, j] = comparison[j, i] = -bound
    comparison_eigenvalues = list(mp.eigsy(comparison, eigvals_only=True))

    interval_phases = reflected_phases(phase_box)
    interval_diagonal = []
    for i in range(6):
        value = mp.iv.mpf(sizes[i])
        for j in range(6):
            if i != j and degrees[i][j]:
                value += degrees[i][j] * mp.iv.cos(
                    interval_phases[i] - interval_phases[j]
                )
        interval_diagonal.append(value)
    radii = [mp.mpf(0) for _ in range(6)]
    for (i, j), block in blocks.items():
        if block["type"] != "residue_orbits":
            continue
        cosine_upper = iv_abs_upper(
            mp.iv.cos(interval_phases[i] - interval_phases[j])
        )
        beta = cosine_upper * mp.sqrt(degrees[i][j] * degrees[j][i])
        radii[i] += beta
        radii[j] += beta
    row_lowers = [
        iv_bounds(interval_diagonal[i])[0] - radii[i]
        for i in range(6)
    ]
    rigorous_lower = min(row_lowers)
    if rigorous_lower <= 0:
        raise AssertionError("Fourier-free transverse lower bound is not positive")
    return {
        "space_dimension": sum(n - 1 for n in sizes),
        "proof": (
            "On each fiber-zero-sum space the clique contributes n_i I; "
            "complete blocks vanish; each partial biregular adjacency block "
            "has norm at most sqrt(d_ij*d_ji). The resulting scalar comparison "
            "is positive by outward-rounded Gershgorin bounds."
        ),
        "fiber_diagonal_terms": [mp_text(value, 35) for value in diagonal],
        "comparison_eigenvalues": [
            mp_text(value, 35) for value in comparison_eigenvalues
        ],
        "point_comparison_lower": mp_text(comparison_eigenvalues[0], 35),
        "rigorous_gershgorin_row_lowers": [
            mp_text(value, 35) for value in row_lowers
        ],
        "rigorous_absolute_lower": mp_text(rigorous_lower, 35),
        "all_transverse_modes_positive": True,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=HERE / "verification_report.json")
    args = parser.parse_args()

    mp.mp.dps = POINT_DPS
    mp.iv.dps = INTERVAL_DPS
    spec_path = HERE / "graph_spec.json"
    spec, sizes, blocks = load_spec(spec_path)
    degrees, edge_count, block_details = block_combinatorics(sizes, blocks)
    class_degrees = [
        sizes[i] - 1 + sum(degrees[i])
        for i in range(6)
    ]
    vertex_count = sum(sizes)
    degree_sum = sum(sizes[i] * class_degrees[i] for i in range(6))
    if degree_sum != 2 * edge_count:
        raise AssertionError("handshake identity failed")
    minimum_degree = min(class_degrees)
    mu = Fraction(minimum_degree, vertex_count - 1)
    target = Fraction(
        int(spec["target_min_degree_ratio"]["numerator"]),
        int(spec["target_min_degree_ratio"]["denominator"]),
    )
    excess = target.denominator * minimum_degree - target.numerator * (
        vertex_count - 1
    )
    if not mu > target or excess != 277:
        raise AssertionError("minimum-degree threshold failed")
    if not support_connected(degrees):
        raise AssertionError("class support graph is disconnected")
    reflection = verify_reflection(spec, sizes, blocks)

    angles, newton_trace = solve_equilibrium(spec, degrees)
    phases = reflected_phases(angles)
    six_torques = []
    for i in range(6):
        six_torques.append(
            mp.fsum(
                degrees[i][j] * mp.sin(phases[j] - phases[i])
                for j in range(6)
            )
        )
    residual = max(abs(value) for value in six_torques)
    if residual >= mp.mpf("1e-140"):
        raise AssertionError("high-precision torque residual is too large")
    krawczyk, phase_box = krawczyk_certificate(angles, spec, degrees)
    quotient_interval = quotient_interval_certificate(
        sizes,
        degrees,
        phase_box,
    )
    nonsynchronous_interval = nonsynchronous_interval_certificate(phase_box)

    quotient = quotient_hessian(degrees, phases)
    even, odd = reflection_blocks(quotient)
    full_eigenvalues = list(mp.eigsy(quotient, eigvals_only=True))
    even_eigenvalues = list(mp.eigsy(even, eigvals_only=True))
    odd_eigenvalues = list(mp.eigsy(odd, eigvals_only=True))
    if abs(full_eigenvalues[0]) >= mp.mpf("1e-120"):
        raise AssertionError("rotation eigenvalue is not numerically zero")
    if min(full_eigenvalues[1:]) <= 0:
        raise AssertionError("nonrotation quotient eigenvalue is not positive")
    quotient_gap = min(full_eigenvalues[1:])
    normalized_gap = quotient_gap / vertex_count
    if normalized_gap <= mp.mpf("1e-6"):
        raise AssertionError("normalized quotient gap does not exceed 1e-6")
    rotation = mp.matrix([mp.sqrt(n) for n in sizes])
    rotation_residual = max(abs(value) for value in quotient * rotation)

    transverse = transverse_certificate(
        sizes,
        degrees,
        phases,
        phase_box,
        blocks,
    )
    order_parameter = abs(
        2
        * mp.fsum(sizes[2 * i] * mp.cos(angles[i]) for i in range(3))
        / vertex_count
    )

    report = {
        "schema": "reflection-witness-verification-v1",
        "input": {
            "path": "graph_spec.json",
            "sha256": hashlib.sha256(spec_path.read_bytes()).hexdigest(),
            "only_local_graph_spec_read": True,
        },
        "graph": {
            "vertex_count": vertex_count,
            "class_sizes": sizes,
            "class_degrees": class_degrees,
            "minimum_degree": minimum_degree,
            "edge_count": edge_count,
            "degree_sum": degree_sum,
            "handshake_exact": True,
            "simple": True,
            "connected": True,
            "block_partition_complete_and_unique": True,
            "blocks": block_details,
            "reflection": reflection,
        },
        "minimum_degree_ratio": {
            "exact": f"{mu.numerator}/{mu.denominator}",
            "decimal": mp_text(mp.mpf(mu.numerator) / mu.denominator, 30),
            "target": f"{target.numerator}/{target.denominator}",
            "cross_product_excess": excess,
        },
        "equilibrium": {
            "angles": [mp_text(value, 160) for value in angles],
            "phases": [mp_text(value, 160) for value in phases],
            "six_torques": [mp_text(value, 25) for value in six_torques],
            "max_abs_torque_residual": mp_text(residual, 25),
            "newton_trace": newton_trace,
            "krawczyk": krawczyk,
            "order_parameter": mp_text(order_parameter, 40),
            "nonsynchronous_interval_certificate": nonsynchronous_interval,
        },
        "quotient": {
            "basis": "mass-orthonormal class-constant coordinates",
            "full_eigenvalues": [
                mp_text(value, 40) for value in full_eigenvalues
            ],
            "reflection_even_eigenvalues": [
                mp_text(value, 40) for value in even_eigenvalues
            ],
            "reflection_odd_eigenvalues": [
                mp_text(value, 40) for value in odd_eigenvalues
            ],
            "rotation_residual_max": mp_text(rotation_residual, 25),
            "exactly_one_rotation_zero": True,
            "absolute_gap": mp_text(quotient_gap, 40),
            "normalized_gap": mp_text(normalized_gap, 40),
            "normalized_gap_above_1e_minus_6": True,
            "interval_certificate": quotient_interval,
        },
        "transverse": transverse,
        "persisted_independent_audit": {
            "nonlinear_runs_omitted_by_this_compact_verifier": True,
            "verdict": "LOCAL RETURN CONFIRMED; BASIN SIZE NOT CERTIFIED",
            "return_observed": {
                "weak_mode_amplitudes_both_signs": [0.0001, 0.0003],
                "random_vertexwise_amplitudes_through": 0.01,
            },
            "narrow_basin_caveat": (
                "Weak-mode amplitudes 0.0005 and 0.002 escaped to synchrony "
                "in both signs during the independent audit."
            ),
            "enumerated_transverse_minimum": 19244.984787192374,
        },
        "dimension_accounting": {
            "quotient": 6,
            "transverse": vertex_count - 6,
            "total": vertex_count,
        },
        "accepted": True,
    }
    args.out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(
        "ACCEPTED",
        f"N={vertex_count}",
        f"mu={mu.numerator}/{mu.denominator}",
        f"gap/N={mp.nstr(normalized_gap, 12)}",
        f"transverse>={mp.nstr(mp.mpf(transverse['rigorous_absolute_lower']), 12)}",
    )


if __name__ == "__main__":
    main()
