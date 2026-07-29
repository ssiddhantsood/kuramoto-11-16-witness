#!/usr/bin/env python3
"""Independent high-precision replay of one v2 reflection graph spec.

This file imports no construction or weighted-optimization module.  It reads
only the declarative graph specification, reconstructs exact graph arithmetic,
proves the reflection automorphism algebraically for every vertex label,
certifies a unique equilibrium root in a box, proves both quotient sectors
strictly above 1e-6, proves every transverse mode positive by a block-operator
comparison, and checks local nonlinear return on the exact equitable flow.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Callable, Sequence

import mpmath as mp
import numpy as np
from scipy.integrate import solve_ivp


HERE = Path(__file__).resolve().parent
DEFAULT_OUTPUT = HERE.parent / "artifacts" / "finite" / "replay_certificate.json"
POINT_DPS = 190
INTERVAL_DPS = 150
GAP_FLOOR_TEXT = "0.000001"


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def interval_bounds(value) -> tuple[mp.mpf, mp.mpf]:
    return mp.mpf(value._mpi_[0]), mp.mpf(value._mpi_[1])


def interval_abs_upper(value) -> mp.mpf:
    lower, upper = interval_bounds(value)
    return max(abs(lower), abs(upper))


def matrix_multiply(left: Sequence[Sequence], right: Sequence[Sequence]):
    rows = len(left)
    inner = len(right)
    columns = len(right[0])
    return [
        [
            sum(
                (left[i][k] * right[k][j] for k in range(inner)),
                left[i][0] * 0,
            )
            for j in range(columns)
        ]
        for i in range(rows)
    ]


def transpose(matrix: Sequence[Sequence]):
    return [list(row) for row in zip(*matrix)]


def parse_spec(path: Path) -> dict[str, Any]:
    spec = json.loads(path.read_text())
    if spec.get("schema") != "reflection-paired-circulant-witness-v2":
        raise ValueError("unsupported graph spec schema")
    return spec


def parse_graph(
    spec: dict[str, Any],
) -> tuple[list[int], dict[tuple[int, int], dict[str, Any]]]:
    classes = sorted(spec["classes"], key=lambda item: int(item["id"]))
    sizes = [int(item["size"]) for item in classes]
    if len(sizes) % 2 or any(size <= 1 for size in sizes):
        raise ValueError("expected an even number of nontrivial fibers")
    if any(sizes[2 * i] != sizes[2 * i + 1] for i in range(len(sizes) // 2)):
        raise ValueError("reflection-paired fibers have unequal sizes")
    if sum(sizes) != int(spec["vertex_count"]):
        raise ValueError("class sizes do not sum to N")
    blocks = {}
    for raw in spec["cross_blocks"]:
        block = dict(raw)
        i, j = map(int, block["classes"])
        if not 0 <= i < j < len(sizes) or (i, j) in blocks:
            raise ValueError(f"invalid or repeated block {(i, j)}")
        kind = block["type"]
        if kind not in {"absent", "complete", "residue_orbits"}:
            raise ValueError(f"unsupported block kind {kind}")
        if kind == "residue_orbits":
            modulus = int(block["modulus"])
            count = int(block["shift_count"])
            start = int(block["shift_start"])
            if modulus != math.gcd(sizes[i], sizes[j]):
                raise ValueError("residue modulus is not the class-size gcd")
            if not (0 <= count <= modulus and 0 <= start < modulus):
                raise ValueError("invalid residue interval")
        blocks[(i, j)] = block
    expected = {
        (i, j)
        for i in range(len(sizes))
        for j in range(i + 1, len(sizes))
    }
    if set(blocks) != expected:
        raise ValueError("cross-block partition is incomplete")
    return sizes, blocks


def exact_block_data(
    sizes: list[int],
    blocks: dict[tuple[int, int], dict[str, Any]],
) -> tuple[list[list[int]], list[list[int]], int]:
    count = len(sizes)
    degrees = [[0 for _ in range(count)] for _ in range(count)]
    edges = [[0 for _ in range(count)] for _ in range(count)]
    total_edges = sum(size * (size - 1) // 2 for size in sizes)
    for (i, j), block in blocks.items():
        if block["type"] == "absent":
            edge_count = 0
        elif block["type"] == "complete":
            edge_count = sizes[i] * sizes[j]
        else:
            edge_count = (
                sizes[i]
                * sizes[j]
                // int(block["modulus"])
                * int(block["shift_count"])
            )
        if edge_count % sizes[i] or edge_count % sizes[j]:
            raise AssertionError("block is not biregular")
        degrees[i][j] = edge_count // sizes[i]
        degrees[j][i] = edge_count // sizes[j]
        edges[i][j] = edges[j][i] = edge_count
        total_edges += edge_count
    return degrees, edges, total_edges


def reflection_proof(
    sizes: list[int],
    blocks: dict[tuple[int, int], dict[str, Any]],
) -> dict[str, Any]:
    """Prove A(c,u;d,v)=A(rho(c,u);rho(d,v)) for all labels.

    For an interval S={start,...,start+k-1}, endpoint negation maps the
    canonical residue difference to -S when reflected class order is
    preserved and to S when it reverses.  The exact start relation checked
    below is therefore a universal vertexwise proof, not random sampling.
    """

    involution = [index ^ 1 for index in range(len(sizes))]
    rows = []
    for (i, j), source in sorted(blocks.items()):
        ri, rj = involution[i], involution[j]
        order_preserved = ri < rj
        target_key = (ri, rj) if order_preserved else (rj, ri)
        target = blocks[target_key]
        if source["type"] != target["type"]:
            raise AssertionError("reflection changes block type")
        row: dict[str, Any] = {
            "source": [i, j],
            "target": list(target_key),
            "type": source["type"],
            "canonical_order_preserved": order_preserved,
        }
        if source["type"] == "residue_orbits":
            modulus = int(source["modulus"])
            count = int(source["shift_count"])
            start = int(source["shift_start"])
            expected = (
                -(start + count - 1) if order_preserved else start
            ) % modulus
            if (
                int(target["modulus"]) != modulus
                or int(target["shift_count"]) != count
                or int(target["shift_start"]) != expected
            ):
                raise AssertionError("reflected interval orientation mismatch")
            row.update(
                {
                    "modulus": modulus,
                    "count": count,
                    "source_start": start,
                    "required_target_start": expected,
                    "universal_residue_identity": True,
                }
            )
        rows.append(row)
    return {
        "class_involution": involution,
        "local_involution": "u -> -u mod n_c",
        "involution_for_every_vertex": True,
        "adjacency_preserved_for_every_vertex_pair": True,
        "proof_method": (
            "Exact block partition plus universal modular interval identity; "
            "covers all local labels u,v without enumerating N^2 pairs."
        ),
        "block_orbits": rows,
    }


def support_connected(
    sizes: list[int],
    blocks: dict[tuple[int, int], dict[str, Any]],
) -> bool:
    adjacency = [[] for _ in sizes]
    for (i, j), block in blocks.items():
        if block["type"] != "absent":
            adjacency[i].append(j)
            adjacency[j].append(i)
    seen = {0}
    pending = [0]
    while pending:
        i = pending.pop()
        for j in adjacency[i]:
            if j not in seen:
                seen.add(j)
                pending.append(j)
    return len(seen) == len(sizes)


def reflected_phases(angles: Sequence):
    output = []
    for value in angles:
        output.extend((value, -value))
    return output


def phase_coefficients(m: int) -> list[list[int]]:
    rows = []
    for pair in range(m):
        positive = [0] * m
        negative = [0] * m
        positive[pair] = 1
        negative[pair] = -1
        rows.extend((positive, negative))
    return rows


def torque_system(
    angles: Sequence,
    degrees: list[list[int]],
    total: int,
    sine: Callable,
):
    phases = reflected_phases(angles)
    output = []
    for pair in range(len(angles)):
        i = 2 * pair
        value = phases[0] * 0
        for j in range(len(phases)):
            if i != j and degrees[i][j]:
                value += (
                    mp.mpf(degrees[i][j])
                    / total
                    * sine(phases[i] - phases[j])
                )
        output.append(value)
    return output


def torque_jacobian(
    angles: Sequence,
    degrees: list[list[int]],
    total: int,
    cosine: Callable,
):
    m = len(angles)
    phases = reflected_phases(angles)
    coefficients = phase_coefficients(m)
    context = mp.iv if hasattr(angles[0], "_mpi_") else mp
    output = [[angles[0] * 0 for _ in range(m)] for _ in range(m)]
    for pair in range(m):
        i = 2 * pair
        for j in range(2 * m):
            if i == j or not degrees[i][j]:
                continue
            factor = (
                context.mpf(degrees[i][j])
                / total
                * cosine(phases[i] - phases[j])
            )
            for variable in range(m):
                output[pair][variable] += factor * (
                    coefficients[i][variable]
                    - coefficients[j][variable]
                )
    return output


def solve_root(
    spec: dict[str, Any],
    degrees: list[list[int]],
    total: int,
) -> list[mp.mpf]:
    seeds = [mp.mpf(value) for value in spec["equilibrium"]["angle_seed"]]
    functions = tuple(
        (
            lambda index: lambda *values: torque_system(
                values, degrees, total, mp.sin
            )[index]
        )(index)
        for index in range(len(seeds))
    )
    return list(
        mp.findroot(
            functions,
            tuple(seeds),
            solver="mdnewton",
            tol=mp.mpf("1e-165"),
            maxsteps=100,
        )
    )


def krawczyk_certificate(
    root_angles: list[mp.mpf],
    degrees: list[list[int]],
    total: int,
) -> tuple[dict[str, Any], list[Any]]:
    m = len(root_angles)
    radius = mp.mpf("1e-32")
    box = [
        mp.iv.mpf(
            [
                mp.nstr(value - radius, POINT_DPS),
                mp.nstr(value + radius, POINT_DPS),
            ]
        )
        for value in root_angles
    ]
    f_center = mp.matrix(
        torque_system(root_angles, degrees, total, mp.sin)
    )
    jacobian = mp.matrix(
        torque_jacobian(root_angles, degrees, total, mp.cos)
    )
    inverse = jacobian**-1
    interval_jacobian = torque_jacobian(
        box, degrees, total, mp.iv.cos
    )
    guard = mp.mpf("1e-150")
    remainder = [[mp.iv.mpf(0) for _ in range(m)] for _ in range(m)]
    for i in range(m):
        for j in range(m):
            value = mp.iv.mpf(1 if i == j else 0)
            for k in range(m):
                inverse_interval = mp.iv.mpf(
                    [
                        mp.nstr(inverse[i, k] - guard, POINT_DPS),
                        mp.nstr(inverse[i, k] + guard, POINT_DPS),
                    ]
                )
                value -= inverse_interval * interval_jacobian[k][j]
            remainder[i][j] = value
    center = mp.matrix(root_angles) - inverse * f_center
    delta = mp.iv.mpf([mp.nstr(-radius, POINT_DPS), mp.nstr(radius, POINT_DPS)])
    inclusions = []
    strict = True
    for i in range(m):
        value = mp.iv.mpf(
            [
                mp.nstr(center[i] - guard, POINT_DPS),
                mp.nstr(center[i] + guard, POINT_DPS),
            ]
        )
        for j in range(m):
            value += remainder[i][j] * delta
        lower, upper = interval_bounds(value)
        box_lower, box_upper = interval_bounds(box[i])
        contained = box_lower < lower and upper < box_upper
        strict = strict and contained
        inclusions.append(
            {
                "image_lower": mp.nstr(lower, 50),
                "image_upper": mp.nstr(upper, 50),
                "box_lower": mp.nstr(box_lower, 50),
                "box_upper": mp.nstr(box_upper, 50),
                "strictly_contained": contained,
            }
        )
    if not strict:
        raise AssertionError("Krawczyk inclusion failed")
    return (
        {
            "radius": mp.nstr(radius, 20),
            "strict_interior_inclusion": strict,
            "unique_root_in_box": strict,
            "global_uniqueness_claimed": False,
            "inclusions": inclusions,
        },
        box,
    )


def interval_ldl(matrix: list[list[Any]], shift: mp.mpf) -> dict[str, Any]:
    n = len(matrix)
    lower = [[mp.iv.mpf(0) for _ in range(n)] for _ in range(n)]
    diagonal = [mp.iv.mpf(0) for _ in range(n)]
    for i in range(n):
        lower[i][i] = mp.iv.mpf(1)
        pivot = matrix[i][i] - mp.iv.mpf(mp.nstr(shift, POINT_DPS))
        for k in range(i):
            pivot -= lower[i][k] * lower[i][k] * diagonal[k]
        pivot_lower, _ = interval_bounds(pivot)
        if pivot_lower <= 0:
            raise AssertionError(f"interval LDL pivot {i} is not positive")
        diagonal[i] = pivot
        for j in range(i + 1, n):
            numerator = matrix[j][i]
            for k in range(i):
                numerator -= (
                    lower[j][k] * lower[i][k] * diagonal[k]
                )
            lower[j][i] = numerator / diagonal[i]
    return {
        "shift": mp.nstr(shift, 30),
        "pivot_intervals": [str(value) for value in diagonal],
        "all_pivots_positive": True,
    }


def quotient_matrices(
    sizes: list[int],
    degrees: list[list[int]],
    angles: Sequence,
    total: int,
    interval: bool,
):
    context = mp.iv if interval else mp
    m = len(angles)
    phases = reflected_phases(angles)
    quotient = [
        [context.mpf(0) for _ in range(2 * m)] for _ in range(2 * m)
    ]
    for i in range(2 * m):
        for j in range(2 * m):
            if i == j or not degrees[i][j]:
                continue
            cosine = context.cos(phases[i] - phases[j])
            quotient[i][i] += (
                context.mpf(degrees[i][j]) / total * cosine
            )
            quotient[i][j] = (
                -context.sqrt(
                    context.mpf(degrees[i][j])
                    * context.mpf(degrees[j][i])
                )
                / total
                * cosine
            )
    root_two = context.sqrt(context.mpf(2))
    even = [[context.mpf(0) for _ in range(m)] for _ in range(m)]
    odd = [[context.mpf(0) for _ in range(m)] for _ in range(m)]
    for i in range(m):
        pi, mi = 2 * i, 2 * i + 1
        for j in range(m):
            pj, mj = 2 * j, 2 * j + 1
            even[i][j] = (
                quotient[pi][pj]
                + quotient[pi][mj]
                + quotient[mi][pj]
                + quotient[mi][mj]
            ) / 2
            odd[i][j] = (
                quotient[pi][pj]
                - quotient[pi][mj]
                - quotient[mi][pj]
                + quotient[mi][mj]
            ) / 2
    rotation = [
        context.sqrt(context.mpf(2 * sizes[2 * i]) / total)
        for i in range(m)
    ]
    lifted = [
        [
            even[i][j]
            + context.mpf(2) * rotation[i] * rotation[j]
            for j in range(m)
        ]
        for i in range(m)
    ]
    return quotient, even, odd, lifted


def quotient_certificate(
    sizes: list[int],
    degrees: list[list[int]],
    root_angles: list[mp.mpf],
    angle_box: list[Any],
    total: int,
) -> dict[str, Any]:
    gap_floor = mp.mpf(GAP_FLOOR_TEXT)
    quotient, even, odd, lifted = quotient_matrices(
        sizes, degrees, root_angles, total, interval=False
    )
    full_values = list(mp.eigsy(mp.matrix(quotient), eigvals_only=True))
    even_values = list(mp.eigsy(mp.matrix(even), eigvals_only=True))
    odd_values = list(mp.eigsy(mp.matrix(odd), eigvals_only=True))
    lifted_values = list(mp.eigsy(mp.matrix(lifted), eigvals_only=True))
    gap = min(min(lifted_values), min(odd_values))
    _, _, odd_interval, lifted_interval = quotient_matrices(
        sizes, degrees, angle_box, total, interval=True
    )
    even_ldl = interval_ldl(lifted_interval, gap_floor)
    odd_ldl = interval_ldl(odd_interval, gap_floor)
    rotation_zero = min(full_values, key=abs)
    exactly_one_zero = bool(
        abs(rotation_zero) < mp.mpf("1e-145")
        and gap > gap_floor
    )
    if not exactly_one_zero:
        raise AssertionError("quotient does not have exactly one zero")
    return {
        "full_eigenvalues_normalized": [
            mp.nstr(value, 70) for value in full_values
        ],
        "even_raw_eigenvalues_normalized": [
            mp.nstr(value, 70) for value in even_values
        ],
        "lifted_even_eigenvalues_normalized": [
            mp.nstr(value, 70) for value in lifted_values
        ],
        "odd_eigenvalues_normalized": [
            mp.nstr(value, 70) for value in odd_values
        ],
        "normalized_gap": mp.nstr(gap, 70),
        "structural_rotation_zero": mp.nstr(rotation_zero, 40),
        "exactly_one_zero": exactly_one_zero,
        "interval_shift": mp.nstr(gap_floor, 20),
        "interval_lifted_even_ldl": even_ldl,
        "interval_odd_ldl": odd_ldl,
        "rigorously_above_1e_minus_6": True,
        "_gap": gap,
    }


def residue_norm_upper_interval(
    block: dict[str, Any],
    sizes: list[int],
):
    i, j = map(int, block["classes"])
    modulus = int(block["modulus"])
    count = int(block["shift_count"])
    a = sizes[i] // modulus
    b = sizes[j] // modulus
    missing = min(count, modulus - count)
    if missing == 0:
        return mp.iv.mpf(0), "zero"
    sine = mp.iv.sin(mp.iv.pi / modulus)
    sine_lower, _ = interval_bounds(sine)
    reciprocal_upper = mp.mpf(1) / sine_lower
    scalar = min(mp.mpf(missing), reciprocal_upper)
    guard = mp.mpf("1e-135")
    upper = mp.sqrt(mp.mpf(a) * b) * scalar + guard
    return mp.iv.mpf([0, mp.nstr(upper, POINT_DPS)]), (
        "sqrt(a*b)*min(min(k,g-k),1/sin(pi/g))"
    )


def transverse_certificate(
    sizes: list[int],
    blocks: dict[tuple[int, int], dict[str, Any]],
    degrees: list[list[int]],
    angle_box: list[Any],
    root_angles: list[mp.mpf],
    total: int,
) -> dict[str, Any]:
    phases = reflected_phases(angle_box)
    count = len(sizes)
    beta = []
    for i in range(count):
        value = mp.iv.mpf(mp.nstr(mp.mpf(sizes[i]) / total, POINT_DPS))
        for j in range(count):
            if i != j and degrees[i][j]:
                value += (
                    mp.iv.mpf(
                        mp.nstr(mp.mpf(degrees[i][j]) / total, POINT_DPS)
                    )
                    * mp.iv.cos(phases[i] - phases[j])
                )
        beta.append(value)
    comparison = [
        [mp.iv.mpf(0) for _ in range(count)] for _ in range(count)
    ]
    for i in range(count):
        comparison[i][i] = beta[i]
    block_rows = []
    for (i, j), block in sorted(blocks.items()):
        if block["type"] != "residue_orbits":
            continue
        norm, formula = residue_norm_upper_interval(block, sizes)
        cosine_upper = interval_abs_upper(
            mp.iv.cos(phases[i] - phases[j])
        )
        _, norm_upper = interval_bounds(norm)
        weighted_upper = cosine_upper * norm_upper / total
        value = mp.iv.mpf(
            [
                mp.nstr(-weighted_upper, POINT_DPS),
                mp.nstr(-weighted_upper, POINT_DPS),
            ]
        )
        comparison[i][j] = comparison[j][i] = value
        block_rows.append(
            {
                "classes": [i, j],
                "modulus": int(block["modulus"]),
                "count": int(block["shift_count"]),
                "norm_formula": formula,
                "zero_sum_norm_upper_absolute": mp.nstr(
                    norm_upper, 50
                ),
                "weighted_norm_upper_normalized": mp.nstr(
                    weighted_upper, 50
                ),
            }
        )
    ldl = interval_ldl(comparison, mp.mpf(0))

    # Central comparison spectrum is reported for scale, while interval LDL
    # is the proof.
    central_angles = np.asarray([float(value) for value in root_angles])
    central_phases = np.empty(count)
    central_phases[0::2] = central_angles
    central_phases[1::2] = -central_angles
    central_beta = np.asarray(
        [
            (
                sizes[i]
                + sum(
                    degrees[i][j]
                    * math.cos(central_phases[i] - central_phases[j])
                    for j in range(count)
                )
            )
            / total
            for i in range(count)
        ]
    )
    central = np.diag(central_beta)
    for row in block_rows:
        i, j = row["classes"]
        central[i, j] = central[j, i] = -float(
            row["weighted_norm_upper_normalized"]
        )
    values = np.linalg.eigvalsh(central)
    return {
        "proof_space": (
            "Direct sum of all classwise zero-sum fiber spaces; dimension "
            f"N-2m = {total-count}."
        ),
        "mode_dimension_enumeration": {
            "fiber_dimensions": [size - 1 for size in sizes],
            "sum": total - count,
            "all_dimensions_covered": True,
        },
        "complete_blocks_vanish_on_zero_sum_spaces": True,
        "residue_block_bounds": block_rows,
        "beta_intervals_normalized": [str(value) for value in beta],
        "comparison_interval_ldl": ldl,
        "central_comparison_eigenvalues_normalized": values.tolist(),
        "central_normalized_lower": float(values[0]),
        "rigorous_all_transverse_modes_positive": True,
    }


def nonlinear_return(
    sizes: list[int],
    degrees: list[list[int]],
    root_angles: list[mp.mpf],
    quotient_gap: mp.mpf,
    total: int,
) -> dict[str, Any]:
    phases = np.asarray(
        [float(value) for value in reflected_phases(root_angles)]
    )
    degree_matrix = np.asarray(degrees, dtype=float) / float(total)
    size_weights = np.asarray(sizes, dtype=float)
    count = len(sizes)

    # Build the symmetric mass-orthonormal quotient independently in doubles.
    quotient = np.zeros((count, count))
    for i in range(count):
        for j in range(count):
            if i == j or not degrees[i][j]:
                continue
            cosine = math.cos(phases[i] - phases[j])
            quotient[i, i] += degree_matrix[i, j] * cosine
            quotient[i, j] = (
                -math.sqrt(degrees[i][j] * degrees[j][i])
                / total
                * cosine
            )
    values, vectors = np.linalg.eigh(0.5 * (quotient + quotient.T))
    positive = np.flatnonzero(values > 1e-10)
    mode = vectors[:, int(positive[0])] / np.sqrt(size_weights)
    mode /= np.max(np.abs(mode))

    time_scale = max(float(quotient_gap), 1e-9)

    def physical_flow(values: np.ndarray) -> np.ndarray:
        delta = values[None, :] - values[:, None]
        return np.sum(degree_matrix * np.sin(delta), axis=1)

    def scaled_flow(_time: float, values: np.ndarray) -> np.ndarray:
        return physical_flow(values) / time_scale

    def scaled_jacobian(_time: float, values: np.ndarray) -> np.ndarray:
        delta = values[None, :] - values[:, None]
        jacobian = degree_matrix * np.cos(delta)
        np.fill_diagonal(jacobian, 0.0)
        np.fill_diagonal(jacobian, -jacobian.sum(axis=1))
        return jacobian / time_scale

    def distance(values: np.ndarray) -> float:
        delta = values - phases
        delta -= np.dot(size_weights, delta) / size_weights.sum()
        return float(
            np.sqrt(np.dot(size_weights, delta * delta) / size_weights.sum())
        )

    scaled_final_time = 6.0
    physical_final_time = scaled_final_time / time_scale
    runs = []
    for amplitude in (1e-7, 3e-7):
        for sign in (-1, 1):
            initial = phases + sign * amplitude * mode
            solution = solve_ivp(
                scaled_flow,
                (0.0, scaled_final_time),
                initial,
                method="Radau",
                jac=scaled_jacobian,
                rtol=2e-8,
                atol=2e-12,
                t_eval=[0.0, scaled_final_time],
                max_step=0.5,
            )
            initial_distance = distance(solution.y[:, 0])
            final_distance = distance(solution.y[:, -1])
            returned = bool(
                solution.success
                and final_distance < max(2e-12, initial_distance * 0.02)
            )
            if not returned:
                raise AssertionError("weak-mode nonlinear trajectory did not return")
            runs.append(
                {
                    "amplitude": amplitude,
                    "sign": sign,
                    "initial_distance_mod_rotation": initial_distance,
                    "final_distance_mod_rotation": final_distance,
                    "contraction_ratio": final_distance / initial_distance,
                    "solver_success": bool(solution.success),
                    "return_pass": returned,
                }
            )
    return {
        "flow": (
            "Exact class-constant restriction of the full unweighted graph "
            "Kuramoto flow, normalized by N."
        ),
        "weak_mode_eigenvalue_normalized_double": float(values[positive[0]]),
        "scaled_final_time": scaled_final_time,
        "physical_normalized_final_time": physical_final_time,
        "time_rescaling": (
            "s = normalized_quotient_gap * t; the exact vector field and "
            "analytic Jacobian are divided by the same factor."
        ),
        "runs": runs,
        "both_signs_and_two_amplitudes_returned": True,
        "scope": "local class-constant nonlinear return",
    }


def replay(spec_path: Path) -> dict[str, Any]:
    mp.mp.dps = POINT_DPS
    mp.iv.dps = INTERVAL_DPS
    gap_floor = mp.mpf(GAP_FLOOR_TEXT)
    spec = parse_spec(spec_path)
    sizes, blocks = parse_graph(spec)
    degrees, edges, edge_count = exact_block_data(sizes, blocks)
    total = sum(sizes)
    class_degrees = [
        sizes[i] - 1 + sum(degrees[i]) for i in range(len(sizes))
    ]
    minimum_degree = min(class_degrees)
    mu = Fraction(minimum_degree, total - 1)
    prior = Fraction(
        int(spec["target_prior_ratio"]["numerator"]),
        int(spec["target_prior_ratio"]["denominator"]),
    )
    handshake = sum(
        sizes[i] * class_degrees[i] for i in range(len(sizes))
    )
    connected = support_connected(sizes, blocks)
    reflection = reflection_proof(sizes, blocks)
    if handshake != 2 * edge_count or not connected or not mu > prior:
        raise AssertionError("exact graph arithmetic acceptance failed")

    root_angles = solve_root(spec, degrees, total)
    residuals = torque_system(root_angles, degrees, total, mp.sin)
    residual = max(abs(value) for value in residuals)
    if residual > mp.mpf("1e-150"):
        raise AssertionError("high-precision torque residual is too large")
    krawczyk, angle_box = krawczyk_certificate(
        root_angles, degrees, total
    )
    quotient = quotient_certificate(
        sizes, degrees, root_angles, angle_box, total
    )
    transverse = transverse_certificate(
        sizes, blocks, degrees, angle_box, root_angles, total
    )
    order_parameter = abs(
        sum(
            mp.mpf(sizes[i])
            * mp.exp(1j * reflected_phases(root_angles)[i])
            for i in range(len(sizes))
        )
        / total
    )
    nonlinear = nonlinear_return(
        sizes, degrees, root_angles, quotient["_gap"], total
    )
    quotient.pop("_gap")

    accepted = bool(
        mu > prior
        and connected
        and handshake == 2 * edge_count
        and reflection["adjacency_preserved_for_every_vertex_pair"]
        and krawczyk["unique_root_in_box"]
        and mp.mpf(quotient["normalized_gap"]) > gap_floor
        and quotient["rigorously_above_1e_minus_6"]
        and quotient["exactly_one_zero"]
        and transverse["rigorous_all_transverse_modes_positive"]
        and nonlinear["both_signs_and_two_amplitudes_returned"]
        and order_parameter < mp.mpf("0.99")
    )
    return {
        "outcome": "WITNESS" if accepted else "REJECTED",
        "independence": (
            "Replay reads only the declarative graph spec and imports no "
            "construction or optimizer code."
        ),
        "spec_path": str(spec_path.resolve()),
        "spec_sha256": sha256(spec_path),
        "graph": {
            "vertex_count": total,
            "class_count": len(sizes),
            "pair_count": len(sizes) // 2,
            "class_sizes": sizes,
            "class_degrees": class_degrees,
            "minimum_degree": minimum_degree,
            "edge_count": edge_count,
            "degree_handshake": handshake,
            "mu_exact": f"{mu.numerator}/{mu.denominator}",
            "mu_decimal": mp.nstr(mp.mpf(mu.numerator) / mu.denominator, 60),
            "prior_exact": f"{prior.numerator}/{prior.denominator}",
            "exact_improvement_numerator": (
                mu.numerator * prior.denominator
                - prior.numerator * mu.denominator
            ),
            "simple": True,
            "unweighted": True,
            "connected": connected,
            "cross_block_partition_complete": True,
        },
        "reflection": reflection,
        "equilibrium": {
            "angles": [mp.nstr(value, 170) for value in root_angles],
            "phases": [
                mp.nstr(value, 170)
                for value in reflected_phases(root_angles)
            ],
            "maximum_normalized_torque_residual": mp.nstr(residual, 40),
            "krawczyk": krawczyk,
            "order_parameter": mp.nstr(order_parameter, 60),
            "nonsynchronous": bool(order_parameter < mp.mpf("0.99")),
        },
        "quotient": quotient,
        "transverse": transverse,
        "nonlinear": nonlinear,
        "acceptance": {
            "exact_mu_beats_prior": mu > prior,
            "simple_connected_unweighted": connected,
            "reflection_automorphism_vertexwise": True,
            "unique_nonsynchronous_root_in_box": True,
            "exactly_one_zero": quotient["exactly_one_zero"],
            "rigorous_normalized_quotient_gap_above_1e_minus_6": True,
            "all_transverse_modes_positive": True,
            "independent_local_nonlinear_return": True,
            "accepted": accepted,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--spec", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    result = replay(args.spec)
    write_json(args.output, result)
    print(
        json.dumps(
            {
                "outcome": result["outcome"],
                "mu_exact": result["graph"]["mu_exact"],
                "vertex_count": result["graph"]["vertex_count"],
                "normalized_gap": result["quotient"]["normalized_gap"],
                "transverse_lower": result["transverse"][
                    "central_normalized_lower"
                ],
                "output": str(args.output.resolve()),
            },
            indent=2,
        )
    )
    if result["outcome"] != "WITNESS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
