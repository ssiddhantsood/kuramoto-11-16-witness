#!/usr/bin/env python3
"""Independent exact/Arb certificate for the six-class Kuramoto witness.

This verifier reads ../graph_spec.json as data. It imports no producer code.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path
from typing import Any, Iterable, Sequence

import flint
from flint import arb, ctx, fmpq


PRECISION_BITS = 512
ROOT_RADIUS = Fraction(1, 10**40)
PRECONDITIONER_DECIMAL_PLACES = 70
QUOTIENT_LOWER_BOUND = Fraction(487, 10)
TRANSVERSE_LOWER_BOUND = Fraction(102_000)
REFLECTION = (1, 0, 3, 2, 5, 4)
REDUCED_TORQUE_CLASSES = (0, 2, 5)
PARAMETER_COEFFICIENTS = (
    (0, 0, 0),   # 0
    (1, 0, 0),   # s
    (0, 1, 0),   # x
    (1, -1, 0),  # s-x
    (1, 0, -1),  # s-y
    (0, 0, 1),   # y
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def fq(value: Fraction | int) -> fmpq:
    value = Fraction(value)
    return fmpq(value.numerator, value.denominator)


def ar(value: Fraction | int) -> arb:
    return arb(fq(Fraction(value)))


def interval(lo: Fraction, hi: Fraction) -> arb:
    require(lo <= hi, "invalid rational interval")
    return ar(lo).union(ar(hi))


def dyadic_value(value: arb) -> Fraction:
    """Convert an exact Arb floating-point value to a Python Fraction."""
    require(value.is_exact(), "Arb endpoint was not exact")
    mantissa, exponent = value.man_exp()
    mantissa = int(mantissa)
    exponent = int(exponent)
    if exponent >= 0:
        return Fraction(mantissa * (2**exponent))
    return Fraction(mantissa, 2 ** (-exponent))


def lower(value: arb) -> Fraction:
    return dyadic_value(value.lower())


def upper(value: arb) -> Fraction:
    return dyadic_value(value.upper())


def midpoint(value: arb) -> Fraction:
    return dyadic_value(value.mid())


def frac_text(value: Fraction | int) -> str:
    value = Fraction(value)
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def interval_json(value: arb) -> dict[str, str]:
    return {"lower": frac_text(lower(value)), "upper": frac_text(upper(value))}


def matrix_interval_json(matrix: Sequence[Sequence[arb]]) -> list[list[dict[str, str]]]:
    return [[interval_json(value) for value in row] for row in matrix]


def fraction_matrix_json(
    matrix: Sequence[Sequence[Fraction]],
) -> list[list[str]]:
    return [[frac_text(value) for value in row] for row in matrix]


def exact_decimal(value: Fraction) -> str:
    """Print a terminating rational decimal without rounding."""
    value = Fraction(value)
    sign = "-" if value < 0 else ""
    numerator = abs(value.numerator)
    denominator = value.denominator
    twos = 0
    fives = 0
    while denominator % 2 == 0:
        denominator //= 2
        twos += 1
    while denominator % 5 == 0:
        denominator //= 5
        fives += 1
    require(denominator == 1, "nonterminating decimal requested")
    places = max(twos, fives)
    numerator *= 2 ** (places - fives)
    numerator *= 5 ** (places - twos)
    if places == 0:
        return f"{sign}{numerator}"
    digits = str(numerator).rjust(places + 1, "0")
    return f"{sign}{digits[:-places]}.{digits[-places:]}"


def nearest_decimal(value: Fraction, places: int) -> Fraction:
    scale = 10**places
    scaled = value * scale
    floor_value = scaled.numerator // scaled.denominator
    remainder = scaled - floor_value
    rounded = floor_value + (1 if remainder >= Fraction(1, 2) else 0)
    return Fraction(rounded, scale)


def sup_abs(value: arb) -> Fraction:
    return max(abs(lower(value)), abs(upper(value)))


def add(*values: arb) -> arb:
    result = arb(0)
    for value in values:
        result += value
    return result


def fraction_matrix_inverse(
    matrix: Sequence[Sequence[Fraction]],
) -> list[list[Fraction]]:
    n = len(matrix)
    require(n > 0 and all(len(row) == n for row in matrix), "matrix is not square")
    augmented = [
        [Fraction(matrix[i][j]) for j in range(n)]
        + [Fraction(int(i == j)) for j in range(n)]
        for i in range(n)
    ]
    for column in range(n):
        pivot_row = next(
            (row for row in range(column, n) if augmented[row][column] != 0),
            None,
        )
        require(pivot_row is not None, "singular rational matrix")
        if pivot_row != column:
            augmented[column], augmented[pivot_row] = (
                augmented[pivot_row],
                augmented[column],
            )
        pivot = augmented[column][column]
        augmented[column] = [value / pivot for value in augmented[column]]
        for row in range(n):
            if row == column:
                continue
            factor = augmented[row][column]
            if factor:
                augmented[row] = [
                    augmented[row][j] - factor * augmented[column][j]
                    for j in range(2 * n)
                ]
    return [row[n:] for row in augmented]


def fraction_matrix_determinant(
    matrix: Sequence[Sequence[Fraction]],
) -> Fraction:
    n = len(matrix)
    work = [[Fraction(value) for value in row] for row in matrix]
    sign = 1
    determinant = Fraction(1)
    for column in range(n):
        pivot_row = next(
            (row for row in range(column, n) if work[row][column] != 0),
            None,
        )
        if pivot_row is None:
            return Fraction(0)
        if pivot_row != column:
            work[column], work[pivot_row] = work[pivot_row], work[column]
            sign *= -1
        pivot = work[column][column]
        determinant *= pivot
        for row in range(column + 1, n):
            factor = work[row][column] / pivot
            for j in range(column + 1, n):
                work[row][j] -= factor * work[column][j]
    return sign * determinant


def interval_ldlt(
    matrix: Sequence[Sequence[arb]],
) -> tuple[list[list[arb]], list[arb]]:
    """Interval LDL^T; positive pivot balls prove every point matrix is PD."""
    n = len(matrix)
    require(n > 0 and all(len(row) == n for row in matrix), "matrix is not square")
    l_matrix = [[arb(0) for _ in range(n)] for _ in range(n)]
    pivots: list[arb] = []
    for k in range(n):
        pivot = matrix[k][k]
        for j in range(k):
            pivot -= l_matrix[k][j] * l_matrix[k][j] * pivots[j]
        require(pivot > arb(0), f"interval LDL pivot {k} is not strictly positive")
        pivots.append(pivot)
        l_matrix[k][k] = arb(1)
        for i in range(k + 1, n):
            numerator = matrix[i][k]
            for j in range(k):
                numerator -= l_matrix[i][j] * l_matrix[k][j] * pivots[j]
            l_matrix[i][k] = numerator / pivot
    return l_matrix, pivots


def phases(parameters: Sequence[arb]) -> list[arb]:
    s, x, y = parameters
    return [arb(0), s, x, s - x, s - y, y]


def reduced_torques(
    parameters: Sequence[arb],
    directed_degrees: Sequence[Sequence[int]],
) -> list[arb]:
    phase = phases(parameters)
    values: list[arb] = []
    for i in REDUCED_TORQUE_CLASSES:
        torque = arb(0)
        for j in range(6):
            if directed_degrees[i][j]:
                torque += directed_degrees[i][j] * (phase[j] - phase[i]).sin()
        values.append(torque)
    return values


def reduced_jacobian(
    parameters: Sequence[arb],
    directed_degrees: Sequence[Sequence[int]],
) -> list[list[arb]]:
    phase = phases(parameters)
    jacobian = [[arb(0) for _ in range(3)] for _ in range(3)]
    for row, i in enumerate(REDUCED_TORQUE_CLASSES):
        for j in range(6):
            degree = directed_degrees[i][j]
            if not degree:
                continue
            cosine = (phase[j] - phase[i]).cos()
            for column in range(3):
                coefficient = (
                    PARAMETER_COEFFICIENTS[j][column]
                    - PARAMETER_COEFFICIENTS[i][column]
                )
                if coefficient:
                    jacobian[row][column] += degree * coefficient * cosine
    return jacobian


def build_quotient(
    parameter_box: Sequence[arb],
    directed_degrees: Sequence[Sequence[int]],
) -> list[list[arb]]:
    phase = phases(parameter_box)
    quotient = [[arb(0) for _ in range(6)] for _ in range(6)]
    for i in range(6):
        for j in range(i + 1, 6):
            degree_ij = directed_degrees[i][j]
            degree_ji = directed_degrees[j][i]
            if not degree_ij:
                continue
            cosine = (phase[i] - phase[j]).cos()
            quotient[i][i] += degree_ij * cosine
            quotient[j][j] += degree_ji * cosine
            singular_bound = ar(Fraction(degree_ij * degree_ji)).sqrt()
            off_diagonal = -singular_bound * cosine
            quotient[i][j] = off_diagonal
            quotient[j][i] = off_diagonal
    return quotient


def build_projection(class_sizes: Sequence[int]) -> list[list[arb]]:
    total = sum(class_sizes)
    projection = [[arb(0) for _ in range(6)] for _ in range(6)]
    for i in range(6):
        for j in range(6):
            radial = ar(Fraction(class_sizes[i] * class_sizes[j])).sqrt() / total
            projection[i][j] = arb(int(i == j)) - radial
    return projection


def build_transverse_comparison(
    parameter_box: Sequence[arb],
    class_sizes: Sequence[int],
    directed_degrees: Sequence[Sequence[int]],
    block_types: dict[tuple[int, int], str],
) -> list[list[arb]]:
    phase = phases(parameter_box)
    comparison = [[arb(0) for _ in range(6)] for _ in range(6)]
    for i in range(6):
        comparison[i][i] = arb(class_sizes[i])
    for i in range(6):
        for j in range(i + 1, 6):
            degree_ij = directed_degrees[i][j]
            degree_ji = directed_degrees[j][i]
            if not degree_ij:
                continue
            cosine = (phase[i] - phase[j]).cos()
            comparison[i][i] += degree_ij * cosine
            comparison[j][j] += degree_ji * cosine
            if block_types[i, j] == "residue_orbits":
                singular_bound = ar(Fraction(degree_ij * degree_ji)).sqrt()
                off_diagonal = -abs(cosine) * singular_bound
                comparison[i][j] = off_diagonal
                comparison[j][i] = off_diagonal
    return comparison


def graph_data(spec: dict[str, Any]) -> dict[str, Any]:
    require(spec.get("schema") == "six-class-circulant-witness-v1", "wrong schema")
    classes = spec["classes"]
    require([entry["id"] for entry in classes] == list(range(6)), "class IDs differ")
    class_sizes = [int(entry["size"]) for entry in classes]
    require(all(size > 1 for size in class_sizes), "clique fiber is too small")
    total_vertices = sum(class_sizes)
    require(total_vertices == int(spec["vertex_count"]), "vertex count mismatch")
    require(
        spec["within_class"] == {"type": "clique", "self_loops": False},
        "within-class rule differs",
    )

    directed_degrees = [[0 for _ in range(6)] for _ in range(6)]
    cross_edge_counts = [[0 for _ in range(6)] for _ in range(6)]
    block_types: dict[tuple[int, int], str] = {}
    block_certificates: list[dict[str, Any]] = []
    seen_pairs: set[tuple[int, int]] = set()
    for block in spec["cross_blocks"]:
        i, j = map(int, block["classes"])
        require(0 <= i < j < 6, "cross-block pair is not ordered")
        require((i, j) not in seen_pairs, "duplicate cross block")
        seen_pairs.add((i, j))
        block_type = block["type"]
        block_types[i, j] = block_type
        if block_type == "absent":
            edge_count = 0
            degree_ij = 0
            degree_ji = 0
        elif block_type == "complete":
            edge_count = class_sizes[i] * class_sizes[j]
            degree_ij = class_sizes[j]
            degree_ji = class_sizes[i]
        elif block_type == "residue_orbits":
            modulus = int(block["modulus"])
            shift_start = int(block["shift_start"])
            shift_count = int(block["shift_count"])
            require(modulus > 0, "nonpositive residue modulus")
            require(class_sizes[i] % modulus == 0, "modulus does not divide class i")
            require(class_sizes[j] % modulus == 0, "modulus does not divide class j")
            require(0 <= shift_start < modulus, "invalid shift start")
            require(0 < shift_count <= modulus, "invalid shift count")
            degree_ij = shift_count * (class_sizes[j] // modulus)
            degree_ji = shift_count * (class_sizes[i] // modulus)
            edge_count = class_sizes[i] * degree_ij
            require(
                edge_count == class_sizes[j] * degree_ji,
                "biregularity identity failed",
            )
        else:
            raise AssertionError(f"unknown block type {block_type}")
        directed_degrees[i][j] = degree_ij
        directed_degrees[j][i] = degree_ji
        cross_edge_counts[i][j] = edge_count
        cross_edge_counts[j][i] = edge_count
        block_certificate: dict[str, Any] = {
            "classes": [i, j],
            "type": block_type,
            "directed_degrees": {f"{i}_to_{j}": degree_ij, f"{j}_to_{i}": degree_ji},
            "edge_count": edge_count,
            "biregularity_identity": f"{class_sizes[i]}*{degree_ij}"
            f"={class_sizes[j]}*{degree_ji}={edge_count}",
        }
        if block_type == "residue_orbits":
            block_certificate["residue_data"] = {
                "modulus": int(block["modulus"]),
                "shift_start": int(block["shift_start"]),
                "shift_count": int(block["shift_count"]),
                "multiplicity_in_first_class": class_sizes[i] // int(block["modulus"]),
                "multiplicity_in_second_class": class_sizes[j] // int(block["modulus"]),
            }
        block_certificates.append(block_certificate)
    require(
        seen_pairs == {(i, j) for i in range(6) for j in range(i + 1, 6)},
        "not all cross blocks were specified exactly once",
    )

    support_adjacency = [set() for _ in range(6)]
    for i in range(6):
        for j in range(i + 1, 6):
            if cross_edge_counts[i][j] > 0:
                support_adjacency[i].add(j)
                support_adjacency[j].add(i)
    reached = {0}
    queue = [0]
    spanning_tree: list[list[int]] = []
    while queue:
        i = queue.pop(0)
        for j in sorted(support_adjacency[i]):
            if j not in reached:
                reached.add(j)
                queue.append(j)
                spanning_tree.append([i, j])
    require(reached == set(range(6)), "class support graph is disconnected")

    total_degrees = [
        class_sizes[i] - 1 + sum(directed_degrees[i]) for i in range(6)
    ]
    minimum_degree = min(total_degrees)
    clique_edges = [size * (size - 1) // 2 for size in class_sizes]
    total_edges = sum(clique_edges) + sum(
        cross_edge_counts[i][j] for i in range(6) for j in range(i + 1, 6)
    )
    threshold_excess = 16 * minimum_degree - 11 * (total_vertices - 1)
    require(total_vertices == 460_800, "unexpected N")
    require(minimum_degree == 316_802, "unexpected minimum degree")
    require(threshold_excess == 43, "threshold excess is not 43")
    require(total_edges == 74_003_856_975, "unexpected edge count")
    for i in range(6):
        for j in range(6):
            require(
                directed_degrees[REFLECTION[i]][REFLECTION[j]]
                == directed_degrees[i][j],
                "directed degrees violate reflection symmetry",
            )

    return {
        "class_sizes": class_sizes,
        "total_vertices": total_vertices,
        "directed_degrees": directed_degrees,
        "cross_edge_counts": cross_edge_counts,
        "block_types": block_types,
        "block_certificates": block_certificates,
        "support_edges": [
            [i, j]
            for i in range(6)
            for j in range(i + 1, 6)
            if cross_edge_counts[i][j] > 0
        ],
        "support_spanning_tree": spanning_tree,
        "total_degrees": total_degrees,
        "minimum_degree": minimum_degree,
        "clique_edge_counts": clique_edges,
        "total_edges": total_edges,
        "threshold_excess": threshold_excess,
    }


def certify(spec_path: Path) -> tuple[dict[str, Any], str]:
    ctx.prec = PRECISION_BITS
    raw_spec = spec_path.read_bytes()
    spec = json.loads(raw_spec)
    graph = graph_data(spec)
    class_sizes: list[int] = graph["class_sizes"]
    directed_degrees: list[list[int]] = graph["directed_degrees"]
    block_types: dict[tuple[int, int], str] = graph["block_types"]

    phase_strings = spec["equilibrium"]["phase_strings"]
    require(len(phase_strings) == 6, "expected six stored phases")
    stored_phases = [Fraction(value) for value in phase_strings]
    require(stored_phases[0] == 0, "stored gauge is not phi_0=0")

    # The rational unwrapped x center is s - stored_phi_3.  Its wrapped
    # representative x+2*pi is the stored positive phase for class 2.
    center = (
        stored_phases[1],
        stored_phases[1] - stored_phases[3],
        stored_phases[5],
    )
    box_bounds = [
        (value - ROOT_RADIUS, value + ROOT_RADIUS) for value in center
    ]
    center_arb = [ar(value) for value in center]
    box_arb = [interval(lo, hi) for lo, hi in box_bounds]

    center_phase = phases(center_arb)
    wrapped_center_phase = list(center_phase)
    wrapped_center_phase[2] += 2 * arb.pi()
    stored_phase_errors = [
        wrapped_center_phase[i] - ar(stored_phases[i]) for i in range(6)
    ]
    for error in stored_phase_errors:
        require(
            sup_abs(error) < Fraction(1, 10**88),
            "reflection center is not within 1e-88 of stored phases",
        )

    f_center = reduced_torques(center_arb, directed_degrees)
    jacobian_box = reduced_jacobian(box_arb, directed_degrees)
    jacobian_center = reduced_jacobian(center_arb, directed_degrees)

    jacobian_midpoint = [
        [midpoint(jacobian_center[i][j]) for j in range(3)] for i in range(3)
    ]
    preconditioner_unrounded = fraction_matrix_inverse(jacobian_midpoint)
    preconditioner = [
        [
            nearest_decimal(value, PRECONDITIONER_DECIMAL_PLACES)
            for value in row
        ]
        for row in preconditioner_unrounded
    ]
    preconditioner_determinant = fraction_matrix_determinant(preconditioner)
    require(preconditioner_determinant != 0, "preconditioner is singular")
    preconditioner_arb = [[ar(value) for value in row] for row in preconditioner]

    defect = [[arb(int(i == j)) for j in range(3)] for i in range(3)]
    for i in range(3):
        for j in range(3):
            product = arb(0)
            for k in range(3):
                product += preconditioner_arb[i][k] * jacobian_box[k][j]
            defect[i][j] -= product

    contraction_rows = [sum(sup_abs(value) for value in row) for row in defect]
    contraction_bound = max(contraction_rows)
    require(
        contraction_bound < Fraction(1, 10**34),
        f"weak contraction bound: {frac_text(contraction_bound)}",
    )

    root_image: list[arb] = []
    inclusion_margins: list[Fraction] = []
    delta = interval(-ROOT_RADIUS, ROOT_RADIUS)
    for i in range(3):
        offset = arb(0)
        for j in range(3):
            offset -= preconditioner_arb[i][j] * f_center[j]
            offset += defect[i][j] * delta
        image = center_arb[i] + offset
        root_image.append(image)
        lo, hi = box_bounds[i]
        margin = min(lower(image) - lo, hi - upper(image))
        require(margin > 0, f"Krawczyk image {i} is not in the open box")
        inclusion_margins.append(margin)
    require(
        min(inclusion_margins) > Fraction(9, 10**41),
        "root inclusion margin unexpectedly small",
    )

    root_phase_box = phases(box_arb)
    wrapped_root_phase_box = list(root_phase_box)
    wrapped_root_phase_box[2] += 2 * arb.pi()
    real_sum = arb(0)
    imaginary_sum = arb(0)
    for size, phase in zip(class_sizes, root_phase_box):
        real_sum += size * phase.cos()
        imaginary_sum += size * phase.sin()
    order_parameter_squared = (
        real_sum * real_sum + imaginary_sum * imaginary_sum
    ) / (graph["total_vertices"] ** 2)
    require(
        upper(order_parameter_squared) < Fraction(1, 900),
        "order parameter is not certified below 1/30",
    )

    quotient = build_quotient(box_arb, directed_degrees)
    projection = build_projection(class_sizes)
    shifted_quotient = [
        [
            quotient[i][j] - ar(QUOTIENT_LOWER_BOUND) * projection[i][j]
            for j in range(6)
        ]
        for i in range(6)
    ]
    quotient_ground_index = 0
    quotient_grounded = [
        [
            shifted_quotient[i][j]
            for j in range(6)
            if j != quotient_ground_index
        ]
        for i in range(6)
        if i != quotient_ground_index
    ]
    _, quotient_pivots = interval_ldlt(quotient_grounded)
    quotient_pivot_simple_lowers = (37_000, 2_500, 1_000, 13, 7)
    for pivot, simple_lower in zip(
        quotient_pivots, quotient_pivot_simple_lowers, strict=True
    ):
        require(pivot > arb(simple_lower), "weak quotient LDL pivot")

    comparison = build_transverse_comparison(
        box_arb,
        class_sizes,
        directed_degrees,
        block_types,
    )
    shifted_comparison = [
        [
            comparison[i][j]
            - (ar(TRANSVERSE_LOWER_BOUND) if i == j else arb(0))
            for j in range(6)
        ]
        for i in range(6)
    ]
    _, transverse_pivots = interval_ldlt(shifted_comparison)
    transverse_pivot_simple_lowers = (16_000, 9_700, 4_600, 1_400, 25_000, 25_000)
    for pivot, simple_lower in zip(
        transverse_pivots, transverse_pivot_simple_lowers, strict=True
    ):
        require(pivot > arb(simple_lower), "weak transverse LDL pivot")

    partial_norms = []
    for (i, j), block_type in sorted(block_types.items()):
        if block_type != "residue_orbits":
            continue
        degree_ij = directed_degrees[i][j]
        degree_ji = directed_degrees[j][i]
        partial_norms.append(
            {
                "classes": [i, j],
                "degree_ij": degree_ij,
                "degree_ji": degree_ji,
                "adjacency_norm_squared_upper": degree_ij * degree_ji,
                "adjacency_norm_upper": {
                    "radicand": degree_ij * degree_ji,
                    "expression": f"sqrt({degree_ij * degree_ji})",
                },
            }
        )

    transverse_dimension = sum(size - 1 for size in class_sizes)
    require(transverse_dimension == 460_794, "wrong transverse dimension")
    require(6 + transverse_dimension == graph["total_vertices"], "dimension mismatch")

    certificate: dict[str, Any] = {
        "schema": "formal-interval-kuramoto-witness-certificate-v1",
        "status": "CONFIRMED_RIGOROUS",
        "independence": {
            "construction_input_only": "../graph_spec.json",
            "imports_external_producer_code": False,
            "comparison_certificate_read_by_verifier": False,
        },
        "rigor": {
            "library": "python-flint",
            "library_version": getattr(flint, "__version__", "unknown"),
            "backend": "Arb midpoint-radius ball arithmetic",
            "precision_bits": PRECISION_BITS,
            "rounding": "Arb operations and lower()/upper() are outward-rounded",
            "interval_endpoints_encoding": "exact reduced rational dyadic bounds",
        },
        "input": {
            "spec_sha256": hashlib.sha256(raw_spec).hexdigest(),
            "vertex_count": graph["total_vertices"],
            "class_sizes": class_sizes,
            "phase_parameterization_unwrapped": "[0,s,x,s-x,s-y,y]",
            "stored_branch_for_class_2": "x+2*pi",
            "parameter_center_rationals": {
                name: frac_text(value) for name, value in zip(("s", "x", "y"), center)
            },
            "parameter_radius_rational": frac_text(ROOT_RADIUS),
            "preconditioner_rational": fraction_matrix_json(preconditioner),
            "quotient_target_rational": frac_text(QUOTIENT_LOWER_BOUND),
            "transverse_target_rational": frac_text(TRANSVERSE_LOWER_BOUND),
        },
        "graph": {
            "simple": True,
            "undirected": True,
            "biregular_cross_blocks": True,
            "connected": True,
            "class_sizes": class_sizes,
            "directed_cross_degrees": directed_degrees,
            "block_certificates": graph["block_certificates"],
            "support_edges": graph["support_edges"],
            "support_spanning_tree": graph["support_spanning_tree"],
            "connectivity_argument": (
                "Each fiber is a clique. The listed support spanning tree has a "
                "nonempty cross block on every tree edge, so clique paths to cross-edge "
                "endpoints lift every support path to a graph path."
            ),
            "within_class_degrees": [size - 1 for size in class_sizes],
            "total_degrees": graph["total_degrees"],
            "minimum_degree": graph["minimum_degree"],
            "edge_count": graph["total_edges"],
            "minimum_degree_ratio": {
                "numerator": graph["minimum_degree"],
                "denominator": graph["total_vertices"] - 1,
                "reduced": frac_text(
                    Fraction(graph["minimum_degree"], graph["total_vertices"] - 1)
                ),
            },
            "threshold_comparison": {
                "target": "11/16",
                "cross_product_excess": graph["threshold_excess"],
                "identity": (
                    f"16*{graph['minimum_degree']}"
                    f"-11*{graph['total_vertices'] - 1}"
                    f"={graph['threshold_excess']}"
                ),
            },
        },
        "root": {
            "reduced_torque_classes": list(REDUCED_TORQUE_CLASSES),
            "parameter_box": {
                name: {
                    "center": frac_text(center_value),
                    "radius": frac_text(ROOT_RADIUS),
                    "lower": frac_text(bounds[0]),
                    "upper": frac_text(bounds[1]),
                }
                for name, center_value, bounds in zip(
                    ("s", "x", "y"), center, box_bounds
                )
            },
            "stored_phase_center_error_intervals": [
                interval_json(value) for value in stored_phase_errors
            ],
            "center_reduced_torque_intervals": [
                interval_json(value) for value in f_center
            ],
            "jacobian_box_intervals": matrix_interval_json(jacobian_box),
            "preconditioner_determinant": frac_text(preconditioner_determinant),
            "krawczyk_defect_intervals": matrix_interval_json(defect),
            "contraction_inf_norm_upper": frac_text(contraction_bound),
            "contraction_simple_bound": "1/10^34",
            "krawczyk_image_intervals": [
                interval_json(value) for value in root_image
            ],
            "strict_inclusion_margins": [
                frac_text(value) for value in inclusion_margins
            ],
            "conclusion": (
                "G(z)=z-CF(z) maps the rational box strictly into itself and "
                "has derivative infinity norm <1e-34 there. Since det(C)!=0, "
                "Banach's theorem gives exactly one F-root in the box."
            ),
            "full_six_torque_residual_at_root": {
                "bounds": [["0", "0"] for _ in range(6)],
                "reason": (
                    "The reduced torques tau_0,tau_2,tau_5 vanish. Exact reflection "
                    "D[R(i),R(j)]=D[i,j] and phi_R(i)=s-phi_i gives "
                    "tau_R(i)=-tau_i."
                ),
            },
            "wrapped_six_phase_box_intervals": [
                interval_json(value) for value in wrapped_root_phase_box
            ],
        },
        "nonsynchrony": {
            "order_parameter_squared_interval": interval_json(order_parameter_squared),
            "certified_order_parameter_upper": "1/30",
            "certified_distance_from_one_lower": "29/30",
        },
        "quotient_hessian": {
            "basis": "e_i=1_{V_i}/sqrt(n_i)",
            "matrix_interval": matrix_interval_json(quotient),
            "rotation_vector": "[sqrt(n_0),...,sqrt(n_5)]",
            "exact_kernel_identity": "Q*r=0 term-by-term",
            "projection": "P=I-r*r^T/N",
            "certified_nonzero_eigenvalue_lower": frac_text(
                QUOTIENT_LOWER_BOUND
            ),
            "grounded_index": quotient_ground_index,
            "grounded_shifted_matrix_interval": matrix_interval_json(
                quotient_grounded
            ),
            "interval_ldlt_pivots": [
                interval_json(value) for value in quotient_pivots
            ],
            "simple_pivot_lower_bounds": list(quotient_pivot_simple_lowers),
            "argument": (
                "The grounded principal matrix of Q-(487/10)P is interval-PD. "
                "That matrix annihilates r; the grounded-minor lemma makes it PSD "
                "with kernel span(r), so Q is >487/10 on r^perp."
            ),
        },
        "transverse_hessian": {
            "dimension": transverse_dimension,
            "partial_block_norm_bounds": partial_norms,
            "comparison_matrix_interval": matrix_interval_json(comparison),
            "certified_eigenvalue_lower": frac_text(TRANSVERSE_LOWER_BOUND),
            "shifted_comparison_matrix_interval": matrix_interval_json(
                shifted_comparison
            ),
            "interval_ldlt_pivots": [
                interval_json(value) for value in transverse_pivots
            ],
            "simple_pivot_lower_bounds": list(transverse_pivot_simple_lowers),
            "argument": (
                "On each class-zero-sum space a clique contributes n_i I; complete "
                "cross adjacencies vanish. Every partial biregular adjacency satisfies "
                "||A_ij||_2<=sqrt(d_ij*d_ji) by the 1/infinity norm inequality. "
                "Cauchy-Schwarz therefore bounds H below by the comparison matrix, "
                "and interval LDL^T proves C-102000 I positive definite."
            ),
        },
        "dimension_split": {
            "quotient": 6,
            "transverse": transverse_dimension,
            "total": graph["total_vertices"],
            "identity": f"6+{transverse_dimension}={graph['total_vertices']}",
        },
        "theorem": {
            "hessian": (
                "positive semidefinite with kernel exactly span(1_N)"
            ),
            "potential": (
                "strict local minimum modulo global phase rotation"
            ),
            "dynamics": (
                "locally asymptotically stable modulo global phase rotation"
            ),
            "simulation_required": False,
        },
    }

    root_lines = []
    for name, value in zip(("s", "x", "y"), center):
        root_lines.append(f"  {name} = {exact_decimal(value)} +/- 1e-40")
    report = "\n".join(
        [
            "CONFIRMED_RIGOROUS",
            f"spec_sha256={certificate['input']['spec_sha256']}",
            f"N={graph['total_vertices']}",
            f"edge_count={graph['total_edges']}",
            f"minimum_degree={graph['minimum_degree']}",
            (
                "mu="
                f"{graph['minimum_degree']}/{graph['total_vertices'] - 1}"
                f" > 11/16 (cross-product excess {graph['threshold_excess']})"
            ),
            "unique_root_box:",
            *root_lines,
            "contraction_inf_norm<1e-34",
            "order_parameter<1/30",
            f"quotient_nonzero_gap>{frac_text(QUOTIENT_LOWER_BOUND)}",
            f"transverse_gap>{frac_text(TRANSVERSE_LOWER_BOUND)}",
            f"dimension_split=6+{transverse_dimension}={graph['total_vertices']}",
            "Hessian PSD with kernel exactly span(1); strict local minimum and",
            "local asymptotic stability modulo rotation.",
        ]
    ) + "\n"
    return certificate, report


def main() -> None:
    directory = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--spec",
        type=Path,
        default=directory.parent / "graph_spec.json",
    )
    args = parser.parse_args()
    certificate, report = certify(args.spec.resolve())
    (directory / "certificate.json").write_text(
        json.dumps(certificate, indent=2) + "\n"
    )
    (directory / "verifier_output.txt").write_text(report)
    print(report, end="")


if __name__ == "__main__":
    main()
