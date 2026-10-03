#!/usr/bin/env python3
"""Independent Arb certificate for the N=80,002 Kuramoto witness.

This verifier reads only the adjacent graph_spec.json.  It does not import the
mpmath verifier, search code, saved reports, or audit outputs.  It reuses the
exact-rational and Arb helper routines from ../../formal/verify.py, whose
original target is the distinct N=460,800 record.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
from fractions import Fraction
from pathlib import Path
from typing import Any, Sequence

import flint
from flint import arb, ctx


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PRECISION_BITS = 512
ROOT_RADIUS = Fraction(1, 10**100)
PRECONDITIONER_DECIMAL_PLACES = 90
QUOTIENT_GAP_LOWER = Fraction(2, 5)
TRANSVERSE_GAP_LOWER = 17_000
PHASE_COEFFICIENTS = (
    (1, 0, 0),
    (-1, 0, 0),
    (0, 1, 0),
    (0, -1, 0),
    (0, 0, 1),
    (0, 0, -1),
)
REDUCED_CLASSES = (0, 2, 4)
REFLECTION = (1, 0, 3, 2, 5, 4)


def load_arb_helpers():
    path = ROOT / "formal" / "verify.py"
    module_spec = importlib.util.spec_from_file_location("arb_helpers", path)
    if module_spec is None or module_spec.loader is None:
        raise RuntimeError("could not load formal Arb helpers")
    module = importlib.util.module_from_spec(module_spec)
    module_spec.loader.exec_module(module)
    return module


helpers = load_arb_helpers()
require = helpers.require
ar = helpers.ar
interval = helpers.interval
lower = helpers.lower
upper = helpers.upper
midpoint = helpers.midpoint
sup_abs = helpers.sup_abs
frac_text = helpers.frac_text
interval_json = helpers.interval_json
matrix_interval_json = helpers.matrix_interval_json
fraction_matrix_json = helpers.fraction_matrix_json
fraction_matrix_inverse = helpers.fraction_matrix_inverse
fraction_matrix_determinant = helpers.fraction_matrix_determinant
nearest_decimal = helpers.nearest_decimal
interval_ldlt = helpers.interval_ldlt


def phases(parameters: Sequence[arb]) -> list[arb]:
    a, b, c = parameters
    return [a, -a, b, -b, c, -c]


def reduced_torques(
    parameters: Sequence[arb], degrees: Sequence[Sequence[int]]
) -> list[arb]:
    phase = phases(parameters)
    result = []
    for i in REDUCED_CLASSES:
        value = arb(0)
        for j in range(6):
            value += degrees[i][j] * (phase[j] - phase[i]).sin()
        result.append(value)
    return result


def reduced_jacobian(
    parameters: Sequence[arb], degrees: Sequence[Sequence[int]]
) -> list[list[arb]]:
    phase = phases(parameters)
    result = [[arb(0) for _ in range(3)] for _ in range(3)]
    for row, i in enumerate(REDUCED_CLASSES):
        for j in range(6):
            if not degrees[i][j]:
                continue
            cosine = (phase[j] - phase[i]).cos()
            for column in range(3):
                coefficient = (
                    PHASE_COEFFICIENTS[j][column]
                    - PHASE_COEFFICIENTS[i][column]
                )
                result[row][column] += degrees[i][j] * coefficient * cosine
    return result


def graph_data(spec: dict[str, Any]) -> dict[str, Any]:
    require(
        spec.get("schema") == "reflection-paired-circulant-witness-v1",
        "wrong graph schema",
    )
    classes = sorted(spec["classes"], key=lambda item: int(item["id"]))
    require([int(item["id"]) for item in classes] == list(range(6)), "bad IDs")
    sizes = [int(item["size"]) for item in classes]
    require(sum(sizes) == int(spec["vertex_count"]), "bad vertex count")
    require(all(size > 1 for size in sizes), "empty or singleton clique fiber")
    require(
        spec["within_class"] == {"type": "clique", "self_loops": False},
        "within-class rule is not loopless clique",
    )

    degrees = [[0 for _ in range(6)] for _ in range(6)]
    block_types: dict[tuple[int, int], str] = {}
    seen: set[tuple[int, int]] = set()
    cross_edges = 0
    for block in spec["cross_blocks"]:
        i, j = map(int, block["classes"])
        require(0 <= i < j < 6 and (i, j) not in seen, "bad block pair")
        seen.add((i, j))
        kind = str(block["type"])
        block_types[i, j] = kind
        if kind == "absent":
            degree_ij = degree_ji = 0
        elif kind == "complete":
            degree_ij, degree_ji = sizes[j], sizes[i]
        elif kind == "residue_orbits":
            modulus = int(block["modulus"])
            count = int(block["shift_count"])
            start = int(block["shift_start"])
            require(modulus == math.gcd(sizes[i], sizes[j]), "wrong modulus")
            require(0 <= start < modulus and 0 < count <= modulus, "bad residues")
            degree_ij = count * (sizes[j] // modulus)
            degree_ji = count * (sizes[i] // modulus)
        else:
            raise AssertionError(f"unknown block type {kind}")
        require(sizes[i] * degree_ij == sizes[j] * degree_ji, "not biregular")
        degrees[i][j], degrees[j][i] = degree_ij, degree_ji
        cross_edges += sizes[i] * degree_ij

    require(
        seen == {(i, j) for i in range(6) for j in range(i + 1, 6)},
        "not all cross blocks specified",
    )
    for i in range(6):
        for j in range(6):
            require(
                degrees[REFLECTION[i]][REFLECTION[j]] == degrees[i][j],
                "degree quotient is not reflection symmetric",
            )

    support = [set() for _ in range(6)]
    for i in range(6):
        for j in range(i + 1, 6):
            if degrees[i][j] > 0:
                support[i].add(j)
                support[j].add(i)
    reached, stack = {0}, [0]
    while stack:
        i = stack.pop()
        for j in support[i]:
            if j not in reached:
                reached.add(j)
                stack.append(j)
    require(reached == set(range(6)), "class support graph is disconnected")

    total_degrees = [sizes[i] - 1 + sum(degrees[i]) for i in range(6)]
    minimum_degree = min(total_degrees)
    total_edges = sum(size * (size - 1) // 2 for size in sizes) + cross_edges
    require(sum(sizes) == 80_002, "unexpected N")
    require(total_degrees == [55_021, 55_021, 55_018, 55_018, 56_684, 56_684],
            "unexpected class degrees")
    require(minimum_degree == 55_018, "unexpected minimum degree")
    require(total_edges == 2_232_211_521, "unexpected edge count")
    require(16 * minimum_degree - 11 * (sum(sizes) - 1) == 277,
            "11/16 comparison failed")
    return {
        "sizes": sizes,
        "degrees": degrees,
        "block_types": block_types,
        "total_degrees": total_degrees,
        "minimum_degree": minimum_degree,
        "total_edges": total_edges,
    }


def quotient_hessian(
    parameter_box: Sequence[arb], sizes: Sequence[int], degrees: Sequence[Sequence[int]]
) -> list[list[arb]]:
    phase = phases(parameter_box)
    matrix = [[arb(0) for _ in range(6)] for _ in range(6)]
    for i in range(6):
        for j in range(i + 1, 6):
            if not degrees[i][j]:
                continue
            cosine = (phase[i] - phase[j]).cos()
            matrix[i][i] += degrees[i][j] * cosine
            matrix[j][j] += degrees[j][i] * cosine
            off = -ar(Fraction(degrees[i][j] * degrees[j][i])).sqrt() * cosine
            matrix[i][j] = matrix[j][i] = off
    return matrix


def transverse_comparison(
    parameter_box: Sequence[arb],
    sizes: Sequence[int],
    degrees: Sequence[Sequence[int]],
    block_types: dict[tuple[int, int], str],
) -> list[list[arb]]:
    phase = phases(parameter_box)
    matrix = [[arb(0) for _ in range(6)] for _ in range(6)]
    for i in range(6):
        matrix[i][i] = arb(sizes[i])
    for i in range(6):
        for j in range(i + 1, 6):
            if not degrees[i][j]:
                continue
            cosine = (phase[i] - phase[j]).cos()
            matrix[i][i] += degrees[i][j] * cosine
            matrix[j][j] += degrees[j][i] * cosine
            if block_types[i, j] == "residue_orbits":
                off = -abs(cosine) * ar(
                    Fraction(degrees[i][j] * degrees[j][i])
                ).sqrt()
                matrix[i][j] = matrix[j][i] = off
    return matrix


def certify(spec_path: Path) -> dict[str, Any]:
    ctx.prec = PRECISION_BITS
    raw = spec_path.read_bytes()
    spec = json.loads(raw)
    graph = graph_data(spec)
    sizes = graph["sizes"]
    degrees = graph["degrees"]

    center = tuple(Fraction(value) for value in spec["equilibrium"]["angle_strings"])
    require(len(center) == 3, "expected three angle centers")
    center_arb = [ar(value) for value in center]
    box_bounds = [(value - ROOT_RADIUS, value + ROOT_RADIUS) for value in center]
    box = [interval(lo, hi) for lo, hi in box_bounds]

    f_center = reduced_torques(center_arb, degrees)
    jac_center = reduced_jacobian(center_arb, degrees)
    jac_box = reduced_jacobian(box, degrees)
    jac_mid = [[midpoint(jac_center[i][j]) for j in range(3)] for i in range(3)]
    inverse = fraction_matrix_inverse(jac_mid)
    preconditioner = [
        [nearest_decimal(value, PRECONDITIONER_DECIMAL_PLACES) for value in row]
        for row in inverse
    ]
    require(fraction_matrix_determinant(preconditioner) != 0, "singular C")
    c_arb = [[ar(value) for value in row] for row in preconditioner]

    defect = [[arb(int(i == j)) for j in range(3)] for i in range(3)]
    for i in range(3):
        for j in range(3):
            for k in range(3):
                defect[i][j] -= c_arb[i][k] * jac_box[k][j]
    contraction = max(sum(sup_abs(value) for value in row) for row in defect)
    require(
        contraction < Fraction(1, 10**80),
        f"weak Krawczyk contraction: {frac_text(contraction)}",
    )

    delta = interval(-ROOT_RADIUS, ROOT_RADIUS)
    images: list[arb] = []
    margins: list[Fraction] = []
    for i in range(3):
        offset = arb(0)
        for j in range(3):
            offset -= c_arb[i][j] * f_center[j]
            offset += defect[i][j] * delta
        image = center_arb[i] + offset
        images.append(image)
        lo, hi = box_bounds[i]
        margin = min(lower(image) - lo, hi - upper(image))
        require(margin > 0, "Krawczyk image is not strictly interior")
        margins.append(margin)

    quotient = quotient_hessian(box, sizes, degrees)
    total = sum(sizes)
    shifted_quotient = [[arb(quotient[i][j]) for j in range(6)] for i in range(6)]
    for i in range(6):
        for j in range(6):
            rotation_projector = ar(Fraction(sizes[i] * sizes[j])).sqrt() / total
            shifted_quotient[i][j] += rotation_projector
        shifted_quotient[i][i] -= ar(QUOTIENT_GAP_LOWER)
    _, quotient_pivots = interval_ldlt(shifted_quotient)

    transverse = transverse_comparison(
        box, sizes, degrees, graph["block_types"]
    )
    shifted_transverse = [
        [
            transverse[i][j]
            - (arb(TRANSVERSE_GAP_LOWER) if i == j else arb(0))
            for j in range(6)
        ]
        for i in range(6)
    ]
    _, transverse_pivots = interval_ldlt(shifted_transverse)

    separation = 2 * box[0]
    require(separation > arb(0) and separation < arb.pi(), "not nonsynchronous")

    return {
        "schema": "n80002-arb-application-certificate-v1",
        "status": "CONFIRMED_RIGOROUS",
        "input": {
            "path": "graph_spec.json",
            "sha256": hashlib.sha256(raw).hexdigest(),
            "only_adjacent_graph_spec_read": True,
        },
        "rigor": {
            "library": "python-flint",
            "version": getattr(flint, "__version__", "unknown"),
            "backend": "Arb outward-rounded ball arithmetic",
            "precision_bits": PRECISION_BITS,
            "independent_of_mpmath_verifier": True,
            "shared_code": "generic exact-rational/Arb helpers from ../../formal/verify.py",
        },
        "graph": {
            "simple_unweighted_undirected_by_construction": True,
            "vertex_count": total,
            "class_sizes": sizes,
            "class_degrees": graph["total_degrees"],
            "minimum_degree": graph["minimum_degree"],
            "edge_count": graph["total_edges"],
            "minimum_degree_ratio": f"{graph['minimum_degree']}/{total - 1}",
            "cross_product_excess_over_11_16": 277,
            "biregularity_checked_exactly": True,
            "reflection_quotient_checked_exactly": True,
            "connected": True,
        },
        "equilibrium": {
            "parameter_center_rationals": [frac_text(value) for value in center],
            "box_radius": frac_text(ROOT_RADIUS),
            "center_torque_intervals": [interval_json(value) for value in f_center],
            "preconditioner_rational": fraction_matrix_json(preconditioner),
            "jacobian_box_intervals": matrix_interval_json(jac_box),
            "krawczyk_defect_intervals": matrix_interval_json(defect),
            "contraction_infinity_norm_upper": frac_text(contraction),
            "krawczyk_image_intervals": [interval_json(value) for value in images],
            "strict_inclusion_margins": [frac_text(value) for value in margins],
            "unique_root_in_box": True,
            "all_six_torques_zero_by_exact_reflection": True,
            "nonsynchronous": True,
        },
        "hessian": {
            "jacobian_relation": "DF=-(K/N)H for K>0",
            "quotient_method": "Arb LDLT on H+R-(2/5)I",
            "quotient_nonrotation_gap_lower": frac_text(QUOTIENT_GAP_LOWER),
            "quotient_ldlt_pivots": [interval_json(value) for value in quotient_pivots],
            "transverse_method": "Arb LDLT on the biregular norm comparison matrix",
            "transverse_gap_lower": str(TRANSVERSE_GAP_LOWER),
            "transverse_ldlt_pivots": [
                interval_json(value) for value in transverse_pivots
            ],
            "positive_semidefinite_with_rotation_only_kernel": True,
        },
        "application": {
            "lean_theorem": "KuramotoBlockLimit.cliqueBlowup_preserves_stable_equilibrium",
            "base_hypotheses_certified": True,
            "conclusion": (
                "For every integer m>=1, the replicated phase on G[K_m] is an "
                "exact nonsynchronous equilibrium and its potential Hessian is "
                "positive semidefinite with kernel exactly global rotation."
            ),
            "blowup_order": "80002*m",
            "blowup_minimum_degree": "55019*m-1",
            "exact_excess_over_11_16": "282*m-5",
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=HERE / "arb_certificate.json")
    args = parser.parse_args()
    certificate = certify(HERE / "graph_spec.json")
    args.out.write_text(json.dumps(certificate, indent=2, sort_keys=True) + "\n")
    print(
        "CONFIRMED_RIGOROUS",
        "N=80002",
        "quotient_gap>2/5",
        "transverse_gap>17000",
        "application=all_m>=1",
    )


if __name__ == "__main__":
    main()
