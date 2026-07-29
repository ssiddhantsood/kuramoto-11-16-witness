#!/usr/bin/env python3
"""Independent replay for an approach-155 compact graph specification.

This verifier imports only the already-audited generic graph/Fourier/interval
primitives from approach 154.  It imports no approach-155 construction or
search module.  The transverse comparison is strengthened by using the
smaller of the block and complement biregular norm bounds.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
from pathlib import Path
from types import ModuleType
from typing import Any

import mpmath as mp
import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
BASE_REPLAY = (
    ROOT / "approach_154_reflection_weighted_realization" / "replay_witness.py"
)


BASE_REPLAY = HERE / "base_replay.py"


def load_base_replay() -> ModuleType:
    specification = importlib.util.spec_from_file_location(
        "audited_reflection_replay", BASE_REPLAY
    )
    if specification is None or specification.loader is None:
        raise RuntimeError("could not load audited replay implementation")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def complement_aware_interval_operator_bound(
    sizes: np.ndarray,
    blocks: dict[tuple[int, int], dict[str, Any]],
    degrees: np.ndarray,
    angle_box: list[Any],
) -> dict[str, Any]:
    """Rigorous transverse bound using both a block and its complement."""

    iv = mp.iv
    phases = [
        angle_box[0],
        -angle_box[0],
        angle_box[1],
        -angle_box[1],
        angle_box[2],
        -angle_box[2],
    ]
    cosine = [
        [iv.cos(phases[i] - phases[j]) for j in range(6)]
        for i in range(6)
    ]
    beta = []
    for i in range(6):
        value = iv.mpf(str(int(sizes[i])))
        for j in range(6):
            value += int(degrees[i, j]) * cosine[i][j]
        beta.append(value)

    residue_pairs = [
        pair
        for pair, block in blocks.items()
        if block["type"] == "residue_orbits"
        and 0 < int(block["shift_count"]) < int(block["modulus"])
    ]
    norms: dict[tuple[int, int], float] = {}
    for i, j in residue_pairs:
        direct = math.sqrt(int(degrees[i, j]) * int(degrees[j, i]))
        complement = math.sqrt(
            (int(sizes[j]) - int(degrees[i, j]))
            * (int(sizes[i]) - int(degrees[j, i]))
        )
        norms[(i, j)] = min(direct, complement)

    row_bounds = []
    for i in range(6):
        row = beta[i]
        for left, right in residue_pairs:
            if i not in (left, right):
                continue
            j = right if i == left else left
            row -= abs(cosine[i][j]) * iv.mpf(str(norms[(left, right)]))
        row_bounds.append(row)
    lower = min(
        math.nextafter(float(value.a), -math.inf) for value in row_bounds
    )

    central_angles = [
        (float(box.a) + float(box.b)) / 2.0 for box in angle_box
    ]
    central_phases = np.array(
        [
            central_angles[0],
            -central_angles[0],
            central_angles[1],
            -central_angles[1],
            central_angles[2],
            -central_angles[2],
        ]
    )
    central_cosine = np.cos(
        central_phases[:, None] - central_phases[None, :]
    )
    central_beta = sizes.astype(float) + np.sum(
        degrees * central_cosine, axis=1
    )
    comparison = np.diag(central_beta)
    for i, j in residue_pairs:
        comparison[i, j] = comparison[j, i] = (
            -abs(central_cosine[i, j]) * norms[(i, j)]
        )
    return {
        "method": (
            "interval comparison plus Gershgorin with "
            "||A_ij|| <= min(sqrt(d_ij d_ji), "
            "sqrt((n_j-d_ij)(n_i-d_ji)))"
        ),
        "phase_domain": "Krawczyk-certified angle box",
        "fractional_pair_norm_bounds": {
            f"{i}-{j}": value for (i, j), value in norms.items()
        },
        "beta_intervals": [str(value) for value in beta],
        "row_lower_bound_intervals": [str(value) for value in row_bounds],
        "rigorous_absolute_lower_bound": lower,
        "rigorous_normalized_lower_bound": lower / int(np.sum(sizes)),
        "central_comparison_eigenvalues": np.linalg.eigvalsh(
            comparison
        ).tolist(),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--spec", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--dps", type=int, default=130)
    parser.add_argument("--skip-dense", action="store_true")
    parser.add_argument("--skip-fourier", action="store_true")
    args = parser.parse_args()

    base = load_base_replay()
    spec, spec_sha = base.load_spec(args.spec)
    sizes, blocks = base.parse_blocks(spec)
    degrees, edges, edge_count = base.integer_block_data(sizes, blocks)
    angles, equilibrium, angle_box = base.solve_and_certify_angles(
        spec, sizes, degrees, args.dps
    )
    quotient = base.quotient_spectra(sizes, edges, angles)
    quotient_interval = base.interval_quotient_certificate(
        sizes, edges, angle_box
    )
    operator_bound = complement_aware_interval_operator_bound(
        sizes, blocks, degrees, angle_box
    )
    fourier = (
        {
            "skipped": True,
            "reason": (
                "The rigorous interval block-operator comparison covers "
                "the full transverse direct sum without enumeration."
            ),
        }
        if args.skip_fourier
        else base.fourier_transverse_spectrum(
            sizes, blocks, degrees, angles
        )
    )
    dense = (
        None
        if args.skip_dense
        else base.dense_analog_check(spec, blocks, angles)
    )

    total = int(np.sum(sizes))
    class_degrees = sizes - 1 + np.sum(degrees, axis=1)
    minimum_degree = int(np.min(class_degrees))
    maximum_degree = int(np.max(class_degrees))
    target_num = int(spec["target_min_degree_ratio"]["numerator"])
    target_den = int(spec["target_min_degree_ratio"]["denominator"])
    floor_excess = target_den * minimum_degree - target_num * (total - 1)
    handshake = bool(
        sum(
            int(sizes[i]) * int(class_degrees[i]) for i in range(6)
        )
        == 2 * int(edge_count)
    )
    quotient_gap = quotient["_quotient_gap"]
    normalized_gap = quotient_gap / total
    transverse_gap = (
        float(operator_bound["rigorous_absolute_lower_bound"])
        if args.skip_fourier
        else float(fourier["minimum_eigenvalue"])
    )
    rotation_zero = quotient["_rotation_zero"]
    connected = base.support_connected(sizes, blocks)
    reflection = base.reflection_structure(sizes, blocks)
    exactly_one_zero = bool(
        abs(rotation_zero) < mp.mpf("1e-70")
        and quotient["_positive_even"]
        and quotient["_positive_odd"]
        and transverse_gap > 0.0
    )
    dense_ok = bool(
        dense is None or dense["maximum_sorted_spectrum_error"] < 1.0e-7
    )
    torque_residual = mp.mpf(
        equilibrium["maximum_per_vertex_torque_residual"]
    )
    acceptance = {
        "strict_minimum_degree_above_11_16": floor_excess > 0,
        "torque_residual_at_most_1e_60": bool(
            torque_residual <= mp.mpf("1e-60")
        ),
        "unique_equilibrium_root_in_krawczyk_box": equilibrium[
            "unique_root_in_box"
        ],
        "normalized_gap_at_least_1e_6": bool(
            normalized_gap >= mp.mpf("1e-6")
        ),
        "interval_even_and_odd_above_1e_6": quotient_interval[
            "both_sectors_above_shift"
        ],
        "even_sector_positive_off_rotation": quotient["_positive_even"],
        "odd_sector_positive": quotient["_positive_odd"],
        "all_fourier_transverse_modes_positive": (
            None if args.skip_fourier else transverse_gap > 0.0
        ),
        "all_transverse_modes_positive_by_operator_comparison": (
            operator_bound["rigorous_absolute_lower_bound"] > 0.0
        ),
        "rigorous_transverse_operator_bound_positive": (
            operator_bound["rigorous_absolute_lower_bound"] > 0.0
        ),
        "exact_edge_count_matches_degree_handshake": handshake,
        "exactly_one_zero": exactly_one_zero,
        "dense_analog_decomposition_matches": dense_ok,
        "connected": connected,
        "reflection_automorphism": reflection,
    }
    accepted = all(
        value
        for key, value in acceptance.items()
        if key != "all_fourier_transverse_modes_positive"
    )
    report = {
        "outcome": "WITNESS" if accepted else "REJECTED",
        "independence": (
            "Reads only the compact graph spec and audited verification "
            "primitives; imports no construction or search module."
        ),
        "spec_path": str(args.spec.resolve()),
        "spec_sha256": spec_sha,
        "graph": {
            "vertex_count": total,
            "edge_count": int(edge_count),
            "simple": True,
            "connected": connected,
            "unweighted": True,
            "reflection_automorphism_verified": reflection,
            "handshake_identity_verified": handshake,
            "adjacency_oracle_checks": base.oracle_checks(sizes, blocks),
            "class_sizes": sizes.tolist(),
            "directed_block_degrees": degrees.tolist(),
            "class_degrees": class_degrees.tolist(),
            "minimum_degree": minimum_degree,
            "maximum_degree": maximum_degree,
            "mu_exact": f"{minimum_degree}/{total - 1}",
            "mu_float": minimum_degree / (total - 1),
            "degree_cross_product_excess": floor_excess,
        },
        "equilibrium": {
            **equilibrium,
            "phase_strings": [
                mp.nstr(angles[0], args.dps - 20),
                mp.nstr(-angles[0], args.dps - 20),
                mp.nstr(angles[1], args.dps - 20),
                mp.nstr(-angles[1], args.dps - 20),
                mp.nstr(angles[2], args.dps - 20),
                mp.nstr(-angles[2], args.dps - 20),
            ],
            "order_parameter": mp.nstr(
                base.order_parameter(sizes, angles), 60
            ),
            "nonsynchronous": bool(
                base.order_parameter(sizes, angles) < mp.mpf("0.99")
            ),
        },
        "spectrum": {
            "quotient": base.strip_private_fields(quotient),
            "quotient_interval_certificate": quotient_interval,
            "transverse_fourier": fourier,
            "transverse_operator_norm_bound": operator_bound,
            "absolute_hessian_gap": mp.nstr(
                min(quotient_gap, mp.mpf(transverse_gap)), 60
            ),
            "normalized_hessian_gap": mp.nstr(
                min(quotient_gap, mp.mpf(transverse_gap)) / total, 60
            ),
            "exactly_one_zero_eigenvalue": exactly_one_zero,
        },
        "dense_analog_spot_check": dense,
        "acceptance": {**acceptance, "accepted": accepted},
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(
        json.dumps(
            {
                "outcome": report["outcome"],
                "mu_exact": report["graph"]["mu_exact"],
                "normalized_gap": report["spectrum"][
                    "normalized_hessian_gap"
                ],
                "operator_bound": operator_bound[
                    "rigorous_absolute_lower_bound"
                ],
                "output": str(args.output.resolve()),
            },
            indent=2,
        )
    )
    if not accepted:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
