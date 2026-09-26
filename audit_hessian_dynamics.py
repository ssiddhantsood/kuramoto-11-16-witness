#!/usr/bin/env python3
"""Independent, readable audit of the Kuramoto Hessian/stability claim.

This script intentionally uses only the Python standard library.  It reads the
compact graph specification, reconstructs the equitable six-class dynamics,
solves for the displayed reflection-symmetric equilibrium, and checks that the
matrix called the Hessian is exactly the negative Jacobian of the stated
Kuramoto flow (up to the optional positive factor K/N).

It also checks the six class-constant modes and the repository's conservative
lower bound for every class-zero-sum (transverse) mode.  The generated Markdown
report is meant to be read alongside the machine-readable JSON output.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import subprocess
import sys
from pathlib import Path
from typing import Any, Iterable, Sequence


HERE = Path(__file__).resolve().parent
DEFAULT_RECORD = HERE / "records" / "n80002"
PHASE_COEFFICIENTS = (
    (1.0, 0.0, 0.0),
    (-1.0, 0.0, 0.0),
    (0.0, 1.0, 0.0),
    (0.0, -1.0, 0.0),
    (0.0, 0.0, 1.0),
    (0.0, 0.0, -1.0),
)
POSITIVE_CLASSES = (0, 2, 4)


def mat_vec(matrix: Sequence[Sequence[float]], vector: Sequence[float]) -> list[float]:
    return [sum(a * b for a, b in zip(row, vector)) for row in matrix]


def max_abs(values: Iterable[float]) -> float:
    return max(abs(value) for value in values)


def solve_linear(matrix: Sequence[Sequence[float]], rhs: Sequence[float]) -> list[float]:
    """Solve a small dense system using Gaussian elimination with pivoting."""
    n = len(rhs)
    augmented = [list(matrix[i]) + [rhs[i]] for i in range(n)]
    for column in range(n):
        pivot = max(range(column, n), key=lambda row: abs(augmented[row][column]))
        if abs(augmented[pivot][column]) < 1e-15:
            raise ArithmeticError("singular Newton matrix")
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        scale = augmented[column][column]
        augmented[column] = [value / scale for value in augmented[column]]
        for row in range(n):
            if row == column:
                continue
            factor = augmented[row][column]
            augmented[row] = [
                augmented[row][entry] - factor * augmented[column][entry]
                for entry in range(n + 1)
            ]
    return [augmented[row][-1] for row in range(n)]


def jacobi_eigenvalues(matrix: Sequence[Sequence[float]]) -> list[float]:
    """Eigenvalues of a real symmetric 6x6 matrix via Jacobi rotations."""
    a = [list(row) for row in matrix]
    n = len(a)
    scale = max(1.0, max_abs(value for row in a for value in row))
    for _ in range(100 * n * n):
        p, q = max(
            ((i, j) for i in range(n) for j in range(i + 1, n)),
            key=lambda pair: abs(a[pair[0]][pair[1]]),
        )
        if abs(a[p][q]) <= 1e-14 * scale:
            break
        angle = 0.5 * math.atan2(2.0 * a[p][q], a[q][q] - a[p][p])
        cosine, sine = math.cos(angle), math.sin(angle)
        app, aqq, apq = a[p][p], a[q][q], a[p][q]
        for k in range(n):
            if k in (p, q):
                continue
            akp, akq = a[k][p], a[k][q]
            a[k][p] = a[p][k] = cosine * akp - sine * akq
            a[k][q] = a[q][k] = sine * akp + cosine * akq
        a[p][p] = cosine * cosine * app - 2 * sine * cosine * apq + sine * sine * aqq
        a[q][q] = sine * sine * app + 2 * sine * cosine * apq + cosine * cosine * aqq
        a[p][q] = a[q][p] = 0.0
    return sorted(a[i][i] for i in range(n))


def reflected_phases(angles: Sequence[float]) -> list[float]:
    return [angles[0], -angles[0], angles[1], -angles[1], angles[2], -angles[2]]


def block_data(spec: dict[str, Any]) -> tuple[list[int], list[list[int]], dict[tuple[int, int], str]]:
    classes = sorted(spec["classes"], key=lambda item: int(item["id"]))
    sizes = [int(item["size"]) for item in classes]
    degrees = [[0] * 6 for _ in range(6)]
    block_types: dict[tuple[int, int], str] = {}
    for block in spec["cross_blocks"]:
        i, j = map(int, block["classes"])
        kind = str(block["type"])
        block_types[i, j] = kind
        if kind == "absent":
            edges = 0
        elif kind == "complete":
            edges = sizes[i] * sizes[j]
        elif kind == "residue_orbits":
            modulus = int(block["modulus"])
            if modulus != math.gcd(sizes[i], sizes[j]):
                raise ValueError(f"block {(i, j)} has the wrong modulus")
            edges = sizes[i] * sizes[j] * int(block["shift_count"]) // modulus
        else:
            raise ValueError(f"unknown block type {kind!r}")
        if edges % sizes[i] or edges % sizes[j]:
            raise ValueError(f"block {(i, j)} is not biregular")
        degrees[i][j] = edges // sizes[i]
        degrees[j][i] = edges // sizes[j]
    if len(block_types) != 15:
        raise ValueError("all 15 unordered class pairs must be specified")
    for i in range(6):
        for j in range(6):
            if sizes[i] * degrees[i][j] != sizes[j] * degrees[j][i]:
                raise ValueError(f"edge balance fails for classes {i}, {j}")
    return sizes, degrees, block_types


def class_flow(phases: Sequence[float], degrees: Sequence[Sequence[int]]) -> list[float]:
    """F_i(theta) = sum_j d_ij sin(theta_j-theta_i), without K/N."""
    return [
        sum(degrees[i][j] * math.sin(phases[j] - phases[i]) for j in range(6))
        for i in range(6)
    ]


def reduced_torque(angles: Sequence[float], degrees: Sequence[Sequence[int]]) -> list[float]:
    flow = class_flow(reflected_phases(angles), degrees)
    return [flow[i] for i in POSITIVE_CLASSES]


def reduced_jacobian(angles: Sequence[float], degrees: Sequence[Sequence[int]]) -> list[list[float]]:
    phases = reflected_phases(angles)
    result: list[list[float]] = []
    for i in POSITIVE_CLASSES:
        row = []
        for variable in range(3):
            row.append(
                sum(
                    degrees[i][j]
                    * math.cos(phases[j] - phases[i])
                    * (PHASE_COEFFICIENTS[j][variable] - PHASE_COEFFICIENTS[i][variable])
                    for j in range(6)
                )
            )
        result.append(row)
    return result


def solve_equilibrium(seed: Sequence[float], degrees: Sequence[Sequence[int]]) -> tuple[list[float], int]:
    angles = list(seed)
    for iteration in range(20):
        residual = reduced_torque(angles, degrees)
        if max_abs(residual) < 1e-9:
            return angles, iteration
        correction = solve_linear(reduced_jacobian(angles, degrees), [-value for value in residual])
        angles = [angle + delta for angle, delta in zip(angles, correction)]
    raise ArithmeticError("Newton solve did not reach a small residual")


def plain_jacobian(phases: Sequence[float], degrees: Sequence[Sequence[int]]) -> list[list[float]]:
    """Jacobian of the actual six-class phase flow in ordinary coordinates."""
    matrix = [[0.0] * 6 for _ in range(6)]
    for i in range(6):
        for j in range(6):
            if i == j:
                continue
            value = degrees[i][j] * math.cos(phases[j] - phases[i])
            matrix[i][j] = value
            matrix[i][i] -= value
    return matrix


def mass_jacobian(
    jacobian: Sequence[Sequence[float]], sizes: Sequence[int]
) -> list[list[float]]:
    """Jacobian in q_i=sqrt(n_i)*delta-theta_i orthonormal coordinates."""
    return [
        [math.sqrt(sizes[i] / sizes[j]) * jacobian[i][j] for j in range(6)]
        for i in range(6)
    ]


def quotient_hessian(
    phases: Sequence[float], degrees: Sequence[Sequence[int]]
) -> list[list[float]]:
    """Cosine-weighted full-graph Hessian on class-constant directions."""
    matrix = [[0.0] * 6 for _ in range(6)]
    for i in range(6):
        for j in range(i + 1, 6):
            cosine = math.cos(phases[i] - phases[j])
            matrix[i][i] += degrees[i][j] * cosine
            matrix[j][j] += degrees[j][i] * cosine
            value = -math.sqrt(degrees[i][j] * degrees[j][i]) * cosine
            matrix[i][j] = matrix[j][i] = value
    return matrix


def finite_difference_error(
    phases: Sequence[float],
    sizes: Sequence[int],
    degrees: Sequence[Sequence[int]],
    hessian: Sequence[Sequence[float]],
) -> float:
    direction = [0.17, -0.31, 0.29, 0.07, -0.11, 0.23]
    epsilon = 1e-3

    def projected_flow(signed_epsilon: float) -> list[float]:
        shifted = [
            phases[i] + signed_epsilon * direction[i] / math.sqrt(sizes[i])
            for i in range(6)
        ]
        flow = class_flow(shifted, degrees)
        return [math.sqrt(sizes[i]) * flow[i] for i in range(6)]

    plus, minus = projected_flow(epsilon), projected_flow(-epsilon)
    difference = [(plus[i] - minus[i]) / (2 * epsilon) for i in range(6)]
    expected = [-value for value in mat_vec(hessian, direction)]
    return max_abs(difference[i] - expected[i] for i in range(6)) / max_abs(expected)


def transverse_lower_bound(
    phases: Sequence[float],
    sizes: Sequence[int],
    degrees: Sequence[Sequence[int]],
    block_types: dict[tuple[int, int], str],
) -> tuple[list[float], float]:
    """Conservative point bound for all N-6 fiber-zero-sum modes."""
    diagonal = [
        sizes[i]
        + sum(degrees[i][j] * math.cos(phases[i] - phases[j]) for j in range(6))
        for i in range(6)
    ]
    radii = [0.0] * 6
    for (i, j), kind in block_types.items():
        if kind != "residue_orbits":
            continue
        bound = abs(math.cos(phases[i] - phases[j])) * math.sqrt(
            degrees[i][j] * degrees[j][i]
        )
        radii[i] += bound
        radii[j] += bound
    row_lowers = [diagonal[i] - radii[i] for i in range(6)]
    return row_lowers, min(row_lowers)


def fmt(value: float) -> str:
    return f"{value:.12g}"


def run_rigorous_verifier(record: Path) -> dict[str, Any]:
    """Regenerate the rigorous report so the readable layer cannot use stale data."""
    verifier = record / "verify.py"
    report = record / "verification_report.json"
    venv_python = HERE / ".venv" / "bin" / "python"
    python = venv_python if venv_python.exists() else Path(sys.executable)
    completed = subprocess.run(
        [str(python), str(verifier), "--out", str(report)],
        cwd=HERE,
        check=True,
        capture_output=True,
        text=True,
    )
    return {
        "command": f"{python} {verifier} --out {report}",
        "stdout": completed.stdout.strip(),
        "fresh_report_generated": "ACCEPTED" in completed.stdout,
    }


def audit(record: Path, verifier_run: dict[str, Any]) -> dict[str, Any]:
    spec_path = record / "graph_spec.json"
    spec = json.loads(spec_path.read_text())
    spec_sha256 = hashlib.sha256(spec_path.read_bytes()).hexdigest()
    rigorous_report_path = record / "verification_report.json"
    rigorous_report = json.loads(rigorous_report_path.read_text())
    rigorous_input_matches = rigorous_report["input"]["sha256"] == spec_sha256
    root_certificate = rigorous_report["equilibrium"]["krawczyk"]
    nonsynchronous_certificate = rigorous_report["equilibrium"][
        "nonsynchronous_interval_certificate"
    ]
    quotient_certificate = rigorous_report["quotient"]["interval_certificate"]
    quotient_pivot_lowers = [
        float(pivot["lower"])
        for pivot in quotient_certificate["interval_ldlt_pivots"]
    ]
    rigorous_transverse_lower = float(
        rigorous_report["transverse"]["rigorous_absolute_lower"]
    )
    sizes, degrees, block_types = block_data(spec)
    seed = [float(value) for value in spec["equilibrium"]["angle_seed"]]
    angles, iterations = solve_equilibrium(seed, degrees)
    phases = reflected_phases(angles)
    flow = class_flow(phases, degrees)

    ordinary_jacobian = plain_jacobian(phases, degrees)
    transformed_jacobian = mass_jacobian(ordinary_jacobian, sizes)
    hessian = quotient_hessian(phases, degrees)
    identity_error = max_abs(
        transformed_jacobian[i][j] + hessian[i][j]
        for i in range(6)
        for j in range(6)
    )
    symmetry_error = max_abs(
        hessian[i][j] - hessian[j][i] for i in range(6) for j in range(6)
    )
    rotation = [math.sqrt(size) for size in sizes]
    rotation_error = max_abs(mat_vec(hessian, rotation))
    eigenvalues = jacobi_eigenvalues(hessian)
    transverse_rows, transverse_lower = transverse_lower_bound(
        phases, sizes, degrees, block_types
    )
    fd_relative_error = finite_difference_error(phases, sizes, degrees, hessian)
    vertex_count = sum(sizes)
    order_real = sum(sizes[i] * math.cos(phases[i]) for i in range(6))
    order_imag = sum(sizes[i] * math.sin(phases[i]) for i in range(6))
    order_parameter = math.hypot(order_real, order_imag) / vertex_count
    cross_counts = [sum(degrees[i][j] > 0 for j in range(6) if j != i) for i in range(6)]

    checks = {
        "all_six_torques_near_zero": max_abs(flow) < 1e-8,
        "rigorous_verifier_ran_fresh": verifier_run["fresh_report_generated"],
        "rigorous_report_matches_graph_spec": rigorous_input_matches,
        "exact_equilibrium_exists_in_krawczyk_box": (
            root_certificate["strict_interior_inclusion"]
            and root_certificate["unique_root_within_box"]
        ),
        "hessian_is_symmetric": symmetry_error < 1e-10,
        "jacobian_equals_negative_hessian": identity_error < 1e-8,
        "finite_difference_confirms_derivative": fd_relative_error < 1e-7,
        "rotation_is_only_quotient_zero": abs(eigenvalues[0]) < 1e-7 and eigenvalues[1] > 0,
        "interval_certificate_proves_quotient_gap": (
            quotient_certificate["all_interval_pivot_lowers_positive"]
            and quotient_certificate["uniform_over_certified_root_box"]
            and all(lower > 0 for lower in quotient_pivot_lowers)
            and float(
                quotient_certificate["certified_nonrotation_absolute_gap_lower"]
            )
            > 0
        ),
        "all_transverse_modes_have_positive_bound": transverse_lower > 0,
        "outward_rounded_transverse_bound_is_positive": rigorous_transverse_lower > 0,
        "equilibrium_is_nonsynchronous": (
            order_parameter < 0.99
            and nonsynchronous_certificate["nonsynchronous_throughout_root_box"]
            and float(
                nonsynchronous_certificate["phase_separation_interval"]["lower"]
            )
            > 0
            and float(
                nonsynchronous_certificate["phase_separation_interval"]["upper"]
            )
            < math.pi
        ),
    }
    return {
        "schema": "kuramoto-dynamics-hessian-audit-v1",
        "system_contract": {
            "flow": "dtheta_i/dt = (K/N) sum_j A_ij sin(theta_j-theta_i)",
            "assumptions": [
                "K > 0",
                "A is the undirected unweighted graph in graph_spec.json",
                "identical natural frequencies, removed in a rotating frame",
            ],
            "scaling_note": "The code checks K=1 without 1/N; multiplying by positive K/N changes rates, not stability signs.",
        },
        "input": {
            "graph_spec": str(spec_path.relative_to(HERE)),
            "sha256": spec_sha256,
            "rigorous_report": str(rigorous_report_path.relative_to(HERE)),
        },
        "graph": {
            "vertex_count": vertex_count,
            "class_sizes": sizes,
            "cross_class_neighbors": cross_counts,
            "transverse_dimension": vertex_count - 6,
        },
        "equilibrium": {
            "angles": angles,
            "newton_iterations": iterations,
            "maximum_torque_residual": max_abs(flow),
            "order_parameter": order_parameter,
        },
        "derivative_identity": {
            "statement": "J_mass = -H_quotient",
            "maximum_entry_error": identity_error,
            "finite_difference_relative_error": fd_relative_error,
            "hessian_symmetry_error": symmetry_error,
            "rotation_null_residual": rotation_error,
        },
        "spectrum": {
            "quotient_eigenvalues": eigenvalues,
            "quotient_eigenvalues_divided_by_N": [value / vertex_count for value in eigenvalues],
            "transverse_gershgorin_row_lowers": transverse_rows,
            "transverse_point_lower_bound": transverse_lower,
        },
        "rigorous_certificates": {
            "verifier_run": verifier_run,
            "root": root_certificate,
            "nonsynchronous": nonsynchronous_certificate,
            "quotient": quotient_certificate,
            "transverse_absolute_lower": rigorous_transverse_lower,
        },
        "checks": checks,
        "verdict": "PASS" if all(checks.values()) else "FAIL",
        "scope": "Local asymptotic stability modulo global rotation for the stated system; not global attraction and not a certified basin size.",
    }


def markdown(result: dict[str, Any]) -> str:
    graph = result["graph"]
    equilibrium = result["equilibrium"]
    identity = result["derivative_identity"]
    spectrum = result["spectrum"]
    rigorous = result["rigorous_certificates"]
    rows = "\n".join(
        f"| {name.replace('_', ' ')} | {'PASS' if passed else 'FAIL'} |"
        for name, passed in result["checks"].items()
    )
    eigenvalue_rows = "\n".join(
        f"| {index} | {fmt(value)} | {fmt(spectrum['quotient_eigenvalues_divided_by_N'][index])} |"
        for index, value in enumerate(spectrum["quotient_eigenvalues"])
    )
    return f"""# Dynamics-to-Hessian audit

**Verdict: {result['verdict']}**

This audit answers one narrow question: *is the matrix being called the
Hessian the matrix that governs linear stability of the stated Kuramoto
system?* It first reruns `records/n80002/verify.py`, then uses its freshly
generated interval report together with `graph_spec.json`. The readable layer
itself uses only the Python standard library.

## 1. System being checked

The claim is conditional on this precise homogeneous, first-order model:

```text
d theta_i / dt = (K/N) sum_j A_ij sin(theta_j - theta_i),   K > 0.
```

Natural frequencies are identical and removed in a rotating frame. The audit
sets `K=1` and omits `1/N`. That positive factor rescales time and eigenvalues,
but cannot change stability signs.

For

```text
V(theta) = sum over edges (i,j) in E of [1 - cos(theta_i - theta_j)],
```

the flow is `F=-grad(V)`. Therefore its Jacobian must satisfy `DF=-H`, where

```text
H_ii = sum_j A_ij cos(theta_i-theta_j)
H_ij = -A_ij cos(theta_i-theta_j),  i != j.
```

If the intended physical system differs from the displayed ODE, this report
does **not** transfer automatically.

## 2. Graph and equilibrium reconstruction

- Graph specification SHA-256: `{result['input']['sha256']}`
- Vertices: `{graph['vertex_count']}`
- Six class sizes: `{graph['class_sizes']}`
- Other classes connected to each class: `{graph['cross_class_neighbors']}`
- Newton iterations from the short stored seed: `{equilibrium['newton_iterations']}`
- Point maximum torque residual: `{fmt(equilibrium['maximum_torque_residual'])}`
- Order parameter: `{fmt(equilibrium['order_parameter'])}` (not synchronous)

Biregularity is checked exactly with integer edge counts. It is what permits
the six class equations to represent every one of the `{graph['vertex_count']}`
oscillator equations at a class-constant phase assignment.

The point residual is only a numerical sanity check. Existence of an exact
equilibrium is supplied by an outward-rounded Krawczyk inclusion over a
radius-`{rigorous['root']['box_radius_about_short_seed']}` box. Its interval
contraction bound is `{rigorous['root']['remainder_infinity_norm_upper']}`, and
the Krawczyk image lies strictly inside the box.

Nonsynchrony is also certified throughout the box: the phase separation between
classes 0 and 1 lies in
`[{rigorous['nonsynchronous']['phase_separation_interval']['lower']},
{rigorous['nonsynchronous']['phase_separation_interval']['upper']}]`, strictly
between `0` and `pi`.

## 3. Direct Jacobian-versus-Hessian test

The code first differentiates the actual sine flow in ordinary class phase
coordinates. It then changes to the orthonormal coordinates
`q_i=sqrt(n_i) delta-theta_i`. In this basis the matrix should be symmetric and
equal to the negative quotient Hessian.

- Maximum entry error in `J_mass + H`: `{fmt(identity['maximum_entry_error'])}`
- Separate central finite-difference relative error: `{fmt(identity['finite_difference_relative_error'])}`
- Hessian symmetry error: `{fmt(identity['hessian_symmetry_error'])}`
- Rotation-null residual: `{fmt(identity['rotation_null_residual'])}`

The finite-difference test is deliberately separate from the analytic matrix
construction: it perturbs the phases, reevaluates the sine flow, and compares
the observed derivative with `-H`. It is a sanity check, not the rigorous
spectral certificate.

## 4. Quotient and transverse modes

| Quotient mode | Hessian eigenvalue | Eigenvalue after `/N` normalization |
| ---: | ---: | ---: |
{eigenvalue_rows}

The near-zero mode is global rotation. The other five quotient eigenvalues are
positive, so the corresponding flow eigenvalues are negative.

For rigorous rounding control, the verifier evaluates an interval Hessian over
the *entire certified root box*. It adds the exact rotation projector and proves

```text
H + rotation_projector - 0.4 I  is positive definite
```

by interval `LDL^T`. Every interval pivot has a positive lower endpoint. Thus
every non-rotation quotient eigenvalue at the exact enclosed equilibrium is
strictly greater than `{rigorous['quotient']['certified_nonrotation_absolute_gap_lower']}`
(or `{rigorous['quotient']['certified_nonrotation_normalized_gap_lower']}` after
division by `N`).

The quotient alone is not enough: it covers only six directions. The remaining
`{graph['transverse_dimension']}` fiber-zero-sum directions are bounded using:

1. `n_i I` from the clique inside fiber `i`;
2. zero action from complete cross-blocks on zero-sum vectors; and
3. `||A_ij|| <= sqrt(d_ij d_ji)` for every partial biregular block.

The conservative point lower bound for every transverse Hessian mode is
`{fmt(spectrum['transverse_point_lower_bound'])}`. The outward-rounded bound
over the full equilibrium box is `{fmt(rigorous['transverse_absolute_lower'])}`,
which is positive.

## 5. Machine checks

| Check | Result |
| --- | --- |
{rows}

## Conclusion and limits

For the ODE in section 1, the evidence supports a positive-semidefinite full
Hessian whose only zero direction is global phase rotation. Equivalently, the
flow Jacobian has one symmetry zero and is negative in every physical
direction. This implies **local asymptotic stability modulo rotation**.

It does not prove global attraction or a basin radius. The double-precision
numbers above are for readability; the root, quotient sign, and transverse sign
are backed by the packaged outward-rounded certificates summarized here.

Reproduce with:

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements-current-record.txt
.venv/bin/python records/n80002/verify.py
python3 audit_hessian_dynamics.py
```
"""


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record", type=Path, default=DEFAULT_RECORD)
    parser.add_argument("--report", type=Path, default=HERE / "DYNAMICS_HESSIAN_AUDIT.md")
    parser.add_argument("--json", type=Path, default=HERE / "dynamics_hessian_audit.json")
    args = parser.parse_args()
    record = args.record.resolve()
    verifier_run = run_rigorous_verifier(record)
    result = audit(record, verifier_run)
    args.report.write_text(markdown(result))
    args.json.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(f"{result['verdict']}: wrote {args.report.name} and {args.json.name}")
    if result["verdict"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
