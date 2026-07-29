#!/usr/bin/env python3
"""Independent replay of a compact reflection-paired graph witness.

This verifier reads only ``final_graph_spec.json``.  It does not import the
search or construction modules.  It reconstructs exact graph arithmetic,
resolves and interval-certifies the equilibrium, checks the real even/odd
quotient sectors, enumerates every transverse Fourier mode, and supplies an
independent interval operator-norm bound for the full transverse sector.
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


HERE = Path(__file__).resolve().parent
DEFAULT_SPEC = HERE / "artifacts" / "final_graph_spec.json"
DEFAULT_OUTPUT = HERE / "artifacts" / "witness_certificate.json"


def load_spec(path: Path) -> tuple[dict[str, Any], str]:
    raw = path.read_bytes()
    spec = json.loads(raw)
    if spec.get("schema") != "reflection-paired-circulant-witness-v1":
        raise ValueError("unsupported graph specification")
    return spec, hashlib.sha256(raw).hexdigest()


def parse_blocks(
    spec: dict[str, Any],
) -> tuple[np.ndarray, dict[tuple[int, int], dict[str, Any]]]:
    classes = list(spec["classes"])
    sizes = np.array([int(entry["size"]) for entry in classes], dtype=np.int64)
    if len(sizes) != 6 or np.any(sizes <= 1):
        raise ValueError("certificate requires six clique fibers of size > 1")
    if not (
        sizes[0] == sizes[1]
        and sizes[2] == sizes[3]
        and sizes[4] == sizes[5]
    ):
        raise ValueError("reflected class sizes are not paired")
    if int(sizes.sum()) != int(spec["vertex_count"]):
        raise ValueError("class sizes do not sum to vertex_count")

    blocks: dict[tuple[int, int], dict[str, Any]] = {}
    for raw_block in spec["cross_blocks"]:
        block = dict(raw_block)
        i, j = map(int, block["classes"])
        if not 0 <= i < j < 6:
            raise ValueError(f"invalid cross-block pair {(i, j)}")
        if (i, j) in blocks:
            raise ValueError(f"duplicate cross-block pair {(i, j)}")
        kind = block["type"]
        if kind not in {"absent", "complete", "residue_orbits"}:
            raise ValueError(f"unknown cross-block type {kind!r}")
        if kind == "residue_orbits":
            modulus = int(block["modulus"])
            expected = math.gcd(int(sizes[i]), int(sizes[j]))
            start = int(block["shift_start"])
            count = int(block["shift_count"])
            if modulus != expected:
                raise ValueError(
                    f"block {(i, j)} modulus {modulus} != gcd {expected}"
                )
            if not (0 <= start < modulus and 0 <= count <= modulus):
                raise ValueError(f"invalid residue interval in block {(i, j)}")
        blocks[(i, j)] = block

    expected_pairs = {(i, j) for i in range(6) for j in range(i + 1, 6)}
    if set(blocks) != expected_pairs:
        missing = sorted(expected_pairs - set(blocks))
        extra = sorted(set(blocks) - expected_pairs)
        raise ValueError(f"cross-block coverage mismatch: {missing=}, {extra=}")
    return sizes, blocks


def integer_block_data(
    sizes: np.ndarray, blocks: dict[tuple[int, int], dict[str, Any]]
) -> tuple[np.ndarray, np.ndarray, int]:
    """Return exact directed block degrees, block edges, and total edges."""

    degrees = np.zeros((6, 6), dtype=np.int64)
    edge_totals = np.zeros((6, 6), dtype=object)
    total_edges = sum(
        int(size) * (int(size) - 1) // 2 for size in sizes
    )
    for (i, j), block in blocks.items():
        kind = block["type"]
        if kind == "absent":
            edges = 0
        elif kind == "complete":
            edges = int(sizes[i]) * int(sizes[j])
        else:
            modulus = int(block["modulus"])
            count = int(block["shift_count"])
            edges = int(sizes[i]) * int(sizes[j]) // modulus * count
        if edges % int(sizes[i]) or edges % int(sizes[j]):
            raise ValueError(f"block {(i, j)} is not biregular")
        degrees[i, j] = edges // int(sizes[i])
        degrees[j, i] = edges // int(sizes[j])
        edge_totals[i, j] = edge_totals[j, i] = edges
        total_edges += edges
    return degrees, edge_totals, total_edges


def edge_predicate(
    class_u: int,
    local_u: int,
    class_v: int,
    local_v: int,
    sizes: np.ndarray,
    blocks: dict[tuple[int, int], dict[str, Any]],
) -> bool:
    """Exact adjacency oracle for two class-local vertex labels."""

    if not (0 <= class_u < 6 and 0 <= class_v < 6):
        raise IndexError("class index out of range")
    if not (0 <= local_u < sizes[class_u] and 0 <= local_v < sizes[class_v]):
        raise IndexError("local vertex label out of range")
    if class_u == class_v:
        return local_u != local_v
    if class_u > class_v:
        class_u, class_v = class_v, class_u
        local_u, local_v = local_v, local_u
    block = blocks[(class_u, class_v)]
    if block["type"] == "absent":
        return False
    if block["type"] == "complete":
        return True
    modulus = int(block["modulus"])
    start = int(block["shift_start"])
    count = int(block["shift_count"])
    return (local_v - local_u - start) % modulus < count


def support_connected(
    sizes: np.ndarray, blocks: dict[tuple[int, int], dict[str, Any]]
) -> bool:
    """Connectivity follows from clique fibers and connected class support."""

    support = [[] for _ in range(6)]
    for (i, j), block in blocks.items():
        nonempty = block["type"] == "complete" or (
            block["type"] == "residue_orbits"
            and int(block["shift_count"]) > 0
        )
        if nonempty:
            support[i].append(j)
            support[j].append(i)
    seen = {0}
    pending = [0]
    while pending:
        current = pending.pop()
        for neighbor in support[current]:
            if neighbor not in seen:
                seen.add(neighbor)
                pending.append(neighbor)
    return bool(np.all(sizes > 0) and len(seen) == 6)


def reflection_structure(
    sizes: np.ndarray,
    blocks: dict[tuple[int, int], dict[str, Any]],
) -> bool:
    """Verify rho(c,u)=(pairmate(c),-u) is a graph automorphism."""

    involution = {0: 1, 1: 0, 2: 3, 3: 2, 4: 5, 5: 4}
    if any(sizes[i] != sizes[involution[i]] for i in range(6)):
        return False

    for i in range(6):
        for j in range(i + 1, 6):
            source = blocks[(i, j)]
            ii, jj = involution[i], involution[j]
            order_preserved = ii < jj
            target = blocks[(ii, jj) if order_preserved else (jj, ii)]
            if source["type"] != target["type"]:
                return False
            if source["type"] != "residue_orbits":
                continue
            modulus = int(source["modulus"])
            count = int(source["shift_count"])
            if (
                int(target["modulus"]) != modulus
                or int(target["shift_count"]) != count
            ):
                return False
            start = int(source["shift_start"])
            expected_start = (
                -(start + count - 1) if order_preserved else start
            ) % modulus
            if int(target["shift_start"]) != expected_start:
                return False
    return True


def torque_system(
    angles: list[Any],
    sizes: np.ndarray,
    degrees: np.ndarray,
    context: Any = mp,
) -> list[Any]:
    t1, t2, t3 = angles
    phases = [t1, -t1, t2, -t2, t3, -t3]
    result = []
    for i in (0, 2, 4):
        value = context.mpf("0")
        for j in range(6):
            value += int(degrees[i, j]) * context.sin(phases[i] - phases[j])
        result.append(value)
    return result


def torque_jacobian(
    angles: list[Any],
    sizes: np.ndarray,
    degrees: np.ndarray,
    context: Any = mp,
) -> list[list[Any]]:
    """Differentiate the three positive-class torques directly."""

    del sizes
    t1, t2, t3 = angles
    phases = [t1, -t1, t2, -t2, t3, -t3]
    phase_derivatives = [
        [1, 0, 0],
        [-1, 0, 0],
        [0, 1, 0],
        [0, -1, 0],
        [0, 0, 1],
        [0, 0, -1],
    ]
    out = [[context.mpf("0") for _ in range(3)] for _ in range(3)]
    for row, i in enumerate((0, 2, 4)):
        for j in range(6):
            coefficient = int(degrees[i, j]) * context.cos(
                phases[i] - phases[j]
            )
            for column in range(3):
                out[row][column] += coefficient * (
                    phase_derivatives[i][column]
                    - phase_derivatives[j][column]
                )
    return out


def solve_and_certify_angles(
    spec: dict[str, Any],
    sizes: np.ndarray,
    degrees: np.ndarray,
    dps: int,
) -> tuple[list[mp.mpf], dict[str, Any], list[Any]]:
    """High-precision solve plus a Krawczyk unique-root inclusion."""

    mp.mp.dps = dps
    mp.iv.dps = max(80, dps - 30)
    seed = [mp.mpf(value) for value in spec["equilibrium"]["angle_seed"]]

    def f0(t1: mp.mpf, t2: mp.mpf, t3: mp.mpf) -> mp.mpf:
        return torque_system([t1, t2, t3], sizes, degrees)[0]

    def f1(t1: mp.mpf, t2: mp.mpf, t3: mp.mpf) -> mp.mpf:
        return torque_system([t1, t2, t3], sizes, degrees)[1]

    def f2(t1: mp.mpf, t2: mp.mpf, t3: mp.mpf) -> mp.mpf:
        return torque_system([t1, t2, t3], sizes, degrees)[2]

    angles = list(
        mp.findroot(
            (f0, f1, f2),
            tuple(seed),
            solver="mdnewton",
            tol=mp.mpf(10) ** (-(dps - 20)),
            maxsteps=100,
        )
    )
    residual_values = torque_system(angles, sizes, degrees)
    residual = max(abs(value) for value in residual_values)

    radius = mp.mpf("1e-45")
    angle_box = [
        mp.iv.mpf([str(value - radius), str(value + radius)])
        for value in angles
    ]
    jacobian_center = mp.matrix(
        torque_jacobian(angles, sizes, degrees)
    )
    inverse = jacobian_center**-1
    function_center = mp.matrix(residual_values)
    newton_center = mp.matrix(angles) - inverse * function_center
    jacobian_interval = torque_jacobian(
        angle_box, sizes, degrees, context=mp.iv
    )
    rounding_guard = mp.mpf("1e-100")
    remainder = [[mp.iv.mpf("0") for _ in range(3)] for _ in range(3)]
    for i in range(3):
        for j in range(3):
            value = mp.iv.mpf("1" if i == j else "0")
            for k in range(3):
                inverse_interval = mp.iv.mpf(
                    [
                        str(inverse[i, k] - rounding_guard),
                        str(inverse[i, k] + rounding_guard),
                    ]
                )
                value -= (
                    inverse_interval * jacobian_interval[k][j]
                )
            remainder[i][j] = value
    delta = mp.iv.mpf([str(-radius), str(radius)])
    inclusions = []
    included = True
    for i in range(3):
        value = mp.iv.mpf(
            [
                str(newton_center[i] - rounding_guard),
                str(newton_center[i] + rounding_guard),
            ]
        )
        for j in range(3):
            value += remainder[i][j] * delta
        inner = mp.iv.mpf(
            [
                str(angles[i] - radius / 2),
                str(angles[i] + radius / 2),
            ]
        )
        contained = bool(value in inner)
        included = included and contained
        inclusions.append(
            {
                "krawczyk_interval": str(value),
                "root_box": str(angle_box[i]),
                "inner_half_box": str(inner),
                "strictly_contained": contained,
            }
        )
    if not included:
        raise RuntimeError("Krawczyk root inclusion failed")

    stored = [
        mp.mpf(value) for value in spec["equilibrium"]["angle_strings"]
    ]
    phase_match = max(abs(left - right) for left, right in zip(angles, stored))
    certificate = {
        "angles": [mp.nstr(value, dps - 20) for value in angles],
        "maximum_per_vertex_torque_residual": mp.nstr(residual, 40),
        "maximum_angle_difference_from_stored": mp.nstr(phase_match, 40),
        "krawczyk_radius": str(radius),
        "krawczyk_inclusions": inclusions,
        "unique_root_in_box": included,
    }
    return angles, certificate, angle_box


def quotient_spectra(
    sizes: np.ndarray,
    edge_totals: np.ndarray,
    angles: list[mp.mpf],
) -> dict[str, Any]:
    """Reconstruct full, reflection-even, and reflection-odd quotients."""

    phases = [
        angles[0],
        -angles[0],
        angles[1],
        -angles[1],
        angles[2],
        -angles[2],
    ]
    quotient = mp.matrix(6)
    for i in range(6):
        for j in range(6):
            quotient[i, j] = mp.mpf("0")
    for i in range(6):
        for j in range(i + 1, 6):
            conductance = mp.mpf(int(edge_totals[i, j])) * mp.cos(
                phases[i] - phases[j]
            )
            normalized = conductance / mp.sqrt(
                mp.mpf(int(sizes[i] * sizes[j]))
            )
            quotient[i, i] += conductance / int(sizes[i])
            quotient[j, j] += conductance / int(sizes[j])
            quotient[i, j] = quotient[j, i] = -normalized

    even_transform = mp.matrix(6, 3)
    odd_transform = mp.matrix(6, 3)
    root2 = mp.sqrt(2)
    for row in range(6):
        for column in range(3):
            even_transform[row, column] = mp.mpf("0")
            odd_transform[row, column] = mp.mpf("0")
    for pair in range(3):
        plus, minus = 2 * pair, 2 * pair + 1
        even_transform[plus, pair] = 1 / root2
        even_transform[minus, pair] = 1 / root2
        odd_transform[plus, pair] = 1 / root2
        odd_transform[minus, pair] = -1 / root2
    even = even_transform.T * quotient * even_transform
    odd = odd_transform.T * quotient * odd_transform
    even_raw_values = mp.eigsy(even, eigvals_only=True)
    odd_values = mp.eigsy(odd, eigvals_only=True)
    full_values = mp.eigsy(quotient, eigvals_only=True)

    total = int(sizes.sum())
    rotation = mp.matrix(
        [
            mp.sqrt(mp.mpf(2 * int(sizes[0])) / total),
            mp.sqrt(mp.mpf(2 * int(sizes[2])) / total),
            mp.sqrt(mp.mpf(2 * int(sizes[4])) / total),
        ]
    )
    even_lifted = even + 2 * total * rotation * rotation.T
    even_lifted_values = mp.eigsy(even_lifted, eigvals_only=True)
    even_gap = min(even_lifted_values)
    odd_gap = min(odd_values)
    quotient_gap = min(even_gap, odd_gap)
    return {
        "full_absolute": [mp.nstr(value, 60) for value in full_values],
        "full_normalized": [
            mp.nstr(value / total, 60) for value in full_values
        ],
        "even_raw_absolute": [
            mp.nstr(value, 60) for value in even_raw_values
        ],
        "even_with_rotation_lift_absolute": [
            mp.nstr(value, 60) for value in even_lifted_values
        ],
        "odd_absolute": [mp.nstr(value, 60) for value in odd_values],
        "even_nonrotation_gap_absolute": mp.nstr(even_gap, 60),
        "odd_gap_absolute": mp.nstr(odd_gap, 60),
        "quotient_gap_absolute": mp.nstr(quotient_gap, 60),
        "quotient_gap_normalized": mp.nstr(quotient_gap / total, 60),
        "rotation_zero_absolute": mp.nstr(even_raw_values[0], 40),
        "_quotient_gap": quotient_gap,
        "_rotation_zero": even_raw_values[0],
        "_positive_even": bool(even_gap > 0),
        "_positive_odd": bool(odd_gap > 0),
    }


def interval_ldl(matrix: list[list[Any]], shift: mp.mpf) -> dict[str, Any]:
    """Prove an interval-symmetric matrix exceeds ``shift * I``."""

    iv = mp.iv
    count = len(matrix)
    lower = [[iv.mpf("0") for _ in range(count)] for _ in range(count)]
    diagonal = [iv.mpf("0") for _ in range(count)]
    for i in range(count):
        lower[i][i] = iv.mpf("1")
        pivot = matrix[i][i] - iv.mpf(str(shift))
        for k in range(i):
            pivot -= lower[i][k] * lower[i][k] * diagonal[k]
        if float(pivot.a) <= 0:
            raise RuntimeError(f"nonpositive interval LDL pivot {i}: {pivot}")
        diagonal[i] = pivot
        for j in range(i + 1, count):
            numerator = matrix[j][i]
            for k in range(i):
                numerator -= (
                    lower[j][k] * lower[i][k] * diagonal[k]
                )
            lower[j][i] = numerator / diagonal[i]
    return {
        "shift_absolute": str(shift),
        "pivot_intervals": [str(value) for value in diagonal],
        "all_pivots_positive": True,
    }


def interval_quotient_certificate(
    sizes: np.ndarray,
    edge_totals: np.ndarray,
    angle_box: list[Any],
) -> dict[str, Any]:
    """Interval-LDL proof of both real quotient sectors above 1e-6 N."""

    iv = mp.iv
    phases = [
        angle_box[0],
        -angle_box[0],
        angle_box[1],
        -angle_box[1],
        angle_box[2],
        -angle_box[2],
    ]
    quotient = [[iv.mpf("0") for _ in range(6)] for _ in range(6)]
    for i in range(6):
        for j in range(i + 1, 6):
            conductance = iv.mpf(str(int(edge_totals[i, j]))) * iv.cos(
                phases[i] - phases[j]
            )
            quotient[i][i] += conductance / int(sizes[i])
            quotient[j][j] += conductance / int(sizes[j])
            scale = iv.sqrt(
                iv.mpf(str(int(sizes[i]) * int(sizes[j])))
            )
            quotient[i][j] = quotient[j][i] = -conductance / scale

    even = [[iv.mpf("0") for _ in range(3)] for _ in range(3)]
    odd = [[iv.mpf("0") for _ in range(3)] for _ in range(3)]
    for i in range(3):
        pi, mi = 2 * i, 2 * i + 1
        for j in range(3):
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

    total = int(sizes.sum())
    pair_sizes = [int(sizes[0]), int(sizes[2]), int(sizes[4])]
    rotation = [
        iv.sqrt(iv.mpf(str(2 * size)) / total) for size in pair_sizes
    ]
    even_lifted = [
        [
            even[i][j]
            + iv.mpf(str(2 * total)) * rotation[i] * rotation[j]
            for j in range(3)
        ]
        for i in range(3)
    ]
    shift = mp.mpf(total) * mp.mpf("1e-6")
    return {
        "normalized_shift": "0.000001",
        "even_with_rotation_lift": interval_ldl(even_lifted, shift),
        "odd": interval_ldl(odd, shift),
        "both_sectors_above_shift": True,
    }


def _geometric_sum(angle: float, start: int, count: int) -> complex:
    if count == 0:
        return 0.0j
    denominator = np.expm1(1j * angle)
    if abs(denominator) < 1e-14:
        return complex(count)
    return (
        np.exp(1j * start * angle)
        * np.expm1(1j * count * angle)
        / denominator
    )


def fourier_transverse_spectrum(
    sizes: np.ndarray,
    blocks: dict[tuple[int, int], dict[str, Any]],
    degrees: np.ndarray,
    angles: list[mp.mpf],
    collect: bool = False,
) -> dict[str, Any]:
    """Enumerate all N-6 nonconstant cyclic characters."""

    phases = np.array(
        [
            float(angles[0]),
            -float(angles[0]),
            float(angles[1]),
            -float(angles[1]),
            float(angles[2]),
            -float(angles[2]),
        ]
    )
    cosine = np.cos(phases[:, None] - phases[None, :])
    beta = sizes.astype(float) + np.sum(degrees * cosine, axis=1)
    period = math.lcm(*(int(value) for value in np.unique(sizes)))
    characters = np.unique(
        np.concatenate(
            [
                np.arange(int(size), dtype=np.int64) * (period // int(size))
                for size in np.unique(sizes)
            ]
        )
    )
    characters = characters[characters != 0]
    residue_blocks = [
        (pair, block)
        for pair, block in blocks.items()
        if block["type"] == "residue_orbits"
        and 0 < int(block["shift_count"]) < int(block["modulus"])
    ]

    mode_count = 0
    minimum = math.inf
    minimum_character: int | None = None
    minimum_active: list[int] | None = None
    minimum_values: list[float] | None = None
    all_values: list[float] = []
    for raw_character in characters:
        character = int(raw_character)
        active = np.flatnonzero((character * sizes) % period == 0)
        mode_count += len(active)
        positions = {int(cls): index for index, cls in enumerate(active)}
        matrix = np.diag(beta[active]).astype(np.complex128)
        angle = 2 * np.pi * character / period
        for (i, j), block in residue_blocks:
            if i not in positions or j not in positions:
                continue
            modulus = int(block["modulus"])
            start = int(block["shift_start"])
            count = int(block["shift_count"])
            multiplier = (
                math.sqrt(int(sizes[i]) * int(sizes[j]))
                / modulus
                * _geometric_sum(angle, start, count)
            )
            value = -cosine[i, j] * multiplier
            ii, jj = positions[i], positions[j]
            matrix[ii, jj] = value
            matrix[jj, ii] = value.conjugate()
        values = np.linalg.eigvalsh(matrix)
        if collect:
            all_values.extend(float(value) for value in values)
        if float(values[0]) < minimum:
            minimum = float(values[0])
            minimum_character = character
            minimum_active = [int(value) for value in active]
            minimum_values = [float(value) for value in values]

    expected = int(sizes.sum()) - 6
    if mode_count != expected:
        raise AssertionError(
            f"enumerated {mode_count} transverse modes, expected {expected}"
        )
    out: dict[str, Any] = {
        "common_cyclic_period": period,
        "nonzero_character_count": int(len(characters)),
        "mode_count": mode_count,
        "minimum_eigenvalue": minimum,
        "minimum_normalized_eigenvalue": minimum / int(sizes.sum()),
        "minimum_character_index": minimum_character,
        "minimum_character_angle": (
            2 * math.pi * minimum_character / period
            if minimum_character is not None
            else None
        ),
        "minimum_active_classes": minimum_active,
        "minimum_block_eigenvalues": minimum_values,
        "fiber_diagonal_terms": [float(value) for value in beta],
    }
    if collect:
        out["all_values"] = all_values
    return out


def interval_operator_bound(
    sizes: np.ndarray,
    blocks: dict[tuple[int, int], dict[str, Any]],
    degrees: np.ndarray,
    angle_box: list[Any],
) -> dict[str, Any]:
    """Rigorous block-operator bound over the Krawczyk phase box.

    On the direct sum of classwise zero-sum subspaces, a residue block
    ``A_ij`` has norm at most ``sqrt(d_ij*d_ji)``.  Complete blocks vanish.
    The comparison operator therefore has diagonal beta_i and off-diagonal
    magnitudes at most ``|cos(phi_i-phi_j)| sqrt(d_ij*d_ji)``.  Interval
    Gershgorin lower bounds on this comparison matrix cover every transverse
    mode, independently of the Fourier enumeration.
    """

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
    row_bounds = []
    for i in range(6):
        row = beta[i]
        for left, right in residue_pairs:
            if i not in (left, right):
                continue
            j = right if i == left else left
            norm = iv.sqrt(
                iv.mpf(
                    str(int(degrees[i, j]) * int(degrees[j, i]))
                )
            )
            row -= abs(cosine[i][j]) * norm
        row_bounds.append(row)
    lower_floats = [
        math.nextafter(float(value.a), -math.inf) for value in row_bounds
    ]
    lower = min(lower_floats)

    # Also report the central numerical comparison spectrum.
    central_angles = [
        (float(box.a) + float(box.b)) / 2
        for box in angle_box
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
        comparison[i, j] = comparison[j, i] = -abs(
            central_cosine[i, j]
        ) * math.sqrt(int(degrees[i, j]) * int(degrees[j, i]))
    comparison_values = np.linalg.eigvalsh(comparison)
    return {
        "method": (
            "interval block-operator comparison plus Gershgorin: "
            "||A_ij|| <= sqrt(d_ij*d_ji)"
        ),
        "phase_domain": "Krawczyk-certified angle box",
        "beta_intervals": [str(value) for value in beta],
        "row_lower_bound_intervals": [str(value) for value in row_bounds],
        "rigorous_absolute_lower_bound": lower,
        "rigorous_normalized_lower_bound": lower / int(sizes.sum()),
        "central_comparison_eigenvalues": [
            float(value) for value in comparison_values
        ],
    }


def scaled_analog(
    spec: dict[str, Any],
    blocks: dict[tuple[int, int], dict[str, Any]],
) -> dict[str, Any]:
    sizes = np.array([90, 90, 45, 45, 120, 120], dtype=np.int64)
    out = {
        "schema": spec["schema"],
        "vertex_count": int(sizes.sum()),
        "classes": [
            {"id": index, "size": int(size)}
            for index, size in enumerate(sizes)
        ],
        "cross_blocks": [],
    }
    for (i, j), block in blocks.items():
        copied: dict[str, Any] = {
            "classes": [i, j],
            "type": block["type"],
        }
        if block["type"] == "residue_orbits":
            old_modulus = int(block["modulus"])
            density = int(block["shift_count"]) / old_modulus
            modulus = math.gcd(int(sizes[i]), int(sizes[j]))
            count = min(modulus, max(0, round(density * modulus)))
            copied.update(
                {
                    "modulus": modulus,
                    "shift_start": 0,
                    "shift_count": count,
                }
            )
        out["cross_blocks"].append(copied)
    return out


def numeric_quotient_values(
    sizes: np.ndarray, edges: np.ndarray, phases: np.ndarray
) -> np.ndarray:
    conductance = np.zeros((6, 6))
    for i in range(6):
        for j in range(i + 1, 6):
            value = float(edges[i, j]) * math.cos(phases[i] - phases[j])
            conductance[i, j] = conductance[j, i] = value
    hessian = np.diag(np.sum(conductance, axis=1)) - conductance
    quotient = hessian / np.sqrt(np.outer(sizes, sizes))
    return np.linalg.eigvalsh(0.5 * (quotient + quotient.T))


def dense_analog_check(
    spec: dict[str, Any],
    blocks: dict[tuple[int, int], dict[str, Any]],
    angles: list[mp.mpf],
) -> dict[str, Any]:
    analog = scaled_analog(spec, blocks)
    sizes, analog_blocks = parse_blocks(analog)
    degrees, edges, _ = integer_block_data(sizes, analog_blocks)
    starts = np.concatenate([[0], np.cumsum(sizes)]).astype(int)
    total = int(sizes.sum())
    adjacency = np.zeros((total, total))
    for i, size in enumerate(sizes):
        low, high = starts[i], starts[i + 1]
        adjacency[low:high, low:high] = 1
    np.fill_diagonal(adjacency, 0)
    for (i, j), block in analog_blocks.items():
        ni, nj = int(sizes[i]), int(sizes[j])
        if block["type"] == "absent":
            array = np.zeros((ni, nj))
        elif block["type"] == "complete":
            array = np.ones((ni, nj))
        else:
            u = np.arange(ni)[:, None]
            v = np.arange(nj)[None, :]
            array = (
                (v - u - int(block["shift_start"])) % int(block["modulus"])
                < int(block["shift_count"])
            ).astype(float)
        adjacency[
            starts[i] : starts[i + 1], starts[j] : starts[j + 1]
        ] = array
        adjacency[
            starts[j] : starts[j + 1], starts[i] : starts[i + 1]
        ] = array.T
    if not np.array_equal(adjacency, adjacency.T):
        raise AssertionError("analog adjacency is asymmetric")
    if np.any(np.diag(adjacency)):
        raise AssertionError("analog adjacency has self-loops")

    phases = np.array(
        [
            float(angles[0]),
            -float(angles[0]),
            float(angles[1]),
            -float(angles[1]),
            float(angles[2]),
            -float(angles[2]),
        ]
    )
    lifted = np.repeat(phases, sizes)
    cosine = np.cos(lifted[:, None] - lifted[None, :])
    conductance = adjacency * cosine
    hessian = np.diag(np.sum(conductance, axis=1)) - conductance
    dense_values = np.linalg.eigvalsh(hessian)
    quotient_values = numeric_quotient_values(sizes, edges, phases)
    fourier = fourier_transverse_spectrum(
        sizes, analog_blocks, degrees, angles, collect=True
    )
    decomposed = np.sort(
        np.concatenate(
            [quotient_values, np.array(fourier.pop("all_values"))]
        )
    )
    maximum_error = float(np.max(np.abs(dense_values - decomposed)))
    return {
        "class_sizes": [int(value) for value in sizes],
        "vertex_count": total,
        "maximum_sorted_spectrum_error": maximum_error,
        "dense_smallest_eight": [float(value) for value in dense_values[:8]],
        "decomposed_smallest_eight": [
            float(value) for value in decomposed[:8]
        ],
    }


def oracle_checks(
    sizes: np.ndarray, blocks: dict[tuple[int, int], dict[str, Any]]
) -> int:
    checks = 0
    for (i, j), block in blocks.items():
        ni, nj = int(sizes[i]), int(sizes[j])
        probes = [
            (0, 0),
            (ni - 1, nj - 1),
            (ni // 2, nj // 3),
            (ni // 3, nj // 2),
        ]
        if block["type"] == "residue_orbits":
            modulus = int(block["modulus"])
            start = int(block["shift_start"])
            count = int(block["shift_count"])
            for delta in {0, count - 1, count, modulus - 1}:
                probes.append((0, (start + delta) % modulus))
        for u, v in probes:
            forward = edge_predicate(i, u, j, v, sizes, blocks)
            reverse = edge_predicate(j, v, i, u, sizes, blocks)
            if forward != reverse:
                raise AssertionError(f"asymmetric oracle result in {(i, j)}")
            checks += 1
    return checks


def order_parameter(sizes: np.ndarray, angles: list[mp.mpf]) -> mp.mpf:
    phases = [
        angles[0],
        -angles[0],
        angles[1],
        -angles[1],
        angles[2],
        -angles[2],
    ]
    return abs(
        mp.fsum(
            mp.mpf(int(sizes[i])) * mp.exp(1j * phases[i])
            for i in range(6)
        )
        / int(sizes.sum())
    )


def strip_private_fields(values: dict[str, Any]) -> dict[str, Any]:
    return {key: value for key, value in values.items() if not key.startswith("_")}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--spec", type=Path, default=DEFAULT_SPEC)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--dps", type=int, default=130)
    parser.add_argument("--skip-dense", action="store_true")
    args = parser.parse_args()

    spec, spec_sha = load_spec(args.spec)
    sizes, blocks = parse_blocks(spec)
    degrees, edges, edge_count = integer_block_data(sizes, blocks)
    angles, equilibrium, angle_box = solve_and_certify_angles(
        spec, sizes, degrees, args.dps
    )
    quotient = quotient_spectra(sizes, edges, angles)
    quotient_interval = interval_quotient_certificate(
        sizes, edges, angle_box
    )
    fourier = fourier_transverse_spectrum(sizes, blocks, degrees, angles)
    operator_bound = interval_operator_bound(
        sizes, blocks, degrees, angle_box
    )
    dense = (
        None
        if args.skip_dense
        else dense_analog_check(spec, blocks, angles)
    )

    total = int(sizes.sum())
    class_degrees = sizes - 1 + np.sum(degrees, axis=1)
    handshake_ok = bool(
        sum(
            int(sizes[i]) * int(class_degrees[i]) for i in range(6)
        )
        == 2 * int(edge_count)
    )
    minimum_degree = int(class_degrees.min())
    maximum_degree = int(class_degrees.max())
    target_num = int(spec["target_min_degree_ratio"]["numerator"])
    target_den = int(spec["target_min_degree_ratio"]["denominator"])
    floor_excess = (
        target_den * minimum_degree - target_num * (total - 1)
    )
    quotient_gap = quotient["_quotient_gap"]
    normalized_gap = quotient_gap / total
    transverse_gap = fourier["minimum_eigenvalue"]
    full_gap = min(float(quotient_gap), transverse_gap)
    connected = support_connected(sizes, blocks)
    reflection_ok = reflection_structure(sizes, blocks)
    oracle_count = oracle_checks(sizes, blocks)
    rotation_zero = quotient["_rotation_zero"]
    exactly_one_zero = bool(
        abs(rotation_zero) < mp.mpf("1e-70")
        and quotient["_positive_even"]
        and quotient["_positive_odd"]
        and transverse_gap > 0
    )
    torque_residual = mp.mpf(
        equilibrium["maximum_per_vertex_torque_residual"]
    )
    dense_ok = bool(
        dense is None or dense["maximum_sorted_spectrum_error"] < 1e-7
    )
    accepted = bool(
        floor_excess > 0
        and torque_residual <= mp.mpf("1e-60")
        and equilibrium["unique_root_in_box"]
        and normalized_gap >= mp.mpf("1e-6")
        and quotient_interval["both_sectors_above_shift"]
        and exactly_one_zero
        and connected
        and reflection_ok
        and handshake_ok
        and operator_bound["rigorous_absolute_lower_bound"] > 0
        and dense_ok
    )

    report = {
        "outcome": "WITNESS" if accepted else "REJECTED",
        "independence": (
            "Replay reads only the compact graph spec and imports no "
            "construction/search module."
        ),
        "spec_path": str(args.spec.resolve()),
        "spec_sha256": spec_sha,
        "graph": {
            "vertex_count": total,
            "edge_count": int(edge_count),
            "simple": True,
            "connected": connected,
            "unweighted": True,
            "reflection_automorphism_verified": reflection_ok,
            "handshake_identity_verified": handshake_ok,
            "adjacency_oracle_checks": oracle_count,
            "class_sizes": [int(value) for value in sizes],
            "directed_block_degrees": [
                [int(value) for value in row] for row in degrees
            ],
            "class_degrees": [int(value) for value in class_degrees],
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
            "order_parameter": mp.nstr(order_parameter(sizes, angles), 60),
            "nonsynchronous": bool(order_parameter(sizes, angles) < mp.mpf("0.99")),
        },
        "spectrum": {
            "quotient": strip_private_fields(quotient),
            "quotient_interval_certificate": quotient_interval,
            "transverse_fourier": fourier,
            "transverse_operator_norm_bound": operator_bound,
            "absolute_hessian_gap": mp.nstr(min(quotient_gap, mp.mpf(transverse_gap)), 60),
            "normalized_hessian_gap": mp.nstr(
                min(quotient_gap, mp.mpf(transverse_gap)) / total, 60
            ),
            "full_gap_float": full_gap,
            "exactly_one_zero_eigenvalue": exactly_one_zero,
        },
        "dense_analog_spot_check": dense,
        "acceptance": {
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
            "all_fourier_transverse_modes_positive": transverse_gap > 0,
            "rigorous_transverse_operator_bound_positive": (
                operator_bound["rigorous_absolute_lower_bound"] > 0
            ),
            "exact_edge_count_matches_degree_handshake": handshake_ok,
            "exactly_one_zero": exactly_one_zero,
            "dense_analog_decomposition_matches": dense_ok,
            "accepted": accepted,
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(
        json.dumps(
            {
                "outcome": report["outcome"],
                "vertex_count": total,
                "mu_exact": report["graph"]["mu_exact"],
                "normalized_gap": report["spectrum"]["normalized_hessian_gap"],
                "torque_residual": equilibrium[
                    "maximum_per_vertex_torque_residual"
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
