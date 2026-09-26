# Dynamics-to-Hessian audit

**Verdict: PASS**

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

- Graph specification SHA-256: `5dca79605da2e0eeff844f39f5d1901ed77778c105d6fba1f92687e82f54bf9f`
- Vertices: `80002`
- Six class sizes: `[14133, 14133, 7024, 7024, 18844, 18844]`
- Other classes connected to each class: `[4, 4, 4, 4, 3, 3]`
- Newton iterations from the short stored seed: `0`
- Point maximum torque residual: `3.63797880709e-12`
- Order parameter: `0.0332994842555` (not synchronous)

Biregularity is checked exactly with integer edge counts. It is what permits
the six class equations to represent every one of the `80002`
oscillator equations at a class-constant phase assignment.

The point residual is only a numerical sanity check. Existence of an exact
equilibrium is supplied by an outward-rounded Krawczyk inclusion over a
radius-`1.00000000000e-10` box. Its interval
contraction bound is `0.0000007440931187659487059725447`, and
the Krawczyk image lies strictly inside the box.

Nonsynchrony is also certified throughout the box: the phase separation between
classes 0 and 1 lies in
`[1.41500458424797120000000000000,
1.41500458464797120000000000000]`, strictly
between `0` and `pi`.

## 3. Direct Jacobian-versus-Hessian test

The code first differentiates the actual sine flow in ordinary class phase
coordinates. It then changes to the orthonormal coordinates
`q_i=sqrt(n_i) delta-theta_i`. In this basis the matrix should be symmetric and
equal to the negative quotient Hessian.

- Maximum entry error in `J_mass + H`: `0`
- Separate central finite-difference relative error: `1.10380764291e-10`
- Hessian symmetry error: `0`
- Rotation-null residual: `1.45519152284e-10`

The finite-difference test is deliberately separate from the analytic matrix
construction: it perturbs the phases, reevaluates the sine flow, and compares
the observed derivative with `-H`. It is a sanity check, not the rigorous
spectral certificate.

## 4. Quotient and transverse modes

| Quotient mode | Hessian eigenvalue | Eigenvalue after `/N` normalization |
| ---: | ---: | ---: |
| 0 | -4.09953381684e-13 | -5.12428916383e-18 |
| 1 | 0.426830958998 | 5.33525360613e-06 |
| 2 | 23.2776753613 | 0.000290963667924 |
| 3 | 621.332457962 | 0.00776646156299 |
| 4 | 20522.8961451 | 0.256529788569 |
| 5 | 22942.7299176 | 0.286776954546 |

The near-zero mode is global rotation. The other five quotient eigenvalues are
positive, so the corresponding flow eigenvalues are negative.

For rigorous rounding control, the verifier evaluates an interval Hessian over
the *entire certified root box*. It adds the exact rotation projector and proves

```text
H + rotation_projector - 0.4 I  is positive definite
```

by interval `LDL^T`. Every interval pivot has a positive lower endpoint. Thus
every non-rotation quotient eigenvalue at the exact enclosed equilibrium is
strictly greater than `0.40000000000000000000`
(or `4.9998750031249218770e-6` after
division by `N`).

The quotient alone is not enough: it covers only six directions. The remaining
`79996` fiber-zero-sum directions are bounded using:

1. `n_i I` from the clique inside fiber `i`;
2. zero action from complete cross-blocks on zero-sum vectors; and
3. `||A_ij|| <= sqrt(d_ij d_ji)` for every partial biregular block.

The conservative point lower bound for every transverse Hessian mode is
`17725.8475415`. The outward-rounded bound
over the full equilibrium box is `17725.8475295`,
which is positive.

## 5. Machine checks

| Check | Result |
| --- | --- |
| all six torques near zero | PASS |
| rigorous verifier ran fresh | PASS |
| rigorous report matches graph spec | PASS |
| exact equilibrium exists in krawczyk box | PASS |
| hessian is symmetric | PASS |
| jacobian equals negative hessian | PASS |
| finite difference confirms derivative | PASS |
| rotation is only quotient zero | PASS |
| interval certificate proves quotient gap | PASS |
| all transverse modes have positive bound | PASS |
| outward rounded transverse bound is positive | PASS |
| equilibrium is nonsynchronous | PASS |

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
