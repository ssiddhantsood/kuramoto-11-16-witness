# Independent audit

Verdict: **CONFIRMED**

The audit reconstructed the graph, equilibrium, quotient Hessian, complete
unequal-fiber transverse decomposition, and nonlinear dynamics from
`graph_spec.json` without importing the original search, construction,
spectrum, replay, or simulation code.

## Exact graph checks

- Input specification SHA-256:
  `03ed576a456da292bb86aa9b6c01e3d0fb531e600d3650671b2812593be74c38`
- Class sizes:
  `[80325, 80325, 108675, 108675, 41400, 41400]`
- Class degrees:
  `[316835, 316835, 326096, 326096, 316802, 316802]`
- Edge count: `74003856975`
- Exact connectivity:
  \[
  \mu=\frac{316802}{460799}
  =\frac{11}{16}+\frac{43}{16\cdot460799}.
  \]

All 460,800 loop checks, 60,000 random rule/symmetry queries, and 18 complete
row enumerations passed. Clique fibers and the connected positive-edge class
graph prove connectivity.

## Equilibrium

A fresh gauge-fixed Newton solve at 160 decimal digits converged to the same
phase orbit as the stored decimals within \(4.08\times10^{-90}\) radians.

- Stored-decimal torque residual: \(1.29094\times10^{-85}\)
- Refined torque residual: \(1.60\times10^{-156}\)
- Order parameter: \(0.0321026318990039571\)

The original certificate's \(1.41\times10^{-89}\) residual refers to an
unprinted intermediate refinement. This is a non-material reporting
discrepancy.

## Quotient spectrum

In the orthonormal class-constant basis:

```text
0
48.7439639682460131553749244344726
73.7880140018670535021889912309646
3794.02570606130988371980595924
118081.610793737812860654843733
131873.477910801260229842267359
```

The normalized gap is
`0.0001057811718060894382712129436511992094`.
Outward-rounded interval entries plus Weyl bounds keep the gap strictly
positive. Exact edge balance gives the single rotational zero.

## Complete transverse accounting

For unequal classes, compatible Fourier characters split as:

- `4724 * 4 = 18896` modes shared by classes 0–3;
- `75600 * 2 = 151200` class-0/1 modes;
- `103950 * 2 = 207900` class-2/3 modes;
- `41399 * 2 = 82798` class-4/5 modes.

Therefore

```text
18896 + 151200 + 207900 + 82798 = 460794
6 + 460794 = 460800.
```

Every block was enumerated. The minimum transverse eigenvalue is
`111089.563790368676779...`, at character 4724 (conjugate character 1).
Independent operator bounds are also strictly positive.

## Nonlinear checks

Both signs of weakest-quotient perturbations of RMS \(10^{-4}\) and
\(10^{-3}\) returned with decreasing energy. Independent full-size random
perturbations also decayed toward the equilibrium orbit.

Coherent weak-mode perturbations at RMS \(10^{-2}\) can instead converge to
synchrony, demonstrating a finite basin. The claim is local asymptotic
stability modulo rotation, not global attraction.
