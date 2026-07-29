# Independent adversarial audit of the 0.691120 record

## Verdict

**CONFIRMED.**

The compact specification defines a finite simple connected unweighted graph,
and the independently solved equilibrium has a rigorously positive Hessian on
every direction except global rotation. The small quotient margin survives an
outward-rounded interval certificate:

- exact density
  \[
  \mu=\frac{27,607,015}{39,945,315}
     =\frac{5,521,403}{7,989,063}
     =0.691120222734505911394114679031\ldots;
  \]
- exact excess over \(11/16\):
  \[
  16(27,607,015)-11(39,945,315)=2,313,775,
  \]
  equivalently
  \(\mu-11/16=462755/127825008\);
- exact improvement over the approach-154 record \(55018/80001\):
  \[
  \frac{5,521,403}{7,989,063}-\frac{55,018}{80,001}
  =\frac{725,164,423}{213,044,343,021}
  =0.0034038191895502233402\ldots.
  \]

No material discrepancy was found after comparing the prior artifacts.

## Independence protocol

Before completing the calculations below, the audit read only
`approach_155_six_class_density_frontier/records/mu_27607015_39945315/graph_spec.json`.
It did not import approach-155 construction, search, replay, or report code.
Only after the graph, root, quotient, transverse, operator, Fourier diagnostic,
and nonlinear calculations had completed were the prior reports and the
approach-154 baseline inspected.

The interval calculations use python-flint Arb at 180–220 decimal digits. The
root solve uses a separate 240-decimal Newton implementation initialized only
from the three short seeds in the graph specification.

## 1. Exact compact graph

The class sizes are

```text
[7,652,907, 7,652,907, 2,866,160, 2,866,160, 9,453,591, 9,453,591]
```

and sum exactly to \(N=39,945,316\). All 15 unordered cross-class pairs occur
exactly once. Every residue modulus equals the corresponding gcd, every shift
interval has valid start/count, and exact biregular handshakes hold.

The independently reconstructed directed degree matrix, including clique
degrees on the diagonal, is

```text
7,652,906  6,020,680  2,866,160  2,866,160  8,201,109          0
6,020,680  7,652,906  2,866,160  2,866,160          0  8,201,109
7,652,907  7,652,907  2,866,159          0          0  9,453,591
7,652,907  7,652,907          0  2,866,159  9,453,591          0
6,638,993          0          0  2,866,160  9,453,590  9,453,591
        0  6,638,993  2,866,160          0  9,453,591  9,453,590
```

Thus the class degrees are

```text
[27,607,015, 27,607,015, 27,625,564,
 27,625,564, 28,412,334, 28,412,334].
```

The exact handshake gives \(559,051,789,848,239\) undirected edges. The
six-class support is connected, and each class is a nontrivial clique, so the
full graph is connected.

Reflection was checked with the class involution
`[1,0,3,2,5,4]` and local map \(u\mapsto-u\). In particular, reversing a stored
residue interval changes its start to

```text
(-shift_start - shift_count + 1) mod modulus.
```

All 8,553,249 residues across the three fractional blocks were exhaustively
checked. The reflected interval for block 0–4 has start 0 because
`59643 + 390529 - 1 = 450171`; it matches block 1–5 exactly. There were no
mismatches in 200,000 deterministic random adjacency, symmetry, or reflection
queries. Six representative rows were also enumerated against all \(N\)
targets; their counts and packed-bit hashes are in `graph_audit.json`.

## 2. Equilibrium and interval root

Starting only from

```text
(-0.7298280527779402, 0.9562664566002271, -2.3588991933226247)
```

the independent solve converged to

```text
theta1 = -0.7298280527765897229697129801179898060357449160706467174414658360662728516421171864560293797789015254...
theta2 =  0.9562664565985826186444191292220588927360317239017627527188665517775061686058749353996018110354844152...
theta3 = -2.3588991933212602771956656878482580781737733343514834633984263092394449420181618577115259312482157758...
```

The six phases are
`[theta1,-theta1,theta2,-theta2,theta3,-theta3]`. All three phase-reflection
equations and all three paired torque-reflection equations vanish identically.
The six independently evaluated per-vertex torque residuals have maximum
absolute value

```text
1.25803686906194009926309306706e-234.
```

The order parameter is

```text
0.0327161548791958039801716763194271908900980497078415655751404,
```

with zero imaginary part by reflection, so the root is nonsynchronous.

For each coordinate, the certified local box is the displayed center plus or
minus `1e-30` (with Arb also enclosing decimal-conversion roundoff). A
220-decimal outward-rounded Krawczyk computation has strict interior inclusion
in this box. Its infinity-norm contraction bound is

```text
4.005572036354451769158588502670102855586e-25.
```

Therefore a root exists and is unique in the stated box. Full outward box and
Krawczyk-image endpoints are in `equilibrium_certificate.json`.

## 3. Quotient Hessian and the close acceptance margin

The full symmetric 6×6 quotient was reconstructed on the orthonormal
class-indicator basis and independently split into 3×3 reflection-even and
reflection-odd blocks. The full numerical spectrum is

```text
0 (structural rotation)
45.6585929825134962401175155113207645956323725
215.514274299024421706616201671772350364549362
205,538.577711938996723267265931050426576964243
10,238,016.716570021137031503084082492648645477
11,470,110.4241268750242756228246965772085339456
```

The even block has the exact null vector

```text
(sqrt(2*7,652,907), sqrt(2*2,866,160), sqrt(2*9,453,591)).
```

Its two nonzero eigenvalues were enclosed from the exact rank-two trace
polynomial, avoiding any numerical decision about the structural zero. The odd
block was certified by outward interval eigenvalue bounds and shifted
Sylvester tests. The weakest eigenvalue lies in

```text
[45.65859298251349624011717980759663310298848583997537671,
 45.65859298251349624011785121504407673463041035405039999].
```

Normalizing by \(N=39,945,316\) gives

```text
[1.143027457399848739214309377539950694169e-6,
 1.143027457399848739214326185704578647835e-6].
```

The rigorous lower margin above `1e-6` is therefore

```text
1.4302745739984873921430937753995069e-7,
```

or 14.3027% of the acceptance floor. Equivalently, the absolute lower bound
exceeds the required `39.945316` by more than `5.7132769825`. This directly
addresses the close-margin risk.

Differentiating the phase dynamics independently gives \(J=-H\) edge by edge.
On class-constant coordinates,
`diag(sqrt(n_i)) J_class diag(1/sqrt(n_i)) = -Q`; the high-precision entry
check is below `6e-95`.

## 4. All 39,945,310 transverse dimensions

Let each class component have zero sum. A clique adjacency then acts as
\(-I\), while every complete cross block vanishes exactly. For a fractional
biregular block \(A_{ij}\), its complement satisfies \(A_{ij}=-A^c_{ij}\) on
this subspace, hence

```text
||A_ij|| <= min(
    sqrt(d_ij*d_ji),
    sqrt((n_j-d_ij)*(n_i-d_ji))
).
```

This is Fourier-free and applies to every transverse vector simultaneously.
Using outward phase intervals produces a 6×6 scalar comparison matrix. Its
reflection-even and reflection-odd blocks both pass shifted Sylvester tests at
the integer shift

```text
9,589,733.
```

Thus the full transverse Hessian is rigorously greater than `9,589,733 I`.
This is stronger than, and therefore confirms, the claimed prior lower bound
`9,532,144.352250697`.

As a non-rigorous error-detection cross-check, exact character accounting gives

```text
14,405,472  A-only dimensions
 1,800,680  shared A/C dimensions
18,006,840  C-only dimensions
 5,732,318  B dimensions
-----------
39,945,310  transverse dimensions.
```

Among 15,968 deterministic diagnostic characters, the smallest sampled
eigenvalue was

```text
9,591,967.7375068178927685793043216899
```

at \(q=\pm1/450171\), with active classes 0,1,4,5. Interval/complement
coefficient magnitudes agreed to at least 70 decimal places.

An independently assembled 620-vertex unequal-fiber dense analog was compared
against a separately coded residue-aggregation Hessian matvec. Across 24
random vectors, the maximum relative error was `8.58e-16`; 12 class-mean-zero
trials also preserved the complement to roundoff.

Combining quotient positivity, transverse positivity, and the exact rotation
null vector proves that the full Hessian has **exactly one zero eigenvalue**.
The full Jacobian has exactly one zero and all other eigenvalues are negative.

## 5. Nonlinear local return

The weakest quotient eigenvector was perturbed with both signs at amplitudes
`1e-8` and `1e-7`. All four radial flow projections were negative, all energy
changes were positive, and the measured restoring curvature approached
`45.6585929825...`. Another 160 deterministic random compact perturbations,
which assign phases to all \(N\) vertices classwise, also returned.

There is also a conservative analytic full-graph certificate. The spectral
work proves \(\lambda_2(H)>45\). If every vertex moves by at most \(r\), each
edge cosine changes by at most \(2r\), and the Hessian perturbation norm is at
most \(4d_{\max}r\). At \(r=10^{-7}\),

```text
4*d_max*r = 11.3649336,
```

so every rotation-orthogonal perturbation on all \(N\) vertices has strong
radial return coefficient greater than

```text
45 - 11.3649336 = 33.6350664.
```

Finite-basin caveat: this certifies only the stated local sup-norm ball. It
does not establish global attraction or characterize the full basin.

## 6. Post-audit comparison

After the independent calculations, the prior approach-155 replay was read.
Graph counts, degrees, edges, density, root, quotient gap, exactly-one-zero
conclusion, and nonlinear verdict agree.

Two non-material differences were found:

1. The graph spec's construction-time binary64 screening gap is
   `1.1430274655129156e-6`, while the resolved-root value is
   `1.1430274573998487e-6`, a difference of `8.1131e-15`. The prior interval
   replay and this audit agree on the latter.
2. The prior report stated the conservative Gershgorin transverse bound
   `9,532,144.352...`; this audit certifies the stronger comparison-matrix
   eigenvalue bound `>9,589,733`.

Neither affects acceptance.

The other five exact approach-155 specs were spot-checked independently for
class sums, gcd moduli, reflection interval orientation, directed degrees,
handshake/edge count, and exact reduced density. Their replay summaries agree
and report accepted witnesses with one zero and positive transverse bounds:

```text
22115608/32008211
14133442/20480013
641528/930169
688720/999961
880544/1279921
```

## Artifacts

- `verdict.json` — machine-readable final verdict.
- `graph_audit.json` — exact block, degree, reflection, row, and adjacency audit.
- `equilibrium_certificate.json` — 240-decimal solve and Arb/Krawczyk certificate.
- `spectral_certificate.json` — full/even/odd quotient and transverse certificates.
- `fourier_mode_crosscheck.json` — exact mode accounting and sampled characters.
- `operator_validation.json` — independent dense/compact matvec validation.
- `nonlinear_return.json` — weak-mode probes and rigorous local-return radius.
- `other_records_spotcheck.json` — checks of the five lower exact records.
- `artifact_comparison.json` — post-independence comparison and baseline improvement.
- `audit_common.py`, `audit_graph.py`, `certify_equilibrium.py`,
  `certify_spectrum.py`, `audit_fourier_modes.py`, `validate_operators.py`,
  `nonlinear_return.py`, `spotcheck_other_records.py`, and
  `compare_prior_artifacts.py` — independent audit implementations.
