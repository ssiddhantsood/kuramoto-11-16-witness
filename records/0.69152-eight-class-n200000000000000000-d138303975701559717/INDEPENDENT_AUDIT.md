# Adversarial independent audit: 0.6915198785 record

## Verdict

**CONFIRMED**, with scope notes on the nonlinear basin and the producer's
certificate packaging.

The central near-threshold assertion survives an independent 768-bit Arb
calculation:

```text
normalized quotient gap in
[0.0000010000000156152773495392052108047429893159771133581723801621732323,
 0.0000010000000156152773495392052108047429893159771133581723801621732324]
```

Thus the gap is rigorously greater than \(10^{-6}\). The normalized margin is
strictly greater than

```text
0.000000000000015615277349539205210804742989315977113358172380162173232392690
```

Equivalently, the exact absolute shift is
\(N/10^6=200000000000\), while the certified eigenvalue is approximately
`200000003123.055469907841...`; its absolute excess is greater than
`3123.0554699078410421609485`.

The audit used only the supplied graph JSON through completion of its own exact
graph reconstruction, high-precision equilibrium solve, Krawczyk proof, and
spectral calculations. Approach 156 certificates and reports were opened only
after those independent results had been written.

## 1. Exact graph reconstruction

I interpreted a residue block `(i,j,m,a,s)` in the standard consecutive cyclic
way: vertices have local integer labels, and `(i,x)` is adjacent to `(j,y)` iff
`((y-x-a) mod m) < s`. This is also the interpretation later found in the
producer verifier.

The eight sizes sum exactly to `200000000000000000`. All 28 unordered class
pairs occur exactly once. The pair partition is:

- residue blocks: `01, 04, 26, 37, 45`;
- complete blocks: `02, 03, 06, 07, 12, 13, 15, 23, 24, 34, 46, 47, 56, 57, 67`;
- absent blocks: `05, 14, 16, 17, 25, 27, 35, 36`.

Every modulus divides both incident class sizes, every start is in range, and
every count is at most its shift capacity (the modulus). The nontrivial block
data are:

| Pair | Modulus / shift capacity | Forward start | Count | Reverse start |
|---|---:|---:|---:|---:|
| 0–1 | 2000000000 | 0 | 1736289442 | 263710559 |
| 0–4 | 38303976000000000 | 0 | 30195442942996172 | 8108533057003829 |
| 2–6 | 7457675000000000 | 0 | 6829047 | 7457674993170954 |
| 3–7 | 6893523000000000 | 0 | 1637747 | 6893522998362254 |
| 4–5 | 2000000000 | 0 | 1736289442 | 263710559 |

For each of all 28 blocks, `reports/graph_exact.json` records the pair capacity
`n_i*n_j`, shift capacity, starts, count, both directed degrees, and exact edge
count. For a residue block,

```text
k[i,j] = (n_j/m)*s,   k[j,i] = (n_i/m)*s,
|E_ij| = n_i*k[i,j] = n_j*k[j,i].
```

The independently reconstructed directed cross-degree rows are:

```text
0: [0,41102160758563546,7457675000000000,6893523000000000,30195442942996172,0,7457675000000000,6893523000000000]
1: [33253394557710696,0,7457675000000000,6893523000000000,0,47344826000000000,0,0]
2: [38303976000000000,47344826000000000,0,6893523000000000,38303976000000000,0,6829047,0]
3: [38303976000000000,47344826000000000,7457675000000000,0,38303976000000000,0,0,1637747]
4: [30195442942996172,0,7457675000000000,6893523000000000,0,41102160758563546,7457675000000000,6893523000000000]
5: [0,47344826000000000,0,0,33253394557710696,0,7457675000000000,6893523000000000]
6: [38303976000000000,0,6829047,0,38303976000000000,47344826000000000,0,6893523000000000]
7: [38303976000000000,0,0,1637747,38303976000000000,47344826000000000,7457675000000000,0]
```

Adding the clique degree `n_i-1` gives:

| Class | Size | Total degree |
|---:|---:|---:|
| 0 | 38303976000000000 | 138303975701559717 |
| 1 | 47344826000000000 | 142294244557710695 |
| 2 | 7457675000000000 | 138303976006829046 |
| 3 | 6893523000000000 | 138303976001637746 |
| 4 | 38303976000000000 | 138303975701559717 |
| 5 | 47344826000000000 | 142294244557710695 |
| 6 | 7457675000000000 | 138303976006829046 |
| 7 | 6893523000000000 | 138303976001637746 |

The minimum is attained exactly on classes 0 and 4. The exact edge count is

```text
14019316159188853020996070000000000
```

and the degree handshake sum is exactly twice this value. The positive
class-support graph is connected. Cliques have no loops; cross blocks join
different classes; the pair partition and residue predicate give no duplicate
edges. The graph is therefore finite, simple, undirected, and exactly
biregular on every cross block.

An explicit label-level reflection is

```text
R(c,x) = (c+4 mod 8, -x+b_c mod n_c)
b = (0,1736289441,0,0,0,1736289441,0,0).
```

It is involutory. On `01 ↔ 45`, it sends a shift `t` to
`1736289441-t`, reversing the allowed interval. On self-paired residue blocks
`04`, `26`, and `37`, the oriented shift is preserved after the endpoints are
reordered. All complete and absent blocks map to blocks of the same type.

Twenty deterministic modular boundary samples checked the first and last
included shift and the adjacent excluded shifts in both orientations, without
materializing the graph. Their SHA-256 is
`b2589f7808b8d3b7927d9da398861756fe8909c70a33a624accc9cc622a9426e`.

The exact density is

```text
mu = 138303975701559717/199999999999999999
   = 0.6915198785077986...
```

and

```text
mu - 11/16 = 12863611224955483/3199999999999999984.
```

Exact improvements over the prior records are:

```text
over approach 154 (55018/80001):
  60856360100478974735/16000199999999999919999
  = 0.003803474962842900404...

over approach 155 (5521403/7989063):
  638575030229782896574/1597812599999999992010937
  = 0.000399655773292677063...
```

Both cross-product excesses are positive.

## 2. Independent equilibrium and interval uniqueness

Starting only from the four short phase seeds, a 260-decimal Newton solve
converged in five iterations. With
`theta=(x0,x1,x2,x3,-x0,-x1,-x2,-x3)`, the independently obtained centers are:

```text
x0 = 0.729558778480953506520991965332611684009029146181460523415641697663019221438076012231668377942167531160576841162017604603964712583051762265130453562439015256093713151812500811553610754126364056390853483866920038
x1 = 2.358746108012727693153241587473315254157358613178879374913254676050679400980303954682875603454978206993501349017414467192263054558040411684242948592946514431993863722308863902659855299661289343127759869368811224
x2 = 0.956640994361590213471092208356671377489813449353931084846813920488998079579715898132895649887740007759741679516627971824882792216043176540913924942117092228230893642282501893564616724214369065456074738011883871
x3 = 0.956640994450157511720289930946293217057755914277296407851752110609490285985418882638270155661950401476556053977533578012934671753346587966460737617158548832032435262624787737191900411057751030823870331009072954
```

The explicit rational box is `x_i = center_i ± 10^-190`. A genuine
python-flint/Arb Krawczyk calculation at 768 bits used an exact rational,
nonsingular preconditioner. Its image offsets are only about `4e-211`, strictly
inside the rational box. Hence there is exactly one reflection-antisymmetric
root in this box. The eight scaled torque intervals all contain zero, with
maximum radius below `5.8e-191`.

The mass-orthonormal Hessian has no kernel beyond rotation (Sections 3–4), so
the equilibrium is also isolated modulo rotation in the full phase space.
This is local uniqueness; no global uniqueness among all equilibria is claimed.

The order parameter is rigorously enclosed by

```text
[0.0326139110233434853562408862029204230408922888440323093903249675290108317428
 +/- 1.51e-77]
```

with zero imaginary part exactly. The equilibrium is nonsynchronous.

## 3. Quotient spectrum and the near-threshold gap

For class-constant vectors, the independently derived normalized
mass-orthonormal quotient is

```text
Q_ii = (1/N) sum_j k_ij cos(theta_j-theta_i)
Q_ij = -(1/N) sqrt(k_ij*k_ji) cos(theta_j-theta_i).
```

The complete 8×8 interval matrix is in
`reports/quotient_spectrum_arb.json`. Reflection gives the even and odd
4×4 blocks `E=A+C` and `O=A-C`. Their certified normalized spectra are:

```text
even: 0 (exact),
      0.00000119776475066004021518294326042786639...,
      0.25629544709348496189225755384144709336...,
      0.27605019125679198064035159772445913066...

odd:  0.00000100000001561527734953920521080474299...,
      0.00515803939882596627944760044911132922...,
      0.27605019124291765511651131889681572014...,
      0.28720555849923163896452460801635730655....
```

The strict threshold was not inferred from rounded eigenvalues. At 768 bits:

1. `O - 10^-6 I` has all four outward-rounded leading principal minors
   strictly positive.
2. `E - 10^-6 I`, restricted by an interval basis to the exact
   mass-orthogonal complement of rotation, has all three leading principal
   minors strictly positive.
3. `Q*sqrt(n)=0` holds algebraically term by term from
   `n_i*k_ij=n_j*k_ji`.

Therefore the normalized quotient has one exact rotation zero and seven
strictly positive eigenvalues, and its smallest positive eigenvalue is
strictly above `1e-6`.

Rounding sensitivity was checked adversarially. Evaluating at the 30-digit
short seed changes the gap by only about `-1.62e-30`. Converting the large
directed degrees to binary64 changes some integers by `-2`; a binary64
eigensolve still reports a positive margin near `1.5613e-14`, but that
calculation cannot establish the strict inequality. The Arb and shifted
Sylvester tests are the certification.

## 4. All transverse dimensions

Let `W_i` be the zero-sum subspace inside class `i`. Its dimension is
`n_i-1`, and the full transverse space is

```text
W = direct_sum_i W_i,
dim(W) = sum_i(n_i-1) = N-8 = 199999999999999992.
```

The proof does not enumerate Fourier modes. On `W_i`, a clique contributes
`n_i I`, and a complete cross block vanishes between class-zero-sum spaces.
For a residue block with modulus `m`, fiber multiplicities
`q_i=n_i/m`, `q_j=n_j/m`, and `s` present shifts, aggregate each residue
fiber. On the nonconstant residue space,

```text
C_S = -C_(S complement),
||C_S|| <= min(s,m-s),
||A_ij restricted|| <= sqrt(q_i*q_j)*min(s,m-s).
```

This is a complement-aware regular-operator bound and uses no character
enumeration. Applying Cauchy–Schwarz to the five residue blocks reduces the
full transverse quadratic form to an outward-rounded 8×8 comparison matrix.
The smallest normalized comparison eigenvalue is rigorously in

```text
[0.24007629108526466650363936535515075386137367423869322394160335639,
 0.24007629108526466650363936535515075386137367423869322394160335640].
```

All eight leading principal minors after the exact shift `0.24 I` are
strictly positive. Thus all `N-8` transverse dimensions are positive with
normalized lower bound strictly greater than `0.24`. Sample character
formula evaluations were made only as spot checks; all lie below the
Fourier-free bound and are not used in the proof.

Combining quotient and transverse spaces gives:

```text
positive dimensions: 7 + (N-8) = N-1 = 199999999999999999
zero dimensions:     1 (global rotation)
```

## 5. Nonlinear return and basin limitations

The exact equitable quotient was integrated in normalized time
`tau=N*t` along the weakest quotient eigenvector. Independent tests returned
for signed maximum phase amplitudes

```text
-1e-6, +1e-6, -1e-5, +1e-5, -5e-5, +1e-4, +2e-4.
```

At `tau=5e6`, the smallest tests contracted by factors close to the linear
prediction `exp(-gap*tau) ≈ 0.00673795`.

The basin is narrow and asymmetric. Tests at `-1e-4` and `+2.5e-4` escaped
the local basin and converged numerically to other equilibria. Those boundary
locations are numerical, not interval-certified.

There is also an analytical local statement. If all quotient class phases
remain within sup-radius `r`, the normalized Hessian perturbation is at most
`4*r*d_max/N`. Choosing

```text
r = 1.7569228093580840579633828257e-7
```

leaves at least half the certified quotient gap. Monotonicity of the weighted
distance then gives a conservative guaranteed weak-mode amplitude below
`3.7445493608132229181e-8`. Along a class-constant quotient trajectory the
transverse space remains invariant. While every class phase obeys the same
tube condition,

```text
lambda_trans/N >= 0.240076291085264... - 4*r*d_max/N.
```

This stays positive for `r < 0.0843591010414677...`; for example, at
`r=0.05` the analytical transverse lower bound remains about `0.0977820`.
For arbitrary non-equitable perturbations, this sector invariance is absent;
the statement is then only a local Hessian-curvature bound, not a decoupled
transverse simulation.

Accordingly, the nonlinear verdict is **CONFIRMED_LOCAL_ONLY**. No global,
macroscopic, or full-graph numerical basin claim is made.

## 6. Post-independence comparison

After the independent outputs were fixed:

- all exact graph fields, degree data, edge count, `mu`, dimension split, and
  zero multiplicity matched approach 156;
- the producer gap and transverse decimals agree with the independent
  intervals to every digit the producer reports;
- the producer's stored phases agree with the independent root to every stored
  digit (they are truncated at about 160 decimal places and therefore are not,
  as exact decimals, inside the much tighter `1e-190` audit box);
- every cited graph-spec SHA-256 for the highest record matches
  `caa32a330b61ab4907f77a11773d84b7c74087f91d52c844d426235c542bf774`.

Three scope/package discrepancies were found:

1. The producer 512-bit formal certificate proves the quotient threshold
   `>1e-6` but only proves the transverse bound `>1/5`; its `>0.240` value was
   numerical. This audit independently supplies a 768-bit formal
   shifted-`0.24` proof.
2. The per-record numerical certificate contains
   `nonlinear_local_return: null` while its acceptance block says
   `nonlinear_return: true`. The separate highest-record certificate does
   contain eight small-amplitude return tests.
3. Return is genuinely local; larger weak-direction amplitudes can escape.

None changes the record verdict.

As a secondary check, exact graph/mu reconstruction passed for all six
approach 156 ledger records, all six approach 155 records, and the approach
154 record. Referenced approach 156 certificate hashes and approach 154/155
artifact hashes are recorded in `reports/ledger_spotchecks.json`.

## Reproduction and artifacts

Run from this directory:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python independent_verify.py
.venv/bin/python nonlinear_return.py
.venv/bin/python compare_artifacts.py
.venv/bin/python ledger_spotcheck.py
```

Primary outputs:

- `verdict.json`
- `reports/independent_results.json`
- `reports/graph_exact.json`
- `reports/root_arb.json`
- `reports/quotient_spectrum_arb.json`
- `reports/transverse_arb.json`
- `reports/nonlinear_return.json`
- `reports/artifact_comparison.json`
- `reports/ledger_spotchecks.json`

The full graph was never materialized.
