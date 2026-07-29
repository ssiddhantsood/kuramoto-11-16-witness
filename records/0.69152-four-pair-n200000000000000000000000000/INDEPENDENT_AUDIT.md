# Adversarial independent audit of the 0.6915227109103000046 record

## Verdict

**CONFIRMED.**

The declarative graph specification defines a finite, connected, simple,
undirected, unweighted graph with the claimed reflection automorphism. An
independent 280-decimal-digit solve and 1024-bit Arb audit certify a unique
reflected equilibrium in the stated local box, exactly one Hessian zero from
global rotation, a normalized quotient gap strictly above \(10^{-6}\), and
positivity of all \(N-8\) transverse modes. Independent class-level nonlinear
tests return from both signs of the weak mode at five amplitudes.

The critical lower bound is

\[
\lambda_{\rm quotient}>
\frac{1000047662728726701}{10^{24}}
=0.000001000047662728726701,
\]

so the rigorously certified excess over the acceptance floor is

\[
\frac{47662728726701}{10^{24}}
=4.7662728726701\times 10^{-11}.
\]

That margin is genuinely small, but it is positive by a wide interval-arithmetic
margin relative to the 1024-bit rounding radii.

## Independence protocol

The core calculation read only:

`approach_158_reflection_pair_expansion/artifacts/finite/m4_gamma_1e-06_lift_buffer_1p0005e-6_N200000000000000000000000000_candidate0074_ceil_graph_spec.json`

No approach 158 code, report, certificate, or optimizer output was imported.
The core result was written before prior/producer artifacts were opened. The
post-audit comparison is a separate script.

- Graph-spec SHA-256:
  `50e073d06ec0384f455c9c75e24adeaf7974b0db296b9d1b40a5728f1fbc8c2a`
- Core arithmetic: mpmath at 280 decimal digits.
- Rigorous arithmetic: python-flint Arb at 1024 bits.
- Krawczyk radius: \(10^{-230}\) nominal; printed enclosing radius at most
  \(1.01\times10^{-230}\).
- Producer/prior code imported: none.

## 1. Exact graph reconstruction

### Classes and total size

The class sizes reconstructed as integers are:

- \(n_0=n_1=47344974154690000000000000\)
- \(n_2=n_3=14350483663250000000000000\)
- \(n_4=n_5=38300542182060000000000000\)
- \(n_6=n_7=4000000000000000000000\)

Their exact sum is

\[
N=200000000000000000000000000.
\]

Every class is a loop-free clique.

### All 28 cross blocks

For a block \(i<j\), write `R(g,k,s)` for the residue rule with modulus
\(g=\gcd(n_i,n_j)\), count \(k\), and start \(s\). Write `C(g)` for a
complete block (\(k=g\)) and `A(g)` for an absent block (\(k=0\)).

```text
(0,1) R(47344974154690000000000000,47344974154689987773202432,0)
(0,2) C(10000000000000)
(0,3) A(10000000000000)
(0,4) R(10000000000000,8681551781677,1318448218324)
(0,5) A(10000000000000)
(0,6) R(10000000000000,8516949755081,1483050244920)
(0,7) A(10000000000000)

(1,2) A(10000000000000)
(1,3) C(10000000000000)
(1,4) A(10000000000000)
(1,5) R(10000000000000,8681551781677,0)
(1,6) A(10000000000000)
(1,7) R(10000000000000,8516949755081,0)

(2,3) R(14350483663250000000000000,915026159,0)
(2,4) C(10000000000000)
(2,5) C(10000000000000)
(2,6) C(250000000000000)
(2,7) C(250000000000000)

(3,4) C(10000000000000)
(3,5) C(10000000000000)
(3,6) C(250000000000000)
(3,7) C(250000000000000)

(4,5) R(38300542182060000000000000,30193013201710114515451904,0)
(4,6) C(20000000000000)
(4,7) R(20000000000000,16174995918920,3825004081081)
(5,6) R(20000000000000,16174995918920,0)
(5,7) C(20000000000000)
(6,7) R(4000000000000000000000,396032970297458880,0)
```

For a residue block, a vertex in class \(i\) has exactly
\(k n_j/g\) neighbors in class \(j\), while a vertex in \(j\) has exactly
\(k n_i/g\) neighbors in \(i\). Thus

\[
n_i(k n_j/g)=n_j(k n_i/g),
\]

which proves exact biregularity and gives the cross-edge count without
enumerating vertices. Every count lies in \([0,g]\), so no multiple edge can
arise. The 28 blocks partition all unordered class pairs, and the relation is
used as an undirected adjacency relation.

The machine-readable core artifact records, for every block, both endpoint
degrees, its exact edge count, and the biregular product identity.

### Degrees, handshake, support, and reflection

The exact class degrees are:

```text
142294652774517828803813893
142294652774517828803813893
138304542182060000915026158
138304542182060000915026158
138304542182063947166013416
138304542182063947166013416
138304542182065487930135628
138304542182065487930135628
```

Hence

\[
\delta(G)=138304542182060000915026158.
\]

The exact edge and handshake counts are:

```text
|E|       = 14019365901080421922041506310039040918630000000000000
sum_v d_v = 28038731802160843844083012620078081837260000000000000
```

The handshake identity \(\sum_vd_v=2|E|\) holds exactly. The nonzero
class-support graph reaches all eight classes, and each class is a clique, so
the full graph is connected.

The reflection is

\[
(0\,1)(2\,3)(4\,5)(6\,7),\qquad u\mapsto -u\pmod {n_i}.
\]

For an interval \(S=\{s,\ldots,s+k-1\}\pmod g\):

- if reflected endpoint order is preserved, \(v-u\mapsto -(v-u)\), so the
  target start must be \(1-s-k\pmod g\);
- if canonical endpoint order reverses, the reordering supplies a second
  minus sign and the start remains \(s\).

All block orbits satisfy the applicable identity exactly. This proves
adjacency preservation for every local label \(u,v\), not just sampled
vertices. Clique adjacency and nonadjacency are also preserved.

### Exact density and prior improvements

The exact minimum-degree ratio is

\[
\mu=
\frac{138304542182060000915026158}
     {199999999999999999999999999}
=0.691522710910300004575130793457613554551500022875653967288067\ldots
\]

The cross-product test against \(11/16\) is

\[
16\delta(G)-11(N-1)
=12872674912960014640418539>0.
\]

After the independent audit completed, every graph spec in approaches 154,
155, and 156 was reconstructed exactly. The selected comparator is the maximum
within each family:

- Approach 154 maximum:
  \(55018/80001\). Exact improvement:
  `60901679106982133203007721176/16000199999999999999999999919999`
  \(=0.003806307365344316521231467199\ldots\).
- Approach 155 maximum:
  \(27607015/39945315=5521403/7989063\). Exact improvement:
  `643100678634817090201628431357/1597812599999999999999999992010937`
  \(=0.000402488175794093181016114427\ldots\).
- Approach 156 maximum:
  \(138303975701559717/199999999999999999\). Exact improvement:
  `113296100056644700689556243974786533559/39999999999999999799999999800000000000000001`
  \(=0.00000283240250141611753140091862\ldots\).

All 18 approach 158 graph specs were also reconstructed. The best other spec is
candidate 0073, whose minimum degree is exactly one smaller. Candidate 0074
therefore improves on every other producer spec by exactly

\[
\frac{1}{199999999999999999999999999}.
\]

## 2. Independent equilibrium and Krawczyk certificate

Starting only from the four short decimal seeds in the graph spec, Newton's
method converged in seven iterations at 280 decimal digits to:

```text
theta_1 = 0.7828438541688917951158131010732577071178649662104798882182382937046197380322793876104060183450045965...
theta_2 = 2.1849559284471583623000363349052103459802385925906175273311689023306393238888176282725056944327454234...
theta_3 = 2.4120224564927292000579521882853590401531322772914432981422497133351274978752084076108083165890391341...
theta_4 = 2.4420224577353025591124460648884278500283668601914747950343756856210073701524743353447143884115195097...
```

The full 270-digit centers and outward lower/upper endpoints are serialized in
`verdict.json` and `artifacts/independent_core_audit.json`. Each coordinate box
has nominal radius \(10^{-230}\), enclosed by an Arb radius no larger than
\(1.01\times10^{-230}\).

For the reflected phases
\((\theta_1,-\theta_1,\ldots,\theta_4,-\theta_4)\), the independently evaluated
maximum normalized torque residual at the 280-digit center is

\[
5.87183870453646174646614156108297236748\times10^{-282}.
\]

All eight point torques are listed in the core artifact. Evaluating all eight
torques over the complete root box gives absolute endpoints no larger than

\[
6.16477092517325835842653809427490164\times10^{-231}.
\]

The real signed order parameter is enclosed around

\[
-0.032610437885513306036624544539621538588703509601443691931359\ldots,
\]

so its magnitude is
\(0.0326104378855133060366245445396215385887\ldots\), and its imaginary
part is exactly zero by reflection.

### Genuine outward Krawczyk test

With exact binary midpoint \(x_0\), interval box \(X\), and a point
preconditioner \(C\), the audit evaluated

\[
K(X)=x_0-CF(x_0)+(I-CJ(X))(X-x_0)
\]

entirely with 1024-bit outward-rounded Arb operations.

- Arb proves \(\det C\) excludes zero.
- Every coordinate of \(K(X)\) is strictly inside the corresponding
  coordinate of \(X\).
- The largest displayed Krawczyk-image radius is
  \(3.54\times10^{-245}\), far below the \(10^{-230}\) box radius.

Therefore the box contains exactly one reflected torque root. This is a local
existence-and-uniqueness statement; no global uniqueness claim is made.

## 3. Quotient Hessian and the critical gap

Let \(r_{ij}\) be the degree from one class-\(i\) vertex into class \(j\), and
let \(E_{ij}=n_i r_{ij}=n_jr_{ji}\). In the orthonormal class basis
\(\mathbf1_i/\sqrt{n_i}\), the normalized quotient Hessian is

\[
\frac{Q_{ii}}N
=\frac1N\sum_{j\ne i}r_{ij}\cos(\phi_j-\phi_i),
\]

\[
\frac{Q_{ij}}N
=-\frac{E_{ij}}{N\sqrt{n_in_j}}\cos(\phi_j-\phi_i)
\quad(i\ne j).
\]

Within-class clique terms cancel exactly on class-constant vectors. Applying
the reflection-even basis
\((e_{2a}+e_{2a+1})/\sqrt2\) and reflection-odd basis
\((e_{2a}-e_{2a+1})/\sqrt2\) gives independent \(4\times4\) even and odd
blocks. The full \(8\times8\), even \(4\times4\), and odd \(4\times4\)
matrices are independently derived and printed in
`artifacts/independent_matrices.txt`.

The numerical normalized eigenvalues at the certified root are:

```text
0  (structural rotation; numerical residual about 2.7e-282)
0.000001000047662728726701958565362407914413768881923595119907...
0.000001000957177412321747046702900097749527296185388732473173...
0.005157559662474237448435672107690078044606691146414170200410...
0.256293734525421363061574795523621911955271137671050573183428...
0.258777170256315674599535764194312708031948086005821390387126...
0.258782822389926055897920456900876998675222826586194620349099...
0.287202775667496579028647306598044896011837272585206208050233...
```

For \(s_i=\sqrt{n_i}\), every cross pair cancels in \(Qs\) exactly:

\[
r_{ij}\cos(\Delta)s_i
-\frac{E_{ij}\cos(\Delta)}{\sqrt{n_in_j}}s_j=0.
\]

Thus global rotation gives one exact zero. In the even block, the audit adds
the rank-one lift \(2pp^\mathsf T\), where
\(p_a=\sqrt{2n_{2a}/N}\), and applies interval LDL to
\(Q_{\rm even}/N+2pp^\mathsf T-\lambda I\). It separately applies interval
LDL to \(Q_{\rm odd}/N-\lambda I\).

At

```text
lambda = 1000047662728726701 / 1000000000000000000000000
```

all outward interval pivots are strictly positive. The smallest lifted-even
pivot lower bound is approximately \(1.9965\times10^{-20}\), while its Arb
radius is approximately \(1.46\times10^{-140}\). This certifies the critical
strict inequality, not merely a floating-point estimate. Since the lifted
even and odd blocks are positive above this threshold and the unlifted even
block has the exact rotation zero, the quotient has exactly one zero.

## 4. All transverse modes

The full vertex space splits as the eight-dimensional class-constant quotient
plus classwise zero-sum spaces \(W_i\). Their exact dimensions are:

```text
47344974154689999999999999
47344974154689999999999999
14350483663249999999999999
14350483663249999999999999
38300542182059999999999999
38300542182059999999999999
3999999999999999999999
3999999999999999999999
```

The sum is exactly

\[
199999999999999999999999992=N-8.
\]

On \(W_i\), the clique contribution is \(n_iI\), so the normalized diagonal
term is

\[
\beta_i=\frac{n_i+\sum_{j\ne i}r_{ij}\cos(\phi_j-\phi_i)}N.
\]

Complete cross blocks annihilate classwise zero-sum vectors. For a residue
block with \(n_i=ga\), \(n_j=gb\), and a consecutive interval of \(k\)
residues, the nonconstant character \(q\ne0\) has singular value

\[
\sqrt{ab}\left|\sum_{h=0}^{k-1}
e^{2\pi iq(s+h)/g}\right|.
\]

The start changes only the complex phase. The audit used the rigorous
complement-aware bound

\[
\left|\sum_{h=0}^{k-1}e^{2\pi iqh/g}\right|
\le
\min\!\left(k,g-k,\frac1{\sin(\pi/g)}\right).
\]

This was evaluated outward with the certified phase boxes for all ten residue
blocks. The artifact also records \(q=1\) and \(q=2\) geometric-sum spot
checks, the \(g-1\) transverse characters, and the exact left/right null
dimensions for each block. No vertex or character enumeration was used.

Replacing each off-diagonal block by the negative of its certified
\(|\cos\Delta|\,\|A_{ij}\|/N\) upper bound gives an \(8\times8\) scalar
comparison matrix. A 1024-bit interval LDL factorization proves

\[
\frac{H|_{\oplus W_i}}N
>
\frac{9603117443111}{40000000000000}I
=0.240077936077775\,I.
\]

Therefore all \(N-8\) transverse modes are strictly positive.

## 5. Independent nonlinear weak-mode return

The audit integrated the exact eight-class restriction

\[
\dot\phi_i=\sum_j\frac{r_{ij}}N\sin(\phi_j-\phi_i)
\]

for normalized time \(8\times10^6\). The initial perturbation was the weakest
mass-orthonormal quotient eigenvector, converted to physical class phases.
Both signs were tested at weighted-RMS amplitudes:

```text
1e-7, 3e-7, 1e-6, 3e-6, 1e-5.
```

Every one of the ten runs had the restoring initial sign, decreased the
class-level energy, and returned to a rotation-aligned distance between
approximately \(3.35\times10^{-4}\) and \(3.42\times10^{-4}\) of its initial
distance. Thus the demonstrated finite symmetric basin along the tested weak
direction is at least weighted-RMS radius \(10^{-5}\).

This is an independent high-accuracy nonlinear numerical stress test. It is
not presented as an interval proof for every point in that basin.

## 6. Post-audit producer comparison and record spot checks

Only after the independent core artifact existed, the producer certificate was
opened. The following matched:

- graph-spec SHA-256;
- exact minimum degree, edge count, and \(\mu\);
- all four angles to more than 160 decimal places;
- quotient-gap estimate to more than 70 decimal places;
- refined root, quotient spectrum, and transverse comparison.

One metadata discrepancy is worth flagging. The graph spec's embedded
`construction_screen.quotient.normalized_gap` is
`1.0000471140778006e-6`, whereas the independently refined equilibrium gives
`1.000047662728726701958565...e-6`. The former is a preliminary
floating-point screen value at the short/unrefined phase data; it must not be
used as the final certificate. The refined producer certificate agrees with
the independent result.

The other five certified record specs were reconstructed exactly. For each,
the actual spec hash equals the hash embedded in its certificate, and the
reconstructed \(N\), minimum degree, edge count, and exact \(\mu\) agree.
Actual certificate hashes are:

- Certificate 00:
  `a625204b35fb852416a019529c999156fd8dd99cbd00a9df6b336678a0364137`;
  spec `79ffcfedb8bc51445d13c4d31ff88833ad548302e137b13e72dcec781d3d1a93`.
- Certificate 01:
  `68c5a8c81a36855c9a46edbcd73ed3e2a3aab9fce4e342b3f436cf46e0ae1fba`;
  spec `c846b1bde72f7dca5da0f1ced5ba633668975148fb2874e9be800ab344bf91c4`.
- Certificate 02:
  `b453c2a06ecf254dd9a6c9613ac398e289a818e7ac27ed4a3c9ff6d7d84d1450`;
  spec `0acbe760b5750f98723de8710900df924243cfd7ee99e9b5a6d11aeb94686dd3`.
- Certificate 03:
  `13a514fa0cc12e5d4457d1e358fa447e227f437ed5bc6b01c35c5a81ee1cc4da`;
  spec `a10ef4dee3a9060eed662296e17e9145ab1b225f719b70b4ed61f746f8d1a9e5`.
- Certificate 05:
  `906593e44b666a99069dc47f636f778d76128f0bbf49b438991518d6b21c838a`;
  spec `0a9fa225159442b260206a524a2ef40db9b613b7eb2693ff652844560350002e`.

The current certificate hash is
`9c0bed7ffe17248807460cad4fd80d59ce5edd61270f9129281a39d81b72704a`.

## Scope and caveats

- The root uniqueness certificate is local to the reported box.
- The exact one-zero statement concerns the full Hessian: one quotient
  rotation zero, seven positive quotient directions, and \(N-8\) positive
  transverse directions.
- The quotient and transverse lower bounds are rigorous outward interval
  statements.
- The nonlinear basin statement is numerical and class-level, as requested;
  it does not claim a rigorous full-dimensional basin radius.
- The preliminary gap in `construction_screen` is not the acceptance
  certificate.

## Reproduction and artifact paths

From `approach_162_independent_069152271_audit/`:

```sh
uv venv
uv pip install -r requirements.txt
.venv/bin/python independent_audit.py \
  ../approach_158_reflection_pair_expansion/artifacts/finite/m4_gamma_1e-06_lift_buffer_1p0005e-6_N200000000000000000000000000_candidate0074_ceil_graph_spec.json \
  --output artifacts
.venv/bin/python compare_artifacts.py
.venv/bin/python finalize_verdict.py
```

Primary outputs:

- `AUDIT_REPORT.md`
- `verdict.json`
- `artifacts/independent_core_audit.json`
- `artifacts/independent_matrices.txt`
- `artifacts/artifact_comparison.json`
- `independent_audit.py`
- `compare_artifacts.py`
- `finalize_verdict.py`
