# N=80,002 reflection witness

This directory contains the current finite Kuramoto record witness. It is a
simple connected unweighted graph with a nonsynchronous equilibrium that is a
strict local minimum of the Kuramoto potential modulo global rotation.

## Exact graph

The six clique fibers are ordered as

```text
P1, M1, P2, M2, P3, M3
```

with sizes

```text
14133, 14133, 7024, 7024, 18844, 18844.
```

Every unordered pair of distinct fibers is specified exactly once in
[`graph_spec.json`](graph_spec.json) as complete, absent, or a consecutive
residue-orbit block. For a residue block between classes \(i<j\), local labels
\(u\) and \(v\) are adjacent exactly when

```text
(v - u - shift_start) mod gcd(n_i, n_j) < shift_count.
```

The resulting class degrees are

```text
55021, 55021, 55018, 55018, 56684, 56684.
```

There are exactly \(2{,}232{,}211{,}521\) edges. The minimum degree is
\(\delta=55{,}018\), so

\[
\mu=\frac{\delta}{N-1}
=\frac{55018}{80001}
=0.687716403544955688\ldots
>\frac{11}{16}.
\]

The comparison is exact:

\[
16(55018)-11(80001)=277>0.
\]

## Reflection rule

The graph involution swaps each positive/negative fiber pair,

```text
P1 <-> M1, P2 <-> M2, P3 <-> M3,
```

and maps local label \(u\) to \(-u\) modulo the fiber size.

The only orientation-sensitive case is the unequal-size P1-P3 block. Its
accepted residues modulo 4711 are

```text
{0} union {721, ..., 4710}.
```

Negation maps this set exactly to `{0, ..., 3990}`, which is the M1-M3 block
with `shift_start=0` and `shift_count=3991`. The standalone verifier checks
every residue difference in every reflected block.

## Equilibrium

The phase assignment is

\[
(\theta_1,-\theta_1,\theta_2,-\theta_2,\theta_3,-\theta_3),
\]

where

```text
theta1 = 0.707502292224047977399735004837712187588879040424...
theta2 = 0.969543644405493611756828635893715869390864151021...
theta3 = 2.360517597335587843038275804522053753812010721930...
```

Starting only from the short seeds in the graph specification, the verifier
solves the three positive-class torque equations at 180 decimal digits. It
also performs an outward-rounded Krawczyk check in a radius-\(10^{-10}\) box
about those seeds, proving existence and uniqueness within that box. This is
local uniqueness; the unrestricted trigonometric system has other roots,
including synchronous ones.

The independent audit found maximum six-torque residual below
\(8.5\times10^{-227}\). The order parameter is approximately
\(0.03329948\), so the equilibrium is nonsynchronous.

## Quotient spectrum

In mass-orthonormal class-constant coordinates, the six quotient eigenvalues
are approximately

```text
0,
0.4268309591,
23.2776753608,
621.3324579626,
20522.8961450791,
22942.7299176051.
```

Reflection splits this matrix into even and odd \(3\times3\) blocks. The
rotation zero lies in the even block and is the only zero eigenvalue. The
point normalized gap is

```text
5.3352536078706716e-6.
```

The independent interval audit gives the rigorous lower bound

```text
5.335072026813738e-6 > 1e-6.
```

## Complete transverse proof

The equitable partition splits the Hessian into the six-dimensional quotient
and \(N-6=79{,}996\) class-zero-sum dimensions.

On the transverse space:

- the clique on fiber \(i\) contributes \(n_i I\);
- complete cross blocks vanish between fiber-zero-sum vectors;
- a partial \((d_{ij},d_{ji})\)-biregular block has adjacency operator norm at
  most \(\sqrt{d_{ij}d_{ji}}\).

The verifier forms the resulting six-class comparison matrix and applies
outward-rounded Gershgorin bounds. Its rigorous lower bound is

```text
17725.84752950048.
```

This proves positivity of all 79,996 transverse modes without enumerating
them. The independent character audit separately enumerated every mode and
found minimum \(19244.984787192374\).

## Nonlinear qualification

The compact verifier intentionally omits the expensive nonlinear runs. Their
persisted independent verdict is:

```text
LOCAL RETURN CONFIRMED; BASIN SIZE NOT CERTIFIED
```

Weakest-mode perturbations of \(10^{-4}\) and \(3\times10^{-4}\) radians
returned in both signs. Perturbations of \(5\times10^{-4}\) and
\(2\times10^{-3}\) radians escaped to synchrony in both signs. Random
full-graph perturbations with vertexwise amplitudes through \(10^{-2}\)
returned because their projection onto the weak quotient direction was much
smaller. No global-attraction claim or explicit basin radius is made.

## Verify

From the repository root:

```bash
python3 records/n80002/verify.py
```

The command reads only the adjacent `graph_spec.json`, writes
`verification_report.json` in this directory, and exits nonzero if any check
fails. A successful run ends with:

```text
ACCEPTED N=80002 mu=55018/80001 ...
```

Machine-readable results are in
[`verification_report.json`](verification_report.json), and the independent
audit decision is in [`audit_verdict.json`](audit_verdict.json).
