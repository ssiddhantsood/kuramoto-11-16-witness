# Applying the machine-checked clique-blow-up theorem

This note records the exact bridge between the `N=80,002` certificate and the
Lean theorem
[`cliqueBlowup_preserves_stable_equilibrium`](https://github.com/ssiddhantsood/kuramoto-block-limit/blob/main/KuramotoBlockLimit/CliqueBlowup.lean).
It separates what is checked by interval arithmetic from what is checked by
Lean.

## Base hypotheses

Run:

```bash
python records/n80002/verify_arb.py
```

The resulting [`arb_certificate.json`](arb_certificate.json) establishes:

1. The compact residue specification defines a finite simple, unweighted,
   undirected graph with `N=80002` and minimum degree `55018`.
2. An outward-rounded Arb Krawczyk calculation encloses an exact
   nonsynchronous class-constant equilibrium.
3. The mass-orthonormal quotient Hessian has rotation as its only zero mode;
   every nonrotation quotient eigenvalue is greater than `2/5`.
4. The clique-fiber/partial-biregular comparison proves every one of the
   `79996` transverse Hessian modes is greater than `17000`.

Items 3 and 4 together say that the full base Hessian is positive
semidefinite with kernel exactly the global rotation vector. Its row sums are
zero by the cosine-Laplacian definition.

## Lean implication

For an arbitrary positive integer `m`, the Lean theorem takes precisely the
following hypotheses:

- the base coupling has no self-loops;
- the base phase is an exact equilibrium;
- the base Hessian quadratic form is nonnegative;
- the base Hessian kernel contains only constant vectors; and
- the Hessian row sums are zero.

Lean constructs the lexicographic clique coupling, proves that its matrix is
the actual Kuramoto potential Hessian at the replicated phase, and proves:

- the replicated phase is an exact equilibrium;
- the blown-up Hessian is positive semidefinite; and
- its kernel consists exactly of constant vectors.

Replication also preserves two distinct phase values, so the equilibrium
remains nonsynchronous.

## Exact family

For the base values `N=80002` and `delta=55018`, the graph `G[K_m]` has

```text
N_m     = 80002 m,
delta_m = 55019 m - 1.
```

Its exact excess over the `11/16` threshold is

```text
16 delta_m - 11(N_m-1) = 282 m - 5 > 0
```

for every integer `m>=1`. Thus the checked implication produces an infinite
family of stable nonsynchronous equilibria above `11/16`.

## Verification boundary

This is not a single end-to-end Lean proof: Arb certifies the concrete base
hypotheses, and Lean proves the universal blow-up implication. The short
bridge matching the certified properties to the Lean hypotheses is ordinary
mathematical reasoning. The final implication from a negative Jacobian on the
rotation quotient to local asymptotic stability modulo rotation is also a
standard finite-dimensional ODE result, not formalized here.
