# Kuramoto system contract

The stability certificate in this repository applies to the following system
and no broader class of models without a separate derivation.

## Certified dynamics

For the simple undirected graph with adjacency matrix `A` specified by
`records/n80002/graph_spec.json`, the phases obey

```text
d theta_i / dt = (K/N) sum_j A_ij sin(theta_j - theta_i),   K > 0.
```

The natural frequencies are identical. Their common value has been removed by
passing to a uniformly rotating frame. The factor `K/N` is positive, so omitting
it while constructing the Hessian changes only the time scale and not any
stability sign.

The associated potential is

```text
V(theta) = (K/N) sum over edges (i,j) of [1 - cos(theta_i - theta_j)].
```

Consequently the dynamics are `d theta/dt = -grad(V)`, and the flow Jacobian at
an equilibrium is the negative Hessian of `V`.

## Not covered by this certificate

The current certificate does not automatically apply to:

- heterogeneous natural frequencies that cannot be removed in one rotating
  frame;
- directed or asymmetric coupling;
- degree-normalized coupling;
- negative coupling strength;
- inertial or other second-order oscillator equations;
- amplitude dynamics; or
- a six-node weighted surrogate in place of the full `N=80,002` graph.

Those systems may have related stability tests, but their Jacobians must be
derived separately.

## Interpretation of the conclusion

For the certified dynamics, the exact enclosed nonsynchronous equilibrium has
a positive Hessian in every direction except global phase rotation. Therefore
it is locally asymptotically (and locally exponentially, after quotienting the
rotation symmetry) stable. This is not a global-attraction claim and does not
certify a basin radius.

External model confirmation remains a human/interface requirement: the person
requesting the check must confirm that the equations above are the “required
system.” Code can verify consequences of an equation, but cannot infer an
unstated physical model.
