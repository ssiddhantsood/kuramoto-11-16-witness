# Rigorous interval certificate

Status: **CONFIRMED_RIGOROUS**.

`verify.py` reads only `../graph_spec.json`. It uses exact integer/rational
arithmetic and outward-rounded Arb intervals from `python-flint==0.9.0` at
512-bit precision. It imports no producer or repository verifier code.

Run from the repository root:

```sh
python3 -m venv /tmp/kuramoto-formal-venv
/tmp/kuramoto-formal-venv/bin/pip install -r requirements.txt
/tmp/kuramoto-formal-venv/bin/python formal/verify.py
```

The verifier always writes generated artifacts to `formal/certificate.json`
and `formal/verifier_output.txt`.

## Certified result

- `N=460800`, with class sizes
  `(80325,80325,108675,108675,41400,41400)`.
- The graph is simple, unweighted, undirected, and connected.
- `|E|=74003856975`.
- Class degrees are
  `(316835,316835,326096,326096,316802,316802)`.
- The exact minimum-degree ratio is
  `316802/460799 > 11/16`, since
  `16*316802-11*460799=43`.
- A unique nonsynchronous equilibrium exists in the rational box

```text
s =  1.40739194443017506833155213230129895323516278019010739034754292834847506656553860411541657 +/- 1e-40
x = -1.65671857734101766114395332154205876503218108062580714800181493665341901327426827868627599 +/- 1e-40
y = -0.268726453317482380461178982598739664045080336777632085737417667463948780917056917786449521 +/- 1e-40
```

  for `phi=(0,s,x,s-x,s-y,y)`, with class 2 wrapped as `x+2*pi`.
- The weighted order parameter satisfies `r<1/30`.
- The mass-orthonormal quotient Hessian has one exact rotation zero and
  nonzero gap `>487/10=48.7`.
- Every transverse eigenvalue is `>102000`.
- Therefore the full Hessian is positive semidefinite with kernel exactly
  `span(1)`. The equilibrium is a strict local minimum and locally
  asymptotically stable modulo global rotation.

## Proof outline

For a residue block with modulus `m` and `k` selected residues, divisibility
of both class sizes by `m` gives directed degrees
`d_ij=k*n_j/m`, `d_ji=k*n_i/m`, and
`n_i*d_ij=n_j*d_ji`. Exact support connectivity follows from clique fibers
and the support spanning tree `(0,1),(0,2),(0,4),(0,5),(1,3)`.

Reflection reduces equilibrium to the torques of classes `0,2,5`. A fixed
rational preconditioner gives an interval contraction
`||I-CJ(X)||_infinity<1e-34`, and its centered Krawczyk image lies strictly
inside `X`. Banach's theorem proves existence and uniqueness; reflection
forces all six torques to vanish.

For the quotient, Arb interval `LDL^T` proves the grounded principal matrix of
`Q-(487/10)P` positive definite, where `P` projects off the exact rotation
vector. For transverse vectors, clique fibers contribute `n_i I`, complete
cross blocks vanish on zero-sum vectors, and every partial biregular adjacency
satisfies
`||A_ij||_2<=sqrt(d_ij*d_ji)`. Interval `LDL^T` then proves the resulting
comparison matrix exceeds `102000 I`. No simulation or Fourier completeness
assumption is used.

Exact rational inputs, interval endpoints, contraction bounds, matrices, and
`LDL^T` pivots are recorded in `certificate.json`.
