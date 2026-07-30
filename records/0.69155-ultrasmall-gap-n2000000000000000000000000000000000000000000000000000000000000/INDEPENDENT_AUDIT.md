# Adversarial independent audit of the `gamma=1e-12` six-class record

## Verdict

**CONFIRMED.**

The graph specification defines a finite simple connected unweighted graph,
the stated reflection is an exact automorphism, the independently resolved
equilibrium is locally unique and nonsynchronous, the normalized Hessian has
exactly one zero, and every other mode is strictly positive.  The exact
density beats the certified approach-158 and approach-156 records.

The independent phase read only
`approach_167_ultrasmall_gap_continuation/finite/gamma_1em12_graph_spec.json`.
It produced `independent_core.json` before any producer certificate, report,
script, or prior record was read.  Producer and prior-record comparisons were
performed afterward by the separate `posthoc_compare.py`.

## 1. Exact graph reconstruction

For a residue block `(i,j)` with modulus `m`, start `s`, and count `k`, the
independent reconstruction uses the declarative edge predicate

```text
(v mod m - u mod m) mod m in {s,s+1,...,s+k-1} mod m.
```

All integer arithmetic was performed with arbitrary-size integers.  The
result is:

```text
N = 2000000000000000000000000000000000000000000000000000000000000

class sizes =
383107521134247498152220474771000000000000000000000000000000 (twice)
473466121294882045140591845053000000000000000000000000000000 (twice)
143426357570870456707187680176000000000000000000000000000000 (twice)

class degrees =
1383107521134247498152220474771356221878069981885091586880410 (twice)
1422991521438334283841520057649123224781458256311876526628430 (twice)
1383107521134247498152220474770999999999999999999999999999999 (twice)

minimum degree =
1383107521134247498152220474770999999999999999999999999999999

edge count =
1401991244059947364597948557895742860151012077301487991307104091431014405180239156876112724000000000000000000000000000000
```

The handshake identity holds exactly.  Each residue block is exactly
biregular, every modulus divides both incident class sizes, every shift set
has no repeated residue, and all 15 unordered class pairs occur exactly once.
Thus the construction is simple and unweighted.  The nonempty class-support
graph is connected; because every class is a clique and every support edge
contains graph edges, the full graph is connected.

The exact density is

```text
1383107521134247498152220474770999999999999999999999999999999
/
1999999999999999999999999999999999999999999999999999999999999

= 0.69155376056712374907611023738549999999999999999999999999999984577688028356187453805511869275...
```

For reflection, the local map is `u -> -u mod class_size`.  Complete and
absent blocks map to their declared mates.  On residue blocks, endpoint
reversal preserves the shift set, while preserved endpoint order negates it.
In particular,

```text
-(131750480144869275378358570540 ... cyclic count 868249519855130724621641429461)
= (0 ... cyclic count 868249519855130724621641429461) mod 10^30.
```

The exact modular identity passes for every block, with no vertex sampling.

## 2. Independent root and Krawczyk certificate

An independent 420-decimal-digit Newton solve reached normalized torque
residual at most `1.81e-422`.  The three reflection-centered angles are:

```text
-0.729628015634143611366870682220297220138016422425426703493546765252914190456679017425543574827687553700711004582438147074997892429418964675508515713454643514456894053345690705003601244221529399908889021600346992211681460363453265751640608153750839659751337315187885531957983856143745182360347982483663496016378571841945915146673701706701985579833801758107359920

-2.35874477038387634610503260043145267139026514540585371601666726341128960565577607793628542208982595704434808325671767759057069919800939054681281552808297595882975364418125866874699516956066048895097543524819581342379754763801987649320941355872481709327555832827722959799844306855465053493755849447331717599445254427253596862611238993044530915140464264400406987

 0.956552492186795164818884573003263023660329791146266272258713715323159113215929185929190694187856697156966852649677178841955836870724508166891166372628566059075915013746418059187582522799719731272826863889215697105259177025731332381004620754126940662795827141271147129551144234054622529515432701245845300102819561301904173284605717808328117488325951542643131622
```

Each center has certified outer radius `1.00000001e-330`.  At 2048 bits,
Arb evaluates a fixed rational-preconditioner Krawczyk map.  Its image radii
are at most `4.78e-380`, and all three image intervals lie strictly inside
the root box.  The rational preconditioner has nonzero exact determinant.
This proves existence and uniqueness of the reflected torque root in the
box.

All three pair masses are exactly distinct.  Arb proves all 15 quantities
`|sin((phi_i-phi_j)/2)|` are positive; the smallest lower value is about
`0.1132189488485467`.  The order parameter is real by reflection and lies in

```text
0.03259556072406864644480342924782100478927426017418894793656143095256...
```

with a `2.81e-102` displayed enclosure radius.  Hence all six phases are
distinct modulo `2*pi` and the equilibrium is nonsynchronous.

The quotient Hessian certificate below, together with the nonlinear
strong-monotonicity ball, also rules out non-reflected nearby quotient roots
modulo global rotation.

## 3. Even, odd, and full quotient

The normalization is `H=-J` divided by `N`.  The mass-symmetric form is
constructed from exact detailed-balance edge densities before any square
roots are taken.  Arb verifies `J=-H`, exact quotient row sums zero, and
commutation with reflection.  The exact reflection decomposition gives:

```text
even:
0 exactly
9.99999992474819970240407568405408540558935507498057392303224797481938...e-13
0.256292574625478300455437949262054281685369990530820027131238589080906...

odd:
1.00000000760698683793082703521121608394699446337419381065449160362255...e-12
0.005149318842718618051010461324596537746401529347283768652878609223812...
0.287180293040303007560506511496740483689945423078804858900281616293823...
```

The full certified spectrum is the union of these parity spectra.  The even
stiffness has the rotation vector in its kernel exactly.  After restricting
to the exact mass-orthogonal tangent basis, Arb/Sylvester proves positivity
after the explicit rational shift below.  The odd block passes the same
shifted Sylvester test.  Therefore the zero multiplicity is exactly one.

The actual positive gap is enclosed by the explicit rational decimal bounds

```text
0.0000000000009999999924748199702404075684054085405589355074980573923032247974819380015549601457103430093402332321092745682542329765408164545
<
lambda_gap
<
0.00000000000099999999247481997024040756840540854055893550749805739230322479748193800155496014571034300934023323210927456825423297654081645453.
```

A shorter explicit strict rational lower bound is

```text
lambda_gap >
1999999984949639940480815136810817081
/
2000000000000000000000000000000000000000000000000

= 0.0000000000009999999924748199702404075684054085405.
```

This is a bound against zero, not against a fixed `1e-6` threshold.

Roundoff is not close to the sign decision:

- Arb precision: 2048 bits (at least the requested 1536 bits);
- certified root radius: at most `1.00000001e-330`;
- resulting gap-ball radius: about `1.825e-330`;
- gap scale: about `1e-12`;
- slack above the stated rational lower: about `5.894e-50`.

Thus the positive sign has roughly 280 decimal orders of interval safety
relative to the computed ball radius.

## 4. Fourier-free transverse certificate

Let each class component have zero coordinate sum.  Complete bipartite
blocks annihilate such components.  For a residue adjacency block `A_S`,
the complete block also annihilates the input, so

```text
A_S = -A_(S complement)
```

on the relevant space.  A row/column Schur bound therefore gives, without
Fourier characters,

```text
||A_S|| <= min(k,m-k) * sqrt(copies_i * copies_j).
```

After multiplying by the outward-rounded cosine weights, the three nonzero
residue-block norm bounds are approximately

```text
0.00451063484999762715570112971275763814077689882237
0.00163531524823937028409293161519317456456806302294
0.00163531524823937028409293161519317456456806302294.
```

The resulting six-by-six scalar comparison matrix has certified smallest
eigenvalue

```text
0.240080898804678670560426900235102625903539942773935438480517054399609...
```

and a shifted six-minor Sylvester certificate proves the explicit bound

```text
lambda_transverse >
1200404494023393352802134501175513129
/
5000000000000000000000000000000000000

= 0.2400808988046786705604269002351026258.
```

This covers every transverse mode.  The exact dimension split is:

```text
quotient dimension                 = 6
transverse dimension               = N-6
                                   = 1999999999999999999999999999999999999999999999999999999999994
reflection-even transverse         = N/2-3
                                   = 999999999999999999999999999999999999999999999999999999999997
reflection-odd transverse          = N/2-3
                                   = 999999999999999999999999999999999999999999999999999999999997
```

The six classwise zero-sum component dimensions are listed exactly in
`transverse_certificate.json`.

## 5. Local nonlinear return

This audit uses a nonlinear analytic certificate rather than relying only on
trajectory plots.  In weighted quotient norm, the graph-Laplacian Hessian
variation obeys

```text
||Delta H|| <= L ||delta||_p,
L < 10.627541787759557 < 11.
```

Combining `L<11` with the strict quotient-gap lower proves uniform positive
curvature throughout the deliberately tiny ball

```text
||delta||_p <= 1e-14.
```

The corresponding phase sup-norm radius is below
`3.734225266858031e-14`.  Strong monotonicity makes this ball forward
invariant and proves convergence to the root modulo rotation.  For both
signed gap-mode initial conditions

```text
delta_i(0) = +/- 1e-15 * v_i/sqrt(p_i),
```

the weighted norm is `1e-15`.  With gap-scaled time
`tau = lambda_gap * t`, the rigorous decay estimate is

```text
||delta(tau)||_p / ||delta(0)||_p
<= exp(-0.8899999991722301904973279928171854464123... * tau).
```

Consequently the ratio is at most `1.8601939574878001e-8` at `tau=20`.
Both signs return.  No macroscopic basin is claimed.

## 6. Post-independence comparisons

The exact improvements are:

```text
over approach 158:
12419862729497800391777571092478865752501847779525228999938304542182060000915026159
/
399999999999999999999999997999999999999999999999999999999999800000000000000000000000001
= 0.0000310496568237445009794439278864454484999771243460...

over approach 156:
13552823730064247336573819952501847779525228999999999999938303975701559718
/
399999999999999997999999999999999999999999999999999999999999800000000000000001
= 0.0000338820593251606185108448465070577120030373050353...
```

Against the direct strict-branch zero-gap endpoint

```text
0.6915537605979708727851504507416565918940086772395606714859289344844874,
```

the finite record remains lower by

```text
0.0000000000308471237090402133561565918940086772395606714859290887076071...
```

as expected for a strictly positive-gap point.

The producer graph hash, exact degree data, edge count, root, parity spectra,
rigorous gap lower, transverse lower, and acceptance fields agree.  The
producer's 132-place displayed point estimate for the gap sits
`1.835455e-133` above the audit's 140-place enclosure, but this is below its
display ULP `1e-132`; the two enclosures overlap when that display precision
is respected.  This is not a certificate discrepancy.

The other four approach-167 records (`5e-7`, `1e-7`, `1e-8`, and `1e-10`)
were spot-checked after the independent target audit.  All eight recomputed
spec/certificate SHA-256 hashes match the bundle, exact graph arithmetic
matches each certificate, reflection and connectivity checks pass, and the
five exact densities form a strictly increasing sequence.

## Artifacts

- `independent_audit.py` — spec-only exact and 2048-bit Arb verifier.
- `graph_arithmetic.json` — exact block, degree, edge, reflection, and density arithmetic.
- `root_krawczyk.json` — 420-digit solve and 2048-bit Krawczyk/phase certificate.
- `quotient_certificate.json` — even, odd, full quotient and strict gap certificate.
- `transverse_certificate.json` — Fourier-free all-mode comparison and dimensions.
- `nonlinear_return.json` — rigorous tiny-basin nonlinear return.
- `independent_core.json` — persisted handoff from the independent phase.
- `posthoc_compare.py` — producer/prior comparison, run only afterward.
- `producer_comparison.json` — target producer comparison.
- `prior_record_comparison.json` — exact approach-158/156 improvements and branch headroom.
- `hash_spot_checks.json` — four-record hash and exact-arithmetic spot checks.
- `requirements.txt` — pinned independent numerical dependencies.
- `verdict.json` — machine-readable final verdict.
