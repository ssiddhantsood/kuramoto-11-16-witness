# 0.691519878268 — eight-class witness

- Status: **FORMAL_INTERVAL_CERTIFIED**
- Audit: `producer_arb_and_independent_replay`
- Source approach: `156`
- Vertices: `40000000000000000`
- Minimum degree: `27660795130736401`
- Exact connectivity: `27660795130736401/39999999999999999`
- Decimal connectivity: `0.6915198782684100`
- Exact 11/16 cross-product excess: `2572722091782427`

## Contents

- `graph_spec.json` — compact exact graph and phase specification.
- `verify.py` — replay verifier.
- `record.json` — standardized metadata and status.
- Certificate, audit, and nonlinear JSON files when available.

## Verify

```bash
python3 verify.py --spec graph_spec.json --out verification_report.json
```
