# 0.691519878032 — eight-class witness

- Status: **FORMAL_INTERVAL_CERTIFIED**
- Audit: `producer_arb_and_independent_replay`
- Source approach: `156`
- Vertices: `100000000000000000`
- Minimum degree: `69151987803168581`
- Exact connectivity: `69151987803168581/99999999999999999`
- Decimal connectivity: `0.6915198780316858`
- Exact 11/16 cross-product excess: `6431804850697307`

## Contents

- `graph_spec.json` — compact exact graph and phase specification.
- `verify.py` — replay verifier.
- `record.json` — standardized metadata and status.
- Certificate, audit, and nonlinear JSON files when available.

## Verify

```bash
python3 verify.py --spec graph_spec.json --out verification_report.json
```
