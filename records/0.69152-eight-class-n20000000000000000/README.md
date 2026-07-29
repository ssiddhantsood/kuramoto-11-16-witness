# 0.691519878268 — eight-class witness

- Status: **FORMAL_INTERVAL_CERTIFIED**
- Audit: `producer_arb_and_independent_replay`
- Source approach: `156`
- Vertices: `20000000000000000`
- Minimum degree: `13830397565368202`
- Exact connectivity: `1975771080766886/2857142857142857`
- Decimal connectivity: `0.6915198782684101`
- Exact 11/16 cross-product excess: `1286361045891243`

## Contents

- `graph_spec.json` — compact exact graph and phase specification.
- `verify.py` — replay verifier.
- `record.json` — standardized metadata and status.
- Certificate, audit, and nonlinear JSON files when available.

## Verify

```bash
python3 verify.py --spec graph_spec.json --out verification_report.json
```
