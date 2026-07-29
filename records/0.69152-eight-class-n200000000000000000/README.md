# 0.691519878387 — eight-class witness

- Status: **FORMAL_INTERVAL_CERTIFIED**
- Audit: `producer_arb_and_independent_replay`
- Source approach: `156`
- Vertices: `200000000000000000`
- Minimum degree: `138303975677354405`
- Exact connectivity: `4769102609563945/6896551724137931`
- Decimal connectivity: `0.6915198783867720`
- Exact 11/16 cross-product excess: `12863610837670491`

## Contents

- `graph_spec.json` — compact exact graph and phase specification.
- `verify.py` — replay verifier.
- `record.json` — standardized metadata and status.
- Certificate, audit, and nonlinear JSON files when available.

## Verify

```bash
python3 verify.py --spec graph_spec.json --out verification_report.json
```
