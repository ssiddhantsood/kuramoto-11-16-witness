# 0.691519878508 — eight-class witness

- Status: **FORMAL_INTERVAL_CERTIFIED**
- Audit: `producer_arb_and_independent_replay`
- Source approach: `156`
- Vertices: `200000000000000000`
- Minimum degree: `138303975701559707`
- Exact connectivity: `138303975701559707/199999999999999999`
- Decimal connectivity: `0.6915198785077985`
- Exact 11/16 cross-product excess: `12863611224955323`

## Contents

- `graph_spec.json` — compact exact graph and phase specification.
- `verify.py` — replay verifier.
- `record.json` — standardized metadata and status.
- Certificate, audit, and nonlinear JSON files when available.

## Verify

```bash
python3 verify.py --spec graph_spec.json --out verification_report.json
```
