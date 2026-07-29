# 0.691522710910 — four-pair witness

- Status: **CONFIRMED**
- Audit: `independent_1024bit_arb_audit`
- Source approach: `158`
- Vertices: `200000000000000000000000000`
- Minimum degree: `138304542182060000915026158`
- Exact connectivity: `138304542182060000915026158/199999999999999999999999999`
- Decimal connectivity: `0.6915227109103000`
- Exact 11/16 cross-product excess: `12872674912960014640418539`

## Contents

- `graph_spec.json` — compact exact graph and phase specification.
- `verify.py` — replay verifier.
- `record.json` — standardized metadata and status.
- Certificate, audit, and nonlinear JSON files when available.

## Verify

```bash
python3 verify.py --spec graph_spec.json --output verification_report.json
```
