# 0.691522473075 — five-pair witness

- Status: **CERTIFIED_REPLAY**
- Audit: `producer_interval_replay`
- Source approach: `158`
- Vertices: `200000000000000000000000000`
- Minimum degree: `138304494614954012746623412`
- Exact connectivity: `138304494614954012746623412/199999999999999999999999999`
- Decimal connectivity: `0.6915224730747701`
- Exact 11/16 cross-product excess: `12871913839264203945974603`

## Contents

- `graph_spec.json` — compact exact graph and phase specification.
- `verify.py` — replay verifier.
- `record.json` — standardized metadata and status.
- Certificate, audit, and nonlinear JSON files when available.

## Verify

```bash
python3 verify.py --spec graph_spec.json --output verification_report.json
```
