# 0.691249367157 — four-pair witness

- Status: **CERTIFIED_REPLAY**
- Audit: `producer_interval_replay`
- Source approach: `158`
- Vertices: `20000000000000000`
- Minimum degree: `13824987343133142`
- Exact connectivity: `13824987343133142/19999999999999999`
- Decimal connectivity: `0.6912493671566571`
- Exact 11/16 cross-product excess: `1199797490130283`

## Contents

- `graph_spec.json` — compact exact graph and phase specification.
- `verify.py` — replay verifier.
- `record.json` — standardized metadata and status.
- Certificate, audit, and nonlinear JSON files when available.

## Verify

```bash
python3 verify.py --spec graph_spec.json --output verification_report.json
```
