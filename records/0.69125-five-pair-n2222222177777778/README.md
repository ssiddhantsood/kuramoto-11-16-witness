# 0.691249139012 — five-pair witness

- Status: **CERTIFIED_REPLAY**
- Audit: `producer_interval_replay`
- Source approach: `158`
- Vertices: `2222222177777778`
- Minimum degree: `1536109167082421`
- Exact connectivity: `1536109167082421/2222222177777777`
- Decimal connectivity: `0.6912491390120725`
- Exact 11/16 cross-product excess: `133302717763189`

## Contents

- `graph_spec.json` — compact exact graph and phase specification.
- `verify.py` — replay verifier.
- `record.json` — standardized metadata and status.
- Certificate, audit, and nonlinear JSON files when available.

## Verify

```bash
python3 verify.py --spec graph_spec.json --output verification_report.json
```
