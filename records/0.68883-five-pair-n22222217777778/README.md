# 0.688828353185 — five-pair witness

- Status: **CERTIFIED_REPLAY**
- Audit: `producer_interval_replay`
- Source approach: `158`
- Vertices: `22222217777778`
- Minimum degree: `15307293675990`
- Exact connectivity: `15307293675990/22222217777777`
- Decimal connectivity: `0.6888283531852447`
- Exact 11/16 cross-product excess: `472303260293`

## Contents

- `graph_spec.json` — compact exact graph and phase specification.
- `verify.py` — replay verifier.
- `record.json` — standardized metadata and status.
- Certificate, audit, and nonlinear JSON files when available.

## Verify

```bash
python3 verify.py --spec graph_spec.json --output verification_report.json
```
