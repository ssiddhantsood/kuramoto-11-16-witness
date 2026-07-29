# 0.690109034599 — reflection witness

- Status: **CERTIFIED_REPLAY**
- Audit: `producer_interval_replay`
- Source approach: `155`
- Vertices: `20480014`
- Minimum degree: `14133442`
- Exact connectivity: `14133442/20480013`
- Decimal connectivity: `0.6901090345987574`
- Exact 11/16 cross-product excess: `854929`

## Contents

- `graph_spec.json` — compact exact graph and phase specification.
- `verify.py` — replay verifier.
- `record.json` — standardized metadata and status.
- Certificate, audit, and nonlinear JSON files when available.

## Verify

```bash
python3 verify.py --spec graph_spec.json --output verification_report.json --skip-fourier
```
