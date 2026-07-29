# 0.688746861128 — reflection witness

- Status: **CERTIFIED_REPLAY**
- Audit: `producer_interval_replay`
- Source approach: `155`
- Vertices: `999962`
- Minimum degree: `688720`
- Exact connectivity: `688720/999961`
- Decimal connectivity: `0.6887468611275840`
- Exact 11/16 cross-product excess: `19949`

## Contents

- `graph_spec.json` — compact exact graph and phase specification.
- `verify.py` — replay verifier.
- `record.json` — standardized metadata and status.
- Certificate, audit, and nonlinear JSON files when available.

## Verify

```bash
python3 verify.py --spec graph_spec.json --output verification_report.json --skip-fourier
```
