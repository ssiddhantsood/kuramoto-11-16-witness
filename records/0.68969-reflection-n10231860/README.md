# 0.689689723050 — reflection witness

- Status: **CERTIFIED_REPLAY**
- Audit: `producer_interval_replay`
- Source approach: `155`
- Vertices: `10231860`
- Minimum degree: `7056808`
- Exact connectivity: `641528/930169`
- Decimal connectivity: `0.6896897230503274`
- Exact 11/16 cross-product excess: `358479`

## Contents

- `graph_spec.json` — compact exact graph and phase specification.
- `verify.py` — replay verifier.
- `record.json` — standardized metadata and status.
- Certificate, audit, and nonlinear JSON files when available.

## Verify

```bash
python3 verify.py --spec graph_spec.json --output verification_report.json --skip-fourier
```
