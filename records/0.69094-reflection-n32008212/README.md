# 0.690935460279 — reflection witness

- Status: **CERTIFIED_REPLAY**
- Audit: `producer_interval_replay`
- Source approach: `155`
- Vertices: `32008212`
- Minimum degree: `22115608`
- Exact connectivity: `22115608/32008211`
- Decimal connectivity: `0.6909354602792390`
- Exact 11/16 cross-product excess: `1759407`

## Contents

- `graph_spec.json` — compact exact graph and phase specification.
- `verify.py` — replay verifier.
- `record.json` — standardized metadata and status.
- Certificate, audit, and nonlinear JSON files when available.

## Verify

```bash
python3 verify.py --spec graph_spec.json --output verification_report.json --skip-fourier
```
