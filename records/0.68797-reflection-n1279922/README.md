# 0.687967460492 — reflection witness

- Status: **CERTIFIED_REPLAY**
- Audit: `producer_interval_replay`
- Source approach: `155`
- Vertices: `1279922`
- Minimum degree: `880544`
- Exact connectivity: `880544/1279921`
- Decimal connectivity: `0.6879674604917022`
- Exact 11/16 cross-product excess: `9573`

## Contents

- `graph_spec.json` — compact exact graph and phase specification.
- `verify.py` — replay verifier.
- `record.json` — standardized metadata and status.
- Certificate, audit, and nonlinear JSON files when available.

## Verify

```bash
python3 verify.py --spec graph_spec.json --output verification_report.json --skip-fourier
```
