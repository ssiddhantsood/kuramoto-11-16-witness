# 0.691120222735 — reflection witness

- Status: **CONFIRMED**
- Audit: `independent_adversarial_audit`
- Source approach: `155`
- Vertices: `39945316`
- Minimum degree: `27607015`
- Exact connectivity: `5521403/7989063`
- Decimal connectivity: `0.6911202227345059`
- Exact 11/16 cross-product excess: `2313775`

## Contents

- `graph_spec.json` — compact exact graph and phase specification.
- `verify.py` — replay verifier.
- `record.json` — standardized metadata and status.
- Certificate, audit, and nonlinear JSON files when available.

## Verify

```bash
python3 verify.py --spec graph_spec.json --output verification_report.json --skip-fourier
```
