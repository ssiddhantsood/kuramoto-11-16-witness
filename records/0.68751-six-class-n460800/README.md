# 0.687505832261 — six-class witness

- Status: **CONFIRMED_RIGOROUS**
- Audit: `independent_audit_and_arb_certificate`
- Source approach: `147`
- Vertices: `460800`
- Minimum degree: `316802`
- Exact connectivity: `316802/460799`
- Decimal connectivity: `0.6875058322609207`
- Exact 11/16 cross-product excess: `43`

## Contents

- `graph_spec.json` — compact exact graph and phase specification.
- `verify.py` — replay verifier.
- `record.json` — standardized metadata and status.
- Certificate, audit, and nonlinear JSON files when available.

## Verify

```bash
python3 verify.py --spec graph_spec.json --out verification_report.json
```
