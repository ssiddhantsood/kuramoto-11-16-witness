# 0.691553760567 — ultrasmall-gap witness

- Status: **CONFIRMED**
- Audit: `independent_2048bit_arb_audit`
- Source approach: `167`
- Vertices: `2000000000000000000000000000000000000000000000000000000000000`
- Minimum degree: `1383107521134247498152220474770999999999999999999999999999999`
- Exact connectivity: `1383107521134247498152220474770999999999999999999999999999999/1999999999999999999999999999999999999999999999999999999999999`
- Decimal connectivity: `0.6915537605671237`
- Exact 11/16 cross-product excess: `129720338147959970435527596335999999999999999999999999999995`

## Contents

- `graph_spec.json` — compact exact graph and phase specification.
- `verify.py` — replay verifier.
- `record.json` — standardized metadata and status.
- Certificate, audit, and nonlinear JSON files when available.

## Verify

```bash
python3 verify.py --spec graph_spec.json
```
