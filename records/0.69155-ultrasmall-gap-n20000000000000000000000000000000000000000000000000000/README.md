# 0.691553757513 — ultrasmall-gap witness

- Status: **FORMAL_INTERVAL_CERTIFIED**
- Audit: `producer_1024bit_arb_replay`
- Source approach: `167`
- Vertices: `20000000000000000000000000000000000000000000000000000`
- Minimum degree: `13831075150265178757390842399999999999999999999999999`
- Exact connectivity: `13831075150265178757390842399999999999999999999999999/19999999999999999999999999999999999999999999999999999`
- Decimal connectivity: `0.6915537575132590`
- Exact 11/16 cross-product excess: `1297202404242860118253478399999999999999999999999995`

## Contents

- `graph_spec.json` — compact exact graph and phase specification.
- `verify.py` — replay verifier.
- `record.json` — standardized metadata and status.
- Certificate, audit, and nonlinear JSON files when available.

## Verify

```bash
python3 verify.py --spec graph_spec.json
```
