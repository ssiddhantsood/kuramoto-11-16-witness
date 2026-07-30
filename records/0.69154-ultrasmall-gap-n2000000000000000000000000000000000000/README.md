# 0.691538348033 — ultrasmall-gap witness

- Status: **FORMAL_INTERVAL_CERTIFIED**
- Audit: `producer_1024bit_arb_replay`
- Source approach: `167`
- Vertices: `2000000000000000000000000000000000000`
- Minimum degree: `1383076696065173424880843046495688336`
- Exact connectivity: `1383076696065173424880843046495688336/1999999999999999999999999999999999999`
- Decimal connectivity: `0.6915383480325867`
- Exact 11/16 cross-product excess: `129227137042774798093488743931013387`

## Contents

- `graph_spec.json` — compact exact graph and phase specification.
- `verify.py` — replay verifier.
- `record.json` — standardized metadata and status.
- Certificate, audit, and nonlinear JSON files when available.

## Verify

```bash
python3 verify.py --spec graph_spec.json
```
