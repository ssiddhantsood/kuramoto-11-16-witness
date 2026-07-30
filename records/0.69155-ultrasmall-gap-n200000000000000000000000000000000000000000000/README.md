# 0.691553452131 — ultrasmall-gap witness

- Status: **FORMAL_INTERVAL_CERTIFIED**
- Audit: `producer_1024bit_arb_replay`
- Source approach: `167`
- Vertices: `200000000000000000000000000000000000000000000`
- Minimum degree: `138310690426227431945615663904064929607280523`
- Exact connectivity: `138310690426227431945615663904064929607280523/199999999999999999999999999999999999999999999`
- Decimal connectivity: `0.6915534521311372`
- Exact 11/16 cross-product excess: `12971046819638911129850622465038873716488379`

## Contents

- `graph_spec.json` — compact exact graph and phase specification.
- `verify.py` — replay verifier.
- `record.json` — standardized metadata and status.
- Certificate, audit, and nonlinear JSON files when available.

## Verify

```bash
python3 verify.py --spec graph_spec.json
```
