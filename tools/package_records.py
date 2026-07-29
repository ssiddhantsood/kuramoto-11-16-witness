#!/usr/bin/env python3
"""Package every finite witness into a uniform, sorted record catalog."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from fractions import Fraction
from pathlib import Path


REPO = Path(__file__).resolve().parent.parent


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def copy_if(source: Path | None, target: Path) -> None:
    if source is not None and source.exists():
        shutil.copy2(source, target)


def read_json(path: Path) -> dict:
    return json.loads(path.read_text())


def make_slug(mu: Fraction, architecture: str, n: int, degree: int, used: set) -> str:
    base = f"{float(mu):.5f}-{architecture}-n{n}"
    slug = base if base not in used else f"{base}-d{degree}"
    used.add(slug)
    return slug


def first_records(source_root: Path) -> list[dict]:
    records = []
    records.append(
        {
            "source_approach": 147,
            "architecture": "six-class",
            "n": 460800,
            "degree": 316802,
            "mu": Fraction(316802, 460799),
            "status": "CONFIRMED_RIGOROUS",
            "audit_status": "independent_audit_and_arb_certificate",
            "spec": REPO / "graph_spec.json",
            "verifier_kind": "first",
            "report": REPO / "verification_report.json",
            "audit": REPO / "audit_verdict.json",
            "formal": REPO / "formal" / "certificate.json",
            "nonlinear": None,
        }
    )
    records.append(
        {
            "source_approach": 154,
            "architecture": "reflection",
            "n": 80002,
            "degree": 55018,
            "mu": Fraction(55018, 80001),
            "status": "CONFIRMED",
            "audit_status": "independent_adversarial_audit",
            "spec": REPO / "records" / "n80002" / "graph_spec.json",
            "verifier_kind": "reflection",
            "report": REPO / "records" / "n80002" / "verification_report.json",
            "audit": REPO / "records" / "n80002" / "audit_verdict.json",
            "formal": None,
            "nonlinear": (
                source_root
                / "approach_154_reflection_weighted_realization"
                / "artifacts"
                / "nonlinear_simulations.json"
            ),
        }
    )
    return records


def approach_155_records(source_root: Path) -> list[dict]:
    root = source_root / "approach_155_six_class_density_frontier"
    audit_root = source_root / "approach_159_independent_069112_audit"
    ledger = read_json(root / "artifacts" / "record_ledger.json")
    records = []
    for row in ledger["records"]:
        spec = Path(row["spec_path"])
        record_dir = spec.parent
        independently_confirmed = (
            row["weighted_record"] == "gamma_1e-06"
            and (audit_root / "verdict.json").exists()
        )
        records.append(
            {
                "source_approach": 155,
                "architecture": "reflection",
                "n": int(row["vertex_count"]),
                "degree": int(row["mu_numerator"]),
                "mu": Fraction(int(row["mu_numerator"]), int(row["mu_denominator"])),
                "status": (
                    "CONFIRMED" if independently_confirmed else "CERTIFIED_REPLAY"
                ),
                "audit_status": (
                    "independent_adversarial_audit"
                    if independently_confirmed
                    else "producer_interval_replay"
                ),
                "spec": spec,
                "verifier_kind": "six-class",
                "report": record_dir / "independent_replay.json",
                "audit": (
                    audit_root / "verdict.json"
                    if independently_confirmed
                    else None
                ),
                "audit_report": (
                    audit_root / "AUDIT_REPORT.md"
                    if independently_confirmed
                    else None
                ),
                "formal": None,
                "nonlinear": record_dir / "nonlinear_return.json",
            }
        )
    return records


def approach_156_records(source_root: Path) -> list[dict]:
    root = source_root / "approach_156_bottleneck_fiber_splitting"
    audit_root = source_root / "approach_160_independent_0691519_audit"
    ledger = read_json(root / "artifacts" / "record_ledger.json")
    records = []
    seen_specs = set()
    for row in ledger["records"]:
        spec = root / row["graph_spec"]
        digest = sha256(spec)
        if digest in seen_specs:
            continue
        seen_specs.add(digest)
        independently_confirmed = (
            row["name"].endswith("d138303975701559717")
            and (audit_root / "verdict.json").exists()
        )
        records.append(
            {
                "source_approach": 156,
                "architecture": "eight-class",
                "n": int(row["N"]),
                "degree": int(row["minimum_degree"]),
                "mu": Fraction(int(row["mu_numerator"]), int(row["mu_denominator"])),
                "status": (
                    "CONFIRMED"
                    if independently_confirmed
                    else "FORMAL_INTERVAL_CERTIFIED"
                ),
                "audit_status": (
                    "independent_768bit_arb_audit"
                    if independently_confirmed
                    else "producer_arb_and_independent_replay"
                ),
                "spec": spec,
                "verifier_kind": "eight-class",
                "report": (
                    root / row["independent_replay"]
                    if row.get("independent_replay")
                    else None
                ),
                "audit": (
                    audit_root / "verdict.json"
                    if independently_confirmed
                    else None
                ),
                "audit_report": (
                    audit_root / "AUDIT_REPORT.md"
                    if independently_confirmed
                    else None
                ),
                "formal": root / row["formal_interval_certificate"],
                "formal_output": root / row["formal_verifier_output"],
                "numerical": root / row["numerical_certificate"],
                "nonlinear": (
                    root / row["nonlinear_certificate"]
                    if row.get("nonlinear_certificate")
                    else None
                ),
            }
        )
    return records


def package_verifier(record: dict, target: Path, source_root: Path) -> str:
    kind = record["verifier_kind"]
    if kind == "first":
        shutil.copy2(REPO / "verify_witness.py", target / "verify.py")
        return "python3 verify.py --spec graph_spec.json --out verification_report.json"
    if kind == "reflection":
        shutil.copy2(REPO / "records" / "n80002" / "verify.py", target / "verify.py")
        return "python3 verify.py"
    if kind == "six-class":
        base = (
            source_root
            / "approach_154_reflection_weighted_realization"
            / "replay_witness.py"
        )
        wrapper = (
            source_root
            / "approach_155_six_class_density_frontier"
            / "scripts"
            / "replay_record.py"
        )
        shutil.copy2(base, target / "base_replay.py")
        text = wrapper.read_text()
        marker = "def load_base_replay() -> ModuleType:"
        text = text.replace(marker, 'BASE_REPLAY = HERE / "base_replay.py"\n\n\n' + marker)
        (target / "verify.py").write_text(text)
        return (
            "python3 verify.py --spec graph_spec.json "
            "--output verification_report.json --skip-fourier"
        )
    if kind == "eight-class":
        source = (
            source_root
            / "approach_156_bottleneck_fiber_splitting"
            / "scripts"
            / "replay_record.py"
        )
        shutil.copy2(source, target / "verify.py")
        return "python3 verify.py --spec graph_spec.json --out verification_report.json"
    raise ValueError(kind)


def package_record(record: dict, target: Path, source_root: Path, slug: str) -> dict:
    target.mkdir(parents=True, exist_ok=True)
    shutil.copy2(record["spec"], target / "graph_spec.json")
    command = package_verifier(record, target, source_root)
    copy_if(record.get("report"), target / "verification_report.json")
    copy_if(record.get("audit"), target / "audit_verdict.json")
    copy_if(record.get("audit_report"), target / "INDEPENDENT_AUDIT.md")
    copy_if(record.get("formal"), target / "formal_certificate.json")
    copy_if(record.get("formal_output"), target / "formal_verifier_output.txt")
    copy_if(record.get("numerical"), target / "numerical_certificate.json")
    copy_if(record.get("nonlinear"), target / "nonlinear_results.json")

    mu = record["mu"]
    metadata = {
        "schema": "kuramoto-finite-record-v1",
        "slug": slug,
        "source_approach": record["source_approach"],
        "architecture": record["architecture"],
        "status": record["status"],
        "audit_status": record["audit_status"],
        "vertex_count": record["n"],
        "minimum_degree": record["degree"],
        "mu_exact": f"{mu.numerator}/{mu.denominator}",
        "mu_decimal": float(mu),
        "excess_over_11_16_cross_product": 16 * record["degree"] - 11 * (record["n"] - 1),
        "graph_spec_sha256": sha256(target / "graph_spec.json"),
        "verify_command": command,
    }
    (target / "record.json").write_text(json.dumps(metadata, indent=2) + "\n")
    (target / "README.md").write_text(
        f"# {float(mu):.12f} — {record['architecture']} witness\n\n"
        f"- Status: **{record['status']}**\n"
        f"- Audit: `{record['audit_status']}`\n"
        f"- Source approach: `{record['source_approach']}`\n"
        f"- Vertices: `{record['n']}`\n"
        f"- Minimum degree: `{record['degree']}`\n"
        f"- Exact connectivity: `{mu.numerator}/{mu.denominator}`\n"
        f"- Decimal connectivity: `{float(mu):.16f}`\n"
        f"- Exact 11/16 cross-product excess: "
        f"`{metadata['excess_over_11_16_cross_product']}`\n\n"
        "## Contents\n\n"
        "- `graph_spec.json` — compact exact graph and phase specification.\n"
        "- `verify.py` — replay verifier.\n"
        "- `record.json` — standardized metadata and status.\n"
        "- Certificate, audit, and nonlinear JSON files when available.\n\n"
        "## Verify\n\n"
        f"```bash\n{command}\n```\n"
    )
    return metadata


def write_index(rows: list[dict]) -> None:
    lines = [
        "# Finite witness catalog",
        "",
        "Every folder below defines a finite simple unweighted graph with a stable "
        "nonsynchronous equilibrium and exact minimum-degree ratio above `11/16`.",
        "",
        "Statuses distinguish independently confirmed records from producer-side "
        "interval/replay certificates whose adversarial audit is pending.",
        "",
        "| μ | N | architecture | status | folder |",
        "|---:|---:|---|---|---|",
    ]
    for row in sorted(
        rows,
        key=lambda item: Fraction(
            int(item["mu_exact"].split("/")[0]),
            int(item["mu_exact"].split("/")[1]),
        ),
    ):
        lines.append(
            f"| `{row['mu_decimal']:.12f}` | `{row['vertex_count']}` | "
            f"{row['architecture']} | `{row['status']}` / `{row['audit_status']}` | "
            f"[`{row['slug']}/`]({row['slug']}/README.md) |"
        )
    (REPO / "records" / "README.md").write_text("\n".join(lines) + "\n")
    (REPO / "records" / "index.json").write_text(
        json.dumps(
            {
                "schema": "kuramoto-record-catalog-v1",
                "threshold": "11/16",
                "records": sorted(rows, key=lambda item: item["mu_decimal"]),
            },
            indent=2,
        )
        + "\n"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--source-root",
        type=Path,
        default=REPO.parent,
        help="workspace containing approach_147, approach_154, approach_155, and approach_156",
    )
    args = parser.parse_args()
    source_root = args.source_root.resolve()
    records_root = REPO / "records"
    records_root.mkdir(exist_ok=True)

    records = (
        first_records(source_root)
        + approach_155_records(source_root)
        + approach_156_records(source_root)
    )
    records.sort(key=lambda row: row["mu"])
    used = set()
    packaged = []
    for record in records:
        slug = make_slug(
            record["mu"],
            record["architecture"],
            record["n"],
            record["degree"],
            used,
        )
        packaged.append(
            package_record(record, records_root / slug, source_root, slug)
        )
    write_index(packaged)
    print(f"packaged {len(packaged)} finite witnesses")


if __name__ == "__main__":
    main()
