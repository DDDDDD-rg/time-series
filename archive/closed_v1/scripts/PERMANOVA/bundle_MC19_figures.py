"""Combine existing MC19 figure pages without rerunning statistics or plotting."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import tempfile

from pypdf import PdfReader, PdfWriter

ENVIRONMENT = ("Nitrate_mgN_L", "Nitrite_mgN_L", "Ammonia_mgN_L", "DIN_mgN_L", "DIP_mgP_L")
GROUPS = {
    "00_1_总体结果与NMDS.pdf": ["environment_R2_BH", "NMDS_TS", *["NMDS_environment_" + e for e in ENVIRONMENT], "host_highlight_descriptive"],
    "00_2_环境与MC时间变化.pdf": ["environment_time_series", "MC_individual_samples_by_TS", "MC_sample_abundance_heatmap"],
    "00_3_环境梯度下MC组成.pdf": ["MC_composition_" + e for e in ENVIRONMENT],
}


def source_pdf(figure_root: Path, stem: str) -> Path:
    """Prefer the organized layout; read legacy flat directories as a fallback."""
    organized = figure_root / "supporting_files" / (stem + ".pdf")
    return organized if organized.is_file() else figure_root / (stem + ".pdf")


def bundle(figure_root: Path, allow_partial: bool = False, environments=ENVIRONMENT) -> dict:
    figure_root = figure_root.resolve(strict=True)
    environments = tuple(environments)
    if not environments or len(set(environments)) != len(environments) or not set(environments) <= set(ENVIRONMENT) | {"SST_C"}:
        raise ValueError("Invalid environmental panel list")
    groups = {name: list(stems) for name, stems in GROUPS.items()}
    for name, stems in groups.items():
        if any(stem.startswith("NMDS_environment_") for stem in stems):
            groups[name] = ["environment_R2_BH", "NMDS_TS", *["NMDS_environment_" + e for e in environments], "host_highlight_descriptive"]
        elif any(stem.startswith("MC_composition_") for stem in stems):
            groups[name] = ["MC_composition_" + e for e in environments]
    sources = {name: [source_pdf(figure_root, stem) for stem in stems] for name, stems in groups.items()}
    missing = {name: [str(p) for p in paths if not p.is_file()] for name, paths in sources.items()}
    missing = {name: paths for name, paths in missing.items() if paths}
    if missing and not allow_partial:
        raise ValueError("Missing PDF sources: " + json.dumps(missing, ensure_ascii=False))
    staged = []
    try:
        for name, paths in sources.items():
            if name in missing:
                continue
            writer = PdfWriter()
            pages = 0
            for path in paths:
                reader = PdfReader(path)
                if reader.is_encrypted:
                    raise ValueError("Encrypted PDF source: " + str(path))
                pages += len(reader.pages)
                writer.append(reader, outline_item=path.stem.replace("_", " "))
            handle, filename = tempfile.mkstemp(prefix=".mc19-bundle-", suffix=".pdf", dir=figure_root)
            temporary = Path(filename)
            staged.append({"name": name, "temporary": temporary, "pages": pages})
            with os.fdopen(handle, "wb") as stream:
                writer.write(stream)
                stream.flush()
                os.fsync(stream.fileno())
            writer.close()
            with temporary.open("rb") as stream:
                if len(PdfReader(stream).pages) != pages:
                    raise ValueError("Unexpected combined page count: " + name)
        # Validate every requested group before replacing any existing combined PDF.
        for item in staged:
            os.replace(item["temporary"], figure_root / item["name"])
        return {"figure_root": str(figure_root), "bundles": [{"file": i["name"], "pages": i["pages"]} for i in staged],
                "skipped_incomplete_groups": list(missing)}
    finally:
        for item in staged:
            item["temporary"].unlink(missing_ok=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--figure-dir", type=Path, required=True)
    parser.add_argument("--allow-partial", action="store_true", help="Refresh complete groups when only a plotting branch was requested")
    parser.add_argument("--environments", nargs="+", default=ENVIRONMENT)
    args = parser.parse_args()
    print(json.dumps(bundle(args.figure_dir, args.allow_partial, args.environments), ensure_ascii=False, indent=2))
