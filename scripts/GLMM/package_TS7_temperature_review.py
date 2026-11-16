"""Package exported saved-result figures; never fit models or change a deck."""
from pathlib import Path
import csv
import hashlib
import json
import argparse

from PIL import Image
from pypdf import PdfReader, PdfWriter
from reportlab.pdfgen import canvas
from reportlab.lib.utils import ImageReader


def digest(path):
    return hashlib.md5(Path(path).read_bytes()).hexdigest()


def images_pdf(images, target):
    out = canvas.Canvas(str(target))
    out.setTitle(target.stem)
    for path in images:
        with Image.open(path) as im:
            width = 12.7 * 72
            height = width * im.height / im.width
            out.setPageSize((width, height))
            out.drawImage(ImageReader(im), 0, 0, width=width, height=height)
            out.showPage()
    out.save()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parents[2] / "figures" / "temperature_review_v1")
    args = parser.parse_args()
    root = args.output.resolve()
    support = root / "supporting_files"
    assert (support / "export_status.txt").read_text().startswith("complete;")
    with (support / "protected_result_and_PPT_md5.csv").open(encoding="utf-8-sig", newline="") as f:
        protected = list(csv.DictReader(f))
    for row in protected:
        assert digest(row["path"]) == row["md5"], row["path"]
    with (support / "unchanged_original_figures.csv").open(encoding="utf-8-sig", newline="") as f:
        originals = list(csv.DictReader(f))
    for row in originals:
        assert digest(row["source"]) == digest(row["copy"]) == row["md5"]
    main_files = [
        ("02_GLMM_heatmap_forest.pdf", "Five-indicator GLMM overview"),
        ("03_GLMM_nutrients_PPT_original.pdf", "Original PPT nutrient GLMM curves"),
        ("04_GLMM_temperature_curves.pdf", "Temperature GLMM curves [R]"),
        ("05_environment_by_TS.pdf", "Four nutrients and temperature by TS"),
        ("06_two_part_nutrients_PPT_original.pdf", "Original PPT nutrient two-part forest"),
        ("08_two_part_temperature_effects.pdf", "Temperature two-part effects [R]"),
    ]
    # Component curves and diagnostics stay as separate supporting figures.
    supporting_only = ["07_two_part_nutrients_curves_original.pdf", "09_two_part_temperature_curves.pdf",
                       "10_temperature_component_residuals.pdf", "11_temperature_model_review.pdf"]
    merge_sources = []
    writer = PdfWriter()
    for scope in ("coral", "ecosystem"):
        folder = root / scope
        pages = sorted(folder.glob("03_GLMM_nutrients_PPT_page_*.png"))
        assert len(pages) == (3 if scope == "coral" else 4)
        if not (folder / "03_GLMM_nutrients_PPT_original.pdf").exists():
            images_pdf(pages, folder / "03_GLMM_nutrients_PPT_original.pdf")
        if not (folder / "06_two_part_nutrients_PPT_original.pdf").exists():
            images_pdf([folder / "06_two_part_nutrients_PPT_original.png"], folder / "06_two_part_nutrients_PPT_original.pdf")
        scope_bookmark = writer.add_outline_item(scope.title(), len(writer.pages))
        for filename, label in main_files:
            page_number = len(writer.pages)
            writer.append(str(folder / filename), import_outline=False)
            writer.add_outline_item(label, page_number, parent=scope_bookmark)
            merge_sources.append({"scope": scope, "source": filename, "first_page": page_number + 1,
                                  "last_page": len(writer.pages)})
    assert len(writer.pages) == 18, len(writer.pages)
    writer.add_metadata({"/Title": "Seven-TS temperature figure review", "/Subject": "Saved results only; PPT style; no PERMANOVA or NMDS"})
    with (root / "00_main_figures.pdf").open("wb") as f:
        writer.write(f)
    writer.close()
    readme = """# 七TS温度补图检查包

## 查看方式

先看 `00_main_figures.pdf`，已精简为18页，按 coral、ecosystem 分组并加书签。主文件只合并02、03、04、05、06、08六类主要结果图。two-part逐项曲线（07、09）、组件残差（10）及模型诊断（11）单独保留，不放进主合并版。

## 图形和结果口径

- 样式依据实际PPT `presentations/time_series_2_text_v6.pptx`。新热图、森林图和温度曲线沿用其颜色、字号、线型及分离的放大框；环境时序图加入温度。旧营养盐GLMM曲线和two-part森林图直接复制PPT内嵌图像，字节一致；旧two-part曲线PDF保留。
- 五指标总览使用当前GLMM的BH q（每范围95项）。03的PPT原曲线页保留旧76项族的q，显著组合没有改变。温度新增：coral MC25正向；ecosystem MC25正向、MC2s1负向。
- 温度two-part使用当前18项（coral）/32项（ecosystem）统一BH：coral MC25两组件均未达0.05；ecosystem MC25两组件均显著，MC2s1仅检出组件显著。
- 为保留原来的营养盐结论和PPT图，营养盐two-part历史页仍保留旧16/28项检验族的q和显著标记。当前完整结果中，coral MC98—总氨氮检出组件q由旧0.04817变为0.05419；该历史页的星号未改。温度补充页单独列出，避免把历史页当成更新后的统一BH结果。
- `[R]`仅标注温度父Beta-binomial GLMM的诊断复核；不在营养盐上加review标记。该标记不表示two-part组件另做了同样的正式诊断。全部19个温度MC诊断摘要及显著组合数据见supporting_files。

## 内容与文件

02：五指标GLMM热图＋显著指标森林图；03：PPT原营养盐GLMM曲线；04：显著温度MC曲线；05：四营养盐＋温度时序；06：PPT原营养盐two-part森林图；07：原营养盐two-part曲线；08：温度two-part森林图；09：温度two-part曲线；10：温度组件残差；11：显著温度MC的原保存诊断图。

未生成PERMANOVA、NMDS和整体组成相关图。未重新拟合模型、模拟或置换；未修改原结果、原图、PPT或联合GLMM。MD5完整性核对、图源和范围检查保存在supporting_files。2026-11-16按用户要求将two-part逐项辅助曲线从主合并版移出，共减少25页；原独立文件保留。

重现：先运行 `scripts/GLMM/export_TS7_temperature_review.R`，再运行同目录 `package_TS7_temperature_review.py`。R入口默认拒绝覆盖已有检查包；明确重绘时可设置 `options(mc19.review.redraw=TRUE)`，只写检查包。需要本地已有结果、PPT图源快照与compact绘图模板。
"""
    (root / "00_README.md").write_text(readme, encoding="utf-8")
    with (support / "main_merge_sources.csv").open("w", encoding="utf-8-sig", newline="") as f:
        out = csv.DictWriter(f, fieldnames=merge_sources[0].keys())
        out.writeheader()
        out.writerows(merge_sources)
    inventory = []
    for path in sorted(root.rglob("*.pdf")):
        reader = PdfReader(path)
        assert len(reader.pages) > 0
        inventory.append({"path": str(path.relative_to(root)), "pages": len(reader.pages), "bytes": path.stat().st_size, "md5": digest(path)})
    assert len(inventory) == 21, len(inventory)
    with (support / "PDF_inventory.csv").open("w", encoding="utf-8-sig", newline="") as f:
        out = csv.DictWriter(f, fieldnames=inventory[0].keys())
        out.writeheader()
        out.writerows(inventory)
    for row in protected:
        assert digest(row["path"]) == row["md5"]
    report = {"pdf_count": len(inventory), "main_pdf_pages": len(PdfReader(root / "00_main_figures.pdf").pages),
              "protected_result_and_PPT_files": len(protected), "protected_hashes_unchanged": True,
              "old_figure_copies_verified": len(originals), "statistical_analysis_run": False,
              "PPT_changed": False, "PERMANOVA_or_NMDS_exported": False,
              "supporting_only_excluded_from_main": supporting_only}
    (support / "validation.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report))


if __name__ == "__main__":
    main()
