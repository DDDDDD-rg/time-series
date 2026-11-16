"""Assemble the six-indicator primary figure pack, keeping auxiliary checks separate."""
from pathlib import Path
import argparse
import csv
import hashlib
import json
from pypdf import PdfReader, PdfWriter


def md5(path):
    return hashlib.md5(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path(__file__).resolve().parents[2] / 'figures/TS7_env6_PPT_v1')
    root = parser.parse_args().output.resolve()
    support = root / 'supporting_files'
    assert (support / 'export_status.txt').read_text().startswith('complete;')
    with (support / 'protected_result_PPT_hashes.csv').open(encoding='utf-8-sig') as f:
        protected = list(csv.DictReader(f))
    protected = [r for r in protected if not Path(r['path']).name.startswith('~$')]
    assert all(md5(r['path']) == r['md5'] for r in protected)
    writer = PdfWriter()
    sources = []
    for scope in ['coral', 'ecosystem']:
        folder = root / scope
        curve_pages = sorted(folder.glob('03_GLMM_curves_page_*.pdf'))
        combined = folder / '03_GLMM_significant_curves.pdf'
        if curve_pages:
            curves = PdfWriter()
            for path in curve_pages:
                curves.append(path, import_outline=False)
            assert len(curves.pages) == (3 if scope == 'coral' else 5)
            with combined.open('wb') as f:
                curves.write(f)
            curves.close()
            # These are staging PDFs authored by the exporter; per-page PNGs remain.
            for path in curve_pages:
                path.unlink()
        files = [
            '01_environment_by_TS.pdf', '06_NMDS2D_nitrite_ammonia.pdf',
            '05_NMDS2D_Nitrate_mgN_L.pdf', '05_NMDS2D_DIN_mgN_L.pdf',
            '05_NMDS2D_DIP_mgP_L.pdf', '05_NMDS2D_SST_C.pdf',
            '02_GLMM_heatmap_forest.pdf', '03_GLMM_significant_curves.pdf',
            '04_two_part_forest.pdf', '07_PERMANOVA_summary.pdf',
        ]
        parent = writer.add_outline_item(scope.title(), len(writer.pages))
        for name in files:
            start = len(writer.pages)
            writer.append(folder / name, import_outline=False)
            writer.add_outline_item(name.removesuffix('.pdf'), start, parent=parent)
            sources.append(dict(scope=scope, source=name, first_page=start + 1, last_page=len(writer.pages)))
    assert len(writer.pages) == 26
    writer.add_metadata({'/Title': 'Seven TS - six indicators - current result figures',
                         '/Subject': 'PPT graphic style; auxiliary checks excluded; current BH values'})
    with (root / '00_main_figures.pdf').open('wb') as f:
        writer.write(f)
    writer.close()
    with (support / 'main_merge_sources.csv').open('w', encoding='utf-8-sig', newline='') as f:
        out = csv.DictWriter(f, fieldnames=sources[0].keys()); out.writeheader(); out.writerows(sources)
    inventory = []
    for p in sorted(root.rglob('*.pdf')):
        pages = len(PdfReader(p).pages)
        assert pages > 0, p
        inventory.append(dict(path=str(p.relative_to(root)), pages=pages, md5=md5(p)))
    assert len(inventory) == 37, len(inventory)
    with (support / 'PDF_inventory.csv').open('w', encoding='utf-8-sig', newline='') as f:
        out = csv.DictWriter(f, fieldnames=inventory[0].keys()); out.writeheader(); out.writerows(inventory)
    assert all(md5(r['path']) == r['md5'] for r in protected)
    readme = '''# 七TS六指标：PPT样式结果图

先看 `00_main_figures.pdf`（26页，按coral/ecosystem分组并带书签）。独立高清PNG/PDF在对应范围目录；辅助曲线、诊断和3D资料在各自 `auxiliary` 子目录，不并入主文件。

## 主要结果图

- 01：六指标环境时序；02：19MC×六指标热图及显著环境的森林图。
- 03：所有BH显著GLMM组合（coral10、ecosystem18），同一MC的指标同页同行；三指标的MC使用加宽页面，保留PPT的字体/面板尺寸。原始点与完整95%CI保留，局部放大框独立。
- 04：当前two-part森林图，coral20、ecosystem36组件BH。coral MC98—总氨氮检出不再标星，正值仍标#。
- 05：NO3、NO2、总氨氮、DIN、DIP、温度六套独立2D NMDS；06：与原PPT相同的NO2/总氨氮双面板；07：六项PERMANOVA解释率/P/q总览。

当前原始P<0.05的是coral NO2（.04762）、ecosystem NO2（.02857）及总氨氮（.03810）；所有环境BH q均>.05。图中保留raw P及BH q，不把原始P显著标成BH显著。

## 数据和样式

全部读取 `MC19_*_TS7_env6_*_run001` 新结果，不复制旧q。GLMM按114项BH，PERMANOVA按6项，two-part按20/36项。DIN与温度已进入所有相关结果图。

实际PPT `time_series_2_text_v6.pptx` 为样式依据，复用其compact绘图定义及配色/字体/线条；PPT本体保持不变。`[R]`继续只标温度父GLMM的诊断复核，其他环境的标记按此前要求暂不加；所有诊断数据/图另存。

NMDS矩阵、样本顺序和MC列已与原PPT的坐标来源逐项核对相同，因此复用原2D/3D坐标，更新当前六指标的着色和本轮PERMANOVA统计。每个范围六指标共用坐标，切换环境不会改变点的位置。未重算拟合、P/BH、残差模拟、置换或NMDS。

auxiliary内的3D HTML可切换六指标，请连同lib依赖目录一起保留。3D投影PDF、Shepard、two-part逐项辅助曲线和显著GLMM诊断均不进入主合并版。

supporting_files包含图源、保护文件MD5、NMDS复用核对、GLMM完整显示范围及主文件页码表。脚本为 `scripts/GLMM/export_TS7_env6_PPT_figures.R` 和 `package_TS7_env6_PPT_figures.py`；旧五指标检查包保留。
'''
    (root / '00_README.md').write_text(readme, encoding='utf-8')
    result = dict(pdf_count=len(inventory), main_pages=26, protected_files=len(protected),
                  protected_hashes_unchanged=True, auxiliary_in_main=False,
                  GLMM_significant_pairs=dict(coral=10, ecosystem=18), NMDS_2D_environment_plots=12,
                  refit_or_recalculation=False)
    (support / 'package_validation.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
