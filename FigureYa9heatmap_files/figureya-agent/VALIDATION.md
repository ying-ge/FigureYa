# FigureYa agent validation

The end-to-end synthetic suite is in `tests/run_demo_suite.py`; its outputs are in the workspace `outputs/synthetic_suite/` directory. It uses NumPy PCG64 seed `20260921` and does not read the removed single-figure example.

The run passed 5/5 semantic routing checks, 5/5 fixed-commit single-file source fetches with Git blob hashes, 5/5 adapted PDF/PNG renders, and 23/23 unit tests. The selected modules were:

| Synthetic task | Data contract / visual intent | Selected FigureYa module | Artifact |
|---|---|---|---|
| Differential expression | `gene`, `logFC`, `padj`, highlighted candidates | `FigureYa59volcanoV2` | `outputs/synthetic_suite/01_volcano/figureya59volcano_adapted.{pdf,png}` |
| Paired method agreement | two measurements per subject; bias and limits of agreement | `FigureYa176BlandAltman` | `outputs/synthetic_suite/02_bland_altman/figureya176_bland_altman_adapted.{pdf,png}` |
| Expression matrix + group + batch | PCA plus color/shape metadata | `FigureYa101PCA` | `outputs/synthetic_suite/03_pca_batch/figureya101_pca_batch_adapted.{pdf,png}` |
| Binary scores | two prediction scores and a binary truth label | `FigureYa102multipanelROC` | `outputs/synthetic_suite/04_multipanel_roc/figureya102_multipanel_roc_adapted.{pdf,png}` |
| Three categorical stages | observed multi-level paths | `FigureYa25Sankey_update` | `outputs/synthetic_suite/05_sankey/figureya25_sankey_adapted.{pdf,png}` |

The five PNGs were opened at high resolution. Labels, legends, axes, threshold/limit lines, ribbons, and panel boundaries were readable and inside the canvas. All five PDFs were independently rendered with Poppler `pdftoppm` successfully. `summary.json` records significant-point count, Bland–Altman bias and limits, PCA explained variance, ROC AUCs, and Sankey path count.

The source Rmd for each selected module was fetched only after selection and verified against the pinned upstream Git blob hash. The local environment has no `Rscript`, so these are Python adaptations of the selected FigureYa drawing semantics; upstream R execution and package-level equivalence remain `not-run`. This limitation is intentionally recorded rather than presented as an upstream execution pass.

Run the test suite with:

```powershell
$env:PYTHONIOENCODING='utf-8'
python figureya-agent/tests/run_demo_suite.py
python -m unittest discover -s figureya-agent/tests -v
```
