# 04_multipanel_roc synthetic validation

- Random generator: NumPy PCG64 seed 20260921, synthetic only.
- Intent: 多指标 ROC 判别 多面板
- Selected FigureYa module: `FigureYa102multipanelROC`
- Fixed upstream commit: `f627917b79f28558779fb2e3ea2014cede41b89e`
- Top score (retrieval only): 93; candidates compared: 8
- Candidate comparison: FigureYa102multipanelROC (93, data=structural_candidate_only); FigureYa24ROC (73, data=structural_candidate_only); FigureYa200pairwiseAUC (70, data=structural_candidate_only)
- Data check: `structural_candidate_only`; blockers: []
- Adaptation: Python rendering of the selected module’s plotting semantics; upstream Rmd fetched and hash checked, not executed because Rscript is unavailable.
- Source reference: https://github.com/ying-ge/FigureYa/tree/f627917b79f28558779fb2e3ea2014cede41b89e/FigureYa102multipanelROC
