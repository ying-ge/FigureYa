# 03_pca_batch synthetic validation

- Random generator: NumPy PCG64 seed 20260921, synthetic only.
- Intent: PCA 分组 批次 batch 形状
- Selected FigureYa module: `FigureYa101PCA`
- Fixed upstream commit: `f627917b79f28558779fb2e3ea2014cede41b89e`
- Top score (retrieval only): 105; candidates compared: 8
- Candidate comparison: FigureYa101PCA (105, data=structural_candidate_only); FigureYa164PCA3D (79, data=structural_candidate_only); FigureYa38PCA (79, data=structural_candidate_only)
- Data check: `structural_candidate_only`; blockers: []
- Adaptation: Python rendering of the selected module’s plotting semantics; upstream Rmd fetched and hash checked, not executed because Rscript is unavailable.
- Source reference: https://github.com/ying-ge/FigureYa/tree/f627917b79f28558779fb2e3ea2014cede41b89e/FigureYa101PCA
