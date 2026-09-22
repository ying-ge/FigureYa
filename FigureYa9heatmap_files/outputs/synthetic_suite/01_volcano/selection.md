# 01_volcano synthetic validation

- Random generator: NumPy PCG64 seed 20260921, synthetic only.
- Intent: ggplot2 火山图 基因标签 候选基因
- Selected FigureYa module: `FigureYa59volcanoV2`
- Fixed upstream commit: `f627917b79f28558779fb2e3ea2014cede41b89e`
- Top score (retrieval only): 96; candidates compared: 8
- Candidate comparison: FigureYa59volcanoV2 (96, data=structural_candidate_only); FigureYa135multiVolcano (76, data=structural_candidate_only); FigureYa321volcanoSE (76, data=structural_candidate_only)
- Data check: `structural_candidate_only`; blockers: []
- Adaptation: Python rendering of the selected module’s plotting semantics; upstream Rmd fetched and hash checked, not executed because Rscript is unavailable.
- Source reference: https://github.com/ying-ge/FigureYa/tree/f627917b79f28558779fb2e3ea2014cede41b89e/FigureYa59volcanoV2
