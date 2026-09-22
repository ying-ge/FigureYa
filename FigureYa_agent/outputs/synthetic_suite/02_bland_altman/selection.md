# 02_bland_altman synthetic validation

- Random generator: NumPy PCG64 seed 20260921, synthetic only.
- Intent: 测量方法一致性 Bland Altman agreement
- Selected FigureYa module: `FigureYa176BlandAltman`
- Fixed upstream commit: `f627917b79f28558779fb2e3ea2014cede41b89e`
- Top score (retrieval only): 99; candidates compared: 4
- Candidate comparison: FigureYa176BlandAltman (99, data=structural_candidate_only); FigureYa138NiceCalibration (8, data=blocked_or_unknown); FigureYa236circGroup (8, data=blocked_or_unknown)
- Data check: `structural_candidate_only`; blockers: []
- Adaptation: Python rendering of the selected module’s plotting semantics; upstream Rmd fetched and hash checked, not executed because Rscript is unavailable.
- Source reference: https://github.com/ying-ge/FigureYa/tree/f627917b79f28558779fb2e3ea2014cede41b89e/FigureYa176BlandAltman
