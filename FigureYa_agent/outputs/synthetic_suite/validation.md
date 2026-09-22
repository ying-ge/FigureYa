# Synthetic FigureYa agent validation

Status: PASS for template routing, exact source retrieval, synthetic rendering and artifact existence.

- Inputs are seeded synthetic data and are not the local FigureYa9heatmap example.
- Each task profiled its own input, compared candidates, selected an expected semantically specific module, fetched one exact Rmd at the pinned upstream commit, and rendered PDF+PNG.
- Numeric checks: volcano significant count, Bland–Altman bias/limits, PCA explained variance, ROC AUC, Sankey path aggregation are recorded in summary.json.
- Visual inspection PASS: all five generated PNGs opened at high resolution. Volcano labels and threshold lines are inside the canvas; Bland–Altman bias/limits and axes are readable; PCA group colors and batch markers have a readable legend; ROC panels share axes and show AUC; Sankey stages, labels and ribbons are visible without clipping.
- Not run: upstream Rmd execution, because Rscript is unavailable. This is an explicit limitation, not a pass.
