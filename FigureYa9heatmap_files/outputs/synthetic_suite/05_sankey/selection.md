# 05_sankey synthetic validation

- Random generator: NumPy PCG64 seed 20260921, synthetic only.
- Intent: 桑基图 多层分类 流向 alluvial
- Selected FigureYa module: `FigureYa25Sankey_update`
- Fixed upstream commit: `f627917b79f28558779fb2e3ea2014cede41b89e`
- Top score (retrieval only): 93; candidates compared: 4
- Candidate comparison: FigureYa25Sankey_update (93, data=structural_candidate_only); FigureYa40lineage (70, data=structural_candidate_only); FigureYa236circGroup (20, data=blocked_or_unknown)
- Data check: `structural_candidate_only`; blockers: []
- Adaptation: Python rendering of the selected module’s plotting semantics; upstream Rmd fetched and hash checked, not executed because Rscript is unavailable.
- Source reference: https://github.com/ying-ge/FigureYa/tree/f627917b79f28558779fb2e3ea2014cede41b89e/FigureYa25Sankey_update
