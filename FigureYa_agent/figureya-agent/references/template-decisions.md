# 模板辨析：依据输入和表达目标选择

完整目录见 `catalog/INDEX.md`（相对于figureya-agent目录）。下面是经过远程说明复核的区别，不是“推荐永远用这些模板”。精确卡片与提示词在 `catalog/curated.json`；实际源码仍需在选定时读取。

## 容易因名称选错的分支

| 目标或数据特征 | 应比较的完整模块 ID | 决定性的区别 |
|---|---|---|
| 一般组间连续分布 | FigureYa12box / FigureYa162boxViolin / FigureYa227boxdensity / FigureYa298ecdfPvalue | 是否需要原始点、密度还是累计分布；小样本不靠平滑隐藏观测 |
| 三个特征的分面分布 | FigureYa254scViolin / FigureYa162boxViolin | 254也接受普通长表，不能只看sc名称排除 |
| 配对前后变化 | FigureYa195PanPaire / 162的改写绘图层 | 195的原流程含TCGA和评分，普通配对数据只借用连线层；不能沿用独立样本检验 |
| 雨云外观 | FigureYa68friendsV2 / FigureYa227boxdensity / 162绘图层 | 68原流程是GO语义相似，不应为外观强制用户提交基因列表 |
| 两组差异火山 | FigureYa59volcanoV2 / FigureYa75bubble_volcano | ggplot标签/高亮与base多属性编码的区别，按真实变量数量选择 |
| 火山附边缘类别比例 | FigureYa135multiVolcano / 59的分面扩展 | 135的multi是多层信息，不等同多个contrast；不要照搬教程简写 |
| APA差异 | FigureYa321volcanoSE / 59的绘图层 | 321输入PDUI和组别，其横轴不应误标log2FC |
| 已有富集表 | FigureYa11bubble / FigureYa5bubbles / FigureYa75bubble_volcano | Term/P/Count/GeneRatio与表达/显著性的组合编码，不必重跑富集 |
| GSEA曲线 | FigureYa60GSEA_clusterProfilerV2 / FigureYa136fgsea / FigureYa13GSEA_Java_update | 必须有排序和基因集；工具、分析环境和数据阶段不同 |
| 样本相关性 | FigureYa9heatmap / FigureYa37correlationV2_update / FigureYa76corrgram | 相关对象是样本还是变量；不复用表达行Z-score色标 |
| 两组变量间相关 | FigureYa97correlationV3 / FigureYa181multiCorrelation / FigureYa228linkCor | 样本ID对齐；矩阵、散点或矩阵+连线的表达需求 |
| 两种测量方法一致性 | FigureYa176BlandAltman（与相关图对照） | 相关不回答可替换性；只有一个合适基座时不要凑三个 |
| 指定顺序/导出聚类 | FigureYa91cluster_heatmap / 9 / FigureYa243scMarkerGroupHeatmap | 是否允许聚类改变生物学顺序，是否按marker模块分块 |
| 两个对齐指标矩阵 | FigureYa144DiagHeatmap / FigureYa278heatmapPoints | 两半色块或色块+气泡；两套色标和数据含义明确 |
| 效应量矩阵+P值矩阵 | FigureYa149rankHeatmap / FigureYa165heatmapPvalue | 149输入已有矩阵可分档；165从表达和分组做检验，不能混为同一输入 |
| 分类yes/no关联 | FigureYa124AssociationHeatmap / FigureYa125FishertestV2 | 类别输入，不能套连续相关检验 |
| 目标基因相关的表达模式 | FigureYa126CorrelationHeatmap / FigureYa201ClusterCorrelation | 表达热图与相关系数热图不同；共识聚类还有额外分析成本 |
| 药物多证据热图 | FigureYa213customizeHeatmap / FigureYa212drugTargetV2 | 带CMap/药敏等前置流程；“customize”不保证通用热图优先 |
| CNV与共享克隆 | FigureYa307CNVHeatmap / FigureYa320ClontypeHeatmap | GISTIC状态与克隆事件计数是不同数据，不能以相同热图外观替换 |
| PCA分组/批次/边界 | FigureYa38PCA / FigureYa101PCA / FigureYa244PCAPlot | 38椭圆、101双层元数据、244凸包/已有坐标；分离外观不是统计证明 |
| 普通ROC、多指标ROC、生存ROC | FigureYa24ROC / FigureYa102multipanelROC / FigureYa85timeROC | 二分类与删失结局的区别先于多面板样式 |
| 生存分组展示 | FigureYa1survivalCurve_update / FigureYa36nSurvV3 / FigureYa284pairwiseLogrank | 事件与时间、组数及对比；不能默认最佳截点优化 |
| 森林图 | FigureYa47HRtable / FigureYa90subgroup / FigureYa216MetaREM | 表格排版、亚组模型、Meta分析是不同数据阶段；有现成CI优先复用绘图层 |
| 多层属性 | FigureYa25Sankey_update / FigureYa174squareCross / FigureYa259circLink | Alluvial的流量和阶段、任意节点边、弦图的语义不同 |
| 集合关系 | FigureYa112venn / FigureYa112Plus_venn / FigureYa237circVenn | 集合数、是否按区域比例填色、是否真正保留完整交集；必要时自编UpSet |

## 如何避免固定偏好

候选比较记录“为何最合适”也记录“什么变化会改选”：如新增批次字段会倾向101，已有坐标和凸包需求会倾向244；提供双矩阵会考虑144/278；研究问题变为一致性会离开相关图族。相同数据不同问题应可能选不同图；相同问题不同数据结构也应可能选不同基座。不要为了展示花样强行用复杂模块。

部分上游说明带有不严谨的泛化，例如PCA二维重叠即可判断聚类过拟合、固定小样本阈值或把相关性叫影响。保留可用绘图代码，独立审查统计解释。不要让教程宣传语代替对当前数据的验证。
