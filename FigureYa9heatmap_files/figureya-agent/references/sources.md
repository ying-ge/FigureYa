# 来源与研究范围

1. 用户指定的 [AI科研技能库](https://medaibox.com/ai-skills/)：理解技能/Agent/工具的入口关系，仅阅读与本任务有关内容。
2. 用户指定的 [FigureYa专题](https://medaibox.com/ai-skills/figureya.html)：学习外部代码库、标准Rmd、输入、参考图与矢量输出的流程。其200+为页面表述；具体数量采用实际仓库快照。网页的整库clone示例明确不采用。
3. 用户指定的 [FigureYa GitHub](https://github.com/ying-ge/FigureYa)：README、目录树、全部Rmd在内存中浏览提取字段，部分易混淆基座额外人工复核关键输入与正文。

当前索引的固定commit与获取时间在 `catalog/modules.json`（相对于figureya-agent）。GitHub tree读取完整（truncated=false），目录文件只存路径/大小/哈希。所有Rmd内容仅用于在线学习及元数据提取，没有保存上游源码副本。

额外阅读了仓库README_AI_RAG.md、docs/For_clinical.md与chapters.json的内容；没有采用其另行部署云服务的路线。发现说明页可能存在过时模块名和宽泛应用描述，因此所有可获取文件必须以真实tree和commit为准；展示推荐不等于接受其医学/统计推论。

上游README声明 CC BY-NC-SA 4.0，并提供论文引用：Lu X, et al. FigureYa: A Standardized Visualization Framework for Enhancing Biomedical Data Interpretation and Research Efficiency. iMetaMed (2025), e70005. [DOI](https://doi.org/10.1002/imm3.70005)。本地自动提取摘要的来源链接随模块保存；使用上游代码时继续保留归属与许可。

本次不读取当前目录热图示例作为学习来源，也不围绕它选择模板或验证。测试数据均为独立构造的合成表，且不代表真实科研结论。
