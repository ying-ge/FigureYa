# FigureYa 绘图智能体：本工程入口

本工程用于根据用户数据与表达目标，从远程 FigureYa 选择真实适配的基座、修改脚本并交付可复现图件。默认中文沟通。本文件所在目录就是项目根目录；所有命令从根目录运行。

## 必须保留的用户约束

- **禁止完整 clone、镜像、仓库 ZIP、批量模块 ZIP、稀疏/浅克隆或用逐文件下载绕过限制。** 只保留轻量索引；选定具体任务之后才按需取个别源码。参考网站中的克隆建议不适用本工程。
- **原始 FigureYa9heatmap 范例已清理。若用户之后放入新的单图范例，也不能把它当成本智能体的知识来源、默认模板或测试依据。** 用户数据不能被自动当作示例分析。
- 知识入口是用户提供的三个网址：[技能库](https://medaibox.com/ai-skills/)、[FigureYa教程](https://medaibox.com/ai-skills/figureya.html)、[上游仓库](https://github.com/ying-ge/FigureYa)。教程提供工作方式，实际路径、输入和功能以固定版本的上游源码为准。
- 不按最容易写、已经在本地、以前用过或下载最小决定图型。语义和统计设计不匹配时，应更换基座或如实说明库的缺口。

## 项目技能（直接读取，不依赖是否出现在全局技能列表）

1. 数据分析、找模板、比较图型、开始绘图：读 [.agents/skills/figureya-select/SKILL.md](.agents/skills/figureya-select/SKILL.md)。
2. 选定基座后生成/改写/修复 R、Rmd、Python：读 [.agents/skills/figureya-adapt/SKILL.md](.agents/skills/figureya-adapt/SKILL.md)。
3. 出图、修版、最终交付：读 [.agents/skills/figureya-qa/SKILL.md](.agents/skills/figureya-qa/SKILL.md)。

技能位于项目 `.agents/skills`，并通过此入口可靠路由；不需要把它们复制到全局目录或配置外部模型服务。若桌面技能列表尚未刷新，仍能按上述路径使用。

## 最短操作链

```powershell
python figureya-agent/scripts/figureya.py profile "用户数据.csv" --out "outputs/任务名/profile.json"
python figureya-agent/scripts/figureya.py recommend --intent "实际表达目标" --profile "outputs/任务名/profile.json" --out "outputs/任务名/candidates.json"
python figureya-agent/scripts/figureya.py inspect "候选完整模块ID"
```

检索器只召回候选，智能体负责复核：数据意义、观测单位、比较设计、输入契约、规模与可读性；选择后才能按 `fetch` 的精确路径下载个别源码。不要把 retrieval_score 当概率或最终结论。

索引在 `figureya-agent/catalog/`；命令说明和限制在 [figureya-agent/README.md](figureya-agent/README.md)。为每次真实任务创建 `outputs/<任务名>/`，保留原始数据、选型理由、改写脚本、来源和验证记录。用户没有提供真实数据时，交付智能体本身，不把仓库示例当作用户数据。
