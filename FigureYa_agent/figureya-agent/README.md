# FigureYa 专用绘图智能体

在这个工程中向 Codex 描述数据与绘图目标即可。根目录 `AGENTS.md` 会引导它使用项目内三个技能：数据与模板选型、脚本适配、图件验收。不需要部署服务器、配置额外模型API或克隆FigureYa。

**原始 FigureYa9heatmap 单图范例已从工程中清理，也从未作为本智能体的参考基础。** 索引来自用户指定网站及整个远程仓库。这里保存的是知识索引与原创辅助工具，上游源码只在将来的具体绘图任务中按需获取。

## 使用示例

- “读取我的差异结果表，比较合适的FigureYa基座，突出这些候选基因，生成PDF和预览。”
- “这是一组患者治疗前后的配对测量，先判断合适的图，不要按独立组检验。”
- “我有两个对齐的矩阵，想同时表达效应量和显著性，比较不同编码方式后出图。”
- “这是样本表达矩阵和批次表，请同时展示分组和批次，并解释基座选择。”

不必记模板编号；若当前技能列表没有显示新技能，直接说“按本工程AGENTS.md使用FigureYa绘图智能体”。新任务在本目录启动即可加载入口。

## 索引内容与可信度

- [INDEX.md](catalog/INDEX.md)：全模块中文分组导航及固定版本源码链接。
- [modules.json](catalog/modules.json)：远程Rmd的场景、输入说明、依赖、字面文件引用、绘图函数和审查提示；**不是完整源码**。
- [curated.json](catalog/curated.json)：人工辨析卡片，特别标出名称误导与前置分析。
- [routing.json](catalog/routing.json)：用途、输入契约与候选检索分组。
- [upstream-tree.json](catalog/upstream-tree.json)：仅路径、大小、Git blob哈希等元数据，用来核验文件存在和按需下载。

快照数量、构建时间、commit见modules.json；不是实时镜像。`source_metadata_extracted`表示自动提取，**不代表324份脚本都已执行或所有统计方法均被认可**。没有Rmd的应用/其他语言模块标记为inventory_only，使用时另读源码。检索分数只作召回，最终选择由智能体根据数据和脚本比较完成。

## 本地命令

在工程根目录使用Python 3.10或以上；CSV/TSV工具只用标准库。XLSX额外需要openpyxl。绘图时才需要相应R/Python包。若PATH中没有Python，使用Codex提供的依赖运行时。

```powershell
python figureya-agent/scripts/figureya.py profile "输入.csv" --out "outputs/demo/profile.json"
python figureya-agent/scripts/figureya.py profile "输入.tsv" --delimiter tab --encoding utf-8-sig
python figureya-agent/scripts/figureya.py profile "输入.xlsx" --sheet "数据表"
python figureya-agent/scripts/figureya.py recommend --intent "配对测量分布和连线" --profile "outputs/demo/profile.json" --design paired --out "outputs/demo/candidates.json"
python figureya-agent/scripts/figureya.py inspect FigureYa176BlandAltman
python figureya-agent/scripts/figureya.py audit "已选脚本.Rmd"
```

非标准列名可提供`--roles roles.json`。`profile`默认只读取前10000行，可用`--limit`调整至最多100000行，报告明确显示是否抽样。编码错误须按数据实际编码指定（如`--encoding gb18030`）；不静默丢字符。

只有选定具体任务的基座后，才精确取文件：

```powershell
python figureya-agent/scripts/figureya.py fetch FigureYa176BlandAltman --path FigureYa176BlandAltman/FigureYa176BlandAltman.Rmd
```

下载固定commit文件并校验Git blob哈希，默认每次5MiB、全缓存32MiB；禁止存档/数据/HTML和递归获取；不会自动执行任何下载的脚本。辅助源码也必须明确给出路径。需要运行的源码复制到任务目录后改写，不修改缓存原件。

## 更新与验证

仅在确实需要更新时执行：

```powershell
python figureya-agent/scripts/build_catalog.py --refresh
python figureya-agent/scripts/figureya.py export-index
python -m unittest discover -s figureya-agent/tests -v
```

刷新读取目录元数据并在线浏览Rmd，只将提取字段写入索引；不会保存脚本、报告、图片和数据。请求量为一次目录查询加数百个小文本请求，通常约几MB传输；不需要每次绘图都刷新。API限流/离线时使用已有索引，取新源码失败要如实报告。更新后检查路径变化和人工卡片是否过期。

本工程验证记录见 [VALIDATION.md](VALIDATION.md)。这次交付的是由Codex执行的绘图智能体技能和工具，不是已经对用户真实数据完成出图，也不是独立运行的模型服务。
