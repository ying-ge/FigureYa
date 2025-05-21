# 最新版Seurat安装

保证电脑上已经安装了最新的R(3.6)和Rstudio

安装Anaconda，

- 官网地址：https://www.anaconda.com/distribution/#download-section

- 网盘地址: 链接：https://share.weiyun.com/5oy1F2b 密码：b8jlhc

- All Users

- Advanced Options: 将Anaconda加入到系统环境变量中(Add Anaconda to the system PATHenvironment variable)

在Anaconda上安装[UMAP](https://github.com/lmcinnes/umap)

```shell
pip install -i https://pypi.tuna.tsinghua.edu.cn/simple umap-learn
```

测试

```python
import umap
from sklearn.datasets import load_digits

digits = load_digits()

embedding = umap.UMAP().fit_transform(digits.data)
```

打开Rstudio 安装Seurat

```r
install.packages("BiocManager")
chooseBioCmirror()
BiocManager::install("Seurat")
```

测试UMAP

```r
library(Seurat)
pbmc_small <- RunUMAP(object = pbmc_small, dims = 1:5)
DimPlot(pbmc_small)
```


注意事项: 

不要在R里面装Python的包

保证环境变量里有如下位置

- /path/to/Anaconda3
- /path/to/Anaconda3/Library/mingw-w64/bin
- /path/to/Anaconda3/Library/usr/bin
- /path/to/Anaconda3/Library/bin
- /path/to/Anaconda3/Scripts
