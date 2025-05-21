核心函数为如下3个函数：

模型验证模块
repeat_ncv(data,outer_fold,inter_fold,biresult_col,model_col,formula,n_repeat=100)：
作用：该函数的作用为对目标数据进行n_reapeat次的嵌套交叉验证
参数：data为包含建模特征与二分类结果的数据框，其中二分类结果的列名必须为“outcome”，二分类结果必须为0，1，列名一定要为数字，如1~n，n为样本数
           outer_fold为嵌套交叉验证的外层折叠数
           inter_fold为嵌套交叉验证的内层折叠数
           biresult_col为二分类结果所在列
           model_col为建模特征所在列
           formula为建模函数，如 结果~健康状况+血脂+血糖
           n_repeat为嵌套交叉验证的重复次数，默认为100次（文献中重复了100次）
核心返回值：mean_auc为 n_repeat次嵌套交叉验证中所有验证集（一共有outer_fold*n_repeat个）auc面积的平均值
                     mean_roc为 n_repeat次嵌套交叉验证中所有验证集的平均roc曲线
注意：
1）外层折叠数的设计要考虑到二分类结果中较少的一类的数目，外层折叠数不能超过较少一类的数目，否则无法保证所有验证集中
都含有两种类型的二分类结果
2）外层折叠数与内层折叠数的设计要考虑到数据本身的构成结果，由于内层使用分层交叉验证，所以要保证内外折叠数都可以让两
种二分类结果的样本得以平均分配。如15个样本，6个结果为1，9个结果为0，二分类结果中较少的一类为1，有6个，外层折叠数可
以设置为3个，内层折叠数可以设置为2个，此时外层每一折有2个1，3个0，内层折叠每一折有2个1和3个0，这样保证了分层交叉验
证的实现。

绘图模块
roc.plot(model_list,roc_color=NA,model_name=NA)：
作用：该函数可以对多个模型经过repeat_ncv()函数返回值制成的列表进行多个模型平均roc曲线的绘制
参数：model_list为多个模型经过repeat_ncv()函数返回值制成的列表
           roc_color为roc曲线人为指定的颜色向量，向量中每一个颜色依次对应model_list中每一个模型，roc_color若没有指定则使用
           默认配色方案
           model_name为模型的名称向量，向量中每一个名称依次对应model_list中每一个模型，model_name为必须指定的参数
返回值：roc曲线图

多重检验模块（Holm–Bonferroni校正的非配对双尾t检验）
ttest_HB(model_list,model_name)：
作用：该函数可以对model_list中包含的模型进行两两n_repeat个mean_AUROC的Holm–Bonferroni校正的非配对双尾t检验
参数：model_list为多个模型经过repeat_ncv()函数返回值制成的列表
           model_name为模型的名称向量，向量中每一个名称依次对应model_list中每一个模型，model_name为必须指定的参数
返回值：函数ttest_HB()的返回值为一个数据框，model列为进行t检验的两个模型的名称；Pvalue列为检验P值，带"*"代表显著；order列
为Pvalue的降序排列的序号，alpha列比较检验的Holm–Bonferroni显著性，对应行的Pvalue必须要小于这个值才算显著

L2逻辑回归建模模块：
L2log(x,y,nfolds)：
参数：x为数据框中自变量切片
           y为数据框中因变量的切片
           nfolds为交叉检验折叠数
返回值：lambda取lambda.min和lambda.1se两个值时的模型系数
              lambda取lambda.min模型准确率更高
              lambda取lambda.1se模型更精简，但牺牲准确率

具体执行过程见代码测试区域







