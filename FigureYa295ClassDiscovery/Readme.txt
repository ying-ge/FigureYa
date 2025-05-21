算法流程：
1：筛选探针，结果储存在probe.stat对象中
	计算各探针的方差(var)、卡方检验单侧p值(pvalue)和稳健变异系数(rCV，去除最大值和最小值后的变异系数)
	根据pvalue<0.01和rCV<10过滤低质量探针
	根据rCV使用不同分位数阈值(0.6, 0.7, 0.8, 0.9, 0.95, 0.975, 0.99, 0.995)，对探针进行筛选，特征集储存在probe.stat$sets中
2：生成聚类树
	根据8个阈值和3种聚类方法(average, complete, ward)，生成24种聚类树，结果储存在Dends对象中，聚类树图储存在Dendrograms.pdf中
3：树稳定性评估
	扰动：生成均值为0，方差为1.5倍探针中位方差的噪声矩阵noise
	重抽样：将添加噪声矩阵的样本加入到原样本中（进行重新聚类）
	使用RF.dist函数评估聚类树间的相似性
	因为本步骤并未对下游聚类结果产生影响，文中也没有给出树稳定性评估的结果，故此处仅列出扰动、重抽样，相似性得分的计算方法
4：识别robust cluster
	IsCoClass：在指定聚类树(dend)下，一对样本(sampleA, sampleB)在聚为k类时，出现在同一组中，且该组样本量不小于4，则认为本聚类树将该对样本识别为一类
	对不同的k，分别计算聚类共识矩阵（在指定的k下，将一对样本识别为一类的聚类树数目）
	聚类共识矩阵储存在mat.list对象中，共识热图储存在ConsensusMatrix.pdf中
	
	选取合适的k（此处选取k为3），将至少22棵聚类树识别为同类的样本视为一类(cutree, h=24-22)，得到初步的聚类结果，储存在sampleInfo对象的PrimaryCluster列中
	对PrimaryCluster中的每一类，将平均至少有18棵聚类树识别为同类的类别视为一类(cutree, h=24-18)，得到最终聚类结果，储存在sampleInfo对象的FinalCluster列中
5：选取差异基因构建分类器
	由于第4步得到的聚类结果与原文略有出入，此处使用的是原文中rC1和rC2的分组结果
	使用limma对rC1-rC2做差异分析，选取padj<0.001的探针构建分类器
6：使用classpredict包进行分类
	classpredict包是BRB ArrayTools的R版本，由于原文使用的Compound Covariate Predictor(CCP)没有对应的包，所以此处使用classpredict包进行类别预测
	预测结果在resList$predNewSamples中，保存在ClassPrediction.txt中