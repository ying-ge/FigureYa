需求：分析空间转录组测序文件中的微生物组信息
参考文献：Effect of the intratumoral microbiota on spatial and cellular heterogeneity in cancer
	DOI: https://www.nature.com/articles/s41586-022-05435-0
示例数据下载位置：
	https://github.com/FredHutch/Galeano-Nino-Bullman-Intratumoral-Microbiota_2022
输入：
	样本数据：空间转录组样本的数据信息（包括FASTQ，tiff）
	参考信息：spaceranger参考基因组，pathseq参考组
输出：	Seurat对象【Results/*/outs/ST.seu.rds】
	微生物组RNA总量空间分布图【Results/*/outs/*TotalMicrobe.pdf】

注意：该流程对内存和存储要求很高，内存应当在200G以上，可支配存储应在300G以上（150G参考基因组+每个样本30G~100G数据量）

使用方法：
①在InputData下按样本放置输入文件，以样本CRC_16为例，相关文件放在InputData/CRC_16路径下，即
	InputData/CRC_16/CRC_16_*.fastq.gz
	InputData/CRC_16/CRC_16.tiff
	InputData/CRC_16/V10S15-020.gpr
②在InputData/sample_sheet.csv输入样本信息（包括sample,slide,area）
	CRC_16,V10S15-020,D1
③在Codes/S2.pathseq.sh中修改gatk变量，以保证能够正常调用gatk软件
	gatk变量可定义为"gatk"或"java -Xmx100G -jar path/to/gatk.jar"
④运行STpathseq.sh脚本，【在主目录下运行STpathseq.sh】
	或根据STpathseq.sh的用法分步运行Codes/S1.spaceranger.sh、Codes/S2.pathseq.sh和Codes/S3.postprocess.sh

其他说明
1：关于GATK参考组和pathseq参考组的下载
	切换到Reference目录下，运行BuildReference.sh脚本，脚本会自动下载所需参考组，并解压到既定位置
2：关于Codes文件夹下S1~S3.sh脚本的功能和用法介绍 【在主目录下bash *.sh】
	S1.spaceranger.sh	参数sample（样本名称）,slide,area，slide和area需询问测序方
	S2.pathseq.sh		参数sample（样本名称）
	S3.postprocess.sh	参数sample（样本名称）
3：本流程仅在样本支持的情况下能够提取微生物组信息，绝大部分的公开数据集是不能满足该条件的。
	关于如何才能生成符合条件的样本数据，请参看原文献的protocol或咨询测序方
	