ExecuteCNMF <- function (datasets, clusterNum, nrun = 30) 
{
    if (is.list(datasets)) {
        temp = NULL
        for (i in 1:length(datasets)) {
            temp = rbind(temp, datasets[[i]])
        }
    }
    else temp = datasets
    data1 = rbind(pmax(temp, 0), -pmin(temp, 0)) #Kim & Tidor’s trick
    index = which(rowSums(data1) == 0)
    if(length(index)>0) {
      data1 = data1[-index, ] # 源代码的bug，若没有行和为0,则不run此行
    }
    res = nmf(data1, rank = clusterNum, nrun = nrun)
    distanceMatrix = slot(res, "consensus")
    attr(distanceMatrix, "class") = "Similarity"
    group = as.numeric(as.vector(predict(res)))
    result = list(group = group, distanceMatrix = distanceMatrix, 
        originalResult = res)
}