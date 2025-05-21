# 函数目标:
# 读取文件中的指定行和指定列
# 不包括注释行
read_part <- function(file, rows = 1, columns = -1, sep = "\t",
                      stringsAsFactors = FALSE,
                      header = FALSE,
                      check.names = FALSE, 
                      comment.char = "#", ...){
  dfl <- list()
  if (grepl("gz$", file)){
    con <- gzfile(file, open = "rb")
  } else{
    con <- file(file, open = "r")
  }
  
  i <- 0
  j <- 1
  repeat{
    
    rec <- readLines(con, 1)
    if (length(rec) == 0) break
    i <- i + 1
 
    # 当rows = -1时, 会读取所有行 
    # 超过目标行时停止读取
    if (i > max(rows) & rows != -1) break  
    # 不考虑注释行
    if (grepl(comment.char, rec )) next
    if ( ! i %in% rows & rows != -1) next
    
    items <- strsplit(rec, split = sep, fixed = TRUE)[[1]]
    if ( columns == -1){
      select_cols <- items
    } else{
      select_cols <- items[columns]
    }
    #print(select_cols)
    dfl[[j]] <- select_cols
    j <- j + 1
    
    
  }
  close(con) 
  df <- do.call(rbind, dfl)
  return(df)
}
