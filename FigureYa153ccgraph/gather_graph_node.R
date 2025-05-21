#' Gather graph node from data.frame
#'
#' @inheritParams gather_graph
#' @param value
#'
#' @return a tibble
#' @export
#'
#' @examples
#'  library(ccgraph)
#'  data(OTU)
#'  nodes <- gather_graph_node(OTU,index = c("p","c","o"))
#' @importFrom tidyr unite
gather_graph_node <- function(df,index=NULL,value=tail(colnames(df),1),root=NULL){
  require(dplyr)
  if (length(index) < 2){
    stop("please specify at least two index column(s)")
  } else {
    list <- lapply(seq_along(index), function(i){
      dots <- index[1:i]
      df %>%
        group_by(.dots=dots) %>%
        summarise(node.size=sum(.data[[value]]),
                  node.level=index[[i]],
                  node.count=n()) %>%
        mutate(node.short_name=as.character(.data[[ dots[[length(dots)]] ]]),
               node.branch = as.character(.data[[ dots[[1]]]])) %>%
        tidyr::unite(node.name,dots,sep = "/")
    })
    data <- do.call("rbind",list) %>% as_tibble()
    data$node.level <- factor(data$node.level,levels = index)

    if (is.null(root)){
      return(data)
    } else {
      root_data <- data.frame(node.name=root,
                              node.size=sum(df[[value]]),
                              node.level=root,
                              node.count=1,
                              node.short_name=root,
                              node.branch=root,
                              stringsAsFactors = F)
      data <- rbind(root_data,data)
      data$node.level <- factor(data$node.level, levels = c(root,index))
      return(data)
    }
  }
}
