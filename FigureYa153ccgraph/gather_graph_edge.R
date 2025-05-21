#' Gather graph node from data.frame
#'
#' @param df a data.frame
#' @param index  grouping
#' @param root root name
#'
#' @return a tibble
#' @export
#'
#' @examples
#'  ibrary(ccgraph)
#'  data(OTU)
#'  edges <- gather_graph_edge(OTU,index = c("p","c","o"))
#' @importFrom tidyr unite
#'
#' @name gather_graph
gather_graph_edge <- function(df,index=NULL,root=NULL){
  require(dplyr)
  if (length(index) < 2){
    stop("please specify at least two index column(s)")
  } else if (length(index)==2){
    data <- df %>% mutate(from=.data[[index[[1]]]]) %>%
      tidyr::unite(to,index,sep="/") %>%
      select(from,to) %>%
      mutate_at(c("from","to"),as.character)
  } else {
    list <- lapply(seq(2,length(index)), function(i){
      dots <- index[1:i]
      df %>% tidyr::unite(from,dots[-length(dots)],sep = "/",remove = F)  %>%
        tidyr::unite(to,dots,sep="/") %>%
        select(from,to) %>%
        mutate_at(c("from","to"),as.character)
    })
    data <- do.call("rbind",list)
  }
  data <- as_tibble(data)
  if (is.null(root)){
    return(data)
  } else {
    root_data <- df %>% group_by(.dots=index[[1]]) %>%
      summarise(count=n()) %>%
      mutate(from=root,to=as.character(.data[[index[[1]]]] )) %>%
      select(from,to)
    rbind(root_data,data)
  }

}
