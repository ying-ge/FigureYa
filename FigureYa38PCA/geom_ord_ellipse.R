##' add confidence ellipse to ordinary plot produced by ggord
##'
##' 
##' @title geom_ord_ellipse 
##' @param mapping aes mapping 
##' @param ellipse_pro confidence value for the ellipse
##' @param fill color to fill the ellipse, NA by default
##' @param ... additional parameters
##' @return ggplot layer
##' @importFrom ggplot2 aes_
##' @importFrom ggplot2 layer
##' @importFrom utils modifyList
##' @export
##' @author Guangchuang Yu
##' @references \url{http://lchblogs.netlify.com/post/2017-12-22-r-addconfellipselda/}
geom_ord_ellipse <- function(mapping = NULL, ellipse_pro = 0.97, fill = NA, ...) {
    default_aes <- aes_(color = ~Groups, group = ~Groups)
    if (is.null(mapping)) {
        mapping <- default_aes
    } else {
        mapping <- modifyList(default_aes, mapping)
    }
    
    layer(
        geom = "polygon",
        stat = StatOrdEllipse,
        mapping = mapping,
        position = 'identity',
        data = NULL,
        params = list(
            ellipse_pro = ellipse_pro,
            fill = fill,
            ...
        )
    )
}

##' @importFrom ggplot2 ggproto
##' @importFrom ggplot2 Stat
##' @importFrom plyr ddply
##' @importFrom grDevices chull
StatOrdEllipse <- ggproto("StatOrdEllipse", Stat,
                          compute_group = function(self, data, scales, params, ellipse_pro) {
                              names(data)[1:2] <- c('one', 'two')
                              theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
                              circle <- cbind(cos(theta), sin(theta))
                              ell <- ddply(data, .(group), function(x) {
                                  if(nrow(x) <= 2) {
                                      return(NULL)
                                  }
                                  sigma <- var(cbind(x$one, x$two))
                                  mu <- c(mean(x$one), mean(x$two))
                                  ed <- sqrt(qchisq(ellipse_pro, df = 2))
                                  data.frame(sweep(circle %*% chol(sigma) * ed, 2, mu, FUN = '+'))
                              })
                              names(ell)[2:3] <- c('one', 'two')
                              ell <- ddply(ell, .(group), function(x) x[chull(x$one, x$two), ])
                              names(ell) <- c('Groups', 'x', 'y')
                              return(ell)
                          },
                          required_aes = c("x", "y", "group")
                          )


## . function was from plyr package
. <- function (..., .env = parent.frame()) {
    structure(as.list(match.call()[-1]), env = .env, class = "quoted")
}