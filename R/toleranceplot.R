toleranceplot <- function(formula,data,coverage=0.8,horizon,plot=TRUE,xlab)
{
    cox1 <- coxph(formula, data)
    y <- cox1[["y"]]
    p <- ncol(y)
    time <- y[,p-1]
    status <- y[,p]
    x <- cox1$linear.predictors
    xs <- sort(unique(x))
    nx <- length(xs)
    if (missing(horizon)) horizon <- max(time[status==1])*1.05
    ylim <- c(0,horizon)
    if (missing(xlab)) xlab <- "x"
    if (plot) plot(range(xs),c(0,horizon),type="n",xlab=xlab,ylab="Tolerance interval")
    cx <- coxph(Surv(time,status) ~ x, method="breslow")
    res <- matrix(NA,nx,3)
    for (i in 1:nx) {
        xi <- xs[i]
        nd <- data.frame(x=xi)
        sf <- survfit(cx,newdata=nd)
        sf <- data.frame(time=sf$time,surv=sf$surv)
        low <- max(sf$time[sf$surv>1-(1-coverage)/2])
        up <- min(sf$time[sf$surv<(1-coverage)/2])
        if (is.infinite(up)) up <- horizon
        lines(rep(xi,2),c(low,up),type="l")
        res[i,] <- c(xi,low,up)
    }
    res <- as.data.frame(res)
    names(res) <- c("x","lower","upper")
    attr(res, "coverage") <- coverage
    attr(res, "horizon") <- horizon
    return(res)
}
