AUC <- function(formula,data,plot=TRUE)
{
    cox1 <- coxph(formula, data)
    y <- cox1[["y"]]
    p <- ncol(y)
    time <- y[,p-1]
    status <- y[,p]
    x <- cox1$linear.predictors
    tt <- sort(unique(time[status==1]))
    nt <- length(tt)
    AUCt <- rep(NA,nt)
    numsum <- denomsum <- 0
    for (i in 1:nt) {
        ti <- tt[i]
        # risk set
        Y <- sum(time>=ti)
        R <- which(time>ti) # !!! R(ti) is which(time>=ti), but for the "controls", j=i should be excluded, only valid without ties
        xi <- x[time==ti]
        num <- sum(x[R]<xi) + 0.5*sum(x[R]==xi)
        AUCt[i] <- num/(Y-1) # Also only valid in absence of ties
        numsum <- numsum + num
        denomsum <- denomsum + Y-1 # Also only valid in absence of ties
    }
    AUC <- numsum/denomsum
    if (plot) {
        plot(tt,AUCt,xlab="Time t",ylab="AUC(t)")
        lines(lowess(data.frame(tt,AUCt)))
        abline(h=0.5,lty=3)
    }
    return(list(AUCt=data.frame(time=tt,AUC=AUCt),AUC=AUC))
}
