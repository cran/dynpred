scatterplot <- function(formula,data,horizon,plot=TRUE,xlab)
{
    cox1 <- coxph(formula, data)
    y <- cox1[["y"]]
    p <- ncol(y)
    time <- y[,p-1]
    status <- y[,p]
    x <- cox1$linear.predictors
    if (missing(horizon)) horizon <- max(time[status==1])*1.05
    ylim <- c(0,1.05*horizon) # extra space to print R-square
    imputed <- time
    cens <- which(status==0)
    events <- which(status==1)
    ncens <- length(cens)
    if (missing(xlab)) xlab <- "x"
    cx <- coxph(Surv(time,status) ~ x, method="breslow")
    for (i in 1:ncens) {
        xi <- x[cens[i]]
        nd <- data.frame(x=xi)
        sf <- survfit(cx,newdata=nd)
        sf <- data.frame(time=sf$time,surv=sf$surv)
        sf <- sf[!duplicated(sf$surv),]
        if (imputed[cens[i]]>max(sf$time))
            imputed[cens[i]] <- horizon
        if (imputed[cens[i]]<max(sf$time)) {
            sf1 <- sf[sf$time <= imputed[cens[i]],]
            sf2 <- sf[sf$time >= imputed[cens[i]],]
            sf2$surv <- 1-sf2$surv/min(sf1$surv)
            if (sf2$surv[1]==0) sf2 <- sf2[-1,]
            rand <- runif(1)
            survtimes <- c(0,sf2$surv,1)
            idx <- as.numeric(cut(rand,survtimes))
            if (idx>nrow(sf2)) imputed[cens[i]] <- horizon
            if (idx<=nrow(sf2)) imputed[cens[i]] <- sf2$time[idx]
        }
    }
    res <- data.frame(x=x, imputed=imputed)
    attr(res, "horizon") <- horizon
    if (plot) {
        plot(x, imputed, type="n", ylim=ylim, xlab=xlab, ylab="(Imputed) survival time")
        points(x[cens], imputed[cens])
        points(x[events], imputed[events], pch=16)
        abline(lsfit(x,imputed),lty=2)
        lmova <- lm(imputed ~ x)
        r2 <- summary(lmova)$r.squared
        text(max(x),1.05*horizon,paste("R-squared =",round(r2,3)),adj=1)
    }
    return(res)
}
