Fwindow <- function(object, width, variance=TRUE, conf.level=0.95)
{
    if (variance)
        sf <- data.frame(time=object$time,surv=object$surv,varHaz=object$std.err^2)
    else sf <- data.frame(time=object$time,surv=object$surv)
    sf$Haz <- -log(sf$surv)
    tt <- c(0,sf$time) # assuming no event at t=0 or t<0
    ttt <- c(0,tt,tt-width)
    ttt <- ttt[ttt >= 0]
    ttt <- sort(unique(ttt))
    ttt <- unique(ttt)
    H <- outer(c(0,sf$Haz),c(0,sf$Haz),"-")
    dimnames(H) <- list(tt,tt)
    tt <- c(tt,Inf)
    idx1 <- as.numeric(cut(ttt,tt,right=FALSE))
    idx2 <- as.numeric(cut(ttt+width,tt,right=FALSE))
    Fw <- diag(H[idx2,idx1])
    nt <- length(Fw)
    Fw[nt] <- Fw[nt-1]
    if (variance) {
        varH <- outer(c(0,sf$varHaz),c(0,sf$varHaz),"-")
        varFw <- diag(varH[idx2,idx1])
        varFw[nt] <- varFw[nt-1]
        ciwidth <- qnorm(1-(1-conf.level)/2)*sqrt(varFw)
        low <- Fw - ciwidth
        up <- Fw + ciwidth
        low[low<0] <- 0 # negative lower values of cum hazard set to zero
    }
    # Return on probability scale
    Fw <- 1 - exp(-Fw)
    if (variance) {
        low <- 1 - exp(-low)
        up <- 1 - exp(-up)
        res <- data.frame(time=ttt,Fw=Fw,low=low,up=up)
    }
    else res <- data.frame(time=ttt,Fw=Fw)
    attr(res,"width") <- width
    return(res)
}
