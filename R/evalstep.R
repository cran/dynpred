evalstep <- function(time, stepf, newtime, subst=-Inf, to.data.frame=FALSE)
{
    n <- length(time)
    if (is.vector(stepf))
        if (length(stepf) != n)
            stop("arguments 'time' and 'stepf' should have the same length")
    if (is.matrix(stepf) | is.data.frame(stepf))
        if (nrow(stepf) != n)
            stop("argument 'stepf' should have the same number of rows as length of argument 'time'")
    # time should be ordered, not contain duplicates, and not contain +/- Inf
    if (any(!(order(time) == 1:n))) stop("argument 'time' should be ordered")
    if (any(duplicated(time))) stop("argument 'time' should not contain duplicates")
    if (any(is.infinite(time))) stop("(-) infinity not allowed in 'time'")
    idx <- cut(newtime,c(-Inf,time,Inf),right=FALSE)
    idx <- as.numeric(idx)
    if (is.vector(stepf)) res <- c(subst,stepf)[idx]
    if (is.matrix(stepf) | is.data.frame(stepf)) {
        stepf <- rbind(rep(subst,ncol(stepf)),stepf)
        res <- stepf[idx,]
    }
    if (to.data.frame) return(data.frame(newtime=newtime,res=res))
    else return(res)
}
