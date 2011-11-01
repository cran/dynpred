CVcindex <- function(formula, data, type="single", matrix=FALSE)
{
    cox1 <- coxph(formula, data)
    y <- cox1[["y"]]
    p <- ncol(y)
    time <- y[,p-1]
    status <- y[,p]
    # Center covariates
    # First get design matrix
    X <- model.matrix(formula, data = data)
    X <- X[,-1] # remove intercept
    X <- t(t(X) - apply(X,2,mean))
    cfull <- coxph(Surv(time,status) ~ X, data = data, method="breslow")
    n <- nrow(data)
    if (type=="single" | type=="pair") {
        m <- floor(log10(n))+1 # for formatting progress count
        pre <- rep("\b",2*m+1)
        xmat <- matrix(NA,n,n) # in columns (-i)
        # xmat[j,i] will contain PI_{j,(-i)}
        for (i in 1:n) {
            cat(pre,i,"/",n,sep=""); flush.console()
            # leave out i
            cmin1 <- coxph(Surv(time[-i], status[-i]) ~ X[-i,], method="breslow")
            # evaluate at all j except i
            xmat[-i,i] <- cmin1$linear.predictors
            # evaluate at i
            xmat[i,i] <- sum(X[i,] * cmin1$coef)
        }
        cat("\n")
    }
    if (type=="single") {
        formula <- as.formula("Surv(time,status) ~ x")
        ndata <- data.frame(time=time,status=status,x=diag(xmat))
        res <- cindex(formula=formula, data=ndata)
        if (matrix) res <- list(concordant=res$concordant,total=res$total,cindex=res$cindex,matrix=xmat)
    }
    if (type=="pair") {
        n <- length(time) # check if = length(status) and length(x)
        ord <- order(time,-status)
        time <- time[ord]
        status <- status[ord]
        xmat <- xmat[ord,ord]
        # pairs (i,j) for which the smallest observed time is an event time
        wh <- which(status==1)
        total <- concordant <- 0
        for (i in wh) {
            if (i < n) {
                for (j in ((i+1):n)) {
                    if (time[j] > time[i]) { # ties not counted
                        total <- total + 2
                        if (xmat[j,i] < xmat[i,i]) concordant <- concordant + 1
                        if (xmat[j,j] < xmat[i,j]) concordant <- concordant + 1
                    }
                }
            }
        }
        if (matrix) res <- list(concordant=concordant,total=total,cindex=concordant/total,matrix=xmat) else res <- list(concordant=concordant,total=total,cindex=concordant/total)
    }
    if (type=="fullpairs") {
        m <- floor(log10(n*(n-1)/2))+1 # for formatting progress count
        pre <- rep("\b",2*m+1)
        # xmat[i,j] will contain PI_{i,(-i,-j)}; xmat[j,i] will contain PI_{j,(-i,-j)}
        xmat <- matrix(NA,n,n)
        cnt <- 0
        for (i in 1:n) {
            if (i < n) {
                for (j in ((i+1):n)) {
                    cnt <- cnt+1
                    cat(pre,cnt,"/",n*(n-1)/2,sep=""); flush.console()
                    # leave out i and j
                    cmin2 <- coxph(Surv(time[-c(i,j)], status[-c(i,j)]) ~ X[-c(i,j),], method="breslow")
                    # evaluate at i
                    xmat[i,j] <- sum(X[i,] * cmin2$coef)
                    # evaluate at j
                    xmat[j,i] <- sum(X[j,] * cmin2$coef)
                }
            }
        }
        ord <- order(time,-status)
        time <- time[ord]
        status <- status[ord]
        xmat <- xmat[ord,ord]
        # pairs (i,j) for which the smallest observed time is an event time
        wh <- which(status==1)
        total <- concordant <- 0
        for (i in wh) {
            if (i < n) {
                for (j in ((i+1):n)) {
                    if (time[j] > time[i]) {# ties not counted
                        total <- total + 1
                        if (xmat[j,i] < xmat[i,j]) concordant <- concordant + 1
                    }
                }
            }
        }
        if (matrix) res <- list(concordant=concordant,total=total,cindex=concordant/total,matrix=xmat) else res <- list(concordant=concordant,total=total,cindex=concordant/total)
    }
    cat("\n")
    attr(res, "type") <- type
    return(res)
}
