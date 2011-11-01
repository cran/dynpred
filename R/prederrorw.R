pew <-
function (time, status, tsurv, survmat, tcens, censmat, width,
    FUN = c("KL", "Brier"), tout)
{
    ### Data needs to be ordered according to time (asc) and status (desc)
    ord <- order(time, -status)
    time <- time[ord]
    status <- status[ord]
    survmat <- survmat[,ord]
    censmat <- censmat[,ord]
    ### Both tsurv and tcens need to be sorted and need to start with 0
    if (any(!(tsurv==sort(tsurv)))) stop("argument 'tsurv' needs to be sorted")
    if (any(!(tcens==sort(tcens)))) stop("argument 'tcens' needs to be sorted")
    if (min(tsurv)<0) stop("no negative values allowed for 'tsurv'")
    if (min(tcens)<0) stop("no negative values allowed for 'tcens'")
    if (min(tsurv)>0) {
      tsurv <- c(0,tsurv)
      survmat <- rbind(rep(1,ncol(survmat)),survmat)
    }
    if (min(tcens)>0) {
      tcens <- c(0,tcens)
      censmat <- rbind(rep(1,ncol(censmat)),censmat)
    }
    ### Prepare for call to prederrw
    FUN <- match.arg(FUN)
    FUNn <- 2
    if (FUN == "Brier")
        FUNn <- 1
    n <- length(time)
    nsurv <- length(tsurv)
    ncens <- length(tcens)
    nout <- length(tout)
    res <- .C("prederrw", sn = as.integer(n), time = as.double(time),
        status = as.integer(status), snsurv = as.integer(nsurv),
        sncens = as.integer(ncens), snout = as.integer(nout),
        tsurv = as.double(tsurv), survmat = as.double(survmat),
        tcens = as.double(tcens), censmat = as.double(censmat),
        w = as.double(width), tout = as.double(tout), sFUN = as.integer(FUNn),
        err = double(nout), work = double(n))
    res <- data.frame(time = tout, Err = res$err)
    attr(res, "Score") <- FUN
    attr(res, "width") <- width
    return(res)
}

pewcox <- function(formula, censformula, width, data, censdata,
    FUN = c("KL", "Brier"), tout, CV = FALSE, progress = FALSE)
{
    cox1 <- coxph(formula, data)
    y <- cox1[["y"]]
    p <- ncol(y)
    time <- y[, p - 1]
    status <- y[, p]
    n <- length(time)
    if (nrow(data) != n)
        stop("missing data in time or status not allowed")
    ord <- order(time, -status)
    time <- time[ord]
    status <- status[ord]
    tt <- sort(unique(time[status == 1]))
    nt <- length(tt)
    if (!CV) {
        progress <- FALSE
        x <- cox1$linear.predictors[ord]
        cox1 <- coxph(Surv(time, status) ~ x)
        # if no covariates, then use Kalbfleisch-Prentice type in survfit
        if (sum(x^2)==0)
          sf <- survfit(cox1, newdata = data.frame(x = x), type="kalbfl")
        else
          sf <- survfit(cox1, newdata = data.frame(x = x))
        tt <- sf$time
        survmat <- sf$surv
    }
    if (CV) {
        x <- cox1$linear.predictors[ord]
        if (sum(x^2)==0) stop("No cross-validation for null model implemented")
        X <- model.matrix(formula, data = data)
        X <- X[, -1, drop = FALSE]
        X <- X[ord, , drop = FALSE]
        if (progress) {
            m <- floor(log10(n)) + 1
            pre <- rep("\b", 2 * m + 1)
            cat("Calculating cross-validated survival curves:\n")
        }
        survmat <- matrix(NA, nt, n)
        for (i in 1:n) {
            if (progress) {
                cat(pre, i, "/", n, sep = "")
                flush.console()
            }
            cmini <- coxph(Surv(time[-i], status[-i]) ~ X[-i,
                , drop = FALSE], method = "breslow")
            xi <- as.vector(X[-i, , drop = FALSE] %*% cmini$coef)
            ximini <- as.numeric(X[i, ] %*% cmini$coef)
            cmini <- coxph(Surv(time[-i], status[-i]) ~ xi, method = "breslow")
            survi <- survfit(cmini, newdata = data.frame(xi = ximini))
            survmat[, i] <- evalstep(survi$time, survi$surv,
                tt, subst = 1)
        }
        if (progress) {
            cat("\nCalculating prediction error ...")
            flush.console()
        }
    }
    if (tt[1] > 0) {
        tsurv <- c(0, tt)
        survmat <- rbind(rep(1, n), survmat)
    }
    else tsurv <- tt
    nsurv <- length(tsurv)
    ## censoring
    if (missing(censdata)) censdata <- data
    coxcens <- coxph(censformula, censdata)
    ycens <- coxcens[["y"]]
    p <- ncol(ycens)
    tcens <- ycens[, p - 1]
    dcens <- ycens[, p]
    xcens <- coxcens$linear.predictors[ord]
    coxcens <- coxph(Surv(tcens, dcens) ~ xcens)
    # if no covariates, then use Kalbfleisch-Prentice type in survfit
    if (sum(xcens^2)==0)
      sfcens <- survfit(coxcens, newdata = data.frame(xcens = xcens),
        type="kalbfl")
    else
      sfcens <- survfit(coxcens, newdata = data.frame(xcens = xcens))
    tcens <- sfcens$time
    censmat <- sfcens$surv
    if (tcens[1] > 0) {
        tcens <- c(0, tcens)
        censmat <- rbind(rep(1, n), censmat)
    }
    ncens <- length(tcens)
    ## vector at which prediction error is to be calculated
    if (missing(tout)) {
        tout <- unique(c(tsurv, tcens))
        tout <- c(tout, tout - width)
        tout <- tout[tout >= 0]
        if (min(tout) > 0)
            tout <- c(0, tout)
        tout <- sort(unique(tout))
        nout <- length(tout)
        if (nout > 0) {
            tout <- tout[-nout]
            nout <- nout - 1
        }
    }
    else {
        tout <- sort(tout)
        nout <- length(tout)
    }
    res <- pew(time, status, tsurv, survmat, tcens, censmat,
        width, FUN, tout)
    if (progress)
        cat("\n")
    return(res)
}
