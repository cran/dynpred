cindex <- function(formula, data)
{
    cox1 <- coxph(formula, data)
    y <- cox1[["y"]]
    p <- ncol(y)
    time <- y[,p-1]
    status <- y[,p]
    x <- cox1$linear.predictors
    n <- length(time)
    ord <- order(time,-status)
    time <- time[ord]
    status <- status[ord]
    x <- x[ord]
    # pairs (i,j) for which the smallest observed time is an event time
    wh <- which(status==1)
    total <- concordant <- 0
    for (i in wh) {
        for (j in ((i+1):n)) {
            if (time[j] > time[i]) {# ties not counted
                total <- total + 1
                if (x[j] < x[i]) concordant <- concordant + 1
                if (x[j] == x[i]) concordant <- concordant + 0.5
            }
        }
    }
    return(list(concordant=concordant,total=total,cindex=concordant/total))
}
