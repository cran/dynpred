deb <- function(x, method=c("print","cat")) {
  call <- match.call()
  method <- match.arg(method)
  if (method=="print") {
    cat(deparse(call$x), ":\n")
    print(x)
  }
  if (method=="cat") cat(deparse(call$x), ":", x, "\n")
  flush.console()
}
