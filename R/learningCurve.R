
#' Fit learning curve for accuracy matrix
#' @param accMat Matrix of accuracy rate where column indicate different sample size
#' @param n Vector indicates the sample size
#' @param auto_initial whether automatical intialise
#' @param a input the parameter a starting point
#' @param b input the parameter a starting point
#' @param c input the parameter a starting point
#' @param plot indicates whether plot or not
#' @param verbose indicates whether verbose or not
#' @return list of results
#' 
#' @examples
#' set.seed(2019)
#' n <- seq(20, 10000, 100)
#' accMat <- do.call(cbind, lapply(1:length(n), function(i){
#' tmp_n <- rep(n[i], 50)
#' y <- -2/(tmp_n^0.8) + 0.95 + rnorm(length(tmp_n), 0, 0.02)
#' }))
#' res <- learningCurve(accMat = accMat, n)
#' N <- getN(res, acc = 0.9)
#' 
#' @author Yingxin Lin
#' @import minpack.lm ggplot2
#' @importFrom stats quantile coef predict
#' @export

learningCurve <- function(accMat, n, auto_initial = T,
                          a = NULL, b = NULL, c = NULL,
                          plot = T, verbose = T){

  if (ncol(accMat) != length(n)) {
    stop("Number of column doesn't match with the length of n")
  }

  dat <- data.frame(y = colMeans(accMat),
                    y_25 = apply(accMat, 2, function(x) quantile(x, 0.25)),
                    y_75 = apply(accMat, 2, function(x) quantile(x, 0.75)),
                    y_05 = apply(accMat, 2, function(x) quantile(x, 0.05)),
                    y_95 = apply(accMat, 2, function(x) quantile(x, 0.95))
  )

  # Fitting model
  model <- list()
  for (i in seq(ncol(dat))) {
    model[[i]] <- fitLC(dat[,i], n, auto_initial = auto_initial,
                        a = a, b = b, c = c,
                        verbose = verbose)
  }


  # Get fitted value
  dat$n <- n
  new_n <- data.frame(n = seq(min(n), max(n), by = 0.1))
  fit <- lapply(model, function(x) predict(x, newdata = new_n))

  names(model) <- names(fit) <- c("mean", "quantile_25", "quantile_75", "quantile_05", "quantile_75")
  fit[["n"]] <- new_n$n
  dat_fit <- data.frame(do.call(cbind, fit))

  if (plot) {

    cols <- c("Mean" = "#c8133b","Quantile25/75" = "#ea8783" ,"Quantile05/95" = "#feb5a2")
    g <-  ggplot2::ggplot(dat, ggplot2::aes(x = n, y = y))  +
      xlab("N") + ylab("Accuracy Rate") +
      ggplot2::geom_point(alpha = 0.7) +
      ggplot2::geom_line(data = dat_fit, aes(x = n, y = mean, color = "Mean"), linetype = "solid", size = 1) +
      ggplot2::geom_line(data = dat_fit, aes(x = n, y = quantile_25,  color = "Quantile25/75"), linetype = "dashed") +
      ggplot2::geom_line(data = dat_fit, aes(x = n, y = quantile_75,  color = "Quantile25/75"), linetype = "dashed") +
      ggplot2::geom_line(data = dat_fit, aes(x = n, y = quantile_05,  color = "Quantile05/95"), linetype = "dotted") +
      ggplot2::geom_line(data = dat_fit, aes(x = n, y = quantile_75,  color = "Quantile05/95"), linetype = "dotted") +
      ggplot2::scale_color_manual(values = cols) +
      # scale_x_continuous(trans = scales::log_trans()) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "bottom")

    return(list(model = model, fit = fit, plot = g))

  }else{
    return(list(model = model, fit = fit))
  }
}


init_fit <- function(acc, n){
  # y1 <- log(max(acc) + 0.001 - acc)
  # x1 <- log(n)
  fit <- lm(log(n) ~ log(max(acc) + 0.001 - acc))

  c <- -coef(fit)[2]
  a <- -exp(coef(fit)[1])
  return(list(a = a, c = c))
}

# Function to fit the learning curve by inverse power function

fitLC <- function(acc, n, auto_initial = T,
                  a, b, c, verbose = T){
  dat_train <- data.frame(n = n, acc = acc)

  if (auto_initial) {
    init <- init_fit(acc, n)
    learning_curve <-  try({minpack.lm::nlsLM(acc ~ I(1 / n^(c) * a) + b,
                                  data = dat_train,
                                  start = list(a = init$a, c = init$c, b = 1))}, silent = T)


  }else if (!is.null(a) & !is.null(b) & !is.null(c)) {
    # For the case that user supplies starting point
    learning_curve <- stats::nls(acc ~ I(1 / n^(c) * a) + b,
                          data = dat_train,
                          start = list(a = a, c = c, b = b))
  }else{
    # this starting point may not work all the time
    learning_curve <- stats::nls(acc ~ I(1 / n^(c) * a) + b,
                          data = dat_train,
                          start = list(a = -10, c = 1, b = 1))
  }

  if (verbose) {
    print(summary(learning_curve))
  }

  return(learning_curve)
}



#' Function to get the required N given by the accuracy and the learning curve model
#' @param res model results returned by \code{learning_curve} function
#' @param acc accuracy that are quired
#' @return sample size that are required
#'
#' @export

getN <- function(res, acc = 0.9) {
  para <- coef(res$model$mean)
  names(para) <- c("a", "c", "b")
  if (acc > para["b"]) {
    stop("Required accuracy is too high to achieve :(")
  }
  N <- round(exp(1/para["c"]*log(para["a"]/(acc - para["b"]))))
  names(N) <- NULL

  return(N)
}

