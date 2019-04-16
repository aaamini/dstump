Dstump <- function(x, y, s) {
  p <- ncol(x)
  n <- length(y)
  r <- floor(n/2)
  mu <- mean(y)
  impurity <- numeric(p)

  # for (j in 1:p) {
  #   m <- median(x[, j])
  #   impurity[j] <- var(y[x[, j] <= m]) + var(y[x[, j] > m])
  # }

  for (j in 1:p) {
    y_st <- y[order(x[, j])]
    impurity[j] <- -abs(mu - mean(y_st[1:r]))
  }

  order(impurity)[1:s]
}

InfeasibleDstump <- function(x, y, s) {
  p <- ncol(x)
  n <- length(y)
  r <- floor(n/2)
  mu <- mean(y)
  impurity <- numeric(p)
  
  for (j in 1:p) {
    impurity[j] <- -abs(sum(y[x[, j] < 0] - mu) / r)
  }
  
  order(impurity)[1:s]
}

SlidingDstump <- function(x, y, s) {
  p <- ncol(x)
  n <- length(y)
  r <- floor(n/2)
  mu <- mean(y)
  impurity <- numeric(p)
  for (j in 1:p) {
    y_st <- y[order(x[, j])]
    y_st <- c(y_st, y_st[1:r])
    impurity[j] <- -max(abs(mu - diff(cumsum(y_st), r) / r))
    # index1 <- 1:r
    # index2 <- (r+1):n
    # impurity[j] <- min(sapply(1:r, function(i) var(y_st[index1+i-1]) + var(y_st[index2+i-1])))
  }
  order(impurity)[1:s]
}

DstumpLasso <- function(x, y, s) {
  require(lars)
  p <- ncol(x)
  n <- length(y)
  x_dstump <- matrix(0, nrow = n, ncol = p)
  for (j in 1:p) {
    x_dstump[, j] <- as.numeric(x[, j] <= median(x[, j]))
  }
  fit <- lars(x_dstump, y, intercept=TRUE, normalize=FALSE) #TODO try normalize=TRUE
  path <- coef(fit)
  selectionCount <- rowSums(path != 0)
  if (any(selectionCount >= s)) {
    referenceSelection <- which(selectionCount >= s)[1]
  } else {
    referenceSelection <- nrow(path)
  }
  which(path[referenceSelection, ] != 0)
}

Sis <- function(x, y, s) {
  correlation <- abs(crossprod(x, y))
  order(correlation, decreasing = TRUE)[1:s]
}

Lasso <- function(x, y, s) {
  require(lars)
  fit <- lars(x, y, intercept=TRUE, normalize=FALSE) #TODO try normalize=TRUE
  path <- coef(fit)
  selectionCount <- rowSums(path != 0)
  if (any(selectionCount >= s)) {
    referenceSelection <- which(selectionCount >= s)[1]
  } else {
    referenceSelection <- nrow(path)
  }
  which(path[referenceSelection, ] != 0)
}

Spam <- function(x, y, s) {
  require(SAM)
  path <- samQL(x, y)
  if (any(path$df >= s)) {
    referenceSelection <- which(path$df >= s)[1]
  } else {
    referenceSelection <- ncol(path$w)
  }
  nonzero_splines <- which(path$w[, referenceSelection] != 0)
  sort(unique((nonzero_splines - 1) %/% path$p) + 1)
}

TreeWeight <- function(x, y, s) {
  require(randomForest)
  fit <- randomForest(x=x, y=y, maxnodes=4)
  order(importance(fit, type = 2), decreasing=TRUE)[1:s]
}

TreeWeightPermutation <- function(x, y, s) {
  require(randomForest)
  fit <- randomForest(x=x, y=y, maxnodes=4, importance = TRUE)
  order(importance(fit, type = 1), decreasing=TRUE)[1:s]
}

Random <- function(x, y, s) {
  p <- ncol(x)
  sample.int(p, s)
}
