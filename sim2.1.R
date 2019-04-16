lseq <- function(from=1, to=100000, length.out=6) {
  # logarithmic spaced sequence
  # blatantly stolen from library("emdbook"), because need only this
  exp(seq(log(from), log(to), length.out = length.out))
}

CreateSqCovMat <- function(p, rho = 0,
                           type = c('cross-section', 'time-series', 'plus-minus-zero', 'symmetric-range', 'theory-corr-sparse'), alpha = 0.4) {
  require(expm)
  
  # Outputs the square-root of Covariance matrix
  ### rho: correlation between covariates
  ### type: pattern of correlation
  ###
  
  type <- type[1]
  if (type == 'cross-section') {
    sigMat <- matrix(rho, p, p)
    sigMat[cbind(1:p, 1:p)] <- 1
    M <- sqrtm(sigMat)
  } else if (type == 'time-series') {
    rho_pow <- rho^(seq_len(p)-1)
    rho_pow_rev <- rho_pow[p:1]
    sigMat <- matrix(0, p, p)
    for (i in 1:p) {
      sigMat[i, ] <- c(rho_pow_rev[p-i+seq_len(i-1)], rho_pow[seq_len(p-i+1)])
    }
    M <- sqrtm(sigMat)
  } else if (type == 'plus-minus-zero') {
    sigMat <- diag(rep(1.0, p))
    sigMat[lower.tri(sigMat)] <- sample(c(0, -rho, rho), p*(p-1)/2, replace = TRUE)
    sigMat <- crossprod(sigMat)
    M <- sqrtm(sigMat)
  } else if (type == 'symmetric-range') {
    sigMat <- diag(rep(1.0, p))
    sigMat[lower.tri(sigMat)] <- runif(p*(p-1)/2, -rho, rho)
    sigMat <- crossprod(sigMat)
    M <- sqrtm(sigMat)
  } else if (type == 'theory-corr-sparse') {
    M <- diag(rep(1.0, p))
    if (rho > 0 && alpha > 0) {
      multi_bernoulli <- as.integer(cut(runif(p*(p-1)/2), breaks = c(0, alpha/2, alpha, 2), include.lowest = TRUE))
      M[upper.tri(M)] <- ifelse(multi_bernoulli == 1, -rho, ifelse(multi_bernoulli == 2, rho, 0))
    }
  } else {
    stop('type can either be a cross-section or time-series')
  }
  
  # Output sq. root of cov. mat and other info.
  list(mat = M, rho=rho, type=type, alpha=alpha)
}

CreateCovariates <- function(n, sqCovMat){
  # covariates
  p <- dim(sqCovMat)[1]
  unif <- sqrt(3)*matrix(runif(n*p, min=-1, max=1), n, p) # Uniform sample, with unit s.d.
  unif %*% (sqCovMat) #correlated uniform sample
}


Population3 <- function(sqCovMat, s, interact_frac, n, sig, additive = FALSE, vary_sign = FALSE) {
  ### s: size of active covariates (support)
  ### interact_frac: fraction of active interaction terms
  ### n: sample size in each trial
  ### sig: std of noise term
  ### sqCovMat : square root of covariance matrix
  ### Output: x, models
  ###         models[[i]]: s, supp, frac, y
  
  p <- dim(sqCovMat$mat)[1]
  x <- CreateCovariates(n, sqCovMat$mat)
  
  shuff <- sample.int(p)
  
  # noise term
  omega <- rnorm(n, sd=sig)
  
  #popname <- sprintf('p%d_s%d_frac%2.0f%%', p, s, 100*interact_frac)
  #cat(sprintf('Preparing %s\n', popname))
  
  #y <- numeric(n_max)
  active_set <- shuff[1:s]
  # interaction selection
  interact_supp <- matrix(runif(s*s) < interact_frac, s, s)
  interact_supp[upper.tri(interact_supp)] <- FALSE
  single_supp <- active_set[rowSums(interact_supp) + colSums(interact_supp) == 0]
  supp_size <- sum(interact_supp) + length(single_supp)
  theta <- 1/sqrt(supp_size)
  if (vary_sign) {
    sign_pool <- c(-1, 1)  
  } else {
    sign_pool <- c(1)
  }
  
  theta_rnd_sgn <- matrix(sample(sign_pool, length(single_supp), replace = TRUE), ncol = 1)
  y <- theta * x[, single_supp] %*% theta_rnd_sgn
  if (additive && length(single_supp) > 0) {
    Additive <- function(x) exp(x)
    # Additive <- function(x, a, b) ifelse(x < a, 0, ifelse(x < b, 1, 0))
    # Additive <- function(x) exp(3*x) + 2*exp(-3*x)
    y <- y + 2 * matrix(Additive(x[, single_supp]), ncol = length(single_supp)) %*% abs(theta_rnd_sgn)
  }
  # It will speed things up if you vectorize the following two for loops
  for (i in 1:s) {
    for (j in 1:i) {
      if (interact_supp[i, j]) 
        y <- y + sample(sign_pool, 1) * theta * x[, active_set[i]] * x[, active_set[j]]
    }
  }
  y <- y + omega
  
  return(list(x = x, s = s, frac = interact_frac, sqCovMat = sqCovMat, additive = additive, sign = vary_sign, y = y, supp = shuff[1:s]))
}

########################
Simulate3par <- function(gen_model, trial_count, 
                         methods = c("Random", "Sis", "Lasso", "Spam", "Dstump", "TreeWeight", "TreeWeightPermutation"),
                         res_dir = "results") {
  # This first run can be avoided if we don't print
  model <- gen_model()
  n <- nrow(model$x)
  p <- ncol(model$x)
  s <- model$s
  ifrac <- model$frac
  rho <- model$sqCovMat$rho
  sign <- model$sign
  type  <- model$sqCovMat$type
  alpha <- model$sqCovMat$alpha
  additive_str <- substring(as.character(model$additive), 1, 1)
  popname <- sprintf('p%d_s%d_frac%2.0f%%_n%d_rho%2.0f_additive%s_sign%s_type%s_alpha%.2f',
                     p, s, 100*ifrac, n, 100*rho, additive_str, sign, type, 100*alpha)
  cat(sprintf('Using %s\n', popname))
  
  # Actual run
  acc <- matrix(0, nrow=trial_count, ncol=length(methods))
  #colnames(acc) <- methods
  
  M <- length(methods)
  acc_full_rtime <- foreach(tr=1:trial_count, 
                       .combine=rbind, 
                       .export=c("Population3", "Accuracy", "FullRecovery", "CreateSqCovMat", "CreateCovariates", methods)) %dopar% {
                         
                         # cat(sprintf('%3d ', tr))
                         res_fname <- file.path(res_dir, sprintf('%d_%s',tr,popname))
                         if (file.exists(res_fname)){
                           temp <- read.csv(file= res_fname, sep = ",")
                           # temp_acc <- temp[1,1:M]
                           # temp_rtime <- temp[1,M+(1:M)]
                           
                         } else {
                           temp <- list()        
                           for (md in methods) {
                             #cat(sprintf('# %10s ', md))
                             model <- gen_model()
                             stime <- Sys.time()
                             supp_hat <- do.call(md, list(x = model$x, y=model$y, s = model$s))
                             etime <- Sys.time()
                             temp[sprintf('%s_rtime',md)] <- etime - stime # run time
                             temp[md] <- Accuracy(supp_hat, model$supp)  # Acc
                             temp[sprintf('%s_rec',md)] <- FullRecovery(supp_hat, model$supp)  # full recovery
                           }
                           # shenanigan needed to write a vector as row to the file
                           #write.csv(as.matrix(t(temp_acc)), file= res_fname, sep = ",", row.names=FALSE)
                           temp <- data.frame(temp)
                           write.csv(temp, file= res_fname, sep = ",", row.names=FALSE)
                           #write.csv(as.matrix(t(c(temp_acc,temp_rtime))), file= res_fname, sep = ",", row.names=FALSE)
                           
                         }
                         #temp_acc #Equivalent to acc = rbind(acc, temp_acc)
                         temp
                       }
  
  #return(list(p = p, s = s, n = n, frac = ifrac, accuracy = acc))
  #return(list(accuracy = acc, runtime=runtime))
  
  full_rec <- acc_full_rtime[,paste(methods,'_rec',sep='')]
  names(full_rec) <- methods
  list(accuracy = acc_full_rtime[,methods], full_rec = full_rec, rtime=acc_full_rtime[,paste(methods,'_rtime',sep='')])
}

Accuracy <- function(supp_hat, supp) {
  s <- length(supp)
  if (length(supp_hat) > s) {
    supp_hat <- supp_hat[1:s]
  }
  mean(supp %in% supp_hat)
}

FullRecovery <- function(supp_hat, supp) {
  s <- length(supp)
  if (length(supp_hat) > s) {
    supp_hat <- supp_hat[1:s]
  }
  as.numeric(all(supp %in% supp_hat))
}

