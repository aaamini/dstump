rm(list=ls())
REMOVE_RES = T
REMOVE_PLOT = T

library(foreach)
library(doParallel)
source('methods.R')
source('sim2.1.R')

start.time <- Sys.time()
#cat(sprintf("Started at %s\n", start.time)


p = 200
# s_arr <- floor(seq(p/20, p/2, length.out=10))
s_arr <- floor(lseq(5, p/2, length.out=10))
#s_arr <- floor(seq(p/20, p/2, length.out=4))

# # non parallel version
# out <- numeric()
# for (s in s_arr){
#   gen_model <- function() Population2(p=p, s=s, interact_frac = .1, n_max=120, sig=0.1)
#   out <- rbind(out, colMeans(Simulate2(gen_model, trial_count = 3)$accuracy))
# }

#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

TRAIL_COUNT <- 100
res_dir <- "results"
dir.create(file.path(res_dir), showWarnings = FALSE)
if (REMOVE_RES) file.remove(file.path(res_dir, list.files(res_dir)))
plot_dir <- "plots"
dir.create(file.path(plot_dir), showWarnings = FALSE)
if (REMOVE_PLOT) file.remove(file.path(plot_dir, list.files(plot_dir)))

args_df <- rbind(expand.grid(frac = c(0), n = 1024, additive = c(FALSE), sign = c(TRUE), rho = c(0.1), type = c('theory-corr-sparse'), alpha = c(0.08)),
                 expand.grid(frac = c(0, 0.1), n = 1024, additive = c(FALSE, TRUE), sign = c(TRUE), rho = c(0), type = c('theory-corr-sparse'), alpha = c(0)),
                 expand.grid(frac = c(0.1), n = 1024, additive = c(TRUE), sign = c(TRUE), rho = c(0.1), type = c('theory-corr-sparse'), alpha = c(0.08)))
for (a_idx in seq_len(nrow(args_df))) {
  mean_acc <- numeric()
  lower_acc <- numeric()
  upper_acc <- numeric()
  mean_rec <- numeric()
  lower_rec <- numeric()
  upper_rec <- numeric()
  mean_rtime <- numeric()
  for (s in s_arr){
    args <- list(s=s, frac=args_df$frac[a_idx], n=args_df$n[a_idx], additive=args_df$additive[a_idx], sign = args_df$sign[a_idx],
                 sqCovMat = CreateSqCovMat(p = 200, rho = args_df$rho[a_idx], type=args_df$type[a_idx], alpha=args_df$alpha[a_idx]))
    
    clusterExport(cl, "args") # exports the arguments of gen_model
    
    gen_model <- function() Population3(sqCovMat= args$sqCovMat,
                                        s = args$s, 
                                        interact_frac = args$frac, 
                                        n = args$n, sig=0.1, additive = args$additive,
                                        vary_sign = args$sign)
    #out <- Simulate2par( gen_model, trial_count = 3)
    result <- Simulate3par(gen_model, trial_count = TRAIL_COUNT, res_dir=res_dir)
    
    mean_acc <- rbind(mean_acc, colMeans(result$accuracy))
    lower_acc <- rbind(lower_acc, apply(result$accuracy, 2, function(x) quantile(x, probs = c(0.025))))
    upper_acc <- rbind(upper_acc, apply(result$accuracy, 2, function(x) quantile(x, probs = c(0.975))))
    
    mean_rec <- rbind(mean_rec, colMeans(result$full_rec))
    lower_rec <- rbind(lower_rec, apply(result$full_rec, 2, function(x) quantile(x, probs = c(0.025))))
    upper_rec <- rbind(upper_rec, apply(result$full_rec, 2, function(x) quantile(x, probs = c(0.975))))
    
    mean_rtime <- rbind(mean_rtime, colMeans(result$rtime))
  }
  
  # Plots 
  method_count <- ncol(mean_acc)
  
  # accuracy
  
  plotname <- sprintf('acc_frac_%d_n%d_rho%.0f_additive%s_sign%s_type%s_alpha%.2f.png',
                      round(100*args$frac), args$n, 100*args$sqCovMat$rho, args$additive, args$sign, args$sqCovMat$type, 100*args$sqCovMat$alpha)
  if ((args$frac == 0) && (args$additive == FALSE)) {
    model_subtitle <- ""  
  } else if ((args$frac > 0) && (args$additive == FALSE)) {
    model_subtitle <- " + interaction"
  } else if ((args$frac == 0) && (args$additive == TRUE)) {
    model_subtitle <- " + additive"
  } else {
    model_subtitle <- " + interaction + additive"
  }
  
  if (args$sqCovMat$rho == 0) {
    pop_subtitle <- "iid"
  } else {
    pop_subtitle <- "correlated"
  }
  title <- sprintf('n = %d; linear%s, %s', args$n, model_subtitle, pop_subtitle)
  png(file.path(plot_dir, plotname), width = 800, height = 600)
  matplot(s_arr, mean_acc, type = 'b', pch = 16, lwd = 4, ylim = c(0, 1), cex.axis = 2, col=seq_len(method_count),
          main=NA,
          xlab=NA,
          ylab=NA)
  mtext(side = 3, cex = 2.5, line = 2, text = title)
  mtext(side = 1, cex = 2.5, line = 2, text = 's')
  mtext(side = 2, cex = 2.5, line = 2.5, text = 'Fraction of support recovered')
  legend('bottomright', colnames(mean_acc), col=seq_len(method_count), lty = seq_len(method_count), fill=seq_len(method_count), cex = 1.5)
  dev.off()
  
  
  # accuracy with confidence interval
  
  plotname <- sprintf('acc_ci_frac_%d_n%d_rho%.0f_additive%s_sign%s_type%s_alpha%.2f.png',
                      round(100*args$frac), args$n, 100*args$sqCovMat$rho, args$additive, args$sign, args$sqCovMat$type, 100*args$sqCovMat$alpha)
  if ((args$frac == 0) && (args$additive == FALSE)) {
    model_subtitle <- ""  
  } else if ((args$frac > 0) && (args$additive == FALSE)) {
    model_subtitle <- " + interaction"
  } else if ((args$frac == 0) && (args$additive == TRUE)) {
    model_subtitle <- " + additive"
  } else {
    model_subtitle <- " + interaction + additive"
  }
  
  if (args$sqCovMat$rho == 0) {
    pop_subtitle <- "iid"
  } else {
    pop_subtitle <- "correlated"
  }
  title <- sprintf('n = %d; linear%s, %s', args$n, model_subtitle, pop_subtitle)
  png(file.path(plot_dir, plotname), width = 800, height = 600)
  matplot(s_arr, mean_acc, type = 'b', pch = 16, lwd = 4, ylim = c(0, 1), cex.axis = 2, lty = 1, col=seq_len(method_count),
          main=NA,
          xlab=NA,
          ylab=NA)
  matplot(s_arr, lower_acc, type = 'l', pch = 16, lwd = 2, ylim = c(0, 1), cex.axis = 2, lty = 2, col=seq_len(method_count), add = TRUE,
          main=NA,
          xlab=NA,
          ylab=NA)
  matplot(s_arr, upper_acc, type = 'l', pch = 16, lwd = 2, ylim = c(0, 1), cex.axis = 2, lty = 2, add = TRUE,
          main=NA,
          xlab=NA,
          ylab=NA)
  mtext(side = 3, cex = 2.5, line = 2, text = title)
  mtext(side = 1, cex = 2.5, line = 2, text = 's')
  mtext(side = 2, cex = 2.5, line = 2.5, text = 'Fraction of support recovered')
  legend('bottomright', colnames(mean_acc), col=seq_len(method_count), lty = 1, fill=seq_len(method_count), cex = 1.5)
  dev.off()
  
  
  # full recovery
  
  plotname <- sprintf('rec_frac_%d_n%d_rho%.0f_additive%s_sign%s_type%s_alpha%.2f.png',
                      round(100*args$frac), args$n, 100*args$sqCovMat$rho, args$additive, args$sign, args$sqCovMat$type, 100*args$sqCovMat$alpha)
  if ((args$frac == 0) && (args$additive == FALSE)) {
    model_subtitle <- ""  
  } else if ((args$frac > 0) && (args$additive == FALSE)) {
    model_subtitle <- " + interaction"
  } else if ((args$frac == 0) && (args$additive == TRUE)) {
    model_subtitle <- " + additive"
  } else {
    model_subtitle <- " + interaction + additive"
  }
  
  if (args$sqCovMat$rho == 0) {
    pop_subtitle <- "iid"
  } else {
    pop_subtitle <- "correlated"
  }
  title <- sprintf('n = %d; linear%s, %s', args$n, model_subtitle, pop_subtitle)
  png(file.path(plot_dir, plotname), width = 800, height = 600)
  matplot(s_arr, mean_rec, type = 'b', pch = 16, lwd = 4, ylim = c(0, 1), cex.axis = 2, col=seq_len(method_count),
          main=NA,
          xlab=NA,
          ylab=NA)
  mtext(side = 3, cex = 2.5, line = 2, text = title)
  mtext(side = 1, cex = 2.5, line = 2, text = 's')
  mtext(side = 2, cex = 2.5, line = 2.5, text = 'Propability of full support recovery')
  legend('bottomright', colnames(mean_rec), col=seq_len(method_count), lty = seq_len(method_count), fill=seq_len(method_count), cex = 1.5)
  dev.off()
  
  # run time
  
  title <- sprintf('rtime_frac_%d_n%d_rho%2.0f_additive%s_sign%s_type%s.png', 
                   round(100*args$frac), args$n, 100*args$sqCovMat$rho, args$additive, args$sign, args$sqCovMat$type)
  png(file.path(plot_dir, title), width = 800, height = 600)
  matplot(s_arr, mean_rtime, type = 'b', pch = 16, lwd = 4,
          main=title,
          xlab=NA,
          ylab=NA)
  mtext(side = 1, cex = 1.5, line = 2, text = 's')
  mtext(side = 2, cex = 1.5, line = 2.5, text = 'Runtime in seconds')
  legend('bottomright', colnames(mean_acc), col=seq_len(method_count), fill=seq_len(method_count))
  dev.off()
  
  
  title <- sprintf('log_rtime_frac_%d_n%d_rho%2.0f_additive%s_sign%s_type%s.png', 
                   round(100*args$frac), args$n, 100*args$sqCovMat$rho, args$additive, args$sign, args$sqCovMat$type)
  png(file.path(plot_dir, title), width = 800, height = 600)
  matplot(s_arr, log10(mean_rtime), type = 'b', pch = 16, lwd = 4,
          main=title,
          xlab=NA,
          ylab=NA)
  mtext(side = 1, cex = 1.5, line = 2, text = 's')
  mtext(side = 2, cex = 1.5, line = 2.5, text = 'Runtime in seconds (log_10 scale)')
  legend('bottomright', colnames(mean_acc), col=seq_len(method_count), fill=seq_len(method_count))
  dev.off()
}

stopCluster(cl)


cat(sprintf("Total time = %3.3e\n", Sys.time()-start.time))

