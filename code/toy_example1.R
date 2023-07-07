library(ggplot2)
library(tidyverse)
library(purrr)
library(ebpmf)
library(glmpca)
library(PoissonPCA)
library(generalizedPCA)
my_ldf = function(Lhat,Fhat){
  dl = apply(Lhat,2,norm,type='2')
  df = apply(Fhat,2,norm,type='2')
  d = dl*df
  ord = order(d,decreasing = T)
  d = d[ord]
  l = apply(Lhat,2,function(z){z/norm(z,'2')})
  l = l[,ord]
  f = apply(Fhat,2,function(z){z/norm(z,'2')})
  f = f[,ord]
  return(list(l = l,
              f=f,
              d = d))
}

plot.factors.general = function(LL,
                                groups,
                                ylims=c(-1,1),
                                kset = NULL,
                                max.pt.size = 2,
                                title = NULL,
                                nonnegative = FALSE) {

  if(is.null(kset)){
    kset = 1:ncol(LL)
  }
  LL = LL[,kset]

  # # Re-normalize loadings so that factors are equally spread out.
  LL <- t(t(LL) / apply(abs(LL), 2, max))
  #
  # To make it easier to compare factors, flip them to make the largest
  #   loadings positive.
  flip <- 2 * (colSums(LL > 0.75) > colSums(LL < -0.75)) - 1
  LL <- t(t(LL) * flip)

  # Make the size of the point depend on how many of that type there are.
  sizes <- max.pt.size / sqrt(table(groups) / min(table(groups)))

  df <- reshape2::melt(LL, value.name = "loading")
  df$groups <- rep(as.factor(groups), length(kset))
  ggplot(df, aes(x = Var2, y = loading, color = groups)) +
    geom_jitter(position = position_jitter(height=0, width=0.4),
                size = rep(sizes[groups], length(kset))) +
    labs(title = title, x = NULL) +
    lims(y = ylims)
}

plot_factor_1by1_ggplot <- function(LL) {

  # Convert the matrix to a data frame
  LL_df <- as.data.frame(LL)

  # Add a row number column, which will serve as the x-axis
  LL_df$rn <- 1:nrow(LL_df)

  # Convert the data to a 'long' format
  LL_long <- LL_df %>%
    gather(key = "variable", value = "value", -rn)

  # Generate the plots
  ggplot(LL_long, aes(x = rn, y = value)) +
    geom_line(linewidth = 1) +
    facet_wrap(~ variable, scales = "free_y", ncol = 1, strip.position = "bottom") +
    theme(strip.text = element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank())
}

# generate data
set.seed(12345)
N = 300
p = 100
K = 3
Ftrue = matrix(0,nrow=p,ncol=K)
Ftrue[1:20,1] = 1+runif(20,0.5,1)
Ftrue[21:40,2] = -2+runif(20,0.5,1)
Ftrue[41:60,3] = 3+runif(20,0.5,1)
sigma2 = 0
Ltrue = matrix(rnorm(N*K), ncol=K)

Lambda = exp(tcrossprod(Ltrue,Ftrue) + matrix(rnorm(N*p,0,sqrt(sigma2)),nrow=N))
Y = matrix(rpois(N*p,Lambda),nrow=N,ncol=p)
Y = t(Y)
LF = list(L=Ltrue,F=Ftrue)
Ltrue = LF$F
Ftrue = LF$L
clusters = c(rep(1,20),rep(2,20),rep(3,20),rep(4,p-60))

ldf_true = my_ldf(Ltrue,Ftrue)
# plot(ldf_true$f[,1])
# plot(ldf_true$f[,2])
# plot(ldf_true$f[,3])
# ldf_true$d
plot.factors.general(ldf_true$l,clusters,ylims = c(-0.55,1.05),title = 'True Loadings')
#plot_factor_1by1_ggplot(ldf_true$l)

fit = ebpmf_log(Y,l0=0,f0=0,flash_control=list(fix_f0=T,fix_l0=T),
                var_type = 'constant',
                sigma2_control = list(return_sigma2_trace=T),
                general_control = list(maxiter=100,save_latent_M=T,conv_tol=1e-8))
plot(fit$elbo_trace,type='l',ylab = 'ELBO')
plot(fit$K_trace-2)
plot(fit$sigma2_trace,ylab=expression(sigma^2),xlab='iteration',type='l',col=2,lwd=2)
fit$sigma2



ldf_split = my_ldf(fit$fit_flash$L.pm[,-c(1,2)],fit$fit_flash$F.pm[,-c(1,2)])
plot.factors.general(ldf_split$l,clusters,ylims = c(-0.55,1.05),title = 'EBPMF')
#plot_factor_1by1_ggplot(ldf_split$l)

## fit GLMPCA


set.seed(12345)
n_rep = 10
smallest_dev = Inf
for(i in 1:n_rep){
  fit_glmpca = glmpca(Y,L=3,fam='poi',sz=rep(1,N))
  if(fit_glmpca$dev[length(fit_glmpca$dev)]<smallest_dev){
    smallest_dev = fit_glmpca$dev[length(fit_glmpca$dev)]
    fit_glmpca_best = fit_glmpca
  }
}

ldf_glmpca = my_ldf(fit_glmpca_best$loadings,fit_glmpca_best$factors)
plot.factors.general(ldf_glmpca$l,clusters,ylims = c(-0.55,1.05),title = 'GLM-PCA')

## fit Poisson PCA

fit_ppca=PoissonPCA::Poisson_Corrected_PCA(Y,3,transformation = 'log')
ldf_ppca = my_ldf(fit_ppca$scores,fit_ppca$loadings)
#plot(ldf_ppca$f[,1])
#plot(ldf_ppca$f[,2])
#plot(ldf_ppca$f[,3])
plot.factors.general(ldf_ppca$l,clusters,ylims = c(-0.55,1.05),title = 'Poisson PCA')

## fit generalized PCA, very slow
#
# fit_gpca = generalizedPCA(Y,3,family = 'poisson')
# ldf_gpca = my_ldf(fit_gpca$PCs,fit_gpca$U)
# plot.factors.general(ldf_gpca$l,clusters,ylims = c(-0.55,1.05),title = 'Generalized PCA')
#
