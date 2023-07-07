########### data prep ########
######################################################
######################################################
## See: https://stephenslab.github.io/single-cell-topics/prepare_droplet.html
## Save the data as droplet.RData

######################################################
######################################################
########### run ebpmf #######
######################################################
######################################################
library(fastTopics)
library(ebpmf)
library(Matrix)
load('~/droplet.RData')

#counts = counts[,colSums(counts!=0)>10]
print(dim(counts))

res = ebpmf_log(counts,var_type = 'by_col',
                flash_control = list(ebnm.fn=c(ebnm::ebnm_point_exponential,ebnm::ebnm_point_exponential),
                                     factors_sign=1,loadings_sign = 1,Kmax=30),
                sigma2_control = list(return_sigma2_trace = T),
                verbose = TRUE,
                init_control = list(n_cores = 10,init_tol=1e-2),
                general_control = list(batch_size=1000,
                                       printevery=1,
                                       maxiter=100,
                                       save_fit_every = 5,
                                       save_fit_path='~/droplet_iteration_results/',
                                       save_fit_name=paste('ebpmf_droplet_full_nonnegLF'))
)

######################################################
######################################################
########## get plots ############
######################################################
######################################################
library(Matrix)
library(ggplot2)
library(dplyr)
library(ggrepel)

source('code/plot_factors.R')
source('code/plot_factors_general.R')
source('code/structure_plot.R')
source('code/get_loadings_order.R')

colors = c('#a6cee3',
           '#1f78b4',
           '#b2df8a',
           '#33a02c',
           '#fb9a99',
           '#e31a1c',
           '#fdbf6f',
           '#ff7f00',
           '#cab2d6',
           '#6a3d9a',
           '#ffff99',
           '#b15928',
           "#8B8A5F",
           "#D26C6D",
           "#F2F4ED"
)

dim(res$fit_flash$L.pm)
p = structure_plot_general(res$fit_flash$L.pm,res$fit_flash$F.pm,samples$tissue,
                           LD=T,remove_l0f0 = T,title = 'EBPMF',
                           n_samples = 2000,K=7,print_plot = T)

kset=11

ldf = my_ldf(res$fit_flash$L.pm[,-c(1,2)],res$fit_flash$F.pm[,-c(1,2)])
Lhat = ldf$l%*%diag(ldf$d)
rare_idx = which(samples$tissue %in% c("Goblet","Ionocyte", "Neuroendocrine","Tuft"))
Lhat = Lhat[rare_idx,1:kset]
Fhat = matrix(1,nrow=3,ncol=ncol(Lhat))
if(is.null(colnames(Lhat))){
  colnames(Lhat) <- paste0("k",1:ncol(Lhat))
}
fit_list     <- list(L = Lhat,F = Fhat)
class(fit_list) <- c("multinom_topic_model_fit", "list")
rare_cell_type = as.factor(as.character(samples$tissue[rare_idx]))
p = structure_plot(fit_list,grouping = rare_cell_type,
                   loadings_order = 'embed',
                   n = 2000,gap = 5,colors=colors,verbose=F) +
  labs(y = "loading",color = "dim",fill = "dim") +ggtitle("EBPMF, rare cell types")
