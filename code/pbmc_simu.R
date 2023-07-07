source('code/simu_func.R')
######################################################
#################run simulations ##################
######################################################
######################################################
# load real data fit results
fit = readRDS('data/ebpmf_fit_pbmc.rds')
simdata = sim_data_pmf_log(fit$fit_flash$L.pm,fit$fit_flash$F.pm,fit$sigma2,var_type = 'by_col',n_simu=5)
rm(fit)
gc()
res = simu_study_PMF(simdata,n_cores = 3,
                     method_list=c('flash','ebpmf','nmf'),
                     ebnm_function = ebnm::ebnm_point_exponential,
                     loadings_sign = 1,
                     factors_sign = 1,
                     Kmax=10,
                     var_type='by_col',
                     init_tol=1e-2,
                     maxiter=100,
                     tol=1e-5)
saveRDS(res,file='data/pbmc_simu_res.rds')


######################################################
################# get plots ##################
######################################################
######################################################

source('code/plot_factors.R')
source('code/plot_factors_general.R')
source('code/structure_plot.R')
source('code/get_loadings_order.R')

loadings_order = get_loadings_order(res$sim_data$Loading,res$sim_data$Factor,
                                    grouping = pbmc_facs$samples$subpop,n_samples = 5000)

plot0=structure_plot_general(res$sim_data$Loading,res$sim_data$Factor,
                             grouping =pbmc_facs$samples$subpop,title = 'True',print_plot = F,
                             loadings_order = loadings_order)

order_topic = function(ref,est){
  n = nrow(ref)
  p = ncol(est)
  ref = apply(ref,2,function(z){z/norm(z,'2')})
  est = apply(est,2,function(z){z/norm(z,'2')})
  dist_mat = matrix(nrow=p,ncol=p)
  for(i in 1:p){
    for(j in 1:p){
      dist_mat[i,j] = sqrt(mean((est[,i]-ref[,j])^2))
    }
  }
  order_output = c()
  for(i in 1:p){
    order_output[i] = which.min(dist_mat[,i])
    rm_idx = c()
    while(i > 1 & order_output[i]%in%order_output[-i]){
      rm_idx = c(rm_idx,order_output[i])
      temp = dist_mat[,i]
      temp[rm_idx] = Inf
      order_output[i] = which.min(temp)
    }
  }
  order_output
}
for(i in c(1,3,4,5)){
  order_nmf = order_topic(res$sim_data$Loading[,-c(1,2)],res$output[[i]]$fitted_model$nmf$W[,-1])
  plot1 = structure_plot_general((res$output[[i]]$fitted_model$nmf$W[,-1])[,order_nmf],
                                 (t(res$output[[i]]$fitted_model$nmf$H[-1,]))[,order_nmf],
                                 grouping=pbmc_facs$samples$subpop,
                                 title='log transformation+NMF(squared loss)',
                                 print_plot = F,
                                 loadings_order=loadings_order,
                                 remove_l0f0 = F)
  order_ebpmf = order_topic(res$sim_data$Loading,res$output[[i]]$fitted_model$ebpmf$fit_flash$L.pm)
  plot2 = structure_plot_general(res$output[[i]]$fitted_model$ebpmf$fit_flash$L.pm[,order_ebpmf],
                                 res$output[[i]]$fitted_model$ebpmf$fit_flash$F.pm[,order_ebpmf],
                                 grouping =pbmc_facs$samples$subpop,
                                 title='EBPMF',
                                 print_plot = F,
                                 loadings_order=loadings_order)
  order_flash = order_topic(res$sim_data$Loading[,-c(1,2)],res$output[[i]]$fitted_model$flash$L.pm[,-1])
  plot3 = structure_plot_general((res$output[[i]]$fitted_model$flash$L.pm[,-1])[,order_flash],
                                 (res$output[[i]]$fitted_model$flash$L.pm[,-1])[,order_flash],
                                 grouping =pbmc_facs$samples$subpop,
                                 title='log transformation+EBNMF',
                                 remove_l0f0 = F,
                                 print_plot = F,
                                 loadings_order=loadings_order)
  grid.arrange(plot0,plot2,plot3, plot1,nrow=4)
}


