# input: fit, topics, grouping

# poisson2multinom
#
library(fastTopics)
library(ggplot2)
library(gridExtra)
structure_plot_general = function(Lhat,Fhat,grouping,title=NULL,
                                  loadings_order = 'embed',
                                  print_plot=TRUE,
                                  seed=12345,
                                  n_samples = NULL,
                                  gap=40,
                                  LD = TRUE,
                                  remove_l0f0 = TRUE,
                                  topic_model=FALSE,
                                  show_legend=TRUE,
                                  K = NULL,
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
                                    '#b15928')){
  set.seed(seed)
  #s       <- apply(Lhat,2,max)
  #Lhat    <-	t(t(Lhat) / s)

  if(is.null(n_samples)&all(loadings_order == "embed")){
    n_samples = 2000
  }

  if(!topic_model){
    if(remove_l0f0){
      ldf = my_ldf(Lhat[,-c(1,2)],Fhat[,-c(1,2)])
    }else{
      ldf = my_ldf(Lhat,Fhat)
    }
  }else{
    ldf = list(l=Lhat,d=1)
    LD=F
  }

  if(LD){
    Lhat = ldf$l%*%diag(ldf$d)
  }else{
    Lhat = ldf$l
  }
  if(!is.null(K)){
    Lhat = Lhat[,1:K]
    Fhat = Fhat[,1:K]
  }
  Fhat = matrix(1,nrow=3,ncol=ncol(Lhat))
  if(is.null(colnames(Lhat))){
    colnames(Lhat) <- paste0("k",1:ncol(Lhat))
  }
  fit_list     <- list(L = Lhat,F = Fhat)
  class(fit_list) <- c("multinom_topic_model_fit", "list")
  p <- structure_plot(fit_list,grouping = grouping,
                      loadings_order = loadings_order,
                      n = n_samples,gap = gap,colors=colors,verbose=F) +
    labs(y = "loading",color = "dim",fill = "dim") + ggtitle(title)
  if(!show_legend){
    p <- p + theme(legend.position="none")
  }
  if(print_plot){
    print(p)
  }
  return(p)
}

my_ldf = function(Lhat,Fhat){
  dl = apply(Lhat,2,norm,type='2')
  df = apply(Fhat,2,norm,type='2')
  d = df*df
  ord = order(d,decreasing = T)
  return(list(l = apply(Lhat,2,function(z){z/norm(z,'2')})[,ord],
              f = apply(Fhat,2,function(z){z/norm(z,'2')})[,ord],
              d = d[ord]))
}




