get_loadings_order = function(Lhat,Fhat,grouping,
                              seed=12345,
                              n_samples = 2000,
                              LD = TRUE,
                              remove_l0f0 = TRUE,
                              topic_model=FALSE){


  set.seed(seed)

  if(!topic_model){
    ldf = my_ldf(Lhat,Fhat)
  }else{
    ldf = list(l=Lhat,d=1)
  }

  if(LD){
    Lhat = ldf$l%*%diag(ldf$d)
  }else{
    Lhat = ldf$l
  }
  if(remove_l0f0){
    Lhat = Lhat[,-c(1,2),drop=F]
  }
  Fhat = matrix(1,nrow=3,ncol=ncol(Lhat))
  if(is.null(colnames(Lhat))){
    colnames(Lhat) <- paste0("k",1:ncol(Lhat))
  }
  fit_list     <- list(L = Lhat,F = Fhat)
  class(fit_list) <- c("multinom_topic_model_fit", "list")

  if(n_samples < nrow(fit_list$L)){
    rows <- sample(nrow(fit_list$L),n_samples)
    fit_list <- select_loadings(fit_list,rows)
    grouping <- grouping[rows,drop = FALSE]
  }

  loadings_order <- NULL
  for (group in levels(grouping)) {
    i <- which(grouping == group)
    if (length(i) > 0)
      y <- tsne_from_topics(select_loadings(fit_list, i),verbose=F)
    loadings_order <- c(loadings_order, i[order(y)])
  }
  loadings_order
}


