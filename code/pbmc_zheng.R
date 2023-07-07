################# data prep ###################
###############################################

# see https://stephenslab.github.io/single-cell-topics/prepare_purified_pbmc.html

######################################################
######################################################
################# fit ebpmf ################
######################################################
######################################################
library(fastTopics)
library(ebpmf)
library(Matrix)
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(glmpca)
library(scry)

library(peakRAM)
load('~/pbmc_purified.RData')
# filter counts using deviance method
devs <- scry::devianceFeatureSelection(t(counts))
dev_ranked_genes <- colnames(counts)[order(devs, decreasing = TRUE)]
topdev <- head(dev_ranked_genes, 3000)
counts = counts[,colnames(counts)%in%topdev]
counts = counts[,colSums(counts!=0)>10]

rm(devs)
gc()

peakRAM(res <- try(ebpmf_log(counts,
                              flash_control = list(ebnm.fn = c(ebnm::ebnm_point_exponential,ebnm::ebnm_point_exponential),
                                                   factors_sign=1,loadings_sign=1,Kmax=20),
                              sigma2_control = list(return_sigma2_trace=T),
                              verbose = T,
                              init_control = list(n_cores = 1,init_tol=1e-4),
                              general_control = list(garbage_collection_every=1,printevery=1,batch_size = 10000,
                                                     maxiter=100,save_fit_every=5,
                                                     save_fit_path='~/iteration_results/',
                                                     save_fit_name=paste('ebpmf_pbmc_3000gene_nonnegLF')))))

######################################################
######################################################
################# get plot ###################
###############################################

source('code/plot_factors.R')
source('code/plot_factors_general.R')
source('code/structure_plot.R')
source('code/get_loadings_order.R')
structure_plot_general(res$fit_flash$L.pm,res$fit_flash$F.pm,cell_names,LD=T,remove_l0f0 = T,title = "EBPMF",n_samples = 5000)
devs <- scry::devianceFeatureSelection(t(counts))
dev_ranked_genes <- colnames(counts)[order(devs, decreasing = TRUE)]
topdev <- head(dev_ranked_genes, 3000)
gene_idx = colnames(counts)%in%topdev
gene_names = genes$symbol[gene_idx]


cols = rep('grey80',length(gene_names))
cols[which(gene_names%in%gene_names[order(res$fit_flash$F.pm[,3],decreasing = T)][1:10])] = 'red'
plot(res$fit_flash$F.pm[,3],xlab = '',
     pch=20,ylab='',col=cols,
     axes = F)
axis(2,at=c(0,1,2,3,4,5),labels = c(0,1,2,3,4,5))
text(which(gene_names=='GNLY')+100,res$fit_flash$F.pm[which(gene_names=='GNLY'),3],label='GNLY',cex=0.8)
text(which(gene_names=='CCL5')+100,res$fit_flash$F.pm[which(gene_names=='CCL5'),3],label='CCL5',cex=0.8)
text(which(gene_names=='NKG7')+100,res$fit_flash$F.pm[which(gene_names=='NKG7'),3],label='NKG7',cex=0.8)
text(which(gene_names=='GZMB')+100,res$fit_flash$F.pm[which(gene_names=='GZMB'),3],label='GZMB',cex=0.8)
text(which(gene_names=='TYROBP')+100,res$fit_flash$F.pm[which(gene_names=='TYROBP'),3],label='TYROBP',cex=0.8)
text(which(gene_names=='FGFBP2'),res$fit_flash$F.pm[which(gene_names=='FGFBP2'),3]+0.2,label='FGFBP2',cex=0.8)
text(which(gene_names=='GZMA')+100,res$fit_flash$F.pm[which(gene_names=='GZMA'),3],label='GZMA',cex=0.8)
text(which(gene_names=='GZMH')+100,res$fit_flash$F.pm[which(gene_names=='GZMH'),3],label='GZMH',cex=0.8)
text(which(gene_names=='CST7')+100,res$fit_flash$F.pm[which(gene_names=='CST7'),3],label='CST7',cex=0.8)
text(which(gene_names=='CLIC3')+100,res$fit_flash$F.pm[which(gene_names=='CLIC3'),3],label='CLIC3',cex=0.8)




