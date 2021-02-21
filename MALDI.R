
# Data Import --------------------------------------------

#Read in Annotations
anno <- read.csv(" add path ...Metaspace_annotations_943peaks.csv", sep=";")

## MALDI data
library(Cardinal)
setwd("add path to imzML files")
data2 <- readImzML("imzML_folder")
image_data <- pixelData(data2) %>% as.data.frame()
image_data$Breaks=1
count=1
for(i in 2:nrow(image_data)){

  if(c(image_data$y[i]-image_data$y[i-1])>=2 | c(image_data$y[i]-image_data$y[i-1])<=c(-2) ){
    count=count+1
    image_data$Breaks[i]=count
    message(count)
  }else{image_data$Breaks[i]=count}

}

mz <- data2@featureData@mz
spectra <- as.matrix(spectraData(data2)[[1]])
rownames(spectra) <- mz
annotation <- data.frame(mz=as.numeric(rownames(spectra)), anno=0)


# Analysis All Data ------------------------------------------------------


#Batch Effectes

#Benchmark Batch effect removal 

spectra.batch <- limma::removeBatchEffect(spectra, image_data$Breaks)
spectra.batch <- sva::ComBat(spectra, image_data$Breaks)

#-> low performance 

# Dimensional reduction and Clustering:

#select 1000 most variable
var <- data.frame(row=rownames(spectra.batch), var=apply(spectra.batch,1,var)) %>% top_n(., 1000, var)
most.var.spectra <- spectra.batch[var$row, ]

#Cluster and DimRed

#PCA
pca <- irlba::prcomp_irlba(most.var.spectra, n=30)$rotation %>% as.data.frame()

pca.plot <- pca %>% mutate(barcodes=colnames(spectra.batch),samples=as.factor(image_data$Breaks))
ggplot(data=pca.plot, aes(x=PC1+PC2, y=PC3, color=samples))+
  geom_point(size=0.2, alpha=0.5)+
  SPATA::scale_color_add_on(variable = "character",clrp = "milo")+ 
  theme_classic()


#Graph Clustering

nearest <- RANN::nn2(pca, k = NN, treetype = "bd", searchtype = "priority")
nearest$nn.idx <- nearest$nn.idx
nearest$nn.dists <- nearest$nn.dists
edges = reshape::melt(t(nearest$nn.idx[, 1:NN]))
base::colnames(edges) = c("B", "A", "C")
edges = edges[, c("A", "B", "C")]
edges$B = edges$C
edges$C = 1
edges = base::unique(base::transform(edges, A = pmin(A, B), B = pmax(A, B)))
names(edges) <- c("V1", "V2", "weight")
edges$V1 <- base::rownames(pca)[edges$V1]
edges$V2 <- base::rownames(pca)[edges$V2]
g <- igraph::graph.data.frame(edges, directed = F)
graph.out = igraph::cluster_louvain(g)
clust.assign = base::factor(graph.out$membership, levels = base::sort(base::unique(graph.out$membership)))
Cluster_out=base::data.frame(ID=base::rownames(pca), Cluster=clust.assign, samples=image_data$Breaks); base::rownames(Cluster_out)=Cluster_out$ID

# Add cluster to plot df
pca.plot <- 
  pca %>% 
  mutate(barcodes=colnames(spectra.batch),
         samples=as.factor(image_data$Breaks),
         Cluster=clust.assign
  )

#UMAP
umap=umap::umap(pca)
UMAP <- base::data.frame(barcodes=base::rownames(pca),samples=image_data$Breaks,cluster=Cluster_out$Cluster,umap1=umap$layout[,1], umap2=umap$layout[,2])

#Run TSNE
TSNE=Rtsne::Rtsne(pca, perplexity=50)
TSNE <- base::data.frame(barcodes=base::rownames(pca),
                         samples=image_data$Breaks,
                         #cluster=Cluster_out$Cluster,
                         umap1=TSNE$Y[,1], umap2=TSNE$Y[,2])

#Plot Dim Red

library(ggplot2)
col_use <- colorRampPalette(confuns::clrp_milo)(23)

confuns::all_colorpanels()
ggplot(TSNE, aes(x=umap1, y=umap2, color=as.character(Cluster_out$samples)))+geom_point(size=0.1, alpha=0.5)+theme_classic()+xlim(-10,12)
ggplot(TSNE, aes(x=umap1, y=umap2, color=as.character(TSNE$samples)))+
  geom_point(size=0.2, alpha=0.5)+theme_classic()+
  scale_color_manual(values=col_use)


# Better Performance in Batch Effect removal using Monocle3

library(monocle3)

#Create a cds file

fdata <- image_data %>% tibble::rownames_to_column("row_ID") %>% mutate(sample=paste0("ID_S",Breaks,"_ID",row_ID))
expression_matrix=spectra.batch
class(expression_matrix)
dim(expression_matrix)
colnames(expression_matrix) <- fdata$sample
gene_metadata = data.frame(gene_short_name=rownames(expression_matrix));rownames(gene_metadata)=rownames(expression_matrix)
cell_metadata=data.frame(barcodes=fdata$sample, sample=paste0("sample_", fdata$Breaks)); rownames(cell_metadata)=fdata$sample
class(cell_metadata)
cell_metadata$sample=as.factor(cell_metadata$sample)
table(rownames(cell_metadata) %in% colnames(expression_matrix))

cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_metadata)
cds=cds[,Matrix::colSums(exprs(cds)) != 0]

#Process and Plot TSNE

cds <- estimate_size_factors(cds)
cds <- preprocess_cds(cds, preprocess_method="PCA", num_dim = 15)
cds <- align_cds(cds, alignment_group = "sample")
cds <- preprocess_cds(cds, preprocess_method="PCA", num_dim = 15)
cds <- reduce_dimension(cds,reduction_method="tSNE", cores=4)
cds <- reduce_dimension(cds)
cds <- cluster_cells(cds,cluster_method="louvain")
col_use <- colorRampPalette(confuns::clrp_milo)(24)
plot_cells(cds, reduction_method="tSNE", color_cells_by="sample", label_cell_groups=F, cell_size=0.5)+scale_color_manual(values=col_use)


# Seperate Patients and Create SPATA file -------------------------------------------------------

#-------------------------Annotate Peaks to Metabolites

#-> Run all samples for SPATA object

for(xx in 1:6){
  

samples_MALDI <- c("248_T", "259_T", "260_T", "262_T", "269_T", "275_T")
Seperate <- image_data %>%
  tibble::rownames_to_column("row_ID") %>%
  dplyr::filter(Breaks==xx) %>%
  dplyr::mutate(sample=paste0("ID_S",Breaks,"_ID",row_ID))
spectra_S1 <- as.data.frame(spectra.batch[,as.numeric(Seperate$row_ID)])
spectra_S1 <- na.omit(spectra_S1)
spectra_S1 <- spectra_S1[rowSums(spectra_S1)>0, ]




#Function to replicate rows with similar peak annotations
row_rep <- function(df, n) {
  df[rep(1:nrow(df), times = n),]
}

#Annotate
library(MALDIquant)
match <- match.closest(anno$mz, as.numeric(rownames(spectra_S1))) 

data_all_Spectra <- as.data.frame(do.call(rbind,pbapply::pblapply(1:nrow(anno), 
                                                                  function(i){
  MET=unlist(stringr::str_split(anno[i,]$moleculeIds, pattern=", "))
  data_out <-
       as.data.frame(spectra_S1[match[i], ,drop=F]) %>%
      'names<-'(colnames(spectra_S1)) %>%
      colSums() %>%
      as.data.frame() %>%
      t() %>%
      as.data.frame() %>%
      row_rep(., n=length(MET)) %>%
      mutate(HMDH=MET ) %>%
      select(HMDH, names(as.data.frame(spectra_S1[r.ow, ,drop=F])))
  dim(data_out)
 return(data_out)
})))
dim(data_all_Spectra)
data_all_Spectra <- data_all_Spectra[rowSums(data_all_Spectra[,2:ncol(data_all_Spectra)])>0, ]
dim(data_all_Spectra)

#Summarize Annotations
duplicated <- data_all_Spectra[duplicated(data_all_Spectra$HMDH),"HMDH"]
#remove=rownames((data_all_Spectra %>% filter(HMDH %in% duplicated)))
data_all <- as.data.frame(do.call(rbind, pbapply::pblapply(1:length(duplicated), function(z){
  out <- data.frame(HMDH=duplicated[z], t(as.data.frame(colMeans(data_all_Spectra %>% filter(HMDH==duplicated[z]) %>% select(!HMDH) ))) )
  rownames(out)=NULL
  return(out)
})))
data_all.dt <- rbind(data_all, data_all_Spectra %>% filter(!(HMDH %in% duplicated)))


rows <- data_all.dt$HMDH
data_all.dt$HMDH <- NULL
rownames(data_all.dt) <- rows
data_all <- data_all.dt %>% as.matrix()

#Coordinate File
S1 <- image_data %>% tibble::rownames_to_column("row_ID") %>% dplyr::filter(Breaks==xx) %>% mutate(sample=paste0("ID_S",Breaks,"_ID",row_ID))
coordinates <- data.frame(barcodes=S1$sample,sample=paste0("Sample_",S1$Breaks), x=as.numeric(S1$X3DPositionX), y=as.numeric(S1$X3DPositionY))

dim(pw_meta)
pw_meta.x <-  pw_meta %>% filter(gene %in% rownames(data_all))

library(SPATA)
S1_object <-  SPATA::initiateSpataObject_MALDI(coordinates=coordinates,
                                        intensity_matrix=data_all,
                                        geneSets = pw_meta.x,
                                        file_name = paste0("MALDI_S",samples_MALDI[xx],".RDS"))

}





















