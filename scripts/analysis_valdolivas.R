
#valdolivas analysis 
library(Seurat)
library(SeuratObject)


load("~/Desktop/transferdata/wanglong/GSMobjfinal.Robj")
load("~/Desktop/transferdata/wanglong/signfinal.Robj")
sign$EMR=NULL


# Custom subsetting function for Visium objects
custom_subset_visium <- function(seurat_object, cells) {
  # Subset the counts matrix
  subset_counts <- seurat_object@assays$Spatial@layers$counts[, cells, drop = FALSE]
  # Create a new Seurat object with the subsetted counts matrix and metadata
  new_seurat_object <- CreateSeuratObject(counts = subset_counts, meta.data = seurat_object@meta.data[cells, ])
  
  
  # Subset the image data if present
  if (!is.null(seurat_object@images)) {
    new_seurat_object@images <- seurat_object@images
    for (image_name in names(seurat_object@images)) {
      new_seurat_object@images[[image_name]]@boundaries$centroids@coords <- seurat_object@images[[image_name]]@boundaries$centroids@coords[cells, , drop = FALSE]
    }
  }
  
  # Subset the graph data if present
  if (!is.null(seurat_object@graphs)) {
    graph_object <- as.Graph(subset_counts)
    new_seurat_object@graphs <- list(graph = graph_object)
  }
  
  return(new_seurat_object)
}

library(Matrix)
allspatial<-CreateSeuratObject(counts = Matrix(cbind(GSM7058756_C1@assays$Spatial@layers$counts,
                                                     GSM7058757_C2@assays$Spatial@layers$counts,
                                                     GSM7058758_C3@assays$Spatial@layers$counts,
                                                     GSM7058759_C4@assays$Spatial@layers$counts),sparse = T),
                               meta.data = as.data.frame(rbind(GSM7058756_C1@meta.data,
                                                               GSM7058757_C2@meta.data,
                                                               GSM7058758_C3@meta.data,
                                                               GSM7058759_C4@meta.data)))


allspatial1<-allspatial@meta.data
pdf('allspatial_wl.pdf')
ggplot(allspatial1,aes(x=log10(Fibroblasts_C2L),y=Fibroblasts_C2L_percentage,colour = orig.ident))+geom_point()
dev.off()


spotstoconsider<-intersect(rownames(allspatial@meta.data)[allspatial$Fibroblasts_C2L>quantile(allspatial$Fibroblasts_C2L,.85)],
                           rownames(allspatial@meta.data)[allspatial$Fibroblasts_C2L_percentage>20])

GSM7058759_C4_fibroenriched<-custom_subset_visium(GSM7058759_C4,intersect(rownames(GSM7058759_C4@meta.data),spotstoconsider))
GSM7058756_C1_fibroenriched<-custom_subset_visium(GSM7058756_C1,intersect(rownames(GSM7058756_C1@meta.data),spotstoconsider))
GSM7058757_C2_fibroenriched<-custom_subset_visium(GSM7058757_C2,intersect(rownames(GSM7058757_C2@meta.data),spotstoconsider))
GSM7058758_C3_fibroenriched<-custom_subset_visium(GSM7058758_C3,intersect(rownames(GSM7058758_C3@meta.data),spotstoconsider))
rm(allspatial)

V573_1<-Load10X_Spatial('/home/carlo/Desktop/transferdata/valdolivas/SN048_A121573_Rep1')
V573_2<-Load10X_Spatial('/home/carlo/Desktop/transferdata/valdolivas/SN048_A121573_Rep2')
V371_1<-Load10X_Spatial('/home/carlo/Desktop/transferdata/valdolivas/SN048_A416371_Rep1')
V371_2<-Load10X_Spatial('/home/carlo/Desktop/transferdata/valdolivas/SN048_A416371_Rep2')
V838_1<-Load10X_Spatial('/home/carlo/Desktop/transferdata/valdolivas/SN84_A120838_Rep1')
V838_2<-Load10X_Spatial('/home/carlo/Desktop/transferdata/valdolivas/SN84_A120838_Rep2')
V763_1<-Load10X_Spatial('/home/carlo/Desktop/transferdata/valdolivas/SN123_A551763_Rep1')
V763_2<-Load10X_Spatial('/home/carlo/Desktop/transferdata/valdolivas/SN124_A551763_Rep2')
V688_1<-Load10X_Spatial('/home/carlo/Desktop/transferdata/valdolivas/SN123_A595688_Rep1')
V688_2<-Load10X_Spatial('/home/carlo/Desktop/transferdata/valdolivas/SN124_A595688_Rep2')
V015_1<-Load10X_Spatial('/home/carlo/Desktop/transferdata/valdolivas/SN123_A798015_Rep1')
V015_2<-Load10X_Spatial('/home/carlo/Desktop/transferdata/valdolivas/SN124_A798015_Rep2')
V797_1<-Load10X_Spatial('/home/carlo/Desktop/transferdata/valdolivas/SN123_A938797_Rep1_X')
V797_2<-Load10X_Spatial('/home/carlo/Desktop/transferdata/valdolivas/SN124_A938797_Rep2')

deconvKorean<-read.csv('/home/carlo/Desktop/transferdata/valdolivas/DeconvKoreanCohort/W_cell_density_q05.csv')
deconvBelgian<-read.csv('/home/carlo/Desktop/transferdata/valdolivas/DeconvBelgianCohort/W_cell_density_q05.csv')

colnames(deconvKorean)=gsub('q05_spot_factors','',colnames(deconvKorean))
colnames(deconvBelgian)=gsub('q05_spot_factors','',colnames(deconvBelgian))

id_korean <- c()
barcode_korean <- c()
for (i in 1:length(deconvKorean$spot_id)) {
  elements <- strsplit(deconvKorean$spot_id[i], '_')[[1]]
  if (elements[1] == "Count") {
    id_korean <- c(id_korean, paste(elements[1:4], collapse = '_'))
    barcode_korean <- c(barcode_korean, elements[5])
  } else {
    id_korean <- c(id_korean, paste(elements[1:3], collapse = '_'))
    barcode_korean <- c(barcode_korean, elements[4])
  }
}

id_belgian <- c()
barcode_belgian <- c()
for (i in 1:length(deconvBelgian$spot_id)) {
  elements <- strsplit(deconvBelgian$spot_id[i], '_')[[1]]
  if (elements[1] == "Count") {
    id_belgian <- c(id_belgian, paste(elements[1:4], collapse = '_'))
    barcode_belgian <- c(barcode_belgian, elements[5])
  } else {
    id_belgian <- c(id_belgian, paste(elements[1:3], collapse = '_'))
    barcode_belgian <- c(barcode_belgian, elements[4])
  }
}

deconvKorean$ST.exp<-id_korean
deconvBelgian$ST.exp<-id_belgian

deconvKorean$ST.barcode<-barcode_korean
deconvBelgian$ST.barcode<-barcode_belgian


library(dplyr)
v573_1_MD<-deconvKorean%>%
  filter(ST.exp=='Count_SN048_A121573_Rep1')
v573_1_MD<-v573_1_MD[,2:39]
v573_1_MD<-v573_1_MD[v573_1_MD$ST.barcode%in%rownames(V573_1@meta.data),]
V573_1<-V573_1[,v573_1_MD$ST.barcode]
rownames(v573_1_MD)=rownames(V573_1@meta.data)
colnames(v573_1_MD)=gsub('$','_KoreanDeconv',colnames(v573_1_MD))
V573_1@meta.data<-cbind(V573_1@meta.data,v573_1_MD)


library(dplyr)
v573_2_MD<-deconvKorean%>%
  filter(ST.exp=='Count_SN048_A121573_Rep2')
v573_2_MD<-v573_2_MD[,2:39]
v573_2_MD<-v573_2_MD[v573_2_MD$ST.barcode%in%rownames(V573_2@meta.data),]
V573_2<-V573_2[,v573_2_MD$ST.barcode]
rownames(v573_2_MD)=rownames(V573_2@meta.data)
colnames(v573_2_MD)=gsub('$','_KoreanDeconv',colnames(v573_2_MD))
V573_2@meta.data<-cbind(V573_2@meta.data,v573_2_MD)


library(dplyr)
v371_1_MD<-deconvKorean%>%
  filter(ST.exp=='Count_SN048_A416371_Rep1')
v371_1_MD<-v371_1_MD[,2:39]
v371_1_MD<-v371_1_MD[v371_1_MD$ST.barcode%in%rownames(V371_1@meta.data),]
V371_1<-V371_1[,v371_1_MD$ST.barcode]
rownames(v371_1_MD)=rownames(V371_1@meta.data)
colnames(v371_1_MD)=gsub('$','_KoreanDeconv',colnames(v371_1_MD))
V371_1@meta.data<-cbind(V371_1@meta.data,v371_1_MD)


library(dplyr)
v371_2_MD<-deconvKorean%>%
  filter(ST.exp=='Count_SN048_A416371_Rep2')
v371_2_MD<-v371_2_MD[,2:39]
v371_2_MD<-v371_1_MD[v371_2_MD$ST.barcode%in%rownames(V371_2@meta.data),]
V371_2<-V371_1[,v371_2_MD$ST.barcode]
rownames(v371_2_MD)=rownames(V371_2@meta.data)
colnames(v371_2_MD)=gsub('$','_KoreanDeconv',colnames(v371_2_MD))
V371_2@meta.data<-cbind(V371_2@meta.data,v371_2_MD)


library(dplyr)
v763_1_MD<-deconvKorean%>%
  filter(ST.exp=='Count_SN123_A551763_Rep1')
v763_1_MD<-v763_1_MD[,2:39]
v763_1_MD<-v763_1_MD[v763_1_MD$ST.barcode%in%rownames(V763_1@meta.data),]
V763_1<-V763_1[,v763_1_MD$ST.barcode]
rownames(v763_1_MD)=rownames(V763_1@meta.data)
colnames(v763_1_MD)=gsub('$','_KoreanDeconv',colnames(v763_1_MD))
V763_1@meta.data<-cbind(V763_1@meta.data,v763_1_MD)


library(dplyr)
v763_2_MD<-deconvKorean%>%
  filter(ST.exp=='Count_SN124_A551763_Rep2')
v763_2_MD<-v763_2_MD[,2:39]
v763_2_MD<-v763_2_MD[v763_2_MD$ST.barcode%in%rownames(V763_2@meta.data),]
V763_2<-V763_2[,v763_2_MD$ST.barcode]
rownames(v763_2_MD)=rownames(V763_2@meta.data)
colnames(v763_2_MD)=gsub('$','_KoreanDeconv',colnames(v763_2_MD))
V763_2@meta.data<-cbind(V763_2@meta.data,v763_2_MD)


library(dplyr)
v688_1_MD<-deconvKorean%>%
  filter(ST.exp=='Count_SN123_A595688_Rep1')
v688_1_MD<-v688_1_MD[,2:39]
v688_1_MD<-v688_1_MD[v688_1_MD$ST.barcode%in%rownames(V688_1@meta.data),]
V688_1<-V688_1[,v688_1_MD$ST.barcode]
rownames(v688_1_MD)=rownames(V688_1@meta.data)
colnames(v688_1_MD)=gsub('$','_KoreanDeconv',colnames(v688_1_MD))
V688_1@meta.data<-cbind(V688_1@meta.data,v688_1_MD)


library(dplyr)
v688_2_MD<-deconvKorean%>%
  filter(ST.exp=='Count_SN124_A595688_Rep2')
v688_2_MD<-v688_2_MD[,2:39]
v688_2_MD<-v688_2_MD[v688_2_MD$ST.barcode%in%rownames(V688_2@meta.data),]
V688_2<-V688_2[,v688_2_MD$ST.barcode]
rownames(v688_2_MD)=rownames(V688_2@meta.data)
colnames(v688_2_MD)=gsub('$','_KoreanDeconv',colnames(v688_2_MD))
V688_2@meta.data<-cbind(V688_2@meta.data,v688_2_MD)


library(dplyr)
v015_1_MD<-deconvKorean%>%
  filter(ST.exp=='Count_SN123_A798015_Rep1')
v015_1_MD<-v015_1_MD[,2:39]
v015_1_MD<-v015_1_MD[v015_1_MD$ST.barcode%in%rownames(V015_1@meta.data),]
V015_1<-V015_1[,v015_1_MD$ST.barcode]
rownames(v015_1_MD)=rownames(V015_1@meta.data)
colnames(v015_1_MD)=gsub('$','_KoreanDeconv',colnames(v015_1_MD))
V015_1@meta.data<-cbind(V015_1@meta.data,v015_1_MD)


library(dplyr)
v015_2_MD<-deconvKorean%>%
  filter(ST.exp=='Count_SN124_A798015_Rep2')
v015_2_MD<-v015_2_MD[,2:39]
v015_2_MD<-v015_2_MD[v015_2_MD$ST.barcode%in%rownames(V015_2@meta.data),]
V015_2<-V015_2[,v015_2_MD$ST.barcode]
rownames(v015_2_MD)=rownames(V015_2@meta.data)
colnames(v015_2_MD)=gsub('$','_KoreanDeconv',colnames(v015_2_MD))
V015_2@meta.data<-cbind(V015_2@meta.data,v015_2_MD)

library(dplyr)
v797_1_MD<-deconvKorean%>%
  filter(ST.exp=='SN123_A938797_Rep1')
v797_1_MD<-v797_1_MD[,2:39]
v797_1_MD<-v797_1_MD[v797_1_MD$ST.barcode%in%rownames(V797_1@meta.data),]
V797_1<-V797_1[,v797_1_MD$ST.barcode]
rownames(v797_1_MD)=rownames(V797_1@meta.data)
colnames(v797_1_MD)=gsub('$','_KoreanDeconv',colnames(v797_1_MD))
V797_1@meta.data<-cbind(V797_1@meta.data,v797_1_MD)


library(dplyr)
v797_2_MD<-deconvKorean%>%
  filter(ST.exp=='Count_SN124_A938797_Rep2')
v797_2_MD<-v797_2_MD[,2:39]
v797_2_MD<-v797_2_MD[v797_2_MD$ST.barcode%in%rownames(V797_2@meta.data),]
V797_2<-V797_2[,v797_2_MD$ST.barcode]
rownames(v797_2_MD)=rownames(V797_2@meta.data)
colnames(v797_2_MD)=gsub('$','_KoreanDeconv',colnames(v797_2_MD))
V797_2@meta.data<-cbind(V797_2@meta.data,v797_2_MD)


library(dplyr)
v838_1_MD<-deconvKorean%>%
  filter(ST.exp=='Count_SN84_A120838_Rep1')
v838_1_MD<-v838_1_MD[,2:39]
v838_1_MD<-v838_1_MD[v838_1_MD$ST.barcode%in%rownames(V838_1@meta.data),]
V838_1<-V838_1[,v838_1_MD$ST.barcode]
rownames(v838_1_MD)=rownames(V838_1@meta.data)
colnames(v838_1_MD)=gsub('$','_KoreanDeconv',colnames(v838_1_MD))
V838_1@meta.data<-cbind(V838_1@meta.data,v838_1_MD)

library(dplyr)
v838_2_MD<-deconvKorean%>%
  filter(ST.exp=='Count_SN84_A120838_Rep2')
v838_2_MD<-v838_2_MD[,2:39]
v838_2_MD<-v838_2_MD[v838_2_MD$ST.barcode%in%rownames(V838_2@meta.data),]
V838_2<-V838_2[,v838_2_MD$ST.barcode]
rownames(v838_2_MD)=rownames(V838_2@meta.data)
colnames(v838_2_MD)=gsub('$','_KoreanDeconv',colnames(v838_2_MD))
V838_2@meta.data<-cbind(V838_2@meta.data,v838_2_MD)


{
  library(dplyr)
  v573_1_MD<-deconvBelgian%>%
    filter(ST.exp=='Count_SN048_A121573_Rep1')
  v573_1_MD<-v573_1_MD[,2:43]
  v573_1_MD<-v573_1_MD[v573_1_MD$ST.barcode%in%rownames(V573_1@meta.data),]
  V573_1<-V573_1[,v573_1_MD$ST.barcode]
  rownames(v573_1_MD)=rownames(V573_1@meta.data)
  colnames(v573_1_MD)=gsub('$','_BelgianDeconv',colnames(v573_1_MD))
  V573_1@meta.data<-cbind(V573_1@meta.data,v573_1_MD)
  
  
  library(dplyr)
  v573_2_MD<-deconvBelgian%>%
    filter(ST.exp=='Count_SN048_A121573_Rep2')
  v573_2_MD<-v573_2_MD[,2:43]
  v573_2_MD<-v573_2_MD[v573_2_MD$ST.barcode%in%rownames(V573_2@meta.data),]
  V573_2<-V573_2[,v573_2_MD$ST.barcode]
  rownames(v573_2_MD)=rownames(V573_2@meta.data)
  colnames(v573_2_MD)=gsub('$','_BelgianDeconv',colnames(v573_2_MD))
  V573_2@meta.data<-cbind(V573_2@meta.data,v573_2_MD)
  
  
  library(dplyr)
  v371_1_MD<-deconvBelgian%>%
    filter(ST.exp=='Count_SN048_A416371_Rep1')
  v371_1_MD<-v371_1_MD[,2:43]
  v371_1_MD<-v371_1_MD[v371_1_MD$ST.barcode%in%rownames(V371_1@meta.data),]
  V371_1<-V371_1[,v371_1_MD$ST.barcode]
  rownames(v371_1_MD)=rownames(V371_1@meta.data)
  colnames(v371_1_MD)=gsub('$','_BelgianDeconv',colnames(v371_1_MD))
  V371_1@meta.data<-cbind(V371_1@meta.data,v371_1_MD)
  
  
  library(dplyr)
  v371_2_MD<-deconvBelgian%>%
    filter(ST.exp=='Count_SN048_A416371_Rep2')
  v371_2_MD<-v371_2_MD[,2:43]
  v371_2_MD<-v371_1_MD[v371_2_MD$ST.barcode%in%rownames(V371_2@meta.data),]
  V371_2<-V371_1[,v371_2_MD$ST.barcode]
  rownames(v371_2_MD)=rownames(V371_2@meta.data)
  colnames(v371_2_MD)=gsub('$','_BelgianDeconv',colnames(v371_2_MD))
  V371_2@meta.data<-cbind(V371_2@meta.data,v371_2_MD)
  
  
  library(dplyr)
  v763_1_MD<-deconvBelgian%>%
    filter(ST.exp=='Count_SN123_A551763_Rep1')
  v763_1_MD<-v763_1_MD[,2:43]
  v763_1_MD<-v763_1_MD[v763_1_MD$ST.barcode%in%rownames(V763_1@meta.data),]
  V763_1<-V763_1[,v763_1_MD$ST.barcode]
  rownames(v763_1_MD)=rownames(V763_1@meta.data)
  colnames(v763_1_MD)=gsub('$','_BelgianDeconv',colnames(v763_1_MD))
  V763_1@meta.data<-cbind(V763_1@meta.data,v763_1_MD)
  
  
  library(dplyr)
  v763_2_MD<-deconvBelgian%>%
    filter(ST.exp=='Count_SN124_A551763_Rep2')
  v763_2_MD<-v763_2_MD[,2:43]
  v763_2_MD<-v763_2_MD[v763_2_MD$ST.barcode%in%rownames(V763_2@meta.data),]
  V763_2<-V763_2[,v763_2_MD$ST.barcode]
  rownames(v763_2_MD)=rownames(V763_2@meta.data)
  colnames(v763_2_MD)=gsub('$','_BelgianDeconv',colnames(v763_2_MD))
  V763_2@meta.data<-cbind(V763_2@meta.data,v763_2_MD)
  
  
  library(dplyr)
  v688_1_MD<-deconvBelgian%>%
    filter(ST.exp=='Count_SN123_A595688_Rep1')
  v688_1_MD<-v688_1_MD[,2:43]
  v688_1_MD<-v688_1_MD[v688_1_MD$ST.barcode%in%rownames(V688_1@meta.data),]
  V688_1<-V688_1[,v688_1_MD$ST.barcode]
  rownames(v688_1_MD)=rownames(V688_1@meta.data)
  colnames(v688_1_MD)=gsub('$','_BelgianDeconv',colnames(v688_1_MD))
  V688_1@meta.data<-cbind(V688_1@meta.data,v688_1_MD)
  
  
  library(dplyr)
  v688_2_MD<-deconvBelgian%>%
    filter(ST.exp=='Count_SN124_A595688_Rep2')
  v688_2_MD<-v688_2_MD[,2:43]
  v688_2_MD<-v688_2_MD[v688_2_MD$ST.barcode%in%rownames(V688_2@meta.data),]
  V688_2<-V688_2[,v688_2_MD$ST.barcode]
  rownames(v688_2_MD)=rownames(V688_2@meta.data)
  colnames(v688_2_MD)=gsub('$','_BelgianDeconv',colnames(v688_2_MD))
  V688_2@meta.data<-cbind(V688_2@meta.data,v688_2_MD)
  
  
  library(dplyr)
  v015_1_MD<-deconvBelgian%>%
    filter(ST.exp=='Count_SN123_A798015_Rep1')
  v015_1_MD<-v015_1_MD[,2:43]
  v015_1_MD<-v015_1_MD[v015_1_MD$ST.barcode%in%rownames(V015_1@meta.data),]
  V015_1<-V015_1[,v015_1_MD$ST.barcode]
  rownames(v015_1_MD)=rownames(V015_1@meta.data)
  colnames(v015_1_MD)=gsub('$','_BelgianDeconv',colnames(v015_1_MD))
  V015_1@meta.data<-cbind(V015_1@meta.data,v015_1_MD)
  
  
  library(dplyr)
  v015_2_MD<-deconvBelgian%>%
    filter(ST.exp=='Count_SN124_A798015_Rep2')
  v015_2_MD<-v015_2_MD[,2:43]
  v015_2_MD<-v015_2_MD[v015_2_MD$ST.barcode%in%rownames(V015_2@meta.data),]
  V015_2<-V015_2[,v015_2_MD$ST.barcode]
  rownames(v015_2_MD)=rownames(V015_2@meta.data)
  colnames(v015_2_MD)=gsub('$','_BelgianDeconv',colnames(v015_2_MD))
  V015_2@meta.data<-cbind(V015_2@meta.data,v015_2_MD)
  
  library(dplyr)
  v797_1_MD<-deconvBelgian%>%
    filter(ST.exp=='SN123_A938797_Rep1')
  v797_1_MD<-v797_1_MD[,2:43]
  v797_1_MD<-v797_1_MD[v797_1_MD$ST.barcode%in%rownames(V797_1@meta.data),]
  V797_1<-V797_1[,v797_1_MD$ST.barcode]
  rownames(v797_1_MD)=rownames(V797_1@meta.data)
  colnames(v797_1_MD)=gsub('$','_BelgianDeconv',colnames(v797_1_MD))
  V797_1@meta.data<-cbind(V797_1@meta.data,v797_1_MD)
  
  
  library(dplyr)
  v797_2_MD<-deconvBelgian%>%
    filter(ST.exp=='Count_SN124_A938797_Rep2')
  v797_2_MD<-v797_2_MD[,2:43]
  v797_2_MD<-v797_2_MD[v797_2_MD$ST.barcode%in%rownames(V797_2@meta.data),]
  V797_2<-V797_2[,v797_2_MD$ST.barcode]
  rownames(v797_2_MD)=rownames(V797_2@meta.data)
  colnames(v797_2_MD)=gsub('$','_BelgianDeconv',colnames(v797_2_MD))
  V797_2@meta.data<-cbind(V797_2@meta.data,v797_2_MD)
  
  
  library(dplyr)
  v838_1_MD<-deconvBelgian%>%
    filter(ST.exp=='Count_SN84_A120838_Rep1')
  v838_1_MD<-v838_1_MD[,2:43]
  v838_1_MD<-v838_1_MD[v838_1_MD$ST.barcode%in%rownames(V838_1@meta.data),]
  V838_1<-V838_1[,v838_1_MD$ST.barcode]
  rownames(v838_1_MD)=rownames(V838_1@meta.data)
  colnames(v838_1_MD)=gsub('$','_BelgianDeconv',colnames(v838_1_MD))
  V838_1@meta.data<-cbind(V838_1@meta.data,v838_1_MD)
  
  library(dplyr)
  v838_2_MD<-deconvBelgian%>%
    filter(ST.exp=='Count_SN84_A120838_Rep2')
  v838_2_MD<-v838_2_MD[,2:43]
  v838_2_MD<-v838_2_MD[v838_2_MD$ST.barcode%in%rownames(V838_2@meta.data),]
  V838_2<-V838_2[,v838_2_MD$ST.barcode]
  rownames(v838_2_MD)=rownames(V838_2@meta.data)
  colnames(v838_2_MD)=gsub('$','_BelgianDeconv',colnames(v838_2_MD))
  V838_2@meta.data<-cbind(V838_2@meta.data,v838_2_MD)
}




#to do: perform estimate of the combined estimates using random forests or xgboost, using a target the mean
{
# Install necessary packages if not already installed
if (!require(randomForest)) install.packages("randomForest")
if (!require(xgboost)) install.packages("xgboost")
if (!require(caret)) install.packages("caret")
if (!require(dplyr)) install.packages("dplyr")

# Load libraries
library(randomForest)
library(xgboost)
library(caret)
library(dplyr)


# Example data for multiple datasets
# Assuming you have data frames for each dataset with cell types as columns and samples as rows
datasets <- list(
  data.frame(CellType1 = c(0.1, 0.3, 0.5), CellType2 = c(0.2, 0.4, 0.6)),
  data.frame(CellType1 = c(0.15, 0.25, 0.55), CellType2 = c(0.25, 0.35, 0.65))
)

# Combine estimates into a single data frame
combined_estimates <- do.call(rbind, lapply(datasets, function(df) {
  df %>% rename_all(function(x) paste0(x, "_", which(datasets == df)))
}))

# Prepare data for training
set.seed(123) # For reproducibility

# Split data into training and testing sets
train_index <- createDataPartition(rowMeans(combined_estimates), p = 0.8, list = FALSE)
train_data <- combined_estimates[train_index, ]
test_data <- combined_estimates[-train_index, ]

# Define a function to prepare the data for model training
prepare_data <- function(data) {
  features <- data
  labels <- rowMeans(data)
  list(features = features, labels = labels)
}

# Prepare training and testing datasets
train <- prepare_data(train_data)
test <- prepare_data(test_data)

# Define a control function for model tuning
control <- trainControl(method = "cv", number = 5)

# Train Random Forest model
rf_grid <- expand.grid(mtry = 1:ncol(train$features))
rf_model <- train(train$features, train$labels,
                  method = "rf",
                  trControl = control,
                  tuneGrid = rf_grid)

# Train XGBoost model
xgb_grid <- expand.grid(nrounds = 100,
                        eta = c(0.01, 0.1, 0.3),
                        max_depth = c(3, 6, 9),
                        gamma = 0,
                        colsample_bytree = 0.8,
                        min_child_weight = 1,
                        subsample = 0.8)
xgb_model <- train(train$features, train$labels,
                   method = "xgbTree",
                   trControl = control,
                   tuneGrid = xgb_grid)

# Predict unified estimates on the test set
rf_predictions <- predict(rf_model, test$features)
xgb_predictions <- predict(xgb_model, test$features)

# Evaluate model performance
rf_rmse <- RMSE(rf_predictions, test$labels)
xgb_rmse <- RMSE(xgb_predictions, test$labels)

cat("Random Forest RMSE:", rf_rmse, "\n")
cat("XGBoost RMSE:", xgb_rmse, "\n")

# Create a data frame with the unified estimates
unified_estimates <- data.frame(
  Test_Labels = test$labels,
  RF_Predictions = rf_predictions,
  XGB_Predictions = xgb_predictions
)

print(unified_estimates)

}

#for now, use the deconvKorean data


process_C2L_Korean<-function(seurat_object,df=FALSE){
  if(df==FALSE){
  mt <- seurat_object@meta.data
  cell_type_columns <- mt[, stringr::str_detect(pattern = "^(?!ST\\.exp_)(?!ST\\.barcode_).+_KoreanDeconv$", string = colnames(mt))]
  row_totals <- rowSums(cell_type_columns)
  cell_type_percentages <- sweep(cell_type_columns, 1, row_totals, "/") * 100
  colnames(cell_type_percentages) <- paste0(colnames(cell_type_percentages), "_percentageKorean")
  mt <- cbind(mt, cell_type_percentages)
  seurat_object@meta.data <- mt
  return(seurat_object)
  }else{
  cell_type_columns <- seurat_object[, stringr::str_detect(pattern = "^(?!spot.id_)(?!ST\\.exp_)(?!ST\\.barcode_).+_KoreanDeconv$", string = colnames(seurat_object))]
  row_totals <- rowSums(cell_type_columns)
  cell_type_percentages <- sweep(cell_type_columns, 1, row_totals, "/") * 100
  colnames(cell_type_percentages) <- paste0(colnames(cell_type_percentages), "_percentageKorean")
  seurat_object <- cbind(seurat_object, cell_type_percentages)
  return(seurat_object)
  }
}

V015_1<-process_C2L_Korean(V015_1)
V015_2<-process_C2L_Korean(V015_2)

V371_1<-process_C2L_Korean(V371_1)
V371_2<-process_C2L_Korean(V371_2)

V573_1<-process_C2L_Korean(V573_1)
V573_2<-process_C2L_Korean(V573_2)

V688_1<-process_C2L_Korean(V688_1)
V688_2<-process_C2L_Korean(V688_2)

V763_1<-process_C2L_Korean(V763_1)
V763_2<-process_C2L_Korean(V763_2)

V797_1<-process_C2L_Korean(V797_1)
V797_2<-process_C2L_Korean(V797_2)

V838_1<-process_C2L_Korean(V838_1)
V838_2<-process_C2L_Korean(V838_2)

colnames(deconvKorean)=gsub('$','_KoreanDeconv',colnames(deconvKorean))
deconvKorean<-process_C2L_Korean(deconvKorean,df = T)

  
pdf('allspatial_valdolivas.pdf',width = 20,height = 10)
ggplot(deconvKorean,aes(x=log10(Myofibroblasts_KoreanDeconv+Stromal.1_KoreanDeconv+Stromal.2_KoreanDeconv+Stromal.3_KoreanDeconv),y=Myofibroblasts_KoreanDeconv_percentageKorean,colour = ST.exp_KoreanDeconv))+geom_point(alpha=0.5)+theme_classic()+
  geom_vline(xintercept = quantile(deconvKorean$Myofibroblasts_KoreanDeconv+deconvKorean$Stromal.1_KoreanDeconv+deconvKorean$Stromal.2_KoreanDeconv+deconvKorean$Stromal.3_KoreanDeconv,0.8),color='red',linetype='dashed')+
  geom_hline(yintercept = 20,color='red',linetype='dashed')
dev.off()


V015_1<-NormalizeData(V015_1)
V015_2<-NormalizeData(V015_2)
V371_1<-NormalizeData(V371_1)
V371_2<-NormalizeData(V371_2)
V573_1<-NormalizeData(V573_1)
V573_2<-NormalizeData(V573_2)
V688_1<-NormalizeData(V688_1)
V688_2<-NormalizeData(V688_2)
V763_1<-NormalizeData(V763_1)
V763_2<-NormalizeData(V763_2)
V797_1<-NormalizeData(V797_1)
V797_2<-NormalizeData(V797_2)
V838_1<-NormalizeData(V838_1)
V838_2<-NormalizeData(V838_2)

calculate_corrected_scores <- function(seurat_obj) {
  for (name in names(sign)) {
    genes <- sign[[name]]
    valid_genes <- genes[genes %in% rownames(seurat_obj@assays$Spatial@layers$data)]
    
    if (length(valid_genes) > 0) {
      mean.exp <- colMeans(seurat_obj@assays$Spatial@layers$data[valid_genes, ], na.rm = TRUE)
      
      if (name == "EMR") {
        mean.exp <- (mean.exp)
      } else {
        mean.exp <- (mean.exp)
      }
      
      score_column <- paste0(name, ".score_GM")
      seurat_obj@meta.data[[score_column]] <- mean.exp
    } else {
      warning(paste("No valid genes found for signature:", name))
    }
  }
  return(seurat_obj)
}


rownames(V015_1@assays$Spatial@features)->rownames(V015_1@assays$Spatial@layers$counts)
rownames(V015_1@assays$Spatial@features)->rownames(V015_1@assays$Spatial@layers$data)
rownames(V015_2@assays$Spatial@features)->rownames(V015_2@assays$Spatial@layers$counts)
rownames(V015_2@assays$Spatial@features)->rownames(V015_2@assays$Spatial@layers$data)
rownames(V371_1@assays$Spatial@features)->rownames(V371_1@assays$Spatial@layers$counts)
rownames(V371_1@assays$Spatial@features)->rownames(V371_1@assays$Spatial@layers$data)
rownames(V371_2@assays$Spatial@features)->rownames(V371_2@assays$Spatial@layers$counts)
rownames(V371_2@assays$Spatial@features)->rownames(V371_2@assays$Spatial@layers$data)
rownames(V573_1@assays$Spatial@features)->rownames(V573_1@assays$Spatial@layers$counts)
rownames(V573_1@assays$Spatial@features)->rownames(V573_1@assays$Spatial@layers$data)
rownames(V573_2@assays$Spatial@features)->rownames(V573_2@assays$Spatial@layers$counts)
rownames(V573_2@assays$Spatial@features)->rownames(V573_2@assays$Spatial@layers$data)
rownames(V688_1@assays$Spatial@features)->rownames(V688_1@assays$Spatial@layers$counts)
rownames(V688_1@assays$Spatial@features)->rownames(V688_1@assays$Spatial@layers$data)
rownames(V688_2@assays$Spatial@features)->rownames(V688_2@assays$Spatial@layers$counts)
rownames(V688_2@assays$Spatial@features)->rownames(V688_2@assays$Spatial@layers$data)
rownames(V763_1@assays$Spatial@features)->rownames(V763_1@assays$Spatial@layers$counts)
rownames(V763_1@assays$Spatial@features)->rownames(V763_1@assays$Spatial@layers$data)
rownames(V763_2@assays$Spatial@features)->rownames(V763_2@assays$Spatial@layers$counts)
rownames(V763_2@assays$Spatial@features)->rownames(V763_2@assays$Spatial@layers$data)
rownames(V797_1@assays$Spatial@features)->rownames(V797_1@assays$Spatial@layers$counts)
rownames(V797_1@assays$Spatial@features)->rownames(V797_1@assays$Spatial@layers$data)
rownames(V797_2@assays$Spatial@features)->rownames(V797_2@assays$Spatial@layers$counts)
rownames(V797_2@assays$Spatial@features)->rownames(V797_2@assays$Spatial@layers$data)
rownames(V838_1@assays$Spatial@features)->rownames(V838_1@assays$Spatial@layers$counts)
rownames(V838_1@assays$Spatial@features)->rownames(V838_1@assays$Spatial@layers$data)
rownames(V838_2@assays$Spatial@features)->rownames(V838_2@assays$Spatial@layers$counts)
rownames(V838_2@assays$Spatial@features)->rownames(V838_2@assays$Spatial@layers$data)




V015_1<-calculate_corrected_scores(V015_1)
V015_2<-calculate_corrected_scores(V015_2)

V371_1<-calculate_corrected_scores(V371_1)
V371_2<-calculate_corrected_scores(V371_2)

V573_1<-calculate_corrected_scores(V573_1)
V573_2<-calculate_corrected_scores(V573_2)

V688_1<-calculate_corrected_scores(V688_1)
V688_2<-calculate_corrected_scores(V688_2)

V763_1<-calculate_corrected_scores(V763_1)
V763_2<-calculate_corrected_scores(V763_2)

V797_1<-calculate_corrected_scores(V797_1)
V797_2<-calculate_corrected_scores(V797_2)

V838_1<-calculate_corrected_scores(V838_1)
V838_2<-calculate_corrected_scores(V838_2)




pdf('spatialplots_valdolivas.pdf',width = 50,height = 160)
ggarrange(ncol=2,nrow=14,
  SpatialFeaturePlot(V015_1,features = 'EMR.score_GM',alpha = 0),
  SpatialFeaturePlot(V015_1,features = 'mrCAF.score_GM',alpha = 0),
  SpatialFeaturePlot(V015_2,features = 'EMR.score_GM',alpha = 0),
  SpatialFeaturePlot(V015_2,features = 'mrCAF.score_GM',alpha = 0),
  SpatialFeaturePlot(V371_1,features = 'EMR.score_GM',alpha = 0),
  SpatialFeaturePlot(V371_1,features = 'mrCAF.score_GM',alpha = 0),
  SpatialFeaturePlot(V371_2,features = 'EMR.score_GM',alpha = 0),
  SpatialFeaturePlot(V371_2,features = 'mrCAF.score_GM',alpha = 0),
  SpatialFeaturePlot(V573_1,features = 'EMR.score_GM',alpha = 0),
  SpatialFeaturePlot(V573_1,features = 'mrCAF.score_GM',alpha = 0),
  SpatialFeaturePlot(V573_2,features = 'EMR.score_GM',alpha = 0),
  SpatialFeaturePlot(V573_2,features = 'mrCAF.score_GM',alpha = 0),
  SpatialFeaturePlot(V688_1,features = 'EMR.score_GM',alpha = 0),
  SpatialFeaturePlot(V688_1,features = 'mrCAF.score_GM',alpha = 0),
  SpatialFeaturePlot(V688_2,features = 'EMR.score_GM',alpha = 0),
  SpatialFeaturePlot(V688_2,features = 'mrCAF.score_GM',alpha = 0),
  SpatialFeaturePlot(V763_1,features = 'EMR.score_GM',alpha = 0),
  SpatialFeaturePlot(V763_1,features = 'mrCAF.score_GM',alpha = 0),
  SpatialFeaturePlot(V763_2,features = 'EMR.score_GM',alpha = 0),
  SpatialFeaturePlot(V763_2,features = 'mrCAF.score_GM',alpha = 0),
  SpatialFeaturePlot(V797_1,features = 'EMR.score_GM',alpha = 0),
  SpatialFeaturePlot(V797_1,features = 'mrCAF.score_GM',alpha = 0),
  SpatialFeaturePlot(V797_2,features = 'EMR.score_GM',alpha = 0),
  SpatialFeaturePlot(V797_2,features = 'mrCAF.score_GM',alpha = 0),
  SpatialFeaturePlot(V838_1,features = 'EMR.score_GM',alpha = 0),
  SpatialFeaturePlot(V838_1,features = 'mrCAF.score_GM',alpha = 0),
  SpatialFeaturePlot(V838_2,features = 'EMR.score_GM',alpha = 0),
  SpatialFeaturePlot(V838_2,features = 'mrCAF.score_GM',alpha = 0)
)
dev.off()



pdf('EMR_CMS_mrCAF_Fibro_Spatial_valdolivas.pdf',width = 500,height = 300)
ggarrange(ncol=10,nrow=14,
          SpatialFeaturePlot(V015_1,features = 'EMR.score_GM',image.alpha = 0.2),
          SpatialFeaturePlot(V015_1,features = 'CMS1_KoreanDeconv',image.alpha=0.2),
          SpatialFeaturePlot(V015_1,features = 'CMS2_KoreanDeconv',image.alpha=0.2),
          SpatialFeaturePlot(V015_1,features = 'CMS3_KoreanDeconv',image.alpha=0.2),
          SpatialFeaturePlot(V015_1,features = 'CMS4_KoreanDeconv',image.alpha=0.2),
          SpatialFeaturePlot(V015_1,features = 'mrCAF.score_GM',image.alpha = 0.2),
          SpatialFeaturePlot(V015_1,features = 'Myofibroblasts_KoreanDeconv',image.alpha = 0.2),
          SpatialFeaturePlot(V015_1,features = 'Stromal.1_KoreanDeconv',image.alpha = 0.2),
          SpatialFeaturePlot(V015_1,features = 'Stromal.2_KoreanDeconv',image.alpha = 0.2),
          SpatialFeaturePlot(V015_1,features = 'Stromal.3_KoreanDeconv',image.alpha = 0.2),
          SpatialFeaturePlot(V015_2,features = 'EMR.score_GM',image.alpha = 0.2),
          SpatialFeaturePlot(V015_2,features = 'CMS1_KoreanDeconv',image.alpha=0.2),
          SpatialFeaturePlot(V015_2,features = 'CMS2_KoreanDeconv',image.alpha=0.2),
          SpatialFeaturePlot(V015_2,features = 'CMS3_KoreanDeconv',image.alpha=0.2),
          SpatialFeaturePlot(V015_2,features = 'CMS4_KoreanDeconv',image.alpha=0.2),
          SpatialFeaturePlot(V015_2,features = 'mrCAF.score_GM',image.alpha = 0.2),
          SpatialFeaturePlot(V015_2,features = 'Myofibroblasts_KoreanDeconv',image.alpha = 0.2),
          SpatialFeaturePlot(V015_2,features = 'Stromal.1_KoreanDeconv',image.alpha = 0.2),
          SpatialFeaturePlot(V015_2,features = 'Stromal.2_KoreanDeconv',image.alpha = 0.2),
          SpatialFeaturePlot(V015_2,features = 'Stromal.3_KoreanDeconv',image.alpha = 0.2),
          SpatialFeaturePlot(V371_1,features = 'EMR.score_GM',image.alpha = 0.2),
          SpatialFeaturePlot(V371_1,features = 'CMS1_KoreanDeconv',image.alpha=0.2),
          SpatialFeaturePlot(V371_1,features = 'CMS2_KoreanDeconv',image.alpha=0.2),
          SpatialFeaturePlot(V371_1,features = 'CMS3_KoreanDeconv',image.alpha=0.2),
          SpatialFeaturePlot(V371_1,features = 'CMS4_KoreanDeconv',image.alpha=0.2),
          SpatialFeaturePlot(V371_1,features = 'mrCAF.score_GM',image.alpha = 0.2),
          SpatialFeaturePlot(V371_1,features = 'Myofibroblasts_KoreanDeconv',image.alpha = 0.2),
          SpatialFeaturePlot(V371_1,features = 'Stromal.1_KoreanDeconv',image.alpha = 0.2),
          SpatialFeaturePlot(V371_1,features = 'Stromal.2_KoreanDeconv',image.alpha = 0.2),
          SpatialFeaturePlot(V371_1,features = 'Stromal.3_KoreanDeconv',image.alpha = 0.2),
          SpatialFeaturePlot(V371_2,features = 'EMR.score_GM',image.alpha = 0.2),
          SpatialFeaturePlot(V371_2,features = 'CMS1_KoreanDeconv',image.alpha=0.2),
          SpatialFeaturePlot(V371_2,features = 'CMS2_KoreanDeconv',image.alpha=0.2),
          SpatialFeaturePlot(V371_2,features = 'CMS3_KoreanDeconv',image.alpha=0.2),
          SpatialFeaturePlot(V371_2,features = 'CMS4_KoreanDeconv',image.alpha=0.2),
          SpatialFeaturePlot(V371_2,features = 'mrCAF.score_GM',image.alpha = 0.2),
          SpatialFeaturePlot(V371_2,features = 'Myofibroblasts_KoreanDeconv',image.alpha = 0.2),
          SpatialFeaturePlot(V371_2,features = 'Stromal.1_KoreanDeconv',image.alpha = 0.2),
          SpatialFeaturePlot(V371_2,features = 'Stromal.2_KoreanDeconv',image.alpha = 0.2),
          SpatialFeaturePlot(V371_2,features = 'Stromal.3_KoreanDeconv',image.alpha = 0.2),
          SpatialFeaturePlot(V573_1,features = 'EMR.score_GM',image.alpha = 0.2),
          SpatialFeaturePlot(V573_1,features = 'CMS1_KoreanDeconv',image.alpha=0.2),
          SpatialFeaturePlot(V573_1,features = 'CMS2_KoreanDeconv',image.alpha=0.2),
          SpatialFeaturePlot(V573_1,features = 'CMS3_KoreanDeconv',image.alpha=0.2),
          SpatialFeaturePlot(V573_1,features = 'CMS4_KoreanDeconv',image.alpha=0.2),
          SpatialFeaturePlot(V573_1,features = 'mrCAF.score_GM',image.alpha = 0.2),
          SpatialFeaturePlot(V573_1,features = 'Myofibroblasts_KoreanDeconv',image.alpha = 0.2),
          SpatialFeaturePlot(V573_1,features = 'Stromal.1_KoreanDeconv',image.alpha = 0.2),
          SpatialFeaturePlot(V573_1,features = 'Stromal.2_KoreanDeconv',image.alpha = 0.2),
          SpatialFeaturePlot(V573_1,features = 'Stromal.3_KoreanDeconv',image.alpha = 0.2),
          SpatialFeaturePlot(V573_2,features = 'EMR.score_GM',image.alpha = 0.2),
          SpatialFeaturePlot(V573_2,features = 'CMS1_KoreanDeconv',image.alpha=0.2),
          SpatialFeaturePlot(V573_2,features = 'CMS2_KoreanDeconv',image.alpha=0.2),
          SpatialFeaturePlot(V573_2,features = 'CMS3_KoreanDeconv',image.alpha=0.2),
          SpatialFeaturePlot(V573_2,features = 'CMS4_KoreanDeconv',image.alpha=0.2),
          SpatialFeaturePlot(V573_2,features = 'mrCAF.score_GM',image.alpha = 0.2),
          SpatialFeaturePlot(V573_2,features = 'Myofibroblasts_KoreanDeconv',image.alpha = 0.2),
          SpatialFeaturePlot(V573_2,features = 'Stromal.1_KoreanDeconv',image.alpha = 0.2),
          SpatialFeaturePlot(V573_2,features = 'Stromal.2_KoreanDeconv',image.alpha = 0.2),
          SpatialFeaturePlot(V573_2,features = 'Stromal.3_KoreanDeconv',image.alpha = 0.2),
          SpatialFeaturePlot(V688_1,features = 'EMR.score_GM',image.alpha = 0.2),
          SpatialFeaturePlot(V688_1,features = 'CMS1_KoreanDeconv',image.alpha=0.2),
          SpatialFeaturePlot(V688_1,features = 'CMS2_KoreanDeconv',image.alpha=0.2),
          SpatialFeaturePlot(V688_1,features = 'CMS3_KoreanDeconv',image.alpha=0.2),
          SpatialFeaturePlot(V688_1,features = 'CMS4_KoreanDeconv',image.alpha=0.2),
          SpatialFeaturePlot(V688_1,features = 'mrCAF.score_GM',image.alpha = 0.2),
          SpatialFeaturePlot(V688_1,features = 'Myofibroblasts_KoreanDeconv',image.alpha = 0.2),
          SpatialFeaturePlot(V688_1,features = 'Stromal.1_KoreanDeconv',image.alpha = 0.2),
          SpatialFeaturePlot(V688_1,features = 'Stromal.2_KoreanDeconv',image.alpha = 0.2),
          SpatialFeaturePlot(V688_1,features = 'Stromal.3_KoreanDeconv',image.alpha = 0.2),
          SpatialFeaturePlot(V688_2,features = 'EMR.score_GM',image.alpha = 0.2),
          SpatialFeaturePlot(V688_2,features = 'CMS1_KoreanDeconv',image.alpha=0.2),
          SpatialFeaturePlot(V688_2,features = 'CMS2_KoreanDeconv',image.alpha=0.2),
          SpatialFeaturePlot(V688_2,features = 'CMS3_KoreanDeconv',image.alpha=0.2),
          SpatialFeaturePlot(V688_2,features = 'CMS4_KoreanDeconv',image.alpha=0.2),
          SpatialFeaturePlot(V688_2,features = 'mrCAF.score_GM',image.alpha = 0.2),
          SpatialFeaturePlot(V688_2,features = 'Myofibroblasts_KoreanDeconv',image.alpha = 0.2),
          SpatialFeaturePlot(V688_2,features = 'Stromal.1_KoreanDeconv',image.alpha = 0.2),
          SpatialFeaturePlot(V688_2,features = 'Stromal.2_KoreanDeconv',image.alpha = 0.2),
          SpatialFeaturePlot(V688_2,features = 'Stromal.3_KoreanDeconv',image.alpha = 0.2),
          SpatialFeaturePlot(V763_1,features = 'EMR.score_GM',image.alpha = 0.2),
          SpatialFeaturePlot(V763_1,features = 'CMS1_KoreanDeconv',image.alpha=0.2),
          SpatialFeaturePlot(V763_1,features = 'CMS2_KoreanDeconv',image.alpha=0.2),
          SpatialFeaturePlot(V763_1,features = 'CMS3_KoreanDeconv',image.alpha=0.2),
          SpatialFeaturePlot(V763_1,features = 'CMS4_KoreanDeconv',image.alpha=0.2),
          SpatialFeaturePlot(V763_1,features = 'mrCAF.score_GM',image.alpha = 0.2),
          SpatialFeaturePlot(V763_1,features = 'Myofibroblasts_KoreanDeconv',image.alpha = 0.2),
          SpatialFeaturePlot(V763_1,features = 'Stromal.1_KoreanDeconv',image.alpha = 0.2),
          SpatialFeaturePlot(V763_1,features = 'Stromal.2_KoreanDeconv',image.alpha = 0.2),
          SpatialFeaturePlot(V763_1,features = 'Stromal.3_KoreanDeconv',image.alpha = 0.2),
          SpatialFeaturePlot(V763_2,features = 'EMR.score_GM',image.alpha = 0.2),
          SpatialFeaturePlot(V763_2,features = 'CMS1_KoreanDeconv',image.alpha=0.2),
          SpatialFeaturePlot(V763_2,features = 'CMS2_KoreanDeconv',image.alpha=0.2),
          SpatialFeaturePlot(V763_2,features = 'CMS3_KoreanDeconv',image.alpha=0.2),
          SpatialFeaturePlot(V763_2,features = 'CMS4_KoreanDeconv',image.alpha=0.2),
          SpatialFeaturePlot(V763_2,features = 'mrCAF.score_GM',image.alpha = 0.2),
          SpatialFeaturePlot(V763_2,features = 'Myofibroblasts_KoreanDeconv',image.alpha = 0.2),
          SpatialFeaturePlot(V763_2,features = 'Stromal.1_KoreanDeconv',image.alpha = 0.2),
          SpatialFeaturePlot(V763_2,features = 'Stromal.2_KoreanDeconv',image.alpha = 0.2),
          SpatialFeaturePlot(V763_2,features = 'Stromal.3_KoreanDeconv',image.alpha = 0.2),
          SpatialFeaturePlot(V797_1,features = 'EMR.score_GM',image.alpha = 0.2),
          SpatialFeaturePlot(V797_1,features = 'CMS1_KoreanDeconv',image.alpha=0.2),
          SpatialFeaturePlot(V797_1,features = 'CMS2_KoreanDeconv',image.alpha=0.2),
          SpatialFeaturePlot(V797_1,features = 'CMS3_KoreanDeconv',image.alpha=0.2),
          SpatialFeaturePlot(V797_1,features = 'CMS4_KoreanDeconv',image.alpha=0.2),
          SpatialFeaturePlot(V797_1,features = 'mrCAF.score_GM',image.alpha = 0.2),
          SpatialFeaturePlot(V797_1,features = 'Myofibroblasts_KoreanDeconv',image.alpha = 0.2),
          SpatialFeaturePlot(V797_1,features = 'Stromal.1_KoreanDeconv',image.alpha = 0.2),
          SpatialFeaturePlot(V797_1,features = 'Stromal.2_KoreanDeconv',image.alpha = 0.2),
          SpatialFeaturePlot(V797_1,features = 'Stromal.3_KoreanDeconv',image.alpha = 0.2),
          SpatialFeaturePlot(V797_2,features = 'EMR.score_GM',image.alpha = 0.2),
          SpatialFeaturePlot(V797_2,features = 'CMS1_KoreanDeconv',image.alpha=0.2),
          SpatialFeaturePlot(V797_2,features = 'CMS2_KoreanDeconv',image.alpha=0.2),
          SpatialFeaturePlot(V797_2,features = 'CMS3_KoreanDeconv',image.alpha=0.2),
          SpatialFeaturePlot(V797_2,features = 'CMS4_KoreanDeconv',image.alpha=0.2),
          SpatialFeaturePlot(V797_2,features = 'mrCAF.score_GM',image.alpha = 0.2),
          SpatialFeaturePlot(V797_2,features = 'Myofibroblasts_KoreanDeconv',image.alpha = 0.2),
          SpatialFeaturePlot(V797_2,features = 'Stromal.1_KoreanDeconv',image.alpha = 0.2),
          SpatialFeaturePlot(V797_2,features = 'Stromal.2_KoreanDeconv',image.alpha = 0.2),
          SpatialFeaturePlot(V797_2,features = 'Stromal.3_KoreanDeconv',image.alpha = 0.2),
          SpatialFeaturePlot(V838_1,features = 'EMR.score_GM',image.alpha = 0.2),
          SpatialFeaturePlot(V838_1,features = 'CMS1_KoreanDeconv',image.alpha=0.2),
          SpatialFeaturePlot(V838_1,features = 'CMS2_KoreanDeconv',image.alpha=0.2),
          SpatialFeaturePlot(V838_1,features = 'CMS3_KoreanDeconv',image.alpha=0.2),
          SpatialFeaturePlot(V838_1,features = 'CMS4_KoreanDeconv',image.alpha=0.2),
          SpatialFeaturePlot(V838_1,features = 'mrCAF.score_GM',image.alpha = 0.2),
          SpatialFeaturePlot(V838_1,features = 'Myofibroblasts_KoreanDeconv',image.alpha = 0.2),
          SpatialFeaturePlot(V838_1,features = 'Stromal.1_KoreanDeconv',image.alpha = 0.2),
          SpatialFeaturePlot(V838_1,features = 'Stromal.2_KoreanDeconv',image.alpha = 0.2),
          SpatialFeaturePlot(V838_1,features = 'Stromal.3_KoreanDeconv',image.alpha = 0.2),
          SpatialFeaturePlot(V838_2,features = 'EMR.score_GM',image.alpha = 0.2),
          SpatialFeaturePlot(V838_2,features = 'CMS1_KoreanDeconv',image.alpha=0.2),
          SpatialFeaturePlot(V838_2,features = 'CMS2_KoreanDeconv',image.alpha=0.2),
          SpatialFeaturePlot(V838_2,features = 'CMS3_KoreanDeconv',image.alpha=0.2),
          SpatialFeaturePlot(V838_2,features = 'CMS4_KoreanDeconv',image.alpha=0.2),
          SpatialFeaturePlot(V838_2,features = 'mrCAF.score_GM',image.alpha = 0.2),
          SpatialFeaturePlot(V838_2,features = 'Myofibroblasts_KoreanDeconv',image.alpha = 0.2),
          SpatialFeaturePlot(V838_2,features = 'Stromal.1_KoreanDeconv',image.alpha = 0.2),
          SpatialFeaturePlot(V838_2,features = 'Stromal.2_KoreanDeconv',image.alpha = 0.2),
          SpatialFeaturePlot(V838_2,features = 'Stromal.3_KoreanDeconv',image.alpha = 0.2)
)
dev.off()


V015_1<-NormalizeData(V015_1)
V015_2<-NormalizeData(V015_2)
V371_1<-NormalizeData(V371_1)
V371_2<-NormalizeData(V371_2)
V573_1<-NormalizeData(V573_1)
V573_2<-NormalizeData(V573_2)
V688_1<-NormalizeData(V688_1)
V688_2<-NormalizeData(V688_2)
V763_1<-NormalizeData(V763_1)
V763_2<-NormalizeData(V763_2)
V797_1<-NormalizeData(V797_1)
V797_2<-NormalizeData(V797_2)
V838_1<-NormalizeData(V838_1)
V838_2<-NormalizeData(V838_2)


#optionally investigate soglie considerando tutto insieme 

{
library(Matrix)
allspatial_all<-CreateSeuratObject(counts = Matrix(cbind(GSM7058756_C1@assays$Spatial@layers$counts,
                                                     GSM7058757_C2@assays$Spatial@layers$counts,
                                                     GSM7058758_C3@assays$Spatial@layers$counts,
                                                     GSM7058759_C4@assays$Spatial@layers$counts,
                                                     V015_1@assays$Spatial@layers$counts,
                                                     V015_2@assays$Spatial@layers$counts,
                                                     V371_1@assays$Spatial@layers$counts,
                                                     V371_2@assays$Spatial@layers$counts,
                                                     V573_1@assays$Spatial@layers$counts,
                                                     V573_2@assays$Spatial@layers$counts,
                                                     V688_1@assays$Spatial@layers$counts,
                                                     V688_2@assays$Spatial@layers$counts,
                                                     V763_1@assays$Spatial@layers$counts,
                                                     V763_2@assays$Spatial@layers$counts,
                                                     V797_1@assays$Spatial@layers$counts,
                                                     V797_2@assays$Spatial@layers$counts,
                                                     V838_1@assays$Spatial@layers$counts,
                                                     V838_2@assays$Spatial@layers$counts),sparse = T),
                               meta.data = as.data.frame(rbind(GSM7058756_C1@meta.data,
                                                               GSM7058757_C2@meta.data,
                                                               GSM7058758_C3@meta.data,
                                                               GSM7058759_C4@meta.data,
                                                               V015_1@meta.data,
                                                               V015_2@meta.data,
                                                               V371_1@meta.data,
                                                               V371_2@meta.data,
                                                               V573_1@meta.data,
                                                               V573_2@meta.data,
                                                               V688_1@meta.data,
                                                               V688_2@meta.data,
                                                               V763_1@meta.data,
                                                               V763_2@meta.data,
                                                               V797_1@meta.data,
                                                               V797_2@meta.data,
                                                               V838_1@meta.data,
                                                               V838_2@meta.data)))


spotstoconsider<-intersect(rownames(allspatial@meta.data)[allspatial$Fibroblasts_C2L>quantile(allspatial$Fibroblasts_C2L,.80)],
                           rownames(allspatial@meta.data)[allspatial$Fibroblasts_C2L_percentage>quantile(allspatial$Fibroblasts_C2L_percentage,.80)])

}


quantile(deconvKorean$Myofibroblasts_KoreanDeconv_percentageKorean+deconvKorean$Stromal.1_KoreanDeconv_percentageKorean+deconvKorean$Stromal.2_KoreanDeconv_percentageKorean+deconvKorean$Stromal.3_KoreanDeconv_percentageKorean,0.8)

#20% is ok for the percentage, like wanglong
deconvKorean1<-deconvKorean[deconvKorean$Myofibroblasts_KoreanDeconv+deconvKorean$Stromal.1_KoreanDeconv+deconvKorean$Stromal.2_KoreanDeconv+deconvKorean$Stromal.3_KoreanDeconv>quantile(deconvKorean$Myofibroblasts_KoreanDeconv+deconvKorean$Stromal.1_KoreanDeconv+deconvKorean$Stromal.2_KoreanDeconv+deconvKorean$Stromal.3_KoreanDeconv,0.85),]
deconvKorean1<-deconvKorean1[deconvKorean1$Myofibroblasts_KoreanDeconv_percentageKorean+deconvKorean1$Stromal.1_KoreanDeconv_percentageKorean+deconvKorean1$Stromal.2_KoreanDeconv_percentageKorean+deconvKorean1$Stromal.3_KoreanDeconv_percentageKorean>0.2,]


deconvKorean1<-deconvKorean1%>%
  nest_by(ST.exp_KoreanDeconv)


V573_1_fibro<-V573_1[,deconvKorean1$data[[1]]$ST.barcode_KoreanDeconv]
V573_2_fibro<-V573_2[,deconvKorean1$data[[2]]$ST.barcode_KoreanDeconv]


V371_1_fibro<-V371_1[,deconvKorean1$data[[3]]$ST.barcode_KoreanDeconv]
V371_2_fibro<-V371_2[,deconvKorean1$data[[4]]$ST.barcode_KoreanDeconv]


V763_1_fibro<-V763_1[,deconvKorean1$data[[5]]$ST.barcode_KoreanDeconv]
V763_2_fibro<-V763_2[,deconvKorean1$data[[8]]$ST.barcode_KoreanDeconv]


V688_1_fibro<-V688_1[,deconvKorean1$data[[6]]$ST.barcode_KoreanDeconv]
V688_2_fibro<-V688_2[,deconvKorean1$data[[9]]$ST.barcode_KoreanDeconv]


V015_1_fibro<-V015_1[,deconvKorean1$data[[7]]$ST.barcode_KoreanDeconv]
V015_2_fibro<-V015_2[,deconvKorean1$data[[10]]$ST.barcode_KoreanDeconv]


V797_1_fibro<-V797_1[,deconvKorean1$data[[14]]$ST.barcode_KoreanDeconv]
V797_2_fibro<-V797_2[,deconvKorean1$data[[11]]$ST.barcode_KoreanDeconv]


V838_1_fibro<-V838_1[,deconvKorean1$data[[12]]$ST.barcode_KoreanDeconv]
V838_2_fibro<-V838_2[,deconvKorean1$data[[13]]$ST.barcode_KoreanDeconv]



library(reticulate)
use_python('/home/carlo/.local/share/r-miniconda/envs/giotto_env/bin/python', required = TRUE)

library(circlize)
library(data.table)
process_PAGEanalysis1 <- function(seurat_object, signatures_all, selected_signatures,only_fibro=T) {
  if(only_fibro==T){signatures_all <- signatures_all[!names(signatures_all) %in% c("CCM", "EMRM")]}
  raw_exprs <- seurat_object@assays$Spatial@layers$counts
  rownames(raw_exprs) <- rownames(seurat_object@assays$Spatial@features@.Data)
  colnames(raw_exprs) <- rownames(seurat_object@meta.data)
  spatial_locs <- as.data.table(GetTissueCoordinates(seurat_object)[,1:2])
  colnames(spatial_locs) <- c("x", "y")
  
  myGiottoObj <- createGiottoObject(raw_exprs = raw_exprs, spatial_locs = spatial_locs)
  myGiottoObj <- normalizeGiotto(gobject = myGiottoObj)
  
  # Create signature matrix for initial PAGE analysis
  all_signatures <- names(signatures_all)
  signature_matrix_complete <- makeSignMatrixPAGE(
    sign_names = all_signatures,
    sign_list = lapply(all_signatures, function(sign) signatures_all[[sign]])
  )
  
  # Run initial PAGE enrichment analysis on all signatures
  myGiottoObj_initial <- runPAGEEnrich(
    gobject = myGiottoObj,
    sign_matrix = signature_matrix_complete,
    min_overlap_genes = 2,
    output_enrichment = c("original", "zscore")
  )
  
  # Run PAGE enrichment analysis with p-values
  myGiottoObj_pval <- runPAGEEnrich(
    gobject = myGiottoObj,
    sign_matrix = signature_matrix_complete,
    min_overlap_genes = 2,
    p_value = TRUE,
    reverse_log_scale = FALSE
  )
  
  # Extract p-value results
  pval_df <- as.data.frame(myGiottoObj_pval@spatial_enrichment$PAGE)
  zscore_df<-as.data.frame(myGiottoObj_initial@spatial_enrichment$PAGE)
  
  zscore_df <- as.data.frame(zscore_df)
  colnames(zscore_df) <- paste0(colnames(zscore_df), "PAGE")
  seurat_object@meta.data <- cbind(seurat_object@meta.data, zscore_df)
  
  # Identify significantly enriched spots - transformed values corresponds to pval < 0.05
  significant_spots <- pval_df > 1.301
  
  # Store significant enrichment information in Seurat metadata
  significant_df <- as.data.frame(significant_spots)
  colnames(significant_df) <- paste0(colnames(significant_df), "PAGE_significant")
  seurat_object@meta.data <- cbind(seurat_object@meta.data, significant_df)
  
  return(seurat_object)
  
}


library(Giotto)
V573_1_fibro <- process_PAGEanalysis1(seurat_object = V573_1_fibro, signatures_all = sign)
V573_2_fibro <- process_PAGEanalysis1(seurat_object = V573_2_fibro, signatures_all = sign)
V371_1_fibro <- process_PAGEanalysis1(seurat_object = V371_1_fibro, signatures_all = sign)
V371_2_fibro <- process_PAGEanalysis1(seurat_object = V371_2_fibro, signatures_all = sign) #!!
V763_1_fibro <- process_PAGEanalysis1(seurat_object = V763_1_fibro, signatures_all = sign) #!!
V763_2_fibro <- process_PAGEanalysis1(seurat_object = V763_2_fibro, signatures_all = sign)
V688_1_fibro <- process_PAGEanalysis1(seurat_object = V688_1_fibro, signatures_all = sign) #!!
V688_2_fibro <- process_PAGEanalysis1(seurat_object = V688_2_fibro, signatures_all = sign) #!!
V015_1_fibro <- process_PAGEanalysis1(seurat_object = V015_1_fibro, signatures_all = sign) #!!
V015_2_fibro <- process_PAGEanalysis1(seurat_object = V015_2_fibro, signatures_all = sign) #!!
V797_1_fibro <- process_PAGEanalysis1(seurat_object = V797_1_fibro, signatures_all = sign) #!!
V797_2_fibro <- process_PAGEanalysis1(seurat_object = V797_2_fibro, signatures_all = sign) #!!
V838_1_fibro <- process_PAGEanalysis1(seurat_object = V838_1_fibro, signatures_all = sign) #!!
V838_2_fibro <- process_PAGEanalysis1(seurat_object = V838_2_fibro, signatures_all = sign) #!!


GSM7058756_C1_fibroenriched@images$slice1@boundaries$centroids@cells=rownames(GSM7058756_C1_fibroenriched@meta.data)
GSM7058757_C2_fibroenriched@images$slice1@boundaries$centroids@cells=rownames(GSM7058757_C2_fibroenriched@meta.data)
GSM7058758_C3_fibroenriched@images$slice1@boundaries$centroids@cells=rownames(GSM7058758_C3_fibroenriched@meta.data)
GSM7058759_C4_fibroenriched@images$slice1@boundaries$centroids@cells=rownames(GSM7058759_C4_fibroenriched@meta.data)


process_PAGEanalysis1 <- function(seurat_object, signatures_all, selected_signatures,only_fibro=T) {
    if(only_fibro==T){signatures_all <- signatures_all[!names(signatures_all) %in% c("CCM", "EMRM")]}
    raw_exprs <- seurat_object@assays$RNA@layers$counts
    rownames(raw_exprs) <- rownames(seurat_object@assays$RNA@features@.Data)
    colnames(raw_exprs) <- rownames(seurat_object@meta.data)
    spatial_locs <- as.data.table(GetTissueCoordinates(seurat_object)[,1:2])
    colnames(spatial_locs) <- c("x", "y")
    
    myGiottoObj <- createGiottoObject(raw_exprs = raw_exprs, spatial_locs = spatial_locs)
    myGiottoObj <- normalizeGiotto(gobject = myGiottoObj)
    
    # Create signature matrix for initial PAGE analysis
    all_signatures <- names(signatures_all)
    signature_matrix_complete <- makeSignMatrixPAGE(
      sign_names = all_signatures,
      sign_list = lapply(all_signatures, function(sign) signatures_all[[sign]])
    )
    
    # Run initial PAGE enrichment analysis on all signatures
    myGiottoObj_initial <- runPAGEEnrich(
      gobject = myGiottoObj,
      sign_matrix = signature_matrix_complete,
      min_overlap_genes = 2,
      output_enrichment = c("original", "zscore")
    )
    
    # Run PAGE enrichment analysis with p-values
    myGiottoObj_pval <- runPAGEEnrich(
      gobject = myGiottoObj,
      sign_matrix = signature_matrix_complete,
      min_overlap_genes = 2,
      p_value = TRUE,
      reverse_log_scale = FALSE
    )
    
    # Extract p-value results
    pval_df <- as.data.frame(myGiottoObj_pval@spatial_enrichment$PAGE)
    zscore_df<-as.data.frame(myGiottoObj_initial@spatial_enrichment$PAGE)
    
    zscore_df <- as.data.frame(zscore_df)
    colnames(zscore_df) <- paste0(colnames(zscore_df), "PAGE")
    seurat_object@meta.data <- cbind(seurat_object@meta.data, zscore_df)
    
    # Identify significantly enriched spots (p-value < 0.05)
    significant_spots <- pval_df > 1.301
    
    # Store significant enrichment information in Seurat metadata
    significant_df <- as.data.frame(significant_spots)
    colnames(significant_df) <- paste0(colnames(significant_df), "PAGE_significant")
    seurat_object@meta.data <- cbind(seurat_object@meta.data, significant_df)
    
    return(seurat_object)
    
  }
GSM7058756_C1_fibroenriched<-process_PAGEanalysis1(seurat_object = GSM7058756_C1_fibroenriched,signatures_all = sign)
GSM7058757_C2_fibroenriched<-process_PAGEanalysis1(seurat_object = GSM7058757_C2_fibroenriched,signatures_all = sign)
GSM7058758_C3_fibroenriched<-process_PAGEanalysis1(seurat_object = GSM7058758_C3_fibroenriched,signatures_all = sign)
GSM7058759_C4_fibroenriched<-process_PAGEanalysis1(seurat_object = GSM7058759_C4_fibroenriched,signatures_all = sign)
  
  
V763_1_fibro_md<- V763_1_fibro@meta.data[,139:156]
V371_2_fibro_md <- V371_2_fibro@meta.data[,181:198]
V688_1_fibro_md<- V688_1_fibro@meta.data[,139:156]
V688_2_fibro_md<- V688_2_fibro@meta.data[,139:156]
V015_1_fibro_md<- V015_1_fibro@meta.data[,139:156]
V015_2_fibro_md<- V015_2_fibro@meta.data[,139:156]
V797_1_fibro_md<-V797_1_fibro@meta.data[,139:156]
V797_2_fibro_md<-V797_2_fibro@meta.data[,139:156]
V838_1_fibro_md<-V838_1_fibro@meta.data[,139:156]
V838_2_fibro_md<-V838_2_fibro@meta.data[,139:156]

GSM7058756_C1_fibroenriched_md<-GSM7058756_C1_fibroenriched@meta.data[,77:94]
GSM7058757_C2_fibroenriched_md<-GSM7058757_C2_fibroenriched@meta.data[,77:94]
GSM7058759_C4_fibroenriched_md<-GSM7058759_C4_fibroenriched@meta.data[,77:94]



# Define the Fisher's z-transformation and its inverse
fisher_z <- function(r) {
  return(0.5 * log((1 + r) / (1 - r)))
}

inverse_fisher_z <- function(z) {
  return((exp(2 * z) - 1) / (exp(2 * z) + 1))
}

# Example data: list of correlation matrices and corresponding number of spots
correlation_matrices <- list(cor(V763_1_fibro_md),
                             cor(V371_2_fibro_md),
                             cor(V688_1_fibro_md), 
                             cor(V015_1_fibro_md), 
                             cor(V015_2_fibro_md), 
                             cor(V797_1_fibro_md),
                             cor(V797_2_fibro_md),
                             cor(V838_1_fibro_md),
                             cor(V838_2_fibro_md),
                             cor(GSM7058756_C1_fibroenriched_md),
                             cor(GSM7058757_C2_fibroenriched_md),
                             cor(GSM7058759_C4_fibroenriched_md)) # Add your correlation matrices




num_spots <- c(dim(V763_1_fibro_md)[1],
               dim(V371_2_fibro_md)[1],
               dim(V688_1_fibro_md)[1], 
               dim(V015_1_fibro_md)[1], 
               dim(V015_2_fibro_md)[1], 
               dim(V797_1_fibro_md)[1],
               dim(V797_2_fibro_md)[1], 
               dim(V838_1_fibro_md)[1],
               dim(V838_2_fibro_md)[1],
               dim(GSM7058756_C1_fibroenriched_md)[1],
               dim(GSM7058757_C2_fibroenriched_md)[1],
               dim(GSM7058759_C4_fibroenriched_md)[1]
               ) # Number of spots for each experiment

# Fisher's z-transformation to each matrix
z_matrices <- lapply(correlation_matrices, function(mat) {
  apply(mat, c(1, 2), fisher_z)
})

# Step 2: Weight the transformed values by the number of spots
weighted_z_matrices <- mapply(function(z_matrix, n_spots) {
  z_matrix * n_spots
}, z_matrices, num_spots, SIMPLIFY = FALSE)

# Step 3: Compute the sum of weighted z-values and total number of spots
sum_weighted_z_matrix <- Reduce("+", weighted_z_matrices)
total_spots <- sum(num_spots)

# Compute the weighted average of z-values
mean_weighted_z_matrix <- sum_weighted_z_matrix / total_spots

# Step 4: Apply the inverse Fisher's z-transformation to the averaged matrix
merged_correlation_matrix <- apply(mean_weighted_z_matrix, c(1, 2), inverse_fisher_z)

# Ensure the diagonal remains 1 (since it's a correlation matrix)
diag(merged_correlation_matrix) <- 1

# Print or inspect the final merged correlation matrix
print(merged_correlation_matrix)


# Perform hierarchical clustering
row_model <- hclust(dist(merged_correlation_matrix))
col_model <- hclust(t(dist(merged_correlation_matrix)))


library(ComplexHeatmap)
# Create a Heatmap object with row and column dendrograms
ht <- Heatmap(merged_correlation_matrix, 
              name = "Merged\nCorrelation\ncorrected for n.spots",
              col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
              cluster_rows = row_model,
              cluster_columns = col_model,
              show_row_dend = TRUE,
              show_column_dend = TRUE,
              show_row_names = TRUE,
              show_column_names = TRUE,
              row_dend_side = "left",
              column_dend_side = "top",
              row_names_side = "left",
              column_names_side = "bottom",
              heatmap_legend_param = list(title = "Merged\nCorrelation\ncorrected\nfor n.spots"),
              row_dend_width = unit(4, "cm"),
              column_dend_height = unit(4, "cm"))


# Save the heatmap to a PNG file
Cairo::CairoPNG("FinalMergedCorr_matrix_85_fixed20_final.png", width = 1000, height = 1000)
draw(ht, heatmap_legend_side = "right")
dev.off()


V763_1_fibro_pval<- V763_1_fibro@meta.data[,157:173]
V371_2_fibro_pval <- V371_2_fibro@meta.data[,199:215]
V688_1_fibro_pval<- V688_1_fibro@meta.data[,157:173]
V688_2_fibro_pval<- V688_2_fibro@meta.data[,157:173]
V015_1_fibro_pval<- V015_1_fibro@meta.data[,157:173]
V015_2_fibro_pval<- V015_2_fibro@meta.data[,157:173]
V797_1_fibro_pval<-V797_1_fibro@meta.data[,157:173]
V797_2_fibro_pval<-V797_2_fibro@meta.data[,157:173]
V838_1_fibro_pval<-V838_1_fibro@meta.data[,157:173]
V838_2_fibro_pval<-V838_2_fibro@meta.data[,157:173]

GSM7058756_C1_fibroenriched_pval<-GSM7058756_C1_fibroenriched@meta.data[,95:111]
GSM7058757_C2_fibroenriched_pval<-GSM7058757_C2_fibroenriched@meta.data[,95:111]
GSM7058759_C4_fibroenriched_pval<-GSM7058759_C4_fibroenriched@meta.data[,95:111]

calculate_hypergeometric_pvals <- function(df) {
  categories <- colnames(df)
  comb_matrix <- combn(categories, 2, simplify = FALSE)
  
  pval_matrix <- matrix(0, nrow = length(categories), ncol = length(categories))
  colnames(pval_matrix) <- categories
  rownames(pval_matrix) <- categories
  
  N <- nrow(df)
  
  for (comb in comb_matrix) {
    cat1 <- comb[1]
    cat2 <- comb[2]
    
    K <- sum(df[[cat1]])
    n <- sum(df[[cat2]])
    k <- sum(df[[cat1]] & df[[cat2]])
    
    pval <- phyper(k - 1, K, N - K, n, lower.tail = FALSE)
    
    pval_matrix[cat1, cat2] <- pval
    pval_matrix[cat2, cat1] <- pval
  }
  
  diag(pval_matrix) <- 1  # Diagonal elements should be 1 (or another suitable value, as self-comparisons are not meaningful)
  
  return(pval_matrix)
}


hyp_V763_1<-calculate_hypergeometric_pvals(V763_1_fibro_pval)
hyp_V371_1<-calculate_hypergeometric_pvals(V371_2_fibro_pval)
hyp_V688_1<-calculate_hypergeometric_pvals(V688_1_fibro_pval)
hyp_V688_2<-calculate_hypergeometric_pvals(V688_2_fibro_pval)
hyp_V015_1<-calculate_hypergeometric_pvals(V015_1_fibro_pval)
hyp_V015_2<-calculate_hypergeometric_pvals(V015_2_fibro_pval)
hyp_V797_1<-calculate_hypergeometric_pvals(V797_1_fibro_pval)
hyp_V797_2<-calculate_hypergeometric_pvals(V797_2_fibro_pval)
hyp_V838_1<-calculate_hypergeometric_pvals(V838_1_fibro_pval)
hyp_V838_2<-calculate_hypergeometric_pvals(V838_2_fibro_pval)
hyp_C1<-calculate_hypergeometric_pvals(GSM7058756_C1_fibroenriched_pval)
hyp_C2<-calculate_hypergeometric_pvals(GSM7058757_C2_fibroenriched_pval)
hyp_C4<-calculate_hypergeometric_pvals(GSM7058759_C4_fibroenriched_pval)


# Function to combine p-values using Fisher's method
combine_pvals <- function(pval_matrices) {
  categories <- rownames(pval_matrices[[1]])
  combined_pval_matrix <- matrix(0, nrow = length(categories), ncol = length(categories))
  colnames(combined_pval_matrix) <- categories
  rownames(combined_pval_matrix) <- categories
  
  for (i in 1:length(categories)) {
    for (j in 1:length(categories)) {
      pvals <- sapply(pval_matrices, function(x) x[i, j])
      if (all(pvals == 1)) {
        combined_pval <- 1
      } else {
        chisq_stat <- -2 * sum(log(pvals))
        combined_pval <- pchisq(chisq_stat, df = 2 * length(pvals), lower.tail = FALSE)
      }
      combined_pval_matrix[i, j] <- combined_pval
    }
  }
  
  return(combined_pval_matrix)
}


# Combine p-values from all datasets
combined_pval_matrix <- combine_pvals(list(hyp_V763_1,
                                           hyp_V371_1,
                                           hyp_V688_1,
                                           hyp_V688_2,
                                           hyp_V015_1,
                                           hyp_V015_2,
                                           hyp_V797_1,
                                           hyp_V797_2,
                                           hyp_V838_1,
                                           hyp_V838_2,
                                           hyp_C1,
                                           hyp_C2,
                                           hyp_C4))
combined_pval_df <- as.data.frame(combined_pval_matrix)
combined_pval_df



# Perform hierarchical clustering
row_model <- hclust(dist(combined_pval_df))
col_model <- hclust(t(dist(combined_pval_df)))

library(ComplexHeatmap)
# Create a Heatmap object with row and column dendrograms
ht <- Heatmap(-log10(combined_pval_df+5e-324), 
              name = "Combined\nenrichment\nscore\nvia Fisher's\nMethod",
              col = colorRamp2(c(0, 20, 400), c("white","#E94B3C",'#6C1413')),
              cluster_rows = row_model,
              cluster_columns = col_model,
              show_row_dend = TRUE,
              show_column_dend = T,
              show_row_names = TRUE,
              show_column_names = TRUE,
              row_dend_side = "left",
              column_dend_side = "top",
              row_names_side = "left",
              column_names_side = "bottom",
              heatmap_legend_param = list(title = "Combined\nenrichment\nscore\nvia Fisher's\nMethod"),
              row_dend_width = unit(3, "cm"),
              column_dend_height = unit(4, "cm"))

# Save the heatmap to a PNG file
Cairo::CairoPNG("Hypergeometri_p_values_F_cmplt_nomrCAF_completelink.png", width = 1000, height = 1000)
draw(ht, heatmap_legend_side = "right")
dev.off()



### let's work on pathology annotations!!!
path_V573_1<-read.csv('/home/carlo/Desktop/transferdata/valdolivas/Pathology_SpotAnnotations/Pathologist_Annotations_SN048_A121573_Rep1.csv')
path_V573_2<-read.csv('/home/carlo/Desktop/transferdata/valdolivas/Pathology_SpotAnnotations/Pathologist_Annotations_SN048_A121573_Rep2.csv')
path_V371_1<-read.csv('/home/carlo/Desktop/transferdata/valdolivas/Pathology_SpotAnnotations/Pathologist_Annotations_SN048_A416371_Rep1.csv')
path_V371_2<-read.csv('/home/carlo/Desktop/transferdata/valdolivas/Pathology_SpotAnnotations/Pathologist_Annotations_SN048_A416371_Rep2.csv')
path_V838_1<-read.csv('/home/carlo/Desktop/transferdata/valdolivas/Pathology_SpotAnnotations/Pathologist_Annotations_SN84_A120838_Rep1.csv')
path_V838_2<-read.csv('/home/carlo/Desktop/transferdata/valdolivas/Pathology_SpotAnnotations/Pathologist_Annotations_SN84_A120838_Rep2.csv')
path_V763_1<-read.csv('/home/carlo/Desktop/transferdata/valdolivas/Pathology_SpotAnnotations/Pathologist_Annotations_SN123_A551763_Rep1.csv')
path_V688_1<-read.csv('/home/carlo/Desktop/transferdata/valdolivas/Pathology_SpotAnnotations/Pathologist_Annotations_SN123_A595688_Rep1.csv')
path_V015_1<-read.csv('/home/carlo/Desktop/transferdata/valdolivas/Pathology_SpotAnnotations/Pathologist_Annotations_SN123_A798015_Rep1.csv')
path_V797_1<-read.csv('/home/carlo/Desktop/transferdata/valdolivas/Pathology_SpotAnnotations/Pathologist_Annotations_SN123_A938797_Rep1_X.csv')
path_V763_2<-read.csv('/home/carlo/Desktop/transferdata/valdolivas/Pathology_SpotAnnotations/Pathologist_Annotations_SN124_A551763_Rep2.csv')
path_V688_2<-read.csv('/home/carlo/Desktop/transferdata/valdolivas/Pathology_SpotAnnotations/Pathologist_Annotations_SN124_A595688_Rep2.csv')
path_V015_2<-read.csv('/home/carlo/Desktop/transferdata/valdolivas/Pathology_SpotAnnotations/Pathologist_Annotations_SN124_A798015_Rep2.csv')
path_V797_2<-read.csv('/home/carlo/Desktop/transferdata/valdolivas/Pathology_SpotAnnotations/Pathologist_Annotations_SN124_A938797_Rep2.csv')

colnames(path_V573_1)=c("Barcode","Pathologist.Annotations")
colnames(path_V573_2)=c("Barcode","Pathologist.Annotations")
colnames(path_V371_1)=c("Barcode","Pathologist.Annotations")
colnames(path_V371_2)=c("Barcode","Pathologist.Annotations")
colnames(path_V838_1)=c("Barcode","Pathologist.Annotations")
colnames(path_V838_2)=c("Barcode","Pathologist.Annotations")
colnames(path_V763_1)=c("Barcode","Pathologist.Annotations")
colnames(path_V688_1)=c("Barcode","Pathologist.Annotations")
colnames(path_V015_1)=c("Barcode","Pathologist.Annotations")
colnames(path_V797_1)=c("Barcode","Pathologist.Annotations")
colnames(path_V763_2)=c("Barcode","Pathologist.Annotations")
colnames(path_V688_2)=c("Barcode","Pathologist.Annotations")
colnames(path_V015_2)=c("Barcode","Pathologist.Annotations")
colnames(path_V797_2)=c("Barcode","Pathologist.Annotations")




rownames(path_V573_1)=path_V573_1$Barcode
rownames(path_V573_2)=path_V573_2$Barcode
rownames(path_V371_1)=path_V371_1$Barcode
rownames(path_V371_2)=path_V371_2$Barcode
rownames(path_V838_1)=path_V838_1$Barcode
rownames(path_V838_2)=path_V838_2$Barcode
rownames(path_V763_1)=path_V763_1$Barcode
rownames(path_V688_1)=path_V688_1$Barcode
rownames(path_V015_1)=path_V015_1$Barcode
rownames(path_V797_1)=path_V797_1$Barcode
rownames(path_V763_2)=path_V763_2$Barcode
rownames(path_V688_2)=path_V688_2$Barcode
rownames(path_V015_2)=path_V015_2$Barcode
rownames(path_V797_2)=path_V797_2$Barcode


V573_1<-AddMetaData(V573_1,metadata = path_V573_1)
V573_2<-AddMetaData(V573_2,metadata = path_V573_2)
V371_1<-AddMetaData(V371_1,metadata = path_V371_1)
V371_2<-AddMetaData(V371_2,metadata = path_V371_2)
V838_1<-AddMetaData(V838_1,metadata = path_V838_1)
V838_2<-AddMetaData(V838_2,metadata = path_V838_2)
V763_1<-AddMetaData(V763_1,metadata = path_V763_1)
V688_1<-AddMetaData(V688_1,metadata = path_V688_1)
V015_1<-AddMetaData(V015_1,metadata = path_V015_1)
V797_1<-AddMetaData(V797_1,metadata = path_V797_1)
V763_2<-AddMetaData(V763_2,metadata = path_V763_2)
V688_2<-AddMetaData(V688_2,metadata = path_V688_2)
V015_2<-AddMetaData(V015_2,metadata = path_V015_2)
V797_2<-AddMetaData(V797_2,metadata = path_V797_2)


pdf('PathologistAnnotations_valdolivas.pdf',width = 30,height = 100)
ggarrange(nrow = 14,ncol = 2,
SpatialDimPlot(V573_1,group.by = 'Pathologist.Annotations'),
VlnPlot(V573_1,features = c('mrCAF.score_GM','EMR.score_GM'),group.by = 'Pathologist.Annotations'),
SpatialDimPlot(V573_2,group.by = 'Pathologist.Annotations'),
VlnPlot(V573_2,features = c('mrCAF.score_GM','EMR.score_GM'),group.by = 'Pathologist.Annotations'),
SpatialDimPlot(V371_1,group.by = 'Pathologist.Annotations'),
VlnPlot(V371_1,features = c('mrCAF.score_GM','EMR.score_GM'),group.by = 'Pathologist.Annotations'),
SpatialDimPlot(V371_2,group.by = 'Pathologist.Annotations'),
VlnPlot(V573_2,features = c('mrCAF.score_GM','EMR.score_GM'),group.by = 'Pathologist.Annotations'),
SpatialDimPlot(V838_1,group.by = 'Pathologist.Annotations'),
VlnPlot(V838_1,features = c('mrCAF.score_GM','EMR.score_GM'),group.by = 'Pathologist.Annotations'),
SpatialDimPlot(V838_2,group.by = 'Pathologist.Annotations'),
VlnPlot(V838_2,features = c('mrCAF.score_GM','EMR.score_GM'),group.by = 'Pathologist.Annotations'),
SpatialDimPlot(V763_1,group.by = 'Pathologist.Annotations'),
VlnPlot(V763_1,features = c('mrCAF.score_GM','EMR.score_GM'),group.by = 'Pathologist.Annotations'),
SpatialDimPlot(V688_1,group.by = 'Pathologist.Annotations'),
VlnPlot(V688_1,features = c('mrCAF.score_GM','EMR.score_GM'),group.by = 'Pathologist.Annotations'),
SpatialDimPlot(V015_1,group.by = 'Pathologist.Annotations'),
VlnPlot(V015_1,features = c('mrCAF.score_GM','EMR.score_GM'),group.by = 'Pathologist.Annotations'),
SpatialDimPlot(V797_1,group.by = 'Pathologist.Annotations'),
VlnPlot(V797_1,features = c('mrCAF.score_GM','EMR.score_GM'),group.by = 'Pathologist.Annotations'),
SpatialDimPlot(V763_2,group.by = 'Pathologist.Annotations'),
VlnPlot(V763_2,features = c('mrCAF.score_GM','EMR.score_GM'),group.by = 'Pathologist.Annotations'),
SpatialDimPlot(V688_2,group.by = 'Pathologist.Annotations'),
VlnPlot(V688_2,features = c('mrCAF.score_GM','EMR.score_GM'),group.by = 'Pathologist.Annotations'),
SpatialDimPlot(V015_2,group.by = 'Pathologist.Annotations'),
VlnPlot(V015_2,features = c('mrCAF.score_GM','EMR.score_GM'),group.by = 'Pathologist.Annotations'),
SpatialDimPlot(V797_2,group.by = 'Pathologist.Annotations'),
VlnPlot(V797_2,features = c('mrCAF.score_GM','EMR.score_GM'),group.by = 'Pathologist.Annotations'))
dev.off()



pdf('PathologistAnnotations_valdolivas.pdf',width = 30,height = 100)
ggarrange(nrow = 14,ncol = 2,
          SpatialDimPlot(V573_1,group.by = 'Pathologist.Annotations'),
          VlnPlot(V573_1,features = c('mrCAF.score_GM','EMR.score_GM'),group.by = 'Pathologist.Annotations'),
          SpatialDimPlot(V573_2,group.by = 'Pathologist.Annotations'),
          VlnPlot(V573_2,features = c('mrCAF.score_GM','EMR.score_GM'),group.by = 'Pathologist.Annotations'),
          SpatialDimPlot(V371_1,group.by = 'Pathologist.Annotations'),
          VlnPlot(V371_1,features = c('mrCAF.score_GM','EMR.score_GM'),group.by = 'Pathologist.Annotations'),
          SpatialDimPlot(V371_2,group.by = 'Pathologist.Annotations'),
          VlnPlot(V573_2,features = c('mrCAF.score_GM','EMR.score_GM'),group.by = 'Pathologist.Annotations'),
          SpatialDimPlot(V838_1,group.by = 'Pathologist.Annotations'),
          VlnPlot(V838_1,features = c('mrCAF.score_GM','EMR.score_GM'),group.by = 'Pathologist.Annotations'),
          SpatialDimPlot(V838_2,group.by = 'Pathologist.Annotations'),
          VlnPlot(V838_2,features = c('mrCAF.score_GM','EMR.score_GM'),group.by = 'Pathologist.Annotations'),
          SpatialDimPlot(V763_1,group.by = 'Pathologist.Annotations'),
          VlnPlot(V763_1,features = c('mrCAF.score_GM','EMR.score_GM'),group.by = 'Pathologist.Annotations'),
          SpatialDimPlot(V688_1,group.by = 'Pathologist.Annotations'),
          VlnPlot(V688_1,features = c('mrCAF.score_GM','EMR.score_GM'),group.by = 'Pathologist.Annotations'),
          SpatialDimPlot(V015_1,group.by = 'Pathologist.Annotations'),
          VlnPlot(V015_1,features = c('mrCAF.score_GM','EMR.score_GM'),group.by = 'Pathologist.Annotations'),
          SpatialDimPlot(V797_1,group.by = 'Pathologist.Annotations'),
          VlnPlot(V797_1,features = c('mrCAF.score_GM','EMR.score_GM'),group.by = 'Pathologist.Annotations'),
          SpatialDimPlot(V763_2,group.by = 'Pathologist.Annotations'),
          VlnPlot(V763_2,features = c('mrCAF.score_GM','EMR.score_GM'),group.by = 'Pathologist.Annotations'),
          SpatialDimPlot(V688_2,group.by = 'Pathologist.Annotations'),
          VlnPlot(V688_2,features = c('mrCAF.score_GM','EMR.score_GM'),group.by = 'Pathologist.Annotations'),
          SpatialDimPlot(V015_2,group.by = 'Pathologist.Annotations'),
          VlnPlot(V015_2,features = c('mrCAF.score_GM','EMR.score_GM'),group.by = 'Pathologist.Annotations'),
          SpatialDimPlot(V797_2,group.by = 'Pathologist.Annotations'),
          VlnPlot(V797_2,features = c('mrCAF.score_GM','EMR.score_GM'),group.by = 'Pathologist.Annotations'))
dev.off()
















pdf('SpatialLocation_FIbroblastsEnrichedSpots.pdf')


#pipeline for calculating distances in a seurat object 
calculate_distances <- function(seurat_obj) {
  coords <- GetTissueCoordinates(seurat_obj)
  dist_matrix <- as.matrix(dist(coords[,c('x','y')],method = 'euclidean'))
  return(dist_matrix)
}



# Calculate Neighborhood
calculate_neighborhood <- function(seurat_obj, cell_type, radius) {
  dist_matrix <- calculate_distances(seurat_obj)
  spots <- colnames(seurat_obj)
  cell_type_spots <- which(Idents(seurat_obj) == cell_type)
  
  neighborhood_scores <- sapply(1:length(spots), function(i) {
    neighbors <- which(dist_matrix[i, ] <= radius)
    presence <- sum(neighbors %in% cell_type_spots)
    total <- length(neighbors)
    score <- presence / total
    return(score)
  })
  
  names(neighborhood_scores) <- spots
  return(neighborhood_scores)
}

# Calculate Proximity
calculate_proximity <- function(seurat_obj, cell_type) {
  dist_matrix <- calculate_distances(seurat_obj)
  spots <- colnames(seurat_obj)
  cell_type_spots <- which(Idents(seurat_obj) == cell_type)
  
  proximity_scores <- sapply(1:length(spots), function(i) {
    distances <- dist_matrix[i, cell_type_spots]
    proximity <- mean(1 / (1 + distances))
    return(proximity)
  })
  
  names(proximity_scores) <- spots
  return(proximity_scores)
}




GSM7058757_C2$Tissue_Classification=gsub('FALSE','Tumoral',gsub('TRUE','Stromal',GSM7058757_C2$Stromal_C2L_percentage>GSM7058757_C2$Tumoral_C2L_percentage))




md<-as.data.frame(matrix('NA',ncol = 18,nrow = 3305))
md[match(rownames(GSM7058757_C2_fibroenriched@meta.data[,77:94]),rownames(GSM7058757_C2@meta.data)),]<-GSM7058757_C2_fibroenriched@meta.data[,77:94]
rownames(md)=rownames(GSM7058757_C2@meta.data)
colnames(md)=gsub('$','_spatial',colnames(GSM7058757_C2_fibroenriched@meta.data[,77:94]))

GSM7058757_C2<-AddMetaData(GSM7058757_C2,metadata = md)


library(Seurat)
library(ggplot2)
library(grid)
library(gridExtra)


GSM7058757_C2$Tissue_Classification <- factor(GSM7058757_C2$Tissue_Classification, levels = c("Stromal", "Tumoral"))

color_map <- c("Stromal" = "grey", "Tumoral" = "orange")


Idents(GSM7058757_C2)='Tissue_Classification'
p_tumoral <- SpatialDimPlot(GSM7058757_C2, 
                                cols.highlight = 'orange', 
                                image.alpha = 0,alpha = 0.2,crop = T) +
  ggtitle("Tumoral Category") +
  theme(legend.position = "bottom")

# Create spatial feature plots for "mrCAFPAGE_spatial" and "vCAFPAGE_spatial"
p_mrCAFPAGE <- SpatialFeaturePlot(GSM7058757_C2, features = "mrCAFPAGE_spatial",image.alpha = 0,alpha = 0.2,crop = T) +
  ggtitle("mrCAFPAGE_spatial") +
  theme(legend.position = "bottom")

p_vCAFPAGE <- SpatialFeaturePlot(GSM7058757_C2, features = "vCAFPAGE_spatial",image.alpha = 0,alpha = 0.2,crop = T) +
  ggtitle("vCAFPAGE_spatial") +
  theme(legend.position = "bottom")



aligned_plot<-align_plots(p_mrCAFPAGE,p_vCAFPAGE,p_tumoral,align = 'hv',axis = 'tblr')


library(cowplot)
combined_plot <- ggdraw(aligned_plot[[1]])+draw_plot(aligned_plot[[2]])+draw_plot(aligned_plot[[3]])

# Save combined plot as PDF
CairoPDF("combined_spatial_plot.pdf", width = 12, height = 12)  # Adjust dimensions as needed
SpatialFeaturePlot(GSM7058757_C2,features = c('mrCAFPAGE_spatial','vCAFPAGE_spatial'),crop = T)
dev.off()


colnames(GSM7058757_C2@meta.data)=gsub('mrCAFPAGE_spatial','mrCAF_PAGE',colnames(GSM7058757_C2@meta.data))
colnames(GSM7058757_C2@meta.data)=gsub('vCAFPAGE_spatial','vCAF_PAGE',colnames(GSM7058757_C2@meta.data))                                      
                                       
                                       
GSM7058757_C2$vCAFPAGE_spatial='vCAF_PAGE'
mrCAFPAGE_spatial

library(ggpubr)
CairoPDF("combined_spatial_plot.pdf", width = 12, height = 12)  # Adjust dimensions as needed
ggarrange(SpatialFeaturePlot(GSM7058757_C2,features = c('mrCAF_PAGE'),crop = T,image.alpha = 0)+NoAxes(),
          SpatialFeaturePlot(GSM7058757_C2,features = c('vCAF_PAGE'),crop = T,image.alpha = 0),
          SpatialDimPlot(GSM7058757_C2,group.by = 'Tissue_Classification',image.alpha = 0),ncol = 1)
dev.off()









