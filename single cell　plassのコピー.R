getwd()
setwd("./Desktop/single")
library(Seurat)
library(dplyr)
library(patchwork)
#ノーマライズされたデータセットでは、NormalizeDataはスキップする。
#今回のデータセットではノーマライズされている。
Data = read.table(gzfile("dge.txt.gz"),sep="\t") #gz fileの読み込み
Data[1:5,1:5]
Plass <- CreateSeuratObject(counts = Data, project = "Plass", min.cells = 3, min.features = 200)

GetAssayData(object = Plass, slot = "counts")[1:10,1:30]#counts= 生データ
Plass[["RNA"]]@counts[1:10,1:10] #これでもいい
GetAssayData(object = Plass, slot = "data")[1:10,1:30] #objentにした段階でここまで入力されている#data=正規化したデータ
GetAssayData(object = Plass, slot = "scale.data") # 0 xx 0 matrix#scale.data=変動データ
head(Plass$nCount_RNA)
head(Plass$nFeature_RNA)

VlnPlot(Plass,features = c("nCount_RNA","nFeature_RNA"), ncol = 2)
Norm_Data = log(x = Data + 1)
Norm_Data[1:10,1:10]
Norm_Data <- as.matrix(Norm_Data)
Plass <- SetAssayData(object = Plass, slot = "data", new.data = Norm_Data)
Plass[["RNA"]]@data[1:10,1:10]

Plass <- subset(Plass, subset = nFeature_RNA < 2500 )

#FindVariableFeaturesは正規化されたデータに対して分散のでかい遺伝子を見つける。
Plass <- FindVariableFeatures(object = Plass, mean.function = ExpMean , 
                              selection.method = "mvp",
                              dispersion.function = LogVMR, 
                              mean.cutoff = c(0.01, 3),
                              dispersion.cutoff = c(0.4, Inf)) #, nfeatures = 2000)#featureの足切しない？デフォルトだと２０００で切られてる
#メモリを使い切った、と表示された。

top10 <- head(VariableFeatures(Plass),10)
plot1 <- VariableFeaturePlot(Plass)
plot2 <- LabelPoints(plot = plot1, points = top10,repel = TRUE)
CombinePlots(plot = list(plot1,plot2)) #画面が小さくてうめく描画できない
#FindVariableFeaturesによって出てきた遺伝子を元にPCA解析を行う

all.genes <- rownames(Plass)
Plass <- ScaleData(Plass, features = all.genes)
GetAssayData(object = Plass, slot = "scale.data")[1:10,1:10]
#FindVariableFeaturesによって出てきた遺伝子を元にPCA解析を行う
Plass <- RunPCA(Plass, features = VariableFeatures(object = Plass),npcs = 50)
print(Plass[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Plass,dims = 1:3, reduction = "pca") #あんまりきれいに別れてない？？
DimPlot(Plass,  reduction = "pca")
DimHeatmap(Plass, dims=1, cells=500, balanced =TRUE)
DimHeatmap(Plass, dims=1:15, cells=500, balanced =TRUE)
Plass <- JackStraw(Plass, dims = 50, num.replicate = 100)#defolt : dims=20
Plass <- ScoreJackStraw(Plass, dims = 1:50, score.thresh = 1e-05)
JackStrawPlot(Plass,dims = 1:50)
ElbowPlot(Plass,ndims = 50) 

Plass <- FindNeighbors(Plass, dims = 1:50)
#PCの値を元に近隣を分類
Plass <- FindClusters(Plass, resolution = 6)#このresolutionのパラメーターでクラスターの数が変動する
Plass <- RunUMAP(Plass, dims = 1:50)
#上記の計算を元にUMAPを出す
DimPlot(Plass, reduction = "umap")


Plass <- RunTSNE(Plass,dims = 1:50)
DimPlot(Plass, reduction = "tsne")


cluster1.markers <- FindMarkers(Plass, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)
cluster5.markers <- FindMarkers(Plass, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)
#FindMarkers　クラスター
plass.markers <- FindAllMarkers(Plass, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
plass.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

VlnPlot(Plass, features = c("dd-Smed-v6-6208-0", "dd-Smed-v6-11968-0"))

FeaturePlot(Plass, features = c("dd-Smed-v6-6208-0", "dd-Smed-v6-11968-0","dd-Smed-v6-432-0","dd-Smed-v6-19336-0","dd-Smed-v6-74478-0"))

1*2

