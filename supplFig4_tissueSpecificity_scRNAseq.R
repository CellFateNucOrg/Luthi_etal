library(monocle3)
library(ggplot2)
library(cowplot)
library(readxl)
library(dplyr)

theme_set(
  theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          title=element_text(size=9),
          axis.title.y=ggtext::element_markdown(size=9),
          axis.title.x=ggtext::element_markdown(size=9)
    )
)


projectDir="."
fountainsDir=paste0(projectDir,"/fountains")
rnaSeqDir=paste0(projectDir,"/RNAseq_DGE")
rnaTxSeqDir=paste0(projectDir,"/RNAseq_DTE")
publicDataDir=paste0(projectDir,"/publicData")
finalFigDir=paste0(projectDir,"/finalFigures")
if(!dir.exists(finalFigDir)){
  dir.create(finalFigDir)
}

source(paste0(projectDir,"/functions.R"))

#fountains<-readRDS(paste0(fountainsDir,"/fountains_base0_uncorrected_20240125.rds"))
rnaSeqFile<-paste0(rnaSeqDir,"/rds/coh1_noOsc/coh1_noOsc_COH1vsTEVonly_DESeq2_fullResults.rds")

fileList<-data.frame(sampleName=c("COH1cs"),
                     filePath=c(rnaSeqFile))

padjVal=0.05

# file cds_baseline_post_sub.rds downloaded from:
#https://zenodo.org/record/7958249/files/cds_baseline_post_sub.rds?download=1
cds<-readRDS(paste0(publicDataDir,"/Ghaddar2023_cds_baseline_post_sub.rds"))

salmon<-readRDS(file=rnaSeqFile)
salmon<-salmon[!is.na(salmon$chr),]

salmon$upVdown<-"ns"
salmon$upVdown[salmon$padj<padjVal & salmon$log2FoldChange>0]<-"Genes up-regulated upon COH-1 cleavage"
salmon$upVdown[salmon$padj<padjVal & salmon$log2FoldChange<0]<-"Genes down-regulated upon COH-1 cleavage"
sig<-salmon[salmon$upVdown!="ns",c("wormbaseID","upVdown")]

colData(cds)

p1<-plot_cells(cds,color_cells_by="cell_type_group", cell_size=1,
               label_cell_groups = T,group_cells_by="partition",
               group_label_size=1.5) +
  ggtitle("Tissue type")


p2<-plot_cells(cds,color_cells_by="assigned_cell_type", cell_size=1,
               label_cell_groups = T,group_cells_by="cluster",
               group_label_size=1.2) +
  ggtitle("Cell type")


p3<-plot_cells(cds,genes = sig, label_cell_groups = F, cell_size=1) +
  scale_color_gradient(low="white",high="black",name="Percent of gene list") +
  guides(colour = guide_colorbar(title.position = "left")) +
  theme(legend.title = element_text(angle=90,hjust=1),
        legend.key.width=unit(0.02,"npc"))



# tissue types - heatmap -----
sig1<-sig
sig1$upVdown<-factor(sig$upVdown)
levels(sig1$upVdown)<-c("down","up")
sig1$upVdown<-relevel(sig1$upVdown,ref="up")
cell_group_df <- tibble::tibble(cell=row.names(colData(cds)),
                                cell_group=colData(cds)$cell_type_group)

agg_mat <- aggregate_gene_expression(cds, sig1, cell_group_df)

# p4<-pheatmap::pheatmap(t(agg_mat), cluster_rows=T, cluster_cols=F,
#                       scale="column", clustering_method="ward.D2",
#                       fontsize=10, angle_col=45)


hm<-ComplexHeatmap::pheatmap(t(agg_mat), legend = F, cluster_rows = T, cluster_cols = F,
                         heatmap_legend_param = list(title = "Aggregated\ngene\nexpression"),
                         angle_col="45")
gb_heatmap = grid::grid.grabExpr(ComplexHeatmap::draw(hm))
p4 <- gtable::gtable_matrix("hm_gtbl", matrix(list(gb_heatmap)), unit(1, "null"), unit(1, "null"))



# Closer look at Neurons -----
## neuron types ----
neurons_cds <- cds[,grepl("Neurons", colData(cds)$cell_type_group, ignore.case=TRUE)]
cell_group_df <- tibble::tibble(cell=row.names(colData(neurons_cds)),
                                cell_group=colData(neurons_cds)$assigned_cell_type)

agg_mat <- aggregate_gene_expression(neurons_cds, sig1, cell_group_df)


# p<-pheatmap::pheatmap(t(agg_mat), cluster_rows=T, cluster_cols=F,
#                       scale="column", clustering_method="ward.D2",
#                       fontsize=9, angle_col=45)

hm1<-ComplexHeatmap::pheatmap(t(agg_mat), legend = TRUE, cluster_rows = F, cluster_cols = F,
                              heatmap_legend_param = list(title = "Aggregated\ngene\nexpression"),
                              row_order = order(t(agg_mat)[,1],decreasing=T),
                              angle_col="45",fontsize_row=4.5)
hm1

gb_heatmap1 = grid::grid.grabExpr(ComplexHeatmap::draw(hm1))
p5 <- gtable::gtable_matrix("hm_gtbl", matrix(list(gb_heatmap1)), unit(1, "null"), unit(1, "null"))


p<-cowplot::plot_grid(plotlist=list(plot_grid(p1,p2,ncol=2,labels=c("c ","d ")),
                      p3,
                      plot_grid(p4,p5,ncol=2,rel_widths=c(1,1.1))),
                      nrow=3, labels=c("","e ","g "), rel_heights=c(1,1,2))

ggsave(filename=paste0(finalFigDir,"/supplFig4_tissueSpecificity_scRNAseq.png"),p,height=29,width=20,units="cm",device="png")


cds<-readRDS(paste0(publicDataDir,"/Ghaddar2023_cds_baseline_post_sub.rds"))

goi<-c("clec-178","skn-1")

cds_subset<-cds[rowData(cds)$gene_short_name %in% goi]

gene_fits <- fit_models(cds_subset, model_formula_str = "~cell_type_group")
fit_coefs <- coefficient_table(gene_fits)
cell_type_term<-fit_coefs %>%
  select(gene_short_name, term, q_value, estimate)

p1<-plot_genes_violin(cds_subset, group_cells_by="cell_type_group", ncol=2) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + ggtitle("Ghaddar et al. (yAd)")


ASI_cds <- cds[,grepl("ASI", colData(cds)$assigned_cell_type, ignore.case=TRUE)]
intestine_cds <- cds[,grepl("intestine", colData(cds)$cell_type_group, ignore.case=TRUE)]
#coel_cds<-cds[,grepl("coelomocytes", colData(cds)$cell_type_group, ignore.case=TRUE)]
cds1<-combine_cds(list(ASI_cds, intestine_cds))#,coel_cds))

cds_subset<-cds1[rowData(cds1)$gene_short_name %in% goi[1]]
gene_fits <- fit_models(cds_subset, model_formula_str = "~assigned_cell_type")
fit_coefs <- coefficient_table(gene_fits)
cell_type_term<-fit_coefs %>%
  select(gene_short_name, term, q_value, estimate)

p2<-plot_genes_violin(cds_subset, group_cells_by="assigned_cell_type", ncol=2) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + ggtitle("Ghaddar et al. (yAd)")

p<-plot_grid(p1,p2, nrow=2)

asi<-read.csv(paste0(publicDataDir,"/ASI.csv"))
intestine_anterior<-read.csv(paste0(publicDataDir,"/Intestine anterior.csv"))
intestine_middle<-read.csv(paste0(publicDataDir,"/Intestine middle.csv"))
intestine_posterior<-read.csv(paste0(publicDataDir,"/Intestine posterior.csv"))

goi<-c("clec-178", "skn-1")
allTissues<-rbind( asi, intestine_anterior, intestine_middle, intestine_posterior)

allTgoi<-allTissues[allTissues$gene_short_name %in% goi[1],]

ggplot(allTgoi,aes(x=cell_type,y=scaled_TPM,fill=gene_short_name)) +
  geom_bar(stat="identity") + facet_wrap(~gene_short_name)

intestine_average<-mean(c(allTgoi$scaled_TPM[grepl("Intestine", allTgoi$cell_type)]))

# crazy fold change is meaningless as they are all 0.
# statistical test above (gene_fits$model_summary) shows they are not significant
