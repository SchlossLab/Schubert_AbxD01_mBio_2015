#data/process/abxD0.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.subsample.shared
#data/process/abxD0.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.cons.taxonomy
#data/process/abxD1.counts

rm(list=ls())

# read in the metadata file
counts_file <- read.table(file="data/process/abxD1.counts", header=T)
top_dose <- counts_file[counts_file$experiment=="top_dose" | counts_file$abx=="control",]

# read in the shared file and get the relative abundance
shared_file <- read.table(file="data/process/abxD0.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.subsample.shared", header=T, row.names=2)
shared_file <- shared_file[,!(colnames(shared_file) %in% c("label", "numOtus"))]
n_seqs <- apply(shared_file, 1, sum)[1]
rel_abund <- shared_file/n_seqs

# need to figure out which samples made it through the pipeline and look at those
overlap <- rownames(rel_abund)[which(rownames(rel_abund) %in% rownames(top_dose))]
rel_abund <- rel_abund[overlap,]
top_dose <- top_dose[overlap,]

abund <- apply(rel_abund, 2, mean) > 0.001
rel_abund_good <- rel_abund[,abund]


# let's get the taxonomy data
taxonomy_file <- read.table(file="data/process/abxD0.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.cons.taxonomy", header=T, row.names=1)
taxonomy <- taxonomy_file$Taxonomy
names(taxonomy) <- rownames(taxonomy_file)
taxonomy <- gsub("\\(\\d*\\)", "", taxonomy)
taxonomy <- gsub(";[^;]*;$", "", taxonomy)
taxonomy <- gsub(";unclassified", "", taxonomy)
taxonomy <- gsub(".*;", "", taxonomy)

taxonomy <- taxonomy[colnames(rel_abund_good)]


get_CFU_cor <- function(OTU){

    cor.test(top_dose$CFU, rel_abund_good[,OTU], method="spearman")[c("estimate","p.value")]

}

otu_cfu_corrs <- matrix(unlist(sapply(1:ncol(rel_abund_good), get_CFU_cor)), ncol=2, byrow=T)
rownames(otu_cfu_corrs) <- colnames(rel_abund_good)
colnames(otu_cfu_corrs) <- c("rho", "pvalue")
sig <- p.adjust(otu_cfu_corrs[,"pvalue"], method="BH") < 0.05

otu_cfu_corrs <- otu_cfu_corrs[sig,]
taxonomy <- taxonomy[sig]


pdf(file="results/figures/figure2.pdf", width=6.5, height=4.0)
par(mar=c(8,6,0.5,0.5))
stripchart(otu_cfu_corrs[,"rho"]~factor(taxonomy), vertical=T, method="jitter",
            ylim=c(-0.8, 0.8), axes=F, ylab="")
axis(1, at=1:16, label=levels(factor(taxonomy)), las=3, tick=F)
axis(2, at=seq(-0.75,0.75,0.25), label=format(seq(-0.75,0.75,0.25), nsmall=2L), las=1)
mtext(2, text="Spearman Correlation\nCoefficient", line=3.5)
box()
dev.off()

#need to sort by Kingdom -> Family label and then relabel with the family or whatever
