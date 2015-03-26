
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
taxonomy <- gsub("_1", "", taxonomy)
taxonomy <- taxonomy[colnames(rel_abund_good)]


get_CFU_cor <- function(OTU){
    cor.test(top_dose$CFU, rel_abund_good[,OTU], method="spearman")[c("estimate","p.value")]
}

otu_cfu_corrs <- sapply(1:ncol(rel_abund_good), get_CFU_cor)

colnames(otu_cfu_corrs) <- colnames(rel_abund_good)

sig <- p.adjust(otu_cfu_corrs["p.value",], method="BH") < 0.05
otu_cfu_corrs <- otu_cfu_corrs[,sig]

sig_taxonomy <- taxonomy[sig]
sig_corrs <- unlist(otu_cfu_corrs["estimate",])


sig_taxonomy <- gsub("Actinobacteria", "5Actinobacteria", sig_taxonomy)
sig_taxonomy <- gsub("Bacteroidetes", "1Bacteroidetes", sig_taxonomy)
sig_taxonomy <- gsub("Firmicutes", "2Firmicutes", sig_taxonomy)
sig_taxonomy <- gsub("Proteobacteria", "3Proteobacteria", sig_taxonomy)
sig_taxonomy <- gsub("Tenericutes", "6Tenericutes", sig_taxonomy)
sig_taxonomy <- gsub("Verrucomicrobia", "4Verrucomicrobia", sig_taxonomy)
sig_taxonomy <- gsub("Bacteria$", "Bacteria;7", sig_taxonomy)

family_levels <- levels(factor(sig_taxonomy))
family_levels <- gsub(".*;","",family_levels)
family_levels <- gsub("\\d","",family_levels)
family_levels[family_levels==""] <- "Bacteria"

phylum_levels <- levels(factor(sig_taxonomy))
phylum_levels <- gsub("Bacteria;\\d(.).*", "(\\1)", phylum)
phylum_levels <- gsub("Bacteria;\\d", "", phylum_levels)

tax_labels <- paste(family, phylum)

pdf(file="results/figures/figure2.pdf", width=6.5, height=4.0)
    par(mar=c(10,6,0.5,0.5))

    plot(NA, xlim=c(1,length(tax_labels)), ylim=c(-0.8, 0.8), axes=F, xlab="", ylab="")

    abline(v=1:length(tax_labels), col="gray")

    stripchart(sig_corrs~factor(sig_taxonomy), vertical=T, method="jitter",
                ylim=c(-0.8, 0.8), axes=F, ylab="", pch=1, add=T)


    text(1:length(tax_labels)+0.5, par("usr")[3]-0.05, label=(tax_labels), xpd=NA, pos=2, srt=70, cex=1)
    axis(2, at=seq(-0.75,0.75,0.25), label=format(seq(-0.75,0.75,0.25), nsmall=2L), las=1)
    mtext(2, text="Spearman Correlation\nCoefficient", line=3.5)
    box()
dev.off()

write.table(file="data/process/top_dose_corr.tsv", cbind(taxonomy=gsub("\\d", "", sig_taxonomy), sig_corrs), quote=F, row.names=F, sep="\t")
