################################################################################
#
# build_figure5.R
#
# This script builds Figure 5.
#
#
# Dependencies...
#   * data/process/abxD0.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.subsample.shared
#   * data/process/abxD0.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.cons.taxonomy
#   * data/process/abxD1.counts
#
# Output...
#   * results/figures/figure5.pdf
#
################################################################################

load_package <- function(package){
    if(!(package %in% rownames(installed.packages()))){
        install.packages(package)
    }
    library(package, quietly=TRUE, character.only=TRUE)
}

load_package("randomForest")



# read in the metadata file
counts_file <- read.table(file="data/process/abxD1.counts", header=T)

# read in the shared file and get the relative abundance
shared_file <- read.table(file="data/process/abxD0.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.subsample.shared", header=T, row.names=2)
shared_file <- shared_file[,!(colnames(shared_file) %in% c("label", "numOtus"))]
n_seqs <- apply(shared_file, 1, sum)[1]
rel_abund <- 100*shared_file/n_seqs


# need to figure out which samples made it through the pipeline and look at those
overlap <- rownames(rel_abund)[which(rownames(rel_abund) %in% rownames(rel_abund))]
rel_abund <- rel_abund[overlap,]
counts_file <- counts_file[overlap,]


# limit the analysis to those OTUs that have an median relative abundance over
# 1% within each experimental unit
grouping <- paste(counts_file$abx, counts_file$dose, counts_file$experiment, sep="-")

medians <- aggregate(rel_abund, by=list(grouping), median)[,-1]
abund <- apply(medians, 2, max) > 1.0
abund_good <- rel_abund[,abund]

n_trees <- 1e4

rf_full <- randomForest(log(counts_file$CFU+1) ~ ., data=abund_good, importance=TRUE, ntree=n_trees)
counts_file$fit_full <- predict(rf_full, abund_good)
plot(log(counts_file$CFU+1), counts_file$fit_full)

imp_full <- sort(importance(rf)[,1], decreasing=T)


# let's get the taxonomy data so that we have the string from the kingdom to
# the family level name or whatever the next level up is that provided a robust
# classification.
taxonomy_file <- read.table(file="data/process/abxD0.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.cons.taxonomy", header=T, row.names=1)
taxonomy <- taxonomy_file$Taxonomy
names(taxonomy) <- rownames(taxonomy_file)
taxonomy <- gsub("\\(\\d*\\)", "", taxonomy)
taxonomy <- gsub(";unclassified", "", taxonomy)
taxonomy <- gsub("_1", "", taxonomy)
taxonomy <- gsub(";$", "", taxonomy)
taxonomy <- gsub("/.*", "", taxonomy)
taxonomy <- gsub(".*;", "", taxonomy)
taxonomy <- gsub("_sensu_stricto", "", taxonomy)

taxonomy <- taxonomy[names(imp_full)]

corr <- round(top_experiment_corr[sig_otus,"sig_corrs"], digits=2)
corr[is.na(corr)] <- "NS"

otu <- gsub("Otu0*", "", names(taxonomy))

label <- paste0(taxonomy, " (OTU ", otu, ")")


cairo_pdf(file="results/figures/figure5.pdf", width=7.5, height=3.5)

layout(matrix(c(1,2), nrow=1), width=c(1,1.2), height=0.8)



mtext(text="A", line=1.5, side=2, at=23, las=2, cex=2, font=2)





plot_importance(rf_baseline_forest)
mtext(text="B", line=8.5, side=2, at=15, las=2, cex=2, font=2)

dev.off()
