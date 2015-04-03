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

otus <- gsub("tu0*", "TU ", names(taxonomy))
tax_otu_labels <- paste0(taxonomy, " (", otus, ")")
names(tax_otu_labels) <- names(taxonomy)

# limit the analysis to those OTUs that have an median relative abundance over
# 1% within each experimental unit
grouping <- paste(counts_file$abx, counts_file$dose, counts_file$experiment, sep="-")
medians <- aggregate(rel_abund, by=list(grouping), median)[,-1]
abund <- apply(medians, 2, max) > 1.0
abund_good <- rel_abund[,abund]

# see that it actually stabilized by about 2000 trees, but 1e4 sounds more
# impressive :)
n_trees <- 1e4

# build the random forest model using log transformed CFU counts.
set.seed("6201976")
rf_full <- randomForest(log10(counts_file$CFU+1) ~ ., data=abund_good, importance=TRUE, ntree=n_trees)

# sort the features by their effect on % increase in the standard error
importance_sorted <- sort(importance(rf_full)[,1], decreasing=T)

# fit the model back on the datar and tack it to the end of the counts file
counts_file$fit_full <- predict(rf_full, abund_good)

# there appears to be a natural break after OTU 11 (the 12th feature)
n_features <- 12
importance_subset <- importance_sorted[1:n_features]

# let's get Rsq for the subset
set.seed("6201976")
rf_partial <- randomForest(log10(counts_file$CFU+1) ~ ., data=abund_good[,names(importance_sorted)[1:n_features]], importance=TRUE, ntree=n_trees)


# let's build figure 5
cairo_pdf(file="results/figures/figure5.pdf", width=7.5, height=3.5)
    layout(matrix(c(1,2), nrow=1), width=c(1,1), height=0.8)

    par(mar=c(3,9.5,2,0.5))
    plot(NA, yaxt="n", xlab="", ylab="", xlim=c(0, max(importance_subset)),
        ylim=c(1, length(importance_subset)), axes=F)

    abline(h=1:length(importance_subset), lty=3, col="gray")
    points(x=rev(importance_subset), y=1:length(importance_subset), pch=19, cex=0.8)
    axis(1, at=seq(0,100,25), label=c("0", "", "50", "", "100"), cex=0.8)
    box()
    mtext(side=2, line=9, adj=0, at=1:length(importance_subset), text=rev(tax_otu_labels[names(importance_subset)]), las=2, cex=0.7)
    mtext(side=1, text="% Increase in MSE", line=2.0)
    mtext(text="A", line=8, side=2, at=n_features+1.5, las=2, cex=2, font=2)


    par(mar=c(3,4,2,0.5))
    plot(log10(counts_file$CFU+1), counts_file$fit_full, xlim=c(0,9), ylim=c(0,9),
            xlab="", ylab="", cex=0.8, axes=F)
    axis(1)
    axis(2, las=2)
    box()
    mtext(side=1, text="Observed colonization (log CFU)", line=2.0)
    mtext(side=2, text="Predicted colonization (log CFU)", line=2.5)
    mtext(text="B", line=2.5, side=2, at=10, las=2, cex=2, font=2)
dev.off()
