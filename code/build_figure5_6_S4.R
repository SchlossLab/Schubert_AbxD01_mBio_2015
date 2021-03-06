################################################################################
#
# build_figure5_6.R
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
#   * results/figures/figure6.pdf
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
logCFU <- log10(counts_file$CFU+1)

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
names(otus) <- names(taxonomy)

tax_otu_labels <- paste0("italic('", taxonomy, "')~plain('(", otus, ")')")
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
rf_full <- randomForest(logCFU ~ ., data=abund_good, importance=TRUE, ntree=n_trees)

write(paste("full", ncol(abund_good), rf_full$rsq[n_trees], sep="\t"), file="data/process/random_forest.data")

# sort the features by their effect on % increase in the standard error
importance_sorted <- sort(importance(rf_full)[,1], decreasing=T)

# fit the model back on the datar and tack it to the end of the counts file
counts_file$fit_full <- predict(rf_full, abund_good)

# let's plot the top 12 features
n_features <- 12
importance_subset <- importance_sorted[1:n_features]


# let's get Rsq for the subset
set.seed("6201976")
rf_partial <- randomForest(logCFU ~ ., data=abund_good[,names(importance_sorted)[1:n_features]], importance=TRUE, ntree=n_trees)
write(paste("partial", n_features, rf_partial$rsq[n_trees], "\t"), file="data/process/random_forest.data", append=T)



colonized <- rep("yes", length(logCFU))
colonized[logCFU<2] <- "no"
binary_model <- randomForest(x=abund_good, y=as.factor(colonized), importance=TRUE, ntree=n_trees)
write(paste("binary", ncol(abund_good), binary_model$err.rate[n_trees], "\t"), file="data/process/random_forest.data", append=T)

# supplemental figure 4: full feature importance plot
tiff(file="results/figures/figureS4.tiff", width=4, height=5.0, unit="in", res=300)

    par(mar=c(3,8,0.5,0.5))
    plot(NA, yaxt="n", xlab="", ylab="", xlim=c(min(importance_sorted), 100),
        ylim=c(1, length(importance_sorted)), axes=F)

    abline(h=1:length(importance_sorted), lty=3, col="gray")
    points(x=rev(importance_sorted), y=1:length(importance_sorted), pch=19, cex=0.8)
    axis(1, at=seq(0,100,25), label=c("0", "", "50", "", "100"), cex=0.8)
    box()
    mtext(side=2, line=7.5, adj=0, at=1:length(importance_sorted),
            text=parse(text=rev(tax_otu_labels[names(importance_sorted)])),
            las=2, cex=0.6,
            col=c(rep("black", length(importance_sorted)-12), rep("red", 12)))
    mtext(side=1, text="% Increase in MSE", line=2.0)

dev.off()


# figures 5 and 6: color fit by treatment group
clrs <- c(
    "amp-0.5-delay" = "#6600FF", #"red",
    "amp-0.5-top_dose" = "#6600FF", #"red",
    "cef-0.1-dilution" = "#0092FFFF", #"pink",
    "cef-0.3-dilution" = "#0092FFFF", #"pink",
    "cef-0.5-top_dose" = "#0092FFFF", #"pink",
    "cipro-10mg/kg-top_dose" = "#FF0000FF", #"lightgreen",
    "clinda-10mg/kg-top_dose" = "#FFCC00", #"darkgreen",
    "control-NA-NA" = "black",
    "metro-1-delay" = "#FF00DBFF", #"blue",
    "metro-1-top_dose" = "#FF00DBFF",#"blue",
    "strep-0.1-dilution" = "#33FFFF", #"orange",
    "strep-0.5-dilution" = "#33FFFF", #"orange",
    "strep-5-top_dose" = "#33FFFF", #"orange",
    "vanc-0.1-dilution" = "#00CC00", #"purple",
    "vanc-0.3-dilution" = "#00CC00", #"purple",
    "vanc-0.625-top_dose" = "#00CC00" #"purple"
)

pch <- c(
    "amp-0.5-delay" = 1,
    "amp-0.5-top_dose" = 19,
    "cef-0.1-dilution" = 15,
    "cef-0.3-dilution" = 17,
    "cef-0.5-top_dose" = 19,
    "cipro-10mg/kg-top_dose" = 19,
    "clinda-10mg/kg-top_dose" = 19,
    "control-NA-NA" = 19,
    "metro-1-delay" = 1,
    "metro-1-top_dose" = 19,
    "strep-0.1-dilution" = 15,
    "strep-0.5-dilution" = 17,
    "strep-5-top_dose" = 19,
    "vanc-0.1-dilution" = 15,
    "vanc-0.3-dilution" = 17,
    "vanc-0.625-top_dose" = 19
)
make_gray_plot <- function(drug){
    plot_data <- !grepl(drug, grouping)
    plot(logCFU[plot_data], counts_file$fit_full[plot_data], xlim=c(0,9),
            ylim=c(0,9), xlab="", ylab="", cex=0.8, axes=F, col="gray",
            pch=pch[grouping[plot_data]])
    box()
}

make_colored_plot <- function(drug){
    plot_data <- grepl(drug, grouping)
    points(logCFU[plot_data], counts_file$fit_full[plot_data],
            cex=0.8, col=clrs[grepl(drug, names(clrs))],
            pch=pch[grouping[plot_data]])
    axis(1, labels=rep("", 6), at=seq(0,10,2), tick=T)
    axis(2, labels=rep("", 6), at=seq(0,10,2), tick=T)
}



tiff(file="results/figures/figure5.tiff", width=5.0, height=6.0, unit="in", res=300)

    design <- matrix(1:8, nrow=4, byrow=T)
    design <- cbind(c(rep(9,4)), design)
    design <- rbind(design, c(0,10,10))
    layout(design, widths=c(0.2,1,1), heights=c(1,1,1,1,0.3))


    par(mar=c(0.5,1,1.5,0.5))

    #control
    make_gray_plot("control")
    make_colored_plot("control")
    text(x=-0.25, y=9.75, xpd=TRUE, label="No antibiotics", adj=c(0,0), font=2)
    #legend("topleft", pch=19, col="black", bty="n", legend="No antibiotics", cex=0.7)
    axis(2, las=2)

    #cipro
    make_gray_plot("cipro")
    make_colored_plot("cipro")
    text(x=-0.25, y=9.75, xpd=TRUE, label="Ciprofloxacin", adj=c(0,0), font=2)
    #legend("topleft", pch=19, bty="n", col="lightgreen", legend="Ciprofloxacin", cex=0.7)

    #vanc
    make_gray_plot("vanc")
    make_colored_plot("vanc")
    text(x=-0.25, y=9.75, xpd=TRUE, label="Vancomycin", adj=c(0,0), font=2)
    legend("topleft", col=clrs["vanc-0.625-top_dose"], cex=0.7, pch=c(19, 17, 15), bty="n",
        legend=c("0.625 mg/mL", "0.3 mg/mL", "0.1 mg/mL"))
    axis(2, las=2)

    #amp
    make_gray_plot("amp")
    make_colored_plot("amp")
    text(x=-0.25, y=9.75, xpd=TRUE, label="Ampicillin", adj=c(0,0), font=2)
    legend("topleft", col=clrs["amp-0.5-top_dose"], cex=0.7, pch=c(19, 1), bty="n",
        legend=c("1 day recovery", "6 days recovery"))

    #clinda
    make_gray_plot("clinda")
    make_colored_plot("clinda")
    text(x=-0.25, y=9.75, xpd=TRUE, label="Clindamycin", adj=c(0,0), font=2)
    #legend("topleft", pch=19, bty="n", col="darkgreen", legend="Clindamycin", cex=0.7)
    axis(2, las=2)

    #strep
    make_gray_plot("strep")
    make_colored_plot("strep")
    text(x=-0.25, y=9.75, xpd=TRUE, label="Streptomycin", adj=c(0,0), font=2)
    legend("topleft", col=clrs["strep-5-top_dose"], cex=0.7, pch=c(19, 17, 15), bty="n",
                legend=c("5 mg/mL", "0.5 mg/mL", "0.1 mg/mL"))

    #metro
    make_gray_plot("metro")
    make_colored_plot("metro")
    text(x=-0.25, y=9.75, xpd=TRUE, label="Metronidazole", adj=c(0,0), font=2)
    legend("topleft", col=clrs["metro-1-top_dose"], cex=0.7, pch=c(19, 1), bty="n",
        legend=c("1 day recovery", "6 days recovery"))
    axis(1)
    axis(2, las=2)

            #cef
    make_gray_plot("cef")
    make_colored_plot("cef")
    text(x=-0.25, y=9.75, xpd=TRUE, label="Cefoperazone", adj=c(0,0), font=2)
    legend("topleft", col=clrs["cef-0.5-top_dose"], cex=0.7, pch=c(19, 17, 15), bty="n",
        legend=c("0.5 mg/mL", "0.3 mg/mL", "0.1 mg/mL"))
    axis(1)


    plot.new()
    par(mar=rep(0.1,4))
    text(x=0.25,y=0.5, label="Observed colonization (log CFU)", cex=1.2, srt=90)


    plot.new()
    par(mar=rep(0.1,4))
    text(x=0.5,y=0.2, label="Predicted colonization (log CFU)", cex=1.2)
dev.off()



# let's build Figure 6 (w/ color & pch)

correlations <- read.table(file="data/process/correlation_analysis.tsv", header=T, row.names=1)
correlations$sig_corrs <- round(correlations$sig_corrs, 2)

rho <- correlations[names(importance_subset), "sig_corrs"]
rho[is.na(rho)] <- "N.S."

tax_otu_imp_labels <- paste0("bolditalic('", taxonomy[names(importance_subset)],
                        "')~bold('(",
                        otus[names(importance_subset)], "; \u03C1=", rho, ")')")
names(tax_otu_imp_labels) <- names(taxonomy[names(importance_subset)])

tiff(file="results/figures/figure6.tiff", width=6.875, height=7.5, unit="in", res=300)

    #want to jitter the relative abundance for those mice that had no Cdiff
    #colonization
    cd_zeroes <- logCFU == 0
    logCFU[cd_zeroes] <- runif(sum(cd_zeroes),0,1)

    par(mar=c(0.5,0.5,0.5,0.5))

    design <- matrix(1:n_features, nrow=4, byrow=T)
    design <- cbind(c(rep(13,4)), design)
    design <- rbind(design, c(0,14,14,14))
    layout(design, widths=c(0.3,1,1,1), heights=c(1,1,1,1,0.3))

    for(i in 1:n_features){
        #get the row and column number for each spot in the layout
        row <- ceiling(i/3)
        column <- ((i-1) %% 3) + 1

        #extract the relative abundance data for this OTU
        otu_abund <- rel_abund[,names(importance_subset)[i]]

        #want to jitter the number of tumors for those mice that had a zero
        #relative abundance
        ra_zeroes <- otu_abund == 0
        otu_abund[ra_zeroes] <- runif(sum(ra_zeroes),1.0e-2,1.5e-2)

        #plot the relative abundance with the number of cdiff for each animal. plot
        #on consistent log scaled x-axis for all OTUs. will throw errors because it
        #can't plot zeroes on a log scale
        plot(otu_abund[grouping == "control-NA-NA"], logCFU[grouping == "control-NA-NA"],
            log="x", pch=19,
            col="black",
            cex=0.8,
            ylab="", xlab="",
            xlim=c(1e-2, 100), ylim=c(0,9),
            yaxt="n", xaxt="n"
        )

        points(otu_abund[grouping != "control-NA-NA"], logCFU[grouping != "control-NA-NA"],
            pch=pch[grouping[grouping != "control-NA-NA"]],
            col=clrs[grouping[grouping != "control-NA-NA"]],
            cex=0.8
        )

        #create a vertical line to denote the limit of detection
        abline(v=2.2e-2, col="gray")

        #create a horizontal line to denote the limit of detection
        abline(h=1.5, col="gray")

        #put the OTU label in the upper left corner of the plot
        text(x=0.7e-2, y=8.8, label=parse(text=tax_otu_imp_labels[i]), pos=4, cex=0.9)

        #if it's on the bottom row, put a customized axis indicating the % rabund
        if(row == 4){
            axis(1, at=c(1.25e-2, 1e-1,1e0,1e1,1e2),
                    label=c("0", "0.1", "1", "10", "100"),
                    cex.axis=1.5)
        }

        #if it's in the first column turn the axis labels to be horizontal
        if(column == 1){
            axis(2, las=2, cex.axis=1.5)
        }
    }

    plot.new()
    text(x=0.15, y=0.5, label="Observed colonization (log CFU)", cex=1.5, srt=90)

    plot.new()
    text(x=0.5, y=0.2, label="Relative abundance at Day 0 (%)", cex=1.5)

dev.off()
