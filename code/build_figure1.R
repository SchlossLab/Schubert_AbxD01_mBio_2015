################################################################################
#
# build_figure1.R
#
# This script builds Figure 1, which is a barchart of the median relative
# abundance of each genus found in the mice treated with the top dose of
# antibiotics as well as the untreated control mice. Also included is the
# number of CFU per gram of feces.
#
# Dependencies...
#   * data/process/abxD0.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.subsample.shared
#   * data/process/abxD0.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.cons.taxonomy
#   * data/process/abxD1.counts
#
# Output...
#   * results/figures/figure1.pdf
#
################################################################################




# let's get the taxonomy data so that we have the string from the kingdom to
# the genus level name or whatever the next level up is that provided a robust
# classification.
taxonomy_file <- read.table(file="data/process/abxD0.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.cons.taxonomy", header=T, row.names=1)
taxonomy <- taxonomy_file$Taxonomy
names(taxonomy) <- rownames(taxonomy_file)
taxonomy <- gsub("\\(\\d*\\)", "", taxonomy)
taxonomy <- gsub(";unclassified", "", taxonomy)
taxonomy <- gsub("/.*", "", taxonomy)
taxonomy <- gsub(";$", "", taxonomy)
taxonomy <- gsub(".*;", "", taxonomy)
#taxonomy <- taxonomy[sig_otus]



# read in the metadata file
counts_file <- read.table(file="data/process/abxD1.counts", header=T)
top_dose <- counts_file[counts_file$experiment=="top_dose" | counts_file$abx=="control",]


# read in the shared file and get the relative abundance
shared_file <- read.table(file="data/process/abxD0.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.subsample.shared", header=T, row.names=2)
shared_file <- shared_file[,!(colnames(shared_file) %in% c("label", "numOtus"))]
n_seqs <- apply(shared_file, 1, sum)[1]
rel_abund <- 100*shared_file/n_seqs


# need to figure out which samples made it through the pipeline and look at those
overlap <- rownames(rel_abund)[which(rownames(rel_abund) %in% rownames(top_dose))]
rel_abund <- rel_abund[overlap,]
top_dose <- top_dose[overlap,]


# limit the analysis to those OTUs that have an median relative abundance over
# 1% within each antibiotic dose
control_rabund <- rel_abund[top_dose$abx == "control",]
control_metadata <- top_dose[top_dose$abx == "control",]
control_median <- apply(control_rabund, 2, median)
control_otus <- control_median > 1.0

otu_hyp_test <- function(drug){

    drug_rabund <- rel_abund[top_dose$abx == drug,]
    drug_metadata <- top_dose[top_dose$abx == drug,]
    drug_otus <- apply(drug_rabund, 2, median) > 3.0

    combined_otus <- control_otus | drug_otus
    combined_rabund <- rbind(drug_rabund, control_rabund)[combined_otus]

    drugged <- c(rep(TRUE, nrow(drug_rabund)), rep(FALSE, nrow(control_rabund)))
    n_otus <- ncol(combined_rabund)

    p_values <- rep(NA, n_otus)

    for(i in 1:n_otus){
        p_values[i] <- wilcox.test(combined_rabund[,i], g=drugged)$p.value
    }

    adj_p_values <- p.adjust(p_values, method="BH")
    colnames(combined_rabund)[adj_p_values<0.05]

}

amp_sig_otus <- otu_hyp_test("amp")
cef_sig_otus <- otu_hyp_test("cef")
cipro_sig_otus <- otu_hyp_test("cipro")
clinda_sig_otus <- otu_hyp_test("clinda")
metro_sig_otus <- otu_hyp_test("metro")
strep_sig_otus <- otu_hyp_test("strep")
vanc_sig_otus <- otu_hyp_test("vanc")

sig_otus <- sort(unique(c(amp_sig_otus, cef_sig_otus, cipro_sig_otus,
                            clinda_sig_otus, metro_sig_otus, strep_sig_otus,
                            vanc_sig_otus)))
x_max <- length(sig_otus) * 1.2

o <- order(control_median[sig_otus], decreasing=T)
rel_abund_sig <- rel_abund[,sig_otus[o]]


#needed?
#top_dose_corr <- read.table(file="data/process/top_dose_corr.tsv", header=T)
#corr <- round(top_dose_corr[sig_otus,"sig_corrs"], digits=2)
#corr[is.na(corr)] <- "NS"

otu <- gsub("Otu0*", "OTU ", names(taxonomy))
names(otu) <- names(taxonomy)

tax_label <- paste0(taxonomy[sig_otus[o]], " (", otu[sig_otus[o]], ")")

single_drug_bars <- function(drug, drug_sig_otus, drug_label){

    drug_rabund <- rel_abund_sig[top_dose$abx == drug,]
    drug_metadata <- top_dose[top_dose$abx == drug,]

    drug_med <- apply(drug_rabund, 2, median)
    drug_uci <- apply(drug_rabund, 2, function(x){quantile(x, prob=0.75)})
    drug_lci <- apply(drug_rabund, 2, function(x){quantile(x, prob=0.25)})

    z <- barplot(drug_med, names.arg=rep("", length(drug_med)),
                    ylim=c(0,1+max(drug_uci)), xlim=c(0,x_max), axes=F,
                    col="white")

    arrows(x0=z, y0=drug_med, y1=drug_uci, angle=90, length=0.05)
    arrows(x0=z, y0=drug_med, y1=drug_lci, angle=90, length=0.05)

    text(x=z[sig_otus[o] %in% drug_sig_otus], y=-0.05*max(drug_uci),
                                            labels="*", cex=2, xpd=TRUE)

    axis(2, las=1)
    box()
    text(x=par("usr")[1], y=par("usr")[4]*1.175, label=drug_label,
                                adj=c(0,1), cex=1.2, font=2, xpd=TRUE)

    summary_stats <- format(quantile(drug_metadata$CFU,
                            prob=c(0.25, 0.50, 0.75)), scientific=T, digits=2)

    summary_string <- paste0(summary_stats[2], " (",
                        summary_stats[1], "-", summary_stats[3], ")")
    summary_string <- gsub("e\\+0", "x10^", summary_string)

    text(x=par("usr")[2], y=1.05*par("usr")[4], labels=summary_string,
                                adj=c(1,0), pos=2, cex=0.8, xpd=TRUE)
    z
}



cairo_pdf(file="results/figures/figure1.pdf", width=4.5, height=10.0)
    par(cex=1.2)

    layout(matrix(c(
            1,
            2,
            3,
            4,
            5,
            6,
            7,
            8,
            9
        ),
        nrow=9, byrow=T), width=c(1), height=c(1,1,1,1,1,1,1,1,1.5))

    par(mar=c(0.5,5,1.5,0.5))

    z <- single_drug_bars("control", "", "No antibiotics")
    z <- single_drug_bars("amp", amp_sig_otus, "Ampicillin")
    z <- single_drug_bars("cef", cef_sig_otus, "Cefoperazone")
    z <- single_drug_bars("cipro", cipro_sig_otus, "Ciprofloxacin")
    z <- single_drug_bars("clinda", clinda_sig_otus, "Clindamycin")
    mtext(side=2, "Relative abundance (%)", line=3, at=110)

    z <- single_drug_bars("metro", metro_sig_otus, "Metronidazole")
    z <- single_drug_bars("strep", strep_sig_otus, "Streptomycin")
    z <- single_drug_bars("vanc", vanc_sig_otus, "Vancomycin")

    text(x=z+0.7, y=par("usr")[1]-10, xpd=NA, label=tax_label, pos=2, srt=70)
    plot.new()

dev.off()
