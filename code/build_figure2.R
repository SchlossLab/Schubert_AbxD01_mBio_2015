################################################################################
#
# build_figure2.R
#
# This script builds Figure 2.
#
#
# Dependencies...
#   * data/process/abxD0.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.subsample.shared
#   * data/process/abxD0.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.cons.taxonomy
#   * data/process/abxD1.counts
#
# Output...
#   * results/figures/figure2.pdf
#
################################################################################

# read in the metadata file
counts_file <- read.table(file="data/process/abxD1.counts", header=T)
titration <- counts_file[counts_file$abx=="vanc" | counts_file$abx=="cef" | counts_file$abx=="strep",]


# read in the shared file and get the relative abundance
shared_file <- read.table(file="data/process/abxD0.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.subsample.shared", header=T, row.names=2)
shared_file <- shared_file[,!(colnames(shared_file) %in% c("label", "numOtus"))]
n_seqs <- apply(shared_file, 1, sum)[1]
rel_abund <- 100*shared_file/n_seqs


# need to figure out which samples made it through the pipeline and look at those
overlap <- rownames(rel_abund)[which(rownames(rel_abund) %in% rownames(titration))]
rel_abund <- rel_abund[overlap,]
titration <- titration[overlap,]


otu_hyp_test <- function(otu, dose){
    kruskal.test(otu, g=factor(dose))$p.value
}

# limit the analysis to those OTUs that have an median relative abundance over
# 1% within each antibiotic dose

cef <- rel_abund[titration$abx == "cef",]
cef_metadata <- titration[titration$abx == "cef",]

cef_med <- aggregate(cef, by=list(cef_metadata$dose), median)[,-1]
cef_abund <- apply(cef_med, 2, max) > 1.0

cef_abund_good <- cef[,cef_abund]
cef_p_value <- rep(NA, ncol(cef_abund_good))
for(i in 1:ncol(cef_abund_good)){
    cef_p_value[i] <- otu_hyp_test(cef_abund_good[,i], cef_metadata$dose)
}
cef_p_value <- p.adjust(cef_p_value, method="BH")
cef_sig_otus <- colnames(cef_abund_good[,cef_p_value<0.05])




vanc <- rel_abund[titration$abx == "vanc",]
vanc_metadata <- titration[titration$abx == "vanc",]

vanc_med <- aggregate(vanc, by=list(vanc_metadata$dose), median)[,-1]
vanc_abund <- apply(vanc_med, 2, max) > 3.0

vanc_abund_good <- vanc[,vanc_abund]
vanc_p_value <- rep(NA, ncol(vanc_abund_good))
for(i in 1:ncol(vanc_abund_good)){
    vanc_p_value[i] <- otu_hyp_test(vanc_abund_good[,i], vanc_metadata$dose)
}
vanc_p_value <- p.adjust(vanc_p_value, method="BH")
vanc_sig_otus <- colnames(vanc_abund_good[,vanc_p_value<0.05])






strep <- rel_abund[titration$abx == "strep",]
strep_metadata <- titration[titration$abx == "strep",]

strep_med <- aggregate(strep, by=list(strep_metadata$dose), median)[,-1]
strep_abund <- apply(strep_med, 2, max) > 3.0

strep_abund_good <- strep[,strep_abund]
strep_p_value <- rep(NA, ncol(strep_abund_good))
for(i in 1:ncol(strep_abund_good)){
    strep_p_value[i] <- otu_hyp_test(strep_abund_good[,i], strep_metadata$dose)
}
strep_p_value <- p.adjust(strep_p_value, method="BH")
strep_sig_otus <- colnames(strep_abund_good[,strep_p_value<0.05])

sig_otus <- sort(unique(c(strep_sig_otus, cef_sig_otus, vanc_sig_otus)))
x_max <- length(sig_otus)*4-0.5

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
taxonomy <- gsub("_.*", "", taxonomy)
taxonomy <- taxonomy[sig_otus]

otu <- gsub("Otu0*", "", names(taxonomy))

label <- paste0(taxonomy, " (OTU ", otu, ")")

cairo_pdf(file="results/figures/figure2.pdf", width=7.5, height=5.75)
    par(cex=1.2)

    layout(matrix(c(1,5,2,6,3,7,4,8), nrow=4, byrow=T), width=c(1,0.25), height=c(1,1,1,1.2))

    par(mar=c(0.5,5,1.5,0.5))
    sig_cef <- cef[,sig_otus]
    cef_med <- aggregate(sig_cef, by=list(cef_metadata$dose), median)[,-1]
    cef_uci <- aggregate(sig_cef, by=list(cef_metadata$dose), function(x){quantile(x, prob=0.75)})[,-1]
    cef_lci <- aggregate(sig_cef, by=list(cef_metadata$dose), function(x){quantile(x, prob=0.25)})[,-1]
    cef_N <- table(cef_metadata$dose)

    z <- barplot(as.matrix(cef_med), beside=T, names.arg=rep("", ncol(cef_med)), ylim=c(0,23), xlim=c(1.5,x_max), axes=F, col=c("black", "gray", "white"))
    arrows(x0=z, y0=as.matrix(cef_med), y1=as.matrix(cef_uci), angle=90, length=0.02)
    arrows(x0=z, y0=as.matrix(cef_med), y1=as.matrix(cef_lci), angle=90, length=0.02)

    text(x=z[2,sig_otus %in% cef_sig_otus], y=-1, labels="*", cex=2, xpd=TRUE)

    abline(v=seq(4.5, length(sig_otus)*4-0.5, 4), col="gray")
    axis(2, las=1)
    box()
    text(x=0, y=26, label="Cefoperazone", adj=c(0,1), cex=1.2, font=2, xpd=TRUE)

    legend(x=x_max*0.78, y=20,
            legend=paste0(levels(factor(cef_metadata$dose)), " mg/mL (N=", cef_N[levels(factor(cef_metadata$dose))], ")"),
            fill=c("black", "gray", "white"), bg="white")

    sig_strep <- strep[,sig_otus]
    strep_med <- aggregate(sig_strep, by=list(strep_metadata$dose), median)[,-1]
    strep_uci <- aggregate(sig_strep, by=list(strep_metadata$dose), function(x){quantile(x, prob=0.75)})[,-1]
    strep_lci <- aggregate(sig_strep, by=list(strep_metadata$dose), function(x){quantile(x, prob=0.25)})[,-1]
    strep_N <- table(strep_metadata$dose)

    z <- barplot(as.matrix(strep_med), beside=T, names.arg=rep("", ncol(strep_med)), ylim=c(0,65), xlim=c(1.5,x_max), axes=F, col=c("black", "gray", "white"))
    arrows(x0=z, y0=as.matrix(strep_med), y1=as.matrix(strep_uci), angle=90, length=0.02)
    arrows(x0=z, y0=as.matrix(strep_med), y1=as.matrix(strep_lci), angle=90, length=0.02)

    text(x=z[2,sig_otus %in% strep_sig_otus], y=-3, labels="*", cex=2, xpd=TRUE)

    abline(v=seq(4.5, length(sig_otus)*4-0.5, 4), col="gray")
    axis(2, las=1, at=seq(0,60,15))
    mtext(side=2, "Relative abundance (%)", line=3)
    box()
    text(x=0, y=73, label="Streptomycin", adj=c(0,1), cex=1.2, font=2, xpd=T)

    legend(x=x_max*0.78, y=57,
            legend=paste0(levels(factor(strep_metadata$dose)), " mg/mL (N=", strep_N[levels(factor(strep_metadata$dose))], ")"),
            fill=c("black", "gray", "white"), bg="white")




    sig_vanc <- vanc[,sig_otus]
    vanc_med <- aggregate(sig_vanc, by=list(vanc_metadata$dose), median)[,-1]
    vanc_uci <- aggregate(sig_vanc, by=list(vanc_metadata$dose), function(x){quantile(x, prob=0.75)})[,-1]
    vanc_lci <- aggregate(sig_vanc, by=list(vanc_metadata$dose), function(x){quantile(x, prob=0.25)})[,-1]
    vanc_N <- table(vanc_metadata$dose)

    z <- barplot(as.matrix(vanc_med), beside=T, names.arg=rep("", ncol(vanc_med)), ylim=c(0,65), xlim=c(1.5,x_max), axes=F, col=c("black", "gray", "white"))
    arrows(x0=z, y0=as.matrix(vanc_med), y1=as.matrix(vanc_uci), angle=90, length=0.02)
    arrows(x0=z, y0=as.matrix(vanc_med), y1=as.matrix(vanc_lci), angle=90, length=0.02)

    text(x=z[2,sig_otus %in% vanc_sig_otus], y=-3, labels="*", cex=2, xpd=TRUE)

    abline(v=seq(4.5, length(sig_otus)*4-0.5, 4), col="gray")
    axis(2, las=1, at=seq(0,60,15))
    box()
    text(x=0, y=73, label="Vancomycin", adj=c(0,1), cex=1.2, font=2, xpd=TRUE)

    legend(x=x_max*0.77, y=57,
            legend=paste0(levels(factor(vanc_metadata$dose)), " mg/mL (N=", vanc_N[levels(factor(vanc_metadata$dose))], ")"),
            fill=c("black", "gray", "white"), bg="white")

    text(x=apply(z, 2, mean)+1, y=par("usr")[1]-8, xpd=NA, label=label, pos=2, srt=70, cex=1)





    plot.new()

    par(mar=c(0.5,0.5,1.5,5))

    cef_cfu_med <- aggregate(cef_metadata$CFU, by=list(cef_metadata$dose), median)[,-1]+0.1
    cef_cfu_uci <- aggregate(cef_metadata$CFU, by=list(cef_metadata$dose), function(x){quantile(x, prob=0.75)})[,-1]+0.1
    cef_cfu_lci <- aggregate(cef_metadata$CFU, by=list(cef_metadata$dose), function(x){quantile(x, prob=0.25)})[,-1]+0.1

    q <- barplot(as.matrix(cef_cfu_med)+1, beside=T, ylim=c(1, 1e9), log="y", axes=F, col=c("black", "gray", "white"))
    arrows(x0=q, y0=as.matrix(cef_cfu_med), y1=as.matrix(cef_cfu_uci), angle=90, length=0.05)
    arrows(x0=q, y0=as.matrix(cef_cfu_med), y1=as.matrix(cef_cfu_lci), angle=90, length=0.05)
    axis(4, las=1, at=c(1, 1e2, 1e4, 1e6, 1e8), label=c(0, expression(10^2), expression(10^4), expression(10^6), expression(10^8)))
    box()
    #pairwise.wilcox.test(cef_metadata$CFU, cef_metadata$dose)
    text(x=as.vector(q), y=rep(4e8,3), labels=c("a", "a", "b"))


    strep_cfu_med <- aggregate(strep_metadata$CFU, by=list(strep_metadata$dose), median)[,-1]+0.1
    strep_cfu_uci <- aggregate(strep_metadata$CFU, by=list(strep_metadata$dose), function(x){quantile(x, prob=0.75)})[,-1]+0.1
    strep_cfu_lci <- aggregate(strep_metadata$CFU, by=list(strep_metadata$dose), function(x){quantile(x, prob=0.25)})[,-1]+0.1

    q <- barplot(as.matrix(strep_cfu_med)+1, beside=T, ylim=c(1, 1e9), log="y",  axes=F, col=c("black", "gray", "white"))
    arrows(x0=q, y0=as.matrix(strep_cfu_med), y1=as.matrix(strep_cfu_uci), angle=90, length=0.05)
    arrows(x0=q, y0=as.matrix(strep_cfu_med), y1=as.matrix(strep_cfu_lci), angle=90, length=0.05)
    axis(4, las=1, at=c(1, 1e2, 1e4, 1e6, 1e8), label=c(0, expression(10^2), expression(10^4), expression(10^6), expression(10^8)))
    box()
    #pairwise.wilcox.test(strep_metadata$CFU, strep_metadata$dose)
    text(x=as.vector(q), y=rep(4e8,3), labels=c("a", "b", "c"))

    mtext(side=4, "C. difficile colonization (CFU/g)", line=3)


    vanc_cfu_med <- aggregate(vanc_metadata$CFU, by=list(vanc_metadata$dose), median)[,-1]+0.1
    vanc_cfu_uci <- aggregate(vanc_metadata$CFU, by=list(vanc_metadata$dose), function(x){quantile(x, prob=0.75)})[,-1]+0.1
    vanc_cfu_lci <- aggregate(vanc_metadata$CFU, by=list(vanc_metadata$dose), function(x){quantile(x, prob=0.25)})[,-1]+0.1

    q <- barplot(as.matrix(vanc_cfu_med)+1, beside=T, ylim=c(1, 1e9), log="y",  axes=F, col=c("black", "gray", "white"))
    arrows(x0=q, y0=as.matrix(vanc_cfu_med), y1=as.matrix(vanc_cfu_uci), angle=90, length=0.05)
    arrows(x0=q, y0=as.matrix(vanc_cfu_med), y1=as.matrix(vanc_cfu_lci), angle=90, length=0.05)
    axis(4, las=1, at=c(1, 1e2, 1e4, 1e6, 1e8), label=c(0, expression(10^2), expression(10^4), expression(10^6), expression(10^8)))
    box()
    #pairwise.wilcox.test(vanc_metadata$CFU, vanc_metadata$dose)
    text(x=q[2,], y=4e8, labels=c("NS"))

    text(x=q+0.3, y=par("usr")[1]-0.5, xpd=NA, label=c("Low", "Medium","High"), pos=2, srt=70, cex=1.2)
    plot.new()

dev.off()
