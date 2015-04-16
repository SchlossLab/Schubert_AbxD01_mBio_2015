################################################################################
#
# build_figure3.R
#
# This script builds Figure 3.
#
#
# Dependencies...
#   * data/process/abxD0.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.subsample.shared
#   * data/process/abxD0.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.cons.taxonomy
#   * data/process/abxD1.counts
#
# Output...
#   * results/figures/figure3.tiff
#
################################################################################

# read in the metadata file
counts_file <- read.table(file="data/process/abxD1.counts", header=T)
delay <- counts_file[counts_file$abx=="amp" | counts_file$abx=="metro",]


# read in the shared file and get the relative abundance
shared_file <- read.table(file="data/process/abxD0.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.subsample.shared", header=T, row.names=2)
shared_file <- shared_file[,!(colnames(shared_file) %in% c("label", "numOtus"))]
n_seqs <- apply(shared_file, 1, sum)[1]
rel_abund <- 100*shared_file/n_seqs


# need to figure out which samples made it through the pipeline and look at those
overlap <- rownames(rel_abund)[which(rownames(rel_abund) %in% rownames(delay))]
rel_abund <- rel_abund[overlap,]
delay <- delay[overlap,]


otu_hyp_test <- function(otu, delay){
    wilcox.test(otu, g=factor(delay))$p.value
}


# limit the analysis to those OTUs that have an median relative abundance over
# 1% within each antibiotic delay group
amp <- rel_abund[delay$abx == "amp",]
amp_metadata <- delay[delay$abx == "amp",]

amp_med <- aggregate(amp, by=list(amp_metadata$experiment), median)[,-1]
amp_abund <- apply(amp_med, 2, max) > 1.0

amp_abund_good <- amp[,amp_abund]
amp_p_value <- rep(NA, ncol(amp_abund_good))
for(i in 1:ncol(amp_abund_good)){
    amp_p_value[i] <- otu_hyp_test(amp_abund_good[,i], amp_metadata$experiment)
}
amp_p_value <- p.adjust(amp_p_value, method="BH")
amp_sig_otus <- colnames(amp_abund_good[,amp_p_value<0.05])




metro <- rel_abund[delay$abx == "metro",]
metro_metadata <- delay[delay$abx == "metro",]

metro_med <- aggregate(metro, by=list(metro_metadata$experiment), median)[,-1]
metro_abund <- apply(metro_med, 2, max) > 3.0

metro_abund_good <- metro[,metro_abund]
metro_p_value <- rep(NA, ncol(metro_abund_good))
for(i in 1:ncol(metro_abund_good)){
    metro_p_value[i] <- otu_hyp_test(metro_abund_good[,i], metro_metadata$experiment)
}
metro_p_value <- p.adjust(metro_p_value, method="BH")
metro_sig_otus <- colnames(metro_abund_good[,metro_p_value<0.05])



sig_otus <- sort(unique(c(amp_sig_otus, metro_sig_otus)))
x_max <- length(sig_otus)*3-0.5



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
taxonomy <- taxonomy[sig_otus]

otu <- gsub("Otu0*", "", names(taxonomy))

label <- paste0("italic(", taxonomy,")~(OTU~", otu, ")")

tiff(file="results/figures/figure3.tiff", width=6.875, height=4.5, unit="in", res=300)
    par(cex=1.2)

    layout_matrix <- matrix(c(1,4,2,5,3,6), nrow=3, byrow=T)
    layout_matrix <- cbind(c(7,7,0), layout_matrix, c(8,8,0))

    layout(layout_matrix, width=c(0.15, 1, 0.25, 0.15), height=c(1,1,1.1))

    par(mar=c(0.75,0.75,1.5,0.5))

    sig_amp <- amp[,sig_otus]
    amp_metadata$experiment <- factor(amp_metadata$experiment, levels=c("top_dose", "delay"))
    amp_med <- aggregate(sig_amp, by=list(amp_metadata$experiment), median)[,-1]
    amp_uci <- aggregate(sig_amp, by=list(amp_metadata$experiment), function(x){quantile(x, prob=0.75)})[,-1]
    amp_lci <- aggregate(sig_amp, by=list(amp_metadata$experiment), function(x){quantile(x, prob=0.25)})[,-1]

    z <- barplot(as.matrix(amp_med), beside=T, names.arg=rep("", ncol(amp_med)), ylim=c(0,60), xlim=c(1.5,x_max), axes=F, col=c("gray", "white"))
    arrows(x0=z, y0=as.matrix(amp_med), y1=as.matrix(amp_uci), angle=90, length=0.02)
    arrows(x0=z, y0=as.matrix(amp_med), y1=as.matrix(amp_lci), angle=90, length=0.02)

    text(x=z[2,sig_otus %in% amp_sig_otus]-0.5, y=-2.5, labels="*", cex=2, xpd=TRUE)

    abline(v=seq(3.5, length(sig_otus)*3-0.5, 3), col="gray")
    axis(2, las=1, at=seq(0,60,15))
    box()
    text(x=0.5, y=67, label="Ampicillin (0.5 mg/mL)", adj=c(0,1), cex=1.2, font=2, xpd=TRUE)

    sig_metro <- metro[,sig_otus]
    metro_metadata$experiment <- factor(metro_metadata$experiment, levels=c("top_dose", "delay"))
    metro_med <- aggregate(sig_metro, by=list(metro_metadata$experiment), median)[,-1]
    metro_uci <- aggregate(sig_metro, by=list(metro_metadata$experiment), function(x){quantile(x, prob=0.75)})[,-1]
    metro_lci <- aggregate(sig_metro, by=list(metro_metadata$experiment), function(x){quantile(x, prob=0.25)})[,-1]

    z <- barplot(as.matrix(metro_med), beside=T, names.arg=rep("", ncol(metro_med)), ylim=c(0,80), xlim=c(1.5,x_max), axes=F, col=c("gray", "white"))
    arrows(x0=z, y0=as.matrix(metro_med), y1=as.matrix(metro_uci), angle=90, length=0.02)
    arrows(x0=z, y0=as.matrix(metro_med), y1=as.matrix(metro_lci), angle=90, length=0.02)

    text(x=z[2,sig_otus %in% metro_sig_otus]-0.5, y=-3, labels="*", cex=2, xpd=TRUE)

    abline(v=seq(3.5, length(sig_otus)*3-0.5, 3), col="gray")
    axis(2, las=1, at=seq(0,80,20))
    box()
    text(x=0.5, y=89, label="Metronidazole (1 mg/mL)", adj=c(0,1), cex=1.2, font=2, xpd=T)


    text(x=apply(z, 2, mean)+1.25, y=par("usr")[1]-8, xpd=NA, label=parse(text=label), pos=2, srt=70, cex=0.9)


    plot.new()




    par(mar=c(0.5,0.5,1.5,1.5))

    n <- table(delay$abx, delay$experiment)

    amp_cfu_med <- aggregate(amp_metadata$CFU, by=list(amp_metadata$experiment), median)[,-1]+0.1
    amp_cfu_uci <- aggregate(amp_metadata$CFU, by=list(amp_metadata$experiment), function(x){quantile(x, prob=0.75)})[,-1]+0.1
    amp_cfu_lci <- aggregate(amp_metadata$CFU, by=list(amp_metadata$experiment), function(x){quantile(x, prob=0.25)})[,-1]+0.1

    q <- barplot(as.matrix(amp_cfu_med)+1, beside=T, ylim=c(1, 1e9), log="y", axes=F, col=c("gray", "white"))
    arrows(x0=q, y0=as.matrix(amp_cfu_med), y1=as.matrix(amp_cfu_uci), angle=90, length=0.05)
    arrows(x0=q, y0=as.matrix(amp_cfu_med), y1=as.matrix(amp_cfu_lci), angle=90, length=0.05)
    axis(4, las=1, at=c(1, 1e2, 1e4, 1e6, 1e8), label=c(0, expression(10^2), expression(10^4), expression(10^6), expression(10^8)))
    box()
    #wilcox.test(amp_metadata$CFU~amp_metadata$experiment, alternative="less")
    text(x=2, y=4e8, labels=c("*"), cex=2)
    text(x=c(1, 1.5,2.5),y=c(4e9, 4e9, 4e9), labels=c("N=", n["amp", 2:1]), xpd=T)

    metro_cfu_med <- aggregate(metro_metadata$CFU, by=list(metro_metadata$experiment), median)[,-1]+0.1
    metro_cfu_uci <- aggregate(metro_metadata$CFU, by=list(metro_metadata$experiment), function(x){quantile(x, prob=0.75)})[,-1]+0.1
    metro_cfu_lci <- aggregate(metro_metadata$CFU, by=list(metro_metadata$experiment), function(x){quantile(x, prob=0.25)})[,-1]+0.1

    q <- barplot(as.matrix(metro_cfu_med)+1, beside=T, ylim=c(1, 1e9), log="y",  axes=F, col=c("gray", "white"))
    arrows(x0=q, y0=as.matrix(metro_cfu_med), y1=as.matrix(metro_cfu_uci), angle=90, length=0.05)
    arrows(x0=q, y0=as.matrix(metro_cfu_med), y1=as.matrix(metro_cfu_lci), angle=90, length=0.05)
    axis(4, las=1, at=c(1, 1e2, 1e4, 1e6, 1e8), label=c(0, expression(10^2), expression(10^4), expression(10^6), expression(10^8)))
    box()
    #wilcox.test(metro_metadata$CFU~metro_metadata$experiment, alternative="less")
    text(x=2, y=rep(4e8), labels=c("*"), cex=2)
    text(x=c(1, 1.5,2.5),y=c(4e9, 4e9, 4e9), labels=c("N=", n["metro", 2:1]), xpd=T)

    text(x=q+0.3, y=par("usr")[1]-0.5, xpd=NA, label=c("1 day recovery", "6 days recovery"), pos=2, srt=70, cex=1.2)
    plot.new()

    plot.new()
    par(mar=c(0,0,0,0))
    text(x=0.5, y=0.5, "Relative abundance (%)", srt=90, cex=1.5)

    plot.new()
    text(x=0.5, y=0.5, expression(italic(C.~difficile)~colonization~(CFU/g)), srt=-90, cex=1.5)


dev.off()
