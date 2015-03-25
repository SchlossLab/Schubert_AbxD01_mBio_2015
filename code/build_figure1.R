#data/process/abxD0.trim.contigs.good.unique.good.filter.unique.precluster.pick.v4.wang.pick.pick.tx.5.subsample.shared
#data/process/abxD0.trim.contigs.good.unique.good.filter.unique.precluster.pick.v4.wang.pick.pick.tx.5.cons.taxonomy
#data/process/abxD0.counts


# read in the metadata file
counts_file <- read.table(file="data/process/abxD1.counts", header=T)
top_dose <- counts_file[counts_file$experiment=="top_dose" | counts_file$abx=="control",]
top_dose_med <- aggregate(top_dose$CFU, by=list(top_dose$abx), median)
abx_cfu <- format(top_dose_med$x, scientific=T, digits=2)
abx_cfu <- gsub("e\\+0", "E", abx_cfu)
abx_cfu <- gsub("0.0E0", "0", abx_cfu)
names(abx_cfu) <- top_dose_med$Group.1

# read in the shared file and get the relative abundance
shared_file <- read.table(file="data/process/abxD0.trim.contigs.good.unique.good.filter.unique.precluster.pick.v4.wang.pick.pick.tx.5.subsample.shared", header=T, row.names=2)
shared_file <- shared_file[,!(colnames(shared_file) %in% c("label", "numOtus"))]
n_seqs <- apply(shared_top_dose, 1, sum)[1]
rel_abund <- shared_top_dose/n_seqs


# need to figure out which samples made it through the pipeline and look at those
overlap <- which(rownames(rel_abund) %in% rownames(counts_file))
rel_abund <- rel_abund[overlap,]
rel_abund <- 100 * rel_abund
top_dose <- top_dose[overlap,]

# let's get the relative abundances for those phyla that have at least one
# sample where they are more than 10% o the community
good_phyla <- apply(rel_abund, 2, max) > 10
rel_abund_good <- rel_abund[,good_phyla]


# let's get the taxonomy data
taxonomy_file <- read.table(file="data/process/abxD0.trim.contigs.good.unique.good.filter.unique.precluster.pick.v4.wang.pick.pick.tx.5.cons.taxonomy", header=T, row.names=1)
taxonomy <- taxonomy_file$Taxonomy
names(taxonomy) <- rownames(taxonomy_file)
taxonomy <- gsub("\\(\\d*\\)", "", taxonomy)
taxonomy <- gsub("Bacteria;", "", taxonomy)
taxonomy <- gsub(";.*", "", taxonomy)
taxonomy <- gsub("unclassified", "Bacteria", taxonomy)
taxonomy <- taxonomy[colnames(rel_abund_good)]


#let's rename the OTUs
colnames(rel_abund_good) <- taxonomy

#let's get the medians and IQRs
med_ra <- aggregate(rel_abund_good, by=list(top_dose$abx), median)
drugs <- med_ra[,1]
med_ra <- med_ra[,-1]
rownames(med_ra) <- drugs
med_ra <- t(med_ra)

l_qtr <- aggregate(rel_abund_good, by=list(top_dose$abx), function(x){quantile(x, probs=c(0.25))})
l_qtr <- l_qtr[,-1]
rownames(l_qtr) <- drugs
l_qtr <- t(l_qtr)

u_qtr <- aggregate(rel_abund_good, by=list(top_dose$abx), function(x){quantile(x, probs=c(0.75))})
u_qtr <- u_qtr[,-1]
rownames(u_qtr) <- drugs
u_qtr <- t(u_qtr)



pdf(height=9, width=3.75, file="results/figures/figure1.pdf")

    z <- layout(
        matrix( c(  9,1,
                    9,2,
                    9,3,
                    9,4,
                    9,5,
                    9,6,
                    9,7,
                    9,8,
                    0,10), byrow=T, ncol=2), widths=c(0.2, 1), heights=c(rep(1,8), 1)
        )

    par(mar=c(0.5, 0.5, 1.0, 0.5))

    pos <- barplot(as.vector(med_ra[,"control"]), ylim=c(0,100), axes=F)
    arrows(x0=pos, x1=pos, y0=med_ra[,"control"], y1=u_qtr[,"control"], angle=90, length=0.1)
    arrows(x0=pos, x1=pos, y0=med_ra[,"control"], y1=l_qtr[,"control"], angle=90, length=0.1)
    axis(2, at=seq(0,100,25), label=seq(0,100,25), las=1)
    text(x=4.9, y=70, label=paste("No antibiotics\n", abx_cfu["control"] ,"CFU/g"),adj = c(1,0))

    pos <- barplot(as.vector(med_ra[,"cipro"]), ylim=c(0,100), axes=F)
    arrows(x0=pos, x1=pos, y0=med_ra[,"cipro"], y1=u_qtr[,"cipro"], angle=90, length=0.1)
    arrows(x0=pos, x1=pos, y0=med_ra[,"cipro"], y1=l_qtr[,"cipro"], angle=90, length=0.1)
    axis(2, at=seq(0,100,25), label=seq(0,100,25), las=1)
    text(x=4.9, y=70, label=paste("Ciprofloxacin\n", abx_cfu["cipro"] ,"CFU/g"),adj = c(1,0))

    pos <- barplot(as.vector(med_ra[,"vanc"]), ylim=c(0,100), axes=F)
    arrows(x0=pos, x1=pos, y0=med_ra[,"vanc"], y1=u_qtr[,"vanc"], angle=90, length=0.1)
    arrows(x0=pos, x1=pos, y0=med_ra[,"vanc"], y1=l_qtr[,"vanc"], angle=90, length=0.1)
    axis(2, at=seq(0,100,25), label=seq(0,100,25), las=1)
    text(x=4.9, y=70, label=paste("Vancomycin\n", abx_cfu["vanc"],"CFU/g"),adj = c(1,0))

    pos <- barplot(as.vector(med_ra[,"amp"]), ylim=c(0,100), axes=F)
    arrows(x0=pos, x1=pos, y0=med_ra[,"amp"], y1=u_qtr[,"amp"], angle=90, length=0.1)
    arrows(x0=pos, x1=pos, y0=med_ra[,"amp"], y1=l_qtr[,"amp"], angle=90, length=0.1)
    axis(2, at=seq(0,100,25), label=seq(0,100,25), las=1)
    text(x=4.9, y=70, label=paste("Ampicillin\n", abx_cfu["amp"],"CFU/g"),adj = c(1,0))

    pos <- barplot(as.vector(med_ra[,"clinda"]), ylim=c(0,100), axes=F)
    arrows(x0=pos, x1=pos, y0=med_ra[,"clinda"], y1=u_qtr[,"clinda"], angle=90, length=0.1)
    arrows(x0=pos, x1=pos, y0=med_ra[,"clinda"], y1=l_qtr[,"clinda"], angle=90, length=0.1)
    axis(2, at=seq(0,100,25), label=seq(0,100,25), las=1)
    text(x=4.9, y=70, label=paste("Clindamycin\n", abx_cfu["clinda"],"CFU/g"),adj = c(1,0))

    pos <- barplot(as.vector(med_ra[,"strep"]), ylim=c(0,100), axes=F)
    arrows(x0=pos, x1=pos, y0=med_ra[,"strep"], y1=u_qtr[,"strep"], angle=90, length=0.1)
    arrows(x0=pos, x1=pos, y0=med_ra[,"strep"], y1=l_qtr[,"strep"], angle=90, length=0.1)
    axis(2, at=seq(0,100,25), label=seq(0,100,25), las=1)
    text(x=4.9, y=70, label=paste("Streptomycin\n", abx_cfu["strep"], "CFU/g"),adj = c(1,0))

    pos <- pos <- barplot(as.vector(med_ra[,"cef"]), ylim=c(0,100), axes=F)
    arrows(x0=pos, x1=pos, y0=med_ra[,"cef"], y1=u_qtr[,"cef"], angle=90, length=0.1)
    arrows(x0=pos, x1=pos, y0=med_ra[,"cef"], y1=l_qtr[,"cef"], angle=90, length=0.1)
    axis(2, at=seq(0,100,25), label=seq(0,100,25), las=1)
    text(x=4.9, y=70, label=paste("Cefoperazone\n", abx_cfu["cef"], "CFU/g"),adj = c(1,0))

    pos <- barplot(as.vector(med_ra[,"metro"]), ylim=c(0,100), axes=F)
    arrows(x0=pos, x1=pos, y0=med_ra[,"metro"], y1=u_qtr[,"metro"], angle=90, length=0.1)
    arrows(x0=pos, x1=pos, y0=med_ra[,"metro"], y1=l_qtr[,"metro"], angle=90, length=0.1)
    axis(2, at=seq(0,100,25), label=seq(0,100,25), las=1)
    text(x=4.9, y=70, label=paste("Metronidazol\n", abx_cfu["metro"],"CFU/g"),adj = c(1,0))


    plot.new()
    par(mar=c(0, 0, 0, 0))
    text(0.25,0.5,"Relative Abundance (%)", srt=90, xpd=T)

    plot.new()
    text(x=as.vector(pos)/4.75, y=1, labels=rownames(med_ra), srt=70, cex=1, font=1, pos=2, xpd=TRUE)

dev.off()
