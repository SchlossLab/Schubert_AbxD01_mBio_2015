################################################################################
#
# build_figure1.R
#
# This script builds Figure 1, which is a barchart of the median relative
# abundance of each phylum found in the mice treated with the top dose of
# antibiotics as well as the untreated control mice. Also included is the
# number of CFU per gram of feces.
#
# Dependencies...
#   * data/process/abxD0.trim.contigs.good.unique.good.filter.unique.precluster.pick.v4.wang.pick.pick.tx.5.subsample.shared
#   * data/process/abxD0.trim.contigs.good.unique.good.filter.unique.precluster.pick.v4.wang.pick.pick.tx.5.cons.taxonomy
#   * data/process/abxD1.counts
#
# Output...
#   * results/figures/figure1.pdf
#
################################################################################



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
n_seqs <- apply(shared_file, 1, sum)[1]
rel_abund <- shared_file/n_seqs


# need to figure out which samples made it through the pipeline and look at those
overlap <- rownames(rel_abund)[which(rownames(rel_abund) %in% rownames(top_dose))]
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

fulldose_phylum_barplot <- function(drug, label){
    pos <- barplot(as.vector(med_ra[,drug]), ylim=c(0,105), axes=F)
    arrows(x0=pos, x1=pos, y0=med_ra[,drug], y1=u_qtr[,drug], angle=90, length=0.1)
    arrows(x0=pos, x1=pos, y0=med_ra[,drug], y1=l_qtr[,drug], angle=90, length=0.1)
    axis(2, at=seq(0,100,25), label=seq(0,100,25), las=1)
    box()
    text(x=4.9, y=70, label=paste(label, "\n", abx_cfu[drug] ,"CFU/g"),adj = c(1,0))
    pos
}

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

    pos <- fulldose_phylum_barplot("control", "No antibiotics")
    pos <- fulldose_phylum_barplot("cipro", "Ciprofloxacin")
    pos <- fulldose_phylum_barplot("vanc", "Vancomycin")
    pos <- fulldose_phylum_barplot("amp", "Ampicillin")
    pos <- fulldose_phylum_barplot("clinda", "Clindamycin")
    pos <- fulldose_phylum_barplot("strep", "Streptomycin")
    pos <- fulldose_phylum_barplot("cef", "Cefoperazone")
    pos <- fulldose_phylum_barplot("metro", "Metronidazol")


    plot.new()
    par(mar=c(0, 0, 0, 0))
    text(0.25,0.5,"Relative Abundance (%)", srt=90, xpd=T)

    plot.new()
    text(x=as.vector(pos)/4.75, y=1, labels=rownames(med_ra), srt=70, cex=1, font=1, pos=2, xpd=TRUE)

dev.off()
