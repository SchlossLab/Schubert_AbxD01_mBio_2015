################################################################################
#
# build_figure3_extra.R
#
#
################################################################################


# read in the metadata file
counts_file <- read.table(file="data/process/abxD1.counts", header=T)
titration <- counts_file[counts_file$abx=="vanc" | counts_file$abx=="cef" | counts_file$abx=="strep",]

#titration_med <- aggregate(titration$CFU, by=list(titration$abx, titration$dose), median)
#abx_cfu <- format(titration_med$x, scientific=T, digits=2)
#abx_cfu <- gsub("e\\+0", "E", abx_cfu)
#abx_cfu <- gsub("0.0E0", "0", abx_cfu)
#names(abx_cfu) <- titration_med$Group.1


# read in the shared file and get the relative abundance
shared_file <- read.table(file="data/process/abxD0.trim.contigs.good.unique.good.filter.unique.precluster.pick.v4.wang.pick.pick.tx.5.subsample.shared", header=T, row.names=2)
shared_file <- shared_file[,!(colnames(shared_file) %in% c("label", "numOtus"))]
n_seqs <- apply(shared_file, 1, sum)[1]
rel_abund <- shared_file/n_seqs


# need to figure out which samples made it through the pipeline and look at those
overlap <- rownames(rel_abund)[which(rownames(rel_abund) %in% rownames(titration))]
rel_abund <- rel_abund[overlap,]
rel_abund <- 100 * rel_abund
titration <- titration[overlap,]


# let's get the relative abundances for those phyla that have at least one
# sample where they are more than 10% of the community
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


titration_phylum_barplot <- function(drug, label){
    drug_rel_abund <- rel_abund_good[titration$abx == drug,]
    drug_metadata <- titration[titration$abx == drug,]

    med_ra <- aggregate(drug_rel_abund, by=list(drug_metadata$dose), median)
    rownames(med_ra) <- paste(med_ra[,1], " mg/mL")
    med_ra <- med_ra[,-1]

    u_qtr <- aggregate(drug_rel_abund, by=list(drug_metadata$dose),
                                        function(x){quantile(x, prob=0.75)})
    rownames(u_qtr) <- paste(u_qtr[,1], " mg/mL")
    u_qtr <- u_qtr[,-1]

    l_qtr <- aggregate(drug_rel_abund, by=list(drug_metadata$dose),
                                        function(x){quantile(x, prob=0.25)})
    rownames(l_qtr) <- paste(l_qtr[,1], " mg/mL")
    l_qtr <- l_qtr[,-1]

    z <- barplot(as.matrix(med_ra), beside=T, ylim=c(0,105), axes=F,
            col=c("black", "gray", "white"), names.arg=rep("", ncol(med_ra)))

    arrows(x0=z, x1=z, y0=as.matrix(med_ra), y1=as.matrix(u_qtr),
                                                    angle=90, length=0.05)
    arrows(x0=z, x1=z, y0=as.matrix(med_ra), y1=as.matrix(l_qtr),
                                                    angle=90, length=0.05)

    legend("topright", bty="n", fill=c("black", "gray", "white"), legend=rownames(med_ra), cex=0.8)
    box()
    axis(2, at=seq(0,100,25), label=seq(0,100,25), las=1)
    text(x=par("usr")[1], y=par("usr")[4]*1.03, label=label,adj = c(0,0),
                font=2, xpd=TRUE)
    z
}

pdf(height=5, width=3.75, file="results/figures/figure3_phylum.pdf")

    z <- layout(
        matrix( c(  1,
                    2,
                    3,
                    4), byrow=T, ncol=1), widths=c(1.0), heights=c(rep(1,3), 0.75)
        )

    par(mar=c(0.5, 5, 1.25, 0.5))

    pos <- titration_phylum_barplot("cef", "Cefoperazone")
    pos <- titration_phylum_barplot("strep", "Streptomycin")
    mtext(side=2, line=3,"Relative Abundance (%)")
    pos <- titration_phylum_barplot("vanc", "Vancomycin")

    text(x=pos[2,]-0.3, y=par("usr")[3]-5, labels=taxonomy, srt=70, cex=1, xpd=NA, adj=c(1,1))


dev.off()
