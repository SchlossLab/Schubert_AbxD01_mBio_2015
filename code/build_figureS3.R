################################################################################
#
# build_figure4_extra.R
#
#
################################################################################


# read in the metadata file
counts_file <- read.table(file="data/process/abxD1.counts", header=T)
delay <- counts_file[counts_file$abx=="amp" | counts_file$abx=="metro",]

#delay_med <- aggregate(delay$CFU, by=list(delay$abx, delay$dose), median)
#abx_cfu <- format(delay_med$x, scientific=T, digits=2)
#abx_cfu <- gsub("e\\+0", "E", abx_cfu)
#abx_cfu <- gsub("0.0E0", "0", abx_cfu)
#names(abx_cfu) <- delay_med$Group.1


# read in the shared file and get the relative abundance
shared_file <- read.table(file="data/process/abxD0.trim.contigs.good.unique.good.filter.unique.precluster.pick.v4.wang.pick.pick.tx.5.subsample.shared", header=T, row.names=2)
shared_file <- shared_file[,!(colnames(shared_file) %in% c("label", "numOtus"))]
n_seqs <- apply(shared_file, 1, sum)[1]
rel_abund <- shared_file/n_seqs


# need to figure out which samples made it through the pipeline and look at those
overlap <- rownames(rel_abund)[which(rownames(rel_abund) %in% rownames(delay))]
rel_abund <- rel_abund[overlap,]
rel_abund <- 100 * rel_abund
delay <- delay[overlap,]


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


delay_phylum_barplot <- function(drug, label){
    drug_rel_abund <- rel_abund_good[delay$abx == drug,]
    drug_metadata <- delay[delay$abx == drug,]

    drug_metadata$experiment <- factor(drug_metadata$experiment,
                                            levels=c("top_dose", "delay"))

    med_ra <- aggregate(drug_rel_abund, by=list(drug_metadata$experiment), median)
    rownames(med_ra) <- med_ra[,1]
    med_ra <- med_ra[,-1]

    u_qtr <- aggregate(drug_rel_abund, by=list(drug_metadata$experiment),
                                        function(x){quantile(x, prob=0.75)})
    rownames(u_qtr) <- u_qtr[,1]
    u_qtr <- u_qtr[,-1]

    l_qtr <- aggregate(drug_rel_abund, by=list(drug_metadata$experiment),
                                        function(x){quantile(x, prob=0.25)})
    rownames(l_qtr) <- l_qtr[,1]
    l_qtr <- l_qtr[,-1]

    z <- barplot(as.matrix(med_ra), beside=T, ylim=c(0,105), axes=F,
            col=c("gray", "white"), names.arg=rep("", ncol(med_ra)))

    arrows(x0=z, x1=z, y0=as.matrix(med_ra), y1=as.matrix(u_qtr),
                                                    angle=90, length=0.05)
    arrows(x0=z, x1=z, y0=as.matrix(med_ra), y1=as.matrix(l_qtr),
                                                    angle=90, length=0.05)

    if(drug == "amp"){
        legend("topright", bty="n", fill=c("gray", "white"),
                    legend=c("1 day recovery", "6 days recovery"), cex=0.8)
    }

    box()
    axis(2, at=seq(0,100,25), label=seq(0,100,25), las=1)
    text(x=par("usr")[1], y=par("usr")[4]*1.03, label=label,adj = c(0,0),
                font=2, xpd=TRUE)
    z
}

tiff(height=3.25, width=3.75, file="results/figures/figureS3.tiff", res=300, unit="in")

    z <- layout(
        matrix( c(  1,
                    2,
                    3), byrow=T, ncol=1), widths=c(1.0), heights=c(rep(1,2), 0.75)
        )

    par(mar=c(0.5, 5, 1.25, 0.5))

    pos <- delay_phylum_barplot("amp", "Ampicillin")
    pos <- delay_phylum_barplot("metro", "Metronidazole")
    mtext(side=2, at=110, line=3,"Relative Abundance (%)")

    text(x=pos[2,]-0.6, y=par("usr")[3]-5, labels=taxonomy, srt=70, cex=1, xpd=NA, adj=c(1,1), font=3)


dev.off()
