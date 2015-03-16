################################################################################
#
# prep_for_top_dose_amova.R
#
# This script will generate the design and shared files that are needed to
# generate a distance matrix and run amova so that we can see which community
# structures are different between the different antibiotic treatment groups
#
# Dependencies...
# * data/process/abx_cdiff_metadata.tsv
# * data/process/abxD0.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.shared
#
# Output...
# * data/process/abx_topdose.design
# * data/process/abxD0.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.topdose.shared
#
################################################################################


#get the metadata file that was generated earlier...
metadata_file <- "data/process/abx_cdiff_metadata.tsv"
metadata_data <- read.table(file=metadata_file, header=T)


#get the shared file and convert it to a table of counts without the other crap
shared_file <- "data/process/abxD0.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.shared"
shared_data <- read.table(file=shared_file, header=T, row.names=2)
shared_data <- shared_data[,-c(1,2)]

#anticipate which samle will be excluded because there aren't enough reads
enough_seqs <- apply(shared_data, 1, sum) >= 1625
shared_data_enough <- shared_data[enough_seqs,]

#get the sample names so we can pull those samples out of the metadata file
day_zero_metadata <- metadata_data[rownames(shared_data_enough),]

#exclude the stuff from the delayed experiments
no_delay_metadata <- day_zero_metadata[!day_zero_metadata$delayed,]

#exclude the stuff from the titration experiments
top_dose <- which(no_delay_metadata$abx == "amp")
top_dose <- c(top_dose, which(no_delay_metadata$abx == "cipro"))
top_dose <- c(top_dose, which(no_delay_metadata$abx == "clinda"))
top_dose <- c(top_dose, which(no_delay_metadata$abx == "metro"))
top_dose <- c(top_dose, which(no_delay_metadata$abx == "cef" & no_delay_metadata$dose == "0.5"))
top_dose <- c(top_dose, which(no_delay_metadata$abx == "strep" & no_delay_metadata$dose == "5"))
top_dose <- c(top_dose, which(no_delay_metadata$abx == "vanc" & no_delay_metadata$dose == "0.625"))
top_dose_metadata <- no_delay_metadata[top_dose,]

#get those sample names
top_dose_samples <- rownames(top_dose_metadata)


#output a design file where the pre-antibiotic samples are the controls
top_dose_design_file <- "data/process/abx_topdose.design"
treatment <- top_dose_metadata$abx
treatment[top_dose_metadata$day < 0] <- "no_abx"
write.table(cbind(top_dose_samples, treatment), file=top_dose_design_file, quote=F, col.names=F, row.names=F, sep="\t")


#output a new shared file for the correct samples
top_dose_shared_file <- gsub("shared", "topdose.shared", shared_file)
label <- rep(0.03, length(top_dose_samples))
group <- top_dose_samples
n_otus <- rep(ncol(shared_data_enough), length(top_dose_samples))
top_dose_shared <- shared_data_enough[top_dose_samples,]
write.table(cbind(label=label, Group=group, numOtus=n_otus, top_dose_shared), file=top_dose_shared_file, row.names=F, quote=F, sep="\t")
