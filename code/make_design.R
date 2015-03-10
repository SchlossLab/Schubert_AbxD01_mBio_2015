################################################################################
#
# make_design.R
#
# This script will generate the various design files that are needed to run the
# amova command from within mothur
#
# Dependencies...
# * The shared file that will be used in the amova command
# * The table that maps the sample names to the treatments (data/raw/abxD01_IDS.txt)
# * Information about the type of design file to make.
#
# Output...
# * The relevant design file
#
################################################################################

mapping_file <- "data/raw/abxD01_IDS.txt"
mapping_data <- read.table(file=mapping_file, header=T)
day_zero <- mapping_data[mapping_data$day==0,]

shared_file <- "data/process/abxD0.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.subsample.shared"
shared_data <- read.table(file=shared_file, row.names=2, header=T)
shared_data <- shared_data[,-c(1,2)]
samples <- rownames(shared_data)
