################################################################################
#
# make_counts_file.R
#
# This script will build the counts file that will indicate the density of C.
# difficile on the day after challenge
#
# Dependencies...
# * data/process/abxD0.files
# * data/process/abx_cdiff_metadata.tsv
#
# Output...
# * data/process/abxD1.counts
#
################################################################################

#extract the names of the samples that represented the day of challenge
fecal_files <- read.table(file="data/process/abxD0.files", header=F)
fecal_samples <- fecal_files$V1
fecal_samples <- fecal_samples[!grepl("mock", fecal_samples)]
fecal_samples <- unique(fecal_samples)

day0_samples <- fecal_samples[grepl("D0", fecal_samples)]
day1_samples <- gsub("D0", "D1", day0_samples)

all_metadata <- read.table(file="data/process/abx_cdiff_metadata.tsv", header=T)
all_metadata <- all_metadata[,c("CFU", "abx", "dose", "delayed")]
day1_metadata <- all_metadata[day1_samples,]


#which mice received the top dose of antibiotic?
top_dose <- which(day1_metadata$abx == "amp")
top_dose <- c(top_dose, which(day1_metadata$abx == "cipro"))
top_dose <- c(top_dose, which(day1_metadata$abx == "clinda"))
top_dose <- c(top_dose, which(day1_metadata$abx == "metro"))
top_dose <- c(top_dose, which(day1_metadata$abx == "cef" & day1_metadata$dose == "0.5"))
top_dose <- c(top_dose, which(day1_metadata$abx == "strep" & day1_metadata$dose == "5"))
top_dose <- c(top_dose, which(day1_metadata$abx == "vanc" & day1_metadata$dose == "0.625"))


#which mice were in the dilution study?
dilution <- which(day1_metadata$abx == "cef" & day1_metadata$dose != "0.5")
dilution <- c(dilution, which(day1_metadata$abx == "strep" & day1_metadata$dose != "5"))
dilution <- c(dilution, which(day1_metadata$abx == "vanc" & day1_metadata$dose != "0.625"))


#which mice were in the delay study?
delay <- which(day1_metadata$delayed==TRUE)


#remove the delay samples from the top dose samples
top_dose <- top_dose[-(which(top_dose %in% delay))]


#create an experiment column
day1_metadata[top_dose, "experiment"] <- "top_dose"
day1_metadata[dilution, "experiment"] <- "dilution"
day1_metadata[delay, "experiment"] <- "delay"


#get samples that received no antibiotics
control_samples <- fecal_samples[grepl("D-", fecal_samples)]
control_metadata <- data.frame(matrix(rep(NA, ncol(day1_metadata) * length(control_samples)), nrow=length(control_samples)))
rownames(control_metadata) <- control_samples
colnames(control_metadata) <- colnames(day1_metadata)
control_metadata$CFU <- rep(0, length(control_samples))
control_metadata$abx <- rep("control", length(control_samples))


#make new metadata file
new_metadata <- rbind(day1_metadata, control_metadata)
new_metadata <- new_metadata[,!(colnames(new_metadata) %in% "delayed")]
rownames(new_metadata) <- gsub("D1", "D0", rownames(new_metadata))

#output metadata file
write.table(new_metadata, file="data/process/abxD1.counts", quote=FALSE, sep="\t")
