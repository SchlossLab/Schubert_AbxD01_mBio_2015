################################################################################
#
# make_files_file.R
#
# This script will build the files file that will be used in the make.contigs
# command from within mothur.
#
# Dependencies...
# * data/process/abx_cdiff_metadata.tsv
# * fastq files stored in ../../cdiff_fastqs/
#
# Output...
# * data/process/abxD0.files
#
################################################################################

metadata <- read.table(file="data/process/abx_cdiff_metadata.tsv", header=T)

#did the sample come from a mouse that eventually received Cdiff?
cdiff_challenged <- metadata[metadata$cdiff==TRUE,]

#get the day 0 samples from the abx treated mice
abx_treated <- cdiff_challenged[cdiff_challenged$day == 0, ]

#get the pre-abx samples from mice that would receive abx and treat as the
#control group

control_groups <- c(21, 24, 50, 600, 7, 39, 40, 103, 47, 627, 2, 4, 15, 52, 97, 98)

untreated <-  cdiff_challenged[cdiff_challenged$preAbx &
                                cdiff_challenged$group %in% control_groups, ]

#concatenate the sample names together
samples <- c(rownames(abx_treated), rownames(untreated))

#get the list of fastq file names
fastqs <- list.files("../../cdiff_fastqs/", pattern="*.fastq")


#need a function that will build the files file line for each stub
find_files <- function(stub, files=fastqs){
    pattern <- paste0(stub, ".*R1.*fastq")  #use the stub to create a search pattern
    r1_file <- files[grep(pattern, files)]  #find all of the possible R1 files
    r2_file <- gsub("R1", "R2", r1_file)    #find all of the possible R2 files
    label <- rep(stub, length(r1_file))     #make a vector of the sample name

    #this all is necessary because some samples were sequenced multiple times
    #the output of the above will be a series of vectors that here we will bind
    #as columns, and paste across the columns for each row
    apply(cbind(label, r1_file, r2_file), 1, paste, collapse="\t")
}

files_lines <- lapply(samples, find_files)
write(unlist(files_lines), "data/process/abxD0.files")


r1_mock <- fastqs[grep("mock[^17].*R1.*fastq", fastqs)]
r2_mock <- fastqs[grep("mock[^17].*R2.*fastq", fastqs)]
stub_mock <- gsub("(.*)_S.*", "\\1", r1_mock)

lines_mock <- apply(cbind(stub_mock, r1_mock, r2_mock), 1, paste, collapse="\t")
write(lines_mock, "data/process/abxD0.files", append=T)
