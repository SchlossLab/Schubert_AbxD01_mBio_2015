################################################################################
#
# extract_abxD01_IDS_xlsx.R
#
# This script will convert Alyx's original abxD01_IDS.xlsx file into a text file
# where the sample names correspond to the fastq sequence names. We are using
# the file that Alyx sent to Pat on 2015-02-03
#
# Dependencies...
# * data/raw/abx_D01IDS.xlsx
#
# Output...
# * data/process/abx_cdiff_metadata.tsv
#
################################################################################

if(!"openxlsx" %in% rownames(installed.packages())){
    install.packages("openxlsx")
}
library("openxlsx")

spread_sheet <- read.xlsx("data/raw/abxD01_IDS.xlsx", sheet=1, startRow=1, colNames=TRUE)
spread_sheet[spread_sheet=="na"] <- NA
colnames(spread_sheet) <- c("sample", "CFU", "group", "mouse", "day", "abx", "dose", "delayed", "cdiff", "preAbx", "recovDays")

#make the number of CFU a numeric value
spread_sheet$CFU <- as.numeric(spread_sheet$CFU)

#make the cage a character variable instead of a numeric variable
spread_sheet$group <- as.character(spread_sheet$group)

#make the mouse a character variable instead of a numeric variable
spread_sheet$mouse <- as.character(spread_sheet$mouse)

#was the C. difficile challenge delayed?
spread_sheet$delayed <- spread_sheet$delayed == "delayed"

#did the animals get C. difficile 630 delta erm?
spread_sheet$cdiff <- spread_sheet$cdiff == "630dE"

#was the sample from before antibiotics were applied?
spread_sheet$preAbx <- spread_sheet$preAbx == "preAbx"

#make the time into the delay a continuous variable
spread_sheet$recovDays[spread_sheet$recovDays == "no"] <-  NA
spread_sheet$recovDays[spread_sheet$recovDays == "recovD1"] <-  1
spread_sheet$recovDays[spread_sheet$recovDays == "recovD2"] <-  2
spread_sheet$recovDays[spread_sheet$recovDays == "recovD3"] <-  3
spread_sheet$recovDays[spread_sheet$recovDays == "recovD4"] <-  4
spread_sheet$recovDays[spread_sheet$recovDays == "recovD5"] <-  5
spread_sheet$recovDays <- as.numeric(spread_sheet$recovDays)


#there are two duplicate sample names...
spread_sheet$sample[duplicated(spread_sheet$sample)]
#[1] "021-3D-1" "117-4-D4"

#sample name suggests it came from Day -1 when it came from day 1
spread_sheet[spread_sheet$sample=="021-3D-1" & spread_sheet$day==1,"sample"] <- "021-3D1"

#sample name suggests it came from Day 4 when it came from day 5
spread_sheet[spread_sheet$sample=="117-4-D4" & spread_sheet$day==5,"sample"] <- "117-4-D5"

spread_sheet$sample[duplicated(spread_sheet$sample)]
#character(0)

#make rownames the same as the sample names
rownames(spread_sheet) <- spread_sheet$sample
spread_sheet <- spread_sheet[,-1]

# Now we need to make sure that our sample names match the names in the list of
# corrected fastq file names

#the files all start with the correct cage number
sum(grepl("^\\d{1,2}-", rownames(spread_sheet)))
#0

#remove the dash between the animal number and the D that appears in some smaple names
rownames(spread_sheet) <- gsub("^(\\d{3}-\\d)-D", "\\1D", rownames(spread_sheet))



#Read in fastq file names so that we can make sure that we have all of the
#samples and that they are named the same in the metadata (here) and the
#sequence data

fastq <- scan("data/raw/clean_fastqs.txt", what="", quiet=TRUE)
stub <- gsub("(.*)_S.*", "\\1", fastq)
unique_stubs <- unique(stub)
mocks <- grepl("mock", unique_stubs)

sample_stub <- unique_stubs[!mocks]
mock_stub <- unique_stubs[mocks]


#all of the sample names in the spreadsheet have a fastq file
rownames(spread_sheet)[!(rownames(spread_sheet) %in% sample_stub)]
#[1] "021-3D1" "103-3D2" "088-2D4" "0071D1"  "0081D1"


sample_stub[!(sample_stub %in% rownames(spread_sheet))]
#[1] "007-1D1" "008-1D1"


# we can take care of these last two...
rownames(spread_sheet)[rownames(spread_sheet) == "0071D1"] <- "007-1D1"
rownames(spread_sheet)[rownames(spread_sheet) == "0081D1"] <- "008-1D1"

sample_stub[!(sample_stub %in% rownames(spread_sheet))]
#character(0)

rownames(spread_sheet)[!(rownames(spread_sheet) %in% sample_stub)]
#[1] "021-3D1" "103-3D2" "088-2D4"
# these final three samples failed to yield any sequence reads when we sequenced them

write.table(spread_sheet, file="data/process/abx_cdiff_metadata.tsv", sep="\t", quote=F)
