################################################################################
#
# extract_abxD01_IDS_xlsx.R
#
# This script will convert Alyx's original abxD01_IDS.xlsx file into a text file
# where the sample names correspond to the fastq sequence names.
#
# Dependencies...
# * data/raw/abx_D01IDS.xlsx
#
# Output...
# * data/process/abxD01_IDS.txt
#
################################################################################

if(!"openxlsx" %in% rownames(installed.packages())){
    install.packages("openxlsx")
}
library("openxlsx")

spread_sheet <- read.xlsx("data/raw/abxD01_IDS.xlsx", sheet=1, startRow=1, colNames=TRUE)
spread_sheet[spread_sheet=="NA"] <- NA
colnames(spread_sheet) <- c("sample", "CFU", "group", "mouse", "day", "abx", "dose", "delayed", "630dE", "preAbx", "recovDays")

#make the number of CFU a numeric value
spread_sheet$CFU <- as.numeric(spread_sheet$CFU)

#make the cage a character variable instead of a numeric variable
spread_sheet$group <- as.character(spread_sheet$group)

#make the mouse a character variable instead of a numeric variable
spread_sheet$mouse <- as.character(spread_sheet$mouse)

#was the C. difficile challenge delayed?
spread_sheet$delayed <- spread_sheet$delayed == "delayed"

#did the animals get C. difficile 630 delta erm?
spread_sheet$630dE <- spread_sheet$630dE == "630dE"

#was the sample from before antibiotics were applied?
spread_sheet$preAbx <- spread_sheet$preAbx == "preAbx"

#make the time into the delay a continuous varaible
spread_sheet$recovDays[spread_sheet$recovDays == "no"] <-  NA
spread_sheet$recovDays[spread_sheet$recovDays == "recovD1"] <-  1
spread_sheet$recovDays[spread_sheet$recovDays == "recovD2"] <-  2
spread_sheet$recovDays[spread_sheet$recovDays == "recovD3"] <-  3
spread_sheet$recovDays[spread_sheet$recovDays == "recovD4"] <-  4
spread_sheet$recovDays[spread_sheet$recovDays == "recovD5"] <-  5
spread_sheet$recovDays <- as.numeric(spread_sheet$recovDays)
