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
# * data/process/abxD01_IDS.txt
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

spread_sheet[spread_sheet$sample=="021-3D-1" & spread_sheet$day==1,"sample"] <- "021-3D1"
spread_sheet[spread_sheet$sample=="117-4-D4" & spread_sheet$day==5,"sample"] <- "117-4-D5"
