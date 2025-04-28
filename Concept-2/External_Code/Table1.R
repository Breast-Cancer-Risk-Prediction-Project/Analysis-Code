#########################
# Table 1
#########################
#Change filepath to wherever you saved the filepath
source("filepath/Table1_Fun.R")

#Please read in your CLEAN file 
valdata<-read.csv("filepath/bcrpp_valdata.csv")


table1<-table1fun(valdata)
table1

#PLEASE INSERT FILEPATH AND THE APPROPRIATE ACRONYM 

write.table(table1, "filepath/Table1_CTS.csv", col.names=T, row.names=F, append=T, sep=',')
write.csv(table1, "filepath/Table1_CTS.csv")
