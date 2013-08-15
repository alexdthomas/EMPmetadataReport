#restart the metdata analysis using Git in the new directory
#redo first_to_merge_class_issues, 
#this time set characterAsFactor=FALSE when import

library(BiocInstaller)
library(XML)
library(XMLSchema)
library(SSOAP)
library(rols)
library(reshape2)
library(plyr)


setwd("~/EarthMicrobiomeProject/R/EMPmetadataReport")

def.par <- par(no.readonly = TRUE) # save default, for resetting...

#set directory to location of all mapping files
setwd("~/EarthMicrobiomeProject/QIIME_metadata_download/mapping_files")

#list mapping files
map.file.names<-list.files()
#create list of data frames for each mapping file
all.maps<-list()

#import all mapping files, name as individual objects and place in list
#note characterAsFactors=FALSE to reduce class issues
#though this only removes issues of character -> factor and vice versa issues
#does not resolve charcter-> integer, numeric, date etc...
for (i in 1:length(map.file.names)){
	all.maps[[i]]<-assign(paste("map", i, sep=""), read.delim(map.file.names[i], quote="", stringsAsFactors=FALSE))
}

#name studies in list
names(all.maps)<-unlist(lapply(all.maps, function(x) unique(x[,"TITLE"])))

#reset working directory 
setwd("~/EarthMicrobiomeProject/R/EMPmetadataReport")

#check merge
emp.map<-Reduce(function(x, y) merge(x, y, all=TRUE), all.maps)
#well no more warnings, so character/ factor resolved

#search colnames
colnames(emp.map)[grep('contact', colnames(emp.map), ignore.case=TRUE)]

#explore data
#make a list of how many values are NA out of total
unlist(lapply(all.maps, function(x) paste(length(which(is.na(x))), 
																					nrow(x)*ncol(x), sep="/")))

#mean NA?
sort(unlist(lapply(all.maps, function(x) mean(is.na(x)))))
#one of them is 25% NA?
which(unlist(lapply(all.maps, function(x) mean(is.na(x))))>0.24)

#explore that one
head(all.maps[[41]])
sort(apply(all.maps[[41]], 2, function(x){length(which(is.na(x)))}))

#can get study titles
lapply(all.maps, function(x) unique(x[,"TITLE"]))
#great, some studies gave each sample a unique  title... annoying
which(lapply(all.maps, function(x) length(unique(x[,"TITLE"])))>1)
#oh, only 2


#fix titles
sort(all.maps[[9]][,"TITLE"])
#so this is an Excel 'draw down error', 
#started with "Intertidal microbes 16s for 2009 and 2010"
all.maps[[9]][,"TITLE"]<-"Intertidal microbes 16s for 2009 and 2010"

#the other one
sort(all.maps[[16]][,"TITLE"])
#same issue
all.maps[[16]][,"TITLE"]<-"EPOCA_Svalbard2018"

#and titles
unlist(lapply(all.maps, function(x) unique(x[,"TITLE"])))
#ah good

#make a little table...
study_count_NA_table<-
	data.frame(TITLE=unlist(lapply(all.maps, function(x) unique(x[,"TITLE"]))),
						 No_SAMPLES=unlist(lapply(all.maps, nrow)),
						 No_COLUMNS=unlist(lapply(all.maps, ncol)),
						 No_NA=unlist(lapply(all.maps, function(x) length(which(is.na(x))))),
						 No_VALUES=unlist(lapply(all.maps, function(x) nrow(x)*ncol(x))),
						 PERCENT_NA=(round(unlist(lapply(all.maps, function(x) mean(is.na(x)))),2)*100)
	)
#sort by TITLE
study_count_NA_table<-study_count_NA_table[order(as.character(study_count_NA_table$TITLE)),]

#export to subfolder of working directory, remove row names
write.csv(study_count_NA_table, 
					file.path(paste(getwd(), 
													"outputs/study_count_NA_table.csv", sep="/")), 
					row.names=FALSE)


#now a little more depth...
#check out column names
study_column_table<-sort(table(unlist(lapply(all.maps, colnames))), decreasing=TRUE)

#clean up for exporting
study_column_table<-as.data.frame(study_column_table)
#study names are currently row names, add as column
study_column_table$COLUMN_NAMES<-rownames(study_column_table)
#rename count column
colnames(study_column_table)[1]<-"No_STUDIES"
#replace row names
rownames(study_column_table)<-as.character(1:nrow(study_column_table))
#reorder
study_column_table<-study_column_table[,c(2,1)]
#subset to fit on MSword page, place side by side
study_column_table<-cbind(study_column_table[1:38,], study_column_table[39:76,])

#add a little note, how many columns not shown...
study_column_table<-rbind(study_column_table, 
													c(paste("*", 
																	paste(as.character(length(unique(unlist(lapply(all.maps, colnames))))-nrow(study_column_table)*2), 
																				"columns not shown", sep=" "), 
																	sep=""), "", "", ""))

#this is handy...
write.csv(study_column_table, 
					file.path(paste(getwd(), 
													"outputs/study_column_table.csv", sep="/")), 
					row.names=FALSE)

#explore most common column values
com.col<-sort(table(unlist(lapply(all.maps, colnames))), decreasing=TRUE)
#subset columns that are in every study
com.col<-com.col[which(com.col==49)]
com.col<-names(com.col)
lapply(all.maps, function(x) head(x[ ,which(colnames(x) %in% com.col)], 2))

#pick through each study and find which common columns have different classes
map.class.tmp<-list()
map.class<-list()

for(i in 1:length(all.maps)){
	for(j in i:length(all.maps)){
		all.maps.cnames<-colnames(all.maps[[i]]
															[which(unlist(lapply(all.maps[i], colnames)) %in%  
																		 	unlist(lapply(all.maps[j], colnames)))])
		dat1.class<-lapply(all.maps[[i]][all.maps.cnames], class)
		dat2.class<-lapply(all.maps[[j]][all.maps.cnames], class)
		dat1.class<-unlist(dat1.class)
		dat2.class<-unlist(dat2.class)
		map.class.tmp[[j]]<-cbind(dat1.class[which(!dat1.class==dat2.class)], 
													dat2.class[which(!dat2.class==dat1.class)])
	}
	map.class[[i]]<-map.class.tmp
}
map.class[[1]]
#volumous, but works

#now can pick through and see what the differences are...
map.class[[1]][6]
head(all.maps[[1]][c("RUN_DATE", "DEPTH", "EXPERIMENT_DESIGN_DESCRIPTION")])
head(all.maps[[6]][c("RUN_DATE", "DEPTH", "EXPERIMENT_DESIGN_DESCRIPTION")])

#see if can get all the colnames with problems
class(map.class[[1]][6])
rownames(map.class[[1]][[6]])

tmp.rows<-list()
col.class.dif<-list()

for(i in 1:length(map.class)){
	tmp.map<-map.class[[i]]
	for(j in 1:length(map.class[[i]])){
		tmp.rows[[i]]<-list(rownames(tmp.map[[j]]))
	}
	col.class.dif[[i]]<-list(tmp.rows[[i]])
}
unlist(col.class.dif)
unique(unlist(col.class.dif))
sort(table(unlist(col.class.dif)))
#so it is mostly depth and a few others

#look at just these columns across studes
col.class.dif<-unique(unlist(col.class.dif))

l1 <- vector('list', 11)

for(i in 1:length(col.class.dif)){
	
	#get example of values in column for each data frame
	check.val<-lapply(all.maps, function(x) {head(x[colnames(x) %in% 
							unique(unlist(col.class.dif))[i]], 1)})
	#give elements names
	#names(check.val)<-paste("map", 1:49, sep="")
	#replace empty data frames with NA
	check.val[which(lapply(check.val, length)==0)]<-"No Field"
	
	#build table
	l1[[i]]<-unlist(check.val)
	names(l1[i])<-paste(col.class.dif[i], ".ex", sep="")
													 
}


l2 <- vector('list', 11)
for(i in 1:length(col.class.dif)){
	
	#get class of each column for each data frame
	check.class<-lapply(all.maps, function(x) {lapply(x[colnames(x) %in% 
								col.class.dif[i]], class)})
	#convert empty values to NA
	check.class[which(lapply(check.class, length)==0)]<-"No Field"
	
	l2[[i]]<-unlist(check.class)
	names(l2[i])<-paste(col.class.dif[i], ".class", sep="")
	
}

#create table
l1<-data.frame(l1)
colnames(l1)<-paste(col.class.dif, ".ex", sep="")

l2<-data.frame(l2)
colnames(l2)<-paste(col.class.dif, ".class", sep="")

study_class_diff_table<-data.frame(TITLE=unlist(lapply(all.maps, function(x) unique(x[,"TITLE"]))))
study_class_diff_table<-cbind(study_class_diff_table, l1, l2)
	
#sort by column name 
study_class_diff_table<-study_class_diff_table[, c("TITLE", sort(colnames(study_class_diff_table)[2:23]))]
#export
write.csv(study_class_diff_table, 
					file.path(paste(getwd(), 
													"outputs/study_class_diff_table.csv", sep="/")), 
					row.names=FALSE)

#nice thought, not much better...
#study_class_diff_table.melt<-melt(data.frame(TITLE=study_class_diff_table[,"TITLE"], l2), id.vars="TITLE")

#study_class_diff_table.cast<-dcast(study_class_diff_table.melt, variable~ TITLE, value.var="value")
#head(study_class_diff_table.cast)

#investigate all of the most common fields using this script...
col.comm<-c(study_column_table[1:40,1], study_column_table[,3])
col.comm<-col.comm[1:43]
#add contact fields (mostly missing...?)
col.comm<-c(col.comm, "PRINCIPAL_INVESTIGATOR_CONTACT", 
						"LAB_PERSON_CONTACT", "MOST_RECENT_CONTACT")

l1 <- vector('list', 43)

for(i in 1:length(col.comm)){
	
	#get example of values in column for each data frame
	check.val<-lapply(all.maps, function(x) {head(x[colnames(x) %in% 
										col.comm[i]], 1)})
	#give elements names
	#names(check.val)<-paste("map", 1:49, sep="")
	#replace empty data frames with NA
	check.val[which(lapply(check.val, length)==0)]<-"No Field"
	
	#build table
	l1[[i]]<-unlist(check.val)
	names(l1[i])<-paste(col.comm[i], ".ex", sep="")
	
}


l2 <- vector('list', 43)
for(i in 1:length(col.comm)){
	
	#get class of each column for each data frame
	check.class<-lapply(all.maps, function(x) {lapply(x[colnames(x) %in% 
																												col.comm[i]], class)})
	#convert empty values to NA
	check.class[which(lapply(check.class, length)==0)]<-"No Field"
	
	l2[[i]]<-unlist(check.class)
	names(l2[i])<-paste(col.comm[i], ".class", sep="")
	
}

#create table
l1<-data.frame(l1)
colnames(l1)<-paste(col.comm, ".ex", sep="")

l2<-data.frame(l2)
colnames(l2)<-paste(col.comm, ".class", sep="")

study_class_comm_col_table<-data.frame(TITLE=unlist(lapply(all.maps, function(x) unique(x[,"TITLE"]))))
study_class_comm_col_table<-cbind(study_class_comm_col_table, l1, l2)

#sort by column name 
study_class_comm_col_table<-study_class_comm_col_table[, c("TITLE", sort(colnames(study_class_comm_col_table)[2:87]))]
#export
write.csv(study_class_comm_col_table, 
					file.path(paste(getwd(), 
													"outputs/study_class_comm_col_table.csv", sep="/")), 
					row.names=FALSE)

#can create list of studies and key data to export as report?

EMP_metadata_issues<-lapply(all.maps, function(x) data.frame(Value=c(
	paste(unique(x[,"TITLE"]), collapse=" "),
	paste(unique(x[,colnames(x) %in% c("PRINCIPAL_INVESTIGATOR_CONTACT", "LAB_PERSON_CONTACT", "MOST_RECENT_CONTACT")]), collapse="; "),
	ifelse(TRUE %in% (colnames(x) %in% "COLLECTION_DATE"), paste(range(x[,"COLLECTION_DATE"]), collapse=" to "), "No Field"),
	ifelse(TRUE %in% (colnames(x) %in% "COLLECTION_DATE"), class(x[,"COLLECTION_DATE"]), "No Field"),
	ifelse(TRUE %in% (colnames(x) %in% "SAMP_SIZE"), paste(range(x[,"SAMP_SIZE"]), collapse=" to "), "No Field"),
	ifelse(TRUE %in% (colnames(x) %in% "DEPTH"), paste(range(x[,"DEPTH"]), collapse=" to "), "No Field"),
	ifelse(TRUE %in% (colnames(x) %in% "DEPTH"), class(x[,"DEPTH"]), "No Field"),
	ifelse(TRUE %in% (colnames(x) %in% "RUN_DATE"), paste(range(x[,"RUN_DATE"]), collapse=" to "), "No Field"),
	ifelse(TRUE %in% (colnames(x) %in% "RUN_DATE"), class(x[,"RUN_DATE"]), "No Field"),
	ifelse(TRUE %in% (colnames(x) %in% "SEQUENCING_METH"), paste(head(x[,"SEQUENCING_METH"], 1), collapse="  "), "No Field"),
	ifelse(TRUE %in% (colnames(x) %in% "LIBRARY_CONSTRUCTION_PROTOCOL"), paste(head(x[,"LIBRARY_CONSTRUCTION_PROTOCOL"], 1), collapse=" "), "No Field"),
	ifelse(TRUE %in% (colnames(x) %in% "PCR_PRIMERS"), paste(head(x[,"PCR_PRIMERS"], 1), collapse=" "), "No Field"),
	ifelse(TRUE %in% (colnames(x) %in% "PLATFORM"), paste(head(x[,"PLATFORM"], 1), collapse=" "), "No Field"),
	ifelse(TRUE %in% (colnames(x) %in% "RUN_CENTER"), paste(head(x[,"RUN_CENTER"], 1), collapse=" "), "No Field"),
	ifelse(TRUE %in% (colnames(x) %in% "SAMPLE_CENTER"), paste(head(x[,"SAMPLE_CENTER"], 1), collapse=" to "), "No Field"),
	ifelse(TRUE %in% (colnames(x) %in% "SAMPLE_LOCATION"), paste(head(x[,"SAMPLE_LOCATION"], 1), collapse=" "), "No Field"),
	ifelse(TRUE %in% (colnames(x) %in% "TARGET_GENE"), paste(head(x[,"TARGET_GENE"], 1), collapse=" "), "No Field"),
	paste("")),	
	#FIELD=c("TITLE","CONTACTS", "RUN_DATE","RUN_DATE.class", "DEPTH","DEPTH.class", "break"),  																	
	row.names=c("TITLE", 
							"CONTACTS", 
							"COLLECTION_DATE_range", 
							"COLLECTION_DATE_class", 
							"SAMP_SIZE_range",
							"DEPTH_range",
							"DEPTH.class", 
							"RUN_DATE_range",
							"RUN_DATE.class",
							"SEQUENCING_METH", 
							"LIBRARY_CONSTRUCTION", 
					  	"PCR_PRIMERS", 
							"PLATFORM", 
							"RUN_CENTER",  
							"SAMPLE_CENTER", 
							"SAMPLE_LOCATION", 
							"TARGET_GENE", 
							"")  																	
))

names(EMP_metadata_issues)<-unlist(lapply(all.maps, function(x) unique(x[,"TITLE"])))

#test
ifelse(TRUE %in% (colnames(all.maps[[1]]) %in% "RUN_CENTER"), paste(head(all.maps[[1]][,"RUN_CENTER"]), 1), "No Field")
ifelse(TRUE %in% (colnames(all.maps[[1]]) %in% "RUN_CENTER"), paste(range(all.maps[[1]][,"RUN_CENTER"]), collaps= ""), "No Field")

ifelse(TRUE %in% (colnames(all.maps[[1]]) %in% "RUN_CENTER"), class(all.maps[[1]][,"RUN_CENTER"]), "No Field")

#lapply(all.maps, function(x) paste(unique(x[,colnames(x) %in% c("PRINCIPAL_INVESTIGATOR_CONTACT", "LAB_PERSON_CONTACT", "MOST_RECENT_CONTACT")]), collapse="; "))
#all.maps[["EPOCA_Svalbard2018"]][,colnames(all.maps[["EPOCA_Svalbard2018"]]) %in% c("PRINCIPAL_INVESTIGATOR_CONTACT", "LAB_PERSON_CONTACT", "MOST_RECENT_CONTACT")]

#convert to data frame and export
EMP_metadata_issues.df<-do.call("rbind", EMP_metadata_issues)
EMP_metadata_issues.df$Field<-rep(c("TITLE", 
																		"CONTACTS", 
																		"COLLECTION_DATE_range", 
																		"COLLECTION_DATE_class", 
																		"SAMP_SIZE_range",
																		"DEPTH_range",
																		"DEPTH.class", 
																		"RUN_DATE_range",
																		"RUN_DATE.class",
																		"SEQUENCING_METH", 
																		"LIBRARY_CONSTRUCTION", 
																		"PCR_PRIMERS", 
																		"PLATFORM", 
																		"RUN_CENTER",  
																		"SAMPLE_CENTER", 
																		"SAMPLE_LOCATION", 
																		"TARGET_GENE", 
																		""), 49)

EMP_metadata_issues.df<-EMP_metadata_issues.df[,c("Field", "Value")]

write.csv(EMP_metadata_issues.df, 
					file.path(paste(getwd(), 
													"outputs/EMP_metadata_issues.csv", sep="/")), 
					row.names=FALSE)

#this seemed useful, but in the end not, may want later though
#http://stackoverflow.com/questions/13006909/export-a-list-of-matrices-nicely-to-the-same-worksheet-in-excel

#wb = loadWorkbook("matrix.xlsx", create = TRUE)
# Create a new sheet
#createSheet(wb, name = "mysheet")
# cumulative length (rows) of matrices
# +2 = 1 for list names, 1 for header row
#cumlen = cumsum(c(1, head(sapply(EMP_metadata_issues, nrow), n = -1) + 2))
# Write data rows (implicitly vectorized!)
#writeWorksheet(wb, data = EMP_metadata_issues, sheet = "mysheet", startRow = cumlen + 1, header = FALSE, rownames=lapply(EMP_metadata_issues, rownames))
# Write list names
#writeWorksheet(wb, data = as.list(names(EMP_metadata_issues)), sheet = "mysheet", startRow = cumlen, header = FALSE)
#saveWorkbook(wb)

#check out some really bad ones...
#tibetan_plateau_salt_lake_sediment
head(all.maps[["tibetan_plateau_salt_lake_sediment"]])
head(all.maps[[22]])
colnames(all.maps[["tibetan_plateau_salt_lake_sediment"]])
colnames(all.maps[[1]])
all.maps[["tibetan_plateau_salt_lake_sediment"]]["COLLECTION_DATE"]

#Kilauea geothermal soils and biofilms
head(all.maps[["Kilauea geothermal soils and biofilms"]])
colnames(all.maps[["Kilauea geothermal soils and biofilms"]])

dim(all.maps[["Environmental metagenomic interrogation of Thar desert microbial communities"]])
all.maps[["Environmental metagenomic interrogation of Thar desert microbial communities"]]