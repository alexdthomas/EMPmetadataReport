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

#name studies in list
names(all.maps)<-unlist(lapply(all.maps, function(x) unique(x[,"TITLE"])))

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

l1 <- vector('list', length(col.class.dif))

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


l2 <- vector('list', length(col.class.dif))
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
#select fields present in greater than half of the studies
col.comm<-names(which(sort(table(unlist(lapply(all.maps, colnames))), decreasing=TRUE)>length(all.maps)/2))

#add contact fields (mostly missing...?)
col.comm<-c(col.comm, "PRINCIPAL_INVESTIGATOR_CONTACT", 
						"LAB_PERSON_CONTACT", "MOST_RECENT_CONTACT")

l1 <- vector('list', length(col.comm))

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


l2 <- vector('list', length(col.comm))
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

colnames(all.maps[["Fermilab_spatial_study"]])
head(all.maps[["Fermilab_spatial_study"]])

colnames(all.maps[["Great Lake Microbiome"]])
head(all.maps[["Great Lake Microbiome"]])

write.csv(head(all.maps[["tibetan_plateau_salt_lake_sediment"]]), 
					file.path(paste(getwd(), 
													"outputs/tibetan_plateau_salt_lake_sediment_head.csv", sep="/")), 
					row.names=FALSE)

write.csv(head(all.maps[["Kilauea geothermal soils and biofilms"]]), 
					file.path(paste(getwd(), 
													"outputs/Kilauea_geothermal_soils_and_biofilms_head.csv", sep="/")), 
					row.names=FALSE)

write.csv(head(all.maps[["Fermilab_spatial_study"]]), 
					file.path(paste(getwd(), 
													"outputs/Fermilab_spatial_study_head.csv", sep="/")), 
					row.names=FALSE)

write.csv(head(all.maps[["Great Lake Microbiome"]]), 
					file.path(paste(getwd(), 
													"outputs/Great_Lake_Microbiome_head.csv", sep="/")), 
					row.names=FALSE)

#try to evaluate all other columns across studies
col.n.comm<-names(which(sort(table(unlist(lapply(all.maps, colnames))), decreasing=TRUE)<=length(all.maps)/2))

#as all spaces between words are '_' can split all column names
#by '_' and find common words for fields might indicate a pattern to data
col.n.comm.sp<-strsplit(col.n.comm, "_")
head(col.n.comm.sp)
col.n.comm.sp<-unlist(col.n.comm.sp)
col.n.comm.sp<-as.matrix(sort(table(col.n.comm.sp), decreasing=TRUE))
#this is useful...
#frequencies...
table(col.n.comm.sp)
#start with most frequent
col.n.comm.sp[col.n.comm.sp>8, ]
length(col.n.comm.sp[col.n.comm.sp>8, ])

col.n.comm[grep(paste(names(col.n.comm.sp[col.n.comm.sp>8, ]), collapse="|"), col.n.comm, ignore.case=TRUE)]
#99 columns contain the 9 most frequent words

#search colnames
col.n.comm[grep('TOT', col.n.comm, fixed=TRUE)]
#well here's some chemistry
col.chem<-col.n.comm[grep('TOT', col.n.comm, fixed=TRUE)]
col.chem<-col.chem[which(!col.chem %in% c("TOTAL_PHYTOPLANKTON_COUNT", "TOTAL_BACTERIA_COUNT", "DISTANCETOTRANSECTA_CM"))]
#look at metadata with chemistry
chem.maps<-lapply(all.maps, function(x) ifelse(TRUE %in% (colnames(x) %in% col.chem), head(x[,"TITLE"], 1), "No Field"))
chem.maps<-names(chem.maps[which(!chem.maps %in% "No Field")])
lapply(all.maps[chem.maps], function(x) head(x[colnames(x) %in% col.chem], 2))
#not many columns shared across studies, see if can focus
length(chem.maps) #14 studies have any of these columns
table(unlist(lapply(all.maps[chem.maps], function(x) colnames(x[colnames(x) %in% col.chem]))))
#10 have TOT_NITRO, only 7 have TOT_N_METH and only 1 has TOT_NITRO_UNITS
lapply(all.maps[chem.maps], function(x) head(x[colnames(x) %in% c("TOT_NITRO", "TOT_N_METH", "TOT_NITRO_UNITS")], 10))
#well, that instersting, values don't mean much without units and methods...

#export table for EMP chemistry fields
chem.maps.df<-lapply(all.maps[chem.maps], function(x) head(x[colnames(x) %in% c("TITLE", col.chem)], 1))
unlist(chem.maps.df)
melt(chem.maps.df)

chem.maps.df <- lapply(chem.maps.df, unlist)
chem.max <- max(sapply(chem.maps.df, length))
do.call(rbind, lapply(chem.maps.df, function(z)c(z, rep(NA, chem.max-length(z)))))

chem.maps.df<-Reduce(function(x, y) merge(x, y, all=TRUE), chem.maps.df)

chem.maps.df<-vector('list', length(chem.maps))

chem.maps.df<-lapply(all.maps[chem.maps], function(x) data.frame(
	TITLE=paste(unique(x[,"TITLE"]), collapse=" "),
	CONTACTS=paste(unique(x[,colnames(x) %in% c("PRINCIPAL_INVESTIGATOR_CONTACT", "LAB_PERSON_CONTACT", "MOST_RECENT_CONTACT")]), collapse="; "))	
)

chem.maps.df<-lapply(all.maps[chem.maps], "[[", )
	lapply(col.chem, function(i) ifelse(TRUE %in% (colnames(x) %in% col.chem[i]), paste(head(x[,"COLLECTION_DATE"], 1), collapse=" "), "No Field"))
	
	
for(i in 1:length(col.chem)){
			paste(col.chem[i]) = ifelse(TRUE %in% (colnames(x) %in% col.chem[i]), paste(head(x[,"COLLECTION_DATE"], 1), collapse=" "), "No Field")
		)
}
	
#search colnames
col.n.comm[grep('DATE', col.n.comm, fixed=TRUE)]
#yup, dates
col.date<-col.n.comm[grep('DATE', col.n.comm, fixed=TRUE)]
date.maps<-lapply(all.maps, function(x) ifelse(TRUE %in% (colnames(x) %in% col.date), head(x[,"TITLE"], 1), "No Field"))
date.maps<-names(date.maps[which(!date.maps %in% "No Field")])
length(date.maps)
#well, all 49 have some date field (knew that)
sort(table(unlist(lapply(all.maps[date.maps], function(x) colnames(x[colnames(x) %in% col.date])))), decreasing=TRUE)
#well, basically only COLLECTION_DATE and RUN_DATE, most of the rest look like Juan's

#search colnames
col.n.comm[grep('METH', col.n.comm, fixed=TRUE)]
#pretty diverse, already looked at SEQUENCING_METH and chemistry *METH
#some soil ones are ineresting
#WATER_CONTENT_SOIL_METH, SOIL_TYPE_METH, TEXTURE_METH
sort(table(unlist(lapply(all.maps, function(x) colnames(x[colnames(x) %in% c("WATER_CONTENT_SOIL_METH", "SOIL_TYPE_METH", "TEXTURE_METH")])))), decreasing=TRUE)
lapply(all.maps, function(x) head(x[colnames(x) %in% c("WATER_CONTENT_SOIL_METH", "SOIL_TYPE_METH", "TEXTURE_METH")]))
#well it's something...

#search colnames
col.n.comm[grep('SOIL', col.n.comm, fixed=TRUE)]
#I like soil...
col.soil<-col.n.comm[grep('SOIL', col.n.comm, fixed=TRUE)]
soil.maps<-lapply(all.maps, function(x) ifelse(TRUE %in% (colnames(x) %in% col.soil), head(x[,"TITLE"], 1), "No Field"))
soil.maps<-names(soil.maps[which(!soil.maps %in% "No Field")])
length(soil.maps)
#hmm, 9 studies have any fields with soil...
as.matrix(sort(table(unlist(lapply(all.maps[soil.maps], function(x) colnames(x[colnames(x) %in% col.soil])))), decreasing=TRUE))
#and most of them only occur once 
#again most common is WATER_CONTENT_SOIL (7), but only 4 have a method...
lapply(all.maps[soil.maps], function(x) head(x[colnames(x) %in% col.soil], 2))

#search colnames
col.n.comm[grep('CARB', col.n.comm, fixed=TRUE)]
#hmm, maybe TOT_ORG_CAB and ORG_CARB are duplicates
col.carb<-col.n.comm[grep('CARB', col.n.comm, fixed=TRUE)]
carb.maps<-lapply(all.maps, function(x) ifelse(TRUE %in% (colnames(x) %in% col.carb), head(x[,"TITLE"], 1), "No Field"))
carb.maps<-names(carb.maps[which(!carb.maps %in% "No Field")])
length(carb.maps)
#15 studies have fields related to carbon
as.matrix(sort(table(unlist(lapply(all.maps[carb.maps], function(x) colnames(x[colnames(x) %in% col.carb])))), decreasing=TRUE))
#the most common is TOT_ORG_CARB (8), however only 1 TOT_ORG_CARB_UNITS and no methods...

#search colnames
col.n.comm[grep('COUNT', col.n.comm, fixed=FALSE)]
#note "COUNTY" is less frequent, not "COUNTRY"
#includes many microbial counts, maybe interesting, but not right now...

#search colnames
col.n.comm[grep('TEMP', col.n.comm, fixed=FALSE)]
#hey, here are some potential duplicates
col.temp<-col.n.comm[grep('TEMP', col.n.comm, fixed=FALSE)]
temp.maps<-lapply(all.maps, function(x) ifelse(TRUE %in% (colnames(x) %in% col.temp), head(x[,"TITLE"], 1), "No Field"))
temp.maps<-names(temp.maps[which(!temp.maps %in% "No Field")])
length(temp.maps)
#20 studies have temp data
as.matrix(sort(table(unlist(lapply(all.maps[temp.maps], function(x) colnames(x[colnames(x) %in% col.temp])))), decreasing=TRUE))
#what is differenve between "TEMP" and all other "*TEMP*"
#ANNUAL_SEASON_TEMP vs. ANNUAL_TEMP vs. SEASON_TEMP
#SOIL_TEMP vs. MEAN_SOIL_TEMP_DAY
lapply(all.maps[temp.maps], function(x) head(x[colnames(x) %in% col.temp], 2))

#search colnames
col.n.comm[grep('SEASON', col.n.comm, fixed=FALSE)]
#some duplicates from TEMP, and other climate data...
col.season<-col.n.comm[grep('SEASON', col.n.comm, fixed=FALSE)]
season.maps<-lapply(all.maps, function(x) ifelse(TRUE %in% (colnames(x) %in% col.season), head(x[,"TITLE"], 1), "No Field"))
season.maps<-names(season.maps[which(!season.maps %in% "No Field")])
length(season.maps)
#8 studies have season data
as.matrix(sort(table(unlist(lapply(all.maps[season.maps], function(x) colnames(x[colnames(x) %in% col.season])))), decreasing=TRUE))
lapply(all.maps[season.maps], function(x) head(x[colnames(x) %in% col.season], 2))
#wonder what the methods were to find these seasonal values, but interesting

#search colnames
col.n.comm[grep('UNIT', col.n.comm, fixed=FALSE)]
#some duplicates from TEMP, and other climate data...
col.unit<-col.n.comm[grep('UNIT', col.n.comm, fixed=FALSE)]
unit.maps<-lapply(all.maps, function(x) ifelse(TRUE %in% (colnames(x) %in% col.unit), head(x[,"TITLE"], 1), "No Field"))
unit.maps<-names(unit.maps[which(!unit.maps %in% "No Field")])
length(unit.maps)
#8 studies have unit data
as.matrix(sort(table(unlist(lapply(all.maps[unit.maps], function(x) colnames(x[colnames(x) %in% col.unit])))), decreasing=TRUE))
lapply(all.maps[unit.maps], function(x) head(x[colnames(x) %in% col.unit], 2))
#can at least find sudies that have units...



