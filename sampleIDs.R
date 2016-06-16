####CAWR blast data
##goals
##merging IDs with sample data
##convert to presence/absence (abundance counts dont matter)
##determine what was in each sample
##then this should be related to the arthropod data?
library(dplyr)
library(ggplot2)
library(reshape2)
library(devtools)

firstids<-read.csv("/Users/colleennell/Documents/R/CAWR/data/CAWR_OTU_tree.csv")
samples<-read.csv("/Users/colleennell/Documents/R/CAWR/data/CAWR_samples.csv")


##################################################################################################
##cleaning of blast data, adding information

##create new variable that informs of the taxonomic resolution for each taxa ID
blast<-read.csv("/Users/colleennell/Documents/R/CAWR/data/CAWR_OTU_tree_june16.csv")
#blast$idres<-match(blast$Taxa.ID,blast$ORDER)##this returns a variable with '10' indicating a match and NA where not
str(blast)##can only compare factors if have same levels..make all character?
blast$Taxa.ID<-as.character(blast$Taxa.ID)
blast$ORDER<-as.character(blast$ORDER)
blast$FAMILY<-as.character(blast$FAMILY)
blast$GENUS<-as.character(blast$GENUS)
blast$SPECIES<-as.character(blast$SPECIES)
blast$reso<-ifelse(blast$Taxa.ID == blast$SPECIES,"species",
                   ifelse(blast$Taxa.ID == blast$GENUS,"genus",
                          ifelse(blast$Taxa.ID==blast$FAMILY,"family",
                                 ifelse(blast$Taxa.ID==blast$ORDER,"order",
                                        ifelse(blast$Taxa.ID=="unknown","none",
                                               ifelse(blast$Taxa.ID==blast$Taxa.Order,"order",
                                                      ifelse(blast$Taxa.ID==blast$Taxa.Family,"family",NA
                                        )))))))
##manually double check that these are right...what does it do for the unknowns?
check<-blast[,c("OTU","Taxa.ID","reso","Taxa.Family","Taxa.Order","ORDER","FAMILY","GENUS","SPECIES")]
View(check)
##amazing
###new Taxa.ID based on 98% similarity filter from BOLD
blast$Taxa.ID.98<-ifelse(blast$Similarity.BOLD>=97.8 | blast$Sequence.identity....>=97.8,blast$Taxa.ID,
                         ifelse(blast$Similarity.BOLD==NA & blast$Sequence.identity....>=97.8,blast$Taxa.ID,
                                ifelse(blast$Similarity.BOLD==NA | blast$Sequence.identity....==NA,"none", "none")))
View(blast)##I guess this is fine, but is this really the correct filter anyways?
write.csv(blast,"/Users/colleennell/Documents/R/CAWR/CAWR_OTU_finalids.csv")

#####################################################################
##combining with sample data
ids<-read.csv("/Users/colleennell/Documents/R/CAWR/CAWR_OTU_finalids.csv")###this is including bad matches 'unknown' ids
samples<-read.csv("/Users/colleennell/Documents/R/CAWR/data/CAWR_samples.csv")

##need to pull otu id's from ids df
##current format: >denovo127 S8_20500, need denovo127 to then merge with samples df
str(ids) ##needs to be a character vector?
ids$ID<-as.character(ids$ID)
##extract first element from string split (on ' ') into a new col
#ids$OTU<-unique(rapply(strsplit(ids$ID,' '), function(x) head(x, 1)))
##remove '>' at the beginning of each OTU id
#ids$OTU<-gsub('>','',ids$OTU)

##option 2: split the OTU### and assign 'denovo' in front of the numbers
ids$denovo<-gsub('OTU','denovo',ids$OTU)

##join df's
samples<-left_join(samples, ids[,c("denovo","Taxa.ID","Taxa.Order","Taxa.Family","reso","Taxa.ID.98")], by=c("OTU"="denovo"))
View(samples)
##collapse df into orders, look at overall sample composition by orders
##are there multiple OTUs of same Taxa.ID in a single sample? which ones?

##redo 'total' col using function...never trust excel
drops<-c("Totals")
samples<-samples[,!(names(samples)%in% drops)]##remove col
##make new variable...abundances are irrelevant so just do it as a count
##first lets make a vector of all the samples so they can be referenced easily
guts<-select(samples,starts_with("S"))
str(guts)
##and for the TaxaIDs!
taxa<-select(samples,starts_with("Taxa."))
str(taxa)
##total it up
samples$total_abun<-rowSums(guts,na.rm=T)
##count of how many samples present in
samples$total_samples<-rowSums(guts>0)

write.csv(samples,"/Users/colleennell/Documents/R/CAWR_sample_ids.csv")

############################
##starting here is building the rmd 'CAWR_gut.Rmd'

##table of the taxa ID's 
##add variables for genus and species where applicable
str(ids)
ids$Taxa.Genus<-as.factor(ifelse(ids$reso=='genus',paste(ids$Taxa.ID),
                       ifelse(ids$reso=='species',paste(ids$GENUS),NA )))
ids$Taxa.Sp<-as.factor(ifelse(ids$reso=='species',paste(ids$Taxa.ID),NA ))
ids$Taxa.Family<-ifelse(ids$Taxa.Family=='unknown'| ids$Taxa.Family=='None',NA,paste(ids$Taxa.Family))

id.table<-ids%>%
  group_by(Taxa.ID,Taxa.Order,Taxa.Family,Taxa.Genus,Taxa.Sp)%>%
  summarize(n_otu=length(OTU))
View(id.table)

######
#histogram of ids by sample count
guts<-read.csv("/Users/colleennell/Documents/R/CAWR_sample_ids.csv")##this is the above 2 combined with a few more cols
View(guts)
sample.matrix<-select(guts,starts_with("S"))
library(dplyr)

histy<-guts%>%
  group_by(Taxa.Order,Taxa.Family)%>%
  summarize(abun_sample=sum(total_abun))
View(histy)


#################
##lets look at each sample composition
##

##make the id x sample matrix transposed and presence/absence
sample.matrix<-select(guts,starts_with("S"))
order<-guts$Taxa.Order
sample.matrix<-cbind(order,sample.matrix)
View(sample.matrix)

###p/a matrix


library(reshape2)
melted.mat<-melt(sample.matrix)
View(melted.mat)
melted.mat$value<-ifelse(melted.mat$value>0,1,0)
melted.mat<-melted.mat%>%
  group_by(order,variable)%>%
  summarize(n_samps = sum(value))

sample.matris.t<-dcast(melted.mat,variable~order)
View(sample.matris.t)
samps<-sample.matris.t[,1]
rowsum<-rowSums(sample.matris.t[,-1])
prop.mat<-mutate_each(sample.matris.t[,-1],funs(./rowsum*100))
prop.mat<-cbind(samps,prop.mat)
neword<-melt(prop.mat)

samp.comp<-ggplot(neword,aes(x=samps,y=value,fill=variable))+
  geom_bar(stat="identity")+labs(x="Sample",y="Arthropod Composition")+theme(panel.background = element_rect(fill='none'))
                       
samp.comp

##an avg sample

std <- function(x) sd(x)/sqrt(length(x))

newsamp<-neword%>%
  group_by(variable)%>%
  summarize(avg.per = mean(value), se = std(value))


avgplot<-ggplot(newsamp,aes(x=variable,y=avg.per))+
  geom_bar(stat="identity")+
  geom_errorbar(aes(ymin=avg.per-se,ymax=avg.per+se),width=.2)+
  theme_minimal()+
  theme(panel.background = element_rect(fill = "transparent",colour = NA))+
  scale_fill_brewer(name="Order",palette = 'Spectral')+
  labs(x="Order",y="Average Proportion of Diet")
                                                                                           
avgplot




