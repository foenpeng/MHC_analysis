library(data.table)
setwd("/Users/pengfoen/OneDrive - University of Connecticut/MHC")
MHC_meta <- fread("./7. Metadata and MHC data combined/7.0.BayesSummaryMergedMetadata.csv")
# parasiteSWdiversity, RelativeInfectionIntensity
#x_var<-"Parasite.richness"
#x_var<-"parasiteSWdiversity"
x_var<-"RelativeInfectionIntensity"
y_var<-"N_aas"
MHC_meta_cleaned<-MHC_meta
MHC_meta_cleaned[,y:=get(y_var)]
MHC_meta_cleaned[,x:=get(x_var)]
MHC_meta_cleaned[,sqx:=x^2]
# create site_ID from sample ID, becasue site name has many NAs (It turned out these NAs are the ones without metadata)
MHC_meta_cleaned[,site_ID:=unlist(lapply(SAMPLE_ID, function(x) strsplit(x, "(?=[A-Za-z])(?<=[0-9])|(?=[0-9])(?<=[A-Za-z])", perl=TRUE)[[1]][1]))]
# remove the rows that does not contain information for regression
MHC_meta_cleaned<-na.omit(MHC_meta_cleaned, cols=c("x", "y"))
# remove lakes with too few data for regression
MHC_meta_cleaned<-MHC_meta_cleaned[,if(.N>1) .SD, by="site_ID"]



setkey(MHC_meta_cleaned,"site_ID")
fits <- lapply(unique(MHC_meta_cleaned$site_ID),
              # function(z)lm(hp~cyl+factor, data=dt[J(z),], y=T))
              function(z)lm(y~x+sqx,data=MHC_meta_cleaned[J(z),],na.action=na.omit))
# coefficients
res<-sapply(fits, function(x)c(coef=summary(x)$coef[,1],p=round(summary(x)$coef[,4],3), rsq=summary(x)$r.squared))
colnames(res)<-unique(MHC_meta_cleaned$site_ID)
sig<-res[,res["p.sqx",]<0.05]
  
# Q-Q plots
par(mfrow=c(3,2),mar=c(1,1,1,1))
lapply(fits,plot,2)

par(mfrow=c(3,2),mar=c(1,1,1,1))
sapply(unique(MHC_meta_cleaned$site_ID),function(z) plot(y~x,data=MHC_meta_cleaned[J(z),]))

# alternative way to run regression
#MHC_meta[,summary(lm(y~x+sqx,na.action=na.omit))$coefficients[3,"Pr(>|t|)", drop=F],by=site_name]
#MHC_meta_cleaned[,as.list(coef(lm(y ~ x + sqx,na.action=na.omit))), by="site_ID"]
