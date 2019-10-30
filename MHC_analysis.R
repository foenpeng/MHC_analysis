library(data.table)
library(ggplot2)
library(tidyr)
library(cowplot)
setwd("/Users/pengfoen/OneDrive - University of Connecticut/MHC")
MHC_meta <- fread("./7. Metadata and MHC data combined/7.0.BayesSummaryMergedMetadata.csv")

########### regression per lake

x_var<-"N_aas"
MHC_meta_cleaned<-MHC_meta
MHC_meta_cleaned[,x:=get(x_var)]
MHC_meta_cleaned[,sqx:=x^2]
# create site_ID from sample ID, becasue site name has many NAs (It turned out these NAs are the ones without metadata)
MHC_meta_cleaned[,site_ID:=unlist(lapply(SAMPLE_ID, function(x) strsplit(x, "(?=[A-Za-z])(?<=[0-9])|(?=[0-9])(?<=[A-Za-z])", perl=TRUE)[[1]][1]))]
# remove the rows that does not contain information for regression
MHC_meta_cleaned<-na.omit(MHC_meta_cleaned, cols=c("x", "Parasite.richness"))
# remove lakes with too few data for regression
MHC_meta_cleaned<-MHC_meta_cleaned[,if(.N>1) .SD, by="site_ID"][DEPTH_ALLELES>300]

y_list<-c("Parasite.richness","parasiteSWdiversity","RelativeInfectionIntensity")
res_combined<-list()
for(y_var in y_list){
      MHC_meta_cleaned[,y:=get(y_var)]
  
      setkey(MHC_meta_cleaned,"site_ID")
      fits <- lapply(unique(MHC_meta_cleaned$site_ID),
                    function(z)lm(y~x+sqx,data=MHC_meta_cleaned[J(z),],na.action=na.omit))
      
      # coefficients
      res<-sapply(fits, function(x)c(coef=summary(x)$coef[,1],p=round(summary(x)$coef[,4],3), t=round(summary(x)$coef[,3],3),rsq=summary(x)$r.squared))
      colnames(res)<-unique(MHC_meta_cleaned$site_ID)

      res_combined[[y_var]]<-res
}

coef.sqx_all<-rapply(res_combined, how="list",f = function(x) x["coef.sqx",]) 
coef.sqx_all<-as.data.frame(do.call("cbind",coef.sqx_all))
coef.sqx_all[,"site_ID"]<-row.names(coef.sqx_all)

t.sqx_all<-rapply(res_combined, how="list",f = function(x) x["t.sqx",]) 
t.sqx_all<-as.data.frame(do.call("cbind",t.sqx_all))
t.sqx_all[,"site_ID"]<-row.names(t.sqx_all)

p.sqx_all<-rapply(res_combined, how="list",f = function(x) x["p.sqx",]) 
p.sqx_all<-as.data.frame(do.call("cbind",p.sqx_all))
p.sqx_all[,"site_ID"]<-row.names(p.sqx_all)

# Q-Q plots
par(mfrow=c(3,2),mar=c(1,1,1,1))
lapply(fits,plot,2)

par(mfrow=c(3,2),mar=c(1,1,1,1))
sapply(unique(MHC_meta_cleaned$site_ID),function(z) plot(y~x,data=MHC_meta_cleaned[J(z),]))

# alternative way to run regression
#MHC_meta[,summary(lm(y~x+sqx,na.action=na.omit))$coefficients[3,"Pr(>|t|)", drop=F],by=site_name]
#MHC_meta_cleaned[,as.list(coef(lm(y ~ x + sqx,na.action=na.omit))), by="site_ID"]


# quadratic regression for all data
summary(lm(Parasite.richness~x+sqx, data=MHC_meta_cleaned))
ggplot(MHC_meta_cleaned, aes(group=N_aas,x=N_aas, y=Parasite.richness)) + 
  geom_boxplot(outlier.shape=NA) + geom_jitter(shape=16, position=position_jitter(0.2)) 

MHC_meta_cleaned[,pop_mean_rich:=mean(Parasite.richness,na.rm=T),by="site_ID"]
MHC_meta_cleaned[,residual_rich:=Parasite.richness-pop_mean_rich]
summary(lm(residual_rich~x+sqx, data=MHC_meta_cleaned))
ggplot(MHC_meta_cleaned, aes(group=N_aas,x=N_aas, y=residual_rich)) + 
  geom_boxplot(outlier.shape=NA) + geom_jitter(shape=16, position=position_jitter(0.2))

model1 <- lm(Parasite.richness~x+sqx +site_ID + x*site_ID + sqx*site_ID + log_std_length, data=MHC_meta_cleaned)
stepAIC(model1)
anova(model1)

##### regression on lake parameter
lake_info<-MHC_meta_cleaned[,lapply(.SD, function(x) mean(x, na.rm=T)),.SDcols=c(23:26,28,179,116:118,13),by=site_ID]
lake_coef_data<-lake_info[coef.sqx_all,on="site_ID"]
lambda_list<-c("lambda.richness","lambda.SW","lambda.Intensity")
setnames(lake_coef_data,c("i.Parasite.richness","i.parasiteSWdiversity","i.RelativeInfectionIntensity"), lambda_list)

for(i in lambda_list){
  #print(sapply(colnames(lake_coef_data)[2:6],function(x) summary(lm(paste0(i,"~",x),data=lake_coef_data))$coefficients[2,c("Estimate","Pr(>|t|)"), drop=F]))
  lake_fits<-sapply(colnames(lake_coef_data)[2:8],function(x) lm(paste0(i,"~",x),data=lake_coef_data))
  print(sapply(lake_fits, function(x)c(p=round(summary(x)$coef[2,4],3),coef=summary(x)$coef[2,1], t=round(summary(x)$coef[2,3],3),rsq=summary(x)$r.squared)))
}
summary(lm(lambda.richness~elevation_m + maximum_depth_m  +surface_area_ha +distance_from_ocean_km + benthicdiet.score + Parasite.richness, data=lake_coef_data  ))

# examine effect of lake parameter on parasite diversity
lake_fits<-sapply(colnames(lake_coef_data)[2:6],function(x) lm(paste0("Parasite.richness","~",x),data=lake_coef_data))
print(sapply(lake_fits, function(x)c(p=round(summary(x)$coef[2,4],3),coef=summary(x)$coef[2,1], t=round(summary(x)$coef[2,3],3),rsq=summary(x)$r.squared)))

# examine effect of lake parameter on MHC diversity


################## Plot results #####################
data_plot<-setDT(gather(t.sqx_all,y_list,key="parameter_name",value="t"))
p_plot<-setDT(gather(p.sqx_all,y_list,key="parameter_name",value="p_value"))
data_plot<-data_plot[p_plot,on=c("site_ID","parameter_name")]
data_plot[,sig:=ifelse(p_value<0.05,site_ID,NA)]
p<-list()
for(i in 1:length(y_list)){
  p[[i]]<-ggplot(data_plot[parameter_name==y_list[i]], aes(x=parameter_name, y=t)) + 
    geom_boxplot(outlier.shape=NA)+ coord_flip() + geom_jitter(shape=16, position=position_jitter(0.2)) +
    geom_text(aes(label = sig), na.rm = TRUE)
  
}
plot_grid(p[[1]],p[[2]],p[[3]],labels = "AUTO",nrow=3)


