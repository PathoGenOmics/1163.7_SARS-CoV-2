# read data
df=read.csv("vaccine_titers.csv",stringsAsFactors = F)
## convert negative values resulting from background subtraction
## to 1e-5
df$relative_gfp[df$relative_gfp<1e-5]=1e-5
# convert the Ab concentration to the reciprocal log
df$concentration=1/df$concentration
df$concentration=log(df$concentration,10)

library(drc)

## cutoff for ic calculation
cutoff=50 

# make res df
res.ic=data.frame( sera=NA,
                     virus=NA,
                     ic=NA,
                     error.ic=NA,
                     ll2=NA,
                     ll3=NA,
                     ll4=NA,
                     testl2l3=NA,
                     testl3l4=NA,
                     model=NA)


## pdf for figures of fit 
pdf(file = "IC.pdf")
par(mfrow=c(3,4))
# get sera to iterate on
sera=unique(df$sample)

## iterates on sera, then on virus, 
## picks model with lowest # param that fits data
## calculates IC and graphs fit
for (i in sera){ 
  sera.df=df[df$sample==i,]## select  each sera
  
  for (v in unique(sera.df$virus)){
    ## get calculation for each virus
    virus=sera.df[sera.df$virus==v,]
    
    ## calculate different model parameter
    modl.2p=drm(data=virus, relative_gfp~concentration,fct = LL.2(),type = "continuous")
    modl.3p=drm(data=virus, relative_gfp~concentration,fct = LL.3(),type = "continuous")
    modl.4p=drm(data=virus, relative_gfp~concentration,fct = LL.4(),type = "continuous")
    
    # compare to other models
    comp=data.frame(mselect( modl.2p, list( LL.4(), LL.3(), LL.2()), icfct=AIC))
    # get AIC fit for each model
    ll2=comp["LL.2","IC"]
    ll3=comp["LL.3","IC"]
    ll4=comp["LL.4","IC"]
    
    ## anova to test if 3 better than 2
    t23=anova(modl.2p, modl.3p)[2,5]
    ## anova to test if 4 better than 3
    t34=anova(modl.3p, modl.4p)[2,5]
    
    
    if (t23< 0.05/12)  { # if 3 param is better than 2 param model
      mod="LL3";
      final.model=modl.3p;
      eds.final=ED(modl.3p, cutoff)}
    
    if (t23 > 0.05/12)  { # if 3 param is not better than 2 param model
      mod="LL2";
      final.model=modl.2p;
      eds.final=ED(modl.2p, cutoff);
    }
    
    if(t34< 0.05/12)  { # if 4 param is better than 3 param
      mod="LL4"
      final.model=modl.4p;
      eds.final=ED(modl.4p, cutoff)}
    
    ## get IC formatted
    ic=paste0(format(eds.final[1],digits = 3),
                " (",format(eds.final[2],digits = 3),")")
    
    # add results to result dataframe
    results=c(i, v,eds.final[1],eds.final[2],
              ll2,ll3,ll4,t23,t34,mod)
    res.ic=rbind(res.ic,results)
    rm(results)
    
    if (v=="S2-G614-V222") {v="20E"}
    if (v=="S2-G614-V222-Y1163-V1167") {v="VOI1163.7"}
    plot(final.model,
         ylim = c(0,1.5),
         col = "black",
         type="all",
         cex=1,
         xlab=paste0("Reciprocal dilution (Log10)"),
         ylab="GFP \n (% of untreated)",
         legend = F,
         # main=paste0("virus ",d))
         main=paste0(i,"\n"," ", v," ", ic),
         cex.main=.75, cex.lab=.75, cex.axis=0.75)
    
    rm(v,virus,comp,final.model, eds.final,
       modl.2p,modl.3p,modl.4p,ic,ll2,ll3,ll4,mod,
       t23,t34)
    
    
  }
  rm(i,sera.df)
}
dev.off()
#remove empty row
res.ic=res.ic[-1,]
# convet data to correct type
res.ic$ic=as.numeric(res.ic$ic)
res.ic$error.ic=as.numeric(res.ic$error.ic)
# filter out extra columns not needed
res.ic=res.ic[,colnames(res.ic)%in%c("virus","ic","error.ic","sera","model")]
res.ic=res.ic[,c("virus","ic","error.ic","sera","model")]
## write IC results
write.csv(res.ic,paste0("IC",cutoff,"_min_1e5_recip.csv"),row.names = F)

### stats for differences
### t test on IC data (already log transformed)
# test normality # p>0.05, so normal
shapiro.test(res.ic$ic)

t.test(res.ic$ic[res.ic$virus=="S2-G614-V222"],
            res.ic$ic[res.ic$virus=="S2-G614-V222-Y1163-V1167"],
            paired=T)
## just in case, we can look at non-parametric
### non-parametric test on IC data (already log transformed)
wilcox.test(res.ic$ic[res.ic$virus=="S2-G614-V222"],
       res.ic$ic[res.ic$virus=="S2-G614-V222-Y1163-V1167"],
       paired=T)


##### Plots
library(ggplot2)
ggplot(res.ic,aes(x=virus,y=ic))+ geom_point()+
  geom_errorbar(aes(ymin=ic-error.ic,ymax=ic+error.ic),width=.2,
                position=position_dodge(.9))+
  facet_grid(.~sera)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5))+
  ylab("% virus neutralization titer")



########## calculate ratio and means/range for effect size
### undo log for titer to calculate ratio differences in neutralization  
res.ic$titer=10^(res.ic$ic)

## result dataframe
res.ratio=data.frame(cbind(sera=unique(res.ic$sera),ratio=rep(NA,length(unique(res.ic$sera)))))

## iterate on each sera and get ratio
for( x in unique(res.ic$sera)){
  # x=res.ratio[1,1]
  E20=res.ic$titer[which(res.ic$virus=="S2-G614-V222" &
                             res.ic$sera==x)];
  VOI1163.7=res.ic$titer[res.ic$virus=="S2-G614-V222-Y1163-V1167" &
                             res.ic$sera==x];
  res.ratio$ratio[res.ratio$sera==x]= E20/VOI1163.7
  rm(E20,VOI1163.7)}

res.ratio$ratio=as.numeric(res.ratio$ratio)

## summary stats for each 
mean(res.ratio$ratio)
sd(res.ratio$ratio)
range(res.ratio$ratio)

