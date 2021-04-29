library(drc)
# read data
df=read.csv("stability_analysis.csv",stringsAsFactors = F)

## convert negative values resulting from background subtraction
## to 1e-5
df$relative_gfp[df$relative_gfp<1e-5]=1e-5

## fit the 3 parameter log-logistic function
temp_res=drm( relative_gfp~temperature, virus, data = df, fct = LL.3())

## plot data
plot(temp_res, col = c("blue", "red", "black","green"))
# obtain the ED50 (temperature at which 50% reduction is observed)
ED(temp_res,50)

## compare the ED50 to test for statistical significance.
EDcomp(temp_res, c(50,50))


# plot IC50 and error
ed50=data.frame(ED(temp_res,50))
ed50$virus=gsub("e:|:50","",row.names(ed50))


## plot ED50 and error

library(ggplot2)
ggplot(ed50,aes(x=virus,y=Estimate))+ geom_point()+
  geom_errorbar(aes(ymin=Estimate-Std..Error,
                    ymax=Estimate+Std..Error),width=.2,
                position=position_dodge(.9))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5))+
  ylab("50% inhibitory temperature")

