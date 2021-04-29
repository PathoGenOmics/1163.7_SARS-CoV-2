library(tidyverse)
library(broom)
library(drc)
library(modelr)
library(egg) 

set_wd <- function(){
  .expr <- quote({
    print(environment())
  })
  .envir <- sys.frame(which = -1)
  eval(expr = .expr,
       envir = .envir)
}

set_wd()
df=read.csv("stability_analysis.csv",stringsAsFactors = F)
df$relative_gfp[df$relative_gfp<1e-5]=1e-5
df$Mutant <- factor(df$Mutant, levels = c("Wuhan","D614G","20E: A222V_D614G", "20E: A222V_D614G_D1163Y_G1167V"))

# view the data in normal and log scale for temperature
p1 <- df %>% ggplot() + geom_point(aes(temperature, relative_gfp, color = Mutant)) + theme_bw()
p2 <- df %>% ggplot() + geom_point(aes(log(temperature), relative_gfp, color = Mutant)) + theme_bw()
ggarrange(p1, p2)

# define drm function to use with map
drm.func <- function(x) {
  drm(relative_gfp ~ temperature, 
      fct = LL.3(names = c("a", "b", "c")), 
      data = x)
}

predict.fun <- function(x) {
  add_predictions(data.frame(temperature = seq(30,50)), x)
}

coefs.fun <- function(x) {coef(x) %>% tidy}

df2 <- df %>% group_by(Mutant) %>% nest() %>%
  mutate(drmod = map(data, drm.func), 
         pred = map(drmod, predict.fun),
         coefs = map(drmod, coefs.fun))


a<-group_by(df, temperature,Mutant) %>%
  summarise(
    median = median(relative_gfp, na.rm = TRUE),
    IQR = IQR(relative_gfp, na.rm = TRUE),
    mean= mean(relative_gfp, na.rm = TRUE),
    c=n()
  )

dots <-   ggplot(a, aes(x=temperature, y=mean,color=Mutant,fill=Mutant)) + 
  geom_point()

# plot raw data, model and ED50 line
df2 %>% unnest(data) %>% 
  ggplot()  +
  geom_line(aes(temperature, pred, color = Mutant,size=Mutant),size=1, data =
              df2 %>% unnest(pred)) +
  geom_point(aes(temperature, mean, color = Mutant,shape=Mutant,size=Mutant),size=2, data =
              a)+
  geom_vline(aes(xintercept = x, color = Mutant), 
             linetype = 5,
             data = df2 %>% unnest(coefs) %>% filter(names == "ED50:(Intercept)")) +
  theme_minimal()+
  scale_color_manual(values=c("#f9c22e","#e76f51","#40916c","#b00286"))+
  scale_shape_manual(values=c(16, 2, 5,4))+ labs(title="",x="Temperature",y="Relative GFP")+ theme(legend.position = c(0.8, 0.9)) + labs(fill = "Mutant")
