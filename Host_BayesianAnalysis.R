#---
#R script to analyse parasite and temperature effects on host reproduction  
#authors: Andrew Jackson & Charlotte Kunze 
#---

#load packages
library(Hmisc)
library(cowplot)
library(tidyverse)
library(ggplot2)
library(scales)
library(rjags)
library(ggpubr)


#set the working directory in the same folder, our Rscript is saved 
setwd(dirname(rstudioapi::getSourceEditorContext()$path))


#### import host data ####
host_repro_data <- read.csv2("sum_daphnia_noNA.csv", header = TRUE, 
                             stringsAsFactors = FALSE) %>%
  mutate(treat2 = treat)

#check import:
str(host_repro_data)
names(host_repro_data)

#first plot
ggplot(host_repro_data, aes(x = mean_temp, y = sum, group = treat2, color = treat2)) + 
  geom_point() + 
  geom_smooth()+
  theme_bw()

#### Model formulation ####
#Source the beta function we used to model the TPCs which is held in an external `*.R` file as it is called from multiple notebooks.
source("betaFunction.R")


# Temperature T was replaced in the maths above with x when specifying the model.
model_beta_fun_poisson = '
model
{
  # Likelihood
  for (i in 1:n) {
  
  y[i] ~ dpois(lambda[i])
  
  # log10 link function
  lambda[i] <- 10^(mu[i])
  
  mu[i] <-      c_m[treat[i]] *
  ( (x_max[treat[i]] - x[i])  / (x_max[treat[i]] - x_opt[treat[i]]) ) *
  ( (x[i] - x_min[treat[i]])  / (x_opt[treat[i]] - x_min[treat[i]]) ) ^
  ( (x_opt[treat[i]] - x_min[treat[i]]) / (x_max[treat[i]] - x_opt[treat[i]]) )
  }
  
  # Priors
  for (j in 1:M) {
  
  x_min[j] ~ dunif(0, 14)
  x_opt[j] ~ dunif(x_min[j] + 0, x_min[j] + 25)
  x_max[j] ~ dunif(x_opt[j] + 0, 40)
  
  
  # these rather narrow priors worked during model testing
  # x_min[j] ~ dunif(5, 14)
  # x_opt[j] ~ dunif(15, 22)
  # x_max[j] ~ dunif(23, 26)
  
  
  
  c_m[j] ~ dunif(0, 10)
  }
  
  
  sigma ~ dunif(0, 100)
}
'

model_beta_cnx <- textConnection(model_beta_fun_poisson)

# prepare the new two-factor treatment column by combining usinig : notation 
# to denote interactions between treat and inf
host_repro_data <- host_repro_data %>% mutate(treat2 = factor(treat):factor(inf))

model_data_repro <- list(y = host_repro_data$sum, 
                         x = host_repro_data$mean_temp, 
                         n = nrow(host_repro_data),
                         treat = host_repro_data$treat2,
                         treat_levels = levels(host_repro_data$treat2),
                         M = length(levels(host_repro_data$treat2)))


# The jags model is initialised using jags.model()
model_repro <- jags.model(model_beta_cnx, data = model_data_repro, 
                          n.chains = 3, quiet = TRUE)

close(model_beta_cnx) # close the model connection as we dont need it any more


# We then use coda.samples to ask for posterior draws
# which we will use as a reflectin of the posterior.
output_repro <- coda.samples(model = model_repro, 
                             variable.names = c("c_m", "x_max", "x_min", "x_opt"), 
                             n.iter = 10000, thin = 5, quiet = TRUE)

gelman.diag(output_repro)

#summarise the model
post_smry_repro <- summary(output_repro)
print(post_smry_repro)
post_means_repro <- post_smry_repro$statistics[,1]
print(post_smry_repro)

#save model output in a text file
#sink("host_stats.csv")
print(summary(output_repro))
#sink()

#save model information in a dataframe
summary_model_repro<- summary(output_repro)

#sub data frame with only quantiles to add to the plot
repro_CI <- as.data.frame(summary_model_repro$quantiles) %>%
  tibble::rownames_to_column() %>%
  rename(dummy = rowname) %>%
  mutate(treatment =  ifelse(str_detect(dummy, '1'), paste('CS:I'),ifelse(str_detect(dummy, '2'), paste('CS:U'),ifelse(str_detect(dummy, '3'), paste('FLU:I'),ifelse(str_detect(dummy, '4'), paste('FLU:U'), ifelse(str_detect(dummy, '5'), paste('PULSE:I'),'PULSE:U'))) )))%>%
  gather(key = 'Quantile', value = 'value', -dummy, -treatment) %>%
  mutate(transform= ifelse(str_detect(dummy, 'c_m'), 10^value, value)) %>% #backtransform maximum spore data 
  select(-value) %>% #remove raw values to spread data
  spread(key = Quantile, value = transform)


# specify a range of temps over which to evaluate our function
# add in the x_min values + a very small amount to ensure the 
# lines are drawn at the lowest treshold. The resultant
# vector is sorted sequentially.

# find the indices in post_means_repro that match "x_min" and "x_max"
x_min_idx <- grep("x_min", names(post_means_repro))
x_max_idx <- grep("x_max", names(post_means_repro))
x_opt_idx <- grep("x_opt", names(post_means_repro))
c_m_idx   <- grep("c_m"  , names(post_means_repro))

# create the vector
xx <- seq(min(post_means_repro[x_min_idx]), 
          max(post_means_repro[x_max_idx]), length = 10^3)

# and the coresponding estimate based on the means of the posteriors
yy <- matrix(0, nrow = length(xx), ncol = model_data_repro$M)

for (i in 1:model_data_repro$M) {
  yy[,i] <- b_fun(xx, 
                  c_m =   post_means_repro[c_m_idx[i]], 
                  x_max = post_means_repro[x_max_idx[i]],
                  x_min = post_means_repro[x_min_idx[i]], 
                  x_opt = post_means_repro[x_opt_idx[i]])
}

#transform the list into a data frame
host_data <- host_repro_data %>%
  mutate(id = paste(treat2)) %>% #new column 
  separate(treat2, c('treat', 'inf'), ':') %>% #separate treat2 column into two new one
  group_by(treat, inf, mean_temp) %>% #group our data frame after treatment, inf and mean temp
  mutate(mean = mean(sum, na.rm = T), #calculate mean, sd, se for every treat, inf, mean_temp combination
         sd = sd(sum, na.rm =T),
         se = sd/sqrt(n())) %>%
  distinct(treat, inf, mean_temp,id, mean, se ) #remove duplicates, keeps only unique entries for those columns

pred1 <- as.data.frame(cbind(yy, xx) ) #new dataframe containing predictions and temperature (named xx)
names(pred1) <-c('CS:I', 'CS:U', 'FLU:I', 'FLU:U', 'PULSE:I', 'PULSE:U','xx') #name the columns 

pred <- pred1 %>% #dataframe consisting of our predictions
  pivot_longer(-xx, names_to = "treat2", values_to = "sum") %>% #melt my prediction columns together
  separate(treat2, c('treat', 'inf'), ':') %>% #separate column treat2 into two columns named treat and inf
  mutate(trans_pred = 10^sum) #create new column named trans_pred with transformed predictions.


#### plots for MS ####

#### Figure 3 ####
label_treat <- c(CS = ' ', PULSE = '  ', FLU = ' ')
data1 <- filter(host_data, treat == 'CS')

pred_1 <- filter(pred, treat == 'CS')
baby1 <- ggplot(data1, aes(y=mean, x=mean_temp, col = treat))+
  geom_line(data = pred_1, aes(x = xx, y = trans_pred, col = treat, linetype = inf), size = 0.8, alpha = 0.9)+
  geom_point(aes(shape = id),size = 3.5)+
  geom_errorbar(aes(ymin = mean-se, ymax =  mean+se), width =0.8)+
  #t opt
  geom_segment(x = repro_CI[19,3] , xend = repro_CI[19,7], y = 132,  yend = 132,  col = '#003C67FF', linetype = 6, size = 1.3)+ #constant
  geom_segment(x = repro_CI[20,3] , xend = repro_CI[20,7], y = 129.5,  yend = 129.5,  col = '#4A6990FF', size = 1.3)+ #un constant
  scale_y_continuous(limits = c(0,120), breaks = c(0,40, 80,120))+
  scale_x_continuous(breaks = seq(10, 28,3), limits=c(8,30))+
  scale_shape_manual(values=c( 16,1,17,2,15,0))+
  scale_color_manual(values=c('#003C67FF'),labels = c('constant'))+
  scale_linetype_manual(values=c(6,1), labels = c('exposed', 'unexposed'))+
  coord_cartesian(ylim = c(0, 120), clip="off")+ #changes coord system to include the CI arrows in our plot
  labs(x = '', y ='reproductive output', shape = ' ', linetype = 'Exposure',color = 'Treatment')+
 facet_wrap(~treat, labeller = labeller(treat = label_treat))+
  theme( panel.background = element_rect(fill = NA), #removes background
         #panel.grid.major.y = element_line(color='grey', linetype = 'dashed', size=0.2),
         panel.border= element_rect(colour = "black", fill=NA, size=0.5),
         strip.background = element_rect(fill = NA),
         legend.background = element_blank(),
         legend.title = element_text(hjust=0), #moves text to the left
         legend.key = element_blank(),
         legend.position = 'none',
         text = element_text(size=17),
         plot.margin = unit(c(1,0.8,1.5,0), "cm")) 
ggpar(baby1)

##ggplot of offspring in fluctuating treatment 
data2 <- filter(host_data, treat == 'FLU') #- subset
pred_2 <- filter(pred, treat == 'FLU')
baby2 <- ggplot(data2, aes(y=mean, x=mean_temp, col = treat))+
  geom_line(data = pred_2, aes(x = xx, y = trans_pred, col = treat, linetype = inf), size = 0.8, alpha = 0.9)+
  geom_point(aes(shape = id),size = 3.5)+
  geom_errorbar(aes(ymin = mean-se, ymax =  mean+se), width =0.8)+
  #
  #t opt
  geom_segment(x = repro_CI[21,3] , xend = repro_CI[21,7], y = 132,  yend = 132, col = '#EFC000FF', linetype = 6, size = 1.3)+ #flux
  geom_segment(x = repro_CI[22,3] , xend = repro_CI[22,7], y = 129.5,  yend = 129.5,  col = '#c7b514', size = 1.3)+ #un flux
 # scale adjusts
scale_y_continuous(breaks=NULL)+
  scale_x_continuous(breaks = seq(10, 28,3), limits=c(8,30))+
  scale_shape_manual(values=c( 17,2))+
  scale_color_manual(values=c('#EFC000FF','#A73030FF'),labels = c( 'fluctuation', 'heat wave'))+
  scale_linetype_manual(values=c(6,1), labels = c('exposed', 'unexposed'))+
  coord_cartesian(ylim = c(0, 120), clip="off")+ #changes coord system to include the CI arrows in our plot
  labs(x = 'mean Temperature (in °C)', y ='  ', shape = ' ', linetype = 'Exposure',color = 'Treatment')+
  facet_wrap(~treat, labeller = labeller(treat = label_treat))+
  theme( panel.background = element_rect(fill = NA), 
         #panel.grid.major.y = element_line(color='grey', linetype = 'dashed', size=0.2),
         panel.border= element_rect(colour = "black", fill=NA, size=0.5),
         strip.background = element_rect(fill = NA),
         legend.background = element_blank(),
         legend.title = element_text(hjust=0), 
         legend.key = element_blank(),
         legend.position = 'none',
         text = element_text(size=17),
         plot.margin = unit(c(1,0.8,1.5,0), "cm")) 

ggpar(baby2)
#ggsave(plot = last_plot(), file = 'final_babies_per_treatment.tiff', width = 14, height = 10)

#pulse baby
data3 <- filter(host_data, treat == 'PULSE')
pred_3 <- filter(pred, treat == 'PULSE')

baby3 <- ggplot(data3, aes(y=mean, x=mean_temp, col = treat))+
  geom_line(data = pred_3, aes(x = xx, y = trans_pred, col = treat, linetype = inf), size = 0.8, alpha = 0.9)+
  geom_point(aes(shape = inf),size = 3.5)+
  geom_errorbar(aes(ymin = mean-se, ymax =  mean+se), width =0.8)+
 #C max
  geom_segment(x = 31.7, xend = 31.7, y = repro_CI[1,3], yend = repro_CI[1,7], col = '#003C67FF' , size = 1.3)+#constant # inf
  geom_segment(x = 32.5, xend = 32.5, y = repro_CI[3,3],yend =repro_CI[3,7], col ='#EFC000FF', size = 1.3)+#flux #inf
  geom_segment(x = 33, xend = 33, y =  repro_CI[5,3],yend = repro_CI[5,7],col = '#A73030FF', size = 1.3) +#pulse #inf
  geom_segment(x = 31.9 , xend = 31.9, y = repro_CI[2,3],  yend = repro_CI[2,7], col ='#4A6990FF', size = 1.3)+ #constant un
  geom_segment(x = 32.7, xend = 32.7, y = repro_CI[4,3],  yend = repro_CI[4,7], col ='#EFC000FF',size = 1.3)+ #flux un
  geom_segment(x = 33  , xend = 33, y = repro_CI[6,3],  yend = repro_CI[6,7], col = '#CD534CFF' , size = 1.3)+ #pulse un
  #t opt
  geom_segment(x = repro_CI[23,3] , xend = repro_CI[23,7], y = 129.5,  yend = 129.5, col = '#A73030FF',linetype = 6, size = 1.3)+ #pulse
  geom_segment(x = repro_CI[24,3] , xend = repro_CI[24,7], y = 129.5,  yend = 129.5, col = '#CD534CFF', size = 1.3)+ #un pulse
 # scale adjusts
scale_y_continuous(breaks=NULL)+  
  scale_x_continuous(breaks = seq(10, 28,3), limits=c(8,30))+
  scale_shape_manual(values=c(15,0))+
  scale_color_manual(values=c('#A73030FF'),labels = c( 'heat wave'))+
  scale_linetype_manual(values=c(6,1), labels = c('exposed', 'unexposed'))+
  coord_cartesian(ylim = c(0, 120), clip="off")+ #changes coord system to include the CI arrows in our plot
  labs(x = ' ', y =' ', shape = ' ', linetype = 'Exposure',color = 'Treatment')+
  facet_wrap(~treat, labeller = labeller(treat = label_treat))+
  theme( panel.background = element_rect(fill = NA), 
         #panel.grid.major.y = element_line(color='grey', linetype = 'dashed', size=0.2),
         panel.border= element_rect(colour = "black", fill=NA, size=0.5),
         strip.background = element_rect(fill = NA),
         legend.background = element_blank(),
         legend.title = element_text(hjust=0),
         legend.key = element_blank(),
         legend.position = 'none',
         text = element_text(size=17),
         plot.margin = unit(c(1,0.8,1.5,0), "cm")) 

ggpar(baby3)
library(cowplot)
allBAB<-plot_grid(baby1,baby2,baby3,nrow = 1)
ggsave(allBAB, file = 'allBABIES.tiff', width = 12, height = 7)

#### Figure 2c ####
#### only infected animals ###
inf_pred <-pred %>%
  filter(inf == 'I') 
inf_host_data <- host_data%>%
  filter(inf == 'I') 

plot <- ggplot(inf_host_data, aes(y=mean, x=mean_temp, col = treat, shape = treat))+
  geom_line(data = inf_pred, aes(x = xx, y = trans_pred, col = treat), alpha = 0.9,linetype =  'longdash', size = 0.8)+
  geom_point(size = 3.5)+ #change shape of my points 
  geom_errorbar(aes(ymin = mean-se, ymax =  mean+se), width = .2)+
  #C max
  geom_segment(x = 31.8 , xend = 31.8, y = repro_CI[1,3], yend = repro_CI[1,7],  col = '#003C67FF', size = 1.3)+ #constant
  geom_segment(x = 32.2, xend = 32.2, y = repro_CI[3,3],yend =repro_CI[3,7], col ='#EFC000FF', size = 1.3)+ #flux
  geom_segment(x = 32.6  , xend = 32.6, y = repro_CI[5,3],yend = repro_CI[5,7],  col = '#A73030FF',size = 1.3)+ #pulse
  #t opt
  geom_segment(x = repro_CI[21,3] , xend = repro_CI[22,7], y = 135,  yend = 135,  col = '#003C67FF', size = 1.3)+ #constant
  geom_segment(x = repro_CI[22,3] , xend = repro_CI[22,7], y = 133,  yend = 133,  col = '#EFC000FF', size = 1.3)+ #flux
  geom_segment(x = repro_CI[23,3] , xend = repro_CI[23,7], y = 131,  yend = 131,  col = '#A73030FF',size = 1.3)+ #pulse
   scale_x_continuous(breaks = seq(10, 28,3), limits=c(8,30))+
  scale_y_continuous(limits = c(0,120), breaks = seq(0,120, 30))+
  scale_shape_manual(values=c( 16,17,15), labels = c('constant', 'fluctuating', 'heat wave'))+
  scale_color_manual(values=c('#003C67FF','#EFC000FF','#A73030FF'),labels = c('constant', 'fluctuating', 'heat wave'))+
  labs(title = ' ',x = 'mean Temperature (in °C)', y = 'reproductive output', col = ' ', shape = ' ')+
  coord_cartesian(ylim = c(0, 120), clip="off")+ #changes coord system to include the CI arrows in our plot
  theme( panel.background = element_rect(fill = NA), 
         #panel.grid.major.y = element_line(color='grey', linetype = 'dashed', size=0.2),
         panel.border= element_rect(colour = "black", fill=NA, size=0.5),
         strip.background = element_rect(color = 'black', fill = 'grey95'),
         legend.background = element_blank(),
         legend.title = element_text(hjust=0), 
         legend.position  =c(0.85,0.94),
         legend.key = element_blank(),
         text = element_text(size=21),
         plot.margin = unit(c(1,2,1.5,1.5), "cm")) #t, r, b, l bzw. top, right, bottom, left
ggpar(plot)



