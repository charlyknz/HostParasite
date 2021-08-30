#---
#R script to analyse spore data and infectivity of Odospora colligata 
#authors: Andrew Jackson & Charlotte Kunze 
#---

#load packages
library(Hmisc)
library(tidyverse)
library(ggplot2)
library(scales)
library(rjags)
library(ggpubr)
library(cowplot)

#set the working directory in the same folder, our Rscript is saved 
setwd(dirname(rstudioapi::getSourceEditorContext()$path))


#### Import the data ####
spore_data <- read.csv("~/Desktop/Irland/R code/Baysian/SporesNoMaleNA.csv", 
                       header = TRUE, stringsAsFactors = FALSE)
# only exposed animals
spore_data_I <- spore_data %>% filter(exposed == "I")
names(spore_data_I)


#first-plots
g1 <- ggplot(spore_data_I, aes(x = realtemp, y = no_spore, 
                               color = treatment, 
                               group = treatment)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_point() + 
  geom_smooth()

print(g1)



## Try fitting the more complex beta function to the infection rate data
#These data are the column `infect` in the subsetted data object `spore_data_I`.
#In the jags code below i need to force T_max and T_min to be outside the range of the data, else it generates invalid numbers and wont run. 
                                                                                                                                                                                                                                               
#### Model formulation infectivity ####

#Source the beta function we used to model the TPCs which is held in an external `*.R` file as it is called from multiple notebooks.
source("betaFunction.R") 

# add points for the average by temperature
summary_spore <- spore_data_I  %>%  
  group_by(treatment, temperature) %>% 
  summarise(mu = mean(infect), tt = mean(realtemp), 
            Ninfect = sum(infect), Ntrials = length(infect))

#Fitting binominal model, aggregating the data by temperature and treatment.

# Here Temperature T was replaced in the maths above with x when specifying the model.
model_beta_fun_binomial = '
model
{
  # Likelihood
  for (i in 1:n) {
  
  y[i] ~ dbinom(p[i], N[i])
  
  # make sure all values are >= 0
  p[i] <- min(max(mu[i], 0) , 1)
  
  
  mu[i] <-      c_m[treat[i]] *
  ( (x_max[treat[i]] - x[i])  / (x_max[treat[i]] - x_opt[treat[i]]) ) *
  ( (x[i] - x_min[treat[i]])  / (x_opt[treat[i]] - x_min[treat[i]]) ) ^
  ( (x_opt[treat[i]] - x_min[treat[i]]) / (x_max[treat[i]] - x_opt[treat[i]]) )
  }
  
  # Priors
  for (j in 1:M) {
  x_min[j] ~ dunif(5, 15)
  x_opt[j] ~ dunif(x_min[j] + 0, x_min[j] + 20)
  x_max[j] ~ dunif(x_opt[j] + 0, 35)
  
  c_m[j] ~ dunif(0, 1)
  }
  
  
}
'


#Set up and fit the binomial model to the `spore_data_I` object.

model_beta_cnx <- textConnection(model_beta_fun_binomial)

model_data_binomial <- list(y = summary_spore$Ninfect, 
                            x = summary_spore$tt,
                            N = summary_spore$Ntrials,
                            n = nrow(summary_spore), 
                            treat = as.numeric(factor(summary_spore$treatment)),
                            M = length(unique(summary_spore$treatment)))

# The jags model is initialised using jags.model()
model_binomial <- jags.model(model_beta_cnx, data = model_data_binomial,
                             n.chains = 3, quiet = TRUE)

close(model_beta_cnx) # close the model connection as we dont need it any more


# We then use coda.samples to ask for posterior draws
# which we will use as a reflectin of the posterior.
output_binomial <- coda.samples(model = model_binomial, 
                                variable.names = c("c_m", "x_max", "x_min", "x_opt"), 
                                n.iter = 1000, quiet = TRUE)

#BGR test
gelman.diag(output_binomial)

#Summarise the posterior
post_smry_binomial <- summary(output_binomial)
post_means_binomial <- post_smry_binomial$statistics[,1]
print(post_smry_binomial)

#save model output in a text file using sink()
#sink('infectivity.csv')
print(summary(output_binomial))

#save model information in a dataframe
summary_model<- summary(output_binomial)

#sub data frame with only quantiles to add to the plot
summy_model <- as.data.frame(summary_model$quantiles) %>%
  tibble::rownames_to_column() %>%
  rename(dummy = rowname) %>%
  mutate(treatment =  ifelse(str_detect(dummy, '1'), paste('constant'),ifelse(str_detect(dummy, '2'), paste('fluctuation'), 'pulse')))


# specify a range of temps over which to evaluate our function
xx <- seq(10, 30, length = 100)

# and the coresponding estimate based on the means of the posteriors
yy1 <- b_fun(xx, 
             c_m = post_means_binomial["c_m[1]"], 
             x_max = post_means_binomial["x_max[1]"],
             x_min = post_means_binomial["x_min[1]"], 
             x_opt = post_means_binomial["x_opt[1]"])

yy2 <- b_fun(xx, 
             c_m = post_means_binomial["c_m[2]"], 
             x_max = post_means_binomial["x_max[2]"],
             x_min = post_means_binomial["x_min[2]"], 
             x_opt = post_means_binomial["x_opt[2]"])

yy3 <- b_fun(xx, 
             c_m = post_means_binomial["c_m[3]"], 
             x_max = post_means_binomial["x_max[3]"],
             x_min = post_means_binomial["x_min[3]"], 
             x_opt = post_means_binomial["x_opt[3]"])


# write the predictions in a data frame for further use 
inf_cs <- as.data.frame(cbind(xx,yy1)) %>%
  mutate(treatment = paste('constant'),
         prediction = yy1) %>%
  select(-yy1)
inf_flux <- as.data.frame(cbind(xx,yy2))%>%
  mutate(treatment = paste('fluctuation'),
         prediction = yy2)%>%
  select(-yy2)
inf_pulse <- as.data.frame(cbind(xx,yy3))%>%
  mutate(treatment = paste('pulse'),
         prediction = yy3)%>%
  select(-yy3)

#merge predictions and data for single treatments into one df
inf_total <- rbind(inf_cs, inf_flux, inf_pulse)


# save the predicitions in a dataframe to make the values accessable 
df <- as.data.frame(model_data_binomial) %>%
  mutate(treatment = paste(ifelse(treat == 1, 'constant', ifelse(treat == 2, 'fluctuation', 'pulse')))) 

#calculate the CI of the model for the different temperatures
data = binconf(x=df$y, n=df$N, alpha = 0.05) #store in a list
data = as.data.frame(data) #change list to data.frame

#import data with the mean temperature and treatment to bind to the orginial dataset
MeanTemp <- read_csv("MeanTempTreatments.csv") %>%
  rename( x = meanTemp) %>%
  select(-X1)

#create a new total df with temperature information
all_df <- cbind(data, MeanTemp) %>% 
  left_join(df, by = c('x', 'treatment')) %>% #merge with the old one by temperature and treatment
  group_by(x, treatment) %>%
  mutate(sd = (sqrt(n())* (Upper-Lower)/3.92), #calculate se from lower and upper ci
         se = sd/sqrt(n()))


#### infectivity plot Fig. 2a####

# Note: 2.5 % and 97.5 % CI of maximum infectivity, Tmin, Topt and Tmax are added to the plot using geom_segment
# data are stored in the summy_model datafile
inf_plot <- ggplot(all_df, aes(y = I(model_data_binomial$y / model_data_binomial$N), x = model_data_binomial$x, col = treatment, shape = treatment))+
  geom_point(size = 3.5)+
  geom_errorbar(aes(ymin = I(model_data_binomial$y / model_data_binomial$N) -se, ymax = I(model_data_binomial$y / model_data_binomial$N)+se), width = .2)+
   #cm
  geom_segment(x = 32, xend = 32, y =  summy_model[1,2], yend = summy_model[1,6], col ='#003C67FF',size = 1.3)+ #constant
  geom_segment(x = 32.5, xend = 32.5, y = summy_model[2,2],yend =summy_model[2,6], col ='#EFC000FF', size = 1.3)+ #fluctuation
  geom_segment(x = 33, xend = 33, y =  summy_model[3,2],yend = summy_model[3,6],col = '#A73030FF',size = 1.3) +#pulse
  #Topt
  geom_segment(x = summy_model[10,2] , xend = summy_model[10,6], y = 1.29,  yend = 1.29, col = '#003C67FF',size = 1.3)+ #constant
  geom_segment(x = summy_model[11,2], xend = summy_model[11,6], y = 1.27,  yend = 1.27, col = '#EFC000FF',size = 1.3)+ #fluctuation
  geom_segment(x = summy_model[12,2]  , xend = summy_model[12,6], y = 1.25,  yend = 1.25,col = '#A73030FF',size = 1.3)+ #pulse
  # Tmin
  # for better visibility we included only the 97.5% interval, the 2.5% interval below 8 is indicated by arrows
  geom_segment(x = 8 , xend = summy_model[7,6] , y = 1.29,  yend = 1.29, arrow = arrow(ends = 'first', length = unit(0.11,'cm')), col = '#003C67FF',size = 1.3)+#constant
  geom_segment(x = 8, xend = summy_model[8,6], y = 1.27,  yend = 1.27, arrow = arrow(ends = 'first', length = unit(0.11,'cm')), col = '#EFC000FF',size = 1.3)+ #fluctuation
  geom_segment(x = 8, xend =  summy_model[9,6] , y = 1.25,  yend = 1.25, arrow = arrow(ends = 'first',length = unit(0.11,'cm')), col = '#A73030FF',size = 1.3)+#pulse
  #Tmax
  # for better visibility Tmax values over 30 degree are indicated by arrows
  geom_segment(x = summy_model[4,2] , xend =  30, y = 1.29,  yend = 1.29, arrow = arrow(length = unit(0.11,'cm')), col = '#003C67FF',size = 1.3)+#constant
  geom_segment(x = summy_model[5,2] , xend = summy_model[5,6], y = 1.27,  yend = 1.27, col = '#EFC000FF',size = 1.3)+#fluctuation
  geom_segment(x = summy_model[6,2], xend = 30, y = 1.25, yend = 1.25, arrow = arrow(length = unit(0.11,'cm')), col = '#A73030FF',size = 1.3)+#pulse
  geom_line(data = inf_total, aes(x = xx, y = prediction, color = treatment), linetype = 'longdash', size = 0.8)+
  scale_x_continuous(breaks = seq(10, 28,3), limits=c(8,30))+
  scale_y_continuous(limits = c(0,1.15), breaks= seq(0,1,0.2))+
  scale_shape_manual(values=c( 16,17,15), labels = c('constant', 'fluctuating', 'heat wave'))+
  scale_color_manual(values=c('#003C67FF','#EFC000FF','#A73030FF'),labels = c('constant', 'fluctuating', 'heat wave'))+
  labs(title = '',x = 'mean Temperature (in °C)', y = 'infection rate', col = ' ', shape = ' ')+
  coord_cartesian(ylim = c(0, 1.15), clip="off")+ #changes coord system to include the CI arrows in our plot
  theme( panel.background = element_rect(fill = NA), 
         panel.border= element_rect(colour = "black", fill=NA, size=0.5),
         strip.background = element_rect(color = 'black', fill = 'grey95'),
         legend.background = element_blank(),
         legend.title = element_text(hjust=0), 
         legend.position  = 'none',
         legend.key = element_blank(),
         text = element_text(size=21),
         plot.margin = unit(c(1,2,1.5,1.5), "cm")) #t, r, b, l bzw. top, right, bottom, left
ggpar(inf_plot)

#####################################################################################################################

#### Burden Data ####

#Now we can explore and model the data for parasite burden. 
#These are individuals that have been infected *and* have positive burden *and* died within the last X days of the experiment.

burden_data <- spore_data_I %>% filter(no_spore > 0 & lastday == 1)
par <- group_by(burden_data,  treatment) %>%
  count()
#Plot the burden data

g3 <- ggplot(burden_data, aes(x = realtemp, y = no_spore, 
                              color = treatment, 
                              group = treatment)) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_point() + 
  geom_smooth()

print(g3)
 
#### Model formulation burden ####

#Define a JAGS model for gaussian errors with a beta function - useful for the burden data. 

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
  x_min[j] ~ dunif(5, 14)
  # x_opt[j] ~ dunif(15, 22)
  # x_max[j] ~ dunif(23, 26)
  
  x_opt[j] ~ dunif(x_min[j] + 0, x_min[j] + 10)
  x_max[j] ~ dunif(x_opt[j] + 0, 30)
  
  
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

model_data_burden <- list(y = burden_data$no_spore, 
                          x = burden_data$realtemp, 
                          n = nrow(burden_data),
                          treat = as.numeric(factor(burden_data$treatment)),
                          M = length(unique(burden_data$treatment)))

# The jags model is initialised using jags.model()
model_burden <- jags.model(model_beta_cnx, data = model_data_burden, 
                           n.chains = 3, quiet = TRUE)

close(model_beta_cnx) # close the model connection as we dont need it any more


# We then use coda.samples to ask for posterior draws
# which we will use as a reflectin of the posterior.
output_burden <- coda.samples(model = model_burden, 
                              variable.names = c("c_m", "x_max", "x_min", "x_opt"), 
                              n.iter = 1000, quiet = TRUE)


gelman.diag(output_burden) #The BGR test should return values less than 1.1 to indicate convergence.

#summarise the burden model#
post_smry_burden <- summary(output_burden)
print(post_smry_burden)

post_means_burden <- post_smry_burden$statistics[,1]


print(post_smry_burden)

#save stats in a txt file
#sink('parasite_burden.csv')
print(summary(output_burden))
#sink()

# note: C_m has to be transponated with 10^


#save model information in a dataframe
summary_model_burden<- summary(output_burden)

#sub data frame with only quantiles to add to the plot
burden_CI <- as.data.frame(summary_model_burden$quantiles) %>%
  tibble::rownames_to_column() %>%
  rename(dummy = rowname) %>%
  mutate(treatment =  ifelse(str_detect(dummy, '1'), paste('constant'),ifelse(str_detect(dummy, '2'), paste('fluctuation'), 'pulse'))) %>%
  gather(key = 'Quantile', value = 'value', -dummy, -treatment) %>%
  mutate(transform= ifelse(str_detect(dummy, 'c_m'), 10^value, value)) %>% #backtransform maximum spore data 
  select(-value) %>% #remove raw values to spread data
  spread(key = Quantile, value = transform)

# specify a range of temps over which to evaluate our function
# add in the x_min values + a very small amount to ensure the 
# lines are drawn just above the lowest treshold. The resultant
# vector is sorted sequentially.
xx <- sort(c(post_means_burden["x_min[1]"] + 10^-20,
             post_means_burden["x_min[2]"] + 10^-20, 
             post_means_burden["x_min[3]"] + 10^-20, 
             seq(10, 30, length = 10^3) ) )

# and the coresponding estimate based on the means of the posteriors
yy1 <- b_fun(xx, 
             c_m = post_means_burden["c_m[1]"], 
             x_max = post_means_burden["x_max[1]"],
             x_min = post_means_burden["x_min[1]"], 
             x_opt = post_means_burden["x_opt[1]"])

yy2 <- b_fun(xx, 
             c_m = post_means_burden["c_m[2]"], 
             x_max = post_means_burden["x_max[2]"],
             x_min = post_means_burden["x_min[2]"], 
             x_opt = post_means_burden["x_opt[2]"])

yy3 <- b_fun(xx, 
             c_m = post_means_burden["c_m[3]"], 
             x_max = post_means_burden["x_max[3]"],
             x_min = post_means_burden["x_min[3]"], 
             x_opt = post_means_burden["x_opt[3]"])

#store predictions in data frames
spore_cs <- as.data.frame(cbind(yy1,xx))%>%
  mutate(treatment = paste('constant'),
         prediction = yy1) %>%
  select(-yy1)
spore_flux <- as.data.frame(cbind(yy2,xx))%>%
  mutate(treatment = paste('fluctuation'),
         prediction = yy2) %>%
  select(-yy2)
spore_pulse <- as.data.frame(cbind(yy3,xx))%>%
  mutate(treatment = paste('pulse'), #add treatment specification
         prediction = yy3) %>%
  select(-yy3)

#merge the predictions
spore_total <- rbind(spore_cs, spore_flux, spore_pulse) %>%
  mutate(trans_pred = 10^prediction)

#calculate mean, sd,se
data <- as.data.frame(model_data_burden) %>%
  mutate(treatment = paste(ifelse(treat == 1, 'constant', ifelse(treat == 2, 'fluctuation', 'pulse')))) %>%
  group_by(treatment, x) %>%
  mutate(mean = mean(y, na.rm = T),
         sd = sd(y, na.rm = T),
         se = sd/sqrt(n()))

#### Plot burden Fig. 2b####

# Note: CI of maximum spore burden, topt, tmin,tmax are stored in the dataframe burden_CI
## and are added in geom_segment
spore_plot <- ggplot(data, aes(x = x, y = mean, color = treatment, shape = treatment ))+
  geom_line(data = spore_total, aes(x = xx, y = trans_pred, color = treatment), linetype = 'longdash', size = 0.8)+
  geom_point(size = 3.5)+
  geom_errorbar(aes(ymin = mean-se, ymax =  mean+se), width = .2)+
  #Cm
  geom_segment(x = 32, xend = 32, y = burden_CI[1,3], yend =burden_CI[1,7],col ='#EFC000FF', size = 1.3)+#cs
  geom_segment(x = 32.5, xend = 32.5, y =  burden_CI[2,3], yend = burden_CI[2,7], col ='#003C67FF', size = 1.3)+#flux
  geom_segment(x = 33, xend = 33, y =   burden_CI[3,3],yend =  burden_CI[3,7],col = '#A73030FF', size = 1.3) +#pulse
  #Topt
  geom_segment(x = burden_CI[10,3], xend = burden_CI[10,7], y = 950,  yend = 950,  col = '#003C67FF', size = 1.3)+ #constant
  geom_segment(x = burden_CI[11,3], xend = burden_CI[11,7], y = 938,  yend = 938, col = '#EFC000FF', size = 1.3)+ #flux
  geom_segment(x = burden_CI[12,3], xend = burden_CI[12,7], y = 930,  yend = 930, col = '#A73030FF', size = 1.3)+ #pulse
  # Tmin
  geom_segment(x = burden_CI[7,3] , xend = burden_CI[7,7] , y = 950,  yend = 950,  col = '#003C67FF', size = 1.3)+#cs
  geom_segment(x = burden_CI[8,3], xend = burden_CI[8,7], y = 938,  yend = 938,  col = '#EFC000FF', size = 1.3)+#fl
  geom_segment(x = burden_CI[9,3], xend =  burden_CI[9,7] , y = 930,  yend = 930,  col = '#A73030FF', size = 1.3)+#pl
  #Tmax
  geom_segment(x = burden_CI[4,3] , xend =  burden_CI[4,7], y = 950,  yend = 950,  col = '#003C67FF', size = 1.3)+#cs
  geom_segment(x = burden_CI[5,3] , xend = burden_CI[5,7], y = 938,  yend = 938,  col = '#EFC000FF', size = 1.3)+#fl
  geom_segment(x = burden_CI[6,3], xend = burden_CI[6,7], y = 930,  yend = 930, col = '#A73030FF', size = 1.3)+#pl
  scale_y_continuous(limits = c(0,850), breaks= seq(0,800,200))+
  scale_x_continuous(breaks = seq(10, 28,3), limits=c(8,30))+
  scale_shape_manual(values=c( 16,17,15), labels = c('constant', 'fluctuating', 'heat wave'))+
  scale_color_manual(values=c('#003C67FF','#EFC000FF','#A73030FF'),labels = c('constant', 'fluctuating', 'heat wave'))+
  labs(title = '',x = 'mean Temperature (in °C)', y = 'spore burden (no of clusters)', col = ' ', shape = ' ')+
  coord_cartesian(ylim = c(0, 850), clip="off")+ #changes coord system to include the CI arrows in our plot
  theme( panel.background = element_rect(fill = NA), 
         panel.border= element_rect(colour = "black", fill=NA, size=0.5),
         strip.background = element_rect(color = 'black', fill = 'grey95'),
         legend.background = element_blank(),
         legend.title = element_text(hjust=0),
         legend.position  ='none',
         legend.key = element_blank(),
         text = element_text(size=21),
         plot.margin = unit(c(1,2,1.5,1.5), "cm")) #t, r, b, l bzw. top, right, bottom, left
ggpar(spore_plot)


#### Figure 2 ####
#note: plot is the host data using the bayesian_host.R script
plot_grid(inf_plot,spore_plot,plot,nrow = 1)
ggsave(plot = last_plot(), 'host_and_parasite_traits_se_inf.tiff', width = 24, height = 9)


