# R Script to analyse HOBO data
## by Charlotte Kunze 01 May 2019

setwd("~/Desktop/Irland/data")
library(lubridate)
library(tidyverse)
library(scales)

#### import data ####
## constant HOBOs
constant_10 <- read_csv2('~/Desktop/Irland/data/hobo/10_degree_constant_end.csv', col_types = cols(date = col_datetime(format = "%m.%d.%y %I:%M:%S%p")) )%>%
  select(no, date, temp) %>%
  filter(temp != 'NA') 

constant_13 <- read_csv2('~/Desktop/Irland/data/hobo/13_degree_constant_end.csv', col_types = cols(date = col_datetime(format = "%m.%d.%y %I:%M:%S%p")) )%>%
  select(no, date, temp) %>%
  filter(temp != 'NA') 

constant_16 <- read_csv2('~/Desktop/Irland/data/hobo/16_degree_light_constant_end.csv', col_types = cols(date = col_datetime(format = "%m.%d.%y %I:%M:%S%p")) )%>%
  select(no, date, temp, light) %>%
  filter(temp != 'NA') 
constant_16_t <- select(constant_16, -light)
constant_19 <- read_csv2('~/Desktop/Irland/data/hobo/19_degree_constant_end.csv', col_types = cols(date = col_datetime(format = "%m.%d.%y %I:%M:%S%p")) )%>%
  select(no, date, temp) %>%
  filter(temp != 'NA') 

constant_22 <- read_csv2('~/Desktop/Irland/data/hobo/22_degree_constant_end.csv', col_types = cols(date = col_datetime(format = "%m.%d.%y %I:%M:%S%p")) )%>%
  select(no, date, temp)%>%
  filter(temp != 'NA') 

constant_25 <- read_csv2('~/Desktop/Irland/data/hobo/25_degree_constant_end.csv', col_types = cols(date = col_datetime(format = "%m.%d.%y %I:%M:%S%p")) )%>%
  select(no, date, temp)%>%
  filter(temp != 'NA') 

constant_28 <- read_csv2('~/Desktop/Irland/data/hobo/28_degree_constant_end.csv', col_types = cols(date = col_datetime(format = "%m.%d.%y %I:%M:%S%p")) )%>%
  select(no, date, temp)%>%
  filter(temp != 'NA') 



## FLUX HOBOs ### 
flux_mean13 <- read_csv2('~/Desktop/Irland/data/hobo/10-16_degree_flux_end.csv', col_types = cols(date = col_datetime(format = "%m.%d.%y %I:%M:%S%p")) )%>%
  select(no, date, temp)%>%
  filter(temp != 'NA') 
  
flux_mean16 <- read_csv2('~/Desktop/Irland/data/hobo/13-19_degree_flux_end.csv', col_types = cols(date = col_datetime(format = "%m.%d.%y %I:%M:%S%p")) )%>%
  select(no, date, temp)%>%
  filter(temp != 'NA') 

  
flux_mean19 <- read_csv2('~/Desktop/Irland/data/hobo/16-22_degree_flux_end.csv', col_types = cols(date = col_datetime(format = "%m.%d.%y %I:%M:%S%p")) )%>%
  select(no, date, temp)%>%
  filter(temp != 'NA') 


flux_mean22 <- read_csv2('~/Desktop/Irland/data/hobo/19-25_degree_flux_end.csv', col_types = cols(date = col_datetime(format = "%m.%d.%y %I:%M:%S%p")) )%>%
  select(no, date, temp)%>%
  filter(temp != 'NA') 

flux_mean25 <- read_csv2('~/Desktop/Irland/data/hobo/22-28_degree_flux_end.csv', col_types = cols(date = col_datetime(format = "%m.%d.%y %I:%M:%S%p")) )%>%
  select(no, date, temp)%>%
  filter(temp != 'NA') 


## PULSE HOBOs

pulse_13 <- read_csv2('~/Desktop/Irland/data/hobo/13_degree_pulse_+6_end.csv', col_types = cols(date = col_datetime(format = "%m.%d.%y %I:%M:%S%p")) )%>%
  select(no, date, temp)%>%
  filter(temp != 'NA') 

pulse_16 <- read_csv2('~/Desktop/Irland/data/hobo/16_degree_pulse_+6_end.csv', col_types = cols(date = col_datetime(format = "%m.%d.%y %I:%M:%S%p")) )%>%
  select(no, date, temp)%>%
  filter(temp != 'NA')

pulse_19 <- read_csv2('~/Desktop/Irland/data/hobo/19_degree_pulse_+6_end.csv', col_types = cols(date = col_datetime(format = "%m.%d.%y %I:%M:%S%p")) )%>%
  select(no, date, temp)%>%
  filter(temp != 'NA')

pulse_22 <- read_csv2('~/Desktop/Irland/data/hobo/22_degree_light_pulse_+6_end.csv', col_types = cols(date = col_datetime(format = "%m.%d.%y %I:%M:%S%p")) )%>%
  select(no, date, temp, light)%>% # light HOBO
  filter(temp != 'NA')
pulse_22_t <- select(pulse_22, -light)


### --------------------------------------------------------------------------------###
#### add tratment information 
constant_10 <- mutate(constant_10, id = paste('constant 10'),
                      expected_temp = paste(10)) 
constant_13 <- mutate(constant_13, id = paste('constant 13'),
                      expected_temp = paste(13))
constant_16_t <- mutate(constant_16_t, id = paste('constant 16'),
                        expected_temp = paste(16))
constant_19 <- mutate(constant_19, id = paste('constant 19'),
                      expected_temp = paste(19))
constant_22 <- mutate(constant_22, id = paste('constant 22'),
                      expected_temp = paste(22))
constant_25 <- mutate(constant_25, id = paste('constant 25'),
                      expected_temp = paste(25))
constant_28 <- mutate(constant_28, id = paste('constant 28'),
                      expected_temp = paste(28))
flux_mean13 <- mutate(flux_mean13, id = paste('flux 13'),
                      expected_temp = paste(13))
flux_mean16 <- mutate(flux_mean16, id = paste('flux 16'),
                      expected_temp = paste(16))
flux_mean19 <- mutate(flux_mean19, id = paste('flux 19'),
                      expected_temp = paste(19))
flux_mean22 <- mutate(flux_mean22, id = paste('flux 22'),
                      expected_temp = paste(22))
flux_mean25 <- mutate(flux_mean25, id = paste('flux 25'),
                      expected_temp = paste(25))
pulse_13 <- mutate(pulse_13, id = paste('pulse 13'),
                   expected_temp = paste(13))
pulse_16 <- mutate(pulse_16, id = paste('pulse 16'),
                   expected_temp = paste(16))
pulse_19 <- mutate(pulse_19, id = paste('pulse 19'),
                   expected_temp = paste(19))
pulse_22_t <- mutate(pulse_22_t, id = paste('pulse 22'),
                     expected_temp = paste(22))

# bind all data together
all_raw_temp <- rbind(constant_10, constant_13, constant_16_t, constant_19, constant_22, constant_25, constant_28,
                      flux_mean13,flux_mean16,flux_mean19, flux_mean22,flux_mean25, pulse_13, pulse_16,  pulse_19, pulse_22_t) 
#write.csv(all_raw_temp, file = 'all_raw_temperatures_hobo.csv')


### --------------------------------------------------------------------------------###
## --------------------------------------------------------------------------------###
#-----------------------------------------------------------------------------------#
#### calculate mean temperatures####
names(all_raw_temp)
mean_temp <- all_raw_temp %>%
  group_by(id) %>%
  summarise(realtemp = mean(temp)) %>%
  separate(id, into = c('treat', 'expected_temp'), sep = ' ')

# plot the expected temperature against observed temperatures 
ggplot(mean_temp, aes(x = as.numeric(expected_temp), y= realtemp, col = as.factor(treatment)))+
  geom_point(size = 3)+
  geom_abline(intercept = 0, slope = 1)+ #adds a vertical line
  labs( col = 'Treatment', x = 'expected mean temperature (in °C)', y = 'measured mean temperature (in °C)')+
  scale_color_manual(values=c('#003C67FF','#EFC000FF','#A73030FF'),labels = c('constant', 'fluctuating', 'heat wave'))+
  theme( panel.background = element_rect(fill = NA), #loescht den Hintergrund meines Plots/ fuellt ihn mit nichts
         panel.grid.major.y = element_line(color='grey', linetype = 'dashed', size=0.2),
         panel.border= element_rect(colour = "black", fill=NA, size=0.5),
         strip.background = element_rect(color = 'black', fill = 'grey95'),
         legend.background = element_blank(),
         legend.title = element_text(hjust=3), #schiebt text nach links in der Legende
         legend.position  ='bottom',
         legend.key = element_blank(),
         text = element_text(size = 13))
#ggsave('temperature_range.png',plot = last_plot(), width = 8, height = 6) #saves your plot at the given directory

mean_temp <- mean_temp %>%
  select(-expected_temp) %>%
  rename(meanTemp = realtemp, 
         treatment = treat) %>%
  mutate(treatment = paste(ifelse(treatment == 'flux', 'fluctuation', treatment)))
write.csv(MeanTempTreatments.csv)
names(mean_temp)
 
write.csv(mean_temp, 'MeanTempTreatments.csv')


#### temperature plots ####
temperature<-all_raw_temp  %>%
  filter(date > ymd_hms("2019-04-09 18:24:09"))%>%
  filter(date < ymd_hms("2019-05-06 09:24:09"))
#constant
ggplot(subset(temperature, id %in% c('constant 10', 'constant 13', 'constant 16', 'constant 19', 'constant 22', 'constant 25', 'constant 28' )), aes(x = date, y = temp))+
         geom_point()+
  facet_wrap(~id)+
  #scale_y_continuous(breaks = seq(10,28,6), limits = c(5, 30))+
  #scale_x_datetime(breaks = '7 days', date_labels = '%d.%m')+  
  labs(x = 'Date', y = 'Temperature (in °C)')+
  theme_classic()

#other
ggplot(subset(temperature, !id %in% c('pulse 22','constant 10', 'constant 13', 'constant 16', 'constant 19', 'constant 22', 'constant 25', 'constant 28' )), aes(x = date, y = temp))+
  geom_point(size = 0.5)+
  geom_line(size = 1)+
  facet_wrap(~id)+
 labs(x = 'Date', y = 'Temperature (in °C)')+
  theme_classic()
#ggsave(last_plot(), file = 'variable_temperatures.png',height = 8, width = 11)

pulse_22 %>%
  filter(date > ymd_hms("2019-04-09 18:24:09"))%>%
  filter(date < ymd_hms("2019-05-06 09:24:09"))%>%
  ggplot(., aes(x = date, y = temp))+
  geom_point()+
  geom_line()+
  geom_hline(aes(yintercept = mean(temp)), colour="blue")+
  scale_x_datetime(breaks = '3 days', date_labels = "%d.%m")+
  scale_y_continuous(breaks = seq(21, 29, 2), limits = c(20, 29))+
  labs(x = ' ', y = 'Temperature (in °C)')+
  theme_bw()


#light plot
pulse_22 %>%
  filter(date > ymd_hms("2019-04-09 18:24:09"))%>%
  filter(date < ymd_hms("2019-05-06 09:24:09"))%>%
  ggplot(., aes(x = date, y = light))+
  geom_point()+
  geom_line()+
  geom_hline(aes(yintercept = mean(light)), colour="blue")+
  scale_x_datetime(breaks = '3 days', date_labels = "%d.%m")+
  scale_y_continuous(breaks = seq(0, 3500, 500), limits=c(0, 3500))+
  labs(x = ' ', y = 'Light intensity (in lux)')+
  theme_bw()
#ggsave('pulse_22_light.png',plot = last_plot(), width = 8, height = 6) #saves your plot at the given directory

