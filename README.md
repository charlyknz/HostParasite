# HostParasite


README: Alternate patterns of temperature variation bring about very different disease outcomes at different mean temperatures


list of authors: Charlotte Kunze 1,2*, Pepijn Luijckx 2*#, Andrew L. Jackson 2, Ian Donohue 2

1 Institute for Chemistry and Biology of the Marine Environment [ICBM], Carl-von-Ossietzky, University Oldenburg
2 Department of Zoology, School of Natural Sciences, Trinity College Dublin, Dublin, Ireland


These datasets contain observations on host fecundity, infection status and severity of infection collected during an experiment on Daphnia magna and its parasite Odospora colligata at School of Natural Sciences at the Trinity College Dublin in spring 2019. 
In this study we compared two patterns of temperature variation, a diurnal temperature fluctuation (of ±3°C in a 12:12 cycle) and a heat wave (of +6°C with 3 days duration 20 days after exposure to the parasite), with a constant treatment, which had the same mean temperature over the whole experimental period. 
Data were analyzed using R version 3.6.1(R Core Team, 2018) interfacing with JAGS. A Beta Function was fitted to each of our different fitness estimates (that is, host fecundity, parasite infectivity and burden) for each of the three temperature regimes using a baysian framework.



Datasets and R scripts used for the analyses: 

List of RScripts: 

betaFunction.R: Beta Function used to fit to the data.

Host_BayesianAnalysis.R: Bayesian analysis for host data.

Parasite_BayesianAnalysis.R: Bayesian analysis for parasite data.


List of datafiles used for analysis: 

MeanTempTreatments.csv: mean temperature of the different treatments, calculated from temperature logger data that can be found in all_raw_temperatures_hobo.csv.

all_raw_temperatures_hobo.csv: merged data file of recordings made by the temperature loggers within each of the water baths. 

SporesNoMaleNA.csv: observations on infection status and spore burden.

sum_daphnia_noNA.csv: observations on Daphnia fecundity over the duration of the experiment.

