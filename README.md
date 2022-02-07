# HostParasite


README for: 

Alternate patterns of temperature variation bring about very different disease outcomes at different mean temperatures
DOI: 10.7554/eLife.72861

list of authors: Charlotte Kunze 1,2*, Pepijn Luijckx 2*#, Andrew L. Jackson 2, Ian Donohue 2

1 Institute for Chemistry and Biology of the Marine Environment [ICBM], Carl-von-Ossietzky, University Oldenburg
2 Department of Zoology, School of Natural Sciences, Trinity College Dublin, Dublin, Ireland


These datasets contain observations on host fecundity, infection status and severity of infection collected during an experiment on Daphnia magna and its parasite Odospora colligata at School of Natural Sciences at the Trinity College Dublin in spring 2019. 
In this study we compared two patterns of temperature variation, a diurnal temperature fluctuation (of ±3°C in a 12:12 cycle) and a heat wave (of +6°C with 3 days duration 20 days after exposure to the parasite), with a constant treatment, which had the same mean temperature over the whole experimental period. 
Data were analyzed using R version 3.6.1(R Core Team, 2018) interfacing with JAGS. A Beta Function was fitted to each of our different fitness estimates (that is, host fecundity, parasite infectivity and burden) for each of the three temperature regimes using a baysian framework.


##
Datasets and R scripts used for the analyses: 

List of RScripts: 

*betaFunction.R*: Beta Function used to fit to the data.

*Host_BayesianAnalysis.R*: Bayesian analysis for host data.

*Parasite_BayesianAnalysis2.R*: Bayesian analysis for parasite data.

*HOBO.R*: R script to calculate mean temperature and merge all raw temperature files in the Rawtemperature folder.

###

List of datafiles used for analysis: 

*MeanTempTreatments.csv*: mean temperature of the different treatments, calculated from temperature logger data that can be found in all_raw_temperatures_hobo.csv.

*SporesNoMaleNA.csv*: observations on infection status and spore burden.

*sum_daphnia_noNA.csv*: observations on Daphnia fecundity over the duration of the experiment.

Rawtemperature folder: comprises all raw temeprature files for each treatment temperature logger. Here, *all_raw_temperatures_hobo.csv* contains all merged data files of recordings made by the temperature loggers within each of the water baths. 


###
Explanations: 

*MeanTempTreatments.csv*: mean temperature of the different treatments, calculated from temperature logger data that can be found in all_raw_temperatures_hobo.csv.

	meanTemp: mean temperature over the duration of the experiment.
  
	treatment: treatment information, either constant, pulse (or heatwave) or fluctuation corresponding to the used temperature regimes in the experiment.


###
*all_raw_temperatures_hobo.csv*: merged data file of recordings made by the temperature loggers within each of the water baths. 
	
	no: number of the measurement.  
  
	date: date and time of the measurement.

	expected_temp: aimed temperatues

	temp: measured temperature in °C

	id: name of the specific HOBO logger which consist of the treatment and target temperature (e.g., constant 10)


###

*SporesNoMaleNA.csv*: observations on infection status and spore burden

	no: unique sample number from 1 to 516.

	treatment: the temperature regime used, abbreviations, CS = constant, PULSE = heat wave, FLU = fluctuating. 

	exposed: exposure to the parasite, U = unexposed (control), I = exposed.

	temperature: target temperature (10, 13, 16, 19, 22, 25 and 28).  Note: for fluctuating regime this is the lower temperature.

	meantemp: mean temperature, different from target temperature for fluctuation regime which alternates been the target temperature and 6℃ higher.

	realtemp: measured mean temperature over the experimental period provided by temperature loggers see MeanTempTreatments.csv.

	replicate: replicate number from 1 to 18.

	death: date of death (last day 08/05/2019)

	lastday: was the animal alive the last day of the experiment? Either 0 (dead) or 1 (alive). Used for spore analysis.

	infect: infection status, either 0 (no parasite) or 1 (infected).

	no_spore: number of spores observed in the dissected Daphnia gut.

	size: size of the Daphnia measured under the binocular (for the real size this value  must be transformed by the magnification factor). Not used for the paper and is correlated with fecundity.

	comment: gives additional information about the Daphnia during their live span.

###

*sum_daphnia_noNA.csv*: Daphnia data

	no: unique sample number from 1 to 516.

	treat: the temperature regime used, abbreviations, CS = constant, PULSE = heat wave, FLU = fluctuating. 

	inf: exposure to the parasite, U = unexposed (control), I = exposed.

	mean_temp: measured mean temperature over the experimental period provided by temperature loggers see MeanTempTreatments.csv.

	temp: target temperature (10, 13, 16, 19, 22, 25 and 28)

	treatment_id: unique treatment ID comprised of Infection status, temperature, treatment, replicate number

	repl: replicate number from 1 to 18.

	sum: sum of offspring produced over the experimental period.
###


