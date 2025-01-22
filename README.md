# Downing-et-al-2021
Code and data for: Downing PA, Griffin AS, &amp; Cornwallis CK. 2021. Hard-working helpers contribute to long breeder lifespans in cooperative birds. Philosophical Transactions of the Royal Society B, 376: 20190742.

Number of supplementary items: three
1. HWH_R_Code.R
2. HWH_Tables_S1-S5.xlsx
3. HWH_Data_Extraction.txt
4. HWH_Supp_Methods_Figures.pdf


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

File name: HWH_Tables_S1-S5.xlsx

This Excel document contains the following sheets:

- Table S1 (PRISMA workflow)
	+ systematic search for data on survival in cooperatively breeding birds
	+ search term: 'common name OR binomial species name OR synonym' (studies published up to June 2020)
	+ 144 data rows = 144 species
	+ column descriptions:\
		A. common = English name of each species\
		B. animal = latin binomial (matches the Jetz et al. nomenclature)\
		C. n.studies.WofS = number of studies on the species in Web of Science\
		D. n.studies.Scopus = number of studies on the species in Scopus\
		E. additional.studies = other studies on the species (e.g. from the grey literature)\
		F. screened = number of studies whose titles / abstract were examined for suitability\
		G. exluded = number of studies excluded based on the first screen\
		H. full.text = number of studies read in full to identify data\
		I. full.text.exluded = number of studies excluded after reading full text\
		J. n.studies = number of studies with suitable data\
		K. n.effect.sizes = number of effect sizes extracted from studies with suitable data


- Table S2 (data on breeder survival)
	+ annual survival estimates for female and male breeders in pairs and groups in 23 species
	+ 23 data rows = 23 species
	+ column descriptions:\
		A. animal = latin binomial (matches the Jetz et al. nomenclature)\
		B. common = English name of each species\
		C. female pair (%) = % of female breeders in pairs surviving between years\
		D. female group (%) = % of female breeders in groups surviving between years\
		E. male pair (%) = % of male breeders in pairs surviving between years\
		F. male group (%) = % of male breeders in groups surviving between years\
		G. N female pair = sample size used to estimate survival of female breeders in pairs\
		H. N female group = sample size used to estimate survival of female breeders in groups\
		I. N male pair = sample size used to estimate survival of male breeders in pairs\
		J. N male group = sample size used to estimate survival of male breeders in groups\
		K. female p value = test of survival difference between females in pairs and groups (from study)\
		L. male p value = test of survival difference between males in pairs and groups (from study)\
		M. source = where the data were extracted from in the study\
		N. notes = method used to estimate survival in the original study and further details\
		O. reference = study from which data were extracted


- Table S3 (data on breeder and helper feeding rates)
	+ feeding effort of female and male breeders in pairs and groups in 16 / 23 species for which survival data available
	+ helper feeding effort in 10 / 16 species for which breeder feeding data available
	+ 23 data rows = 23 species (note that there are missing data)
	+ column descriptions:\
		A. animal = latin binomial (matches the Jetz et al. nomenclature)\
		B. common = English name of each species\
		C. addresses question = which question in the main text the data address\
		D. units = the units in which each measurement was made\
		E. femalePair = mean feeding effort of breeding females in pairs\
		F. varFP = variance in feeding effort of breeding females in pairs\
		G. n.FemaleP = the number of breeding females in pairs studied\
		H. femaleGroup = mean feeding effort of breeding females in groups\
		I. varFG = variance in feeding effort of breeding females in groups\
		J. n.femaleG = the number of breeding females in groups studied\
		K. malePair = mean feeding effort of breeding males in pairs\
		L. varMP = variance in feeding effort of breeding males in pairs\
		M. n.maleP = the number of breeding males in pairs studied\
		N. maleGroup = mean feeding effort of breeding males in groups\
		O. varMG = variance in feeding effort of breeding males in groups\
		P. n.maleG = the number of breeding males in groups studied\
		Q. helper = mean feeding effort of helpers\
		R. varHelper = variance in helper feeding effort\
		S. n.groups = the number of groups studied to estimate helper feeding effort\
		T. notes =  where the data were extracted from in the study and further details\
		U. study = if the feeding effort data are from the same study as the survival data\
		V. reference = study from which data were extracted


- Table S4 (publication bias tests and heterogeneity estimates)
	+ parameter estimates from publication bias tests (see Downing_et_al_Code.R for further info.)
	+ heterogeneity estimates (see DowningetalRcode.R for further info.)
	+ 8 rows = output from 3 statistical models (one for each of the three effect sizes)
	+ variable names:
		I^2 = between-study variability / total variability


- Table S5 (results)
	+ parameter estimates from statistical models (see Downing_et_al_Code.R for further info.)
	+ 53 rows = output from 6 statistical models
	+ variable names:
		beta = fixed effect parameter estimate from the model
		lwr CI = lower 95% credible interval from the posterior distribution of the model
		upr CI = upper 95% credible interval from the posterior distribution of the model


- longData
	+ the data from Tables S2 and S3 organised in long format and analysed in R (see Downing_et_al_Code.R)
	+ export as a csv.doc and read into R using: data <- read.csv(".../data/longData.csv")
	+ 46 data rows = 23 species (each species has 2 estimates i.e. female and male breeder survival)
	+ column descriptions:\
		A. animal = latin binomial (matches the Jetz et al. nomenclature)\
		B. commonName = English name of each species\
		C. sex = biological sex of the estimates in columns D to M\
		D. sxP = % of breeders in pairs surviving between years\
		E. sxPN = sample size used to estimate survival of breeders in pairs\
		F. sxG = % of breeders in groups surviving between years\
		G. sxGN = sample size used to estimate survival of breeders in groups\
		H. brPFeed = mean feeding effort of breeders in pairs\
		I. brPVar = variance in feeding effort of breeders in pairs\
		J. brPN = the number of breeders in pairs studied\
		K. brGFeed = mean feeding effort of breeders in groups\
		L. brGVar = variance in feeding effort of breeders in groups\
		M. brGN = the number of breeders in groups studied\
		N. hpFeed = mean feeding effort of helpers\
		O. hpVar = variance in helper feeding effort\
		P. hpN = the number of groups studied to estimate helper feeding effort\


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

File name: HWH_R_Code.R

This R script contains all R code needed to replicate the analyses (including packages and functions)
Lettering is consistent with the methods section of the manuscript
Models cited in the manuscript are located in this code

- Data manipulation (lines 47 to 102)
- Publication bias tests and heterogeneity (lines 107 to 151)
- Question A. Does breeder survival increase with group size? (lines 157 to 267)
- Mean estimates for the D effect sizes (lines 273 to 382)
- Question B. Do increases in breeder survival depend on investment in care? (lines 388 to 497)
- Question C. Does breeder investment in parental care  depend on how much care helpers provide? (lines 503 to 612)


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

File name: HWH_Data_Extraction.txt

This plain text document contains details of how each effect size was calculated, organised alphabetically by latin name


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

File name: HWH_Supp_Methods_Figures.pdf

This pdf document contains Figures S1 to S3.


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
