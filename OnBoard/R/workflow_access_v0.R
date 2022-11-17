# This script controls the functions needed to preocess solemon data collected onboard.
# it runs on R 4.0.5 32 bit version due to compatibility issue with RODBC package needed to read from .accdb
# Last update: 14.11.2022
# function code is stored in the 'functions_access_v2_0.R' file. Functions are accompained by explanatory text.

rm(list = ls())
main_wd="C:/Users/e.armelloni/OneDrive/Lavoro/Solemon/github/SoleMon_project/OnBoard"    # paste working directory where xxx is located
setwd(main_wd)
source('R/functions_access_v2_0.R')

# inspect lists of target species and shells
unique(target_species$species_name) # these are the species for which you collect individual length AND individual weight
shells  # these are the species for which you collect total weight and total number 

# set parameters
haul=22 # single haul OR 'all'
db='test' # to be specified only for single hauls
updateID='N'
area_sepia='D'
year=2021
area='TEST'

# single haul application ####
# function1 extract data from access db and format them
hauldata=function1(haul=haul, 
                   db=db,
                   year=year)# extract and format data

# function 2: perform checks
function2(xdat=hauldata, 
          haul=haul)

# function 3: format data to trust format
trustdat=function3(xdat=hauldata[[1]], 
                  haul=haul, 
                  year = year, 
                  weight_not_target = hauldata[[2]],  
                  subsamples_target=hauldata[[3]],
                  catch_sample_disattivati = catch_sample_disattivati) # function 2

# function4: save PDF
function4(trustdat = trustdat, 
          year=year,
          area = area,
          haul=haul)

# multi-haul application: need extra file ####
haul_summary=read_excel("data/haul_order.xlsx")
haul_summary=haul_summary[1:5,]

for(xhaul in 1:nrow(haul_summary)){
  
  
  # loop parameters
  haul=haul_summary[xhaul,]$haul
  db=haul_summary[xhaul,]$DB
  area=haul_summary[xhaul,]$country
  
  cat('processing haul no.', haul, '(', xhaul,'/', nrow(haul_summary),')' )
  
  # function1 extract data from access db and format them
  hauldata=function1(haul=haul, 
                     db=db,
                     year=year)# extract and format data
  
  # function 2: perform checks
  function2(xdat=hauldata, 
            haul=haul)
  
  # function 3: format data to trust format
  trustdat=function3(xdat=hauldata[[1]], 
                     haul=haul, 
                     year = year, 
                     weight_not_target = hauldata[[2]],  
                     subsamples_target=hauldata[[3]],
                     catch_sample_disattivati = catch_sample_disattivati) # function 2
  
  # function4: save PDF
  function4(trustdat = trustdat, 
            year=year,
            area=area,
            haul=haul)


}
