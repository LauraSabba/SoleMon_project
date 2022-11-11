setwd("C:/Users/e.armelloni/OneDrive/Lavoro/Solemon/AtSeaData/2021/code")
rm(list = ls())
library(tidyverse)
library(readxl)
library(RODBC)
'%ni%'=Negate('%in%')

# set parameters
haul_date='haul_sheets'
updateID='Y'
DRIVERINFO <- "Driver={Microsoft Access Driver (*.mdb, *.accdb)};"
area_sepia='D'
#haul='33'

# load data ####
species_list=read_excel("../other_data/species_list.xlsx")
mat_list=read_excel("../other_data/maturity_stages.xlsx")
lw_pars=read_excel("../other_data/lw_pars.xlsx")
haul_order <- read_excel("../haul_order.xlsx")
catch_sample_disattivati <- read_excel("../other_data/catch_sample_disattivati.xlsx")

# set target species ####
  target_species=data.frame(species_name=c("SOLEVUL","SOLEAEG","PLATFLE","SCOHRHO","PSETMAX",
                                           "MERLMER","MULLBAR",
                                           "RAJAAST","RAJACLA","RAJAMIR","TORPMAR",
                                           "SCYOCAN", "SCYIOSTE",
                                           "PAPELON","MELIKER","NEPHNOR","SQUIMAN",
                                           "SEPIOFF",
                                           "PECTJAC","AEQUOPE","CHLAPRO", "CHLAVAR"), target=1)

  shells=c('MUREBRA', 'HEXATRU', 
           'OSTREDU', 'CRASGIG',
           'MYTGALL', 'NATISTE',
           'RAPAVEN', 'GALEECH')


# format species names ####
species_list=species_list[,c('Species','Medits')]
names(species_list)[2]='species_name'
species_list=left_join(species_list, target_species,by='species_name')%>%
  replace(is.na(.),0)


# Run

xfiles=list.files(path=file.path("..", haul_date))

xfiles=data.frame(file=xfiles)%>%
  dplyr::mutate(haul=str_remove(str_remove(file, "cala"), '.xlsx'))%>%
  left_join(haul_order, by='haul')%>%arrange(id)

for(k in 1:nrow(xfiles)){
  cat(k)
  FishIDs=read_excel("../other_data/fishID.xlsx")
  
  ### New part to connect with accdb
  
  MDBPATH <- paste0("C:/Users/e.armelloni/OneDrive/Lavoro/Solemon/AtSeaData/2021/Maschera inserimento SOLEMON_",xfiles[k,]$DB,".accdb")
  PATH <- paste0(DRIVERINFO, "DBQ=", MDBPATH)
  ## Establish connection
  channel <- odbcDriverConnect(PATH)
  
  # look for tables
  #acctables=sqlTables(channel)
  
  ## Load data into R dataframe
  xdat <- sqlQuery(channel,
                 paste0("SELECT *
FROM [", paste0('cala', xfiles[k,]$haul), "]
ORDER BY [ID]"),
                 stringsAsFactors = FALSE)
  
  close(channel)
  rm(channel)
 # xdat2=read_excel(file.path("..", haul_date, xfiles[k,'file']))%>%
 #  arrange(ID)
  xdat=as_tibble(xdat)
  xdat[is.na(xdat$weight_g),'weight_g']=0
  xdat[is.na(xdat$Mat),'Mat']=0

if(length(grep("TAG",names(xdat)))==0){
  
  xdat$TAG=0
}

xdat$total_number=ifelse(is.na(xdat$total_number), 0, xdat$total_number)



 haul=xfiles[k,]$haul
  
  # fill gaps ####
  for(i in 2:nrow(xdat)){
    
    # Fill gaps in gear
    if(is.na(xdat[i,]$gear)){xdat[i,]$gear=xdat[i-1,]$gear}
    # Fill gaps in species name
    if(is.na(xdat[i,]$species_name)){xdat[i,]$species_name=xdat[i-1,]$species_name}
    # Fill gaps in Sex
    if(xdat[i,]$species_name=='MELIKER' & is.na(xdat[i,]$Sex)){xdat[i,]$Sex=xdat[i-1,]$Sex}
    if(xdat[i,]$species_name=='SQUIMAN' & is.na(xdat[i,]$Sex)){xdat[i,]$Sex=xdat[i-1,]$Sex}
    
  }
  
  # attach info on target species
  xdat=left_join(xdat, species_list, by='species_name')
  if(nrow(xdat[is.na(xdat$target),])>0){
    xdat[is.na(xdat$target),]$target=0
  }
  
  
  # oto & gen ID ####
  soleaID=FishIDs[FishIDs$type=='solea',]$fishID+1
  first_soleaID=soleaID
  
  ELASID=FishIDs[FishIDs$type=='ELAS',]$fishID+1
  first_ELASID=ELASID
  
  RHOID=FishIDs[FishIDs$type=='RHO',]$fishID+1
  first_RHOID=RHOID
  
  PSETTAID=FishIDs[FishIDs$type=='PSETTA',]$fishID+1
  first_PSETTAID=PSETTAID
  
  sepiaID=FishIDs[nchar(FishIDs$type)==7 & substr(FishIDs$type,7,7)==tolower(area_sepia),]$fishID+1
  first_sepiaID=sepiaID
  
  
  xdat$Oth=xdat$oto_id
  xdat$Gen.Samp=xdat$genetic_id
  xdat$fish_ID=as.character(NA)
  
  
  # oto 
  for(i in 1:nrow(xdat)){
    if(xdat[i,]$species_name=='SOLEVUL' | xdat[i,]$species_name=='SOLEAEG'){
      
      if(is.na(xdat[i,]$oto_id)==FALSE){
        soleaID=first_soleaID
        xdat[i,]$fish_ID=paste0('SS', soleaID)
        first_soleaID=first_soleaID+1
      }
      
    }
    
    if(xdat[i,]$species_name=='SCOHRHO'){
      
      if(is.na(xdat[i,]$oto_id)==FALSE){
        RHOID=first_RHOID
        xdat[i,]$fish_ID=paste0('SR', RHOID)
        first_RHOID=first_RHOID+1
      }
      
    }
    
    if(xdat[i,]$species_name=='PSETMAX'){
      
      if(is.na(xdat[i,]$oto_id)==FALSE){
        PSETTAID=first_PSETTAID
        xdat[i,]$fish_ID=paste0('PM', PSETTAID)
        first_PSETTAID=first_PSETTAID+1
      }
      
    }
    
  }
  
  # gen
  for(i in 1:nrow(xdat)){
    if(xdat[i,]$species_name=='SEPIOFF'){
      
      if(is.na(xdat[i,]$genetic_id)==FALSE){
        sepiaID=first_sepiaID
        xdat[i,]$fish_ID=paste0(area_sepia, sepiaID)
        first_sepiaID=first_sepiaID+1}
      
    }
    
    if(xdat[i,]$species_name %in% c('RAJAAST', 'RAJACLA', 'RAJAMIR', 'TORPMAR')){
      
      if(is.na(xdat[i,]$genetic_id)==FALSE){
        ELASID=first_ELASID
        xdat[i,]$fish_ID=paste0('Elas', ELASID)
        first_ELASID=first_ELASID+1}
      
    }
    
  }
  
  # update id sheet
  if(first_soleaID != soleaID){ FishIDs[FishIDs$type=='solea',]$fishID=first_soleaID-1 }
  if(first_ELASID != ELASID){ FishIDs[FishIDs$type=='ELAS',]$fishID=first_ELASID-1 }
  if(first_RHOID != RHOID){ FishIDs[FishIDs$type=='RHO',]$fishID=first_RHOID-1 }
  if(first_PSETTAID != PSETTAID){ FishIDs[FishIDs$type=='PSETTA',]$fishID=first_PSETTAID-1 }
  if(first_sepiaID != sepiaID){FishIDs[nchar(FishIDs$type)==7 & substr(FishIDs$type,7,7)==tolower(area_sepia),]$fishID=first_sepiaID-1}
  if(updateID=='Y'){writexl::write_xlsx(FishIDs, "../other_data/fishID.xlsx")}
  
  
  
  # maturity ####
  xdat[xdat$target!=1,]$Sex='I'
  xdat[xdat$target!=1,]$Mat=NA
  xdat[xdat$target==1 & xdat$Sex %ni% c('M', 'F', 'I'),]$Sex='N'
  xdat$Mat=ifelse(xdat$species_name=='MELIKER' & xdat$Sex=='F' & xdat$Mat==0, 1, xdat$Mat)
  xdat$Mat=ifelse(xdat$species_name=='PAPELON' & xdat$Sex=='F' & xdat$Mat==0, 1, xdat$Mat)
  xdat[xdat$Mat%ni%c(1,2,3,4,5,6),]$Mat=NA
  
  xdat$Mat=as.character(xdat$Mat)
  for(i in 1:nrow(xdat)){
    if(xdat[i,]$species_name %in% mat_list$SPECIES){
      
      xmat=mat_list[mat_list$SPECIES==xdat[i,]$species_name,]
      xmat=xmat[xmat$SEX==xdat[i,]$Sex,]
      
      if(nrow(xmat)==1 & xdat[i,]$Mat %in% c(1,2,3,4,5,6)){
        
        xdat[i,]$Mat=paste(xmat$SCALE, xdat[i,]$Mat, sep='-')
        
      }
      
    }
    
  }
  
  
  # weight data for non target ####
  shells_w=xdat[xdat$species_name%in% shells,]
  
  weight_not_target=xdat[xdat$target!=1 & xdat$species_name%ni% shells,]%>%
    dplyr::group_by(gear, species_name)%>%
    dplyr::summarise(w=sum(weight_g)/1000, n=n())
  
  if(nrow(shells_w)>0){
    shells_w=shells_w%>%
      dplyr::group_by(gear, species_name)%>%
      dplyr::summarise(w=sum(weight_g)/1000, n=sum(total_number))
    weight_not_target=rbind(weight_not_target, shells_w)
    
  }
  
  xdat[xdat$target!=1,]$weight_g=-1
  
  # remove shells
  xdat=xdat[xdat$species_name%ni% shells, ]
  
  
  # weight data for sub-samples ####
  subsamples=xdat[xdat$total_number >0,]
  
  if(nrow(subsamples)>0){
    subsamples=subsamples%>%
      dplyr::group_by(gear, species_name)%>%
      dplyr::summarise(w=sum(weight_g)/1000, n=sum(total_number))
  }
 
  
  # remove subsamples
  xdat=xdat[xdat$total_number ==0,]
  
  # Fill gaps for target ####
  for(i in 1:nrow(xdat)){
    #rajas
    if(xdat[i,]$species_name %in% c('RAJAAST', 'RAJACLA', 'RAJAMIR', 'TORPMAR')){xdat[i,]$lenght_mm=as.numeric(strsplit(xdat[i,]$Notes, '/' )[[1]][1])}
    # lw
    if(xdat[i,]$species_name %in% lw_pars$species_name){
      # adjust weight
      if(xdat[i,]$weight_g==0|is.na(xdat[i,]$weight_g)){
        if(xdat[i,]$species_name%in% c('MELIKER', 'MERLMER')){
          a=lw_pars[lw_pars$species_name==xdat[i,]$species_name & lw_pars$sex==xdat[i,]$Sex,]$a
          b=lw_pars[lw_pars$species_name==xdat[i,]$species_name& lw_pars$sex==xdat[i,]$Sex,]$b
        }else{
          a=lw_pars[lw_pars$species_name==xdat[i,]$species_name,]$a
          b=lw_pars[lw_pars$species_name==xdat[i,]$species_name,]$b
        }
        if(xdat[i,]$species_name=='SQUIMAN'|xdat[i,]$species_name=='MELIKER'|xdat[i,]$species_name=='PAPELON'){
          xdat[i,]$weight_g=round((a*xdat[i,]$lenght_mm^b))
        }
        else{
          xdat[i,]$weight_g=round((a*xdat[i,]$lenght_mm^b)/1000)
        }
      }
      # adjust length
      if(xdat[i,]$lenght_mm==0|is.na(xdat[i,]$lenght_mm)){
        if(xdat[i,]$species_name%in% c('MELIKER', 'MERLMER')){
          a=lw_pars[lw_pars$species_name==xdat[i,]$species_name & lw_pars$sex==xdat[i,]$Sex,]$a
          b=lw_pars[lw_pars$species_name==xdat[i,]$species_name& lw_pars$sex==xdat[i,]$Sex,]$b
        }else{
          a=lw_pars[lw_pars$species_name==xdat[i,]$species_name,]$a
          b=lw_pars[lw_pars$species_name==xdat[i,]$species_name,]$b
        }
        if(xdat[i,]$species_name=='SQUIMAN'|xdat[i,]$species_name=='MELIKER'|xdat[i,]$species_name=='PAPELON'){
          xdat[i,]$lenght_mm=round((xdat[i,]$weight_g/(a))^(1/b))
        }
        else{
          xdat[i,]$lenght_mm=round((xdat[i,]$weight_g/(a/1000))^(1/b))
        }
      }
    }
  }
  
  # Format: moduli Trust ####
  
  # bio data
  bio_data_trust=xdat
  bio_data_trust$Survey='SOLEMON2021'
  bio_data_trust$Area=xfiles[k,]$country
  bio_data_trust$Station=haul
  bio_data_trust$Gear=ifelse(bio_data_trust$gear=='A', '1-RAP', '2-RAP')
  bio_data_trust$SampN=1
  bio_data_trust$SpecN=1
  bio_data_trust$Age=NA
  bio_data_trust=bio_data_trust[, c('Survey', 'Area', 'Station', 'Gear', 'species_name', 'SampN', 'SpecN', "lenght_mm", "weight_g", "Sex", 'Mat', 'Oth','Age', 'fish_ID', 'Gen.Samp','Notes', 'TAG')]
  names(bio_data_trust)=c('Survey',	'Area',	'Station'	,'Gear',	'SpecCode',	'SampN'	,'SpecN',	'L(mm)'	,'W(g)',	'Sex'	,'MatStage',	'Oth',	'Age'	,'FishID'	,'Gen.Samp.',	'Notes', 'TAG')
  
  bio_data_trust=bio_data_trust%>%
    dplyr::mutate(MatStage=ifelse(nchar(MatStage)<3,NA,MatStage))
  
 bio_data_trust$Oth=ifelse(is.na(bio_data_trust$Oth), 0,1)
 if(haul==18){
    bio_data_trust[bio_data_trust$SpecCode=='SOLEVUL'&!is.na(bio_data_trust$FishID),]$FishID=NA
  }
  
  writexl::write_xlsx(bio_data_trust, file.path("..", "TrustData","Bio" ,paste0("Bio_Trust_", haul, ".xlsx")))
  
  # catch sample data
  catch_samples=bio_data_trust%>%
    #dplyr::rename('Code'=SpecCode)%>%
    #left_join(species_list[,c('Species', "species_name")] %>%
    #            dplyr::rename('Code'=species_name, 'SpeciesSN'=Species),by='Code')%>%
    #dplyr::distinct(Survey, Area, Station, Gear, SpeciesSN, Code)%>%
    dplyr::distinct(Survey, Area, Station, Gear, SpecCode)%>%
    dplyr::mutate(SampN=1,
                  'W(kg)'=-1,
                  SampMeth='SIMRANDO',
                  InUse='Y',
                  Picts=NA,
                  Notes=NA)
  if(haul %in% catch_sample_disattivati$Station){
    
    catch_sample_disattivati$Gear=ifelse(catch_sample_disattivati$Gear=='A', '1-RAP', '2-RAP')
    catch_samples[ catch_samples$SpecCode ==catch_sample_disattivati[catch_sample_disattivati$Station==haul,]$SpecCode &
                     catch_samples$Gear== catch_sample_disattivati[catch_sample_disattivati$Station==haul,]$Gear, ]$InUse='N'
  }
  

  writexl::write_xlsx(catch_samples, file.path("..", "TrustData","Catch_sample" ,paste0("Catch_sample_Trust_", haul, ".xlsx")))
  
  # catch data
  catch_target=xdat[xdat$target==1,]%>%
    dplyr::group_by(gear, species_name)%>%
    dplyr::summarise(w=sum(weight_g)/1000,
                     n=n())
  
  catch_file_trust=rbind(weight_not_target,catch_target)
  
  if(nrow(subsamples)>0){
    catch_file_trust=rbind(catch_file_trust, subsamples)
    catch_file_trust=catch_file_trust%>%
      dplyr::group_by(gear, species_name)%>%
      dplyr::summarise(w=sum(w), n=sum(n))
    }
  
  
  
  catch_file_trust=catch_file_trust%>%
    left_join(species_list[,c('Species', "species_name")],by = "species_name")
  
  catch_file_trust=catch_file_trust[,c("gear", "Species", "species_name", 'w','n' )]
  catch_file_trust$Survey='SOLEMON2021'
  catch_file_trust$Area=xfiles[k,]$country
  catch_file_trust$Station=haul
  catch_file_trust=catch_file_trust%>%dplyr::select(Survey, Area, Station, gear, Species,
                                   species_name,w,n)
  names(catch_file_trust)=c('Survey' , 'Area', 'Station', 'Gear', 'SpeciesSN', 'Code', 'W(kg)', 'Numb')
  catch_file_trust$RaisF=1
  catch_file_trust$Notes=NA
  catch_file_trust$UserIns=NA
  catch_file_trust$Gear=ifelse(catch_file_trust$Gear=='A', '1-RAP', '2-RAP')
  
  writexl::write_xlsx(catch_file_trust, file.path("..", "TrustData","Catch", paste0("Catch_Trust_", haul, ".xlsx")))
  
  rm(bio_data_trust, catch_file_trust,catch_samples, 
     catch_target, FishIDs, weight_not_target, xdat,xmat,shells_w,
     a,b,ELASID,RHOID,PSETTAID, sepiaID,soleaID,
     first_ELASID,first_RHOID, first_PSETTAID, first_sepiaID, first_soleaID,
     i)

}






