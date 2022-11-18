rm(list = ls())
library(tidyverse)
library(lubridate)
library(readxl)
main_wd=ifelse(Sys.info()[['user']]=="solemon_pc", 'C:/Users/solemon_pc/Desktop/solemon/2022/raccolta_dati',
               ifelse(Sys.info()[['user']]=="e.armelloni", "C:/Users/e.armelloni/OneDrive/Lavoro/Solemon/github/SoleMon_project/OnBoard", 
                      ifelse(Sys.info()[['user']]=="Franc", "C:/Users/Franc/OneDrive/Desktop/solemon/2022/raccolta_dati", NA)))
setwd(main_wd)
source('R/functions_access_v2_0.R')

# set parameters
minilog='24K9475'
dates=1912

# load and format data 
mlog_dat=load_minilog(main_wd, minilog = minilog, dates=dates)

haul_order=read_excel(file.path(main_wd, 'data', 'haul_order.xlsx'))%>%
  dplyr::mutate(s_time=as.POSIXct(paste(day, substr(inizio, 11,19)), tz='UTC'))%>%
  dplyr::mutate(f_time=as.POSIXct(paste(day, substr(fine, 11,19)), tz='UTC'))%>%
  dplyr::select(-inizio, -fine)


# join
xdat_haul=full_join(mlog_dat, haul_order, by='day')%>%
  dplyr::filter(date >=s_time & date <=f_time)


# check
fine=xdat_haul%>%dplyr::group_by(haul)%>%
  slice_max(date)

inizio=xdat_haul%>%dplyr::group_by(haul)%>%
  slice_min(date)

difftime(as.POSIXct(inizio$date), fine$date) # to verify that haul duration is fine


# summarise
haul_summary=xdat_haul%>%dplyr::group_by(haul)%>%dplyr::summarise(m_temp=mean(temp), m_sal=mean(sal),
                                                      min_t=min(temp), max_t=max(temp),
                                                      min_sal=min(sal),
                                                      max_sal=max(sal))

# plot

p1=ggplot(data=xdat_haul)+
  geom_point(aes(x=date, y=temp))+
  facet_wrap(~haul,scales='free')+
  ggtitle(paste(minilog, dates, 'temp', sep='_'))

ggsave(plot=p1,file.path(main_wd, 'output', 'minilog', paste0(paste(minilog, dates, sep='_'), 'temp_profile.png')),
       width=20, height = 15, units='cm')

p2=ggplot(data=xdat_haul)+
  geom_point(aes(x=date, y=sal))+
  facet_wrap(~haul,scales='free')+
  ggtitle(paste(minilog, dates,'sal', sep='_'))

p3=ggplot(data=xdat_haul)+
  geom_point(aes(x=date, y=depth))+
  facet_wrap(~haul,scales='free')+
  ggtitle(paste(minilog, dates, 'depth', sep='_'))


ggsave(plot=p2,file.path(main_wd, 'output', 'minilog', paste0(paste(minilog, dates, sep='_'), 'sal_profile.png')),
       width=20, height = 15, units='cm')
ggsave(plot=p3,file.path(main_wd, 'output', 'minilog', paste0(paste(minilog, dates, sep='_'), 'dep_profile.png')),
       width=20, height = 15, units='cm')

#ggplot(data=xdat_haul%>%dplyr::group_by(haul, depth)%>%
#         dplyr::summarise(sal=mean(sal), depth=-depth))+
#  geom_point(aes(x=sal, y=depth))+
#  facet_wrap(~haul,scales='free')
#
## profilo CTD: considerare anche calo e salpo
#
#haul_order$s_time=haul_order$s_time-(5*60)
#haul_order$f_time=haul_order$f_time+(5*60)
#
#
#xdat_CTD=full_join(xdat, haul_order, by='day')%>%
#  dplyr::filter(date >=s_time & date <=f_time)
#
#
#ggplot(data=xdat_CTD%>%
#         dplyr::filter(depth>=0)%>%dplyr::group_by(haul, depth)%>%
#         dplyr::summarise(sal=mean(sal), depth=-depth))+
#  geom_point(aes(x=sal, y=depth))+
#  facet_wrap(~haul,scales='free')
#