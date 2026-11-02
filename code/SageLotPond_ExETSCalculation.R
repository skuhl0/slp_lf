
#########################Lateral Flux Data Processing############################

# This script uses EXO biogeochemical parameters and water flow data
# to calculate lateral fluxes of DIC and DOC following the ExETS method. DIC 
# and DOC concentrations are calculated using predictive power law equations. 
# Flow data is corrected for overland flow and to balance the water budget. 
# Instantaneous DIC and DOC fluxes are calculated as the product of DIC/DOC 
# concentration and water flow. Tidal cycles are identified and tidal 
# characteristics and fluxes are calculated. Monthly average fluxes are 
# calculated for total carbon flux and marsh-derived carbon flux (corrected
# for groundwater inputs). Flux relationships with air temperature are used
# to estimate fluxes where they can not be calculated. 

library(tidyverse)
library(zoo)

# Sage Lot Pond Creek sensor data with high and low tide extrema and 
# individual tidal cycles designated
data <- read_csv('data/SLP_creek_data.csv')%>%
  mutate(month = month(datetime),year = year(datetime))

# removing flow data during unbalanced tides to avoid bias. These were 
# determined by performing tidal calculations on the whole time series, then 
# determining which tides had a net flow outside the limit (1000 m3). These 
# tides are removed (flow set to NA) for final calculations
rm_tides <- c(90,92,97,167,174,193,338,363,414,416,418,420,429,517,572,585,595,
              612, 638,654,665,688,723,724,725,730,746,985,1042,1055,1057,1094,
              1111,1115,1134,1135,1170,1171,1189,1203,1344,1367,1528,1628,1638,
              1639,1662,1663, 1809,1811,1813,1869,1918,2128,2132,2134,2135,2136,
              2188,2192,2212,2234,2245,2265)

rm_inst <- which(data$tide_num %in% rm_tides)

data$flow[rm_inst] <- NA

####################DIC and DOC Concentration Calculations#######################

# calculating DIC concentration (umol/L) from predictive equation 
data <- mutate(data, temp_adj = ifelse(temp<=0,NA,temp),
               DIC = 10^(3.27035)*sal^(0.60693)*ph^(-0.88417)*
                 domgL^(-0.06894)*temp_adj^(-0.03633))

# removing calculated DIC concentrations where predictors are out of bounds
aa <- which(data$ph<6.8 | data$ph>8.5)
bb <- which(data$sal<19 | data$sal>32)
cc <- which(data$temp_adj<2.5 | data$temp_adj>30)
dd <- which(data$domgL<0.28 | data$domgL>11.85)
dic_filter <- c(aa,bb,cc,dd)
data$DIC[dic_filter] <- NA

# converting DIC conc to g/m3
data$dic_gm3 <- data$DIC*0.012

# calculating DOC concentration (umol/L) from predictive equation 
data$DOC1 <- 10^(3.02976)*data$temp_adj^(0.50412)*data$fdom^(0.20955)*data$ph^(-1.75138)

# removing calculated DOC concentrations where predictors are out of bounds
ee <- which(data$fdom<0.9 | data$fdom>80)
doc_filter_fdom <- c(aa,cc,ee)
data$DOC1[doc_filter_fdom] <- NA

# calculating DOC concentration (umol/L) from second predictive equation when fDOM is missing
data$DOC2 <- 10^(4.68979)*data$temp_adj^(0.63080)*data$ph^(-3.50524)

# removing calculated DOC concentrations where predictors are out of bounds
doc_filter_no_fdom <- c(aa,cc)
data$DOC2[doc_filter_no_fdom] <- NA

# selecting equation 1 value when available and 2 value otherwise
data$DOC <- ifelse(is.na(data$DOC1),data$DOC2,data$DOC1)

# designating equation used as 1 or 2, and 0 for NA values
data$DOCeq <- ifelse(is.na(data$DOC1),2,1)
data$DOCeq[which(is.na(data$DOC))] <- 0

# determining tides where both equations 1 and 2 are used
TF_eq <- group_by(data,tide_num)%>%
  summarize(DOCeq = sum(unique(DOCeq)))%>%
  filter(DOCeq==3)

# replacing some equation 1 values with 2 values so that no individual tidal 
# cycle uses both equations 
data$DOC[which(data$tide_num%in%TF_eq$tide_num)] <- 
  data$DOC2[which(data$tide_num%in%TF_eq$tide_num)]

# updating equation numbers based on replaced data
data$DOCeq[which(data$tide_num%in%TF_eq$tide_num)] <- 2
data$DOCeq[which(is.na(data$DOC))] <- 0

# converting DOC conc to g/m3
data$doc_gm3 <- data$DOC*0.012

########################Water Flow Calculations and Corrections##################

###Overland Flow Corrections###

# create column containing overland flow corrected flow values
data$flow_corr <- data$flow

# selecting values where quadratic overland flow correction 
# is applied. flood correction applied to positive flows, ebb correction 
# applied to negative flows
flood_excess <- which(data$flow>0 & data$h>0.25 & data$h<=0.5)
ebb_excess <- which(data$flow<0 & data$h>0.26 & data$h<=0.5)

# applying ebb and flood overland flow corrections 
data$flow_corr[flood_excess] <- (-15.4*data$h[flood_excess]^2 + 
                                   20.4*data$h[flood_excess]- 3.11)*data$flow[flood_excess]
data$flow_corr[ebb_excess] <- (12.2*data$h[ebb_excess]^2 + 
                                 0.275*data$h[ebb_excess]+ 0.097)*data$flow[ebb_excess]

# selecting values where linear overland flood correction is 
# applied, flood and ebb use same correction for water level > 0.5 m
excess_linear <- which(data$h>0.5)

# applying linear overland flow correction 
data$flow_corr[excess_linear] <- 
  (data$h[excess_linear]*9.170833-1.320587)*data$flow[excess_linear]

####Water Budget Balance/Net Flow Corrections####

# select overland-corrected flow data and salinity for salt balance calculation
Qall<-data$flow_corr
salall<-data$sal

# create vectors for calculation parameters for each month and year
dates <- unique(data[c('month', 'year')])
L <- nrow(dates)
q_net_in <- rep(NA, L)
q_net_out <- rep(NA, L)
s_net_in <- rep(NA, L)
s_net_out <- rep(NA, L)
Qr_TEF1 <- rep(NA, L)

# perform calculation for each month of each year
for (tindex in 1:L){
  
  # position of values from current month in dataset
  dn_range <- which(data$month==dates$month[tindex] & data$year==dates$year[tindex]) 
  
  # make flow and salinity vectors for data in time window
  q_temp <-  Qall[dn_range]
  s_temp <-  salall[dn_range]
  
  # set up parameters for Knudsen calculation
  zin <- which(q_temp>0) # inflowing data positions
  zout <- which(q_temp<0) # outflowing data positions
  tsize <- length(dn_range) # length of current month of data
  dsbin <- 0.1 # salinity bin width
  sbin <- round(seq(0,60,0.1),1) #salinity bin values
  
  # round salinity data to sort into bins
  s_temp_round <- round(s_temp, digits=1)
  
  # create vectors for flow in each bin
  q_temp_in <- rep(0, length(sbin))
  q_temp_out <- rep(0, length(sbin))
  
  # sort flow data into salinity bins
  for(bb in 1:length(sbin)){
    
    # positions of data within each salinity bin 
    bin_range <- which(s_temp_round == sbin[bb])
    
    # split binned data into inflowing and outflowing values
    bin_range_in <- intersect(bin_range, zin)
    bin_range_out <- intersect(bin_range, zout)
    
    # determine number of entries within the bin
    in_length <- length(bin_range_in)
    out_length <- length(bin_range_out)
    
    # create a vector of flow data and zeros equivalent to the total 
    # number of entries in the month
    q_in_zeros <- c(q_temp[bin_range_in],rep(0,(tsize-in_length)))
    q_out_zeros <- c(q_temp[bin_range_out],rep(0,(tsize-out_length)))
    
    # average inflowing and outflowing data within the bin
    q_temp_in[bb] <- mean(q_in_zeros,na.rm=TRUE)
    q_temp_out[bb] <- mean(q_out_zeros,na.rm=TRUE)
    
  }
  
  # sum inflowing and outflowing data in each bin
  q_temp_sum <- q_temp_in+q_temp_out
  
  # calculate overall water exchange in and out
  q_net_in[tindex] = sum(q_temp_sum[which(q_temp_sum > 0)], na.rm = TRUE)  
  q_net_out[tindex] = sum(q_temp_sum[which(q_temp_sum < 0)], na.rm = TRUE)
  
  # calculate salt flux by total exchange flow in and out
  q_temp_sum_in <- q_temp_sum
  q_temp_sum_in[which(q_temp_sum_in<0)] <- 0
  qs_net_in = sum(q_temp_sum_in * sbin, na.rm = TRUE)  
  s_net_in[tindex] <- qs_net_in/q_net_in[tindex] # salt flux in
  
  q_temp_sum_out <- q_temp_sum
  q_temp_sum_out[which(q_temp_sum_out>0)] <- 0
  qs_net_out = sum(q_temp_sum_out * sbin, na.rm = TRUE)
  s_net_out[tindex] <- qs_net_out/q_net_out[tindex] # salt flux out
  
  # calculate net freshwater discharge as a function of salt flux and exchange flow
  Qr_TEF1[tindex]<-(s_net_in[tindex]-s_net_out[tindex])*abs(q_net_out[tindex])/s_net_in[tindex]
  
}

# salt balance calculation results, modifying sign to match convention of 
# creek data
net_flow_calc <- cbind(dates, Qr_TEF1)%>%
  mutate(flow_m3s = -1*Qr_TEF1)

## calculate the data coverage for each month for the net flow calculation
net_flow_cover <- group_by_at(data, c("year","month"))%>%
  summarize(cover = length(which(!is.na(flow)&!is.na(sal)))/n(),
            mean_in = sum(flow_corr[which(flow_corr>0&!is.na(sal))]*300,na.rm=TRUE),
            mean_out = sum(flow_corr[which(flow_corr<0&!is.na(sal))]*300,na.rm=TRUE),
            ratio = mean_in/mean_out)
net_flow_calc$cover <- net_flow_cover$cover
net_flow_calc$ratio <- net_flow_cover$ratio

# remove values where there is bias in flow data coverage
net_flow_calc$flow_m3s[which(net_flow_calc$ratio>-0.9)] <- NA
net_flow_calc$flow_m3s[which(net_flow_calc$cover<0.2)] <- NA

# use mean freshwater discharge value where monthly value could not be calculated
net_flow_calc$flow_m3s[which(is.na(net_flow_calc$flow_m3s))] <-
  mean(net_flow_calc$flow_m3s,na.rm=TRUE)

# create column for flow data with corrected water budget
data$flow_corr_shift <- data$flow_corr

# create dataframe to store net flow correction factor
shift_corr <- tibble(month = rep(NA,nrow(dates)), year = rep(NA,nrow(dates)), 
                     factor = rep(NA,nrow(dates)))

# calculate corrected instantaneous flow values for each month
for(i in 1:nrow(shift_corr)){
  
  # current month and year
  M <- shift_corr$month[i] <- dates$month[i]
  Y <- shift_corr$year[i] <- dates$year[i]
  
  # creek data from current month
  month <- filter(data, month == M, year == Y)
  
  # net flow value from salt balance calculation for current month
  water_budget <- filter(net_flow_calc, year == Y, month == M)
  
  # calculation parameters for net flow correction factor (equation 3)
  ebb_mean <- mean(month$flow_corr[which(month$flow_corr<0)],na.rm=TRUE)
  flood_mean <- mean(month$flow_corr[which(month$flow_corr>0)],na.rm=TRUE)
  tot_mean <- mean(month$flow_corr,na.rm=TRUE)
  ebb_n <- length(which(month$flow_corr<0))
  flood_n <- length(which(month$flow_corr>0))
  tot_n <- ebb_n+flood_n
  
  # calculate correction factor
  shift_corr$factor[i] <- (water_budget$flow_m3s-(flood_n*flood_mean)/tot_n)*
    (tot_n/(ebb_n*ebb_mean))
  
  # apply correction factor to ebb (outflowing) data to determine corrected
  # instantaneous flow values
  ebbs <- which(data$flow_corr<0 & data$month == M & data$year == Y)
  data$flow_corr_shift[ebbs] <- data$flow_corr[ebbs]*shift_corr$factor[i]
  
}

#####################Calculating Instantaneous Carbon Fluxes#####################

# this function takes the "gap check" vector (0 in places where there is flux data,
# and 1 in places where flux = NA) and creates a corresponding vector that contains 
# lengths associated with each run of real values or gap (run of NAs)
gapcheck2length <- function(gap_check){
  
  # creates an object containing the value of each run (0 or 1) and the length of the run
  fluxgaps <- rle(gap_check)
  lengths <- fluxgaps$lengths
  values <- fluxgaps$values
  
  # creates a counter that shows the start position of each gap or series of real values
  lengths_sum <- cumsum(lengths)
  lengths_sum <- c(0,lengths_sum)
  gap_length <- rep(0,length(gap_check))
  
  # assign the length of the sequence to each entry that is part of that run
  for (i in 1:length(lengths)) {
    gap_length[(lengths_sum[i]+1):lengths_sum[i+1]] <- lengths[i]
  }
  
  return(gap_length)
  
}

# instantaneous DIC and DOC fluxes (g/s) calculated as concentration (g/m3)
# multiplied by flow (m3/s)
data$docflux <- data$doc_gm3*data$flow_corr_shift
data$dicflux <- data$dic_gm3*data$flow_corr_shift

###Interpolation of Missing C Flux Data###

# create a marker for gaps in DIC and DOC flux data
data$doc_gap_check <- ifelse(is.na(data$docflux),1,0)
data$dic_gap_check <- ifelse(is.na(data$dicflux),1,0)

# create columns that list the length of the gap each value is part of 
data$doc_gap_length <- gapcheck2length(data$doc_gap_check)
data$dic_gap_length <- gapcheck2length(data$dic_gap_check)

# perform spline interpolation to fill in gaps in flux
data$docinterp <- na.spline(data$docflux)
data$dicinterp <- na.spline(data$dicflux)

# set any interpolated values where the gap length is longer than 10 entries back to NA
data$docinterp[which(data$doc_gap_check == 1 & data$doc_gap_length > 10)] <- NA
data$dicinterp[which(data$dic_gap_check == 1 & data$dic_gap_length > 10)] <- NA

#########################Calculating Tidal Carbon Fluxes#########################

# calculate DIC/DOC fluxes + relevant tidal characteristics for each half tide 
HTfluxes <- group_by(data, halftide_num) %>%
  summarize(tide_num = tide_num[1],
            total_doc_interp = sum(docinterp)*5*60/4132,#total flux during tide g/m2
            total_dic_interp = sum(dicinterp)*5*60/4132,#total flux during tide g/m2
            HL = HL[length(HL)], # high/low extreme tide marker
            start = datetime[1]-hms("0:5:0"),
            end = datetime[length(datetime)],duration = (end-start)/3600,
            N = n(),tot_NAs_DOC = length(which(is.na(docflux))),
            tot_NAs_DIC = length(which(is.na(dicflux)))) %>%
  mutate(in_out = ifelse(HL == "H", "in","out"))

# create markers for half tides that contain good (complete/reliable) flux data
# good tide = 1, bad tide = 0 
HTfluxes <- mutate(HTfluxes, good_doc = 1, good_dic = 1)

# giving "bad" markers to half tides with incomplete C flux data
HTfluxes$good_doc[which(is.na(HTfluxes$total_doc_interp))] <- 0
HTfluxes$good_dic[which(is.na(HTfluxes$total_dic_interp))] <- 0

# giving "bad" markers to half tides where more than 10 entries were missing and 
# interpolated back in
HTfluxes$good_doc[which(HTfluxes$tot_NAs_DOC>10)] <- 0
HTfluxes$good_dic[which(HTfluxes$tot_NAs_DIC>10)] <- 0

# giving "bad" markers to half tides of abnormal duration
HTfluxes$good_doc[which(HTfluxes$duration>7.4)] <- 0
HTfluxes$good_dic[which(HTfluxes$duration<5)] <- 0

# giving "bad" markers to half tides where water level data is missing around
# the low tide marker. Missing data could cause the placement of the low tide
# marker to be incorrect, such that the duration of /data included in the tide
# is not accurate
HL_place <- which(!is.na(data$HL))
HL_num <- data$halftide_num[HL_place]
before_HL <- HL_place - 1
after_HL <- HL_place + 1

# determine half tide numbers where water level data is missing before or after
HT_missing_before <- HL_num[which(is.na(data$h[before_HL]))]+1
HT_missing_after <- HL_num[which(is.na(data$h[after_HL]))]

# create vector containing half tide numbers with missing low tide water level
end_mark <- sort(unique(c(HT_missing_before,HT_missing_before+1,
                          HT_missing_after,HT_missing_after+1)))

# set these as "bad" half tides
HTfluxes$good_doc[which(HTfluxes$halftide_num %in% end_mark)] <- 0
HTfluxes$good_dic[which(HTfluxes$halftide_num %in% end_mark)] <- 0

# calculate DIC/DOC fluxes + relevant tidal characteristics for each tide 
Tfluxes <- group_by(data,tide_num)%>%
  summarize(start_datetime = datetime[1],
            total_doc_interp = sum(docinterp)*5*60/4132,
            total_dic_interp = sum(dicinterp)*5*60/4132,
            avg_doc_interp = mean(docinterp),
            avg_dic_interp = mean(dicinterp),
            net_flow = sum(flow_corr_shift)*300,
            high_h=ifelse(length(h[which(HL=="H")])>0,h[which(HL=="H")],NA),
            conc_diff_dic = mean(dic_gm3[which(flow_corr_shift<0)])-
              mean(dic_gm3[which(flow_corr_shift>0)]),
            conc_diff_doc = mean(doc_gm3[which(flow_corr_shift<0)])-
              mean(doc_gm3[which(flow_corr_shift>0)]),
            HT_in = halftide_num[1],HT_out = halftide_num[length(halftide_num)],
            good_doc = 0,good_dic = 0)

# determine if full tidal cycles are "good" based on whether or not they
# are made up of two "good" half tides. Both half tides are good, value is 2.
Tfluxes$good_doc <- HTfluxes$good_doc[which(HTfluxes$halftide_num %in%
                                              Tfluxes$HT_in)]+
  HTfluxes$good_doc[which(HTfluxes$halftide_num %in% Tfluxes$HT_out)]

Tfluxes$good_dic <- HTfluxes$good_dic[which(HTfluxes$halftide_num %in%
                                              Tfluxes$HT_in)]+
  HTfluxes$good_dic[which(HTfluxes$halftide_num %in% Tfluxes$HT_out)]

# select final values used in analysis
Tfluxes_final <- filter(Tfluxes,good_doc==2|good_dic==2,abs(net_flow)<1000)%>%
  mutate(total_dic_interp = ifelse(good_dic == 2,total_dic_interp,NA),
         total_doc_interp = ifelse(good_doc == 2,total_doc_interp,NA),
         avg_dic_interp = ifelse(good_dic == 2,avg_dic_interp,NA),
         avg_doc_interp = ifelse(good_doc == 2,avg_doc_interp,NA))%>%
  select(start_datetime,tide_num, total_dic_interp,total_doc_interp,
         avg_dic_interp,avg_doc_interp,net_flow,high_h,
         conc_diff_dic,conc_diff_doc)

####################Calculating Monthly/Annual Carbon Fluxes#####################

# this function determines the boundary values for a set number of  bins based 
# on instantaneous water flow values. The bins are symmetrical around zero and 
# each bin is sized to hold approximately the same percentage of all flow data
# based on the input variable "bin_per". We used this function to create 20 flow
# bins each containing ~5% of data
makeBinsPer <- function(flows,bins_per){
  
  # create a vector of all instantaneous flow values and sort in ascending order
  flow_vals <- sort(flows)
  
  # separate incoming (positive) and outgoing (negative) flows
  flow_vals_in <- flow_vals[which(flow_vals > 0)]
  flow_vals_out <- flow_vals[which(flow_vals < 0)]
  
  # determine the number of flows in each group
  l_in <- length(flow_vals_in)
  l_out <- length(flow_vals_out)
  
  # determine the number of entries that should go in each bin (separately for 
  # incoming and outgoing) 
  n_per_bin_in <- bins_per*2/100*l_in
  n_per_bin_out <- bins_per*2/100*l_out
  
  # find the minimum and maximum bin values that will be less than lowest/
  # greater than greatest flow values
  min <- floor(min(flow_vals_out,na.rm=TRUE))
  max <- ceiling(max(flow_vals_in,na.rm=TRUE))
  
  # determine the number of bins based on the input % per bin
  n_bins <- 100/bins_per/2-1
  
  # create a vector to store bin values for bins holding incoming flows
  bins_in <- rep(NA,n_bins)
  
  # determine the numeric flow value associated with each bin based on the 
  # average of the incoming and outgoing values for that bin
  for (i in 1:n_bins) {
    
    val_in <- flow_vals_in[n_per_bin_in*i]
    
    val_out <- flow_vals_out[length(flow_vals_out)-n_per_bin_out*i]
    
    bins_in[i] <- (val_in+abs(val_out))/2
    
  }
  
  # create values for outgoing bins as the negatives of the incoming bin values
  bins_out <- bins_in*-1
  
  # create one vector to hold all bin values, including 0
  bins <- sort(c(min,max,bins_in,bins_out,0))
  
  # create a list to store function output
  out <- list()
  out$bins <- bins
  
  # calculate width of each bin and add to output list
  bin_width <- bins[-c(1,length(bins))]
  bin_width <- c(min(flow_vals_out),bin_width,max(flow_vals_in))
  bin_width_prev <- c(NA,bin_width)
  bin_width_prev <- bin_width_prev[-length(bin_width_prev)]
  bin_width <- bin_width-bin_width_prev
  out$bin_width <- bin_width[which(!is.na(bin_width))]
  
  return(out)
  
}

# this function assigns a bin to each entry of instantaneous flow data. inputs 
# are the data to be sorted and a vector of bin values (calculated with 
# makeBinsPer function).
setBins <- function(data,bins){
  
  #upper bounds for values in each bin
  bin_h <- bins[-1]
  
  #lower bounds for values in each bin
  bin_l <- bins[-length(bins)]
  
  # vector to store bin number 
  bin_num <- rep(NA,length(data))
  
  # sequentially find and label values that fall with in the span of values for 
  # each bin
  for (i in 1:length(bin_h)){
    
    counter <- which(data<=bin_h[i]
                     & data>bin_l[i])
    
    bin_num[counter] <- bin_h[i]
    
  }
  
  return(bin_num)
  
}

# this function is used to calculate the average carbon flux rate in each 
# flow bin in order to calculate a weighted average of carbon fluxes based on 
# flow data
binSum <- function(data){
  
  # separately group bin values in each season and year and calculate 
  # flux averages
  bin_calc <- group_by_at(data, c("year","month","flow_bin"))%>%
    summarize(N=n(),doc_flux = mean(docinterp,na.rm=TRUE),
              N_doc = length(which(!is.na(docinterp))),
              dic_flux = mean(dicinterp,na.rm=TRUE),
              N_dic = length(which(!is.na(dicinterp))))%>%
    filter(!is.na(flow_bin))
  
  bin_calc <- merge(bin_calc, bins_tb, by= c("year","month","flow_bin"), 
                    all = TRUE)
  
  # calculate total number of entries in all flow bins for each season/year
  N_sum <- group_by_at(bin_calc,c("year","month"))%>%summarize(N=sum(N))
  
  # create column for total number of entries in summary data frame
  sums <- rep(N_sum$N,rep(20,nrow(N_sum)))
  bin_calc$Nsum <- sums
  
  # calculate the percent of time (/flows) of the total that is contained in
  # each bin, add to data frame
  per_time <- bin_calc$N/sums
  bin_calc$per_time <- per_time
  
  # calculate products of average flux/conc per bin and time per bin
  bin_calc <- mutate(bin_calc,doc_flux_time = doc_flux*per_time,
                     dic_flux_time = dic_flux*per_time)
  
  return(bin_calc)
  
}

# create 5% bins of instantaneous flow values 
bin_vals <- makeBinsPer(data$flow_corr_shift,5)

bins_tb <- dates[rep(seq_len(nrow(dates)), 20), ]%>%arrange_at(c("year","month"))

bins_tb$flow_bin <- rep(bin_vals$bins[-1],nrow(dates))

# assign each instantaneous flow to a bin
data$flow_bin <- setBins(data$flow_corr_shift,bin_vals$bins)

# calculate summary stats for each flow bin
data_bin_sum <- binSum(data)

# determine flow-weighted flux/conc averages by summing the 
# flux*time columns containing all bins for each month/year
month_sum <- group_by_at(data_bin_sum,c('year','month'))%>%
  summarize(N_flows = sum(N), N_doc_flux = sum(N_doc),N_dic_flux = sum(N_dic),
            doc_flux_avg = sum(doc_flux_time)*31536000/4132, #g/m2/yr
            dic_flux_avg = sum(dic_flux_time)*31536000/4132)

# determine the total number of possible entries in each season/year and 
# calculate net flow, mean tidal exchange volume
N_all <- group_by_at(data,c("year","month"))%>%
  summarize(N=n(), flow_m3yr = mean(flow_corr_shift,na.rm=TRUE)*31536000,
            flow_in_avg = mean(flow_corr_shift[which(flow_corr_shift>0)],
                               na.rm=TRUE),
            flow_out_avg = mean(flow_corr_shift[which(flow_corr_shift<0)],
                                na.rm=TRUE),
            flow_in_out_avg = (flow_in_avg + abs(flow_out_avg))/2,
            MTEV = flow_in_out_avg*60*60*6.2,
            N_sal = length(which(!is.na(sal))),
            N_flow_sal = length(which(!is.na(sal)&!is.na(flow_corr_shift)))) # m3

# use mean values where values couldn't be calculated
N_all$flow_m3yr[which(is.na(N_all$flow_m3yr))] <- mean(N_all$flow_m3yr,na.rm=TRUE)

# calculate average MTEV across all months
MTEV_overall <- mean(N_all$MTEV,na.rm=TRUE)

# calculate the C flux data coverage for each season/year 
month_sum <- merge(month_sum,N_all,by=c("year","month"))%>%
  mutate(doc_cover = N_doc_flux/N,dic_cover = N_dic_flux/N,
         flow_cover = N_flows/N,sal_cover = N_sal/N, flow_sal_cover = N_flow_sal/N)

# removing flux values with low coverage 
month_sum$doc_flux_avg[which(month_sum$doc_cover < 0.2)] <- NA
month_sum$dic_flux_avg[which(month_sum$dic_cover < 0.2)] <- NA

##############Calculate Groundwater Contribution to Monthly Fluxes###############

# incorporate monthly averages of meteorological data calculated from 
# WBNERR data
month_met <- read_csv("data/WBNERR_met_averages.csv")
month_sum <- merge(month_sum,month_met,by=c("year","month"),all=TRUE)

# estimated groundwater DIC and DOC concentrations (reported in umol/L and 
# converted to g/m3)
gw_conc_doc <- 76*0.012
gw_conc_dic <- 1081*0.012 

# calculate mean groundwater contribution
month_sum <- mutate(month_sum, gw_m3yr = flow_m3yr-ET_m3yr-precip_m3yr)

## calculate the average of net flow terms for each month across years
mean_gw <- mean(month_sum$gw_m3yr,na.rm=TRUE)

# calculate estimated groundwater contribution to net freshwater discharge, 
# calculate groundwater DIC and DOC fluxes, calculate marsh-derived DIC 
# and DOC fluxes by subtracting groundwater contributions,
# and calculate normalized DIC and DOC fluxes basted on MTEV
month_sum <- mutate(month_sum, gw_fill = ifelse(gw_m3yr>0,mean_gw,gw_m3yr),
                    gw_doc_flux = gw_fill*gw_conc_doc/4132,# g/m2/yr
                    gw_dic_flux = gw_fill*gw_conc_dic/4132,# g/m2/yr
                    marsh_doc_flux = doc_flux_avg - gw_doc_flux,
                    marsh_dic_flux = dic_flux_avg - gw_dic_flux,
                    norm_doc_flux = marsh_doc_flux/MTEV*MTEV_overall,
                    norm_dic_flux = marsh_dic_flux/MTEV*MTEV_overall)

######################Estimate Missing Monthly Fluxes##########################

# linear models of normalized C fluxes vs air temperature 
doc_model <- summary(lm(norm_doc_flux ~ air_temp_C, data = month_sum))
dic_model <- summary(lm(norm_dic_flux ~ air_temp_C, data = month_sum))

# calculate estimated normalized DIC and DOC fluxes based on the above relationship
# with air temperature, un-normalize values to reflect monthly MTEV, then for 
# monthly fluxes, use measured values when available and estimated value when not
month_sum <- mutate(month_sum,marsh_doc_mod_norm = 
                      ifelse(air_temp_C>0,doc_model$coef[1]+
                               doc_model$coef[2]*air_temp_C,NA),
                    marsh_doc_mod = marsh_doc_mod_norm/MTEV_overall*MTEV,
                    marsh_dic_mod_norm = 
                      ifelse(air_temp_C>0,dic_model$coef[1]+
                               dic_model$coef[2]*air_temp_C,NA),
                    marsh_dic_mod = marsh_dic_mod_norm/MTEV_overall*MTEV,
                    marsh_doc_fill = ifelse(!is.na(marsh_doc_flux),
                                            marsh_doc_flux,marsh_doc_mod),
                    marsh_dic_fill = ifelse(!is.na(marsh_dic_flux),
                                            marsh_dic_flux,marsh_dic_mod))%>%
  arrange(year,month)

# select final values 
month_final <- select(month_sum, year, month, marsh_dic_flux, marsh_doc_flux,
                      marsh_dic_fill, marsh_doc_fill,MTEV,air_temp_C)
