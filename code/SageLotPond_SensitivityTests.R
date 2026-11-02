
###################################################################################

source('code/SageLotPond_ExETSCalculation.R')

DOC_er <- list()
DIC_er <- list()
m_list <- group_by_at(data, c("year","month"))%>%summarize()
error_vals <- read_csv("data/errors_sensitivity_test.csv")

############ error vals

suppressMessages(

for (j in 1:nrow(error_vals)){
  
  # current error values for instantaneous parameters
  EpH <- error_vals$EpH[j]
  Etemp <- error_vals$Etemp[j]
  Edic <- error_vals$Edic[j]
  EdocF <- error_vals$EdocF[j]
  EdocnF <- error_vals$EdocnF[j]
  Esal <- error_vals$Esal[j]
  EDO <- error_vals$EDO[j]
  Efdom <- error_vals$Efdom[j]
  Eflow <- error_vals$Eflow[j]
  
  E <- c(EpH,Etemp,Edic,EdocF,EdocnF,Esal,EDO,Efdom,Eflow)
  
  # matrix to store error values
  range <- matrix(NA,nrow = 101,ncol=length(E))
  
  # potential error values (uniform distribution) for instantaneous values
  for (i in 1:length(E)){
    range[,i] <- seq(-E[i],E[i],E[i]/50)
  }
  
  # groundwater error values
  EGW_val <- error_vals$EGW[j]
  
  if(EGW_val == 0){
    EGW <- 0
  }else{
    EGW <- seq(-1*EGW_val,EGW_val,EGW_val/50)
  }
  
  # overland flow error values
  OFE_val <- error_vals$OFE[j]
  
  if(OFE_val == 0){
    OFE <- 0
  }else{
    OFE <- seq(-1*OFE_val,OFE_val,OFE_val/50)
  }
  
  ############## 
  
  # iterations for each set of error values
  n_runs <- 1000
  
  # matrices to store results
  DOC <- matrix(NA,nrow=nrow(m_list),ncol=n_runs)
  DIC <- matrix(NA,nrow=nrow(m_list),ncol=n_runs)
  
  ############# 
  
  # randomly assign error values to instantaneous parameters
  for (ii in 1:n_runs){
    
    data_run <- data
    
    rand_error <- tibble(.rows = nrow(data_run))
    
    for (i in 1:length(E)){
      
      rand_error[,i] <- sample(range[,i],nrow(data_run),replace = TRUE)
      
    }
    
    colnames(rand_error) <- c("EpH","Etemp","Edic","EdocF","EdocnF","Esal","EDO",
                              "Efdom","Eflow")
    
    # add error values to data frame
    data_run <- cbind(data_run,rand_error)
    
    # modify instantaneous measurements to include error values
    data_run <- mutate(data_run, ph = ph+EpH, temp = temp + Etemp,
                      sal = (Esal*sal)+sal,domgL = domgL + EDO, 
                      fdom = (Efdom*fdom)+fdom,flow = (Eflow*flow)+flow)
    
    
    #################################################################################
    
    # calculating DIC concentration (umol/L) from predictive equation, incorporating error value
    data_run <- mutate(data_run, temp_adj = ifelse(temp<=0,NA,temp),
                   DIC = 10^(3.27035)*sal^(0.60693)*ph^(-0.88417)*
                     domgL^(-0.06894)*temp_adj^(-0.03633)+Edic)
    
    # removing calculated DIC concentrations where predictors are out of bounds
    aa <- which(data_run$ph<6.8 | data_run$ph>8.5)
    bb <- which(data_run$sal<19 | data_run$sal>32)
    cc <- which(data_run$temp_adj<2.5 | data_run$temp_adj>30)
    dd <- which(data_run$domgL<0.28 | data_run$domgL>11.85)
    dic_filter <- c(aa,bb,cc,dd)
    data_run$DIC[dic_filter] <- NA
    
    # converting DIC conc to g/m3
    data_run$dic_gm3 <- data_run$DIC*0.012
    
    # calculating DOC concentration (umol/L) from predictive equation, incorporating error value 
    data_run$DOC1 <- 10^(3.02976)*data_run$temp_adj^(0.50412)*data_run$fdom^(0.20955)*
      data_run$ph^(-1.75138)+data_run$EdocF
    
    # removing calculated DOC concentrations where predictors are out of bounds
    ee <- which(data_run$fdom<0.9 | data_run$fdom>80)
    doc_filter_fdom <- c(aa,cc,ee)
    data_run$DOC1[doc_filter_fdom] <- NA
    
    # calculating DOC concentration (umol/L) from predictive equation w/o fDOM, incorporating error value
    data_run$DOC2 <- 10^(4.68979)*data_run$temp_adj^(0.63080)*data_run$ph^(-3.50524)+
      data_run$EdocnF
    
    # removing calculated DOC concentrations where predictors are out of bounds
    doc_filter_no_fdom <- c(aa,cc)
    data_run$DOC2[doc_filter_no_fdom] <- NA
    
    # selecting equation 1 value when available and 2 value otherwise
    data_run$DOC <- ifelse(is.na(data_run$DOC1),data_run$DOC2,data_run$DOC1)
    
    # designating equation used as 1 or 2, and 0 for NA values
    data_run$DOCeq <- ifelse(is.na(data_run$DOC1),2,1)
    data_run$DOCeq[which(is.na(data_run$DOC))] <- 0
    
    # determining tides where both equations 1 and 2 are used
    TF_eq <- group_by(data_run,tide_num)%>%
      summarize(DOCeq = sum(unique(DOCeq)))%>%
      filter(DOCeq==3)
    
    # replacing some equation 1 values with 2 values so that no individual tidal 
    # cycle uses both equations 
    data_run$DOC[which(data_run$tide_num%in%TF_eq$tide_num)] <- 
      data_run$DOC2[which(data_run$tide_num%in%TF_eq$tide_num)]
    
    # updating equation numbers based on replaced data
    data_run$DOCeq[which(data_run$tide_num%in%TF_eq$tide_num)] <- 2
    data_run$DOCeq[which(is.na(data_run$DOC))] <- 0
    
    # converting DOC conc to g/m3
    data_run$doc_gm3 <- data_run$DOC*0.012
    
    ########################Water Flow Calculations and Corrections##################
    
    ###Overland Flow Corrections###
    
    # create column containing overland flow corrected flow values
    data_run$flow_corr <- data_run$flow
    
    # selecting values where quadratic overland flow correction 
    # is applied. flood correction applied to positive flows, ebb correction 
    # applied to negative flows
    flood_excess <- which(data_run$flow>0 & data_run$h>0.25 & data_run$h<=0.5)
    ebb_excess <- which(data_run$flow<0 & data_run$h>0.26 & data_run$h<=0.5)
    
    # applying ebb and flood overland flow corrections, including error value
    data_run$flow_corr[flood_excess] <- (-15.4*data_run$h[flood_excess]^2 + 
                                       20.4*data_run$h[flood_excess]- 3.11)*data_run$flow[flood_excess]
    data_run$flow_corr[flood_excess] <- data_run$flow_corr[flood_excess]+
      (data_run$flow_corr[flood_excess]*sample(OFE,length(flood_excess),replace=TRUE))
    
    data_run$flow_corr[ebb_excess] <- (12.2*data_run$h[ebb_excess]^2 + 
                                     0.275*data_run$h[ebb_excess]+ 0.097)*data_run$flow[ebb_excess]
    data_run$flow_corr[ebb_excess] <- data_run$flow_corr[ebb_excess]+
      (data_run$flow_corr[ebb_excess]*sample(OFE,length(ebb_excess),replace=TRUE))
    
    # selecting values where linear overland flood correction is 
    # applied, flood and ebb use same correction for water level > 0.5 m
    excess_linear <- which(data_run$h>0.5)
    
    # applying linear overland flow correction, including error value
    data_run$flow_corr[excess_linear] <- 
      (data_run$h[excess_linear]*9.170833-1.320587)*data_run$flow[excess_linear]
    data_run$flow_corr[excess_linear] <- data_run$flow_corr[excess_linear]+
      (data_run$flow_corr[excess_linear]*sample(OFE,length(excess_linear),replace=TRUE))
    
    
    ####Water Budget Balance/Net Flow Corrections####
    
    # net flow values calculated from original data
    net_flow_calc_run <- net_flow_calc
    
    # apply error values to calculated net flow
    net_flow_calc_run$flow_m3s <- sample(EGW,nrow(net_flow_calc_run),replace=TRUE)*
      net_flow_calc_run$flow_m3s+net_flow_calc_run$flow_m3s
    
    # create column for flow data with corrected water budget
    data_run$flow_corr_shift <- data_run$flow_corr
    
    # create dataframe to store net flow correction factor
    shift_corr <- tibble(month = rep(NA,nrow(dates)), year = rep(NA,nrow(dates)), 
                         factor = rep(NA,nrow(dates)))
    
    # calculate corrected instantaneous flow values for each month
    for(i in 1:nrow(shift_corr)){
      
      # current month and year
      M <- shift_corr$month[i] <- dates$month[i]
      Y <- shift_corr$year[i] <- dates$year[i]
      
      # creek data_run from current month
      month <- filter(data_run, month == M, year == Y)
      
      # net flow value from salt balance calculation for current month
      water_budget <- filter(net_flow_calc_run, year == Y, month == M)
      
      # calculation parameters for net flow correction factor
      ebb_mean <- mean(month$flow_corr[which(month$flow_corr<0)],na.rm=TRUE)
      flood_mean <- mean(month$flow_corr[which(month$flow_corr>0)],na.rm=TRUE)
      tot_mean <- mean(month$flow_corr,na.rm=TRUE)
      ebb_n <- length(which(month$flow_corr<0))
      flood_n <- length(which(month$flow_corr>0))
      tot_n <- ebb_n+flood_n
      
      # calculate correction factor
      shift_corr$factor[i] <- (water_budget$flow_m3s-(flood_n*flood_mean)/tot_n)*
        (tot_n/(ebb_n*ebb_mean))
      
      # apply correction factor to ebb (outflowing) data_run to determine corrected
      # instantaneous flow values
      ebbs <- which(data_run$flow_corr<0 & data_run$month == M & data_run$year == Y)
      data_run$flow_corr_shift[ebbs] <- data_run$flow_corr[ebbs]*shift_corr$factor[i]
      
    }
    
    #####################Calculating Instantaneous Carbon Fluxes#####################
    
    # instantaneous DIC and DOC fluxes (g/s) calculated as concentration (g/m3)
    # multiplied by flow (m3/s)
    data_run$docflux <- data_run$doc_gm3*data_run$flow_corr_shift
    data_run$dicflux <- data_run$dic_gm3*data_run$flow_corr_shift
    
    ###Interpolation of Missing C Flux data_run###
    
    # create a marker for gaps in DIC and DOC flux data_run
    data_run$doc_gap_check <- ifelse(is.na(data_run$docflux),1,0)
    data_run$dic_gap_check <- ifelse(is.na(data_run$dicflux),1,0)
    
    # create columns that list the length of the gap each value is part of 
    data_run$doc_gap_length <- gapcheck2length(data_run$doc_gap_check)
    data_run$dic_gap_length <- gapcheck2length(data_run$dic_gap_check)
    
    # perform spline interpolation to fill in gaps in flux
    data_run$docinterp <- na.spline(data_run$docflux)
    data_run$dicinterp <- na.spline(data_run$dicflux)
    
    # set any interpolated values where the gap length is longer than 10 entries back to NA
    data_run$docinterp[which(data_run$doc_gap_check == 1 & data_run$doc_gap_length > 10)] <- NA
    data_run$dicinterp[which(data_run$dic_gap_check == 1 & data_run$dic_gap_length > 10)] <- NA
    
    ####################Calculating Monthly/Annual Carbon Fluxes#####################
    
    # create 5% bins of instantaneous flow values 
    bin_vals <- makeBinsPer(data_run$flow_corr_shift,5)
    
    bins_tb <- dates[rep(seq_len(nrow(dates)), 20), ]%>%arrange_at(c("year","month"))
    
    bins_tb$flow_bin <- rep(bin_vals$bins[-1],nrow(dates))
    
    # assign each instantaneous flow to a bin
    data_run$flow_bin <- setBins(data_run$flow_corr_shift,bin_vals$bins)
    
    # calculate summary stats for each flow bin
    data_run_bin_sum <- binSum(data_run)
    
    # determine flow-weighted flux/conc averages by summing the 
    # flux*time columns containing all bins for each month/year
    month_sum_run <- group_by_at(data_run_bin_sum,c('year','month'))%>%
      summarize(N_flows = sum(N), N_doc_flux = sum(N_doc),N_dic_flux = sum(N_dic),
                doc_flux_avg = sum(doc_flux_time)*31536000/4132, #g/m2/yr
                dic_flux_avg = sum(dic_flux_time)*31536000/4132)
    
    #####################################################################################
    
    # add results to dataframe
    DOC[,ii] <- month_sum_run$doc_flux_avg
    DIC[,ii] <- month_sum_run$dic_flux_avg
    
    if(ii%%10==0){print(paste(ii,j,sep = ", "))}
    if(ii%%50==0){cat("\014")}
    
  }
  
  # add results to list
  DOC_er[[j]] <- DOC
  DIC_er[[j]] <- DIC
  
}

)

# month and year values to assign to recalculated data
month_sum <- arrange_at(month_sum,c("year","month"))

all_doc_sd <- tibble(year=month_sum$year, month = month_sum$month,
                     doc_cover = month_sum$doc_cover)
all_dic_sd <- tibble(year=month_sum$year, month = month_sum$month,
                     dic_cover = month_sum$dic_cover)

for(k in 1:length(DOC_er)){
  
  # calculate stdev for monte carlo results
  sy_doc_sd <- cbind(apply(DOC_er[[k]],1,sd))
  sy_dic_sd <- cbind(apply(DIC_er[[k]],1,sd))
  
  all_doc_sd <- cbind(all_doc_sd,sy_doc_sd)
  all_dic_sd <- cbind(all_dic_sd,sy_dic_sd)
  
}

# colnames for calculated measurement uncertainty values in different error groups
doc_colnames <- c( "year","month","cover", "all", "all_conc",
                   "all_flow","pH","temp","C_eq","sal","fdom",
                   "flow","of","gw1")
dic_colnames <- c( "year","month","cover", "all", "all_conc",
                   "all_flow","pH","temp","C_eq","sal","DO",
                   "flow","of","gw1")

colnames(all_doc_sd) <- doc_colnames
colnames(all_dic_sd) <- dic_colnames
