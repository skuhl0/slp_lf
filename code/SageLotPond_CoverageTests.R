
source('code/SageLotPond_ExETSCalculation.R')

####### logging interval data coverage test

# different timestep options, will select every Xth data point for logging 
# intervals of 10, 15, 30, and 60 minutes 
coverage <- c(2,3,6,12)

# lists to store results
DOC_ts <- list()
DIC_ts <- list()

#################################################################################

# outer loop repeats over each option of logging interval
for(j in 1:length(coverage)){
  
  # matrices to store results
  DOC <- matrix(NA,nrow=nrow(dates),ncol=coverage[j])
  DIC <- matrix(NA,nrow=nrow(dates),ncol=coverage[j])
  
  # inner loop repeats over unique subsets of data with given logging interval
  for(ii in 1:coverage[j]){
    
    # select data for current test
    data_run <- data[seq(ii,nrow(data),coverage[j]),]
    
    # calculate corrected instantaneous flow based on original freshwater discharge
    for(i in 1:nrow(dates)){
      
      # current month and year
      M <- dates$month[i]
      Y <- dates$year[i]
      
      # creek data from current month
      month <- filter(data_run, month == M, year == Y)
      
      # net flow value from salt balance calculation for current month
      water_budget <- filter(net_flow_calc, year == Y, month == M)
      
      # calculation parameters for net flow correction factor
      ebb_mean <- mean(month$flow_corr[which(month$flow_corr<0)],na.rm=TRUE)
      flood_mean <- mean(month$flow_corr[which(month$flow_corr>0)],na.rm=TRUE)
      tot_mean <- mean(month$flow_corr,na.rm=TRUE)
      ebb_n <- length(which(month$flow_corr<0))
      flood_n <- length(which(month$flow_corr>0))
      tot_n <- ebb_n+flood_n
      
      # calculate correction factor
      factor <- (water_budget$flow_m3s-(flood_n*flood_mean)/tot_n)*
        (tot_n/(ebb_n*ebb_mean))
      
      # apply correction factor to ebb (outflowing) data to determine corrected
      # instantaneous flow values
      ebbs <- which(data_run$flow_corr<0 & data_run$month == M & 
                      data_run$year == Y)
      data_run$flow_corr_shift[ebbs] <- data_run$flow_corr[ebbs]*factor
      
    }
    
    # instantaneous DIC and DOC fluxes (g/s) calculated as concentration (g/m3)
    # multiplied by flow (m3/s)
    data_run$docflux <- data_run$doc_gm3*data_run$flow_corr_shift
    data_run$dicflux <- data_run$dic_gm3*data_run$flow_corr_shift
    
    ###Interpolation of Missing C Flux Data###
    
    # create a marker for gaps in DIC and DOC flux data
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
      summarize(N_flows = sum(N), N_doc_flux = sum(N_doc),
                N_dic_flux = sum(N_dic),
                doc_flux_avg = sum(doc_flux_time)*31536000/4132, #g/m2/yr
                dic_flux_avg = sum(dic_flux_time)*31536000/4132)
    
    # arrange data
    month_sum_run <- arrange_at(month_sum_run,c("year","month"))
    
    # store results from current test
    DOC[,ii] <- month_sum_run$doc_flux_avg
    DIC[,ii] <- month_sum_run$dic_flux_avg
    
  }
  
  # store results from current test
  DOC_ts[[j]] <- DOC
  DIC_ts[[j]] <- DIC
  
}

# original data to compare simulations to
doc_ts_m <- tibble(year=month_sum$year, month = month_sum$month,
                   MY = paste(month,year,sep = "-"),
                   avg = month_sum$doc_flux_avg,
                   cover = month_sum$doc_cover)
dic_ts_m <- tibble(year=month_sum$year, month = month_sum$month,
                   MY = paste(month,year,sep = "-"),
                   avg = month_sum$dic_flux_avg,
                   cover = month_sum$dic_cover)

# months with high data coverage to use
high_cov <- doc_ts_m$MY[which(doc_ts_m$cover>0.8&dic_ts_m$cover>0.8)]

# table to store restructured results
doc_long_ts <- tibble()
dic_long_ts <- tibble()

# names for logging intervals tested ("m10" = 10 minutes, etc.)
log_ints <- c("m10","m15","m30","m60")

# extract simulation results and restructure data
for(k in 1:length(DOC_ts)){
  
  # match original measurements to simulation results
  doc_curr <- cbind(doc_ts_m,DOC_ts[[k]])
  dic_curr <- cbind(dic_ts_m,DIC_ts[[k]])
  
  # pivot data 
  doc_curr_long <- pivot_longer(doc_curr, cols = c(6:(ncol(DOC_ts[[k]])+5)),
                                names_to = "run", values_to = "C_recalc")
  
  dic_curr_long <- pivot_longer(dic_curr, cols = c(6:(ncol(DIC_ts[[k]])+5)),
                                names_to = "run", values_to = "C_recalc")
  
  # select only high coverage months
  doc_long_sel <- filter(doc_curr_long,MY%in%high_cov)%>%
    mutate(c_type = "doc",log_int = log_ints[k])
  
  dic_long_sel <- filter(dic_curr_long,MY%in%high_cov)%>%
    mutate(c_type = "dic",log_int = log_ints[k])
  
  # add current group to existing data table
  doc_long_ts <- rbind(doc_long_ts,doc_long_sel)
  
  dic_long_ts <- rbind(dic_long_ts,dic_long_sel)
  
}

# table reporting DIC and DOC logging interval test results
ts_long <- rbind(doc_long_ts,dic_long_ts)

###################################################################################

### random missing data coverage test 

# options for percentage of missing data applied
per_missing <- seq(0.1,0.9,0.1)

# lists to store results
DOC_co <- list()
DIC_co <- list()

# outer loop repeats over possible percent missing data values
for (j in 1:length(per_missing)){
  
  # number of iterations of random data selection
  n_runs <- 100
  
  # matrices to store results
  DOC <- matrix(NA,nrow=nrow(dates),ncol=n_runs)
  DIC <- matrix(NA,nrow=nrow(dates),ncol=n_runs)
  
  # inner loop repeats over prescribed number of iterations for random selection
  for(ii in 1:n_runs){
    
    # select data subset
    data_run <- data[-sample(nrow(data),nrow(data)*per_missing[j]),]
    
    # calculate corrected instantaneous flow values for each month
    for(i in 1:nrow(dates)){
      
      # current month and year
      M <- dates$month[i]
      Y <- dates$year[i]
      
      # creek data from current month
      month <- filter(data_run, month == M, year == Y)
      
      # net flow value from salt balance calculation for current month
      water_budget <- filter(net_flow_calc, year == Y, month == M)
      
      # calculation parameters for net flow correction factor
      ebb_mean <- mean(month$flow_corr[which(month$flow_corr<0)],na.rm=TRUE)
      flood_mean <- mean(month$flow_corr[which(month$flow_corr>0)],na.rm=TRUE)
      tot_mean <- mean(month$flow_corr,na.rm=TRUE)
      ebb_n <- length(which(month$flow_corr<0))
      flood_n <- length(which(month$flow_corr>0))
      tot_n <- ebb_n+flood_n
      
      # calculate correction factor
      factor <- (water_budget$flow_m3s-(flood_n*flood_mean)/tot_n)*
        (tot_n/(ebb_n*ebb_mean))
      
      # apply correction factor to ebb (outflowing) data to determine corrected
      # instantaneous flow values
      ebbs <- which(data_run$flow_corr<0 & data_run$month == M & 
                      data_run$year == Y)
      data_run$flow_corr_shift[ebbs] <- data_run$flow_corr[ebbs]*factor
      
    }
    
    # instantaneous DIC and DOC fluxes (g/s) calculated as concentration (g/m3)
    # multiplied by flow (m3/s)
    data_run$docflux <- data_run$doc_gm3*data_run$flow_corr_shift
    data_run$dicflux <- data_run$dic_gm3*data_run$flow_corr_shift
    
    ###Interpolation of Missing C Flux Data###
    
    # create a marker for gaps in DIC and DOC flux data
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
      summarize(N_flows = sum(N), N_doc_flux = sum(N_doc),
                N_dic_flux = sum(N_dic),
                doc_flux_avg = sum(doc_flux_time)*31536000/4132, #g/m2/yr
                dic_flux_avg = sum(dic_flux_time)*31536000/4132)
    
    # arrange data
    month_sum_run <- arrange_at(month_sum_run,c("year","month"))
    
    # store results from current test
    DOC[,ii] <- month_sum_run$doc_flux_avg
    DIC[,ii] <- month_sum_run$dic_flux_avg
    
  }
  
  # store results from current test
  DOC_co[[j]] <- DOC
  DIC_co[[j]] <- DIC
  
}

# original data to compare simulations to
doc_co_m <- tibble(year=month_sum$year, month = month_sum$month,
                   MY = paste(month,year,sep = "-"),
                   avg = month_sum$doc_flux_avg,
                   cover = month_sum$doc_cover)
dic_co_m <- tibble(year=month_sum$year, month = month_sum$month,
                   MY = paste(month,year,sep = "-"),
                   avg = month_sum$dic_flux_avg,
                   cover = month_sum$dic_cover)

# months with high data coverage to use
high_cov <- doc_co_m$MY[which(doc_co_m$cover>0.8&dic_co_m$cover>0.8)]

# table to store restructured results
doc_long_co <- tibble()
dic_long_co <- tibble()

# names for random data values tested ("m10" = 10% of data missing)
log_ints <- c("m10","m20","m30","m40","m50",
              "m60","m70","m80","m90")

# extract simulation results and restructure data
for(k in 1:length(DOC_co)){
  
  # match original measurements to simulation results
  doc_curr <- cbind(doc_co_m,DOC_co[[k]])
  dic_curr <- cbind(dic_co_m,DIC_co[[k]])
  
  # pivot data 
  doc_curr_long <- pivot_longer(doc_curr, cols = c(6:(ncol(DOC_co[[k]])+5)),
                                names_to = "run", values_to = "C_recalc")
  
  dic_curr_long <- pivot_longer(dic_curr, cols = c(6:(ncol(DIC_co[[k]])+5)),
                                names_to = "run", values_to = "C_recalc")
  
  # select only high coverage months
  doc_long_sel <- filter(doc_curr_long,MY%in%high_cov)%>%
    mutate(c_type = "doc",log_int = log_ints[k])
  
  dic_long_sel <- filter(dic_curr_long,MY%in%high_cov)%>%
    mutate(c_type = "dic",log_int = log_ints[k])
  
  # add current group to existing data table
  doc_long_co <- rbind(doc_long_co,doc_long_sel)
  
  dic_long_co <- rbind(dic_long_co,dic_long_sel)
  
}

# table reporting DIC and DOC coverage test results
co_long <- rbind(doc_long_co,dic_long_co)

######################################################################################

####### decreased deployment length data coverage tests

# determine high data coverage months
high_cover <- filter(month_sum, doc_cover > 0.8, dic_cover > 0.8)%>%
  select(year,month)

# select data from only high coverage months
data_high_cover <- merge(data,high_cover, 
                        by=c("year","month"), all=FALSE)%>%
  arrange(datetime)

# select tidal cycles completely within these months
data_hc_tides <- group_by_at(data_high_cover, c("year","month"))%>%
  summarise(tide_start = min(tide_num)+1,tide_end=max(tide_num)-1)

# different deployment length options in number of tidal cycles
n_tides <- c(2,7,14,21,28,35,42,49,56)

# lists to store results
DOC_dep <- list()
DIC_dep <- list()

# outer loop repeats over number of tidal cycles in deployment
for (j in 1:length(n_tides)){
  
  # matrices to store results
  DOC <- matrix(NA,nrow=nrow(high_cover),ncol=59-n_tides[j])
  DIC <- matrix(NA,nrow=nrow(high_cover),ncol=59-n_tides[j])

  # inner loop repeats over unique subsets of tidal cycles
  for(ii in 0:(58-n_tides[j])){
    
    # determine first and last tidal cycle in deployment for each month
    a <- data_hc_tides$tide_start + ii
    b <- a + n_tides[j]-1
    
    # list of selected tidal cycles for each month
    select_tides_list <- lapply(a, function(x) x:b[b>x][1])
    
    # remove any tidal cycles from lists not completely within given month
    for(k in 1:length(select_tides_list)){
      
      oob <- which(select_tides_list[[k]]>data_hc_tides$tide_end[k])
      
      if(length(oob)>0){select_tides_list[[k]] <-select_tides_list[[k]][-oob]}
      
    }
    
    # restructure
    select_tides <- unlist(select_tides_list)
    
    # select data
    data_run <- filter(data_high_cover,tide_num %in% select_tides)
    
    # calculate corrected instantaneous flow values for each month
    for(i in 1:nrow(high_cover)){
      
      # current month and year
      M <- high_cover$month[i]
      Y <- high_cover$year[i]
      
      # creek data from current month
      month <- filter(data_run, month == M, year == Y)
      
      # net flow value from salt balance calculation for current month
      water_budget <- filter(net_flow_calc, year == Y, month == M)
      
      # calculation parameters for net flow correction factor
      ebb_mean <- mean(month$flow_corr[which(month$flow_corr<0)],na.rm=TRUE)
      flood_mean <- mean(month$flow_corr[which(month$flow_corr>0)],na.rm=TRUE)
      tot_mean <- mean(month$flow_corr,na.rm=TRUE)
      ebb_n <- length(which(month$flow_corr<0))
      flood_n <- length(which(month$flow_corr>0))
      tot_n <- ebb_n+flood_n
      
      # calculate correction factor
      factor <- (water_budget$flow_m3s-(flood_n*flood_mean)/tot_n)*
        (tot_n/(ebb_n*ebb_mean))
      
      # apply correction factor to ebb (outflowing) data to determine corrected
      # instantaneous flow values
      ebbs <- which(data_run$flow_corr<0 & data_run$month == M & 
                      data_run$year == Y)
      data_run$flow_corr_shift[ebbs] <- data_run$flow_corr[ebbs]*factor
      
    }
    
    # instantaneous DIC and DOC fluxes (g/s) calculated as concentration (g/m3)
    # multiplied by flow (m3/s)
    data_run$docflux <- data_run$doc_gm3*data_run$flow_corr_shift
    data_run$dicflux <- data_run$dic_gm3*data_run$flow_corr_shift
    
    ###Interpolation of Missing C Flux Data###
    
    # create a marker for gaps in DIC and DOC flux data
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
      summarize(N_flows = sum(N), N_doc_flux = sum(N_doc),
                N_dic_flux = sum(N_dic),
                doc_flux_avg = sum(doc_flux_time)*31536000/4132, #g/m2/yr
                dic_flux_avg = sum(dic_flux_time)*31536000/4132)
    
    # check number of tidal cycles included in each test (might be fewer if 
    # tides not completely within the month were removed)
    month_n_tides <- group_by_at(data_run,c('year','month'))%>%
      summarize(n_tides = length(unique(tide_num)),
                tn_start = unique(tide_num)[1])
    month_n_tides <- merge(month_n_tides,high_cover,
                           by=c("year","month"), all=TRUE)
    month_n_tides <- arrange_at(month_n_tides,c("year","month"))
    
    # remove runs were there were not enough tidal cycles
    month_sum_run <- merge(month_sum_run,high_cover, by=c("year","month"), all.y=TRUE)
    month_sum_run <- arrange_at(month_sum_run,c("year","month"))
    month_sum_run[which(month_n_tides$n_tides<n_tides[j]),] <- NA
    
    # store results from current test
    DOC[,ii+1] <- month_sum_run$doc_flux_avg
    DIC[,ii+1] <- month_sum_run$dic_flux_avg

  }
  
  #store results from current test
  DOC_dep[[j]] <- DOC
  DIC_dep[[j]] <- DIC
  
}

# original data to compare simulations to
month_sum_sel <- merge(month_sum,high_cover, by=c("year","month"), all=FALSE)
month_sum_sel <- arrange_at(month_sum_sel,c("year","month"))

# months with high data coverage to use
doc_dep_m <- tibble(year=month_sum_sel$year, month = month_sum_sel$month,
                    MY = paste(month,year,sep = "-"),
                    avg = month_sum_sel$doc_flux_avg,
                    cover = month_sum_sel$doc_cover)
dic_dep_m <- tibble(year=month_sum_sel$year, month = month_sum_sel$month,
                    MY = paste(month,year,sep = "-"),
                    avg = month_sum_sel$dic_flux_avg,
                    cover = month_sum_sel$dic_cover)

# months with high data coverage to use
high_cov <- doc_dep_m$MY[which(doc_dep_m$cover>0.8&dic_dep_m$cover>0.8)]

# table to store restructured results
doc_long_co <- tibble()
dic_long_co <- tibble()

# names for deployment lengths tested ("t02" = 2 tidal cycles)
log_ints <- c("t02","t07","t14","t21","t28","t35","t42","t49","t56")

# extract simulation results and restructure data
for(k in 1:length(DOC_dep)){
  
  # match original measurements to simulation results
  doc_curr <- cbind(doc_dep_m,DOC_dep[[k]])
  dic_curr <- cbind(dic_dep_m,DIC_dep[[k]])
  
  # pivot data 
  doc_curr_long <- pivot_longer(doc_curr, cols = c(6:(ncol(DOC_dep[[k]])+5)),
                                names_to = "run", values_to = "C_recalc")
  
  dic_curr_long <- pivot_longer(dic_curr, cols = c(6:(ncol(DIC_dep[[k]])+5)),
                                names_to = "run", values_to = "C_recalc")
  
  # select only high coverage months
  doc_long_sel <- filter(doc_curr_long,MY%in%high_cov)%>%
    mutate(c_type = "doc",log_int = log_ints[k])
  
  dic_long_sel <- filter(dic_curr_long,MY%in%high_cov)%>%
    mutate(c_type = "dic",log_int = log_ints[k])

  # add current group to existing data table
  doc_long_co <- rbind(doc_long_co,doc_long_sel)

  dic_long_co <- rbind(dic_long_co,dic_long_sel)
  
}

# table reporting DIC and DOC deployment length test results
dep_long <- rbind(doc_long_co,dic_long_co)

###############################################################################
#############################################

####### spring-neap sampling tests

# determine high data coverage months
high_cover <- filter(month_sum, doc_cover > 0.8, dic_cover > 0.8)%>%
  select(year,month)

# select data from only high coverage months
data_high_cover <- merge(data,high_cover, 
                         by=c("year","month"), all=FALSE)%>%
  arrange(datetime)

# select tidal cycles completely within these months
data_hc_tides <- group_by_at(data_high_cover, c("year","month"))%>%
  summarise(tide_start = min(tide_num)+1,tide_end=max(tide_num)-1)

# different deployment length options in number of tidal cycles
n_tides <- c(1,2,3,4)

# lists to store results
DOC_bes <- list()
DIC_bes <- list()
tn_start_bes <- list()

# outer loop repeats over number of tidal cycles in deployment
for (j in 1:length(n_tides)){
  
  # matrices to store results
  DOC <- matrix(NA,nrow=nrow(high_cover),ncol=59-n_tides[j])
  DIC <- matrix(NA,nrow=nrow(high_cover),ncol=59-n_tides[j])
  tn_start <- matrix(NA,nrow=nrow(high_cover),ncol=59-n_tides[j])
  
  # inner loop repeats over unique subsets of tidal cycles
  for(ii in 0:(58-n_tides[j])){
    
    # determine first tidal cycle in deployment for each month
    a <- data_hc_tides$tide_start + ii
    select_tides <- a
    
    # if more than one tidal cycle in deployment
    if(n_tides[j]>1){
      
      # determine last tidal cycle in deployment for each month
      b <- a + n_tides[j]-1
      
      # list of selected tidal cycles for each month
      select_tides_list <- lapply(a, function(x) x:b[b>x][1])
      
      # remove any tidal cycles from lists not completely within given month
      for(k in 1:length(select_tides_list)){
        
        oob <- which(select_tides_list[[k]]>data_hc_tides$tide_end[k])
        
        if(length(oob)>0){select_tides_list[[k]] <-select_tides_list[[k]][-oob]}
        
      }
      
      # restructure
      select_tides <- unlist(select_tides_list)
      
    }
    
    # select data
    data_run <- filter(data_high_cover,tide_num %in% select_tides)
    
    # decrease the concentration sampling frequency to 30 min
    conc_exclude <- c(1:nrow(data_run))[!(c(1:nrow(data_run))%in%
                                            seq(1,nrow(data_run),6))]
    data_run$doc_gm3[conc_exclude] <- NA
    data_run$dic_gm3[conc_exclude] <- NA
    
    # calculate corrected instantaneous flow values for each month
    for(i in 1:nrow(high_cover)){
      
      # current month and year
      M <- high_cover$month[i]
      Y <- high_cover$year[i]
      
      # creek data from current month
      month <- filter(data_run, month == M, year == Y)
      
      # net flow value from salt balance calculation for current month
      water_budget <- filter(net_flow_calc, year == Y, month == M)
      
      # calculation parameters for net flow correction factor
      ebb_mean <- mean(month$flow_corr[which(month$flow_corr<0)],na.rm=TRUE)
      flood_mean <- mean(month$flow_corr[which(month$flow_corr>0)],na.rm=TRUE)
      tot_mean <- mean(month$flow_corr,na.rm=TRUE)
      ebb_n <- length(which(month$flow_corr<0))
      flood_n <- length(which(month$flow_corr>0))
      tot_n <- ebb_n+flood_n
      
      # calculate correction factor
      factor <- (0-(flood_n*flood_mean)/tot_n)*
        (tot_n/(ebb_n*ebb_mean))
      
      # apply correction factor to ebb (outflowing) data to determine corrected
      # instantaneous flow values
      ebbs <- which(data_run$flow_corr<0 & data_run$month == M & 
                      data_run$year == Y)
      data_run$flow_corr_shift[ebbs] <- data_run$flow_corr[ebbs]*factor
      
    }
    
    # interpolation of concentration
    data_run$doc_conc_interp <- na.approx(data_run$doc_gm3,na.rm=FALSE)
    data_run$dic_conc_interp <- na.approx(data_run$dic_gm3,na.rm=FALSE)
    
    # instantaneous DIC and DOC fluxes (g/s) calculated as concentration (g/m3)
    # multiplied by flow (m3/s)
    data_run$docflux <- data_run$doc_conc_interp*data_run$flow_corr_shift
    data_run$dicflux <- data_run$dic_conc_interp*data_run$flow_corr_shift
    
    ###Interpolation of Missing C Flux Data###
    
    # create a marker for gaps in DIC and DOC flux data
    data_run$doc_gap_check <- ifelse(is.na(data_run$docflux),1,0)
    data_run$dic_gap_check <- ifelse(is.na(data_run$dicflux),1,0)
    
    # create columns that list the length of the gap each value is part of 
    data_run$doc_gap_length <- gapcheck2length(data_run$doc_gap_check)
    data_run$dic_gap_length <- gapcheck2length(data_run$dic_gap_check)
    
    # perform spline interpolation to fill in gaps in flux
    data_run$docinterp <- na.spline(data_run$docflux,na.rm=FALSE)
    data_run$dicinterp <- na.spline(data_run$dicflux,na.rm=FALSE)
    
    # set any interpolated values where the gap length is longer than 10 entries back to NA
    data_run$docinterp[which(data_run$doc_gap_check == 1 & data_run$doc_gap_length > 10)] <- NA
    data_run$dicinterp[which(data_run$dic_gap_check == 1 & data_run$dic_gap_length > 10)] <- NA
    
    # create 5% bins of instantaneous flow values
    bin_vals <- makeBinsPer(data_run$flow_corr_shift,5)

    bins_tb <- dates[rep(seq_len(nrow(dates)), 20), ]%>%arrange_at(c("year","month"))

    bins_tb$flow_bin <- rep(bin_vals$bins[-1],nrow(dates))

    # assign each instantaneous flow to a bin
    data_run$flow_bin <- setBins(data_run$flow_corr_shift,bin_vals$bins)

    # bin calc for reduced data
    binSumRe <- function(data_run){
      
      # separately group bin values in each season and year and calculate 
      # flux averages
      bin_calc <- group_by_at(data_run, c("year","month","flow_bin"))%>%
        summarize(N=n(),doc_flux = mean(docinterp,na.rm=TRUE),
                  N_doc = length(which(!is.na(docinterp))),
                  dic_flux = mean(dicinterp,na.rm=TRUE),
                  N_dic = length(which(!is.na(dicinterp))))%>%
        filter(!is.na(flow_bin))
      
      bin_calc <- merge(bin_calc, bins_tb, by= c("year","month","flow_bin"), 
                        all = TRUE)
      
      bin_calc$N[which(is.na(bin_calc$N))] <- 0
      bin_calc$N_dic[which(is.na(bin_calc$N_dic))] <- 0
      bin_calc$N_doc[which(is.na(bin_calc$N_doc))] <- 0
      
      # calculate total number of entries in all flow bins for each season/year
      N_sum <- group_by_at(bin_calc,c("year","month"))%>%summarize(N=sum(N))
      
      # create column for total number of entries in summary data_run frame
      sums <- rep(N_sum$N,rep(20,nrow(N_sum)))
      bin_calc$Nsum <- sums
      
      # calculate the percent of time (/flows) of the total that is contained in
      # each bin, add to data_run frame
      per_time <- bin_calc$N/sums
      bin_calc$per_time <- per_time
      
      # calculate products of average flux/conc per bin and time per bin
      bin_calc <- mutate(bin_calc,doc_flux_time = doc_flux*per_time,
                         dic_flux_time = dic_flux*per_time)
      
      return(bin_calc)
      
    }
    
    # calculate summary stats for each flow bin
    data_run_bin_sum <- binSumRe(data_run)

    # determine flow-weighted flux/conc averages by summing the
    # flux*time columns containing all bins for each month/year
    month_sum_run <- group_by_at(data_run_bin_sum,c('year','month'))%>%
      summarize(N_flows = sum(N), N_doc_flux = sum(N_doc),
                N_dic_flux = sum(N_dic),
                doc_flux_avg = sum(doc_flux_time,na.rm=TRUE)*31536000/4132, #g/m2/yr
                dic_flux_avg = sum(dic_flux_time,na.rm=TRUE)*31536000/4132)
    
    # check number of tidal cycles included in each test (might be fewer if 
    # tides not completely within the month were removed)
    month_n_tides <- group_by_at(data_run,c('year','month'))%>%
      summarize(n_tides = length(unique(tide_num)),
                tn_start = unique(tide_num)[1])
    month_n_tides <- merge(month_n_tides,high_cover,
                               by=c("year","month"), all=TRUE)
    month_n_tides <- arrange_at(month_n_tides,c("year","month"))

    # remove runs were there were not enough tidal cycles    
    month_sum_run <- merge(month_sum_run,high_cover, by=c("year","month"), 
                           all.y=TRUE)
    month_sum_run <- arrange_at(month_sum_run,c("year","month"))
    month_sum_run[which(month_n_tides$n_tides<n_tides[j]),] <- NA
    
    
    # store results from current test
    DOC[,ii+1] <- month_sum_run$doc_flux_avg
    DIC[,ii+1] <- month_sum_run$dic_flux_avg
    tn_start[,ii+1] <- month_n_tides$tn_start
    
  }
  
  # store results from current test
  DOC_bes[[j]] <- DOC
  DIC_bes[[j]] <- DIC
  tn_start_bes[[j]] <- tn_start
  
}

# original data to compare simulations to
month_sum_sel <- merge(month_sum,high_cover, by=c("year","month"), all=FALSE)
month_sum_sel <- arrange_at(month_sum_sel,c("year","month"))

# months with high data coverage to use
doc_bes_m <- tibble(year=month_sum_sel$year, month = month_sum_sel$month,
                    MY = paste(month,year,sep = "-"),
                    avg = month_sum_sel$doc_flux_avg,
                    cover = month_sum_sel$doc_cover)
dic_bes_m <- tibble(year=month_sum_sel$year, month = month_sum_sel$month,
                    MY = paste(month,year,sep = "-"),
                    avg = month_sum_sel$dic_flux_avg,
                    cover = month_sum_sel$dic_cover)

# months with high data coverage to use
high_cov <- doc_bes_m$MY[which(doc_bes_m$cover>0.8&dic_bes_m$cover>0.8)]

# table to store restructured results
doc_long_bes <- tibble()
dic_long_bes <- tibble()

# names for deployment lengths tested ("t02" = 2 tidal cycles)
log_ints <- c("t01","t02","t03","t04")

# extract simulation results and restructure data
for(k in 1:length(DOC_bes)){
  
  # match original measurements to simulation results
  doc_curr <- cbind(doc_bes_m,DOC_bes[[k]])
  dic_curr <- cbind(dic_bes_m,DIC_bes[[k]])
  tn_curr <- tn_start_bes[[k]]
  
  # pivot data
  doc_curr_long <- pivot_longer(doc_curr, cols = c(6:(ncol(DOC_bes[[k]])+5)),
                                names_to = "run", values_to = "C_recalc")
  
  dic_curr_long <- pivot_longer(dic_curr, cols = c(6:(ncol(DIC_bes[[k]])+5)),
                                names_to = "run", values_to = "C_recalc")
  
  # add tide number to results to select from later
  tn_start_long <- as.vector(t(tn_start_bes[[k]]))
  doc_curr_long$tn_start <- tn_start_long
  dic_curr_long$tn_start <- tn_start_long
  
  # select only high coverage months
  doc_long_sel <- filter(doc_curr_long,MY%in%high_cov)%>%
    mutate(c_type = "doc",log_int = log_ints[k])
  
  dic_long_sel <- filter(dic_curr_long,MY%in%high_cov)%>%
    mutate(c_type = "dic",log_int = log_ints[k])
  
  # add current group to existing data table
  doc_long_bes <- rbind(doc_long_bes,doc_long_sel)
  
  dic_long_bes <- rbind(dic_long_bes,dic_long_sel)
  
}

# report DIC and DOC results together in one table
bes_long <- rbind(doc_long_bes,dic_long_bes)

# calculate subsequent high tide water levels in original tidal cycle data
Tfluxes <- mutate(Tfluxes, nxt_hh_2 = c(high_h[-1],NA),
                  nxt_hh_3 = c(high_h[-c(1,2)],NA,NA),
                  nxt_hh_4 = c(high_h[-c(1:3)],NA,NA,NA))

# select data from high coverage months
Tfluxes_high_cover <- merge(Tfluxes%>%mutate(month = month(start_datetime),
                                             year = year(start_datetime)),
                            high_cover, 
                            by=c("year","month"), all=FALSE)%>%
  arrange(tide_num)%>%mutate(MY = paste(month,year,sep="-"))

# merge high tide water level data with simulation results
bes_long_hh <- merge(bes_long,Tfluxes_high_cover%>%select(tide_num,high_h,
                                                          nxt_hh_2,nxt_hh_3,
                                                          nxt_hh_4),
                     by.x = "tn_start", by.y = "tide_num",
                     all.x = TRUE)

# calculate bounds for upper and lower third of tidal cycles
mid_tides <- group_by_at(Tfluxes_high_cover, c("year","month"))%>%
  summarize(hh_13 = quantile(high_h, .33,na.rm=TRUE),
            hh_23 = quantile(high_h, .66,na.rm=TRUE))

# add bounds to simulation results
bes_long_hh <- merge(bes_long_hh, mid_tides, by = c("year","month"),
                     all = TRUE)

# names for high tides to average over
hh_names <- c("high_h","hh_avg2","hh_avg3","hh_avg4")

# table to store DIC and DOC spring-neap sampling results
bes_resample <- tibble()

# loop over sampling length options
for(k in 1:length(hh_names)){
  
  # which sampling length
  int <- log_ints[k]
  sel_hh <- hh_names[k]
  
  # select data for given sampling length
  bes_doc_hh <- filter(bes_long_hh,c_type == "doc",log_int == int)%>%
    mutate(c_diff = (-1*avg-(-1*C_recalc)),
           hh_avg2 = (high_h+nxt_hh_2)/2,hh_avg3 = (high_h+nxt_hh_2+nxt_hh_3)/3,
           hh_avg4 = (high_h+nxt_hh_2+nxt_hh_3+nxt_hh_4)/4)
  
  bes_dic_hh <- filter(bes_long_hh,c_type == "dic",log_int == int)%>%
    mutate(c_diff = (-1*avg-(-1*C_recalc)),
           hh_avg2 = (high_h+nxt_hh_2)/2,hh_avg3 = (high_h+nxt_hh_2+nxt_hh_3)/3,
           hh_avg4 = (high_h+nxt_hh_2+nxt_hh_3+nxt_hh_4)/4)
  
  # choose high tide averaging group
  bes_doc_hh$sel_h <- bes_doc_hh[,which(colnames(bes_doc_hh)==sel_hh)]
  
  bes_dic_hh$sel_h <- bes_dic_hh[,which(colnames(bes_dic_hh)==sel_hh)]
  
  # assign label for tide height
  bes_doc_hh <- mutate(bes_doc_hh,qt_class = ifelse(sel_h<hh_13,1,
                                                    ifelse(sel_h<hh_23,2,3)))
  
  bes_dic_hh <- mutate(bes_dic_hh,qt_class = ifelse(sel_h<hh_13,1,
                                                    ifelse(sel_h<hh_23,2,3)))
  
  # filter NAs
  bes_doc_nNA <- filter(bes_doc_hh,!is.na(C_recalc))
  
  bes_dic_nNA <- filter(bes_dic_hh,!is.na(C_recalc))
  
  for(i in 1:100){
    
    ############## calculations for 1 spring-neap sampling
    
    # sample a random set of tidal cycles and select spring and neap
    sample <- group_by_at(bes_doc_nNA,c("year","month","qt_class"))%>%
      slice_sample(n=1)%>%filter(qt_class %in% c(1,3))%>%ungroup()%>%
      group_by_at(c("year","month"))%>%
      summarize(C_recalc = mean(C_recalc),avg = avg[1])%>%
      mutate(c_diff = (-1*avg-(-1*C_recalc))/(-1*avg),log_int = int,
             scheme = paste(int,"hl1",sep="-"),
             c_type = "doc")

    # store data
    bes_resample <- rbind(bes_resample,sample)

    # same for dic 
    sample <- group_by_at(bes_dic_nNA,c("year","month","qt_class"))%>%
      slice_sample(n=1)%>%filter(qt_class %in% c(1,3))%>%ungroup()%>%
      group_by_at(c("year","month"))%>%
      summarize(C_recalc = mean(C_recalc),avg = avg[1])%>%
      mutate(c_diff = (-1*avg-(-1*C_recalc))/(-1*avg),log_int = int,
             scheme = paste(int,"hl1",sep="-"),
             c_type = "dic")

    bes_resample <- rbind(bes_resample,sample)

    ######### calculations for 2 sping-neap sampling events

    # calcs for doc
    sample <- group_by_at(bes_doc_nNA,c("year","month","qt_class"))%>%
      slice_sample(n=2)%>%filter(qt_class %in% c(1,3))%>%ungroup()%>%
      group_by_at(c("year","month"))%>%
      summarize(C_recalc = mean(C_recalc),avg = avg[1])%>%
      mutate(c_diff = (-1*avg-(-1*C_recalc))/(-1*avg),log_int = int,
             scheme = paste(int,"hl2",sep="-"),
             c_type = "doc")

    bes_resample <- rbind(bes_resample,sample)

    # calcs for dic
    sample <- group_by_at(bes_dic_nNA,c("year","month","qt_class"))%>%
      slice_sample(n=2)%>%filter(qt_class %in% c(1,3))%>%ungroup()%>%
      group_by_at(c("year","month"))%>%
      summarize(C_recalc = mean(C_recalc),avg = avg[1])%>%
      mutate(c_diff = (-1*avg-(-1*C_recalc))/(-1*avg),log_int = int,
             scheme = paste(int,"hl2",sep="-"),
             c_type = "dic")

    bes_resample <- rbind(bes_resample,sample)
  
  }
}
