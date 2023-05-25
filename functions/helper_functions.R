#Required Functions

#### action potential extraction ####
AP_extract <- function(loc, df) {
  # loc <- AP_location
  # df <- df
  
  AP <- NULL
  
  for(i in seq_along(loc$sample_loc)){
  
    if(loc$sample_loc[[i]] >= 16 & loc$sample_loc[[i]] <= length(df$MSNA_raw) - 16) {
            
      temp <- cbind(ID = loc$ID[[i]], condition = loc$condition[[i]], ap_time = loc$time[[i]], ap_no = i, sample = c(0:33), ap_v = df$MSNA_raw[(loc$sample_loc[[i]]-16):(loc$sample_loc[[i]]+16)])
      
      AP <- rbind(AP, temp)
      }
   
    incProgress(amount = 1/length(loc$sample_loc)) 
      
  }
  
  AP <- as.data.frame(AP)
  AP$ap_v <- as.numeric(as.character(AP$ap_v)) 
  AP$sample <- as.numeric(as.character(AP$sample))
  
  #browser()
  
  AP <- AP %>% 
    mutate(ID = factor(ID),
           condition = factor(condition, ordered = TRUE),
           ap_time = as.numeric(ap_time),
           ap_no = as.integer(ap_no)
    ) %>% 
    group_by(ap_no, ID, condition, ap_time) %>% 
    nest() %>% 
    group_by(condition) %>%
    mutate(inter_ap_int = (ap_time*1000)-lag(ap_time*1000))
  
  # Determine AP amplitude, minimum negative amplitude, and determine APs set for exclusion.
  # warning: includes those with inter_spike_interval < 2ms)
  AP <- AP %>% 
    mutate(ap_amplitude = as.numeric(map(data, amp)), 
           ap_minimum = as.numeric(map(data, minAP)), 
           ap_include = inter_ap_int > 2 & lead(inter_ap_int > 2)) 
  AP
  
}

#### Amplitude function ####
#determine action potential amplitude
amp <- function(y) {
  ap_amplitude <- max(y$ap_v[11:24])-min(y$ap_v[11:24])
  ap_amplitude
}

#### Minimum function ####
#determine action potential amplitude
minAP <- function(y) {
  ap_minimum <- min(y$ap_v[11:24])
  ap_minimum
}

#### Find Cluster Bins ####
find_cluster_bins <- function(x) {
  
  #browser()
  
  temp <- NULL
  
  for (i in seq_along(unique(x$condition))) {
    h <- x %>% plotly::filter(condition == unique(x$condition)[[i]]) %>%
      .$ap_amplitude %>% hist(breaks = "Scott", plot = FALSE)
    temp <- rbind(temp,cbind(condition = unique(x$condition)[[i]],
                             breaks = h$breaks, mids = h$mids, bin_width = h$breaks - lag(h$breaks))) %>%
      as.data.frame()
  }
  binwidth <- min(as.numeric(temp$bin_width), na.rm = TRUE) # %>% round(2) removed on April 6 2022
  max_bin_mid <- max(as.numeric(temp$mids))
  min_amp <- min(x$ap_amplitude)
  
  bin_centre <- seq(from = max_bin_mid, to = min_amp-(binwidth/2), by = binwidth*-1) %>% rev()
  breaks <- seq(from = min(bin_centre) - binwidth/2, 
                to = max(bin_centre) + binwidth/2, by = binwidth) # added +binwidth to deal with breaks not spanning x range on April 6 2022.
  breaks  
}

#### Find Cluster Bins Version 2 ####
find_cluster_bins2 <- function(x, base_cond) {
  
  #browser()
  
  h <- x %>% plotly::filter(condition == base_cond) %>%
      .$ap_amplitude %>% hist(breaks = "Scott", plot = FALSE)
  
  binwidth <- h$breaks[[2]] - h$breaks[[1]] # %>% round(2) removed on April 6 2022
  max_bin_mid <- max(x$ap_amplitude)
  min_amp <- min(x$ap_amplitude)
  
  bin_centre <- seq(from = max_bin_mid, to = min_amp-(binwidth/2), by = binwidth*-1) %>% rev()
  breaks <- seq(from = min(bin_centre) - binwidth/2, 
                to = max(bin_centre) + binwidth/2, by = binwidth) # added +binwidth to deal with breaks not spanning x range on April 6 2022.
  breaks  
}


#### Generates an index ####
index <- function(x) {
  index <- as.numeric(rownames(x))
  x <- cbind(index = index, x)
}

#### Link APs to beats ####
link_AP <- function(times, ap_beat, APs, signal) {
  
  #browser()
  
  # get beat data from experimental timings
  temp <- NULL
  
  for(i in seq_along(times$condition)) {
    temp1 <- ap_beat %>% filter(beat_time >= times$t_min[[i]] & beat_time <= times$t_max[[i]])
    temp1 <- cbind(ID = times$ID[[i]], condition = times$condition[[i]], temp1)
    temp <- rbind(temp, temp1)
  }
  
  ap_beat <- temp %>% arrange(beat_time) %>% mutate(beat_no = c(1:length(beat_no)),
                                                    ID = factor(ID),
                                                    condition = factor(condition))
  
  # get average burst latency by condition
  ap_beat <- ap_beat %>% group_by(ID, condition) %>% rename(burst_latency = latency) %>%
    mutate(burst_latency_average = mean(burst_latency, na.rm = TRUE),
           tmin = beat_time + (burst_latency_average - (RRI/2)), tmax = beat_time + (burst_latency_average + (RRI/2))) %>% ungroup()
    
    
  # Assign APs to beats/bursts #
  # This loop takes a while # - add progress indicator.
  
  #browser()
  
  temp <- NULL
  temp_noise <- NULL
  
  APs <- APs %>% select(!c(ID, condition))
  
  for(i in seq_along(ap_beat$beat_no)) {
    incProgress(amount = 1/length(ap_beat$beat_no))
    # deal with the case of no APs associated with a beat
    if(APs[APs$ap_time >= ap_beat$tmin[[i]] & APs$ap_time <= ap_beat$tmax[[i]],] %>% nrow() == 0) { 
      
      # extract raw MSNA to determine noise
      temp_df <- signal[signal$time >= ap_beat$tmin[[i]] & signal$time <= ap_beat$tmax[[i]],] %>% 
        select(ID, condition, MSNA_raw) %>% group_by(ID, condition) %>%
        summarize(sd = sd(MSNA_raw, na.rm = TRUE))
      
      temp_noise <- temp_noise %>% 
        rbind(., cbind('beat' = i, 'sd_noise' = temp_df))
      
      temp2 <- APs[APs$ap_time >= ap_beat$tmin[[i]] & APs$ap_time <= ap_beat$tmax[[i]],]
      temp2[1,] <- NA
      temp <- cbind(temp2, ap_beat[i,]) %>% rbind(temp,.)
      
    } else { # deal with case where APs associated with a burst
      
      temp <-  cbind(APs[APs$ap_time >= ap_beat$tmin[[i]] & APs$ap_time <= ap_beat$tmax[[i]],], ap_beat[i,]) %>% rbind(temp, .)
    }
  }
  
  #### remove duplicated APs ####
  df_combined <- temp %>% .[!duplicated(.$ap_time, incomparables = NA),]
  
  #### Calculate SNR ####
  noise <- temp_noise %>% 
    group_by(ID, condition) %>% 
    summarise(beat_num = length(beat), ap_noise = mean(sd, na.rm = TRUE))
  
  df_combined <- df_combined %>% 
    full_join(., noise, by = c("ID", "condition")) %>%
    mutate(SNR = abs(ap_minimum/ap_noise),
           ap_burst_latency = ap_time - burst_latency + (RRI/2) - beat_time,
           ap_latency = ap_time - beat_time,
           ) %>%
    arrange(beat_no, ap_no) %>%
    group_by(beat_no) %>% 
    mutate(inter_ap_int = (ap_time - lag(ap_time))*1000,
           ID = factor(ID),
           condition = factor(condition),
           cluster = factor(cluster)
           ) %>%
    select(!beat_num) %>% relocate(ID, condition) %>%
    arrange(beat_no, ap_no) %>%
    rename(burst_amplitude = amp, burst_area = area) %>%
    relocate(ap_noise, SNR, ap_burst_latency, ap_latency, .before = beat_no)
  df_combined
}




#### LOAD DF and calculate normalize cluster - adapt to fit in app later (i.e. remove read_csv)
# df <- read_csv("df.csv") %>% 
#   group_by(ID, condition, cluster, beat_no, ap_no) %>% 
#   nest(data = c(sample, ap_v)) %>% 
#   relocate(data, .before = ap_include)

#### AP_Summary
Absolute_Summary <- function(df) {

  #browser()
  
  temp <- df %>%
    filter(ap_include == TRUE | is.na(ap_include) & ap_keep == TRUE | is.na(ap_keep)) %>% ungroup() %>%
    group_by(ID, condition) %>%
    nest() %>% filter(!is.na(ID)) %>%
    mutate(
      dt = map(data, ~ {
        max(.x$beat_time) - min(.x$beat_time)
      }),
      no_ap = map(data, ~ {
        length(.x$ap_no[!is.na(.x$ap_no)])
      }),
      no_clusters = map(data, ~ {
        length(unique(.x$cluster[!is.na(.x$cluster)]))
      }),
      no_nclusters = map(data, ~ {
        length(unique(.x$ncluster[!is.na(.x$ncluster)]))
      }),
      no_beats = map(data, ~ {
        length(unique(.x$beat_no[!is.na(.x$beat_no)]))
      }),
      no_bursts = map(data, ~{
        length(unique(.x$burst_time[!is.na(.x$burst_time)]))
      }),
      no_ap_perBurst = map(data, ~ {
        sum(!is.na(.x$ap_no[!is.na(.x$burst_time)])) / length(unique(.x$burst_time[!is.na(.x$burst_time)]))
      }),
      no_synchronous_ap = map(data, ~{
        sum(!is.na(.x$ap_no[!is.na(.x$burst_time)]))
      }),
      percent_synchronous_ap = map(data, ~ {
        sum(!is.na(.x$ap_no[!is.na(.x$burst_time)])) / length(.x$ap_no[!is.na(.x$ap_no)]) * 100
      }),
      percent_asynchronous_ap = map(data, ~ {
        sum(!is.na(.x$ap_no[is.na(.x$burst_time)])) / length(.x$ap_no[!is.na(.x$ap_no)]) * 100
      }),
      no_asynchronous_ap = map(data, ~{
        sum(!is.na(.x$ap_no[is.na(.x$burst_time)]))
      }),
      ap_frequency = map(data, ~ {
        (length(.x$ap_no[!is.na(.x$ap_no)]) / (max(.x$beat_time) - min(.x$beat_time))) * 60
      }),
      ap_incidence = map(data, ~ {
        length(.x$ap_no[!is.na(.x$ap_no)]) / length(unique(.x$beat_no[!is.na(.x$beat_no)])) * 100
      }),
      mean_ap_amplitude = map(data, ~ {
        mean(.x$ap_amplitude, na.rm = TRUE)
      }),
      mean_ap_namplitude = map(data, ~ {
        mean(.x$ap_namplitude, na.rm = TRUE)
      }),
      mean_ap_burstLatency = map(data, ~ {
        mean(.x$ap_burst_latency, na.rm = TRUE)
      }),
      mean_ap_latency = map(data, ~ {
        mean(.x$ap_latency, na.rm = TRUE)
      }),
      mean_inter_ap_interval = map(data, ~ {
        mean(.x$inter_ap_int, na.rm = TRUE)
      }),
      mean_ap_noise = map(data, ~ {
        mean(.x$ap_noise, na.rm = TRUE)
      }),
      mean_ap_SNR = map(data, ~ {
        mean(.x$SNR, na.rm = TRUE)
      }),
      mean_cluster_perBeat = map(data, ~ {
        .x %>% group_by(beat_no) %>% filter(!is.na(cluster)) %>% summarise(n_cluster = length(unique(cluster))) %>% 
          .$n_cluster %>% mean(., na.rm = TRUE)
      }),
      mean_maxCluster_perBurst = map(data, ~ {
        .x %>% filter(!is.na(burst_time) & !is.na(cluster)) %>% group_by(beat_no) %>% summarise(max_cluster = max(as.integer(cluster), na.rm = TRUE)) %>%
          .$max_cluster %>% mean(., na.rm = TRUE)
      }),
      mean_ncluster_perBeat = map(data, ~ {
        .x %>% group_by(beat_no) %>% filter(!is.na(ncluster)) %>% summarise(n_cluster = length(unique(ncluster))) %>% 
          .$n_cluster %>% mean(., na.rm = TRUE)
      }),
      mean_maxnCluster_perBurst = map(data, ~ {
        .x %>% filter(!is.na(burst_time) & !is.na(ncluster)) %>% group_by(beat_no) %>% summarise(max_cluster = max(as.integer(ncluster), na.rm = TRUE)) %>%
          .$max_cluster %>% mean(., na.rm = TRUE)
      }),
      mean_burstAmplitude = map(data, ~{
        .x %>% filter(!is.na(burst_time)) %>% group_by(beat_no) %>% summarise(BurstAmplitude = unique(burst_amplitude)) %>%
          .$BurstAmplitude %>% mean(., na.rm = TRUE)
      }),
      mean_burstNAmplitude = map(data, ~{
        .x %>% filter(!is.na(burst_time)) %>% group_by(beat_no) %>% summarise(BurstNAmplitude = unique(burst_namplitude)) %>%
          .$BurstNAmplitude %>% mean(., na.rm = TRUE)
      }),
      mean_burstArea = map(data, ~{
        .x %>% filter(!is.na(burst_time)) %>% group_by(beat_no) %>% summarise(BurstArea = unique(burst_area)) %>%
          .$BurstArea %>% mean(., na.rm = TRUE)
      }),
      mean_burstNArea = map(data, ~{
        .x %>% filter(!is.na(burst_time)) %>% group_by(beat_no) %>% summarise(BurstNArea = unique(burst_narea)) %>%
          .$BurstNArea %>% mean(., na.rm = TRUE)
      }),
      mean_burstFrequency = map(data, ~ {
        (length(unique(.x$burst_time[!is.na(.x$burst_time)])) /  (max(.x$beat_time) - min(.x$beat_time))) * 60
      }),
      mean_burstIncidence = map(data, ~ {
        length(unique(.x$burst_time[!is.na(.x$burst_time)])) /  length(unique(.x$beat_no[!is.na(.x$beat_no)])) * 100
      }),
      CV_Summary = map(data, ~ {
        .x %>% select(!c(cluster:ap_latency, beat_time, burst_time:tmax)) %>% group_by(beat_no) %>% summarize_all(mean, na.rm = TRUE) %>% ungroup() %>% 
          select(!beat_no) %>% summarize_all(mean, na.rm = TRUE)
      })
    ) %>% unnest(CV_Summary) %>% 
    nest(Absolute_Summary = c(dt:last_col()))
}

#df <- temp %>% unnest(data)

#### Normalized Cluster Summary
nCluster_Summary <- function(df) {
  
  #browser()
  temp <- df %>%
    filter(ap_include == TRUE | is.na(ap_include) & ap_keep == TRUE | is.na(ap_keep)) %>% ungroup() %>%
    unnest(Absolute_Summary, names_sep = "_") %>%
    group_by(ID, condition, ncluster) %>%
    nest() %>% filter(!is.na(ID)) %>%
    mutate(
      no_ap = map(data, ~ {
        length(.x$ap_no[!is.na(.x$ap_no)])
      }),
      no_beats = map(data, ~ {
        length(unique(.x$beat_no[!is.na(.x$beat_no)]))
      }),
      no_bursts = map(data, ~{
        length(unique(.x$burst_time[!is.na(.x$burst_time)]))
      }),
      no_ap_perBurst = map(data, ~ {
        sum(!is.na(.x$ap_no[!is.na(.x$burst_time)])) / length(unique(.x$burst_time[!is.na(.x$burst_time)]))
      }),
      no_synchronous_ap = map(data, ~{
        sum(!is.na(.x$ap_no[!is.na(.x$burst_time)]))
      }),
      percent_synchronous_ap = map(data, ~ {
        sum(!is.na(.x$ap_no[!is.na(.x$burst_time)])) / length(.x$ap_no[!is.na(.x$ap_no)]) * 100
      }),
      percent_asynchronous_ap = map(data, ~ {
        sum(!is.na(.x$ap_no[is.na(.x$burst_time)])) / length(.x$ap_no[!is.na(.x$ap_no)]) * 100
      }),
      no_asynchronous_ap = map(data, ~{
        sum(!is.na(.x$ap_no[is.na(.x$burst_time)]))
      }),
      Proportion_firing_1 = map(data, ~ {
        .x %>% group_by(beat_no) %>% summarize(no_ap = sum(!(is.na(unique(ap_no))))) %>% filter(no_ap == 1) %>% .$no_ap %>% length() /
          length(.x$ap_no[!is.na(.x$ap_no)]) * 100
      }),
      Proportion_firing_2 = map(data, ~ {
        .x %>% group_by(beat_no) %>% summarize(no_ap = sum(!(is.na(unique(ap_no))))) %>% filter(no_ap == 2) %>% .$no_ap %>% length() /
          length(.x$ap_no[!is.na(.x$ap_no)]) * 100
      }),
      Proportion_firing_3 = map(data, ~ {
        .x %>% group_by(beat_no) %>% summarize(no_ap = sum(!(is.na(unique(ap_no))))) %>% filter(no_ap == 3) %>% .$no_ap %>% length() /
          length(.x$ap_no[!is.na(.x$ap_no)]) * 100
      }),
      Proportion_firing_4 = map(data, ~ {
        .x %>% group_by(beat_no) %>% summarize(no_ap = sum(!(is.na(unique(ap_no))))) %>% filter(no_ap >= 4) %>% .$no_ap %>% length() /
          length(.x$ap_no[!is.na(.x$ap_no)]) * 100
      }),
      ap_frequency = map(data, ~ {
        length(.x$ap_no[!is.na(.x$ap_no)]) / mean(unlist(.x$Absolute_Summary_dt)) * 60
      }),
      ap_incidence = map(data, ~ {
        length(.x$ap_no[!is.na(.x$ap_no)]) / mean(unlist(.x$Absolute_Summary_no_beats)) * 100
      }),
      mean_ap_amplitude = map(data, ~ {
        mean(.x$ap_amplitude, na.rm = TRUE)
      }),
      mean_ap_burstLatency = map(data, ~ {
        mean(.x$ap_burst_latency, na.rm = TRUE)
      }),
      mean_ap_latency = map(data, ~ {
        mean(.x$ap_latency, na.rm = TRUE)
      }),
      mean_inter_ap_interval = map(data, ~ {
        mean(.x$inter_ap_int, na.rm = TRUE)
      }),
      mean_ap_noise = map(data, ~ {
        mean(.x$ap_noise, na.rm = TRUE)
      }),
      mean_ap_SNR = map(data, ~ {
        mean(.x$SNR, na.rm = TRUE)
      }),
      mean_burstAmplitude = map(data, ~{
        .x %>% filter(!is.na(burst_time)) %>% group_by(beat_no) %>% summarise(BurstAmplitude = unique(burst_amplitude)) %>%
          .$BurstAmplitude %>% mean(., na.rm = TRUE)
      }),
      mean_burstNAmplitude = map(data, ~{
        .x %>% filter(!is.na(burst_time)) %>% group_by(beat_no) %>% summarise(BurstNAmplitude = unique(burst_namplitude)) %>%
          .$BurstNAmplitude %>% mean(., na.rm = TRUE)
      }),
      mean_burstArea = map(data, ~{
        .x %>% filter(!is.na(burst_time)) %>% group_by(beat_no) %>% summarise(BurstArea = unique(burst_area)) %>%
          .$BurstArea %>% mean(., na.rm = TRUE)
      }),
      mean_burstNArea = map(data, ~{
        .x %>% filter(!is.na(burst_time)) %>% group_by(beat_no) %>% summarise(BurstNArea = unique(burst_narea)) %>%
          .$BurstNArea %>% mean(., na.rm = TRUE)
      }),
      mean_burstFrequency = map(data, ~ {
        (length(unique(.x$burst_time[!is.na(.x$burst_time)])) /  mean(unlist(.x$Absolute_Summary_dt))) * 60
      }),
      mean_burstIncidence = map(data, ~ {
        length(unique(.x$burst_time[!is.na(.x$burst_time)])) /  mean(unlist(.x$Absolute_Summary_no_beats)) * 100
      }),
      nCluster_CV_Summary = map(data, ~ {
        .x %>% select(!c(cluster:ap_latency, beat_time, burst_time:last_col())) %>% group_by(beat_no) %>% summarize_all(mean, na.rm = TRUE) %>% ungroup() %>% 
          select(!beat_no) %>% summarize_all(mean, na.rm = TRUE)
      })
    ) %>% unnest(nCluster_CV_Summary) %>% 
    nest(nCluster_Summary = c(no_ap:last_col())) %>%
    unnest(data) %>% arrange(beat_no, ap_no) %>%
    nest(Absolute_Summary = c(Absolute_Summary_dt:last_col(1)))
}





#### correlogram histogram ####
## x and y are time series.  For autocorrelation input the same time series for both x and y.  For cross-correlation input different time series.  Then plot the vectors of inter stimulus intervals to geom_histogram with a 50ms binwidth and a time of -5 to 5s for standard MSNA vs ECG or MSNA vs RESP or ECG vs ECG or resp vs resp plots.
correlation_histogram <- function(x,y) {
  temp <- NULL
  for (i in 1:length(x$TIME)) {
    
    temp <- (y$TIME - x$TIME[[i]]) %>% c(temp, .)
  }
  return(temp)
}


#### lm analysis of baroreflex curves for ap brs sensitivity ####
ap_brs_model <- function(df) {
  lm(cluster_unique ~ value, data = df)
}

#### lm analysis of baroreflex curves for ap brs threshold ####
threshold_model <- function(df) {
  lm(probability ~ DBP, weights = No_beats, data = df)
}

#### function to set up ap brs threshold ####
ap_brs_threshold <- function(df) {
  Num_beats <- length(unique(df$beat_no))
  DBP <- mean(df$value)
  df <- df %>% 
    select(cluster_norm, ap_per_ncluster) %>% 
      group_by(cluster_norm) %>% 
        summarise(probability = (sum(ap_per_ncluster)/Num_beats)*100) %>% 
          mutate(DBP = DBP, No_beats = Num_beats)
  df  
}

#### function to compile a dataframe of brs data for each ap cluster ####


cluster_compile <- function(x, cluster_list) {
  
  DBP_bins <- unique(x$bins)
  
  temp <- NULL
  for (i in seq_along(cluster_list)) {
    
    temp2 <- NULL
    
    for (j in seq_along(DBP_bins)) {
      #j <- 1
      
      temp3 <- x %>% filter(bins == DBP_bins[[j]])
      
      #deal with case where desired cluster is present
      if (cluster_list[[i]] %in% temp3$cluster_norm) {
        temp2 <- temp3 %>% filter(cluster_norm ==cluster_list[[i]]) %>% rbind(temp2, .)
      } else { #deal with case where desired cluster is not present and set probability to 0
        temp2 <- rbind(temp2, c("brs_stimulus" = "DBP", "bins" = as.character(temp3$bins[[1]]), "cluster_norm" =cluster_list[[i]], "probability" =0, "DBP" = temp3$DBP[[1]], "No_beats" = temp3$No_beats[[1]])) %>% as.data.frame()
      }
      
    }
   
   if (sum(as.numeric(temp2$probability)) > 0) {
     temp <- temp2 %>% rbind(temp, .)
   } 
     
  }
  
  temp$brs_stimulus <- factor(temp$brs_stimulus)
  temp$cluster_norm <- as.numeric(temp$cluster_norm)
  temp$probability <- as.numeric(temp$probability)
  temp$DBP <- as.numeric(temp$DBP)
  temp$No_beats <- as.numeric(temp$No_beats)
  temp <- rename(temp, cluster = cluster_norm)
  temp
}
