#Required Functions


#### action potential extraction ####
AP_extract <- function(loc, df) {
  # loc <- AP_location
  # df <- df
  
  AP <- NULL
  
  for(i in seq_along(loc$sample_loc)){
  
    if(loc$sample_loc[[i]] >= 16 & loc$sample_loc[[i]] <= length(df$MSNA_raw) - 16) {
            
      temp <- cbind(ID = loc$ID[[i]], condition = loc$condition[[i]], time = loc$time[[i]], AP_no = i, sample = c(0:33), AP_v = df$MSNA_raw[(loc$sample_loc[[i]]-16):(loc$sample_loc[[i]]+16)])
      
      AP <- rbind(AP, temp)
      }
   
    incProgress(amount = 1/length(loc$sample_loc)) 
      
  }
  
  AP <- as.data.frame(AP)
  AP$AP_v <- as.numeric(as.character(AP$AP_v)) 
  AP$sample <- as.numeric(as.character(AP$sample))
  
  AP <- AP %>% 
    mutate(ID = factor(ID),
           condition = factor(condition, ordered = TRUE),
           time = as.numeric(time),
           AP_no = as.integer(AP_no)
    ) %>% 
    group_by(AP_no, ID, condition, time) %>% 
    nest() %>% 
    group_by(condition) %>%
    mutate(inter_spike_int = (time*1000)-lag(time*1000))
  
  # Determine AP amplitude, minimum negative amplitude, and determine APs set for exclusion.
  # warning: includes those with inter_spike_interval < 2ms)
  AP <- AP %>% 
    mutate(amplitude = as.numeric(map(data, amp)), 
           minimum = as.numeric(map(data, minAP)), 
           include = inter_spike_int > 3)
  AP
  
}

#### Amplitude function ####
#determine action potential amplitude
amp <- function(y) {
  amplitude <- max(y$AP_v[11:24])-min(y$AP_v[11:24])
  amplitude
}

#### Minimum function ####
#determine action potential amplitude
minAP <- function(y) {
  minimum <- min(y$AP_v[11:24])
  minimum
}





#### Find Cluster Bins ####
find_cluster_bins <- function(x) {
  
  temp <- NULL
  
  for (i in seq_along(unique(x$condition))) {
    h <- x %>% filter(condition == unique(x$condition)[[i]]) %>%
      .$amplitude %>% hist(breaks = "Scott", plot = FALSE)
    temp <- rbind(temp,cbind(condition = unique(x$condition)[[i]],
                             breaks = h$breaks, mids = h$mids, bin_width = h$breaks - lag(h$breaks))) %>%
      as.data.frame()
  }
  binwidth <- min(temp$bin_width, na.rm = TRUE) %>% round(2)
  max_bin_mid <- max(temp$mids)
  min_amp <- min(x$amplitude)
  
  bin_centre <- seq(from = max_bin_mid, to = min_amp-(binwidth/2), by = binwidth*-1) %>% rev()
  breaks <- seq(from = min(bin_centre) - binwidth/2, 
                to = max(bin_centre) + binwidth/2, by = binwidth)
  breaks  
}

#### Generates an index ####
index <- function(x) {
  index <- as.numeric(rownames(x))
  x <- cbind(index = index, x)
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
