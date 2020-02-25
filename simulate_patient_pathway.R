# probability of occurence (truncated gaussian)
proba_occurence <- function(time, mean, sd){
  p1 <- ptruncnorm(q = time-1, a = 0, b = +Inf, mean = mean, sd = sd)
  p2 <- ptruncnorm(q = time, a = 0, b = +Inf, mean = mean, sd = sd)
  return( (p2 - p1) / (1 - p1) )
}

# probability of shift (exponential for now)
proba_shift <- function(time1, time2, mean = 30){
  pexp(q = time2 - time1, rate = 1 / mean)
  # p1 <- pexp(q = time1, rate = 1 / mean)
  # p2 <- pexp(q = time2, rate = 1 / mean)
  # return( (p2 - p1) / (1 - p1) )
}

# simulation parameters
NB_DAYS_SIMULATION <- 365 * 4
SYNDROME <- syndrome_metadata$id[2]

# collect number of symptoms from metadata
NB_SYMPTOMS <- nrow(symptom_metadata)

#
one_simulation <- function(NB_DAYS_SIMULATION = 365, SYNDROME = NULL, SIMULATION_SEED = 643, proba_occurence, proba_shift){
  
  # identify subset of cross metadata relevant to the syndrome
  syndrome_specific_cross_metadata <- cross_metadata[cross_metadata$syndrome_id == SYNDROME,]
  
  # simulation "container"
  sim <- data.frame('time' = 1:NB_DAYS_SIMULATION, 'syndrome' = SYNDROME)
  mat <- array(data = FALSE, dim = c(NB_DAYS_SIMULATION, NB_SYMPTOMS * 3))
  colnames(mat) <- c(
    paste0(symptom_metadata$id, '_is_present'),
    paste0(symptom_metadata$id, '_is_observed'),
    paste0(symptom_metadata$id, '_has_been_observed')
  )
  sim <- cbind.data.frame(sim, mat)
  
  # identify columns of "sim" dataframe for easy manipulation
  colindexes_is_present <- grep(x = colnames(sim), pattern = '_is_present')
  colindexes_has_been_observed <- grep(x = colnames(sim), pattern = '_has_been_observed')
  colindexes_is_observed <- grep(x = colnames(sim), pattern = '_is_observed')
  
  # identify rows of symptom's metadata for each type
  latent_symptoms <- symptom_metadata$type == 'latent'
  chronic_symptoms <- symptom_metadata$type == 'chronic'
  beginning_to_end_symptoms <- symptom_metadata$type == 'beginning_to_end'
  
  # simulation seed
  set.seed(SIMULATION_SEED)
  
  # can be randomized, just uncomment the last line ...
  symptoms_that_shall_be_present <- (runif(n = NB_SYMPTOMS) < syndrome_specific_cross_metadata$proba_symptom_cond_syndrome)
  
  # iterate over time indexes (days)
  for(time_index in 1:NB_DAYS_SIMULATION){
    
    ###
    ### check if new symptom is present
    ###
    
    # which symptoms already are present ?
    if(time_index >= 2){
      check_0 <- sim[time_index-1, colindexes_is_present]
    }else{
      check_0 <- rep(FALSE, NB_SYMPTOMS)
    }
    # test occurence
    proba_0 <- proba_occurence(time = time_index, mean = syndrome_specific_cross_metadata$mean_time_before_symptom_declaration, sd = syndrome_specific_cross_metadata$sd_time_before_symptom_declaration)
    test_0 <- runif(n = NB_SYMPTOMS) < proba_0
    # final
    new_sym <- (symptoms_that_shall_be_present & (check_0 == FALSE) & test_0)
    sim[time_index, colindexes_is_present] <- check_0
    if(any(new_sym)){
      sim[time_index, colindexes_is_present[new_sym]] <- TRUE
    }
    
    ###
    ### now, qualify which symptom are currently / have been observed
    ###
    
    # has been observed ...
    if(time_index >= 2){
      sim[time_index, colindexes_has_been_observed] <- sim[time_index-1, colindexes_has_been_observed]
    }else{
      sim[time_index, colindexes_has_been_observed] <- rep(FALSE, NB_SYMPTOMS)
    }
    non_latent_new_sym <- new_sym & !latent_symptoms
    if(any(non_latent_new_sym)){
      sim[time_index, colindexes_has_been_observed[non_latent_new_sym]] <- TRUE
    }
    
    # currently observed (only beggining-to-end and chronic ones)
    sim[time_index, colindexes_is_observed[beginning_to_end_symptoms]] <- sim[time_index, colindexes_is_present[beginning_to_end_symptoms]]
    if(time_index >= 2){
      last_shift_time <- apply(X = sim[1:(time_index-1), colindexes_is_observed[chronic_symptoms]], MARGIN = 2, FUN = function(x){
        if(any(x)){
          v <- which(x)
          return(v[length(v)])
        }else{
          return(1)
        } 
      })
      proba_of_shift <- proba_shift(time1 = last_shift_time, time2 = time_index, mean = 30)
      shifting <- (runif(n = sum(chronic_symptoms)) < proba_of_shift)
      new_state <- sim[time_index-1, colindexes_is_observed[chronic_symptoms]]
      new_state[shifting] <- (!new_state[shifting])
      sim[time_index, colindexes_is_observed[chronic_symptoms]] <- sim[time_index, colindexes_is_present[chronic_symptoms]] & new_state 
    }else{
      sim[time_index, colindexes_is_observed[chronic_symptoms]] <- sim[time_index, colindexes_is_present[chronic_symptoms]]  
    }
    
    # if we just observed the symptom, we update has_been_observed
    sim[time_index, colindexes_has_been_observed] <- (sim[time_index, colindexes_has_been_observed] | sim[time_index, colindexes_is_observed])
    
  }
  
  return(sim)
  
}

#
NB_SYNDROMES <- nrow(syndrome_metadata)
B <- 100
set.seed(0)
seeds <- ceiling(runif(n = B)*1e06)
list_simulations <- list()
predictor_training_data <- NULL
count <- 0
pb <- txtProgressBar(min = 0, max = NB_SYNDROMES * B)
for(s in 1:NB_SYNDROMES){
  for(b in 1:B){
    count <- count + 1
    setTxtProgressBar(pb = pb, value = count)
    temp <- one_simulation(NB_DAYS_SIMULATION = NB_DAYS_SIMULATION, SYNDROME = syndrome_ids[s], SIMULATION_SEED = seeds[b],
                           proba_occurence = proba_occurence, proba_shift = proba_shift)
    temp$time_since_first_symptom <- 0
    time_first_symptom <- min(apply(X = temp[,colindexes_is_observed], MARGIN = 2, FUN = function(x){ which(x)[1] }), na.rm = TRUE)
    if(is.infinite(time_first_symptom) == FALSE){
      #
      temp$time_since_first_symptom[time_first_symptom:NB_DAYS_SIMULATION] <- 1:(NB_DAYS_SIMULATION - time_first_symptom + 1)
      #
      temp2 <- temp[time_first_symptom:NB_DAYS_SIMULATION,colindexes_has_been_observed]
      temp2$time_since_first_symptom <- 1:nrow(temp2)
      temp2$syndrome <- syndrome_ids[s]
      predictor_training_data <- rbind(predictor_training_data, temp2)
    }
    # 
    list_simulations[[length(list_simulations)+1]] <- temp
  }
}

#
length(list_simulations)
dim(predictor_training_data)
str(predictor_training_data)
## colnames(predictor_training_data)[11] <- 'time_since_first_symptom'

#
library(randomForest)
syndrome_metadata
table(predictor_training_data$syndrome)
predictor_training_data$syndrome_is_rare <- factor(x = predictor_training_data$syndrome %in% syndrome_metadata$id[syndrome_metadata$is_rare], levels = c(TRUE, FALSE))
col_index_syndrome <- which(colnames(predictor_training_data) == 'syndrome')
dim(predictor_training_data)
set.seed(28543)
row_indexes <- sample.int(n = nrow(predictor_training_data), size = 100000, replace = TRUE)
rf0 <- randomForest(formula = syndrome_is_rare ~ ., data = predictor_training_data[row_indexes,-col_index_syndrome])
rf0
varImpPlot(rf0)

#
ls_symptom_observed <- sapply(list_simulations, function(x){
  any(x$time_since_first_symptom > 0)
})
ls_symptom_observed <- which(ls_symptom_observed)

#
threshold_discrete <- seq(from = 0, to = 1, length.out = 100)
observed_cost <- stopping_time <- array(data = NA, dim = c(length(ls_symptom_observed), length(threshold_discrete)))

#
costs <- c(
  'errance_par_jour' = 0.08, # 250 / NB_DAYS_SIMULATION * 4 * 2,
  'medecin_non_spe_MR' = 75,
  'medecin_spe_MR' = 200,
  'nb_moyen_medecin_consultes' = 5,
  'temps_moyen_errance' = NB_DAYS_SIMULATION
)

#
for(j in 1:length(ls_symptom_observed)){
  lindex <- ls_symptom_observed[j]
  list_simulations[[lindex]]$proba_rare_disease <- 
    predict(rf0, list_simulations[[lindex]], type = 'prob')[,1]
  list_simulations[[lindex]]$proba_rare_disease[list_simulations[[lindex]]$time_since_first_symptom == 0] <- NA
  time_since_first_symptom_observed <- which(list_simulations[[lindex]]$time_since_first_symptom == 1)
  for(k in 1:length(threshold_discrete)){
    time_threshold_exceeded <- which(list_simulations[[lindex]]$proba_rare_disease >= threshold_discrete[k])[1]
    syndrome_is_rare <- (list_simulations[[lindex]]$syndrome[1] %in% syndrome_metadata$id[syndrome_metadata$is_rare])
    recommendation_to_RDC <- (is.na(time_threshold_exceeded) == FALSE)
    if(syndrome_is_rare){
      if(recommendation_to_RDC){
        observed_cost[j,k] <- costs['errance_par_jour'] * (time_threshold_exceeded - time_since_first_symptom_observed + 1) + costs['medecin_spe_MR']
      }else{
        observed_cost[j,k] <- costs['errance_par_jour'] * costs['temps_moyen_errance'] + costs['medecin_non_spe_MR'] * costs['nb_moyen_medecin_consultes']
      }
    }else{
      if(recommendation_to_RDC){
        observed_cost[j,k] <- costs['errance_par_jour'] * (time_threshold_exceeded - time_since_first_symptom_observed + 1)  + costs['medecin_spe_MR']
      }else{
        observed_cost[j,k] <- costs['errance_par_jour'] * (NB_DAYS_SIMULATION - time_since_first_symptom_observed + 1) + costs['medecin_non_spe_MR']
      }
    }
    if(is.na(time_threshold_exceeded) == FALSE){
      stopping_time[j,k] <- time_threshold_exceeded - list_simulations[[lindex]]$time_since_first_symptom[time_threshold_exceeded]
    }
  }
}

# compute weight vector
weight_vector <- rep(NA, nrow(observed_cost))
syndrome_sam <- sapply(X = 1:length(ls_symptom_observed), FUN = function(j){
  lindex <- ls_symptom_observed[j]
  list_simulations[[lindex]]$syndrome[1]
})
table(syndrome_sam)
for(k in 1:nrow(syndrome_metadata)){
  weight_vector[syndrome_sam == syndrome_metadata$id[k]] <- syndrome_metadata$proba[k]
} 
weight_vector <- weight_vector / sum(weight_vector)
table(weight_vector)

# compute observed cost average and standard deviation
str(observed_cost)
res <- data.frame(
  'threshold' = threshold_discrete,
  'mean_cost' = apply(observed_cost, 2, function(x) mean(x*weight_vector)),
  'sd_cost' = apply(observed_cost, 2, function(x) sd(x*weight_vector))
)
ggplot(data = res[res$threshold >= 0.03,], mapping = aes(x = threshold, y = (mean_cost - min(mean_cost))/(max(mean_cost) - min(mean_cost)) )) +
  geom_line() + geom_point(col = 'blue', cex = 0.5) + 
  xlim(0,1) + ylim(0,1) +
  xlab('Seuil') + ylab('Coût moyen normalisé sur [0,1]') +
  ggtitle(label = "Coût moyen, normalisé, observé sur les données d'apprentissage", subtitle = "selon le seuil sélectionné.")
ggplot(data = res, mapping = aes(x = threshold, y = mean_cost)) +
  geom_line() +
  geom_errorbar(mapping = aes(ymin = mean_cost, ymax = mean_cost + 2 * sd_cost))

#
sapply(list_simulations, function(x){ max(x$proba_rare_disease, na.rm = TRUE)})
j = 102
plot.ts(list_simulations[[j]]$proba_rare_disease)

# preparing data for proper representation
sim <- list_simulations[[j]]
colindexes_is_observed <- grep(x = colnames(sim), pattern = '_is_observed')
symptoms_that_shall_be_present <- colSums(sim[,colindexes_is_observed]) > 0
is_observed_data <- gather(data = sim[,c(1, which(colnames(sim) == 'proba_rare_disease'), colindexes_is_observed[symptoms_that_shall_be_present])], key = key, value = value, -time)
is_observed_data$type <- 'is_observed'
global_data <- is_observed_data
global_data$key_uni <- sub(pattern = '\\_is\\_present|\\_is_\\observed|\\_has\\_been\\_observed', replacement = '', x = global_data$key)
global_data$key_uni <- factor(global_data$key_uni, levels = levels(symptom_metadata$id))
label_mr <- "P('Maladie Rare')"
global_data$symptom_label <- label_mr 
for(i in 1:nrow(symptom_metadata)){
  global_data$symptom_label[global_data$key_uni == symptom_metadata$id[i]] <- as.character(symptom_metadata$label[i])
}
global_data$symptom_label <- factor(global_data$symptom_label)

global_data$type <- factor(global_data$type)
nb_days <- 200
vect_color <- gg_color_hue(n = length(unique(global_data$symptom_label)))
colname_proba <- "proba_rare_disease"
global_data$my_color <- vect_color[global_data$symptom_label]
g0 <- ggplot(data = global_data[global_data$key != colname_proba & global_data$time <= nb_days,], mapping = aes(x = time, y = as.numeric(value), col = my_color )) +
  geom_line(cex = 2, col = global_data$my_color[global_data$key != colname_proba &  global_data$time <= nb_days]) +
  facet_grid(symptom_label ~ .) +
  scale_y_continuous(breaks = c(0, 1),
                     labels = c("Non", "Oui")) +
  xlab(label = '') +
  ylab(label = 'Observation') +
  ### ggtitle(label = 'Presence and observation of symptoms throughout time.') +
  theme(legend.position = 'none', 
        strip.text.y = element_text(angle = 270), 
        axis.text.x.bottom = element_blank(), 
        axis.ticks.x.bottom = element_blank(),
        axis.title.x = element_blank()) 
g0
g1 <- ggplot(data = global_data[global_data$key == colname_proba &  global_data$time <= nb_days,], 
             mapping = aes(x = time, y = as.numeric(value))) +
  geom_hline(yintercept = 0.5, lty = 2) +
  geom_line(cex = 2, col = global_data$my_color[global_data$key == colname_proba &  global_data$time <= nb_days]) +
  ylim(0,1) + 
  ## scale_y_continuous(breaks = seq(from = 0, to = 1.1, length.out = 15), minor_breaks = NULL) +
  facet_grid(symptom_label ~ .) +
  xlab(label = 'Temps, en jours') +
  ylab(label = 'Prédiction') +
  ## ggtitle(label = '') +
  theme(legend.position = 'none', 
        strip.text.y = element_text(angle = 270), 
        axis.title.x.top = element_blank()) 
  ## scale_x_continuous(breaks = seq(from = 0, to = 200, length.out = 9), minor_breaks = NULL) 
g1
library(gridExtra)
g0
g1
## grid.arrange(p1, arrangeGrob(p2,p3,p4, ncol=3), heights=c(2.5/4, 1.5/4), ncol=1)
val <- 0.35
grid.arrange(g0, g1, nrow = 2, ncol = 1, heights = c(1-val, val))

