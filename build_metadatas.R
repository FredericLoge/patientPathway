# # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#     HANDLING LIBRARIES (to not do it afterwards ...)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# load viz library ggplot2
library(ggplot2)

# truncated normal distribution library
library(truncnorm)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#               CONSTRUCTING METADATAS (I)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Metadata on syndromes considered
# ntr = nothing to report
# rd = rare disease
# cd = common disease
syndrome_metadata <- data.frame(
  'id' = paste0('syndrome_', c('ntr', 1:3)),
  'label' = c('RAS', paste0('Syndrome #', 1:3)),
  'is_rare' = c(F, T, F, F),
  'proba' = c(0.40, 0.01, 0.30, 0.29)
)

# Metadata on symptoms considered
# ntr = nothing to report
set.seed(seed = 896543)
N <- 10
set_symptom_types <- c('beginning_to_end', 'chronic', 'latent')
symptom_metadata <- data.frame(
  'id' = paste0('symptom_', 1:N),
  'label' = paste0('Symptôme ', 1:N),
  ## 'label' = paste0('name_symptom_', 1:N),
  'type' = sample(x = set_symptom_types, size = N, replace = TRUE)
)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#               CONSTRUCTING METADATAS (II)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Defining : 
# - independent probability of occurence 
# - conditional occurence times of symptomes conditional on syndromes
# - for chronical syndromes : periodicity

# cross syndromes and symptoms
cross_metadata <- expand.grid(
  'syndrome_id' = syndrome_metadata$id,
  'symptom_id' = symptom_metadata$id
)

rnormalmixture <- function(n, nu, mu, sd){
  grp <- sample(x = 1:length(nu), size = n, replace = TRUE, prob = nu)
  val <- rep(NA, n)
  for(i in 1:length(nu)){
    val[grp == i] <- rnorm(n = sum(grp == i), mean = mu[i], sd = sd[i])
  }
  return(val)
}
# test <- rnormalmixture(n = 1000, nu = c(0.5, 0.2, 0.3), mu = c(-10, 0, 10), sd = c(3,3,3)/100)
# hist(test, breaks = 30)

# set seed for reproducibility of results
set.seed(seed = 19876543)

# probability of observing symptom conditionally to syndrome
cross_metadata$proba_symptom_cond_syndrome <- rbeta(n = nrow(cross_metadata), shape1 = 1/4, shape2 = 1)

# mean and sd of time between syndrome occurence and (first) symptom occurence
cross_metadata$mean_time_before_symptom_declaration <- rnormalmixture(n = nrow(cross_metadata), nu = c(1/3, 1/3, 1/3), mu = c(30, 60, 90), sd = rep(5, 3))
cross_metadata$sd_time_before_symptom_declaration <- 10

# mean time between symptom occurences
cross_metadata$symptom_periodicity <- rchisq(n = nrow(cross_metadata), df = 30, ncp = 0)
temp_chronic_symptom <- cross_metadata$symptom_id %in% symptom_metadata$id[symptom_metadata$type == 'chronic']
cross_metadata$symptom_periodicity[temp_chronic_symptom == FALSE] <- NA

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#           PLOT GRAPH
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1 + 1:n]
}
symptom_color_code <- gg_color_hue(n = nrow(symptom_metadata))

library(DiagrammeR)
#
cross_selection <- (cross_metadata$proba_symptom_cond_syndrome > 0.15)
table(cross_selection)
syndrome_ok_index <- syndrome_metadata$id %in% cross_metadata$syndrome_id[cross_selection]
symptom_ok_index <- symptom_metadata$id %in% cross_metadata$symptom_id[cross_selection]
#
cmd_syndrome_nodes <- paste0(
  "node [fontname = Helvetica, shape = rectangle, style = filled, sides = black, color = black, fontcolor = white, label = '", 
  syndrome_metadata$label[syndrome_ok_index],
  "'] \n ", 
  syndrome_metadata$id[syndrome_ok_index], 
  collapse = " \n\n ")
cat(cmd_syndrome_nodes)
#
cmd_symptom_nodes <- paste0(
  "node [shape = rectangle, style = filled, color = '",
  symptom_color_code,
  "', fontcolor = white, label = '", 
  symptom_metadata$label[symptom_ok_index],
  "'] \n ", 
  symptom_metadata$id[symptom_ok_index],
  collapse = " \n\n ")
cat(cmd_symptom_nodes)
#
cmd_edges <- paste0(
  "edge [color = black, arrowhead = normal, penwidth = ", 4*cross_metadata$proba_symptom_cond_syndrome[cross_selection], "] \n", 
  cross_metadata$syndrome_id[cross_selection],
  " -> ", 
  cross_metadata$symptom_id[cross_selection],
  collapse = " \n\n ")
cat(cmd_edges)
# circo, 
cmd_diagram <- paste0("digraph mygraph{
  graph [overlap = FALSE, layout = twopi] \n\n ",
                      cmd_syndrome_nodes, "\n\n",
                      cmd_symptom_nodes, "\n\n",
                      cmd_edges, "\n\n",
                      "}")
cat(cmd_diagram)
grViz(diagram = cmd_diagram)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#           PATIENT PATHWAY SIMULATION
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

str(cross_metadata)
library(tidyverse)
cross_metadata %>%
  filter(proba_symptom_cond_syndrome > 0.002) %>%
  arrange(syndrome_id, desc(proba_symptom_cond_syndrome)) %>%
  View()

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#           VISUALIZE SYMPTOM OCCURENCE DISTRIBUTION
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# list of syndrome IDs
syndrome_ids <- syndrome_metadata$id

# iterate over IDs
global_mat <- NULL
for(syndrome_id in syndrome_ids){
  
  # collect syndrome metadata  
  temp <- cross_metadata[cross_metadata$syndrome_id == syndrome_id,]
  
  # generate time indices to represent distribution
  max_time <- pmax(100, round(max(temp$mean_time_before_symptom_declaration)*2,-2))
  time_ref <- seq(from = 0, to = max_time, length.out = 1000)
  
  # collect (properly-shaped) : (syndrome, symptom, time, proba global occurence symptom, proba occurence at that time, max proba)
  for(row_index_temp in 1:nrow(temp)){
    mat <- cbind.data.frame(
      'V1' = temp$syndrome_id[row_index_temp],
      'V2' = temp$symptom_id[row_index_temp],
      'V3' = time_ref,
      'V4' = temp$proba_symptom_cond_syndrome[row_index_temp],
      ## dexp(x = time_ref, rate = 1 / temp$mean_time_before_symptom_declaration[row_index_temp]),
      'V5' = dtruncnorm(x = time_ref, a = 0, b = +Inf, mean = temp$mean_time_before_symptom_declaration[row_index_temp], sd = temp$sd_time_before_symptom_declaration[row_index_temp]),
      'V6' = FALSE
    )
    
    # identify max proba time point and flag if symptom is typical (> 15%)
    if(temp$proba_symptom_cond_syndrome[row_index_temp] > 0.15){
      row_index_max <- which.max(mat[,5])
      mat[row_index_max,6] <- TRUE
    }
    
    # add to previous dataset
    global_mat <- rbind.data.frame(global_mat, mat)
    
  }
  
}
str(global_mat)
colnames(global_mat) <- c('syndrome_id', 'symptom_id', 'time', 'proba_symptom_cond_syndrome', 'proba_occurence', 'flag_max_proba_occurence')

#
global_mat$syndrome_label <- syndrome_metadata$label[1]
for(i in 1:nrow(syndrome_metadata)){
  global_mat$syndrome_label[global_mat$syndrome_id == syndrome_metadata$id[i]] <- syndrome_metadata$label[i]
}
table(global_mat$syndrome_label)
global_mat$symptom_label <- symptom_metadata$label[1]
for(i in 1:nrow(symptom_metadata)){
  global_mat$symptom_label[global_mat$symptom_id == symptom_metadata$id[i]] <- symptom_metadata$label[i]
}
table(global_mat$symptom_label)

# represent proba_occurence
ggplot(data = global_mat) +
  geom_line(mapping = aes(x = time, y = proba_occurence, col = symptom_label)) +
  facet_wrap(.~syndrome_label) +
  labs(col = 'Symptom') +
  ylab(label = 'Proba of occurence') +
  xlab(label = 'Time in days')

# represent proba_occurence * proba_symptom_cond_syndrome
x_offset_heuristic <- 0 # 200 / 8
y_offset_heuristic <- with(global_mat, max(proba_occurence * proba_symptom_cond_syndrome)) / (6*2)
ggplot(data = global_mat, mapping = aes(x = time, y = proba_occurence * proba_symptom_cond_syndrome, col = symptom_label)) +
  geom_line() +
  geom_label(data = global_mat[global_mat$flag_max_proba_occurence,], mapping = aes(label = symptom_label), nudge_x = x_offset_heuristic, nudge_y = y_offset_heuristic, size = 3) +
  facet_wrap(.~syndrome_label) +
  theme(legend.position = 'none') +
  ylab(label = "Probabilité") +
  xlab(label = 'Temps, en jours') +
  ggtitle(label = "Probabilité d'occurence des symptômes") +
  scale_colour_hue()

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#           PATIENT PATHWAY SIMULATION
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# !!!!!!!! !!!!!!!! !!!!!!!! !!!!!!!! !!!!!!!! !!!!!!!! !!!!!!!! 
#
#                       important note
#
# - *is_present* indicates that a symptom has developped linked to
# the disease, but it is maybe not yet observed by the patient (e.g. 
# latent symptom, requiring specific medical tests)
# - *is_observed* indicates whether at current day we have confirmation
# of the current presence of the symptom (difference arises in chronic
# cases)
# - *has_been_observed* indicates whether the symptom has ever been observed before :D
#
# !!!!!!!! !!!!!!!! !!!!!!!! !!!!!!!! !!!!!!!! !!!!!!!! !!!!!!!! 

# probability of occurence (truncated gaussian)
proba_occurence <- function(time, mean, sd){
  p1 <- ptruncnorm(q = time-1, a = 0, b = +Inf, mean = mean, sd = sd)
  p2 <- ptruncnorm(q = time, a = 0, b = +Inf, mean = mean, sd = sd)
  return( (p2 - p1) / (1 - p1) )
}

# probability of shift (exponential for now)
proba_shift <- function(time1, time2, mean = 30){
  p1 <- pexp(q = time1, rate = 1 / mean)
  p2 <- pexp(q = time2, rate = 1 / mean)
  return( (p2 - p1) / (1 - p1) )
}

# simulation parameters
NB_DAYS_SIMULATION <- 365
SYNDROME <- syndrome_metadata$id[2]

# collect number of symptoms from metadata
NB_SYMPTOMS <- nrow(symptom_metadata)

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
str(sim)

# identify columns of "sim" dataframe for easy manipulation
colindexes_is_present <- grep(x = colnames(sim), pattern = '_is_present')
colindexes_has_been_observed <- grep(x = colnames(sim), pattern = '_has_been_observed')
colindexes_is_observed <- grep(x = colnames(sim), pattern = '_is_observed')

# identify rows of symptom's metadata for each type
latent_symptoms <- symptom_metadata$type == 'latent'
chronic_symptoms <- symptom_metadata$type == 'chronic'
beginning_to_end_symptoms <- symptom_metadata$type == 'beginning_to_end'

# simulation seed
set.seed(9876543)

# can be randomized, just uncomment the last line ...
symptoms_that_shall_be_present <- rep(FALSE, NB_SYMPTOMS)
symptoms_that_shall_be_present[c(7,4,2)] <- TRUE
### symptoms_that_shall_be_present <- (runif(n = NB_SYMPTOMS) < syndrome_specific_cross_metadata$proba_symptom_cond_syndrome)

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

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#             SIMULATION ANALYSIS
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# syndrome of patient
SYNDROME

# symptoms which will appear
symptom_metadata[symptoms_that_shall_be_present,]

# relevant metadata, sorted by probability of symptom occurence
cross_metadata %>%
  filter(syndrome_id == SYNDROME) %>% 
  arrange(desc(proba_symptom_cond_syndrome))

# preparing data for proper representation
is_present_data <- gather(data = sim[,c(1,colindexes_is_present[symptoms_that_shall_be_present])], key = key, value = value, -time)
is_present_data$type <- 'is_present'
is_observed_data <- gather(data = sim[,c(1,colindexes_is_observed[symptoms_that_shall_be_present])], key = key, value = value, -time)
is_observed_data$type <- 'is_observed'
has_been_observed_data <- gather(data = sim[,c(1,colindexes_has_been_observed[symptoms_that_shall_be_present])], key = key, value = value, -time)
has_been_observed_data$type <- 'has_been_observed'
global_data <- rbind.data.frame(is_present_data, is_observed_data, has_been_observed_data)
global_data$key_uni <- sub(pattern = '\\_is\\_present|\\_is_\\observed|\\_has\\_been\\_observed', replacement = '', x = global_data$key)
global_data$key_uni <- factor(global_data$key_uni)
global_data$type <- factor(global_data$type)
ggplot(data = global_data, mapping = aes(x = time, y = as.numeric(value), col = key_uni)) +
  geom_line() +
  facet_grid(key_uni~type) +
  scale_y_continuous(breaks=c(0, 1),
                     labels=c("No", "Yes")) +
  xlab(label = 'Time, in days') +
  ylab(label = 'Symptom presence and observation, binary') +
  labs(col = 'Symptom') +
  ggtitle(label = 'Presence and observation of symptoms throughout time.')
