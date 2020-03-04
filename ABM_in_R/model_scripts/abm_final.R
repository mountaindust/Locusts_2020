library(tidyverse)
library(grid)
library(data.table)
library(gganimate)
library(zoo)


# To Do -------------------------------------------------------------------

# make resources gif

# Description -------------------------------------------------------------

# This model simulates the movement of a band of Australian Plague Locust hoppers moving through a field of resources. The model is 1-dimensional, with locusts marching in one direction. Individual locusts are either in a stationary or moving state, and the transition between the two is modeled as a Markov process with transition rates dependent on resource levels at the locust's spatial location. Stationary locusts also consume resources at a rate dependent on resource level- # as locusts deplete the resources at a point in space, their searching becomes less efficient, decreasing the rate at which they feed.

# From these simple rules, a traveling pulse of locusts forms, moving through space in a relatively stationary wave form.

# Simulation Parameters ----------------------------------------------------------

# every t timesteps, n agents, and x spatial location's resources to keep. A value of 1 means no thinning will be done, your results will have every simulation data point. This may exhaust R's memory limit on big runs
t_thin <- 1
n_thin <- 1
x_thin <- 1

# do you want to make a directory to save the results of your simulation run?
# this will include a suite of figures, an RDA file containing all the model parameters, and gifs if you want.
save_results <- T

# Should you make the gifs at the end? Heads up, it will take a while
gifs <- T

# the agents are often piled up at super high densities at the start of a simulation, and this throws the histogram scale off. this will set it such that the y axis of the histogram is scaled to the heights of bars after "hist_plot_after" proportion of the model has run. 0.05 works well
hist_plot_after <- 0.05

# in the histograms, how wide should the bins be? 2 works well
binwidth <- 2

# Setting Model Parameters -------------------------------------------------------

# Probability function parameters

R_plus <- 1.5 # initial resources
N_agents <- 1000
T_end <- 35
delta_t <- 0.1
v <-1
a = 1 #alpha
b = 0 #eta
r_sm = 1 #gamma
c = 0.5 #beta
d = 1.3 #theta
r_ms = 2 #delta
lambda = 0.01 # rate of resource proportion decrease per stationary agent

# #new set
# T_end <- 10000 # How long to run
# delta_t <- 1 # length of a tstep
# N_agents <- 5000
# R_plus <- 200
# v <- 0.03
# lambda <- 5.0e-6
# a = 0.005 #alpha
# b = 0.0025 #eta
# r_sm = 0.01 #gamma
# c = 0.0025 #beta
# d = 0.07 #theta
# r_ms = 0.005 #delta


# Initialize a spatial domain width and delta_x based on delta_t
# this will ensure that during each timestep, a moving locust advances 1 spatial unit
domain_length <- ceiling(T_end/delta_t)
delta_x <- v * delta_t

# number of simulation timesteps that need to be run
T_steps_end <- ceiling(T_end/delta_t)


# Create Directory --------------------------------------------------------

if(save_results == T){
  
  # create a directory named by this unique set of parameters. If it already exists, that directory will just be used and any plots/gifs stored there will be overwritten
sim_name <- paste0("simulation_Nagents=", N_agents, "_Rplus=", R_plus, "_Tend=", T_end,
                   "_deltaT=", delta_t, "_v=", v, "_alpha=", a, "_eta=", b, 
                   "_gamma=", r_sm, "_beta=", c, "_theta=", d, 
                   "_delta=", r_ms, "_lambda=", lambda, 
                   "_tthin=", t_thin, "_nthin=", n_thin, "_xthin=", x_thin)

sim_dir <- paste0("simulations/", sim_name)

if(!dir.exists("simulations")){
  dir.create("simulations")
}

if(!dir.exists(sim_dir)){
  dir.create(sim_dir)
}

save(file = paste0(sim_dir, "/params.rda"), list = c("N_agents", "R_plus", "T_end", "delta_t", "v", "a", "b", "r_sm", "c", "d", "r_ms", "lambda", "t_thin", "n_thin", "x_thin"))

}

# Transition + Resource Functions ------------------------------------------------------------

# Probability function of going from S -> M per delta_t
p_SM = function(R) {
  (b - (b-a) * exp(-r_sm * R)) * delta_t
}

# Probability function of going from M -> S per delta_t
p_MS = function(R){
  (d - (d-c) * exp(-r_ms * R)) * delta_t
}

# Exp function of S to reduce R by. R(t+1) = R*R_exp(S)
# per unit time in a given location
R_exp = function(S){
  exp(-lambda * S * (delta_t/delta_x)) #  needs to depend on delta_t and delta_x!
}

# check what the curves look like

# generate sequence of resources to check the transition rates along, then plot
test_resources <- seq(from = 0, to = R_plus, length = 100)

tibble(resources = seq(from = 0, to = R_plus, length = 100), 
       psm = p_SM(resources),
       pms = p_MS(resources)) %>% 
  gather(key = type, value = value, -resources) %>% 
  ggplot(aes(x = resources, y = value, group = type)) +
  geom_line(aes(color = type), size = 1.5) +
  theme_minimal() +
  ylab("probability of transition per delta_t") +
  scale_color_viridis_d(labels = c("moving to stationary", "stationary to moving"))

if(save_results == T){ggsave(paste0(sim_dir, "/transitions.jpg"))}

# now generate sequence of stationary locust counts, and check resource reduction

tibble(num_stationary = 1:1000, r_exp_res = R_exp(1:1000)) %>% 
  ggplot(aes(x = num_stationary, y = r_exp_res)) +
  geom_line() +
  theme_minimal() +
  ylab("proportion of resources remaining after 1 timestep") +
  xlab("number of stationary locusts at a location")

if(save_results == T){ggsave(paste0(sim_dir, "/locusts_feed.jpg"))}

# now we're looking at how a single locust reduces the resources through time
tibble(tstep = 1:1000, resources = R_plus*(R_exp(1)^tstep)) %>% 
  ggplot(aes(x = tstep, y = resources)) +
  geom_line() +
  theme_minimal() +
  xlab("time") + ylab("resources with 1 stationary locust")

if(save_results == T){ggsave(paste0(sim_dir, "/one_locust_feeds.jpg"))}




# Initialization -----------------------------------------------------

agents <- matrix(NA, N_agents, 2)
resources <- rep(R_plus, T_steps_end)

save_agents <- array(data = NA, dim = c(T_steps_end/t_thin, N_agents/n_thin, 2))

save_resources <- matrix(NA, T_steps_end/t_thin, domain_length/x_thin)

# start all locusts in the first 5 spatial points, equally spaced, stationary
ifelse(N_agents >=5,
       agents[,1] <- rep(1:5, out = N_agents),
       agents[,1] <- rep(1, out = N_agents)
)
# set the initial states to 0
agents[,2] <- 0 

agent_save_ids <- sample(1:N_agents, N_agents/n_thin, replace = F)
x_resource_save_ids <- seq(from = x_thin, by = x_thin, length.out = domain_length/x_thin)


# Main Loop -----------------------------------------------------

# start at 2nd timestep since we've initialized with the first one
for (tstep in 2:(T_steps_end)){
  
  ### Step 1: State Change
  
  #  this returns a logical: are you moving?
  are_you_moving <- agents[,2] == 1
  
  # returns a logical for each agent: do you transition to stationary? your probability of doing this is based on a) are you moving to begin with and b) the pMS for the resources where you are
  do_MS <- rbinom(N_agents, 1, (are_you_moving * p_MS(resources[agents[,1]])))
  
  # this returns a logical for every agent: are you stationary?
  are_you_stationary <- agents[,2] == 0
  
  # returns a logical for each agent: do you transition to moving? your probability of doing this is based on a) are you stationary to begin with and b) the pMS for the resources where you are
  do_SM <- rbinom(N_agents, 1, (are_you_stationary * p_SM(resources[agents[,1]])))
  
  # now, given the S to M and M to S instructions, change state to moving or not
  agents[,2] <- agents[,2] - do_MS + do_SM
  
  ### Step 2: Move
  
  # move the moving locusts by an integer, which represents v*delta_x 
  agents[,1] <- agents[,1] + agents[,2]
  
  ### Step 3: Eat
  
  # get x-values of stationary agents, and tabulate into unique values (first column) and counts (second column)
  stationary_agents_locations <- agents[agents[,2] == 0, 1]
  
  counts_at_x <- as.data.frame(table(stationary_agents_locations))
  
  # have to do this because table() generates factors
  counts_at_x$stationary_agents_locations <- as.numeric(as.character(counts_at_x$stationary_agents_locations))
  
  # for each spatial location, calculate the amount of resource reduction to be done and store it in the dataframe
  counts_at_x$r_exp_result <- R_exp(counts_at_x$Freq)
  
  # now do the reduction at the locations where there are currently stationary locusts
  resources[counts_at_x$stationary_agents_locations] <- resources[counts_at_x$stationary_agents_locations] * counts_at_x$r_exp_result

  ### 4: Storage
  if (tstep %% t_thin == 0) {
    save_agents[tstep/t_thin,,] <- agents[agent_save_ids,]
    save_resources[tstep/t_thin,] <- resources[x_resource_save_ids]
  }
}



# Tidy Data -----------------------------------------------------

result <- as.data.table(save_agents) %>% 
  as_tibble() %>% 
  spread(key = V3, value = value) %>% 
  rename(tstep = V1, agent_id = V2, location = `1`, state = `2`) %>% 
  mutate(tstep = tstep*t_thin)

rownames(save_resources) <- seq(from = t_thin, by = t_thin, length.out = T_steps_end/t_thin)

resources <- as.data.table(save_resources, keep.rownames = T) %>% 
  as_tibble() %>% 
  gather(key = x, value = value, -rn) %>% 
  rename(tstep = rn) %>% 
  mutate(x = as.numeric(str_remove_all(x, "V")) * x_thin, tstep = as.numeric(tstep))





# Plot --------------------------------------------------------------------

# a histogram looking at an early time, a quarter through the run, halfway through the run, and at the end of the run. Sort of a proxy for a gif that shows the traveling pulse

plot_steps <- quantile(unique(result$tstep), 
                       probs = seq(from = 0, to = 1, length.out = 4), 
                       na.rm = F, names = T, type = 3)

# we calculate the histogram bar heights after the model has run "hist_plot_after" % of the way through, with hist_plot_after being set up above in the Simulation parameters
counts <- result %>% 
  filter(tstep > hist_plot_after*T_steps_end) %>% 
  group_by(tstep, location) %>% 
  count()

# then we find the tallest one
max_height <- max(counts$n)

result %>% 
  filter(tstep %in% plot_steps) %>% 
  ggplot(aes(x = location, fill = as.factor(tstep))) +
  geom_histogram(binwidth = 2, position = "identity", alpha = 0.5) +
  theme_minimal() +
  scale_fill_viridis_d() +
  coord_cartesian(ylim = c(0,max_height*1.5))

if(save_results == T){ggsave(paste0(sim_dir, "/result_hist.jpg"), width = 8, height = 6)}


# if the x axis is super spread out, this will make a more condensed plot
plot_steps <- quantile(unique(result$tstep), 
                       probs = seq(from = 0, to = 1, length.out = 51), 
                       na.rm = F, names = T, type = 3)
plot_steps <- unique(plot_steps[47:51])

result %>% 
  filter(tstep %in% plot_steps) %>% 
  ggplot(aes(x = location, fill = as.factor(tstep))) +
  geom_histogram(binwidth = 2, position = "identity", alpha = 0.5) +
  theme_minimal() +
  scale_fill_viridis_d() +
  coord_cartesian(ylim = c(0,max_height*1.5))

if(save_results == T){ggsave(paste0(sim_dir, "/result_hist_endrun.jpg"), width = 8, height = 6)}


# plot resources at 100 (roughly) evenly spaced time steps

plot_steps <- quantile(unique(result$tstep), 
                       probs = seq(from = 0, to = 1, length.out = 100), 
                       na.rm = F, names = T, type = 3)
resources %>% 
  filter(tstep %in% plot_steps) %>% 
  group_by(x) %>% 
  mutate(min_val = min(value)) %>% 
  ungroup() %>% 
  filter(min_val != R_plus) %>% 
  ggplot(aes(x = x, y = value, group = tstep)) +
  geom_line(aes(color = tstep), alpha = 0.8) +
  scale_color_viridis_c() +
  theme_minimal() +
  ylab("resources") + xlab("location")

if(save_results == T){ggsave(paste0(sim_dir, "/result_resources.jpg"), width = 8, height = 6)}

# plot resources at 8 (roughly) evenly spaced time steps, and facet wrap them, which is useful if they're bunched up

plot_steps <- quantile(unique(result$tstep), 
                       probs = seq(from = 0, to = 1, length.out = 8), 
                       na.rm = F, names = T, type = 3)
resources %>% 
  filter(tstep %in% plot_steps) %>% 
  group_by(x) %>% 
  mutate(min_val = min(value)) %>% 
  ungroup() %>% 
  filter(min_val != R_plus) %>% 
  ggplot(aes(x = x, y = value, group = tstep)) +
  geom_line(aes(color = tstep)) +
  scale_color_viridis_c() +
  theme_minimal() +
  facet_wrap(~tstep) +
  ylab("resources") + xlab("location")

if(save_results == T){ggsave(paste0(sim_dir, "/result_resources_facet.jpg"), width = 8, height = 6)}


# same thing, but at the end if the x axis is way spread out

plot_steps <- quantile(unique(result$tstep), 
                       probs = seq(from = 0, to = 1, length.out = 101), 
                       na.rm = F, names = T, type = 3)

plot_steps <- unique(plot_steps[97:101])


resources %>% 
  filter(tstep %in% plot_steps) %>% 
  group_by(x) %>% 
  mutate(min_val = min(value), max_val = max(value)) %>% 
  ungroup() %>% 
  filter(min_val != R_plus & max_val > 0.01*R_plus) %>% 
  ggplot(aes(x = x, y = value, group = tstep)) +
  geom_line(aes(color = tstep)) +
  scale_color_viridis_c() +
  theme_minimal() +
  facet_wrap(~tstep) +
  ylab("resources") + xlab("location")

if(save_results == T){ggsave(paste0(sim_dir, "/result_resources_facet_endrun.jpg"), width = 8, height = 6)}


pick <- function(condition){ function(d) d %>% filter(!!enquo(condition))}

# now let's plot all of the agents to look like a band of individual locusts, at 90% of the way through the simulation

result %>% 
  group_by(tstep) %>% 
  mutate(mean_location = mean(location)) %>% 
  ungroup() %>% 
  filter(tstep == plot_steps[1]) %>% 
  ggplot(aes(x = location, y = agent_id)) +
  geom_vline(data = pick(agent_id == 1), aes(xintercept = mean_location), color = "grey80") +
  geom_point(alpha = 0.6, aes(color = as.factor(state))) +
  scale_color_viridis_d(begin = 0, end = 0.8, labels = c("stationary", "moving"), name = "state") +
  theme_minimal() +
  ggtitle(paste("time =", plot_steps[1]))

if(save_results == T){ggsave(paste0(sim_dir, "/result_points.jpg"), width = 8, height = 6)}

result %>% 
  group_by(tstep) %>% 
  summarise(stand_dev = sd(location)) %>% 
  ggplot(aes(x = tstep, y = stand_dev)) +
  geom_line() +
  theme_minimal() +
  ylab("standard deviation of locust location")

if(save_results == T){ggsave(paste0(sim_dir, "/result_stdev.jpg"), width = 8, height = 6)}

# calculate the velocity of the mean locust location through time
vel_calc <- result %>% 
  group_by(tstep) %>% 
  summarise(mean_loc = mean(location))

vel_calc <- vel_calc %>% 
  mutate(velocity = rollapply(data = .$mean_loc, width = 2, FUN = diff, fill = NA, align = "right")/.$tstep[1])

vel_calc %>% 
  ggplot(aes(x = tstep, y = velocity)) +
  geom_line() +
  theme_minimal()

if(save_results == T){ggsave(paste0(sim_dir, "/result_velocity.jpg"), width = 8, height = 6)}


# Gifs --------------------------------------------------------------------

# we will make some gifs to visualize the locust movement

if (gifs == T) {
  
  a <- result %>% 
    ggplot(aes(x = location)) +
    labs(title = 'Time: {frame_time}') +
    theme(axis.text=element_text(size=14),
          axis.text.x = element_blank(), 
          axis.title.x = element_blank(), 
          axis.ticks.x = element_blank()) +
    transition_time(tstep) +
    ease_aes('linear') +
    theme_minimal()
  
  gif_hist <- animate(a + geom_histogram(binwidth = binwidth) +
                        coord_cartesian(ylim = c(0,max_height*3)), 
                      nframes = min(c(length(unique(result$tstep)), 500)), 
                      width = 800, height = 400)
  
  gif_hist
  
  if(save_results == T){anim_save(paste0(sim_dir, "/hist.gif"))}
  
  gif_dens <- animate(a + stat_density(geom="line") +
                        coord_cartesian(ylim = c(0,max_height/N_agents*1.5)), 
                      nframes =  min(c(length(unique(result$tstep)), 500)), 
                      width = 800, height = 400)
  gif_dens
  
  if(save_results == T){anim_save(paste0(sim_dir, "/dens.gif"))}
  
  # this function lets you pick a specific subset of the data to use with a single geom
  pick <- function(condition){ function(d) d %>% filter(!!enquo(condition))}
  
  # now let's plot a bunch of the agents to look like a band of individuals
  # heads up, this will look pretty bad on longer, thinned runs. Much easier to see what's happening on shorter, thinned runs
  
  b <- result %>% 
    group_by(tstep) %>% 
    mutate(mean_location = mean(location)) %>% 
    ungroup() %>% 
    filter(agent_id %in% 1:1000) %>% 
    ggplot(aes(x = location, y = agent_id)) +
    geom_vline(data = pick(agent_id == 1), aes(xintercept = mean_location), color = "grey80") +
    geom_point(alpha = 0.6, aes(color = as.factor(state))) +
    scale_color_viridis_d(begin = 0, end = 0.8) +
    theme(axis.text=element_text(size=20)) +
    labs(title = paste('Time: {frame_time}')) +
    transition_time(tstep) + 
    theme_minimal()
  
  gif_points <- animate(b, nframes =  round(min(c(length(unique(result$tstep)), 500))/2), width = 800, height = 400, res = 100)
  gif_points
  
  if(save_results == T){anim_save(paste0(sim_dir, "/points.gif"))}
  
  a <- resources %>% 
    ggplot(aes(x = x, y = value, group = tstep)) +
    geom_line() +
    theme_minimal() +
    ylab("resources") + xlab("location") +
    labs(title = 'Time: {frame_time}') +
    theme(axis.text=element_text(size=14),
          axis.text.x = element_blank(), 
          axis.title.x = element_blank(), 
          axis.ticks.x = element_blank()) +
    transition_time(tstep) +
    ease_aes('linear') +
    theme_minimal() +
    xlim(0, max(result$location))
  
  gif_resource <- animate(a, 
                      nframes = min(c(length(unique(resources$tstep)), 500)), 
                      width = 800, height = 400)
  
  gif_resource
  
  if(save_results == T){anim_save(paste0(sim_dir, "/resource.gif"))}
}
