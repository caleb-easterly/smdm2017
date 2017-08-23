require(scales)
require(RColorBrewer)
require(lattice)
require(gridExtra)
require(ggplot2)
require(sfceHPV)
require(extraDistr)

# use age_specific_rates to get the mixing matrices
grps <- seq(12, 60, by = 2)
n_age <- length(grps)

group_data <- all_parameters(grps)
all_ind <- group_data$all_ind

out <- data_matrix_longform(all_ind)
Mdat <- out$MOME
Fdat <- out$FOME
Mcounts <- out$Mcounts
Fcounts <- out$Fcounts

# main Likelihood stuff is now in loglikes_varmod_comp.R
# set up A-P mixing matrix
M_APmat <- F_APmat <- matrix(0, nrow = n_age, ncol = n_age)
rho_temp <- lambda_all(group_data$init_vec, group_data, 2, 0)
rhoM_temp <- rho_temp$rhoM
rhoF_temp <- rho_temp$rhoF
j <- 1
for (i in seq(3, n_age*3, by = 3)){
  M_APmat[j, ] <- colSums(rhoM_temp[[i]])
  F_APmat[j, ] <- colSums(rhoF_temp[[i]])
  j <- j + 1
}
AP <- list(AP= list(m_AP = M_APmat, f_AP = F_APmat))

# all distributions
all_out <- pt_choice_all_choose_variance_model(all_ind, group_data$mean_ages,
                                               agevec = grps,
                                               variance_model = "sqrt")
all_out <- append(all_out, AP)

all_mat <- unlist(all_out, recursive = FALSE)

# create data frame for plotting
chsage <- (rep(rep(grps, n_age), 16))
ptage <- (rep(rep(grps, each = n_age), 16))
prob <- as.vector(c(Mdat, Fdat, unlist(all_mat))) # 22 matrices
type_names <- c("Data","gamEmp", "gamConst", "lapEmp", "lapConst", "normEmp", "normConst", "A-P")
type <- factor(rep(type_names, each = 2*n_age^2), levels = c("Data", "A-P", "gamEmp", "gamConst", "lapEmp", "lapConst", "normEmp", "normConst")) #two matrices, each with n_age^2 entries
sex <- rep(rep(c("M", "F"), each = n_age^2), 8)

df_all <- data.frame(chsage = chsage, ptage = ptage, prob = prob, type = type, sex = sex)

# probability of 60-year-old mixing with 16-26 year old
## ap
sum(subset(df_all, chsage == 60 & type == "A-P" & sex == "F" & ptage >= 16 & ptage <= 26)$prob) * 100

## pref
sum(subset(df_all, chsage == 60 & type == "normEmp" & sex == "F" & ptage >= 16 & ptage <= 26)$prob) * 100

## data
sum(subset(df_all, chsage == 60 & type == "Data" & sex == "F" & ptage >= 16 & ptage <= 26)$prob) * 100

# journal article plot, normal only 
df_normal_only <- subset(df_all, chsage > 12 & ptage > 12 & 
                        (type == "normEmp" | type == "A-P" | type == "Data"))
df_normal_only$type <- as.character(df_normal_only$type)
df_normal_only$type[which(df_normal_only$type == "normEmp")] <- "Gaussian Regression\n (Square Root Variance)"
df_normal_only$type[which(df_normal_only$type == "AP")] <- "Assortative-Proportionate"
df_normal_only$type[which(df_normal_only$type == "Data")] <- "Natsal-3 Data"
df_normal_only$sex <- c("Male", "Female")[(df_normal_only$sex == "F")*1 + 1]
df_normal_only$type <- as.factor(df_normal_only$type)
ggplot(df_normal_only) + 
  geom_raster(aes(x = chsage, y = ptage, fill = prob)) +
  facet_grid(sex ~ type) + 
  scale_fill_gradient2(name = "Probability", low = "white", high = "black") +
  scale_x_continuous(breaks = seq(10, 60, by = 10), expand = c(0,0)) + 
  scale_y_continuous(breaks = seq(10, 60, by = 10), expand = c(0,0)) + 
  coord_fixed() +
  labs(x="Chooser's Age", y = "Partner's Age") + 
  theme_bw(base_size = 14)
ggsave("Plots/comparison_levelplots_normal_only.tiff")


