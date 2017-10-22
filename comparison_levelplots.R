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

parms <- all_parameters(grps)
all_ind <- parms$all_ind
char_age <- parms$char_age

out <- data_matrix_longform(all_ind)
Mdat <- out$MOME
Fdat <- out$FOME
Mcounts <- out$Mcounts
Fcounts <- out$Fcounts

# all distributions
all_out <- pt_choice_all_choose_variance_model(all_ind, parms$mean_ages,
                                               agevec = grps,
                                               variance_model = "linear")
AP <- list("AP" = list("AP_M" = parms$AP_agemix_M, "AP_F" = parms$AP_agemix_F))
all_out <- append(all_out, AP)

all_mat <- unlist(all_out, recursive = FALSE)

# all plots, for appendix
chsage <- (rep(rep(grps, n_age), 10))
ptage <- (rep(rep(grps, each = n_age), 10))
prob <- as.vector(c(Mdat, Fdat, unlist(all_mat))) # 22 matrices
type_names <- c("Data","gamEmp", "lapEmp", "normEmp", "A-P")
type <- factor(rep(type_names, each = 2*n_age^2), levels = c("Data", "A-P", "gamEmp", "lapEmp", "normEmp")) #two matrices, each with n_age^2 entries
sex <- rep(rep(c("M", "F"), each = n_age^2), 5)

df_all <- data.frame(chsage = chsage, ptage = ptage, prob = prob, type = type, sex = sex)
ggplot(df_all) + geom_raster(aes(x = chsage, y = ptage, fill = prob)) +
    facet_grid(sex ~ type) + 
    scale_fill_gradient2(name = "Prob.", low = "white", high = "black") +
    scale_x_continuous(breaks = seq(10, 60, by = 10), expand = c(0,0)) + 
    scale_y_continuous(breaks = seq(10, 60, by = 10), expand = c(0,0)) + 
    coord_fixed() +
    labs(x="Chooser's Age", y = "Partner's Age")
# ggsave("Plots/AppendixPlots/comparison_levelplots_all.png", width = 10, height = 4, units = "in")

# probability of 60-year-old mixing with 16-27 year old
## ap
sum(subset(df_all, chsage == 60 & type == "AP" & sex == "F" & ptage >= 16 & ptage <= 26)$prob) * 100

## pref
sum(subset(df_all, chsage == 60 & type == "normEmp" & sex == "F" & ptage >= 16 & ptage <= 26)$prob) * 100

## data
sum(subset(df_all, chsage == 60 & type == "Data" & sex == "F" & ptage >= 16 & ptage <= 26)$prob) * 100


# journal article plot, normal only 
df_normal_only <- subset(df_all, chsage > 12 & ptage > 12 & 
                             (type == "lapEmp" | type == "A-P" | type == "Data"))
df_normal_only$type <- as.character(df_normal_only$type)
df_normal_only$type[which(df_normal_only$type == "lapEmp")] <- "Laplace Regression\n (Linear Variance)"
df_normal_only$type[which(df_normal_only$type == "A-P")] <- "Assortative-Proportionate"
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
ggsave("Plots/comparison_levelplots_normal_only.png", height = 6, width = 10, units = "in")

# journal article plot, two age groups.

ggplot(subset(df_normal_only, chsage == 50)) + 
    geom_area(aes(x = ptage, y = prob, fill = sex), alpha = 0.6) + 
    facet_grid(sex ~ type) + 
    scale_fill_manual(name="Sex",
                      values = c("dodgerblue", "darkorange"))+ 
    labs(x = "Partner Age", y = "Probability") + 
    theme_bw(base_size=16)
ggsave("Plots/comparison_age50.png", height=6, width=10, units="in")
