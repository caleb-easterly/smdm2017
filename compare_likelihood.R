library(sfceHPV)
grps <- c(12, 13, seq(14, 60, by = 2))
group_data <- all_parameters(grps)

# A-P 
all_out_AP <- list("AP" = list(m_AP = group_data$AP_agemix_M, f_AP = group_data$AP_agemix_F))

all_ind <- group_data$all_ind
age_prop_full <- group_data$age_prop_single_year

# linear
all_out_linear <- pt_choice_all_choose_variance_model(all_ind, group_data$mean_ages, grps, "linear")
all_out_sqrt <- pt_choice_all_choose_variance_model(all_ind, group_data$mean_ages, grps, "sqrt")
all_out_log <- pt_choice_all_choose_variance_model(all_ind, group_data$mean_ages, grps, "log")

loglikes_AP <- sapply(all_out_AP, mixing_matrix_loglikelihood, group_data)  
loglikes_linear <- sapply(all_out_linear, mixing_matrix_loglikelihood, group_data)
loglikes_sqrt <- sapply(all_out_sqrt, mixing_matrix_loglikelihood, group_data)
loglikes_log <- sapply(all_out_log, mixing_matrix_loglikelihood, group_data)
all_loglikes_df <- data.frame("VarianceModel" = c(rep(c("Linear", "Constant"), 3),
                                                  rep(c("Square Root", "Constant"), 3),
                                                  rep(c("Log", "Constant"), 3)), 
                              "Distribution" = rep(c("Gamma", "Laplace", "Normal"), 3),
                              "RegOrSh" = rep(c("Regression", "Shared"), 9),
                              "NegLogLike" = -c(loglikes_linear, loglikes_sqrt, loglikes_log))
require(ggplot2)
ggplot(all_loglikes_df) + 
    # scale_y_continuous(limits = c(45000, 65000)) + 
    geom_point(aes(x = RegOrSh, y = NegLogLike, color = VarianceModel), 
               stat = "identity", size = 3) + 
    labs(x = "Method for Mean Partner Age",
         y = "Negative Log Likelihood") +
    facet_grid(.~Distribution) + 
    geom_hline(aes(yintercept = -loglikes_AP), linetype = 2, size = 1,
               color = "blue") + 
    scale_y_continuous(breaks = seq(42000, 56000, by = 2000)) + 
    theme_bw(base_size = 12)
ggsave("Plots/negLogLikelihoodAll.png",
       height = 4, width = 6.5, units = "in")

