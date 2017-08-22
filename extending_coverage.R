#' Compare extending vaccination under different models

# age groups
library(sfceHPV)
library(ggplot2)
library(extraDistr)
library(deSolve)
library(rootSolve)
library(beepr)
library(microbenchmark)

analyze_extension <- function(sigma = 1/100, plot_name = "Plots/temp.png") {
    
    agevec <- c(12, 13, seq(14, 38, by = 2), c(40, 45, 50, 55, 60))
    parms <- all_parameters(agevec, sigma = sigma) # middle scenario from Choi 2010
    init_vec <- parms$init_vec
    
    define_journal_article_vaccination_strategies(agevec)
    
    #' base case: female vaccination, with catchup to 26
    prev_base_case <- list("emp" = estimate_steady_state(parms, "pref", vacc_strategy = JAstratG_base),
                           "ap" = estimate_steady_state(parms, "fact", vacc_strategy = JAstratG_base))
    prev_extra_catch <- list("emp" = estimate_steady_state(parms, "pref", vacc_strategy = JAstratG_comp),
                             "ap" = estimate_steady_state(parms, "fact", vacc_strategy = JAstratG_comp))
    
    #' make data frame for plotting
    relred_extra_catch <- c(
        (prev_base_case$emp - prev_extra_catch$emp)/prev_base_case$emp,
        (prev_base_case$ap - prev_extra_catch$ap)/prev_base_case$ap
    )
    
    sex <- rep(rep(c("M", "F"), parms$n_age), 2)
    age <- rep(rep(agevec, each = 2), 2)
    mixing <- rep(c("Emp.", "A-P"), each = parms$n_age*2)
    coverage <- rep("75% Strategy", parms$n_age * 2)
    
    relred_extra_catch_df <- data.frame(relred_extra_catch = relred_extra_catch,
                                        Sex = sex,
                                        Age = age,
                                        Mixing = mixing,
                                        Vacc = coverage)
    
    older_ind <- which(agevec >= 40)
    avg_red_extra_catch_100yr_old <- with(parms, 
                                          c("M_ap" = t(age_prop[older_ind]) %*% 
                                                subset(relred_extra_catch_df, mixing == "A-P" &
                                                           sex == "M" & Age >= 40)$relred_extra_catch / 
                                                sum(age_prop[older_ind]),
                                            "F_ap" = t(age_prop[older_ind]) %*% 
                                                subset(relred_extra_catch_df, mixing == "A-P" &
                                                           sex == "F" & Age >= 40)$relred_extra_catch / 
                                                sum(age_prop[older_ind]),
                                            "M_emp" = t(age_prop[older_ind]) %*% 
                                                subset(relred_extra_catch_df, mixing == "Emp." &
                                                           sex == "M" & Age >= 40)$relred_extra_catch / 
                                                sum(age_prop[older_ind]),
                                            "F_emp" = t(age_prop[older_ind]) %*% 
                                                subset(relred_extra_catch_df, mixing == "Emp." &
                                                           sex == "F" & Age >= 40)$relred_extra_catch / 
                                                sum(age_prop[older_ind])
                                          )
    )
    # estimate average benefit ratio for older 
    benefit_ratio <- avg_red_extra_catch_100yr_old['F_ap'] / avg_red_extra_catch_100yr_old['F_emp']
    
    ggplot(subset(relred_extra_catch_df, Age >= 14)) +
        geom_line(aes(x = Age, y = 100*relred_extra_catch, color = Mixing), size = 1.2) +
        theme_bw(base_size = 14) +
        facet_grid(.~Sex) +
        scale_x_continuous(breaks = seq(10, 60, by = 5)) +
        labs(y = expression(APRAV[list(s, a)]),
             title =
                 paste("Catch-up Vaccination of 26-40 y/o Females\nsigma = ", sigma))
    ggsave(plot_name)
    
    return(benefit_ratio)
}

benefit_ratio_100yr <- analyze_extension(sigma = 1/100, plot_name = "Plots/100yr_vacc.png")
benefit_ratio_10yr <- analyze_extension(sigma = 1/10, plot_name = "Plots/10yr_vacc.png")

