#' Compare extending vaccination under different models

# age groups
library(sfceHPV)
library(ggplot2)
library(extraDistr)
library(deSolve)
library(rootSolve)
library(beepr)
library(microbenchmark)

analyze_extension <- function(sigma = 1/100, comparison = c("male_catchup75", "no_male_catchup75")) {
    
    agevec <- c(12, 13, seq(14, 38, by = 2), c(40, 45, 50, 55, 60))
    parms <- all_parameters(agevec, sigma = sigma) # middle scenario from Choi 2010
    init_vec <- parms$init_vec
    
    define_journal_article_vaccination_strategies(agevec)
    
    comparison <- match.arg(comparison)
    if (comparison == "male_catchup75"){
        vacc_strat_base <- JAstratG_base
        vacc_strat_comp <- JAstratG_comp
    } else if (comparison == "no_male_catchup75") {
        vacc_strat_base <- JAstratH_base
        vacc_strat_comp <- JAstratH_comp
    }
    
    #' base case: female vaccination, with catchup to 26
    prev_base_case <- list("emp" = estimate_steady_state(parms, "pref", vacc_strategy = vacc_strat_base),
                           "ap" = estimate_steady_state(parms, "fact", vacc_strategy = vacc_strat_base))
    prev_extra_catch <- list("emp" = estimate_steady_state(parms, "pref", vacc_strategy = vacc_strat_comp),
                             "ap" = estimate_steady_state(parms, "fact", vacc_strategy = vacc_strat_comp))
    
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
                                        Vacc = coverage,
                                        Comparison = comparison,
                                        Sigma = paste("Sigma =", as.character(1/sigma), "years"))
    
    
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
    
    return(list("benefit_ratio" = benefit_ratio, "df" = relred_extra_catch_df))
}

# with male catchup
sigma100yr <- analyze_extension(sigma = 1/100,  
                                         "male_catchup75")
sigma10yr <- analyze_extension(sigma = 1/10, 
                                        "male_catchup75")

# without male catchup
sigma100yr_noMcatch <- analyze_extension(sigma = 1/100,
                                                  "no_male_catchup75")
sigma10yr_noMcatch <- analyze_extension(sigma = 1/10,
                                                 "no_male_catchup75")

benefit_ratios <- c("sigma100yr" = sigma100yr$benefit_ratio, 
                    "sigma10yr" = sigma10yr$benefit_ratio,
                    "sigma100yr_noMcatch" = sigma100yr_noMcatch$benefit_ratio,
                    "sigma10yr_noMcatch" = sigma100yr_noMcatch$benefit_ratio)

dfs <- rbind(sigma100yr$df, 
             sigma10yr$df,
             sigma100yr_noMcatch$df,
             sigma10yr_noMcatch$df)

ggplot(subset(dfs, Age >= 14)) +
    geom_line(aes(x = Age,
                  y = 100*relred_extra_catch,
                  color = Mixing,
                  linetype = Comparison), size = 1.2) +
    theme_bw(base_size = 14) +
    facet_grid(Sigma~Sex) +
    scale_x_continuous(breaks = seq(10, 60, by = 5)) +
    scale_y_continuous(limits = c(0, 15)) + 
    labs(y = expression(R[list(k, s, a)]))
ggsave("Plots/relative_reduction_extending_coverage.png")
