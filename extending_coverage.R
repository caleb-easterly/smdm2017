#' Compare extending vaccination under different models

# age groups
library(sfceHPV)
library(ggplot2)

analyze_extension <- function(sigma = 1/100, comparison = c("male_catchup75", "no_male_catchup75")) {
    
    agevec <- c(12, 13, seq(14, 60, by = 2))
    parms <- all_parameters(agevec, sigma = sigma, variance_model = "linear") # middle scenario from Choi 2010
    init_vec <- parms$init_vec
    age_prop <- parms$age_prop
    comparison <- match.arg(comparison)
    if (comparison == "male_catchup75"){
        vacc_strat_base <- parms$JAstratG_base
        vacc_strat_comp <- parms$JAstratG_comp
    } else if (comparison == "no_male_catchup75") {
        vacc_strat_base <- parms$JAstratH_base
        vacc_strat_comp <- parms$JAstratH_comp
    }
    
    #' base case: female vaccination, with catchup to 26
    prev_base_case <- list("emp" = estimate_steady_state(parms, "emp", vacc_strategy = vacc_strat_base),
                           "ap" = estimate_steady_state(parms, "ap", vacc_strategy = vacc_strat_base))
    prev_extra_catch <- list("emp" = estimate_steady_state(parms, "emp", vacc_strategy = vacc_strat_comp),
                             "ap" = estimate_steady_state(parms, "ap", vacc_strategy = vacc_strat_comp))
    
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
              
    
    younger_ind <- which(agevec <= 25 & agevec >= 14)
    avg_red_extra_catch_100yr_young <- with(parms, 
                                          c("M_ap" = t(age_prop[younger_ind]) %*% 
                                              subset(relred_extra_catch_df, mixing == "A-P" &
                                                       sex == "M" & Age <= 25 & Age >= 14)$relred_extra_catch / 
                                              sum(age_prop[younger_ind]),
                                            "F_ap" = t(age_prop[younger_ind]) %*% 
                                              subset(relred_extra_catch_df, mixing == "A-P" &
                                                       sex == "F" & Age <= 25 & Age >= 14)$relred_extra_catch / 
                                              sum(age_prop[younger_ind]),
                                            "M_emp" = t(age_prop[younger_ind]) %*% 
                                              subset(relred_extra_catch_df, mixing == "Emp." &
                                                       sex == "M" & Age <= 25 & Age >= 14)$relred_extra_catch / 
                                              sum(age_prop[younger_ind]),
                                            "F_emp" = t(age_prop[younger_ind]) %*% 
                                              subset(relred_extra_catch_df, mixing == "Emp." &
                                                       sex == "F" & Age <= 25 & Age >= 14)$relred_extra_catch / 
                                              sum(age_prop[younger_ind])
                                          )
    )
    # estimate average benefit ratio for older 
    benefit_ratio_young <- c(avg_red_extra_catch_100yr_young['M_ap'] / avg_red_extra_catch_100yr_young['M_emp'],
                             avg_red_extra_catch_100yr_young['F_ap'] / avg_red_extra_catch_100yr_young['F_emp'])
    benefit_ratio_old <- c(avg_red_extra_catch_100yr_old['M_ap'] / avg_red_extra_catch_100yr_old['M_emp'],
                           avg_red_extra_catch_100yr_old['F_ap'] / avg_red_extra_catch_100yr_old['F_emp'])
    
    benefit_ratios <- data.frame(
      "Sex" = c("M", "F"),
      "<=25" = benefit_ratio_young * 100, # convert to percent
      ">40" = benefit_ratio_old * 100
    )
    colnames(benefit_ratios) <- c("Sex","<26 y/o", ">41 y/o ")
    rownames(benefit_ratios) <- NULL
    return(list("benefit_ratios" = benefit_ratios,
                "df" = relred_extra_catch_df))
}

# with male catchup
sigma100yr <- analyze_extension(sigma = 1/100,  
                                         "male_catchup75")

# without male catchup
sigma100yr_noMcatch <- analyze_extension(sigma = 1/100,
                                                  "no_male_catchup75")

benefit_ratios <- sigma100yr$benefit_ratios


 # write to csv for using in poster
write.csv(format(benefit_ratios, nsmall = 1, digits = 3),
          file = "Poster/reduction_table.csv", quote = FALSE,
          row.names = FALSE)

dfs <- rbind(sigma100yr$df, 
             sigma100yr_noMcatch$df)
dfs$Comparison <- c("Male Catchup\n75% uptake", "No Male Catchup\n75% uptake")[as.numeric(dfs$Comparison)]
    
ggplot(subset(dfs, Age >= 14)) +
    geom_line(aes(x = Age,
                  y = 100*relred_extra_catch,
                  color = Mixing), size = 1.2) +
    theme_bw(base_size = 18) +
    facet_grid(Comparison~Sex) +
    scale_x_continuous(breaks = seq(10, 60, by = 5)) +
    scale_y_continuous(limits = c(0, 15)) + 
    labs(y = "Percent Reduction in Prevalence")
ggsave("Plots/relative_reduction_extending_coverage.png", width = 10, height=6, units ="in")
