#-----------------------------------------
# Open Libraries
#-----------------------------------------
library(tidyverse)
library(openxlsx)
#library(gridExtra)
#-----------------------------------------

#Functions

#---------------------------
# Function to add back-calculated values to Caltable
addDeviations <- function (fit, caltable) {
  caltable %>%
    mutate(BackCalc_Conc = 10^ ((Area - fit$coef[[1]])/fit$coef[[2]]),
           Deviation = (signif(100*(BackCalc_Conc-Conc_orig)/Conc_orig,3)))
}

#---------------------------


# Function to optimize the cal table by iteratively removing the most deviating calpoint
OptimizeCaltable2 <- function (fit, caltable) {
  devs <- addDeviations(fit, caltable)
  print("Initial Deviations:")
  print(devs)
  
  if (max(abs(devs$Deviation), na.rm = TRUE) <= 20) {
    print("All deviations within 20%. Returning original table.")
    
    final_table <- left_join(caltable, devs) %>%
      mutate(excluded = is.na(Deviation),,
             slope = coef(fit)[2],
             intercept = coef(fit)[1],
             adj_R_squared = summary(fit)$adj.r.squared)
    
    # Add LLOQ and ULOQ info
    lowest_included_conc <- final_table %>%
      filter(!excluded) %>%
      summarise(LLOQ = min(Conc_orig, na.rm = TRUE)) %>%
      pull(LLOQ)
    
    highest_included_conc <- final_table %>%
      filter(!excluded) %>%
      summarise(ULOQ = max(Conc_orig, na.rm = TRUE)) %>%
      pull(ULOQ)
    
    # Calculate % of included points
    percent_included <- final_table %>%
      summarise(percent_included = mean(!excluded) * 100) %>%
      pull(percent_included)
    
    # Add to final_table
    final_table <- final_table %>%
      mutate(LLOQ = lowest_included_conc,
             ULOQ = highest_included_conc,
             percent_included = round(percent_included, 1))
    
    final_fit <- list(
      Intercept = coef(fit)[1],
      Slope = coef(fit)[2],
      Adj.R.Squared = summary(fit)$adj.r.squared
    )
    
    return(final_table)
  }
  
  iteration <- 1
  repeat {
    max_dev <- max(abs(devs$Deviation))
    worst_row <- which.max(abs(devs$Deviation))
    message(sprintf("Iteration %d: Max deviation = %.2f%% (removing row %d)", 
                    iteration, max_dev, worst_row))
    
    revised <- devs[-worst_row, ]  # Remove the row with the largest deviation
    lmcalcurve_rev <- lm(Area ~ Conc, data = revised)
    devs <- addDeviations(lmcalcurve_rev, revised)
    
    print(sprintf("Updated coefficients: Intercept = %.4f, Slope = %.4f",
                  coef(lmcalcurve_rev)[1], coef(lmcalcurve_rev)[2]))
    
    if (max(abs(devs$Deviation), na.rm = TRUE) <= 20) {
      message("All deviations now within 20%. Optimization complete.")
      break
    }
    
    iteration <- iteration + 1
  }
  
  # Final output
  final_table <- left_join(caltable, devs) %>%
    mutate(excluded = is.na(Deviation),
           slope = coef(lmcalcurve_rev)[2],
           intercept = coef(lmcalcurve_rev)[1],
           adj_R_squared = summary(lmcalcurve_rev)$adj.r.squared)
  
  # Add LLOQ and ULOQ info
  lowest_included_conc <- final_table %>%
    filter(!excluded) %>%
    summarise(LLOQ = min(Conc_orig, na.rm = TRUE)) %>%
    pull(LLOQ)
  
  highest_included_conc <- final_table %>%
    filter(!excluded) %>%
    summarise(ULOQ = max(Conc_orig, na.rm = TRUE)) %>%
    pull(ULOQ)
  
  # Calculate % of included points
  percent_included <- final_table %>%
    summarise(percent_included = mean(!excluded) * 100) %>%
    pull(percent_included)
  
  # Add to final_table
  final_table <- final_table %>%
    mutate(LLOQ = lowest_included_conc,
           ULOQ = highest_included_conc,
           percent_included = round(percent_included, 1))
  
  final_fit <- data.frame(
    Intercept = coef(lmcalcurve_rev)[1],
    Slope = coef(lmcalcurve_rev)[2],
    Adj.R.Squared = summary(lmcalcurve_rev)$adj.r.squared
  )
  
  return(final_table)
}

#---------------------------
# Now create a function to apply this for multiple compounds
OptimizeCaltableByCompound <- function(fit_function, caltable) {
  if (!"compound" %in% names(caltable)) {
    stop("The 'compound' column is required in the caltable.")
  }
  
  caltable_nested <- caltable %>%
    group_by(compound) %>%
    group_split()
  
  # Apply OptimizeCaltable2 per compound, and get only the final_table
  results <- map_dfr(caltable_nested, function(tbl) {
    compound_name <- unique(tbl$compound)
    
    # Initial fit
    fit <- lm(Area ~ Conc, data = tbl)
    
    # Apply your updated OptimizeCaltable2() which returns a final_table
    final_table <- fit_function(fit, tbl)
    
    return(final_table)  # Already includes compound column
  })
  
  return(results)
}


# Create a function to plot the results
plot_calibration_exclusions <- function(final_table) {
  # Plot with ggplot
  ggplot(final_table, aes(x = Conc, y = Area)) +
    geom_point(aes(color = excluded), size = 2.5, alpha = 0.8) +
    geom_smooth(
      data = final_table %>% filter(!excluded),
      method = "lm", se = FALSE, color = "black"
    ) +
    facet_wrap(~ compound, scales = "free") +
    scale_color_manual(values = c("FALSE" = "blue", "TRUE" = "red")) +
    labs(
      title = paste(validation_plate, "calibration curve points: Included vs. Excluded"),
      x = "Log(Concentration in ng/mL))",
      y = "Log(Area)",
      color = "Excluded"
    ) +
    theme_minimal(base_size = 14)
}


#---------------------------

# Load the data
valdata_VPA_POS <- read.table("data/Apolar_validation/p2024_051_POS_VPA.txt", header = TRUE, sep = "\t")
valdata_VPB_POS <- read.table("data/Apolar_validation/p2024_051_POS_VPB.txt", header = TRUE, sep = "\t")
valdata_VPCD_POS <- read.table("data/Apolar_validation/p2024_051_POS_VPC_VPD.txt", header = TRUE, sep = "\t")

# Define which validation plate to check, so run 1 of the 3 lines below
validation_plate <- "Validation plate A"
validation_plate <- "Validation plate B"
validation_plate <- "Validation plate CD"

if (validation_plate == "Validation plate A") {
  valdata <- valdata_VPA_POS
} else if (validation_plate == "Validation plate B") {
  valdata <- valdata_VPB_POS
} else if (validation_plate == "Validation plate CD") {
  valdata <- valdata_VPCD_POS
}

#Load the targets
targets_POS <- read.table("data/Apolar_validation/20250315_valtargets_XAB_positive.csv", header = TRUE, sep = ",")
targets <- targets_POS %>%
  select(compound, MAC_mix)

# Load the concentrations
CAL_conc <- read.table("data/Apolar_validation/CAL_level_concentrations.csv", header = TRUE, sep = ",") %>%
  mutate(MAC_mix = as.character(MAC_mix),
         cal_level = as.character(cal_level),
         Conc = concentration_ngpmL)

# Clean up the data table
CAL_data <- valdata %>%
  filter(grepl("PLCAL", Sample.Name)) %>%
  filter(!grepl("Batch", Sample.Name)) %>%
  mutate(Area = as.numeric(Area),
         compound = Component.Name) %>%
  select(Sample.Name, compound, Area) %>%
  filter(!grepl("CAL0", Sample.Name)) %>%
  separate(Sample.Name, into = c("Project", "Name", "Injection"), sep = "_", remove = FALSE) %>%
  mutate(cal_level = str_sub(Name, nchar(Name), nchar(Name))) %>%
  select(Name, Injection, compound, cal_level, Area)


# Add the MAC-mix number to the table and concentrations
CAL_data_ext <- left_join(left_join(CAL_data, targets), CAL_conc)

# Filter compounds with concentration levels
CAL_data_filtered <- CAL_data_ext %>%
  filter(!is.na(Conc)) %>%
 # filter(compound == compound[20]) %>%
  mutate(Area_orig = Area,
         Area = log10(Area),
         Conc_orig = Conc,
         Conc = log10(Conc))



# Create the results
final_results <- OptimizeCaltableByCompound(OptimizeCaltable2, CAL_data_filtered)


# Make a plot with all calcurves
plot_calibration_exclusions(final_results %>% filter(MAC_mix == 51 | MAC_mix == 52 | MAC_mix == 53 | MAC_mix == 0))

plot_calibration_exclusions(final_results %>% filter(MAC_mix == 41 | MAC_mix == 42 | MAC_mix == 36))
# Create a table with calibration curve parameters
cal_curve_table <- final_results %>%
  select(compound, slope, intercept, adj_R_squared, LLOQ, ULOQ, percent_included) %>%
  unique() %>%
  mutate(calline_accepted = percent_included > 66)


# Export the cal_curve_table
write.csv(cal_curve_table, paste("output/Apolar_validation/cal_curve_table", validation_plate, ".csv"), row.names = FALSE)

          