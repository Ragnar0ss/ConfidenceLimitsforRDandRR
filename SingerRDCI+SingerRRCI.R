
Singer_CI_RD <- function(
  r1_nprev,         # Sample size for risk 1
  r1_kprev,         # Frequency of positive diagnoses in sample of size r1_nprev
  r1_nsens,         # Sample size for sensitivity of risk 1
  r1_ksens,         # Frequency of positive diagnoses in sample of size r1_nsens
  r1_nspec,         # Sample size for specificity of risk 1
  r1_kspec,         # Frequency of negative diagnoses in sample of size r1_nspec
  r2_nprev,         # Sample size for risk 2
  r2_kprev,         # Frequency of positive diagnoses in sample of size r2_nprev
  r2_nsens,         # Sample size for sensitivity of risk 2
  r2_ksens,         # Frequency of positive diagnoses in sample of size r2_nsens
  r2_nspec,         # Sample size for specificity of risk 2
  r2_kspec,         # Frequency of negative diagnoses in sample of size r2_nspec
  conflevel = .95)  # Confidence level
{
  zcrit <- qnorm((1+conflevel)/2)
  
  # Improved estimates
  
  r1_imp <- (r1_kprev+zcrit^2/2) / (r1_nprev + zcrit^2)
  r2_imp <- (r2_kprev+zcrit^2/2) / (r2_nprev + zcrit^2)
  
  se1_imp <- (r1_ksens + 1) / (r1_nsens + 2)
  se2_imp <- (r2_ksens + 1) / (r2_nsens + 2)
  
  sp1_imp <- (r1_kspec + 1) / (r1_nspec + 2)
  sp2_imp <- (r2_kspec + 1) / (r2_nspec + 2)
  
  # Rogen-Gladen estimates based on the improved estimates
  
  rg_r1 <- max(min((r1_imp+sp1_imp-1)/(se1_imp+sp1_imp-1),1),0)
  rg_r2 <- max(min((r2_imp+sp2_imp-1)/(se2_imp+sp2_imp-1),1),0)
  
  # Calculation of the confidence intervals of the risks
  
  varp_r1 <- ((r1_imp*(1-r1_imp)/(r1_nprev+zcrit^2))+(rg_r1^2*se1_imp*(1-se1_imp)/(r1_nsens+2))+((1-rg_r1)^2*sp1_imp*(1-sp1_imp)/(r1_nspec+2)))/(se1_imp+sp1_imp-1)^2
  varp_r2 <- ((r2_imp*(1-r2_imp)/(r2_nprev+zcrit^2))+(rg_r2^2*se2_imp*(1-se2_imp)/(r2_nsens+2))+((1-rg_r2)^2*sp2_imp*(1-sp2_imp)/(r2_nspec+2)))/(se2_imp+sp2_imp-1)^2
  
  dp_1 <- 2*zcrit^2*((rg_r1 * se1_imp * (1 - se1_imp) / (r1_nsens + 2)) - ((1 - rg_r1) * sp1_imp * (1 - sp1_imp) / (r1_nspec + 2)))
  dp_2 <- 2*zcrit^2*((rg_r2 * se2_imp * (1 - se2_imp) / (r2_nsens + 2)) - ((1 - rg_r2) * sp2_imp * (1 - sp2_imp) / (r2_nspec + 2)))
  
  low1 <- rg_r1+dp_1-zcrit*sqrt(varp_r1)
  upp1 <- rg_r1+dp_1+zcrit*sqrt(varp_r1)
  
  low2 <- rg_r2+dp_2-zcrit*sqrt(varp_r2)
  upp2 <- rg_r2+dp_2+zcrit*sqrt(varp_r2)
  
  # Calculation of the confidence intervals of the RD
  
  CI_lower <- rg_r1-rg_r2-sqrt((rg_r1-low1)^2+(upp2-rg_r2)^2)
  CI_upper <- rg_r1-rg_r2+sqrt((rg_r2-low2)^2+(upp1-rg_r1)^2)
  
  ci = c(CI_lower, CI_upper)
  return(ci)
}

Singer_CI_RD(r1_nprev=509, r1_kprev=3,
             r2_nprev=513, r2_kprev=25,
             r1_nsens=168, r1_ksens=93,
             r2_nsens=168, r2_ksens=93,
             r1_nspec=169, r1_kspec=162,
             r2_nspec=169, r2_kspec=162)

Singer_CI_RD_adj <- function(
    r1_nprev,         # Sample size for risk 1
    r1_kprev,         # Frequency of positive diagnoses in sample of size r1_nprev
    r1_nsens,         # Sample size for sensitivity of risk 1
    r1_ksens,         # Frequency of positive diagnoses in sample of size r1_nsens
    r1_nspec,         # Sample size for specificity of risk 1
    r1_kspec,         # Frequency of negative diagnoses in sample of size r1_nspec
    r2_nprev,         # Sample size for risk 2
    r2_kprev,         # Frequency of positive diagnoses in sample of size r2_nprev
    r2_nsens,         # Sample size for sensitivity of risk 2
    r2_ksens,         # Frequency of positive diagnoses in sample of size r2_nsens
    r2_nspec,         # Sample size for specificity of risk 2
    r2_kspec,         # Frequency of negative diagnoses in sample of size r2_nspec
    conflevel = .95)  # Confidence level
{
  zcrit <- qnorm((1+conflevel)/2)
  
  # Improved estimates
  
  r1_imp <- (r1_kprev+zcrit^2/2) / (r1_nprev + zcrit^2)
  r2_imp <- (r2_kprev+zcrit^2/2) / (r2_nprev + zcrit^2)
  
  se1_imp <- (r1_ksens + 1) / (r1_nsens + 2)
  se2_imp <- (r2_ksens + 1) / (r2_nsens + 2)
  
  sp1_imp <- (r1_kspec + 1) / (r1_nspec + 2)
  sp2_imp <- (r2_kspec + 1) / (r2_nspec + 2)
  
  # Rogen-Gladen estimates based on the improved estimates
  
  rg_r1 <- max(min((r1_imp+sp1_imp-1)/(se1_imp+sp1_imp-1),1),0)
  rg_r2 <- max(min((r2_imp+sp2_imp-1)/(se2_imp+sp2_imp-1),1),0)
  
  # Calculation of the confidence intervals of the risks
  
  varp_r1 <- ((r1_imp*(1-r1_imp)/(r1_nprev+zcrit^2))+(rg_r1^2*se1_imp*(1-se1_imp)/(r1_nsens+2))+((1-rg_r1)^2*sp1_imp*(1-sp1_imp)/(r1_nspec+2)))/(se1_imp+sp1_imp-1)^2
  varp_r2 <- ((r2_imp*(1-r2_imp)/(r2_nprev+zcrit^2))+(rg_r2^2*se2_imp*(1-se2_imp)/(r2_nsens+2))+((1-rg_r2)^2*sp2_imp*(1-sp2_imp)/(r2_nspec+2)))/(se2_imp+sp2_imp-1)^2
  #print(varp_r1)
  #print(varp_r2)
  
  dp_1 <- 2*zcrit^2*((rg_r1 * se1_imp * (1 - se1_imp) / (r1_nsens + 2)) - ((1 - rg_r1) * sp1_imp * (1 - sp1_imp) / (r1_nspec + 2)))
  dp_2 <- 2*zcrit^2*((rg_r2 * se2_imp * (1 - se2_imp) / (r2_nsens + 2)) - ((1 - rg_r2) * sp2_imp * (1 - sp2_imp) / (r2_nspec + 2)))
  #print(dp_1)
  #print(dp_2)
  
  low1 <- rg_r1+dp_1-zcrit*sqrt(varp_r1)
  upp1 <- rg_r1+dp_1+zcrit*sqrt(varp_r1)
  #print(low1)
  #print(upp1)
  
  low2 <- rg_r2+dp_2-zcrit*sqrt(varp_r2)
  upp2 <- rg_r2+dp_2+zcrit*sqrt(varp_r2)
  #print(low2)
  #print(upp2)
  
  # Adjustments based on Zou and Donner
  
  pi1 <- min((1-1/2/r1_nprev), max(rg_r1, 1/2/r1_nprev))
  pi2 <- min((1-1/2/r2_nprev), max(rg_r2, 1/2/r2_nprev))
  #print(pi1)
  #print(pi2)
  
  low1_adj = min(max(low1,1/2/r1_nprev),1-1/2/r1_nprev)
  upp1_adj = min(max(upp1,1/2/r1_nprev),1-1/2/r1_nprev)
  low2_adj = min(max(low2,1/2/r2_nprev),1-1/2/r2_nprev)
  upp2_adj = min(max(upp2,1/2/r2_nprev),1-1/2/r2_nprev)
  #print(low1_adj)
  #print(upp1_adj)
  #print(low2_adj)
  #print(upp2_adj)
  
  # Calculation of the confidence intervals of the RD
  
  CI_lower <- pi1-pi2-sqrt((pi1-low1_adj)^2+(upp2_adj-pi2)^2)
  CI_upper <- pi1-pi2+sqrt((pi2-low2_adj)^2+(upp1_adj-pi1)^2)
  
  ci = c(CI_lower, CI_upper)
  return(ci)
}


Singer_CI_RD_adj(r1_nprev=100, r1_kprev=28,
             r2_nprev=100, r2_kprev=17,
             r1_nsens=24, r1_ksens=22,
             r2_nsens=24, r2_ksens=22,
             r1_nspec=300, r1_kspec=297,
             r2_nspec=300, r2_kspec=297)

Singer_CI_RR <- function(
  r1_nprev,         # Sample size for risk 1
  r1_kprev,         # Frequency of positive diagnoses in sample of size r1_nprev
  r1_nsens,         # Sample size for sensitivity of risk 1
  r1_ksens,         # Frequency of positive diagnoses in sample of size r1_nsens
  r1_nspec,         # Sample size for specificity of risk 1
  r1_kspec,         # Frequency of negative diagnoses in sample of size r1_nspec
  r2_nprev,         # Sample size for risk 2
  r2_kprev,         # Frequency of positive diagnoses in sample of size r2_nprev
  r2_nsens,         # Sample size for sensitivity of risk 2
  r2_ksens,         # Frequency of positive diagnoses in sample of size r2_nsens
  r2_nspec,         # Sample size for specificity of risk 2
  r2_kspec,         # Frequency of negative diagnoses in sample of size r2_nspec
  conflevel = .95)  # Confidence level
{
  
  zcrit <- qnorm((1+conflevel)/2)
  # Improved estimates
  
  r1_imp <- (r1_kprev+zcrit^2/2) / (r1_nprev + zcrit^2)
  r2_imp <- (r2_kprev+zcrit^2/2) / (r2_nprev + zcrit^2)
  #print(r1_imp)
  #print(r2_imp)
  
  se1_imp <- (r1_ksens + 1) / (r1_nsens + 2)
  se2_imp <- (r2_ksens + 1) / (r2_nsens + 2)
  #print(se1_imp)
  #print(se2_imp)
  
  sp1_imp <- (r1_kspec + 1) / (r1_nspec + 2)
  sp2_imp <- (r2_kspec + 1) / (r2_nspec + 2)
  #print(sp1_imp)
  #print(sp2_imp)
  
  # Rogen-Gladen estimates based on the improved estimates
  
  rg_r1 <- max(min((r1_imp+sp1_imp-1)/(se1_imp+sp1_imp-1),1),0)
  rg_r2 <- max(min((r2_imp+sp2_imp-1)/(se2_imp+sp2_imp-1),1),0)
  #print(rg_r1)
  #print(rg_r2)
  
  # Calculation of the confidence intervals of the risks
  
  varp_r1 <- ((r1_imp*(1-r1_imp)/(r1_nprev+zcrit^2))+(rg_r1^2*se1_imp*(1-se1_imp)/(r1_nsens+2))+((1-rg_r1)^2*sp1_imp*(1-sp1_imp)/(r1_nspec+2)))/(se1_imp+sp1_imp-1)^2
  varp_r2 <- ((r2_imp*(1-r2_imp)/(r2_nprev+zcrit^2))+(rg_r2^2*se2_imp*(1-se2_imp)/(r2_nsens+2))+((1-rg_r2)^2*sp2_imp*(1-sp2_imp)/(r2_nspec+2)))/(se2_imp+sp2_imp-1)^2
  #print(varp_r1)
  #print(varp_r2)
  
  dp_1 <- 2*zcrit^2*((rg_r1 * se1_imp * (1 - se1_imp) / (r1_nsens + 2)) - ((1 - rg_r1) * sp1_imp * (1 - sp1_imp) / (r1_nspec + 2)))
  dp_2 <- 2*zcrit^2*((rg_r2 * se2_imp * (1 - se2_imp) / (r2_nsens + 2)) - ((1 - rg_r2) * sp2_imp * (1 - sp2_imp) / (r2_nspec + 2)))
  #print(dp_1)
  #print(dp_2)
  
  low1 <- rg_r1+dp_1-zcrit*sqrt(varp_r1)
  upp1 <- rg_r1+dp_1+zcrit*sqrt(varp_r1)
  
  low2 <- rg_r2+dp_2-zcrit*sqrt(varp_r2)
  upp2 <- rg_r2+dp_2+zcrit*sqrt(varp_r2)

  # Calculation of the confidence intervals of the RR
  
  pi1 <- min((1-1/2/r1_nprev), max(rg_r1, 1/2/r1_nprev))
  pi2 <- min((1-1/2/r2_nprev), max(rg_r2, 1/2/r2_nprev))
  
  CI_lower <- exp(log(pi1)-log(pi2)-sqrt((log(pi1)-log(low1))^2+(log(upp2)-log(pi2))^2))
  CI_upper <- exp(log(pi1)-log(pi2)+sqrt((log(pi2)-log(low2))^2+(log(upp1)-log(pi1))^2))
  
  ci = c(CI_lower, CI_upper)
  return(ci)
}

Singer_CI_RR(r1_nprev=67, r1_kprev=18,
             r2_nprev=252, r2_kprev=90,
             r1_nsens=33, r1_ksens=32,
             r2_nsens=33, r2_ksens=32,
             r1_nspec=20, r1_kspec=20,
             r2_nspec=20, r2_kspec=20)




