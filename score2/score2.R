box::use(janitor[...])

get_score2_coefs = function(sex) {
  SCORE2_men = c(
    "b_age" = 0.3742,
    "b_tc" = 0.1458,
    "b_hdlc" = -0.2698,
    "b_sbp" = 0.2777,
    "b_smoking" = 0.6012,
    "b_smo_age" = -0.0755,
    "b_sbp_age" = -0.0255,
    "b_tc_age" = -0.0281,
    "b_hdl_age" = 0.0426,
    "basesurv" = 0.9605 
  )
  SCORE2_women = c(
    "b_age" = 0.4648,
    "b_tc" = 0.1002,
    "b_hdlc" = -0.2606,
    "b_sbp" = 0.3131,
    "b_smoking" = 0.7744,
    "b_smo_age" = -0.1088,
    "b_sbp_age" = -0.0277,
    "b_tc_age" = -0.0226,
    "b_hdl_age" = 0.0613,
    "basesurv" = 0.9776
  )
  
  if (sex == 0) {
    return(SCORE2_women)
  }
  if (sex == 1) {
    return(SCORE2_men)
  }
  return(-1)
}


get_score2_scalling_factors = function(sex) {
  men = list(
    "low" = c(-0.5699, 0.7476),
    "moderate" = c(-0.1565, 0.8009),
    "high" = c(0.3207, 0.9360),
    "v_high" = c(0.5836, 0.8294)
  )
  women = list(
    "low" = c(-0.7380, 0.7019),
    "moderate" = c(-0.3143, 0.7701),
    "high" = c(0.5710, 0.9369),
    "v_high" = c(0.9412, 0.8329)
  )
  
  
  if (sex == 0) {
    return(women)
  }
  if (sex == 1) {
    return(men)
  }
  return(-1)
}

#this function is not vectorised and therefore needs to be applied using mapply
get_score2_score = function(sex, age, tc, hdlc, sbp, smoking, region) {
  
  SCORE2 = get_score2_coefs(sex)
  
  age = (as.integer(age)-60)/5
  sbp = (as.integer(sbp)-120)/20
  tc = (tc-6)/1
  hdlc = (hdlc-1.3)/0.5
  
  L=SCORE2["b_age"]*age+SCORE2["b_tc"]*tc+SCORE2["b_hdlc"]*hdlc+
    SCORE2["b_sbp"]*sbp+SCORE2["b_smoking"]*smoking+SCORE2["b_smo_age"]*smoking*age+
    SCORE2["b_sbp_age"]*sbp*age+SCORE2["b_tc_age"]*tc*age+SCORE2["b_hdl_age"]*hdlc*age
  
  A=as.numeric(exp(L))
  risk_un = 1-SCORE2["basesurv"]^A
  
  factors = get_score2_scalling_factors(sex)
  ln_risk = log(-log(1-risk_un))
  scaled_risk = 1-exp(-exp(factors[region][[1]][1]+factors[region][[1]][2]*ln_risk))
    
  return(as.numeric(100*scaled_risk))
} 

# > get_assign_score(0, 63, 0, 6.2, 2.5, 134, 0, 0, 2, 3.693)
# [1] 8.797802
# > get_assign_score(1, 56, 0, 5.2, 1.5, 137, 0, 0, 4, 2.974)
# [1] 11.53216