box::use(janitor[...])

get_assign_coefs = function(sex) {
  ASSIGN_men = c(
                  "b_age" = 0.05698,
                  "b_tc" = 0.22286,
                  "b_hdlc" = -0.53684,
                  "b_sbp" = 0.01183,
                  "b_diabetes" = 0.81558,
                  "b_family" = 0.27500,
                  "b_cpd" = 0.02005,
                  "b_SIMDSC10" = 0.06296
                )
  ASSIGN_women = c(
                  "b_age" = 0.07203,
                  "b_tc" = 0.12720,
                  "b_hdlc" = -0.55836,
                  "b_sbp" = 0.01064,
                  "b_diabetes" = 0.97727,
                  "b_family" = 0.49159,
                  "b_cpd" = 0.02724,
                  "b_SIMDSC10" = 0.09386
                 )

  if (sex == 0) {
    return(ASSIGN_women)
  }
  if (sex == 1) {
    return(ASSIGN_men)
  }
  return(-1)
}

#this function is not vectorised and therefore needs to be applied using mapply
get_assign_score = function(sex, age, ra, tc, hdlc, sbp, diabetes, family, cpd, simd10, bp_medicines) {
  
  ASSIGN = get_assign_coefs(sex)
  
  if (ra == 1) {
    cpd = cpd+10
  }
  
  if (bp_medicines == 1) {
    sbp = sbp+20
  }
  
  L=ASSIGN["b_age"]*age+ASSIGN["b_tc"]*tc+ASSIGN["b_hdlc"]*hdlc+
    ASSIGN["b_sbp"]*sbp+ASSIGN["b_diabetes"]*diabetes+ASSIGN["b_family"]*family+
    ASSIGN["b_cpd"]*cpd+ASSIGN["b_SIMDSC10"]*simd10
  
  if (sex == 0) {
    Lbar=ASSIGN["b_age"]*48.7959+ASSIGN["b_tc"]*6.40706+ASSIGN["b_hdlc"]*1.62837+
         ASSIGN["b_sbp"]*130.115+ASSIGN["b_diabetes"]*0.0127275+ASSIGN["b_family"]*0.326328+ 
         ASSIGN["b_cpd"]*6.44058+ASSIGN["b_SIMDSC10"]*2.82470
  }
  if (sex == 1) {
    Lbar=ASSIGN["b_age"]*48.8706+ASSIGN["b_tc"]*6.22520+ASSIGN["b_hdlc"]*1.35042+
         ASSIGN["b_sbp"]*133.810+ASSIGN["b_diabetes"]*0.0152905+ASSIGN["b_family"]*0.263762+
         ASSIGN["b_cpd"]*7.95841+ASSIGN["b_SIMDSC10"]*2.74038
  }
  
  A=L-Lbar
  B=as.numeric(exp(A))

  if (sex == 0) {
    P=100*(1-(0.9365^B))
  }
  if (sex == 1) {
    P=100*(1-(0.8831^B))
  }
  
  return(round_half_up(P))
} 

get_assign_score(0, 63, 0, 6.2, 2.5, 134, 0, 0, 2, 3.693, 0)
get_assign_score(1, 56, 0, 5.2, 1.5, 137, 0, 0, 4, 2.974, 0)

# > get_assign_score(0, 63, 0, 6.2, 2.5, 134, 0, 0, 2, 3.693)
# [1] 8.797802
# > get_assign_score(1, 56, 0, 5.2, 1.5, 137, 0, 0, 4, 2.974)
# [1] 11.53216