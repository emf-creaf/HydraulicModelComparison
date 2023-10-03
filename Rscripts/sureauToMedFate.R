library(medfate)
# Implicit time integration function on small time step dt

# Notes : np1 means "n+1", nph means "n+1/2" (for n plus half)
PLC.comp <- function(Pmin, slope , P50) {
  PLC  = 100 / (1 + exp(slope / 25 * (Pmin - P50)))
  return(PLC)
}
PLCPrime.comp <- function(PLC , slope) {
  return(- slope/25 * PLC/100 * (1 - PLC/100))
}
implicit.temporal.integration.atnp1 <- function(WBveg, WBsoil, dt, opt) {
  
  # 1. Initializing current time step according to computation options (FP)
  
  dbxmin = 1e-100 # FP minimal double to avoid 0/0
  
  Psi_LApo_n = WBveg$Psi_LApo
  
  Psi_SApo_n = WBveg$Psi_SApo
  
  Psi_LSym_n = WBveg$Psi_LSym
  
  Psi_SSym_n = WBveg$Psi_SSym
  
  Psi_LApo_cav = WBveg$Psi_LApo_cav
  
  Psi_SApo_cav = WBveg$Psi_SApo_cav
  
  
  
  K_LSym = opt$Lsym * WBveg$k_LSym   #
  
  K_SSym = opt$Ssym * WBveg$k_SSym   #
  
  C_LSym = opt$Lsym * WBveg$C_LSym   #
  
  C_SSym = opt$Ssym * WBveg$C_SSym   #
  
  
  
  C_LApo = opt$CLapo * WBveg$C_LApo #
  
  C_SApo = opt$CTapo * WBveg$C_SApo #
  
  
  
  
  
  K_SL = WBveg$k_SLApo #
  
  E_nph = WBveg$Elim # COMPUTED WITH WBveg$Psi_LSym AND INTERPOLATE CLIMATE
  
  Eprime_nph = opt$Eord * WBveg$Eprime
  
  
  
  Emin_L_nph = WBveg$Emin
  
  Emin_S_nph = WBveg$Emin_S
  
  
  
  
  
  
  
  #Compute K_L_Cav et K_S_Cav
  
  
  
  PLC_prime_L = PLCPrime.comp(WBveg$PLC_Leaf,WBveg$params$slope_VC_Leaf)
  
  #K_L_Cav = -opt$Lcav * WBveg$Q_LApo_sat_mmol * PLC_prime_L / dt  # avec WBveg$Q_LSym_sat en l/m2 sol
  
  K_L_Cav = -opt$Lcav * WBveg$Q_LApo_sat_mmol_perLeafArea * PLC_prime_L / dt  # avec WBveg$Q_LSym_sat en l/m2 sol # changed by NM (25/10/2021)
  
  PLC_prime_S = PLCPrime.comp(WBveg$PLC_Stem,WBveg$params$slope_VC_Stem)
  
  #K_S_Cav = -opt$Scav * WBveg$Q_SApo_sat_mmol * PLC_prime_S / dt  #opt$Scav * WBveg$K_S_Cav #FP corrected a bug sign here
  
  K_S_Cav = -opt$Scav * WBveg$Q_SApo_sat_mmol_perLeafArea * PLC_prime_S / dt  #opt$Scav * WBveg$K_S_Cav #FP corrected a bug sign herehanged by NM (25/10/2021)
  
  # print(c(PLC_prime_L, K_L_Cav, PLC_prime_S, K_S_Cav))
  
  
  
  # 2. While loop in order to decide if cavitation or not :
  
  # In order to account for the cavitation that occurs only when potentials go below their lowest value "cav" (formerly called "mem" in an earlier version)
  
  # the following computations are done trying sequentially the resolutions of LApo and TApo eventually activating
  
  # the appropriate cavitation events when needed (starting assuming no cavit at all...)
  
  # in case of computational problem, the last case assume no cavitation flux
  
  LcavitWellComputed=F #initialized to false
  
  ScavitWellComputed=F
  
  if ((opt$Lcav==0)&(opt$Scav==0)) { # no cavitation flux computed
    
    delta_L_cavs=c(0)
    
    delta_S_cavs=c(0)
    
  } else if ((opt$Lcav==0)&(opt$Scav==1)) {# Scav only
    
    delta_L_cavs=c(0,0,0)
    
    delta_S_cavs=c(0,1,0)
    
  } else if ((opt$Lcav==0)&(opt$Scav==1)) {# Lcav only
    
    delta_L_cavs=c(0,1,0)
    
    delta_S_cavs=c(0,0,0)
    
  } else { #Lcav=1 and Scav=1 # Mode normal all options are open
    
    delta_L_cavs=c(0,1,0,1,0) # the fifth case is here in case no solution with others...
    
    delta_S_cavs=c(0,0,1,1,0)
    
  }
  
  if (opt$numericalScheme=="Explicit"){ # in case of explicit we do the computation only once
    
    delta_L_cavs=c(opt$Lcav)
    
    delta_S_cavs=c(opt$Scav)
    
  }
  
  nwhilecomp = 0 # count the number of step in while loop (if more than 4 no solution and warning)
  
  while (((!LcavitWellComputed)|(!ScavitWellComputed))&(nwhilecomp<length(delta_L_cavs))) {
    
    #browser()
    
    nwhilecomp = nwhilecomp + 1
    
    #print(paste0('nwhilecomp=',nwhilecomp)) for debugging
    
    delta_L_cav = delta_L_cavs[nwhilecomp]
    
    delta_S_cav = delta_S_cavs[nwhilecomp]
    
    
    
    if (opt$numericalScheme=="Implicit") {
      
      # 2.1 Compute intermediates
      
      
      
      # compute Psi_L_tilda and K_L_tilda and E_L_tilda
      
      klsym = kseriesum(K_LSym, C_LSym/dt+0.5 * Eprime_nph) # for Psi_LSym_n
      
      K_L_td =  C_LApo/dt + klsym + delta_L_cav*K_L_Cav
      
      Psi_L_td = (C_LApo/dt*Psi_LApo_n + klsym*Psi_LSym_n + delta_L_cav*K_L_Cav*Psi_LApo_cav)/(K_L_td + dbxmin) # dbxmin to avoid 0/0
      
      E_L_tilda = (E_nph + Emin_L_nph)/(1+(C_LSym/dt+0.5 * Eprime_nph+dbxmin)/K_LSym)
      
      
      
      # compute Psi_S_tilda and K_S_tilda and E_S_tilda
      
      kssym = kseriesum(K_SSym, C_SSym/dt) # for Psi_SSym_n
      
      K_S_td = C_SApo/dt + kssym + sum(WBveg$k_SoilToStem)  + delta_S_cav*K_S_Cav
      
      Psi_S_Td = (C_SApo/dt*Psi_SApo_n + kssym*Psi_SSym_n + sum(WBveg$k_SoilToStem * WBsoil$PsiSoil) + delta_S_cav*K_S_Cav*Psi_SApo_cav) / (K_S_td + dbxmin) # dbxmin to avoid 0/0
      
      E_S_tilda = K_SL/(K_SL + K_S_td) * Emin_S_nph/(1+(C_SSym/dt + dbxmin)/K_SSym)  # dbxmin to avoid 0/0 
      
      
      
      # 2.2 Compute Psi_LApo_np1
      
      Psi_LApo_np1_Num = kseriesum(K_SL , K_S_td + dbxmin)*Psi_S_Td + K_L_td*Psi_L_td - (E_L_tilda + E_S_tilda)
      
      Psi_LApo_np1_Denom = kseriesum(K_SL, K_S_td + dbxmin) + K_L_td + dbxmin # dbxmin to avoid 0/0
      
      Psi_LApo_np1 = Psi_LApo_np1_Num/Psi_LApo_np1_Denom
      
      
      
      # 2.3 Compute Psi_SApo_np1
      
      Psi_SApo_np1 = ((K_L_td + K_SL)*Psi_LApo_np1 - K_L_td*Psi_L_td + E_L_tilda)/(K_SL+ dbxmin)
      
    } else if (opt$numericalScheme=="Semi-Implicit") {
      
      # 2.1 LApo
      
      alpha = exp(-(K_SL+K_LSym+delta_L_cav*K_L_Cav)/C_LApo*dt)
      
      Psi_td = (K_SL*Psi_SApo_n + K_LSym*Psi_LSym_n + delta_L_cav*K_L_Cav*Psi_LApo_cav)/(K_SL + K_LSym+delta_L_cav*K_L_Cav + dbxmin) # dbxmin to avoid 0/0
      
      Psi_LApo_np1 = alpha * Psi_LApo_n +(1-alpha) * Psi_td
      
      
      
      # 2.2. SApo
      
      alpha = exp(-(K_SL+K_SSym + sum(WBveg$k_SoilToStem)+delta_S_cav*K_S_Cav)/C_SApo*dt)
      
      Psi_td = (K_SL*Psi_LApo_n + K_SSym*Psi_SSym_n + sum(WBveg$k_SoilToStem * WBsoil$PsiSoil)+ delta_S_cav*K_S_Cav*Psi_SApo_cav)/(K_SL + K_SSym+sum(WBveg$k_SoilToStem)+delta_S_cav*K_S_Cav + dbxmin) # dbxmin to avoid 0/0
      
      # print(c(sum(WBveg$k_SoilToStem * WBsoil$PsiSoil), sum(WBveg$k_SoilToStem), alpha, Psi_td))
      Psi_SApo_np1 = alpha * Psi_SApo_n +(1-alpha) * Psi_td
      
      # print(c(Psi_LApo_np1, Psi_LApo_cav, Psi_SApo_np1, Psi_SApo_cav))
      
      
    } else if (opt$numericalScheme=="Explicit"){
      
      # NB for cavitation we use the min(cav-current) because here cav is at time n-1
      
      # psi L_apo
      
      # if (Psi_LApo_cav > Psi_LApo_n) {
      
      #   print(paste0("Lcav",opt$Lcav))
      
      # }
      
      # if (Psi_SApo_cav > Psi_SApo_n) {
      
      #   print(paste0("Scav",opt$Scav))
      
      # }
      
      Psi_LApo_np1 = Psi_LApo_n + (dt/C_LApo) * (K_SL * (Psi_SApo_n - Psi_LApo_n) + K_LSym * (Psi_LSym_n - Psi_LApo_n) + delta_L_cav*K_L_Cav * max(Psi_LApo_cav - Psi_LApo_n,0))
      
      # psi T_apo
      
      Psi_SApo_np1 = Psi_SApo_n + (dt/C_SApo) * (K_SL * (Psi_LApo_n - Psi_SApo_n) + K_SSym * (Psi_SSym_n - Psi_SApo_n) + delta_S_cav*K_S_Cav * max(Psi_SApo_cav - Psi_SApo_n,0) + sum(WBveg$k_SoilToStem * (WBsoil$PsiSoil - Psi_SApo_n)))
      
      
      
      # determine the cfl for each cell
      
      cfl_LApo = C_LApo/(2*max(K_SL,K_LSym,K_L_Cav))
      
      #cfl_LApo = C_LApo/(2*max(K_SL,K_LSym))
      
      cfl_SApo = C_SApo/(2*max(K_SL,K_SSym,K_S_Cav,max(WBveg$k_SoilToStem)))
      
      
      
      
      
    }
    
    
    
    # 2.4 check if cavitation is well computed according to delta_cav, np1 and "cav"
    
    LcavitWellComputed = (delta_L_cav==(Psi_LApo_np1 < Psi_LApo_cav))|(opt$Lcav==0)
    
    ScavitWellComputed = (delta_S_cav==(Psi_SApo_np1 < Psi_SApo_cav))|(opt$Scav==0)
    
    if ((length(delta_L_cavs)>1)&(nwhilecomp==length(delta_L_cavs))) { # we tried the normal cases and the computation is still not ok so we have done a last one desactivating cavitation water source (delta_cav=0)
      
      warning(paste0("water flux due to Cavitation ignored with time step, no solution from the implicit solver=",dt))
      
    }
    
  } # end of the while loop with check on cavitation options
  
  WBveg$Diag_nwhile_cavit = nwhilecomp  # Diagnostic step to track cavit event and eventual errors (corresponding to nwhilecomp==5)
  
  
  
  # 3. Compute Psi_Symp_np1 (L and S)
  
  if (opt$numericalScheme=="Implicit") {
    
    klsym = C_LSym/dt+0.5 * Eprime_nph # for Psi_LSym_n
    
    Psi_LSym_np1 = (K_LSym*Psi_LApo_np1 + klsym*Psi_LSym_n - (E_nph + Emin_L_nph)) / (K_LSym + klsym + dbxmin) # dbxmin to avoid 0/0
    
    Psi_SSym_np1 = (K_SSym*Psi_SApo_np1 + C_SSym/dt*Psi_SSym_n - Emin_S_nph) / (K_SSym + C_SSym/dt + dbxmin) # dbxmin to avoid 0/0
    
  } else if (opt$numericalScheme=="Semi-Implicit") {
    
    alpha = exp(-K_LSym/C_LSym*dt)
    
    Psi_td = (K_LSym*Psi_LApo_n - (E_nph + Emin_L_nph))/(K_LSym + dbxmin) # dbxmin to avoid 0/0
    
    Psi_LSym_np1 = alpha * Psi_LSym_n +(1-alpha) * Psi_td
    
    alpha = exp(-K_SSym/C_SSym*dt)
    
    Psi_td = (K_SSym*Psi_SApo_n - Emin_S_nph)/(K_SSym + dbxmin) # dbxmin to avoid 0/0
    
    Psi_SSym_np1 = alpha * Psi_SSym_n +(1-alpha) * Psi_td
    
  } else if (opt$numericalScheme=="Explicit") {
    
    # psi L_Sym_np1 and cfl
    
    
    
    if(C_LSym==0){ Psi_LSym_np1 = ((K_LSym * Psi_LApo_n) - E_nph - Emin_L_nph)/(K_LSym ) # cas stationnaire lorsque C_Lsymp=0
    
    cfl_LSym = 10
    
    }else {
      
      Psi_LSym_np1 =  Psi_LSym_n + (dt/(C_LSym)) * (K_LSym * (Psi_LApo_n - Psi_LSym_n) - E_nph - Emin_L_nph)
      
      cfl_LSym = C_LSym/(2*K_LSym)}
    
    
    
    # psi T_Sym_np1 and cfls
    
    if (C_SSym==0){Psi_SSym_np1 = ((K_SSym * Psi_SApo_n) - Emin_S_nph)/(K_SSym )# cas stationnaire lorsque C_Lsymp=0
    
    cfl_SSym = 10 # cfl non limitante
    
    }else{Psi_SSym_np1 =  Psi_SSym_n + (dt/(C_SSym)) * (K_SSym * (Psi_SApo_n - Psi_SSym_n) - Emin_S_nph)
    
    cfl_SSym = C_SSym/(2*K_SSym)}
    
    
    
    
    
    # print(paste0("cfl_LSym = ",cfl_LSym))
    
    # print(paste0("cfl_SSym = ",cfl_SSym))
    
    # print(paste0("cfl_LApo = ",cfl_LApo))
    
    # print(paste0("cfl_SApo = ",cfl_SApo))
    
    
    
    cfl_all = min(cfl_LSym,cfl_SSym,cfl_LApo,cfl_SApo)
    
    #print(paste0("cfl_all = ",cfl_all))
    
    # if (cfl_all<=dt)
    
    #   {
    
    #   stop('execution stopped because cfl<dt')
    
    #   }
    
  }
  
  
  
  
  
  #Step 4 : set computed values in WBveg and update Psi_cav, PLC and Psi_AllSoil
  
  WBveg$Psi_LApo<-min(-0.00001,Psi_LApo_np1)
  
  WBveg$Psi_SApo<-min(-0.00001,Psi_SApo_np1)
  
  WBveg$Psi_LSym<-min(-0.00001,Psi_LSym_np1)
  
  WBveg$Psi_SSym<-min(-0.00001,Psi_SSym_np1)
  
  
  
  # Cavitation
  
  if (opt$numericalScheme=="Explicit"){
    
    psiref = min(-0.00001,Psi_LApo_n)  # the reference is at previous time step
    
  } else {
    
    psiref = WBveg$Psi_LApo  # the reference is at current time step for other modes
    
  }
  
  if (psiref < Psi_LApo_cav) {
    
    WBveg$Psi_LApo_cav = psiref
    
    WBveg$PLC_Leaf = PLC.comp(Pmin = psiref, slope = WBveg$params$slope_VC_Leaf, P50 = WBveg$params$P50_VC_Leaf)
    
  }
  
  if (opt$numericalScheme=="Explicit"){
    
    psiref = min(-0.00001,Psi_SApo_n)  # the reference is at previous time step
    
  } else {
    
    psiref = WBveg$Psi_SApo  # the reference is at current time step for other modes
    
  }
  
  if (psiref < Psi_SApo_cav) {
    
    WBveg$Psi_SApo_cav = psiref
    
    WBveg$PLC_Stem = PLC.comp(Pmin = psiref, slope = WBveg$params$slope_VC_Stem, P50 = WBveg$params$P50_VC_Stem)
    
  }
  
  
  
  # browser()
  
  WBveg$Psi_AllSoil <- sum (WBveg$k_SoilToStem * WBsoil$PsiSoil)/sum (WBveg$k_SoilToStem)
  
  return(WBveg)
  
}
semi_implicit_may <-function(WBveg, WBsoil, dt, opt) {
  # Step 1. Initializing current time step according to computation options (FP)
  dbxmin = 1e-100 # FP minimal double to avoid 0/0
  Psi_LApo_n = WBveg$Psi_LApo
  Psi_SApo_n = WBveg$Psi_SApo
  Psi_LSym_n = WBveg$Psi_LSym
  Psi_SSym_n = WBveg$Psi_SSym
  Psi_LApo_cav = WBveg$Psi_LApo_cav
  Psi_SApo_cav = WBveg$Psi_SApo_cav
  
  K_LSym = opt$Lsym * WBveg$k_LSym   #
  K_SSym = opt$Ssym * WBveg$k_SSym   #
  C_LSym = opt$Lsym * WBveg$C_LSym   #
  C_SSym = opt$Ssym * WBveg$C_SSym   #
  
  C_LApo = opt$CLapo * WBveg$C_LApo #
  C_SApo = opt$CTapo * WBveg$C_SApo #
  
  
  K_SL = WBveg$k_SLApo #
  E_nph = WBveg$Elim # COMPUTED WITH WBveg$Psi_LSym AND INTERPOLATE CLIMATE
  Eprime_nph = opt$Eord * WBveg$Eprime
  
  Emin_L_nph = WBveg$Emin
  Emin_S_nph = WBveg$Emin_S
  
  
  
  #Compute K_L_Cav et K_S_Cav
  PLC_prime_L = PLCPrime.comp(WBveg$PLC_Leaf,WBveg$params$slope_VC_Leaf)
  K_L_Cav = -opt$Lcav * WBveg$Q_LApo_sat_mmol_perLeafArea * PLC_prime_L / dt  # avec WBveg$Q_LSym_sat en l/m2 sol # changed by NM (25/10/2021)
  PLC_prime_S = PLCPrime.comp(WBveg$PLC_Stem,WBveg$params$slope_VC_Stem)
  K_S_Cav = -opt$Scav * WBveg$Q_SApo_sat_mmol_perLeafArea * PLC_prime_S / dt  #opt$Scav * WBveg$K_S_Cav #FP corrected a bug sign herehanged by NM (25/10/2021)
  
  
  # Step 2. While loop in order to decide if cavitation or not :
  # In order to account for the cavitation that occurs only when potentials go below their lowest value "cav" (formerly called "mem" in an earlier version)
  # the following computations are done trying sequentially the resolutions of LApo and TApo eventually activating
  # the appropriate cavitation events when needed (starting assuming no cavit at all...)
  # in case of computational problem, the last case assume no cavitation flux
  LcavitWellComputed=F #initialized to false
  ScavitWellComputed=F
  if ((opt$Lcav==0)&(opt$Scav==0)) { # no cavitation flux computed
    delta_L_cavs=c(0)
    delta_S_cavs=c(0)
  } else if ((opt$Lcav==0)&(opt$Scav==1)) {# Scav only
    delta_L_cavs=c(0,0,0)
    delta_S_cavs=c(0,1,0)
  } else if ((opt$Lcav==0)&(opt$Scav==1)) {# Lcav only
    delta_L_cavs=c(0,1,0)
    delta_S_cavs=c(0,0,0)
  } else { #Lcav=1 and Scav=1
    delta_L_cavs=c(0,1,0,1,0) # the fifth case is here in case no solution with others...
    delta_S_cavs=c(0,0,1,1,0)
  }
  
  
  nwhilecomp = 0 # count the number of step in while loop (if more than 4 no solution and warning)
  while (((!LcavitWellComputed)|(!ScavitWellComputed))&(nwhilecomp<length(delta_L_cavs)))
  {
    
    nwhilecomp = nwhilecomp + 1
    delta_L_cav = delta_L_cavs[nwhilecomp]
    delta_S_cav = delta_S_cavs[nwhilecomp]
    
    # 2.1 LApo
    alpha = exp(-(K_SL+K_LSym+delta_L_cav*K_L_Cav)/C_LApo*dt)
    Psi_td = (K_SL*Psi_SApo_n + K_LSym*Psi_LSym_n + delta_L_cav*K_L_Cav*Psi_LApo_cav)/(K_SL + K_LSym+delta_L_cav*K_L_Cav + dbxmin) # dbxmin to avoid 0/0
    Psi_LApo_np1 = alpha * Psi_LApo_n + (1-alpha) * Psi_td
    
    # 2.2. SApo
    alpha = exp(-(K_SL+K_SSym + sum(WBveg$k_SoilToStem)+delta_S_cav*K_S_Cav)/C_SApo*dt)
    Psi_td = (K_SL*Psi_LApo_n + K_SSym*Psi_SSym_n + sum(WBveg$k_SoilToStem * WBsoil$PsiSoil)+ delta_S_cav*K_S_Cav*Psi_SApo_cav)/(K_SL + K_SSym+sum(WBveg$k_SoilToStem)+delta_S_cav*K_S_Cav + dbxmin) # dbxmin to avoid 0/0
    Psi_SApo_np1 = alpha * Psi_SApo_n + (1-alpha) * Psi_td
    
    # 2.3 Compute Psi_SApo_np1
    Psi_SApo_np1 = ((K_L_td + K_SL)*Psi_LApo_np1 - K_L_td*Psi_L_td + E_L_tilda)/(K_SL+ dbxmin)
    
    # 2.4 check if cavitation is well computed according to delta_cav, np1 and "cav"
    LcavitWellComputed = (delta_L_cav==(Psi_LApo_np1 < Psi_LApo_cav))|(opt$Lcav==0)
    ScavitWellComputed = (delta_S_cav==(Psi_SApo_np1 < Psi_SApo_cav))|(opt$Scav==0)
    if ((length(delta_L_cavs)>1)&(nwhilecomp==length(delta_L_cavs))) { # we tried the normal cases and the computation is still not ok so we have done a last one desactivating cavitation water source (delta_cav=0)
      warning(paste0("water flux due to Cavitation ignored with time step, no solution from the implicit solver=",dt))
    }
  } # end of the while loop with check on cavitation options
  WBveg$Diag_nwhile_cavit = nwhilecomp  # Diagnostic step to track cavit event and eventual errors (corresponding to nwhilecomp==5)
  
  # Step 3. Compute Psi_Symp_np1 (L and S)
  alpha = exp(-K_LSym/C_LSym*dt)
  Psi_td = (K_LSym*Psi_LApo_n - (E_nph + Emin_L_nph))/(K_LSym + dbxmin) # dbxmin to avoid 0/0
  Psi_LSym_np1 = alpha * Psi_LSym_n +(1-alpha) * Psi_td
  alpha = exp(-K_SSym/C_SSym*dt)
  Psi_td = (K_SSym*Psi_SApo_n - Emin_S_nph)/(K_SSym + dbxmin) # dbxmin to avoid 0/0
  Psi_SSym_np1 = alpha * Psi_SSym_n +(1-alpha) * Psi_td
  
  #Step 4 : set computed values in WBveg and update Psi_cav, PLC and Psi_AllSoil
  WBveg$Psi_LApo<-min(-0.00001,Psi_LApo_np1)
  WBveg$Psi_SApo<-min(-0.00001,Psi_SApo_np1)
  WBveg$Psi_LSym<-min(-0.00001,Psi_LSym_np1)
  WBveg$Psi_SSym<-min(-0.00001,Psi_SSym_np1)
  
  # Cavitation
  psiref = WBveg$Psi_LApo  # the reference is at current time step for other modes  (implicit, explicit)
  
  if (psiref < Psi_LApo_cav) {
    WBveg$Psi_LApo_cav = psiref
    WBveg$PLC_Leaf = PLC.comp(Pmin = psiref, slope = WBveg$params$slope_VC_Leaf, P50 = WBveg$params$P50_VC_Leaf)
  }
  
  psiref = WBveg$Psi_SApo  # The reference is at current time step for other modes (implicit, explicit)
  
  if (psiref < Psi_SApo_cav) {
    WBveg$Psi_SApo_cav = psiref
    WBveg$PLC_Stem = PLC.comp(Pmin = psiref, slope = WBveg$params$slope_VC_Stem, P50 = WBveg$params$P50_VC_Stem)
  }
  
  
  WBveg$Psi_AllSoil <- sum (WBveg$k_SoilToStem * WBsoil$PsiSoil)/sum (WBveg$k_SoilToStem)
  return(WBveg)
}
transform_sureauToNetwork <- function(network, sureau_Veg, sureau_Soil) {
  # network$params$TPhase_gmin = sureau_Veg$params$TPhase_gmin
  # network$params$Q10_1_gmin = sureau_Veg$params$Q10_1_gmin
  # network$params$Q10_2_gmin = sureau_Veg$params$Q10_2_gmin
  # network$params$turgorPressureAtGsMax = NULL
  # network$params$JarvisPAR = sureau_Veg$params$JarvisPAR 
  # network$params$gsNight = sureau_Veg$params$gsNight
  # network$params$Tgs_optim = sureau_Veg$params$Tgs_optim
  # network$params$Tgs_sens = sureau_Veg$params$Tgs_sens
  # network$params$gCrown0 = sureau_Veg$params$gCrown0
  # network$params$gmin_S = sureau_Veg$params$gmin_S
  # network$params$fTRBToLeaf = sureau_Veg$params$fTRBToLeaf
  # network$params$gmin20 = sureau_Veg$params$gmin20
  # network$params$gsMax = sureau_Veg$params$gsMax
  # network$params$VCroot_P50 = NULL #Not used in sureau
  # network$params$VCroot_slope = NULL #Not used in sureau
  # network$params$k_SLApoInit = sureau_Veg$params$k_SLApoInit
  # network$params$k_RSApoInit = sureau_Veg$params$k_RSApoInit #Beware four values in Medfate but only three in Sureau (but still only 3 water potential in Medfate : BUG !!?)
  # network$params$PiFullTurgor_Leaf = sureau_Veg$params$PiFullTurgor_Leaf
  # network$params$epsilonSym_Leaf = sureau_Veg$params$epsilonSym_Leaf
  # network$params$PiFullTurgor_Stem = sureau_Veg$params$PiFullTurgor_Stem
  # network$params$epsilonSym_Stem = sureau_Veg$params$epsilonSym_Stem
  # network$params$C_SApoInit = sureau_Veg$params$C_SApoInit
  # network$params$C_LApoInit = sureau_Veg$params$C_LApoInit
  
  
  network$params$VCleaf_P50 = sureau_Veg$params$P50_VC_Leaf
  network$params$VCleaf_slope = sureau_Veg$params$slope_VC_Leaf
  network$params$VCstem_P50 = sureau_Veg$params$P50_VC_Stem
  network$params$VCstem_slope = sureau_Veg$params$slope_VC_Stem
  
  network$PsiSoil = sureau_Soil$PsiSoil
  network$k_Soil = sureau_Soil$kSoil

    
  network$LAI = sureau_Veg$LAI
  
  network$Psi_LApo = sureau_Veg$Psi_LApo
  network$Psi_LSym = sureau_Veg$Psi_LSym
  network$Psi_SSym = sureau_Veg$Psi_SSym
  network$Psi_SApo = sureau_Veg$Psi_SApo
  network$Psi_SApo_cav = sureau_Veg$Psi_SApo_cav
  network$Psi_LApo_cav = sureau_Veg$Psi_LApo_cav
  
  network$PLC_Stem = sureau_Veg$PLC_Stem
  network$PLC_Leaf = sureau_Veg$PLC_Leaf
  
  network$C_SApo = sureau_Veg$C_SApo
  network$C_LApo = sureau_Veg$C_LApo
  network$C_SSym = sureau_Veg$C_SSym
  network$C_LSym = sureau_Veg$C_LSym
  
  network$k_Plant = sureau_Veg$k_Plant
  network$k_SSym = sureau_Veg$k_SSym
  network$k_LSym = sureau_Veg$k_LSym
  network$k_SLApo = sureau_Veg$k_SLApo
  network$k_RSApo = sureau_Veg$k_RSApo
  network$k_SoilToStem = sureau_Veg$k_SoilToStem
  
  network$Q_LApo_sat_mmol_perLeafArea = sureau_Veg$Q_LApo_sat_mmol_perLeafArea
  network$Q_SApo_sat_mmol_perLeafArea = sureau_Veg$Q_SApo_sat_mmol_perLeafArea
  network$Q_LSym_sat_mmol_perLeafArea = sureau_Veg$Q_LSym_sat_mmol_perLeafArea
  network$Q_SSym_sat_mmol_perLeafArea = sureau_Veg$Q_SSym_sat_mmol_perLeafArea
  
  
  network$Elim =  sureau_Veg$Elim
  
  network$Emin_L =  sureau_Veg$Emin #Emin is Emin leaf in sureau_veg (i.e. WBVeg object)
  network$Emin_S = sureau_Veg$Emin_S
  
  # network$Einst = sureau_Veg$Einst #NULL : Does not exist in sureau_veg (i.e. WBVeg object)
  # network$Elim_SH = sureau_Veg$Elim_SH #NULL : Does not exist in sureau_veg (i.e. WBVeg object)
  # network$Elim_SL = sureau_Veg$Elim_SL #NULL : Does not exist in sureau_veg (i.e. WBVeg object)
  # network$Emin_L_SH = sureau_Veg$Emin_L_SH #NULL : Does not exist in sureau_veg (i.e. WBVeg object)
  # network$Emin_L_SL = sureau_Veg$Emin_L_SL #NULL : Does not exist in sureau_veg (i.e. WBVeg object)
  # network$Diag_nwhile_cavit = sureau_Veg$Diag_nwhile_cavit #NULL : Does not exist in sureau_veg (i.e. WBVeg object)
  # network$Diag_deltaRegulMax = sureau_Veg$Diag_deltaRegulMax #NULL : Does not exist in sureau_veg (i.e. WBVeg object)
  # network$Diag_deltaPLCMax = sureau_Veg$Diag_deltaPLCMax #NULL : Does not exist in sureau_veg (i.e. WBVeg object)
  # network$Diag_timeStepInSeconds = sureau_Veg$Diag_timeStepInSeconds #NULL : Does not exist in sureau_veg (i.e. WBVeg object)
  
  
  return(network)
  
}


# Forest/topo preparation ------------------------------------------------------
# Example with one single cohort of Q. ilex
data("exampleforestMED2")
forest <- exampleforestMED2
forest$treeData <- forest$treeData[2, ,drop = FALSE]
forest$treeData$LAI <- 2.5
forest$shrubData <- forest$shrubData[numeric(0), ,drop = FALSE]
forest$herbHeight <- NA
forest$herbLAI <- NA

latitude <- 41.82592
elevation <- 100
slope <- 0
aspect <- 0


# Soil definition ---------------------------------------------------------
df_soil <- defaultSoilParams(3)
df_soil$widths <- c(300,1200,2500)
df_soil$rfc <- c(40,65,97)
soil <- soil(df_soil)


# Network -----------------------------------------------------------------

#Initialize control parameters
control <- defaultControl("Cochard")
control$subdailyResults <- TRUE
control$cavitationRefill <- "none"
#Initialize input
x2 <- forest2spwbInput(forest, soil, SpParamsMED, control)

network <- initCochardNetworks(x2)[[1]]

veg <- readRDS("Rdata/veg.rds")
soil <- readRDS("Rdata/soil.rds")

network_2 <- transform_sureauToNetwork(network, veg, soil)

opt = c("Lsym" = 1.0,"Ssym" = 1.0,"Eord" = 1.0,"Lcav" = 1.0, "Scav" = 1.0,
        "CLapo" = 1.0,"CTapo" = 1.0)
opt_list = list("Lsym" = 1.0,"Ssym" = 1.0,"Eord" = 1.0,"Lcav" = 1.0, "Scav" = 1.0,
        "CLapo" = 1.0,"CTapo" = 1.0, numericalScheme="Semi-Implicit")

dt = 600
semi_implicit_integration(network_2, dt = dt, opt = opt, cavitationRefill = "none")


sol <- readRDS("Rdata/outoutVegObject.rds")
sol2 <- implicit.temporal.integration.atnp1(veg, soil, dt = dt, opt = opt_list)

network_2$Psi_LSym
sol$Psi_LSym
sol2$Psi_LSym

network_2$Psi_SSym
sol$Psi_SSym
sol2$Psi_SSym

network_2$Psi_LApo
sol$Psi_LApo
sol2$Psi_LApo

network_2$Psi_SApo
sol$Psi_SApo
sol2$Psi_SApo

network_2$Diag_nwhile_cavit
sol$Diag_nwhile_cavit
sol2$Diag_nwhile_cavit

