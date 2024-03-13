# Tester une variance est nulle dans un modèle à 2 effets aléatoires (M2)
# Tom Guédon
rm(list=ls())


library(lme4)
library(stats)
library(dplyr)
library(MASS)
library(varTestnlme)
library(xtable)
library(doParallel)   

path <- paste("Out_",Sys.Date(),format(Sys.time(), "_%Hh_%Mm_%Ss"),sep="")
dir.create(path)
setwd(path)

# fonctions utiles

createResultatList_power <- function(vliste, covliste , B){
  # le resultat : 
  resultat = list()
  nv = length(vliste)
  ncov = length(covliste)
  #1) les niveaux empiriques
  level_simu_b = data.frame(matrix(data = 0, nrow = nv * ncov, ncol = 1))
  level_simu_asym = data.frame(matrix(data = 0, nrow = nv * ncov, ncol = 1))
  k = 0

  rnam = rep(" ", nv*ncov)
for (v in vliste){for (co in covliste){k = k+1
  rnam[k] = paste("v2 = ", v, ", v12=", co)
  }}
  
  rownames(level_simu_b) = rnam
  rownames(level_simu_asym) = rnam
  
  #2) les lrt bootstrap et les lrt 
  lrt_sauv = 0
  #lrt_b_sauv = matrix(0, nrow = 1, ncol = B)
  
  #3) les paramètres pour générer les echantillons bootstrap 
  param_sauv = matrix(0, nrow = 1, ncol = 4) # beta de dim 2, gamma1, gamma2, gamma12 et sigma
  
  for (nam in rnam){
    
    resultat[[nam]]$parameters = param_sauv
    resultat[[nam]]$lrt_init = lrt_sauv
    }
  
  resultat$level_boot = level_simu_b
  resultat$level_asym = level_simu_asym 
  resultat}

data2effLin <- function(N, J, beta, cov, sbruit){ ## Simulation des donn?es
  
  # On souhaite simuler y = X * beta + Z1 * b1 + epsilon
  x = cbind(rep(1,N*J),rep(1:J, N))
  z1= matrix(0, nrow = N*J, ncol = N)
  z2= matrix(0, nrow = N*J, ncol = N)
  
  
  for (i in 1:N){z1[((i - 1) * J + 1) : (i * J), i] = rep(1 , J)}
  for (i in 1:N){z2[((i - 1) * J + 1) : (i * J), i] = 1:J}
  
  # On construit un vecteur "indiv" pour les individus, afin de l'utiliser dans le package lme4
  indiv = rep('0', N * J)
  for (i in 1:N){indiv[((i - 1) * J + 1) : (i * J)] = as.character(i)}
  
  #on simule les residus et les effets aléatoires
  epsilon = rnorm(N*J, mean = 0, sd = sbruit)
  rdm_effects = mvrnorm(N, c(0,0), cov)
  
  # On construit la variable reponse
  y = x %*% beta + z1 %*% rdm_effects[,1] + z2 %*% rdm_effects[,2] +  epsilon
  
  #On construit le dataset qui sera utilisé pour ajuster le modèle
  data = cbind(y,data.frame(x[,1], x[,2], indiv))
  colnames(data) = c('y', 'intercept', 'x1', 'indiv')
  data}

## paramètres de simulations 

# K nombre de répétitions de l'expérience
# B nombre d'échantillons bootstrap simulés à chaque expérience
K = 2500
B = 500


## paramètres du modèle 

# N_grid vecteur contenant les differents nombres d'individus qui seront utilisés
# J nombre de mesures répétées par individu
# beta vecteur des effets fixes
## s1 écart type des effets aléatoires sous H1 (simulations pour les puissances) ( a approfondir car normalement lpusieurs valeurs )
# s0 écart type des effets aleatoires sous H0 (simulations pour les niveaux)
# sbruit ecart type du bruit homoscedastique
# matrice de covariance des effets aleatoires

N = 30 
J = 5 
beta = c(0, 7)
s0 = c(sqrt(1.3),0)

sbruit = sqrt(1.5)
cov0 = diag(s0**2)

cov1 = c( 0.005, 0.02, 0.05)
cov12 = c(0, 0.2, 0.8)


level_simu_b = data.frame(matrix(data = 0, nrow = length(cov1)*length(cov12), ncol = 1))
level_simu_asym = data.frame(matrix(data = 0, nrow = length(cov1)*length(cov12), ncol = 1))

ncore = 40
grappe = makeCluster(ncore)
registerDoParallel(grappe)

stime <- system.time({
  resTotal <- foreach(k=1:K, .packages = c("lme4", "stats", "dplyr", "MASS", "varTestnlme")) %dopar% {
    resultat = createResultatList_power(cov1, cov12, B)
    for (c1 in cov1){
      for (c12  in cov12){
      set.seed(k*10000)
      cov0 = matrix(c(1.3, c12*sqrt(1.3*c1) , c12*sqrt(1.3*c1), c1),nrow =2, byrow = TRUE)
      d0 = data2effLin(N, J, beta, cov0, sbruit)
      # On ajuste le model sous H0...
      m0 = suppressMessages(lme4::lmer(formula = 'y ~ 1 + x1 +  (1  | indiv)', data = d0, REML = F ))
      # On calcule la statistique de test sur le jeu de donn?es initiales
      m1 = suppressMessages(lme4::lmer(formula = 'y ~ 1 + x1 + (1 + x1 | indiv)', data = d0, REML = F))
      
      #... et r?cup?ration des param?tres n?cessaires au processus de g?n?ration de donn?es 
      gam_hat = as.data.frame(VarCorr(m1))$vcov[1]
      #On projette sur H0
      shat = c(gam_hat, 0)
      covhat = diag(shat)
      #gam_tilde = gam_hat * (gam_hat > cn)
      beta_hat = as.vector(fixef(m1))
      sigma_hat = as.data.frame(VarCorr(m1))$vcov[4]
      sbruit_hat = sigma_hat**0.5
      theta = c(beta_hat, gam_hat,sigma_hat)
      resultat[[paste("v2 = ", c1, ", v12=", c12)]]$parameters = theta
      
      
      l1 = logLik(m1)
      l0 = logLik(m0)
      lrt = -2 * (l0 - l1)
      resultat[[paste("v2 = ", c1, ", v12=", c12)]]$lrt_init = lrt
      
      # lancement de la proc?dure Bootstrap
      
      # On va stocker les statistiques de test
      lrt_B = vector(mode = 'numeric', B)
      
      for (b in 1:B){
        set.seed((k+9999)*100*b)
        # On g?n?re selon le processus de g?n?ration de donn?es ?tabli (restricted parametric bootstrap)
        db = data2effLin(N, J, beta_hat, covhat, sbruit_hat)
        
        # On calcule la statistique de test sur l'?chantillon Bootstrap
        m1_b = suppressMessages(lme4::lmer(formula = 'y ~ 1 + x1 + (1+ x1  | indiv)', data = db, REML = F))
        m0_b = suppressMessages(lme4::lmer(formula = 'y ~ 1 + x1 +  ( 1 | indiv)', data = db, REML = F))
        l1_b = logLik(m1_b)
        l0_b = logLik(m0_b)
        lrt_B[b] = -2 * (l0_b - l1_b)
        
      }
      

      p_b = mean(lrt_B > lrt)
      resultat$level_boot[paste("v2 = ", c1, ", v12=", c12),] = p_b < 0.05
      test_asym = suppressMessages(varTestnlme::varCompTest(m1,m0))
      p_asym = test_asym$p.value[[3]]
      resultat$level_asym[paste("v2 = ", c1, ", v12=", c12),] = p_asym < 0.05
      
    }}

    resultat
  }
})  
stopCluster(grappe)


stime


# stockage des resultats

level_simu_b = data.frame(matrix(data = 0, nrow = length(cov1)*length(cov12), ncol = 1))
level_simu_asym = data.frame(matrix(data = 0, nrow = length(cov1)*length(cov12), ncol = 1))

resultat = list()

k = 0
nv = length(cov1)
ncov = length(cov12)
rnam = rep(" ", nv*ncov)
for (v in cov1){for (co in cov12){k = k+1
rnam[k] = paste("v2 = ", v, ", v12=", co)
}}

rownames(level_simu_b) = rnam
rownames(level_simu_asym) = rnam

#2) les lrt bootstrap et les lrt 
lrt_sauv = rep(0, K)

#3) les paramètres pour générer les echantillons bootstrap 
param_sauv = matrix(0, nrow = K, ncol = 4) # beta de dim 2, gamma1, gamma2, gamma12 et sigma


for (r in rnam){
  
  resultat[[r]]$parameters = param_sauv
  resultat[[r]]$lrt_init = lrt_sauv
}

resultat$level_boot = level_simu_b
resultat$level_asym = level_simu_asym 

for (k in 1:K){
  for (r in  rnam){
    resultat[[r]]$parameters[k,] = resTotal[[k]][[r]]$parameters
    resultat[[r]]$lrt_init[k] = resTotal[[k]][[r]]$lrt_init
  }
  resultat$level_boot = resultat$level_boot + resTotal[[k]]$level_boot / K
  resultat$level_asym = resultat$level_asym + resTotal[[k]]$level_asym / K
}

print(xtable(resultat$level_boot, type = "latex"), file = "puissance_non_lin_boot.tex")
print(xtable(resultat$level_asym, type = "latex"), file = "puissance_non_lin_asym.tex")
saveRDS(resultat, "power_m2_out.rda")
