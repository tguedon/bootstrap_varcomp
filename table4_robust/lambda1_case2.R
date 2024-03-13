rm(list=ls())
require(lme4)
require(stats)
require(dplyr)
library(mvtnorm)

require(MASS)
library(xtable)
library(doParallel) 
set.seed(12345)

dir.create( "resultats")
setwd(dir = "resultats")
data3effLin <- function(N, J, beta, cov, sbruit){ ## Simulation des donn?es
  
  # On souhaite simuler y = X * beta + Z1 * b1 + epsilon
  x = cbind(rep(1,N*J),rep(1:J, N), rep(1:J,N)**2)
  z1= matrix(0, nrow = N*J, ncol = N)
  z2= matrix(0, nrow = N*J, ncol = N)
  z3= matrix(0, nrow = N*J, ncol = N)
  
  
  for (i in 1:N){z1[((i - 1) * J + 1) : (i * J), i] = rep(1 , J)}
  for (i in 1:N){z2[((i - 1) * J + 1) : (i * J), i] = 1:J}
  for (i in 1:N){z3[((i - 1) * J + 1) : (i * J), i] = (1:J)**2}
  
  # On construit un vecteur "indiv" pour les individus, afin de l'utiliser dans le package lme4
  indiv = rep('0', N * J)
  for (i in 1:N){indiv[((i - 1) * J + 1) : (i * J)] = as.character(i)}
  
  #on simule les residus et les effets aléatoires
  epsilon = rnorm(N*J, mean = 0, sd = sbruit)
  rdm_effects = mvrnorm(N, c(0,0,0), cov)
  
  # On construit la variable reponse
  y = x %*% beta + z1 %*% rdm_effects[,1] + z2 %*% rdm_effects[,2] +  epsilon
  
  #On construit le dataset qui sera utilisé pour ajuster le modèle
  data = cbind(y,data.frame(x[,1], x[,2], x[,3], indiv))
  colnames(data) = c('y', 'intercept', 'x1', 'x2', 'indiv')
  data}

ncore = min(64,detectCores())
print(ncore)

grappe = makeCluster(ncore-1) 
registerDoParallel(grappe)
N = 30 
J = 6
alpha = 0.05
#tous les effets ont la m eme variance gama
sigma = 1.5
# nombre d'effets al  eatoires
M = 3
v_ind = 0
vlist = c(0.001,0.01,0.1)

propshrink = matrix(0,nrow = length(vlist), ncol = M)

niv = rep(0,length(vlist))
beta = c(0,7,3)

for(v in vlist){
  
  v_ind = v_ind+1
  var = c(1.3**2,v,0)
  var[M] = 0
  cN = v
  
  # H0
  
  
  
  
  
  
  # var non nuls
  
  #var[1:(M-1)] = seq(from = 0.5, to = 0, length.out = M-1)
  
  gama = diag(var)
  
  # matrice indicatrice des individus
  z= matrix(0, nrow = N*J, ncol = N)
  for (i in 1:N){z[((i - 1) * J + 1) : (i * J), i] = rep(1 , J)}
  
  
  # seuils 
  
  
  K=2500
  B = 500
  data_multip = function(N, J, sigma, gama, M, s, ind_shrink){
    # cr  eation de la matrice de covariables
    X = list()
    
    for (i in 1:N){
      X[[i]] = matrix(0, nrow = M, ncol = J)
      for (j in 1:J){
        
        for (m in 1 : M){
          X[[i]][m,j] = rnorm(1, mean = 2, sd= 0.5)
        }
      }
    }
    #construction de la matrice de covariable   e partir des covariables
    covariables = matrix(0,nrow = N*J, ncol = M)
    for(i in 1:N){
      for (j in 1:J){
        for (m in 1:M){
          covariables[((i-1)*J+j),m] = X[[i]][m,j]
        }
      }
    }
    indiv = rep('0', N * J)
    # matrice des effets al  eatoires
    #rdm_effects = matrix(rnorm(N*M), nrow = M, ncol = N)# ICI mettre un mvnorm(N, rep(0,M), gama) comme d'habitude 
    rdm_effects = rmvnorm(N, rep(0,M), gama)
    # indices des effets al  eatoires param  etres de nuisance
    #rdm_effects[ind_shrink,] = 0 ca ?a va dans ind_shrink
    
    #rdm_effects[1,] = 0 ca aussi 
    
    #simulation des donn  ees selon le mod  ele
    y = rep(0, N*J)
    epsilon = rnorm(N*J, mean = 0, sd = sigma)
    for(i in 1:N){
      res = rep(0,J)
      for (m in 1:M){
        res = res + X[[i]][m,] * (1+rdm_effects[i,m])} #beta =1
      y[((i-1)*J + 1 ): (i*J)] = res
    }
    
    y = y + epsilon 
    
    for (i in 1:N){indiv[((i - 1) * J + 1) : (i * J)] = as.character(i)}
    D = as.data.frame(cbind(y, covariables, indiv))
    colname = c("y",paste0("x",(1:M)), "indiv")
    colnames(D) = colname
    for (j in 1:dim(D)[2]){
      D[,j] = as.numeric(D[,j])}
    
    D
  }
  
  #####
  #sauvegarder les resultats : 
  
  niveaux = 0
  
  t1 = Sys.time()
  #####
  
  for (k in 1:K){
    #set.seed(k+1000)
    if(k%%10==0){print(paste("var",v,"k=",k))}
    
    d0 = data3effLin(N, J, beta, gama, sqrt(sigma))
    
    # On ajuste le model sous H0...
    m0 = suppressMessages(lme4::lmer(formula = 'y~1+x1+x2+(1+x1 ||indiv)', data = d0, REML = F))
    # On calcule la statistique de test sur le jeu de donn?es initiales
    
    m1 = suppressMessages(lme4::lmer(formula = 'y~1+x1+x2+(1+x1+x2  ||indiv)', data = d0, REML = F))
    
    #... et r?cup?ration des param?tres n?cessaires au processus de g?n?ration de donn?es 
    
    d_gam = as.data.frame(VarCorr(m1))$vcov[1:M]
    d_gam[M] = 0
    #indice_shrink = d_gam<cN
    #propshrink[v_ind,] = propshrink[v_ind,] + indice_shrink
    d_gam_shrink = d_gam
    d_gam_shrink[2] = 0
    gam_hat = diag(d_gam_shrink) 
    beta_hat = as.vector(fixef(m1))
    sigma_hat = as.data.frame(VarCorr(m1))$vcov[M+1]
    
    
    
    sbruit_hat = sigma_hat**0.5
    
    l1 = logLik(m1)
    l0 = logLik(m0)
    lrt = -2 * (l0 - l1)
    # lancement de la proc?dure Bootstrap
    
    # On va stocker les statistiques de test
    lrt_Bs = vector(mode = 'numeric', B)
    
    
    
    
    resTotal <- foreach(b=1:B, .packages = c("lme4", "stats", "dplyr", "MASS","mvtnorm")) %dopar% {
      
      # On g?n?re selon le processus de g?n?ration de donn?es ?tabli (restricted parametric bootstrap)
      dbs = data3effLin(N, J, beta_hat, gam_hat, sbruit_hat)
      
      
      
      m1_bs = suppressMessages(lme4::lmer(formula = 'y~1+x1+x2+(1+x1+x2  ||indiv)', data = dbs, REML = F))
      m0_bs = suppressMessages(lme4::lmer(formula = 'y~1+x1+x2+(1+x1 ||indiv)', data = dbs, REML = F))
      l1_bs = logLik(m1_bs)
      l0_bs = logLik(m0_bs)
      lrt_bs = -2 * (l0_bs - l1_bs)
      return(list("lrt_bs" = lrt_bs))
    }
    
    for (b in 1:B){
      lrt_Bs[b] = resTotal[[b]]$lrt_bs}
    
    p_bs = mean(lrt_Bs > lrt)
    niveaux = niveaux+ (p_bs < alpha)
    
    #resultat
    
    #write.table(niveaux/k, file = paste0("niveaux",".csv"))
    #write.table(paste("k=",k), file = paste0("suivi",".csv"))
    
    
  }
  t2 = Sys.time()
  print(paste(round((t2-t1), 2), "min"))
  niv[v_ind] = niveaux
  write.csv2(niveaux/K, file = paste0("levels_",v,".csv"))
}

stopCluster(grappe)
write.csv2(niv,file = "niveaux.csv")
write.csv2(propshrink, file = "propshrink.csv")
