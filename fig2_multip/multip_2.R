rm(list=ls())
require(lme4)
require(stats)
require(dplyr)
library(mvtnorm)

require(MASS)
library(xtable)
library(doParallel) 




path <- paste0("Out_",Sys.Date())
dir.create(path)
setwd(path)
set.seed(12345)



#param  etres du mod  ele
N_grid = c(40)
N = 40 
J = 9
alpha_grid = c(0.01,0.05,0.1)


#tous les effets ont la m eme variance gama
sigma = 1


# nombre d'effets al  eatoires
M = 8
var = rep(1,M)
seuil = 0.5 * 40**{-0.2}


for(cN in c(0, seuil, 0.9)){
  for (s in 2:5){
    # nombre de param  etres de nuisances
    
    
    var[1 : (s+1)] = 0
    gama = diag(var)
    
    # matrice indicatrice des individus
    z= matrix(0, nrow = N*J, ncol = N)
    for (i in 1:N){z[((i - 1) * J + 1) : (i * J), i] = rep(1 , J)}
    
    
    
    
    K=2500
    B = 300
    
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
      rdm_effects = rmvnorm(N, rep(0,M), gama)
      
      
      #simulation des donn  ees selon le mod  ele
      y = rep(0, N*J)
      epsilon = rnorm(N*J, mean = 0, sd = sigma)
      for(i in 1:N){
        res = rep(0,J)
        for (m in 1:M){
          res = res + X[[i]][m,] * rdm_effects[i,m]}
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
    
    # Bootstrap 

    level_simu_bs = rep(0,3)
    
    
    #stime
    resultat = list()
    level_simu_b_shrink = as.data.frame(matrix(data = 0, nrow = 3, ncol = length(N_grid)))
    colnames(level_simu_b_shrink) = paste('N =', N_grid)
    rownames(level_simu_b_shrink) = paste('alpha =', alpha_grid)

    resultat$level_boot_shrink = level_simu_b_shrink
    
    stime <- system.time({
      
      for (k in 1:K){
        #set.seed(k+1000)
        
        
        for (N in N_grid){
          d0 = data_multip(N, J, sigma, gama, M, s, ind_shrink)
          # On ajuste le model sous H0...
          m0 = suppressMessages(lme4::lmer(formula = 'y~(0+x2+x3+x4+x5+x6+x7+x8  ||indiv)', data = d0, REML = F))
          # On calcule la statistique de test sur le jeu de donn?es initiales
          
          m1 = suppressMessages(lme4::lmer(formula = 'y~(0+x1+x2+x3+x4+x5+x6+x7+x8  ||indiv)', data = d0, REML = F))
          
          #... et r?cup?ration des param?tres n?cessaires au processus de g?n?ration de donn?es 
          d_gam = as.data.frame(VarCorr(m1))$vcov[1:8]
          # projection su H0
          d_gam[1]=0
          
          
          d_gam[d_gam<cN] = 0#seuillage des avriances   
          gam_hat = diag(d_gam) 
          
          sigma_hat = as.data.frame(VarCorr(m1))$vcov[9]

          
          
          sbruit_hat = sigma_hat**0.5

          l1 = logLik(m1)
          l0 = logLik(m0)
          lrt = -2 * (l0 - l1)
          resultat$lrt[k] = lrt
          # lancement de la proc?dure Bootstrap
          
          # On va stocker les statistiques de test
          lrt_Bs = vector(mode = 'numeric', B)
          
          
          ncore = 75
          grappe = makeCluster(ncore-1) 
          registerDoParallel(grappe)
          
          resTotal <- foreach(b=1:B, .packages = c("lme4", "stats", "dplyr", "MASS","mvtnorm")) %dopar% {
            
            # On g?n?re selon le processus de g?n?ration de donn?es ?tabli (restricted parametric bootstrap)
            dbs = data_multip(N, J, sigma_hat, gam_hat, M, s, ind_shrink)


            
            
            m1_bs = suppressMessages(lme4::lmer(formula = 'y~(0+x1+x2+x3+x4+x5+x6+x7+x8  ||indiv)', data = dbs, REML = F))
            m0_bs = suppressMessages(lme4::lmer(formula = 'y~(0+x2+x3+x4+x5+x6+x7+x8  ||indiv)', data = dbs, REML = F))
            l1_bs = logLik(m1_bs)
            l0_bs = logLik(m0_bs)
            lrt_bs = -2 * (l0_bs - l1_bs)
            return(list("lrt_bs" = lrt_bs))
          }
          stopCluster(grappe)
          for (b in 1:B){
            lrt_Bs[b] = resTotal[[b]]$lrt_bs}
          
          p_bs = mean(lrt_Bs > lrt)
          resultat$p_Bs[k] = p_bs
          for (alpha in alpha_grid){
            
            resultat$level_boot_shrink[paste('alpha =', alpha), paste('N =', N)] = resultat$level_boot_shrink[paste('alpha =', alpha), paste('N =', N)] + (p_bs < alpha)
          }
        }
        #resultat
      
        write.table(100 *resultat$level_boot_shrink/k, file = paste0("boot_shrink_int","s",s,"cN",cN,".csv"))
        write.table(paste("k=",k), file = paste0("suivi","s",s,"cN",cN,".csv"))
        
        
      }
    })  
    stime
    

    boots_table = resultat$level_boot_shrink
    
    tableau = matrix(0, nrow = length(alpha_grid), ncol = 2*dim(boots_table)[2]) 
    #colnames(tableau) = paste("N=", N_grid)
    rownames(tableau) = paste ("alpha=", alpha_grid)
    
    for (c in 1:dim(boots_table)[2]){
      tableau[,2*c] = paste0(round(boots_table[,c],3)*100, sprintf("(%.3f)", sqrt(boots_table[,c] * (1-boots_table[,c])/K)))
    }
    
    print(xtable(tableau, type = "latex"), file = paste0("s",s,"cN",cN,"niveaux.tex"))
    
    saveRDS(resultat, file = paste0("ResTotal","s",s,"cN",cN,".rda"))
    
    write.csv2(boots_table *100 , file = paste0("s",s,"cN",cN,"niveaux.csv"))}}

