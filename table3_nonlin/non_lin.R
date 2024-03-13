rm(list = ls())
library(lme4)
library(stats)
library(dplyr)
library(MASS)
library(varTestnlme)
library(saemix)
library(ggplot2)
library(doParallel) 
library(xtable)#

path <- paste0("Out_",Sys.Date())
dir.create(path)
setwd(path)

growth_function<-function(psi,id,xidep) {
  x<-xidep[,1]
  psi1<-psi[id,1]
  psi2<-psi[id,2]
  psi3<-psi[id,3]
  
  f <- psi1/(1+exp(-(x-psi2)/psi3))
  return(f)
}

createDatas3eff <- function(N,beta,gamma,sbruit,name){
  require(mvtnorm)
  phi <- rmvnorm(N,beta,gamma)
  d <- data.frame()
  
  for (i in 1:N)
  {  
    y <- growth_function(phi,i,tj) + rnorm(length(tj),0,sbruit)
    d <- rbind(d,data.frame(y = as.vector(y), days = tj, subject = rep(i,length(y))))
    
  }
  d
}

createResultatList <- function(N_grid, B){
  # le resultat : 
  resultat = list()
  alpha_grid = c(0.01, 0.05, 0.1)
  
  #1) les niveaux empiriques
  level_simu_b = data.frame(matrix(data = 0, nrow = 3, ncol = length(N_grid)))
  level_simu_b_shrink = data.frame(matrix(data = 0, nrow = 3, ncol = length(N_grid)))
  level_simu_asym = data.frame(matrix(data = 0, nrow = 3, ncol = length(N_grid)))
  
  colnames(level_simu_b) = paste('N =', N_grid)
  rownames(level_simu_b) = paste('alpha =', alpha_grid)
  colnames(level_simu_b_shrink) = paste('N =', N_grid)
  rownames(level_simu_b_shrink) = paste('alpha =', alpha_grid)
  colnames(level_simu_asym) = paste('N =', N_grid)
  rownames(level_simu_asym) = paste('alpha =', alpha_grid)
  
  #2) les lrt bootstrap et les lrt 
  lrt_sauv = 0
  lrt_b_sauv = matrix(0, nrow = 1, ncol = B)
  
  #3) les paramètres pour générer les echantillons bootstrap 
  param_sauv = matrix(0, nrow = 1, ncol = 7) # beta de dim 3, gamma1 gamma2 gamma12 et sigma
  
  for (N in N_grid){
    
    resultat[[paste0('N',N)]]$parameters = list()
    resultat[[paste0('N',N)]]$lrt_init = list()
    resultat[[paste0('N',N)]]$lrt_boot = list()
    resultat[[paste0('N',N)]]$lrt_boot_shrink = list()}
  
  resultat$level_boot = level_simu_b
  resultat$level_boot_shrink = level_simu_b_shrink
  resultat$level_asym = level_simu_asym 
  resultat}



K = 1000
B = 300
alpha_grid = c(0.01, 0.05, 0.1)
N_grid = c( 40)
J= 5
tj <- matrix(c(seq(50,1000,length.out = J),1100,1200,1300,1400,1500),ncol=1)
j = length(tj)
beta <- c(200,500, 150)
gamas <- diag(c(10**2,10**2,0))  
sbruit <- 5

latex.out.file ="comp_boot_boots.tex"
result.list.file = "comp_boot_boots.rda"



saemix.modelH0 <- saemixModel(model=growth_function,description="Logistic growth",
                              psi0=matrix(beta,ncol=3,byrow=TRUE,
                                          dimnames=list(NULL,c("psi1","psi2","psi3"))),
                              covariance.model=matrix(c(1,0,0,
                                                        0,1,0,
                                                        0,0,0),ncol=3,byrow=TRUE), verbose = FALSE)

saemix.modelH1 <- saemixModel(model=growth_function,description="Logistic growth",
                              psi0=matrix(beta,ncol=3,byrow=TRUE,
                                          dimnames=list(NULL,c("psi1","psi2", "psi3"))),
                              covariance.model=matrix(c(1,0,0,
                                                        0,1,0,
                                                        0,0,1),ncol=3,byrow=TRUE), verbose = FALSE)
saemix.options <- saemixControl(map = FALSE, print = FALSE, 
                                fim =FALSE ,save=FALSE,save.graphs=FALSE, 
                                nmc.is = 30000, nbiter.burn = 100, 
                                nbiter.saemix = c(600, 500), print.is = FALSE,
                                nb.chains = 2) 
gam0s = matrix(0, nrow = K, ncol = 2)
gam1s = matrix(0, nrow = K, ncol = 3)

level_simu = data.frame(matrix(data = 0, nrow = 3, ncol = length(N_grid)))
colnames(level_simu) = paste('N =', N_grid)
rownames(level_simu) = paste('alpha =', alpha_grid)
boot_table = level_simu
boots_table = level_simu
asym_table =  level_simu




stime <- system.time({
  pdf(NULL)
  
  #ON RENTRE DANS LA BOUCLE k
  resultat = createResultatList(N_grid, B)
  
  for (k in 1:K){     
    for (N in N_grid){
      pdf(NULL)
      #on simule les donn?es
      graine = k+10000
      set.seed(graine)
      d0 = createDatas3eff(N,beta,gamas,sbruit,"d0")
      saemix.dataH0 <- saemixData(name.data = d0,header=TRUE,name.group=c("subject"),
                                  name.predictors=c("days"),name.response=c("y"),
                                  units=list(x="days",y="kg"), verbose = FALSE)
      
      
      fity0h0 <- saemix(saemix.modelH0,saemix.dataH0, saemix.options)
      fity0h1 <- saemix(saemix.modelH1,saemix.dataH0, saemix.options)
      
      gam0s[k,] = diag(fity0h0@results@omega)[1:2]
      gam1s[k,] = diag(fity0h1@results@omega)
      
      #param?tres 
      llH1_is_n <- fity0h1@results@ll.is
      llH0_is_n <- fity0h0@results@ll.is
      lrt = -2 * (llH0_is_n - llH1_is_n)
      if (lrt>0){
      
        betaChap1_n <- fity0h1@results@betas[,1]
        omegaChap1_n <- fity0h1@results@omega
        
        gama_hat = as.matrix(omegaChap1_n)
        gama_hat[3,] = 0
        gama_hat[,3] = 0
        
        
        
        sigma1_n <- fity0h1@results@respar[1]
        
        
        
        param <- list(llH1_is=llH1_is_n,llH0_is=llH0_is_n,betaChap1=betaChap1_n,omegaChap1=omegaChap1_n,sigma1=sigma1_n)
        resultat$parameters[[k]] = param
        lrt_B = vector(mode = 'numeric', B)
        lrt_Bs = vector(mode = 'numeric', B)
        rm(d0)
        ncore = 50
        grappe = makeCluster(ncore-1) 
        registerDoParallel(grappe)
        
        Bootstrap = foreach(b=1:B, .packages = c("lme4", "stats", "dplyr", "MASS", "varTestnlme", "saemix")) %dopar% {
          pdf(NULL)
          gr = (b*(k+300000))
          set.seed(gr)
          # On g?n?re selon le processus de g?n?ration de donn?es ?tabli (restricted parametric bootstrap)
          db = createDatas3eff(N,betaChap1_n,gama_hat,sigma1_n,"db")
          write.table(db, file = "db.rds")
          saemix.databootH0 <- saemixData(name.data = db,header=TRUE,name.group=c("subject"),
                                          name.predictors=c("days"),name.response=c("y"),
                                          units=list(x="days",y="kg"), verbose = FALSE)
         
          # On calcule la statistique de test sur l'?chantillon Bootstrap
          m0_b <- saemix(saemix.modelH0,saemix.databootH0, saemix.options)
          m1_b <- saemix(saemix.modelH1,saemix.databootH0, saemix.options)
          l1_b = m1_b@results@ll.is
          l0_b = m0_b@results@ll.is
          lrt_B = -2 * (l0_b - l1_b)
          if (lrt_B<0){lrt_B = 0}
          
          
          rm(list=c("db"))#rm(list=c("db","dbs"))
          dev.off()
          res_boot = list("lrt_b" = lrt_B)#res_boot = list("lrt_b" = lrt_B, "lrt_bs" = lrt_Bs)
          res_boot
        }
        stopCluster(grappe)
        for (b in 1:B){
          lrt_B[b] = Bootstrap[[b]]$lrt_b
        }
        resultat[[paste0('N',N)]]$lrt_boot[[k]] = lrt_B

        
        p_b = mean(lrt_B > lrt)

        for (alpha in alpha_grid){
          
          resultat$level_boot[paste('alpha =', alpha), paste('N =', N)] = resultat$level_boot[paste('alpha =', alpha), paste('N =', N)] + ( p_b < alpha)/K

          test_asym = suppressMessages(varTestnlme::varCompTest(fity0h1,fity0h0))
          p_asym = test_asym$p.value[[3]]
          resultat$level_asym[paste('alpha =', alpha), paste('N =', N)] =  resultat$level_asym[paste('alpha =', alpha), paste('N =', N)] + (p_asym < alpha)/K
      }}
      
      boot_table = resultat$level_boot 
      asym_table = resultat$level_asym
      
      write.table(boot_table *(K/k),file = "boot_int.csv")
      write.table(asym_table*(K/k), file = "asym_int.csv")
      write.table(paste("k = ", k), file = "suivi.csv")
      
    } #fin dela boucle K
    pdf(NULL)
    while (!is.null(dev.list())){dev.off(dev.list()[1])}
    
    
  }#fin de dopar
  
  while (!is.null(dev.list())){dev.off(dev.list()[1])}
  
})
stime

pdf(NULL)
# stockage des resultats




tableau = matrix(0, nrow = length(alpha_grid), ncol = 2*dim(boot_table)[2]) 
rownames(tableau) = paste ("alpha=", alpha_grid)



tableau2 = matrix(0, nrow = length(alpha_grid), ncol = 2*dim(boot_table)[2]) 
rownames(tableau) = paste ("alpha=", alpha_grid)

for (c in 1:dim(boot_table)[2]){
  tableau2[,(2*c-1)] = paste0(round(boot_table[,c],3), sprintf("(%.3f)", sqrt(boot_table[,c] * (1-boot_table[,c])/K)))
  tableau2[,2*c] = paste0(round(asym_table[,c],3), sprintf("(%.3f)", sqrt(asym_table[,c] * (1-asym_table[,c])/K)))
  
  
  pdf(NULL)
  
print(xtable(tableau2, type = "latex"), file = "niveaux_non_lin_asym_boot.tex")}

saveRDS(resultat, "ResTotal.rda")
write.table(boot_table,file = "boot.csv")
write.table(asym_table, file = "asym.csv")
write.table(gam0s, file = "gam0s.csv")
write.table(gam1s, file = 'gam1s.csv')

print(paste("K=",K,"B=",B, "J=" , j,"sbruit=",sbruit))
print(paste("beta=", beta))
print(paste("gamas=", gamas))