##### Article 

rm(list=ls())
set.seed(156354)
library(lme4)
library(stats)
library(dplyr)
library(MASS)
library(varTestnlme)
library(saemix)
library(ggplot2)
library(doParallel) 
library(xtable)


coucal <- read.table("growth_coucal.csv",sep=",",header = TRUE)


path <- paste0("donnees_reelles_",Sys.Date())
dir.create(path)
setwd(path)



growth_function<-function(psi,id,xidep) {
  #x<-xidep[,1]
  x<-xidep[,1]
  psi1<-psi[id,1]
  psi2<-psi[id,2]
  psi3<-psi[id,3]
  
  f <- psi1/(1+exp(-(x-psi2)/psi3))
  return(f)
}

growth_function2<-function(psi,id,xidep) {
 
  x<-xidep
  psi1<-psi[id,1]
  psi2<-psi[id,2]
  psi3<-psi[id,3]
  
  f <- psi1/(1+exp(-(x-psi2)/psi3))
  return(f)
}

createDatas3eff_rd <- function(N,beta,gamma,sbruit,name){
  require(mvtnorm)
  phi <- rmvnorm(N,beta,gamma)
  d <- data.frame()
  
  for (i in 1:N)
  { id = list_indiv[i]
    tj = id_indiv[id_indiv$nestling_ID2_cat==id,2][!is.na(id_indiv[id_indiv$nestling_ID2_cat==id,2])]
    y <- growth_function2(phi,i,tj) + rnorm(length(tj),0,sbruit)
    d <- rbind(d,data.frame(y = as.vector(y), days = tj, subject = rep(id,length(y))))
    
  }
  #write.table(d,paste(name,"_",N,"_",".rds",sep=""))
  d
}

bc <- coucal[coucal$species_num==0,]
vect = !is.na(bc$weight)
bc = bc[vect,]
#bc <- coucal[coucal$species_num==1,]
id_indiv = bc[,c("nestling_ID2_cat","age")]
list_indiv = levels(as.factor(id_indiv$nestling_ID2_cat))
vect = !is.na(bc$age)
bc_clean = bc[vect,]
B = 1000
#p <- ggplot(data=coucal,aes(x=age,y=weight,group=nestling_ID2_cat)) + geom_point() + geom_line() + xlab("Age (days)") + ylab("Weight (g)")
#ggsave("coucal.pdf",p, width=10, units="cm")

#p <- ggplot(data=coucal[coucal$species_num==0,],aes(x=age,y=weight,group=nestling_ID2_cat)) + geom_point() + geom_line() + xlab("Age (days)") + ylab("Weight (g)")
#ggsave("coucal_BC.pdf",p, width=10, units="cm")

#p <- ggplot(data=coucal[coucal$species_num==1,],aes(x=age,y=weight,group=nestling_ID2_cat)) + geom_point() + geom_line() + xlab("Age (days)") + ylab("Weight (g)")
#ggsave("coucal_WBC.pdf",p, width=10, units="cm")

# choix valeurs initiales
x <- seq(0,20)
y <- 130/(1+exp(-(x-10)/4))
dinit <- data.frame(age=x,weight=y,nestling_ID2_cat="mean")


beta_init <- c(120,10,4)

N = length(list_indiv)


saemix.modelH1 <- saemixModel(model=growth_function,description="Logistic growth",
                              psi0=matrix(beta_init,ncol=3,byrow=TRUE,
                                          dimnames=list(NULL,c("psi1","psi2","psi3"))),transform.par=c(0,0,0),fixed.estim=c(1,1,1),
                              covariance.model=matrix(c(1,0,0,
                                                        0,1,0,
                                                        0,0,1),ncol=3,byrow=TRUE),
                              omega.init=diag(0.1*beta_init),error.model="constant")



##### model .1 test var1,var2 = 0 contre var1, var2>0

saemix.modelH0.1 <- saemixModel(model=growth_function,description="Logistic growth",
                              psi0=matrix(beta_init,ncol=3,byrow=TRUE,
                                          dimnames=list(NULL,c("psi1","psi2","psi3"))),transform.par=c(0,0,0),fixed.estim=c(1,1,1),
                              covariance.model=matrix(c(1,0,0,
                                                        0,0,0,
                                                        0,0,0),ncol=3,byrow=TRUE),
                              omega.init=diag(0.1*beta_init),error.model="constant")

##### model .2 test var 3 = 0 contre var 3 > 0 

saemix.modelH0.2 <- saemixModel(model=growth_function,description="Logistic growth",
                                psi0=matrix(beta_init,ncol=3,byrow=TRUE,
                                            dimnames=list(NULL,c("psi1","psi2","psi3"))),transform.par=c(0,0,0),fixed.estim=c(1,1,1),
                                covariance.model=matrix(c(1,0,0,
                                                          0,1,0,
                                                          0,0,0),ncol=3,byrow=TRUE),
                                omega.init=diag(0.1*beta_init),error.model="constant")



saemix.modelH0.3 <- saemixModel(model=growth_function,description="Logistic growth",
                                psi0=matrix(beta_init,ncol=3,byrow=TRUE,
                                            dimnames=list(NULL,c("psi1","psi2","psi3"))),transform.par=c(0,0,0),fixed.estim=c(1,1,1),
                                covariance.model=matrix(c(1,0,0,
                                                          0,0,0,
                                                          0,0,1),ncol=3,byrow=TRUE),
                                omega.init=diag(0.1*beta_init),error.model="constant")





saemix.options <- saemixControl(nbiter.saemix=c(1000,1000),nb.chains=5,print.is = FALSE,print=FALSE, save=FALSE,save.graphs=FALSE, displayProgress = FALSE, fix.seed = TRUE, seed = 123456789, nmc.is = 50000, nu.is = 4)
#saemix.options <- saemixControl(nbiter.saemix=c(100,100),nb.chains=2,print.is = FALSE,print=FALSE, save=FALSE,save.graphs=FALSE, displayProgress = FALSE, fix.seed = TRUE, seed = 123456789, nmc.is = 2000, nu.is = 4)


saemix.data <- saemixData(name.data=bc,header=TRUE,name.group=c("nestling_ID2_cat"),
                          name.predictors=c("age"),name.response=c("weight"),
                          units=list(x="days",y="g"))

fity0h1 <- saemix(saemix.modelH1,saemix.data,saemix.options)

fity0h0.1 <- saemix(saemix.modelH0.1,saemix.data,saemix.options)

fity0h0.2 <- saemix(saemix.modelH0.2,saemix.data,saemix.options)

fity0h0.3 <- saemix(saemix.modelH0.3,saemix.data,saemix.options)

# vraisemblance

llH1_is_n <- fity0h1@results@ll.is

llH0_is_n.1 <- fity0h0.1@results@ll.is
lrt.1 = -2 * (llH0_is_n.1 - llH1_is_n)

llH0_is_n.2 <- fity0h0.2@results@ll.is
lrt.2 = -2 * (llH0_is_n.2 - llH1_is_n)

llH0_is_n.3 <- fity0h0.3@results@ll.is
lrt.3 = -2 * (llH0_is_n.3 - llH1_is_n)

#param?tres 


gam1s = diag(fity0h1@results@omega)

gam0s.1 = diag(fity0h0.1@results@omega)[1:2]

gam0s.2 = diag(fity0h0.2@results@omega)[1:2]

gam0s.3 = diag(fity0h0.3@results@omega)[1:2]

betaChap1_n <- fity0h1@results@betas[,1]

omegaChap1_n <- fity0h1@results@omega
  
gama_hat.1 = as.matrix(omegaChap1_n)
gama_hat.2 = as.matrix(omegaChap1_n) 
gama_hat.3 = as.matrix(omegaChap1_n) 

sigma1_n <- fity0h1@results@respar[1]
  
# Projection sur Theta_0
  
gama_hat.1[2:3,] = 0

gama_hat.1[,2:3] = 0

gama_hat.2[3,] = 0

gama_hat.2[,3] = 0

gama_hat.3[2,] = 0

gama_hat.3[,2] = 0
  
# Seuillage ou non seuillage

gama_hat.2.1 = gama_hat.2 # pas de seuillage

gama_hat.2.2 = gama_hat.2
gama_hat.2.2[2,] = 0 # seuillage
gama_hat.2.2[,2] = 0 # seuillage

gama_hat.3.1 = gama_hat.3 # pas de seuillage 

gama_hat.3.2 = gama_hat.3 
gama_hat.3.2[3,] = 0 # seuillage
gama_hat.3.2[,3] = 0 # seuillage

# stockage des lrt bootstrap 

lrt_B.1 = vector(mode = 'numeric', B)
lrt_B.2.1 = vector(mode = 'numeric', B)
lrt_B.2.2 = vector(mode = 'numeric', B)
lrt_B.3.1 = vector(mode = 'numeric', B)
lrt_B.3.2 = vector(mode = 'numeric', B)


# parallèlisation 
  
ncore = 50
grappe = makeCluster(ncore-1) 
registerDoParallel(grappe)

Bootstrap = foreach(b=1:B, .packages = c("lme4", "stats", "dplyr", "MASS", "varTestnlme", "saemix")) %dopar% {
  pdf(NULL)
  
  gr = (b*300000)
  set.seed(gr)
  
  # On g?n?re selon le processus de g?n?ration de donn?es ?tabli (unrestricted parametric bootstrap)
  db.1 = createDatas3eff_rd(N,betaChap1_n,gama_hat.1,sigma1_n,"db")
  db.2.1 = createDatas3eff_rd(N,betaChap1_n,gama_hat.2.1,sigma1_n,"db")
  db.2.2 = createDatas3eff_rd(N,betaChap1_n,gama_hat.2.2,sigma1_n,"db")
  db.3.1 = createDatas3eff_rd(N,betaChap1_n,gama_hat.3.1,sigma1_n,"db")
  db.3.2 = createDatas3eff_rd(N,betaChap1_n,gama_hat.3.2,sigma1_n,"db")
  

  
  saemix.databootH0.1 <- saemixData(name.data = db.1,header=TRUE,name.group=c("subject"),name.predictors=c("days"),
                                  name.response=c("y"), units=list(x="days",y="kg"), verbose = FALSE)
  saemix.databootH0.2.1 <- saemixData(name.data = db.2.1,header=TRUE,name.group=c("subject"),name.predictors=c("days"),
                                    name.response=c("y"), units=list(x="days",y="kg"), verbose = FALSE)
  saemix.databootH0.2.2 <- saemixData(name.data = db.2.2,header=TRUE,name.group=c("subject"),name.predictors=c("days"),
                                    name.response=c("y"), units=list(x="days",y="kg"), verbose = FALSE)
  saemix.databootH0.3.1 <- saemixData(name.data = db.3.1,header=TRUE,name.group=c("subject"),name.predictors=c("days"),
                                    name.response=c("y"), units=list(x="days",y="kg"), verbose = FALSE)
  saemix.databootH0.3.2 <- saemixData(name.data = db.3.2,header=TRUE,name.group=c("subject"),name.predictors=c("days"),
                                    name.response=c("y"), units=list(x="days",y="kg"), verbose = FALSE)
  

  # On calcule la statistique de test sur l'?chantillon Bootstrap
  
  m0_b.1 <- saemix(saemix.modelH0.1,saemix.databootH0.1, saemix.options)
  m1_b.1 <- saemix(saemix.modelH1,saemix.databootH0.1, saemix.options)
  
  
  m0_b.2.1 <- saemix(saemix.modelH0.2,saemix.databootH0.2.1, saemix.options)
  m1_b.2.1 <- saemix(saemix.modelH1,saemix.databootH0.2.1, saemix.options)
  
  m0_b.2.2 <- saemix(saemix.modelH0.2,saemix.databootH0.2.2, saemix.options)
  m1_b.2.2 <- saemix(saemix.modelH1,saemix.databootH0.2.2, saemix.options)
  
  m0_b.3.1 <- saemix(saemix.modelH0.3,saemix.databootH0.3.1, saemix.options)
  m1_b.3.1 <- saemix(saemix.modelH1,saemix.databootH0.3.1, saemix.options)
  
  m0_b.3.2 <- saemix(saemix.modelH0.3,saemix.databootH0.3.2, saemix.options)
  m1_b.3.2 <- saemix(saemix.modelH1,saemix.databootH0.3.2, saemix.options)
  
  
  
  
  l0_b.1 = m0_b.1@results@ll.is
  l1_b.1 = m1_b.1@results@ll.is
  
  l0_b.2.1 = m0_b.2.1@results@ll.is
  l1_b.2.1 = m1_b.2.1@results@ll.is
  
  l0_b.2.2 = m0_b.2.2@results@ll.is
  l1_b.2.2 = m1_b.2.2@results@ll.is
  
  l0_b.3.1 = m0_b.3.1@results@ll.is
  l1_b.3.1 = m1_b.3.1@results@ll.is
  
  l0_b.3.2 = m0_b.3.2@results@ll.is
  l1_b.3.2 = m1_b.3.2@results@ll.is
  
  lrt_B.1 = -2 * (l0_b.1 - l1_b.1)  
  lrt_B.2.1 = -2 * (l0_b.2.1 - l1_b.2.1)
  lrt_B.2.2 = -2 * (l0_b.2.2 - l1_b.2.2)
  lrt_B.3.1 = -2 * (l0_b.3.1 - l1_b.3.1)
  lrt_B.3.2 = -2 * (l0_b.3.2 - l1_b.3.2)
  
  if (lrt_B.1<0){lrt_B.1 = 0}
  if (lrt_B.2.1<0){lrt_B.2.1 = 0}
  if (lrt_B.2.2<0){lrt_B.2.2 = 0}
  if (lrt_B.3.1<0){lrt_B.3.1 = 0}
  if (lrt_B.3.2<0){lrt_B.3.2 = 0}
  
  # suppression des données pour mémoire
  rm(list=c("db.1","db.2.1","db.2.2","db.3.1","db.3.2"))
  
  dev.off()
  res_boot = list("model.1" =lrt_B.1, 
                  "model.2.1" = lrt_B.2.1, 
                  "model.2.2" =  lrt_B.2.2,
                  "model.3.1" =  lrt_B.3.1,
                  "model.3.2" =  lrt_B.3.2)

  res_boot
}

stopCluster(grappe)

for (b in 1:B){
  
  lrt_B.1[b] = Bootstrap[[b]]$model.1
  lrt_B.2.1[b] = Bootstrap[[b]]$model.2.1
  lrt_B.2.2[b] = Bootstrap[[b]]$model.2.2
  lrt_B.3.1[b] = Bootstrap[[b]]$model.3.1
  lrt_B.3.2[b] = Bootstrap[[b]]$model.3.2
  
}


p_b.1 = mean(lrt_B.1 > lrt.1)
p_b.2.1 = mean(lrt_B.2.1 > lrt.2)
p_b.2.2 = mean(lrt_B.2.2 > lrt.2)
p_b.3.1 = mean(lrt_B.3.1 > lrt.3)
p_b.3.2 = mean(lrt_B.3.2 > lrt.3)

p.value = c(p_b.1, p_b.2.1,p_b.2.2, p_b.3.1, p_b.3.2)

table.results = data.frame("p.value" = p.value*100, "sd" = sqrt(p.value*(1-p.value)/B)*100)

list.result = list("mod.tot" = fity0h1,
                   "mod.1" = fity0h0.1, 
                   "mod.2" = fity0h0.2,
                   "mod.3" = fity0h0.3)
print(xtable(table.results,caption = "comparaison model .1, 2.1, 2.2 , 3.1, 3.2, B = 1000", type = "latex"), file = "donnees reelles.tex")

saveRDS(list.result, "donnees_reelles_models.rda")
LRT_B_tot = cbind(lrt_B.1,lrt_B.2.1,lrt_B.2.2,lrt_B.3.1,lrt_B.3.2  )

write.table(LRT_B_tot, file = "LRT_B.csv")
write.table(table.results, file = "p_values.csv")

p_values_ok = c(mean(LRT_B_tot$lrt_B.1>lrt.1)*100,
                mean(LRT_B_tot$lrt_B.2.1>lrt.2)*100,
                mean(LRT_B_tot$lrt_B.2.2>lrt.2)*100,
                mean(LRT_B_tot$lrt_B.3.1>lrt.3)*100,
                mean(LRT_B_tot$lrt_B.3.2>lrt.3)*100)
write.table(p_values_ok, file = "p_values_percentage.csv")
