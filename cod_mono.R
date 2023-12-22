library(AER)
library(VGAM)
library(pscl)
library(MASS)
library(gamlss)
library(xtable)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(plyr)
library(lattice)
library(psych)
library(TeachingDemos)
library(plotrix)
library(e1071)
library(readr)
library(magrittr)
library(pscl)
library(hnp)

# Descritiva
data(bioChemists)

y <- art <- bioChemists[,1]
gen <- bioChemists[,2]
estcvl <- bioChemists[,3]
flhs <- bioChemists[,4]
phd <- bioChemists[,5]
ment.art <- bioChemists[,6]

med_resu <- ddply(bioChemists,.(flhs),
                  summarise,
                  Mín = min(art),
                  Med = median(art),
                  Média = mean(art),
                  Máx = max(art),
                  DP = sqrt(var(art)),
                  Var = var(art),
                  CV = 100*((sqrt(var(art))/mean(art))),
                  Assimetria = skewness(art),
                  Curtose = kurtosi(art),
                  n = length(art))

xtable(med_resu)

# Criar um gráfico de distribuição de frequência
ggplot(data = bioChemists, aes(x = art)) +
  geom_histogram(binwidth = 0.5, fill = "blue") +
  labs(x = "Número de Artigos",
       y = "Frequência")+
  theme_light()

table(ment.art)
table(art)/915*100

table(gen)
table(gen)/915*100
xtable(table(gen))

table(flhs)
table(flhs)/915*100
xtable(table(flhs))

xtable(table(flhs,art))
xtable(table(gen,art))

# Modelos
eta <- art~gen+flhs+ment.art

# poisson
mpo = glm(eta, family = "poisson"(link = "log"))
summary(mpo)
coef.po = summary(mpo)$coefficients[,c(1,2)]

# poisson - Gama
mpg = glm.nb(eta, link = log)
summary(mpg)
coef.pg = summary(mpg)$coefficients[,c(1,2)]


# poisson - Inversa Gaussiana
mpig = gamlss(eta,
              family = PIG(mu.link = "log"),
              trace = FALSE ) # Global Deviance : -2* loglik
summary(mpig)
coef.pig = summary(mpig)[c(1:4),c(1,2)]

# Log - v e r o s s i m i l h a n a PLG
L = function (y,x,par){
  k <- length (par)
  theta <- par[1]
  beta <- as.matrix(par[2:k])
  ee = function(aux){
    sol = TRUE
    for (i in 1:length(aux)) {
      if (sol & aux[i]) {
      } else {
        sol = FALSE
      }
    }
    return(sol)
  }
  if (theta > 0) {
    al =(exp(x%*%beta)*theta*(theta+1)-1)/(theta+1)
    if (ee(al>0) & max(al)<100) {
      A = sum(log(gamma(y+al))-(log(factorial(y)) + log(gamma(1+al)))+
                (al+1)*log(theta)-(y+al+2)*log(theta+1)+log((al*(theta+1)+
                                                               y+al)))
      return(A)
    } else {
      return ( - Inf )
    }
  } else {
    return ( - Inf )
  }
}

# - logverossimilhanca
Lest <- function (par) {
  return (-L (y,x,par))
}


############################### poisson - Lindley Generalizada
x = model.matrix(mpg)

# ---------- poisson
i = 0.7
j = 1
n = length(y)
p = ncol(x)

matplg.po = matrix(0,ncol = 9, length(seq(0.7 ,1.8, by =0.1))-1)
colnames(matplg.po) = c("Valor inicial", "Theta Estimado",
                        "Beta_0", "Beta_1", "Beta_2", "Beta_3",
                        "AIC", "BIC", "LogV")

while (i <= 1.8) {
  e3 = optim(c(i,coef.po[,1]),Lest,method = "SANN")
  predito = exp(x%*%e3$par[2:5])
  AICplg = -2*L(y,x,e3$par)+2*p
  BICplg = -2*L(y,x,e3$par)+log(n)*p
  Logvero = Lest(e3$par)
  
  matplg.po[j,1] = i
  matplg.po[j,2:6] = unlist(e3$par)
  matplg.po[j,7] = AICplg
  matplg.po[j,8] = BICplg
  matplg.po[j,9] = Logvero
  
  j = j+1
  i = i+0.1
}
matplg.po


# ---------- poisson - Gama
i = 0.7
j = 1
n = length(y)
p = ncol(x)

matplg.pg = matrix(0 , ncol = 9 , length(seq(0.7 ,1.8 , by =0.1) ) -1)
colnames(matplg.pg) = c("Valor inicial", "Theta Estimado",
                        "Beta_0", "Beta_1", "Beta_2", "Beta_3",
                        "AIC", "BIC", "LogV")
while (i <= 1.8) {
  e3 = optim(c(i,coef.pg[,1]),Lest,method = "SANN")
  predito = exp(x%*%e3$par[2:5])
  AICplg = -2*L(y,x,e3$par)+2*p
  BICplg = -2*L(y,x,e3$par)+log(n)*p
  Logvero = Lest(e3$par)
  
  matplg.pg[j,1] = i
  matplg.pg[j,2:6] = unlist(e3$par)
  matplg.pg[j,7] = AICplg
  matplg.pg[j,8] = BICplg
  matplg.pg[j,9] = Logvero
  
  j = j+1
  i = i+0.1
}
matplg.pg

# ---------- poisson - Inversa Guassiana
i = 0.7
j = 1
n = length(y)
p = ncol(x)

matplg.pig = matrix(0 , ncol = 9 , length(seq(0.7, 1.8 , by=0.1))-1)
colnames(matplg.pig) = c("Valor inicial", "Theta Estimado",
                         "Beta_0", "Beta_1", "Beta_2", "Beta_3",
                         "AIC", "BIC", "LogV")

while ( i <= 1.8) {
  e3 = optim(c(i,coef.pig[,1]),Lest,method = "SANN")
  predito = exp(x%*%e3$par[2:5])
  AICplg = -2*L(y,x,e3$par)+2*p
  BICplg = -2*L(y,x,e3$par)+log(n)*p
  Logvero = Lest(e3$par)
  
  matplg.pig[j,1] = i
  matplg.pig[j,2:6] = unlist(e3$par)
  matplg.pig[j,7] = AICplg
  matplg.pig[j,8] = BICplg
  matplg.pig[j,9] = Logvero
  
  j = j+1
  i = i+0.1
}
matplg.pig

# ---------- Modelo poisson - Lindley Generalizada
# Selecionando com melhor AIC / BIC
plg.po = 
  matplg.po[which.min(apply(cbind(matplg.po[,7],matplg.po[,8]),1,sum)),]

plg.pg = 
  matplg.pg[which.min(apply(cbind(matplg.pg[,7],matplg.pg[,8]),1,sum)),]

plg.pig = 
  matplg.pig[which.min(apply(cbind(matplg.pig[,7],matplg.pig[,8]),1,sum)),]

m.plg = 
  matrix(c(plg.po,plg.pg,plg.pig),nrow = 3, byrow = TRUE)
colnames(m.plg) = c("Valor inicial", "Theta Estimado",
                    "(Intercept)","genWomen","flhs","ment.art",
                    "AIC", "BIC", "LogV")

final.plg = m.plg[which.min(apply(cbind(m.plg[,7],m.plg[,8]),1,sum)),]
coef.plg = c(unlist(final.plg[2:6])) # Estimativas
predito.plg = exp(x%*%coef.plg[2:5]) #valore estimados

dseg = function(y,x,par) {
  k = length(par)
  theta = par[1]
  beta = par[2:k]
  al = (exp(x%*%beta)*theta*(theta+1)-1)/(theta+1)
  vdc = exp(x%*%beta)
  res = sum(theta*vdc*(vdc*theta*trigamma(y + al) +
                         digamma(y + al ) + vdc * theta * trigamma(1+ al ) + 
                         digamma(1+ al ) - log( theta /( theta +1) ) +
                         ((theta +1) *( theta +2) *( y *( theta +1) -
                         (theta +2)))/((vdc*theta*(theta +1)-1)*(theta +2) 
                                           + y*(theta +1))^2)) *( t(x)%*%x)
  
  return(res)
}

## sabemos que hat ( beta ) ~ N ( beta , invertida )
n = 915
inform = (1/n)*dseg(y,x,coef.plg);inform # matriz de informacao
invertida = ginv(inform); invertida # inverso da matrz de informacao
epb0 = sqrt(invertida[1 ,1]); epb0 # erro padrao beta0
epb1 = sqrt(invertida[2 ,2]); epb1 # erro padrao beta1
epb2 = sqrt(invertida[3 ,3]); epb2 # erro padrao beta2
epb3 = sqrt(invertida[4 ,4]); epb3 # erro padrao beta3

t0cal = coef.plg[2]/ epb0 ; t0cal # estatistica t calculada
t1cal = coef.plg[3]/ epb1 ; t1cal
t2cal = coef.plg[4]/ epb2 ; t2cal
t3cal = coef.plg[5]/ epb3 ; t3cal

p0 = pt(abs(t0cal),n-1, lower.tail = F ); p0 # valor - p associado
p1 = pt(abs( t1cal ),n-1 , lower.tail = F ) ; p1
p2 = pt(abs( t2cal ),n-1 , lower.tail = F ) ; p2
p3 = pt(abs( t3cal ),n-1 , lower.tail = F ) ; p3

# Modelo Hurdle poisson, com zeros: Binomial
mhpb <- hurdle(eta,
               dist = "poisson",
               zero.dist = "binomial")
summary(mhpb)
coef.hpb.c = summary(mhpb)$coefficients$count[,c(1,2)]
coef.hpb.z = summary(mhpb)$coefficients$zero[,c(1,2)]


# Modelo Hurdle poisson, com zeros: poisson
mhpp <- hurdle(eta,
               dist = "poisson",
               zero.dist = "poisson")
summary(mhpp)
coef.hpp.c = summary(mhpp)$coefficients$count[,c(1,2)]
coef.hpp.z = summary(mhpp)$coefficients$zero[,c(1,2)]

# Modelo Hurdle poisson, com zeros: Binomial Negativo
mhpbn <- hurdle(eta,
                dist = "poisson",
                zero.dist = "negbin")
summary(mhpbn)
coef.hpbn.c = summary(mhpbn)$coefficients$count[,c(1,2)]
coef.hpbn.z = summary(mhpbn)$coefficients$zero[,c(1,2)]

# Modelo Hurdle Binomial Negativo, com zeros: Binomial
mhbnb <- hurdle(eta,
                dist = "negbin",
                zero.dist = "binomial")
summary(mhbnb)
coef.hbnb.c = summary(mhbnb)$coefficients$count[,c(1,2)]
coef.hbnb.z = summary(mhbnb)$coefficients$zero[,c(1,2)]


# Modelo Hurdle Binomial Negativo, com zeros: poisson
mhbnp <- hurdle(eta,
                dist = "negbin",
                zero.dist = "poisson")
summary(mhbnp)
coef.hbnp.c = summary(mhbnp)$coefficients$count[,c(1,2)]
coef.hbnp.z = summary(mhbnp)$coefficients$zero[,c(1,2)]

# Modelo Hurdle Binomial Negativo, com zeros: Binomial Negativa
mhbnbn <- hurdle(eta,
                 dist = "negbin",
                 zero.dist = "negbin")
summary(mhbnbn)
coef.hbnbn.c = summary(mhbnbn)$coefficients$count[,c(1,2)]
coef.hbnbn.z = summary(mhbnbn)$coefficients$zero[,c(1,2)]

#### Estimativas e erros dos parametros
coef.plg = c(unlist(final.plg[3:6]))
coef.plg <- data.frame(coef.plg)
coef.plg <- coef.plg %>%
  mutate(Std.Error = c(epb0,epb1,epb2,epb3))


est_hi <- round(rbind(coef.po[,1],coef.pg[,1]
                      ,coef.pig[,1],coef.plg[,1]),4)
xtable(est_hi, digits = 4)
err_hi <- rbind(round(rbind(coef.po[,2],
                            coef.pg[,2],coef.pig[,2],coef.plg[,2]),4))
xtable(err_hi, digits = 4)


est_huc <- round(rbind(coef.hpb.c[c(1:4),1],
                       coef.hpp.c[c(1:4),1],coef.hpbn.c[c(1:4),1],
                       coef.hbnb.c[c(1:4),1],
                       coef.hbnp.c[c(1:4),1],coef.hbnbn.c[c(1:4),1]),4)
xtable(est_huc, digits = 4)
err_huc <- rbind(round(rbind(coef.hpb.c[c(1:4),2],
                             coef.hpp.c[c(1:4),2],coef.hpbn.c[c(1:4),2],
                             coef.hbnb.c[c(1:4),2],
                             coef.hbnp.c[c(1:4),2],coef.hbnbn.c[c(1:4),2]),4))
xtable(err_huc, digits = 4)

est_huz <- round(rbind(coef.hpb.z[c(1:4),1],
                       coef.hpp.z[c(1:4),1],coef.hpbn.z[c(1:4),1],
                       coef.hbnb.z[c(1:4),1],
                       coef.hbnp.z[c(1:4),1],coef.hbnbn.z[c(1:4),1]),4)
xtable(est_huz, digits = 4)
err_huz <- rbind(round(rbind(coef.hpb.z[c(1:4),2],
                             coef.hpp.z[c(1:4),2],coef.hpbn.z[c(1:4),2],
                             coef.hbnb.z[c(1:4),2],
                             coef.hbnp.z[c(1:4),2],coef.hbnbn.z[c(1:4),2]),4))
xtable(err_huz, digits = 4)

#### AIC e BIC de cada modelo

abmhivsmhu <- rbind(cbind(rbind(AIC(mpo),AIC(mpg),AIC(mpig),
                                final.plg[7]),
                          rbind(BIC(mpo),BIC(mpg),BIC(mpig),
                                final.plg[8])),
                    cbind(rbind(AIC(mhpb),AIC(mhpp),AIC(mhpbn),
                                AIC(mhbnb),AIC(mhbnp),AIC(mhbnbn)),
                          rbind(BIC(mhpb),BIC(mhpp),BIC(mhpbn),
                                BIC(mhbnb),BIC(mhbnp),BIC(mhbnbn))))
dim(abmhivsmhu)
dimnames(abmhivsmhu) <- list(c("MP","MPG","MPIG","MPLG",
                               "MHP-B","MHP-P","MHP-BN",
                               "MHBN-B","MHBN-P","MHBN-BN"),
                             c("AIC","BIC"));abmhivsmhu
xtable(abmhivsmhu, digits = 4)

# ---------- Poisson
predito.po = exp(x %*% coef.po[,1])
var.po = exp(x %*% coef.po[,1])

# ---------- Poisson - Gama
predito.pg = exp(x %*% coef.pg[,1])
var.pg = ((exp(x %*% coef.pg[,1])^2) / mpg$theta) + exp(x %*% coef.pg[,1])

# ---------- Poisson - Inversa Gaussiana
predito.pig = exp(x %*% coef.pig[,1])
var.pig = ((exp(x %*% coef.pig[,1])^3) / mpig$sigma.coefficients) +
  exp(x %*% coef.pig[,1])

# ---------- Poisson Lindley Generalizada
predito.plg
alpha = (exp(x %*% coef.plg[2:5]) * coef.plg[1] *
           (coef.plg[1]+1) - 1) / (coef.plg[1]+1)
var.plg = exp(x %*% coef.plg[2:5]) + 
  (((alpha + 1) * (alpha * (coef.plg[1]+1) + 2)) / (coef.plg[1]^2 *
                     (coef.plg[1]+1))) - exp(x %*% coef.plg[2:5])^2

set.seed(475050)
# ---------- Poisson
Acumu.po = ppois(y - 1, predito.po)
Prob.po = dpois(y, predito.po)
quantilico.po = qnorm(Acumu.po + runif(915) * Prob.po)
wp(resid = quantilico.po,main = "Modelo Hurdle Binominal Negativa-Binominal")

# ---------- Poisson - Gama
Acumu.pg = pnbinom(y - 1, size = mpg$theta,
                   prob = mpg$theta / (mpg$theta + predito.pg))
Prob.pg = dnbinom(y, size = mpg$theta,
                  prob = mpg$theta / (mpg$theta + predito.pg))
quantilico.pg = qnorm(Acumu.pg + runif(915) * Prob.pg)
wp(resid = quantilico.pg)

# ---------- Poisson - Inversa Gaussiana
wp(gamlss(eta, data = bioChemists, family = PIG(mu.link = "log"),
          trace = FALSE))

# ---------- Poisson - Lindley Generalizada
p = coef.plg[1] / (1 + coef.plg[1])
Acumu.plg = (coef.plg[1] / (1 + coef.plg[1])) *
  pnbinom(y - 1, alpha, p) + (1 / (1 + coef.plg[1])) *
  pnbinom(y - 1, alpha + 1, p)
Prob.plg = (coef.plg[1] / (1 + coef.plg[1])) *
  dnbinom(y, alpha, p) + (1 / (1 + coef.plg[1])) * dnbinom(y, alpha + 1, p)
quantilico.plg = qnorm(Acumu.plg + runif(915) * Prob.plg)
wp(resid = quantilico.plg)

set.seed(475050)
par(mfrow = c(2,2))
hnp(mpo, ylab = "Resíduos", xlab = "Quantis teôricos",
    main = "Modelo Poisson",resid.type = "pearson",
    pch = 20,how.many.out=TRUE, paint.out = TRUE)

hnp(mhpb, ylab = "Resíduos", xlab = "Quantis teôricos",
    main = "Modelo Hurdle Poisson-Binomial",
    resid.type = "pearson",pch = 20,how.many.out=TRUE, paint.out = TRUE)

hnp(mpg, ylab = "Resíduos", xlab = "Quantis teôricos",
    main = "Modelo Poisson-Gama",
    resid.type = "pearson",pch = 20,how.many.out=TRUE, paint.out = TRUE)

hnp(mhbnb, ylab = "Resíduos", xlab = "Quantis teôricos",
    main = "Modelo Hurdle Binomial Negativa-Binomial",
    resid.type = "pearson",pch = 20,how.many.out=TRUE,
    paint.out = TRUE )

