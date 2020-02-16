### R code from vignette source 'simposio_salvador.Rnw'
#alterado

# Editado por Andre Chalom 

options(width=60, continue=" ")
# Numero de simulacoes por ponto do espaco de param:
NSIM = 20
# Tempo final da simulacao
TEMPO = 5
SIZE=400

###################################################
### Lê o arquivo com taxas basais e de morte
###################################################
mccoy <- read.csv2("mccoy_mortality_rates.csv", dec=".")
mccoy.m <- mccoy[mccoy$group=="Mammal",]
mccoy.m <- aggregate(mccoy.m[,c(3,5)], by=list(specie=mccoy.m$species), mean)
## species between 10g and 1.5 ton body mass
mccoy.m <- mccoy.m[mccoy.m$dry.mass.g>10&mccoy.m$dry.mass.g<1.5e6,]
## exluding outliers of mortality rate
mccoy.m <- mccoy.m[mccoy.m$mortality.year<=1,]
## Mean and sd of log values
(taxa.morte.mu <- mean(log(mccoy.m$mortality.year)))
(taxa.morte.sd <- sd(log(mccoy.m$mortality.year)))


###################################################
### code chunk number 3: calc-basal
###################################################
## Calculo da taxa intrinseca maxima anual com equacao de Blueweiss et al
mccoy.m$rmax.year <- (10^(-0.26*log10(mccoy.m$dry.mass.g)+log10(0.025)))*365
## Taxa basal como proporcao da taxa de mortalidade
mccoy.m$basal.prop <- (mccoy.m$rmax.year-mccoy.m$mortality.year)/mccoy.m$mortality.year
## media e desvio-padrão dos logs das proporcoes
basal.p.mu <- mean(log(mccoy.m$basal.prop[mccoy.m$basal.prop>=1]))
basal.p.sd <- sd(log(mccoy.m$basal.prop[mccoy.m$basal.prop>=1]))


###################################################
### code chunk number 4: calc-raio
###################################################
## Densidade por km2, equacao Damuth:
mccoy.m$dens.km2 <- 10^(-0.75*log10(mccoy.m$dry.mass.g)+4.23)
## Raio da area individual em metros
mccoy.m$raio.m <- sqrt(1/(mccoy.m$dens.km2*pi))*1000
## Media e sd do log do raio
raio.mu <- mean(log(mccoy.m$raio.m))
raio.sd <- sd(log(mccoy.m$raio.m))


###################################################
### code chunk number 5: create-dynamic-library-for-C-routines
###################################################
Sys.setenv("PKG_CPPFLAGS" = "-fopenmp")
system ("R CMD SHLIB TWoLife.cpp")
dyn.load("TWoLife.so")
TWoLife <- function (tamanho, raio, N, AngVis, passo, move, densi.max, 
                       taxa.basal, taxa.morte, fragmentacao, numero.manchas,
                       tempo) {
	saida.C <- .C("TWoLife", as.double(tamanho), as.double(raio),
				  as.integer(N),as.double(AngVis), as.double(passo),
				  as.double(move), as.double(densi.max), 
				  as.double(taxa.basal), as.double(taxa.morte),
				  as.integer(fragmentacao), as.double(numero.manchas),
				  as.double(tempo), as.integer(0))
	return(saida.C[[13]])
}

###################################################
### code chunk number 7: R-function-for-simulations
###################################################
ProbExt2 <- function (Nsim, N, AngVis, passo.rel, move.rel, ln.taxa.morte, 
                     ln.basal.prop, ln.raio, tempo) {
  taxa.morte <- exp(ln.taxa.morte)
  basal.prop <- exp(ln.basal.prop)
  raio <- exp(ln.raio)
  densi.max <- 1.1/(pi*raio^2)# 1,1 individuo por area de vizinhanca
  tamanho <- sqrt(20*pi*raio^2)
  taxa.basal <- basal.prop*taxa.morte
  passo <- passo.rel*raio
  move <- move.rel*(taxa.basal+taxa.morte)
  ## tempo sera o esperado para a ocorrencia de tempo.rel mortes per capita
  #p.morte <- taxa.morte/(taxa.morte+taxa.basal+move)
  #tempo <- (tempo.rel/p.morte)/(taxa.morte+taxa.basal+move)
  vivo = c()
  for (i in 1:Nsim) {
    vivo[i] <- VivoMorto(tamanho, raio, N, AngVis, passo, move, 
                         densi.max, taxa.basal, taxa.morte,
                         fragmentacao=0, numero.manchas=0, tempo)
  }
  return (sum(vivo==0)/Nsim)
}
modelProbExt2 <- function(dados) {
  return(mapply(ProbExt2, Nsim=NSIM, N=dados[,1], AngVis=dados[,2],
                passo.rel=dados[,3], move.rel=dados[,4], ln.taxa.morte=dados[,5],
                ln.basal.prop=dados[,6], ln.raio=dados[,7], tempo=TEMPO))
}

###################################################
### code chunk number 9: definicao-espaco-parametros
###################################################
factors2 <- c("N", "AngVis", "passo.rel", "move.rel", "ln.taxa.morte", "ln.basal.prop", "ln.raio")
q2 <- rep(c("qunif","qnorm"),c(4,3))
q.arg2 <- list ( 
  list(min=1, max=20),
  list(min=1, max=360),
  list(min=0.1, max=3),
  list(min=0.1, max=20),
  list(mean=taxa.morte.mu, sd=taxa.morte.sd),
  list(mean=basal.p.mu, sd=basal.p.sd),
  list(mean=raio.mu, sd=raio.sd)
  )

###################################################
### code chunk number 16: pse-sem-correlacoes (eval = FALSE)
###################################################
library(pse) #carrega as funcoes de pse
LHSnull <- LHS(model=NULL, factors2, SIZE, q2, q.arg2)
#save(LHSnull, file="LHSnull.Rdata");
### Rodar no cluster:
# load("LHSnull.Rdata")
system.time(res <- modelProbExt2(LHSnull$data))
# save(res, file="modelres.Rdata");
#load("modelres.Rdata");
fullLHS <- tell(LHSnull, res);
save(fullLHS, file="vivomorto.Rdata")
corPlot(fullLHS)
plotprcc(fullLHS) # Mostra N e ln.basal.prop negativos como as unicas vars significativas

###################################################
### code chunk number 17: lhs-sem-cor-dataframe
###################################################
lhs.sc <- cbind(fullLHS$data,fullLHS$res) 
names(lhs.sc)[8] <- "prop.ext"
lhs.sc$taxa.morte <- exp(lhs.sc$ln.taxa.morte)
lhs.sc$taxa.nasc <- exp(lhs.sc$ln.basal.prop)*lhs.sc$taxa.morte
lhs.sc$raio <- exp(lhs.sc$ln.raio)
lhs.sc$D <- with(lhs.sc, sqrt(passo.rel*raio)/
  (2/(move.rel*(taxa.morte+taxa.nasc)))) #coef difusao
lhs.sc$n.ext <- lhs.sc$prop.ext*SIZE
lhs.sc$n.sobr <- (1-lhs.sc$prop.ext)*SIZE

###################################################
### code chunk number 18: glm-partial-plots-lhs-sem-cor
###################################################
library(car)
m1 <- glm(cbind(n.ext,n.sobr)~ N + AngVis + D +  taxa.morte + taxa.nasc + raio, data=lhs.sc, family=binomial)

plot(pcc(lhs.sc[,c(1,2,12,9,10,11)],lhs.sc[,13], rank=TRUE, nboot=100)) # mais coerente, mas ainda mostra AngVis e D como NS
###################################################
### code chunk number 19: partial-plots-lhs-sem-cor
###################################################
lhs.sc.av <- avPlots(m1)
par(mfrow=c(2,3))
plot(lhs.sc.av$N[,2]~I(lhs.sc.av$N[,1]/500))
plot(lhs.sc.av$AngVis[,2]~I(lhs.sc.av$AngVis[,1]/500), ylim=c(0,700))
plot(lhs.sc.av$D[,2]~I(lhs.sc.av$D[,1]/500), ylim=c(0,700))
plot(lhs.sc.av$taxa.morte[,2]~I(lhs.sc.av$taxa.morte[,1]/500), ylim=c(0,700))
plot(lhs.sc.av$taxa.nasc[,2]~I(lhs.sc.av$taxa.nasc[,1]/500), ylim=c(0,1500))
plot(lhs.sc.av$raio[,2]~I(lhs.sc.av$raio[,1]/500), ylim=c(0,700))


###################################################
### code chunk number 20: lhs-sem-cor-model-selection
###################################################
library(MASS)
m1.step <- stepAIC(glm(cbind(n.ext,n.sobr)~ N + AngVis + D +  taxa.morte + taxa.nasc + raio +
                       N:AngVis + N:D + N:taxa.morte + N:taxa.nasc + N:raio +
                       AngVis:D + AngVis:taxa.morte + AngVis:taxa.nasc + AngVis:raio + 
                       D:taxa.morte + D:taxa.nasc + D:raio +
                       taxa.morte:taxa.nasc + taxa.morte:raio +
                       taxa.nasc:raio, 
                       data=lhs.sc, family=binomial), direction="both")
summary(m1.step)
m1.step.b <- update(m1.step, . ~. -AngVis)
AIC(m1.step,m1.step.b)


###################################################
### code chunk number 21: correlation-matrix
###################################################
pse.cor.mat <- matrix(0, ncol=length(factors2), nrow=length(factors2),
                      dimnames=list(factors2, factors2))
## diagonais igual a um
diag(pse.cor.mat) <- 1
## substituindo correlacoes por valores estimados de Mccoy et al
pse.cor.mat[5:7,5:7] <- cor(log(mccoy.m[,c(3,5,7)]))


###################################################
### code chunk number 22: pse-com-correlacoes (eval = FALSE)
###################################################
## LHS.com.cor.200 <- LHS(modelProbExt2, factors2, 200, q2, q.arg2, COR=pse.cor.mat)


###################################################
### code chunk number 23: lhs-com-cor-dataframe
###################################################
lhs.cc <- cbind(fullLHS$data,fullLHS$res) 
names(lhs.sc)[8] <-"prop.ext"
lhs.cc$taxa.morte <- exp(lhs.cc$ln.taxa.morte)
lhs.cc$taxa.nasc <- exp(lhs.cc$ln.basal.prop)*lhs.cc$taxa.morte
lhs.cc$raio <- exp(lhs.cc$ln.raio)
lhs.cc$D <- with(lhs.cc, sqrt(passo.rel*raio)/
  (2/(move.rel*(taxa.morte+taxa.nasc)))) #coef difusao
lhs.cc$n.ext <- lhs.cc$prop.ext*500
lhs.cc$n.sobr <- (1-lhs.cc$prop.ext)*500


###################################################
### code chunk number 24: glm-partial-plots-lsh-com-cor
###################################################
library(car)
m1.cc <- glm(cbind(n.ext,n.sobr)~ N + AngVis + D +  taxa.morte + taxa.nasc + raio, data=lhs.cc, family=binomial)


###################################################
### code chunk number 25: partial-plots-lsh-com-cor
###################################################
lhs.cc.av <- avPlots(m1.cc)
par(mfrow=c(2,3))
plot(lhs.cc.av$N[,2]~I(lhs.cc.av$N[,1]/500))
plot(lhs.cc.av$AngVis[,2]~I(lhs.cc.av$AngVis[,1]/500), ylim=c(0,700))
plot(lhs.cc.av$D[,2]~I(lhs.cc.av$D[,1]/500), ylim=c(0,700))
plot(lhs.cc.av$taxa.morte[,2]~I(lhs.cc.av$taxa.morte[,1]/500), ylim=c(0,700))
plot(lhs.cc.av$taxa.nasc[,2]~I(lhs.cc.av$taxa.nasc[,1]/500), ylim=c(0,1500))
plot(lhs.cc.av$raio[,2]~I(lhs.cc.av$raio[,1]/500), ylim=c(0,700))


###################################################
### code chunk number 26: lsh-com-cor-model-selection
###################################################
library(MASS)
m1.cc.step <- stepAIC(glm(cbind(n.ext,n.sobr)~ N + AngVis + D +  taxa.morte + taxa.nasc + raio +
                          N:AngVis + N:D + N:taxa.morte + N:taxa.nasc + N:raio +
                          AngVis:D + AngVis:taxa.morte + AngVis:taxa.nasc + AngVis:raio + 
                          D:taxa.morte + D:taxa.nasc + D:raio +
                          taxa.morte:taxa.nasc + taxa.morte:raio +
                          taxa.nasc:raio,
                          data=lhs.cc, family=binomial), direction="both")
summary(m1.cc.step)
coef(m1.cc.step)
## Outras tentativas
m1.cc.step.b <- update(m1.cc.step, . ~. -taxa.nasc - D:taxa.nasc - N:taxa.nasc - taxa.morte:taxa.nasc - taxa.nasc:raio)
m2 <- glm(cbind(n.ext,n.sobr)~ N + D +  AngVis + taxa.morte + taxa.nasc +
            N:D + N:taxa.morte + N:taxa.nasc, data=lhs.cc, family=binomial)
summary(m1.cc.step.b)
AIC(m1.cc.step,m1.cc.step.b, m2)


###################################################
### code chunk number 27: data-frames-for-predicted-values
###################################################
plot.lhs <- function(model= m1.cc.step, y="D", cond="taxa.morte" , N=2, df=lhs.cc[,c(1:2,9:12)],...){
  pred.df <- as.data.frame(matrix(apply(df, 2 ,median),nrow=200, ncol=6, byrow=TRUE,
                                  dimnames=list(NULL,names(df))))
  pred.df$N <- N
  i <- which(names(pred.df)==y)
  j <- which(names(pred.df)==cond)
  pred.df[,j] <- rep(range(df[,j]),each=100)
  pred.df[,i] <- rep(seq(min(df[,i]),max(df[,i]),length=100),2)
  m1.cc.fit <- predict(model,newdata=pred.df, type="response")
  plot(m1.cc.fit[1:100]~pred.df[1:100,i], type="l", ylim=range(m1.cc.fit), col="red", ...)
  lines(m1.cc.fit[101:200]~pred.df[101:200,i], type="l", col="blue")
}


###################################################
### code chunk number 28: plot-fitted-values
###################################################
par(mfrow=c(2,2))
plot.lhs(y="taxa.nasc", cond="N")
plot.lhs(y="taxa.morte", cond="N")
plot.lhs(y="D", cond="N")
plot.lhs(y="raio", cond="N")
par(mfrow=c(1,1))


###################################################
### code chunk number 29: plot-fitted-values
###################################################
par(mfrow=c(1,2))
plot.lhs(y="taxa.nasc", cond="taxa.morte", ylab="P extinção", xlab="Taxa nascimentos")
plot.lhs(y="D", cond="taxa.morte", ylab="", xlab="Coeficiente de difusão")
legend("topright", c(paste("Mortalidade=",round(min(lhs.cc$taxa.morte),2)),paste("Mortalidade=",round(max(lhs.cc$taxa.morte),2))),
       lty=1, col=c("red", "blue"),bty="n", cex=1.25)
par(mfrow=c(1,1))


