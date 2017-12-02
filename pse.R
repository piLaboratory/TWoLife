library("pse")

#Cálculo dos parâmetros
L <- 0.03*512 #Em km
mulamb0 <- 1.26
sigmalamb0 <- 0.44
mut <- 1.89
sigmat <- 0.21
muN <- 3.24
sigmaN <- 0.86
mud <- 1.31
sigmad <- 1.92
rholamb0t <- -0.11
rholamb0N <- 0
rholamb0d <- -0.18
rhotN <- 0.3
rhotd <- -0.254
rhoNd <- 0.79

muu <- log(2)/2
sigmau <- log(2)/(3.64*sqrt(2))
muK <- muu + muN
sigmaK <- sqrt(sigmau^2 + sigmaN^2)
mumu0 <- -mut
sigmamu0 <- sigmat
alpha <- 0.4957
beta <- exp(-0.6796)
muR <- log(L) + log(beta) - alpha*muK
sigmaR <- alpha*sigmaK
muv <- log(5)/2
sigmav <- log(5)/(3.64*sqrt(2))
mus <- muv + muR
sigmas <- sqrt(sigmav^2 + sigmaR^2)
mudelta <- log(4/pi) - 2*log(L) - 2*log(beta) + 2*mud + 2*alpha*muu + 2*alpha*muN - 2*muv - mut
sigmadelta <- sqrt( 4*sigmad^2 + 4*alpha^2*sigmau^2 + 4*alpha^2*sigmaN^2 + 4*sigmav^2 + sigmat^2 + 8*alpha*rhoNd*sigmad*sigmaN - 4*rhotd*sigmad*sigmat - 4*alpha*rhotN*sigmaN*sigmat)
muw <- log(50)/2
sigmaw <- log(50)/(3.64*sqrt(2))

rholamb0R <- -alpha*rholamb0N*sigmaN/sigmaR
rhomu0R <- alpha*rhotN**sigmaN*sigmat/(sigmaR*sigmamu0)

mub <- muw + mulamb0 + 2*muR + log(pi/50)
sigmab <- sqrt( sigmaw^2 + sigmalamb0^2 + 4*sigmaR^2 + 4*rholamb0R*sigmalamb0*sigmaR )
mum <- muw + mumu0 + 2*muR + log(pi/50)
sigmam <- sqrt( sigmaw^2 + sigmamu0^2 + 4*sigmaR^2 + 4*rhomu0R*sigmamu0*sigmaR )

rhosR <- sigmaR/sigmas
rhodeltaR <- 2*alpha*rhoNd*sigmaN*sigmad/(sigmaR*sigmadelta) - 2*sigmaR/sigmadelta + alpha*rhotN*sigmaN*sigmat/(sigmaR*sigmadelta)
rhobR <- rholamb0R*sigmalamb0/sigmab + 2*sigmaR/sigmab
rhomR <- rhomu0R*sigmamu0/sigmam + 2*sigmaR/sigmam
rhodeltas <- rhodeltaR*sigmaR/sigmas
rholamb0s <- rholamb0R*sigmaR/sigmas
rhomu0s <- rhomu0R*sigmaR/sigmas
rhobs <- rhobR*sigmaR/sigmas
rhoms <- rhomR*sigmaR/sigmas
rhodeltalamb0 <- 2*rholamb0d*sigmad/sigmadelta - 2*rholamb0s*sigmas/sigmadelta - rholamb0t*sigmat/sigmadelta
rhodeltamu0 <- -2*rhotd*sigmad*sigmat/(sigmadelta*sigmamu0) - 2*rhomu0s*sigmas/sigmadelta + sigmat^2/(sigmadelta*sigmamu0)
rhodeltab <- 2*rholamb0d*sigmad*sigmalamb0/(sigmadelta*sigmab) - 4*alpha*rhoNd*sigmad*sigmaN/(sigmadelta*sigmab) - 2*rholamb0s*sigmas*sigmalamb0/(sigmadelta*sigmab) - 4*rhosR*sigmas*sigmaR/(sigmadelta*sigmab) - rholamb0t*sigmat*sigmalamb0/(sigmadelta*sigmab) + 2*rhomu0R*sigmaR*sigmamu0/(sigmadelta*sigmab)
rhodeltam <- -2*rhotd*sigmad*sigmat/(sigmadelta*sigmam) - 4*alpha*rhoNd*sigmad*sigmaN/(sigmadelta*sigmam) - 2*rhomu0s*sigmas*sigmamu0/(sigmadelta*sigmam) - 4*rhosR*sigmas*sigmaR/(sigmadelta*sigmam) + sigmat^2/(sigmadelta*sigmam) + 2*rhomu0R*sigmaR*sigmamu0/(sigmadelta*sigmab)
rhomu0lamb0 <- -rholamb0t*sigmat/sigmamu0
rholamb0b <- sigmalamb0/sigmab + 2*rholamb0R*sigmaR/sigmab
rholamb0m <- -rholamb0t*sigmat/sigmam + 2*rholamb0R*sigmaR/sigmam
rhomu0b <- rhomu0lamb0*sigmalamb0/sigmab + 2*rhomu0R*sigmaR/sigmab
rhomu0m <- sigmamu0/sigmab + 2*rhomu0R*sigmaR/sigmab
rhobm <- rholamb0m*sigmalamb0/sigmab + 2*rhomR*sigmaR/sigmab

factors <- c("R", "s", "delta", "lamb0", "mu0", "b", "d", "dm")
q <- c("qnorm", "qnorm", "qnorm", "qnorm", "qnorm", "qnorm", "qnorm", "qunif")
q.arg <- list( list(mean=muR, sd=sigmaR ), list(mean=mus, sd=sigmas ), list(mean=mudelta, sd=sigmadelta ), list(mean=mulamb0, sd=sigmalamb0 ), list(mean=mumu0, sd=sigmamu0 ), list(mean=mub, sd=sigmab ), list(mean=mum, sd=sigmam ), list(min=1, max=10) )

COR <- matrix(c(1, rhosR, rhodeltaR, rholamb0R, rhomu0R, rhobR, rhomR, 0,
								rhosR, 1, rhodeltas, rholamb0s, rhomu0s, rhobs, rhoms, 0,
								rhodeltaR, rhodeltas, 1, rhodeltalamb0, rhodeltamu0, rhodeltab, rhodeltam, 0,
								rholamb0R, rholamb0s, rhodeltalamb0, 1, rhomu0lamb0, rholamb0b, rholamb0m, 0,
								rhomu0R, rhomu0s, rhodeltamu0, rhomu0lamb0, 1, rhomu0b, rhomu0m, 0,
								rhobR, rhobs, rhodeltab, rholamb0b, rhomu0b, 1, rhobm, 0,
								rhomR, rhoms, rhodeltam, rholamb0m, rhomu0m, rhobm, 1, 0,
								0, 0, 0, 0, 0, 0, 0, 1), nrow = 8, ncol = 8)

opts <- list(COR)
myLHS <- LHS(model=NULL, factors, N=10, q=q, q.arg=q.arg, opts = opts, nboot=0)
write.csv(get.data(myLHS), file="Hipercubo.csv")
