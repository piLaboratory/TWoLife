mprop <- c()
for(i in seq(5, 80, 5))
{

  if(i == 5)
    prop <- "05"
  else
    prop <- paste(i, "", sep = "")
  print(prop)  
  mpais <- c()  
  for(j in 1:5)
  {
    conf <- paste(j, "", sep = "")
    mrep <- c()
    for(k in 1:5)
    {
      rep <- paste(k, "", sep = "")
      arquivo <- paste("abundancia-output-0", prop, conf, rep, ".txt", sep="")
      dados = file(arquivo, "r")
      lin = readLines(dados, n=1)
      line <- unlist(lin)
      pop <- strtoi(line[1]) 
      mrep <- c(mrep, pop)
      close(dados)
    }
    mpais <- c(mpais, mean(mrep))
  }
  mprop <- c(mprop, mpais)
}


x <- c()
for(i in seq(5, 80, 5))
  x<- c(x, rep(i, 5))

num <- 10

png(paste("exist-nt-", num, ".png", sep = ""))
plot(x, mprop, ylab = "Abundance", xlab = "Habitat Proportion", main= paste("LHS-", num, sep = ""))
dev.off()

#y <- mprop
#x2 <- x[y>1]
#y <- y[y>1]
#ly <- log(y)
#lx <- log(x2)
#linreg <- lm(ly~lx)
#ty <- exp(coef(linreg)[1])*mprop^coef(linreg)[2]

#png(paste("exist-t-", num, ".png", sep = ""))
#plot(x, ty, ylab = "Abundance", xlab = "Habitat Proportion", main= paste("LHS-", num, sep = ""))
#dev.off()

#lin.mod <- lm(ty~x)
lin.mod <- lm(mprop~x)

png(paste("exist-ntlin-", num, ".png", sep = ""))
plot(x, mprop, xlab = "Habitat Proportion", ylab = "Abundance", main= paste("LHS-", num, sep = ""))
abline(lin.mod)
dev.off()

library("segmented")
segmented.mod <- segmented(lin.mod, seg.Z = ~x,psi=37)

#png(paste("exist-tpw-", num, ".png", sep = ""))
png(paste("exist-ntpw-", num, ".png", sep = ""))
#plot(x, ty, xlab = "Habitat Proportion", ylab = "Abundance", main= paste("LHS-", num, sep = ""))
plot(x, mprop, xlab = "Habitat Proportion", ylab = "Abundance", main= paste("LHS-", num, sep = ""))
plot.segmented(segmented.mod, add=TRUE)
dev.off()


#poly.mod <- lm(ty ~ poly(x, 2))
poly.mod <- lm(mprop ~ poly(x, 2))

#png(paste("exist-tpoly-", num, ".png", sep = ""))
png(paste("exist-ntpoly-", num, ".png", sep = ""))
#plot(x, ty, xlab = "Habitat Proportion", ylab = "Abundance", main= paste("LHS-", num, sep = ""))
plot(x, mprop, xlab = "Habitat Proportion", ylab = "Abundance", main= paste("LHS-", num, sep = ""))
lines(sort(x), fitted(poly.mod)[order(x)])
dev.off()

library("AICcmodavg")

linAICc <- AICc(lin.mod)
polyAICc <- AICc(poly.mod)
segAICc <- AICc(segmented.mod)
dAICcPS <- polyAICc - segAICc
dAICcLS <- linAICc - segAICc

file.create("linreg.txt")
file.create("segreg.txt")
file.create("polyreg.txt")
file.create("AICc.txt")

cat(coef(lin.mod), file = "linreg.txt", append = TRUE, sep = "\n")
cat(coef(segmented.mod), file = "segreg.txt", append = TRUE, sep = "\n")
cat(coef(poly.mod), file = "polyreg.txt", append = TRUE, sep = "\n")
cat(paste("lin vs seg: ", dAICcLS, sep =""), file = "AICc.txt", append = TRUE, sep = "\n")
cat(paste("poly vs seg: ", dAICcPS, sep =""), file = "AICc.txt", append = TRUE, sep = "\n")
cat(paste("AICc seg: ",segAICc, sep =""), file = "AICc.txt", append = TRUE, sep = "\n")
cat(paste("AICc poly: ",polyAICc, sep =""), file = "AICc.txt", append = TRUE, sep = "\n")
cat(paste("AICc lin: ",linAICc, sep =""), file = "AICc.txt", append = TRUE, sep = "\n")



