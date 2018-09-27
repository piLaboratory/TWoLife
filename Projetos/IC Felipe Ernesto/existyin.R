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
for(i in seq(.05, .8, .05))
  x<- c(x, rep(i, 5))
x <- x/0.8

mprop <- mprop/max(mprop)
num <- 10

png(paste("exist-yin-", num, ".png", sep = ""))
plot(x, mprop, ylab = "Abundance", xlab = "Relative Habitat Proportion", main= paste("LHS-", num, sep = ""))
#dev.off()

lin.mod <- lm(mprop~x)

library("segmented")
segmented.mod <- segmented(lin.mod, seg.Z = ~x,psi=.37)

bp <- summary(segmented.mod)$psi[1,2] 

x2 <- x[x>bp]
y2 <- mprop[x>bp]
#x2 <- x
#y2 <- mprop


poly.mod <- lm(y2 ~ poly(x2, 2, raw = TRUE))
lines(sort(x2), fitted(poly.mod)[order(x2)])
dev.off()

xr1 <- (1-as.double(coef(poly.mod)[2])) / (2*as.double(coef(poly.mod)[3]))
xr2 <- (-1-as.double(coef(poly.mod)[2])) / (2*as.double(coef(poly.mod)[3]))

file.create("yin.txt")
cat(xr1, file = "yin.txt", append = TRUE, sep = "\n")
cat(xr2, file = "yin.txt", append = TRUE, sep = "\n")
