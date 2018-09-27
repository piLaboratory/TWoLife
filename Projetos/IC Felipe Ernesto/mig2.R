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
      arquivo <- paste("migracao2-output-0", prop, conf, rep, ".txt", sep="")
      dados = file(arquivo, "r")
      lin = readLines(dados, n=1)
      lin = readLines(dados, n=1)
      time <- as.double(unlist(strsplit(system(paste("tail -n 2 ../output/output-0", prop, conf, rep, ".txt | head -n 1", sep=""), intern = TRUE), " "))[1])
      patches <- c()
      while(lin != "EOF")
      {
      	lin<-strsplit(lin, " ")
      	line <- unlist(lin)
      	patches <- c(patches, strtoi(line[2])/time )
        lin = readLines(dados, n=1)
      }
      mrep <- c(mrep, mean(patches))
    }
    mpais <- c(mpais, mean(mrep))
  }
  mprop <- c(mprop, mpais)
}

num <- 3

#print(mprop)
x <- c()
for(i in seq(5, 80, 5))
  x<- c(x, rep(i, 5))
#print(x)
#x<-seq(5,80,5)

y <- mprop[x<=60]
x2 <- x[x<=60]

reg <- lm(y~x2)

png("mig2.png")
plot(x, mprop, ylab = "Migration rate", xlab = "Habitat Proportion", main= paste("LHS-", num, sep = ""))
abline(reg)
dev.off()

file.create("mig2summary.txt")
cat(paste("p-value: ", summary(reg)$coefficients[8], sep = ""), file = "mig2summary.txt", append = TRUE, sep = "\n")
