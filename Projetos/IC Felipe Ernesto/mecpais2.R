mprop <- c()
count <- 0
cont <- 0
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
      arquivo <- paste("taxacol-output-0", prop, conf, rep, ".txt", sep="")
      dados = file(arquivo, "r")
      patches <- c()
      lin = readLines(dados, n=1)
      
      if( strtoi(unlist(strsplit(lin, " "))[1]) == 0)
        lin = readLines(dados, n=1)
      while(lin != "EOF")
      {
      	lin<-strsplit(lin, " ")
      	line <- unlist(lin)
        ptimes <- as.double(line[2:length(line)])
        tax <- 1/(mean(ptimes))
        cont <- cont + length(tax)
        count <- count + length(tax[tax>=8000])
        tax <- tax[tax<=4000]
      	patches <- c(patches, tax )
        lin = readLines(dados, n=1)
      }
      mrep <- c(mrep, mean(patches)) 
        
    }
    mpais <- c(mpais, mean(mrep))
  }
  mprop <- c(mprop, mpais)
}

num <- 10

print(mprop)

x <- c()
for(i in seq(5, 80, 5))
  x<- c(x, rep(i, 5))
print(x)

y <- mprop[x<=60]
x2 <- x[x<=60]
reg <- lm(y~x2)

png(paste0("col-", num, ".png"))
plot(x, mprop, xlab = "Habitat Proportion", ylab = "Colonization rate", main= paste("LHS-", num, sep = ""))
abline(reg)
dev.off()

file.create("colsummary.txt")
cat(paste("p-value: ", summary(reg)$coefficients[8], sep = ""), file = "colsummary.txt", append = TRUE, sep = "\n")

