mprop <- c()
count <- 0
cont <- 0
png("MecanismosFragmentosExt.png")
par(mfrow=c(4,4))
for(i in seq(5, 80, 5))
{

  if(i == 5)
    prop <- "05"
  else
    prop <- paste(i, "", sep = "")
  print(prop)  
  mpais <- c()
  xarea <- c()
  for(j in 1:5)
  {
    conf <- paste(j, "", sep = "")
    nfrag <- strtoi(unlist(strsplit(system(paste("head -n 4 ../output/output-0", prop, conf, "1.txt | tail -n 1", sep=""), intern = TRUE), " "))[3])
    patches <- rep(0, nfrag)
    mrep <- c()
    for(k in 1:5)
    {
      rep <- paste(k, "", sep = "")
      arquivo <- paste("taxaext-output-0", prop, conf, rep, ".txt", sep="")
      dados = file(arquivo, "r")
      lin = readLines(dados, n=1)
      lin = readLines(dados, n=1)
      time <- as.double(unlist(strsplit(system(paste("tail -n 2 ../output/output-0", prop, conf, rep, ".txt | head -n 1", sep=""), intern = TRUE), " "))[1])
      while(lin != "EOF")
      {
      	lin<-strsplit(lin, " ")
      	line <- unlist(lin)
        ptimes <- as.double(line[2:length(line)])
        tax <- 1/(mean(ptimes))
        cont <- cont + length(tax)
        count <- count + length(tax[tax>=8000])
        tax <- tax[tax<=4000]
        if(length(tax)==1)
          patches[strtoi(line[1])] <- patches[strtoi(line[1])]  +  tax
        lin = readLines(dados, n=1)
      }
      #mrep <- c(mrep, mean(patches)) 
    }
    area <- unlist(strsplit(system(paste("head -n 5 ../output/output-0", prop, conf, "1.txt | tail -n 1", sep=""), intern = TRUE), " "))
    area <- area[3:length(area)]
    area <- as.double(area)
    xarea <- c(xarea, area)
    tam <- length(patches)-1
    mpais <- c(mpais, patches[1:tam]/5)

  }
  plot(xarea, mpais, main = prop, xlab = "Area dos Fragmentos", ylab = "Taxa de extincao")
}

dev.off()
