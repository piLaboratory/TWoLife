mprop <- c()
#png("MecanismosFragmentos.png")
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
      arquivo <- paste("migracao2-output-0", prop, conf, rep, ".txt", sep="")
      dados = file(arquivo, "r")
      lin = readLines(dados, n=1)
      lin = readLines(dados, n=1)
      time <- as.double(unlist(strsplit(system(paste("tail -n 2 ../output/output-0", prop, conf, rep, ".txt | head -n 1", sep=""), intern = TRUE), " "))[1])
      while(lin != "EOF")
      {
      	lin<-strsplit(lin, " ")
      	line <- unlist(lin)
      	patches[strtoi(line[1])] <- patches[strtoi(line[1])]  +  strtoi(line[2])/time
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
  plot(xarea, mpais, main = prop, xlab = "Area dos Fragmentos", ylab = "Taxa de migracao")
}

#dev.off()