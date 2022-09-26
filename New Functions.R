indspec<-function(N, TNW , BIC , Method = 2 ){
  
  
  
  WIC<- TNW - BIC
  
  
  ind<- matrix(nrow = N, ncol = 2) 
  colnames(ind) <- c("genotype_Mean", "width_sd") 
  
  
  
  x<- runif(1,0,1-TNW)
  
  
  #Random
  
  if(Method ==  0 || Method == "random" || Method == "Random"){
    
    for (i in 1:N) {
      
      for (i in 1:N) {
        ind[i,1]<- runif(1,(x+WIC/2),(x+TNW-WIC/2))
      }
    }
    
  }
  
  
  #Sistematical
  
  if(Method == 1 || Method =="sistematical"  || Method == "Sistematical"){
    
    ind[,1] <-seq(from = x+WIC/2, to= x+TNW-WIC/2, length.out = N)
    
  }
  
  #Normal
  
  if(Method ==  2 || Method == "Normal" ||  Method == "Normal"){
    
    for (i in 1:N) {
      
      ind[i]<- rnorm(1,(x+TNW/2),BIC/2)
      while (ind[i,1]<x+WIC/2 || ind[i,1]>x+TNW-WIC/2) {
        ind[i]<- rnorm(1,(x+TNW/2),BIC/2)
      }
    }
  }
  
  ind[,2]<-WIC/2
  
  return(ind)
}

show_ind_fitness_curves<-function(nInd,TNW,BIC, Method=2){
  
  
  inds<- indspec(nInd,TNW,BIC, Method)
  
  x<- seq(from= 0, to= 1, length.out = 1000) # Generates the z_scores
  y = matrix(nrow=1000, ncol=nInd) # creates a matrix for storing the Distributions
  
  for (i in 1:nrow(inds)) {
    
    y[,i] = dnorm(x,inds[i,1], inds[i,2])
  }
  
  Y=apply(y, 1, sum)
  
  plot(Y, col="darkblue", type = "line")
  lines(nInd*dnorm(x,mean = mean(inds[,1]),TNW/2),col="grey") # Plots the Distribution of individual optimum trait values
  lines(Y/nInd, type = "line", col="purple")
  for (i in 1:nrow(inds)) {
    lines(y[,i])
  }
  
  info<-list()
  
  info[[1]]<-inds
  info[[2]]<-Y/nInd
  
  names(info)<-c("inds","pop")
  
  return(info)
}
