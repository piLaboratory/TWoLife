multiRun=function(PATH,nrep=10,b0,d0,m0,inc.b,inc.d,step,radius,dens.t,config,N0,lands=land,tm=100)
{
  oldpath <- getwd()
  dir.create(PATH, recursive=TRUE)
  setwd(PATH)
  for(i in 1:nrep)
  {
    TWoLife(taxa.basal=b0,
            taxa.morte=d0,
            move=m0,
            incl.birth=inc.b,
            incl.death=inc.d,
            passo=step,            
            raio=radius,
            density.type=dens.t,
            ini.config=config,
            #####################
            N=N0,
            AngVis=360,
            death.mat=1,
            landscape=lands,
            tempo=tm,            
            out.code=i)
  }
  #Nt.reps=Nt.data(wdir=PATH,intervals=1,tmax=tm)
  #plot.Nt(data=Nt.reps,growth=b0-d0,sum.incl=inc.b+inc.d)
  #return(Nt.reps)
  setwd(oldpath)
}

# land <- Landscape(cover=1,type="b",cell.size=100)
# Dcoef=10000*0.8/4
# sim1=multiRun(PATH="~/Desktop/Comb-0023",
#               nrep=20,
#               b0=0.4,
#               d0=0.1,
#               m0=0.8,
#               inc.b=100450,
#               inc.d=0,
#               step=100,
#               radius=1460,
#               dens.t=1,
#               config=0,
#               N0=50,
#               lands=land,
#               tm=2000
# )


#########################
#  Ancillary functions 1
#########################
Support <- function(taxa.basal=0.6, taxa.morte=0.1, incl.birth=0.5/0.01, 
                    incl.death=0, numb.cells=200, cell.size=2) {
  densi.max = (taxa.basal-taxa.morte)/(incl.birth+incl.death)
  return ((numb.cells*cell.size)^2 * densi.max)
}

####################
#  Get N(t) data
####################
Nt.data=function(wdir,intervals=1,tmax){
  oldpath <- getwd()
  setwd(wdir)
  files=dir()
  Nt.arr=array(0,c(tmax+1,2,length(files)))
  for(i in 1:length(files))
  {
    t1=read.table(files[i],sep=" ")
    Nt=as.data.frame(table(t1[,1]))
    if(is.na(t1[dim(t1)[1],2])) # condicao que indica se houve extincao.
    {
      Nt[dim(Nt)[1],2]=0 # Corrije o data.frame Nt para que N = 0 quando t (coluna 1) = t de extincao. Isso se deve ao fato
      # do table levar em conta o tempo de extincao como um fator representado uma vez na planilha da dados.  
    }
    Nt[,1]=as.double(levels(Nt[,1]))
    if(max(Nt[,1])<=tmax){miss=data.frame(Var1=seq(0,tmax,by=intervals))} else{miss=data.frame(Var1=seq(0,floor(max(Nt[,1])),by=intervals))}
    Nt=merge(miss,Nt,all.x=T,all.y=T)
    
    index=sort(which(is.na(Nt[,2])),decreasing=T)
    if(length(index>0)){
      if(index[1]==dim(Nt)[1])
      {
        Nt[index[1],2]=0
        for(j in index[-1])
        {
          Nt[j,2]=Nt[j+1,2]
        }  
      } else{
        for(j in index)
        {
          Nt[j,2]=Nt[j+1,2]
        }
      }
    }
  Nt.arr[,1,i]=0:tmax
  Nt.arr[,2,i]=Nt[1:(tmax+1),2]
  }
return(Nt.arr)
setwd(oldpath)  
}
# Nt.test=Nt.data("~/Desktop/Comb-0004",tmax=50)

####################
#  Plot N(t) data
####################
plot.Nt=function(dataSim=sim1,growth,sum.incl=0,land.area=10^8,name="N(t)",radius=0)
{
  plot(dataSim[,1,1],dataSim[,2,1],type="n",xlab="t",ylab="N (population size)",ylim=c(0,max(c(max(dataSim[,2,])+10),growth*land.area/sum.incl)),
       cex=1.5, main = name)
  for(i in 1:dim(dataSim)[3])
  {
    lines(dataSim[,1,i],dataSim[,2,i],col="gray50")
  }
  if(sum.incl==0)
  {
    curve(dataSim[1,2,1]*exp(growth*x),add=T,col="red",lwd=4)
    legend(dim(dataSim)[1]/4,max(dataSim[,2,]),legend=c("Observed","Predicted","Replicates"),lty=1,lwd=2,col=c(1,2,"gray50"),bty="n")
    lines(0:(dim(dataSim)[1]-1),apply(dataSim[,2,which(dataSim[dim(dataSim)[1],2,]!=0)],1,mean),col=1,lwd=4)
  } else 
  {
    K=growth*land.area/sum.incl
    lines(0:(dim(dataSim)[1]-1),apply(dataSim[,2,which(dataSim[dim(dataSim)[1],2,]!=0)],1,mean),col=1,lwd=4)
    curve(K/(1+((K/dataSim[1,2,1])-1)*exp(-growth*x)),add=T,col="red",lwd=4)
    legend(dim(dataSim)[1]/4,max(dataSim[,2,])/2,legend=c("Observed","Predicted","Replicates"),lty=1,lwd=2,col=c(1,2,"gray50"),bty="n")
  }
  
}
# plot.Nt(Nt.test,growth=0,0,10^8)
# plot.Nt(growth=0.3,sum.incl=100450)
# sim1=Nt.data(wdir="~/Desktop/Comb-0022",tmax=2000)

#########################
#  Ancillary functions 2
#########################
displacement=function(dado)
{
  disp=sqrt(dado[,3]^2+dado[,4]^2)
  return(disp)
}
dist.front=function(data,ind=10)
{
  vec2=sort(data,decreasing=T)
  if(length(data)>=ind){front=vec2[ind]} else {front=vec2[length(data)]}
  return(front)
}

#########################
#  Get velocity(t) data
#########################
velocity=function(work.dir="~/Desktop/Comb-0003",tmax=50,nrep=20,intervals=1)
{
  oldpath <- getwd()
  setwd(work.dir)
  files=dir()
  vels.arr=matrix(0,tmax,(nrep+1))
  vels.arr[,1]=1:tmax
  for(i in 1:nrep)
  {
    t2=read.table(files[i],sep=" ")
    if(is.na(t2[dim(t2)[1],2])) # condicao que indica se houve extincao.
    {
      t2=t2[-(dim(t2)[1]),] # Corrije o data.frame Nt para que N = 0 quando t (coluna 1) = t de extincao. Isso se deve ao fato
      # do table levar em conta o tempo de extincao como um fator representado uma vez na planilha da dados.  
    }
    vec1=displacement(t2)
    new.t2=data.frame(t2[,1],vec1)
    front=aggregate(new.t2[,2],by = list(new.t2[,1]),FUN=dist.front,ind=10)
    if(dim(front)[1]<tmax+1) # 
    {
      if(dim(front)[1]<(max(front[,1]+1))) # completa os valores faltantes no data frame de distancia da frente (antes da ocorrencia de possiveis extincoes)
      {
        miss=data.frame(Group.1=seq(0,max(front[,1]),by=intervals))# else{miss=data.frame(Var1=seq(0,floor(max(Nt[,1])),by=intervals))}
        front=merge(miss,front,all.x=T,all.y=T)
        index=sort(which(is.na(front[,2])),decreasing=T)
        for(j in index)
        {
          front[j,2]=front[j+1,2]
        }
      }
      
      vels.arr[1:(dim(front)[1]-1),i+1]= diff(front[,2])# preenche a matriz de velocidades antes de um possivel evento de extincao
      #vels.arr[dim(front)[1]:dim(vels.arr)[1],i+1]=0 # preenche o resto da matriz
    } else {vels.arr[,i+1]= diff(front[,2])}# front[-1,2]/front[-1,1]
  }
  return(vels.arr)
  setwd(oldpath)
}

####################
#  Plot velocity(t)
####################
plot.vel=function(object,nrep=20,Dcoef=D,b0=0,d0=0,name="V(t)",yrange=NULL)
{
  growth=b0-d0
  plot(object[,1],object[,2],type="n",main=name,xlab="t",ylab="velocity (distance/time)",ylim=yrange,cex.main=0.7)
  for(i in 1:nrep)
  {
    lines(object[,1],object[,i+1],col="gray50") # replicates
  }
  if(nrep>1){lines(object[,1],apply(object[,-1],1,mean),lwd=4)} # mean replicates (observed)
  if(b0==0 && d0==0){curve(2*sqrt(Dcoef/x),col=2,add=T,lwd=3)} else {abline(h=2*sqrt(growth*Dcoef),col=2,lwd=3)} # expected
}

# teste=velocity("~/Desktop/Comb-0014",nrep=100,tmax=50)
# Dcoef=10000*0.8/4
# plot.vel(teste,Dcoef=Dcoef,nrep=20,growth=0.03)

####################################
#  Plot XY density distribution data
####################################
plot.XYdistrib=function(FILE="output-00001.txt",Dcoef,n.tests=5,tmax=50)
{
  times=sort(sample(1:tmax,n.tests,replace=F))
  t3=read.table(FILE,sep=" ")
  
  par(mfrow=c(n.tests,2))
  for(i in times)
  {
  #x
  hist(t3[which(t3[,1]==i),3],freq=F,breaks=13)
  curve(exp(-(x^2)/(4*Dcoef*i))/sqrt(4*pi*Dcoef*i),col=2,add=T)
  #y
  hist(t3[which(t3[,1]==i),4],freq=F,breaks=13)
  curve(exp(-(x^2)/(4*Dcoef*i))/sqrt(4*pi*Dcoef*i),col=2,add=T)
  }
}

############################
#  Plot sd(x) & sd(y) 
############################
plot.sdXY=function(FILE="output-00001.txt",Dcoef,tmax=50,name="Standard Deviation")
{
  ## Teste desvio padrao
  t3=read.table(FILE,sep=" ")
  #par(mfrow=c(2,1))  
  sd.x=aggregate(t3[,3],list(t3[,1]),sd)
  sd.y=aggregate(t3[,4],list(t3[,1]),sd)
  # sd x
  plot(sd.x[,1],sd.x[,2],xlim=c(0,tmax),type="n",xlab="Time",ylab= "Std. Deviation of X positions",main=name)
  lines(sd.x[,1],sd.x[,2],lwd=3)
  curve(sqrt(2*Dcoef*x),col=2,add=T,lwd=3)
  #sd y
  plot(sd.y[,1],sd.y[,2],xlim=c(0,tmax),type="n",xlab="Time",ylab= "Std. Deviation of Y positions",main=name)
  lines(sd.y[,1],sd.y[,2],lwd=3)
  curve(sqrt(2*Dcoef*x),col=2,add=T,lwd=3)
  #par(mfrow=c(1,1))
}

# Dcoef=10000*0.8/4
# plot.XYdistrib(Dcoef=Dcoef)

################
#  Plot XY data
################

# t3=read.table("output-00001.txt",sep=" ")
# for(i in 1:20){
# plot(t3[which(t3[,1]==i),3],t3[which(t3[,1]==i),4],pch=16)}
# 
# head(t3)
# t3
# TWoPlot <- function(pop, land, col1="gray20", col2="gray70") {
#   n = land$numb.cells
#   s <- seq(-n*land$cell.size/2, n*land$cell.size/2, length=n) # creates the x- and y- sequences for image
#   if (sum(land$scape) == n*n) { 
#     color = col1
#   } else {
#     color = c(col2, col1)
#   }
#   image(s, s, matrix(land$scape,ncol=n), col=color)
#   points(pop[,],pop[,], pch=4, col=2)
# }


###########################
#  Critical Patch Size Test
###########################

#params=c(0.007,0.004,0.8,0,0,500,0,0,0,20,0,1000)
multi.run.cP=function(PATH="~/Desktop/TWoTests/Critical_Patch/Comb-01",nrep=20,b0=0.007,d0=0.004,m0=0.8,inc.b=0,inc.d=0,
                      step=500,radius=0,dens.t=0,config=0,N0=20,tm=1000,lands.range=c(-5000,5000,1000))
{
  oldpath <- getwd()  
  dir.create(PATH)
  setwd(PATH)
  Dcoef=m0*(step^2)/4
  A.crit=2*(pi^2)*Dcoef/(b0-d0)
  land.sides=ceiling(sqrt(A.crit))+seq(lands.range[1],lands.range[2],by=lands.range[3])
  
  for(j in 1:length(land.sides))
  {
    lands=Landscape(cover=1,type="b",cell.size=land.sides[j]/10,numb.cells = 10)
    dir.create(paste(j,land.sides[j],sep="_"))
    setwd(paste(PATH,land.sides[j],sep="/"))          
    
    for(i in 1:nrep)
    {
      TWoLife(taxa.basal=b0,
              taxa.morte=d0,
              move=m0,
              incl.birth=inc.b,
              incl.death=inc.d,
              passo=step,            
              raio=radius,
              density.type=dens.t,
              ini.config=config,
              #####################
              N=N0,
              AngVis=360,
              death.mat=1,
              landscape=lands,
              tempo=tm,            
              out.code=i)
    }          
    setwd(PATH)
  }
  return(land.sides)
  setwd(oldpath)
}

##multi.run.cP(nrep=50)
##multi.run.cP(PATH="TWoTests/Critical_Patch/Comb-03",lands.range=c(2000,17000,1000),nrep=50)

data.Acrit=function(PATH="~/Desktop/TWoTests/Critical_Patch/Comb-01",nreps=50,tmaxi=1000)
{
  directories=paste(PATH,dir(PATH),sep="/")
  Nt.end=matrix(0,nreps,length(directories))
  for(i in 1:length(directories))
  {
    teste=Nt.data(wdir=directories[i],tmax=tmaxi)
    Nt.end[,i]=teste[tmaxi+1,2,]  
  }
  return(list(Nt.end,teste))
}
#testes=data.Acrit("~/Desktop/TWoTests/Critical_Patch/Comb-02")

plot.Acrit=function(object,move=0.8,step=500,growth=0.003,N0=20,range=c(-17000,17000,1000),tmaxi=1000,Ncrit=0)
{
  Dcoef=move*(step^2)/4
  A.crit=2*(pi^2)*Dcoef/growth
  sides=ceiling(sqrt(A.crit)+seq(range[1],range[2],by=range[3]))
  areas=sides^2
  
  par(mfrow=c(1,2))
  plot(rep(areas,each=dim(object)[1]),as.numeric(object),pch=1,xlab="Patch Size",ylab="Population size (Nt =1000)",
       main=paste("Reaction-diffusion: Exponential Growth \n growth = ",growth,", D = ", Dcoef, ", N0 = ", N0,sep=""),
       cex.main=0.8)
  lines(areas,apply(Nt.end,2,mean),col=2)
  abline(v=areas[ceiling(length(areas)/2)],h=N0*exp(growth*tmaxi),lty=c(1,3),col=c(4,1))
  
  
  plot(rep(areas,each=dim(object)[1]),as.numeric(object),type="n",ylim=c(0,1.1),xlab="Patch Size",
       ylab="Extinction Probability for t = 1000",
       main=paste("Reaction-diffusion: Exponential Growth \n growth = ",growth,", D = ", Dcoef, ", N0 = ", N0,sep=""),
       cex.main=0.8)
  for(i in Ncrit)
  {
    indexs=which(object<=i)
    occurence=object
    occurence[indexs]=0
    occurence[-indexs]=1
    points(areas,1-apply(occurence,2,sum)/dim(occurence)[1],col=i+1)
  }
  abline(v=areas[ceiling(length(areas)/2)],h=1,lty=c(1,3),col=c(2,1))
  legend(x =areas[ceiling(length(areas)/2)+2] ,y =0.6 ,"Critical Patch Size",bty="n",lty=3,col=1,cex=0.8)
}

#plot.Acrit(teste[[1]])

#################################
###### Fragmented Landscape Tests
#################################
