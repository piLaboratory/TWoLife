## source("TWoLife.R") ## calls the C++ sources with parallel.
source("Analytical-pckg.R")
##########################################
########## Set of Parameteres Combinations
##########################################
b0s=c(rep(0.04,6),# only growth tests: exponential
      rep(rep(c(0.4,0.7,0.4),each=2),2),#only growth: logistic - global density
      rep(0.4,6),rep(0.7,2),#only growth: logistic - local density
      rep(0,3),# only diffusion
      rep(c(0.04,0.07,0.04),2), # reaction-diffusion: Skellam model (exponential growth)
      rep(c(rep(0.4,3),0.7),2), # reaction-diffusion: FIsher-Kolmogorov (logistic growth, global density)
      rep(c(rep(0.4,3),0.7),2)) # reaction-diffusion: FIsher-Kolmogorov (logistic growth, local density)
d0s=c(rep(0.01,6),# only growth tests: exponential
      rep(rep(c(0.1,0.4,0.1),each=2),2),#only growth: logistic - global density
      rep(0.1,6),rep(0.4,2),#only growth: logistic - local density
      rep(0,3),# only diffusion
      rep(c(0.01,0.04,0.01),2), # reaction-diffusion: Skellam model (exponential growth)
      rep(c(rep(0.1,3),0.4),2), # reaction-diffusion: FIsher-Kolmogorov (logistic growth, global density)
      rep(c(rep(0.1,3),0.4),2)) # reaction-diffusion: FIsher-Kolmogorov (logistic growth, local density)
m0s=c(rep(0,26), # only growth
     0.8,0.8,0.4, # only diffusion
     rep(c(rep(0.8,2),1.6),2), # reaction-diffusion: Skellam model (exponential growth)
     rep(0.8,16))
incl.bs=c(rep(0,6), # only growth tests: exponential
          rep(c(rep(10^5,4),rep(0.5*10^5,2)),2), # only growth: logistic - global density
          rep(100500,2),rep(1750,2),rep(100500,4), # only growth: logistic - local density
          0,0,0, # only diffusion
          rep(0,6), # reaction-diffusion: Skellam model (exponential growth)
          rep(c(rep(100500,2),50250,100500),2), # reaction-diffusion: FIsher-Kolmogorov (logistic growth, global density)
          rep(c(rep(100500,2),50250,100500),2)) # reaction-diffusion: FIsher-Kolmogorov (logistic growth, local density)
incl.ds=c(rep(0,6), # only growth tests: exponential
          rep(c(rep(0,4),rep(0.5*10^5,2)),2), # only growth: logistic - global density
          rep(0,8), # only growth: logistic - local density
          0,0,0, # only diffusion
          rep(0,6), # reaction-diffusion: Skellam model (exponential growth)
          rep(c(rep(0,2),50250,0),2), # reaction-diffusion: FIsher-Kolmogorov (logistic growth, global density)
          rep(c(rep(0,2),50250,0),2)) # reaction-diffusion: FIsher-Kolmogorov (logistic growth, local density)
steps=c(rep(0,26), # only growth
        rep(100,3), # only diffusion
        rep(100,6), # reaction-diffusion: Skellam model (exponential growth)
        rep(100,16)) # reaction-diffusion: Fisher-Kolmogorov model (logistic growth; global and local density)
radii=c(rep(0,18), # only growth: exponential + logistic global density
        750,1500,rep(750,6), # only growth: logistic local density
        rep(0,3), # only diffusion
        rep(0,6), # reaction-diffusion: Skellam model (exponential growth)
        rep(0,8), # reaction-diffusion: Fisher-Kolmogorov model (logistic growth; global density)
        rep(c(750,1500,750,750),2)) # reaction-diffusion: Fisher-Kolmogorov model (logistic growth; local density)
d.types=c(rep(0,18), # only growth: exponential + logistic global density
          rep(1,8), # only growth: logistic local density
          rep(0,3), # only diffusion
          rep(0,6), # reaction-diffusion: Skellam model (exponential growth)
          rep(0,8), # reaction-diffusion: Fisher-Kolmogorov model (logistic growth; global density)
          rep(1,8))
i.configs=c(rep(0:2,2),
            rep(c(0,1),each=6),
            rep(0,4),rep(1,2),0,0,
            0,2,0,
            rep(c(0,1),each=3),
            rep(c(0,1),each=4),
            rep(c(0,1),each=4))
N0s=c(rep(c(50,10),each=3),
      rep(c(50,500),6),
      rep(50,3),500,rep(50,4),
      rep(500,3),
      rep(50,6),
      rep(50,16))
d.matrix=rep(0,51)
times=c(rep(c(50,100),each=3),
        rep(2000,20),
        rep(100,3),
        rep(50,6),
        rep(2000,16))
#n.replicates=c()

cb1=paste("0",1:9,sep="")
cb2=paste(10:51)
CODE=c(cb1,cb2)

# Parameter combinations table for testing TwoLife
params=data.frame(b0=b0s,d0=d0s,movem.rate=m0s,incl.birth=incl.bs,incl.death=incl.ds,step.length=steps,radius=radii,
                  dens.type=d.types,initial.config=i.configs,N0=N0s,death.matrix=d.matrix,tmax=times)
# params$b0[51]
# dim(params)
## Create a directory for the output files, here named TWotests
land <- Landscape(cover=1,type="b",cell.size=100)
combs=1:51
for (i in combs)
{
  multiRun(PATH=paste("./TWoTests/Comb-",CODE[i],sep=""),
                nrep=20,
                b0=params$b0[i],
                d0=params$d0[i],
                m0=params$movem.rate[i],
                inc.b=params$incl.birth[i],
                inc.d=params$incl.death[i],
                step=params$step.length[i],
                radius=params$radius[i],
                dens.t=params$dens.type[i],
                config=params$initial.config[i],
                N0=params$N0[i],
                lands=land,
                tm=params$tmax[i])
}

###################################
###################################
###################################



######################
##### Plot main titles
######################
titles.exp=paste("Simple Birth-Death (COMB ", CODE[1:6],") \n", expression(r),
                 paste(" = ", params[1:6,1]-params[1:6,2],", ini.config = ",params[1:6,9],sep=""), sep="")
titles.log1=paste("General Birth-Death - Logistic Growth (COMB ", CODE[7:18],") \n", expression(b0),
                 paste(" =", params[7:18,1],", "), expression(slope(b)),paste(" =", params[7:18,4],", "),
                 expression(d0), paste(" =", params[7:18,2],", "),expression(slope(d)),paste(" =", params[7:18,5],"\n"),
                 paste("ini.config = ",params[7:18,9],", "),expression(N0),paste(" =", params[7:18,10]),sep="")
titles.log2= paste("General Birth-Death - Logistic Growth (COMB ", CODE[19:26],") \n",
                      expression(b0), paste(" =", params[19:26,1],", "), expression(slope(b)),paste(" =", params[19:26,4],", "),
                      expression(d0), paste(" =", params[19:26,2],", "), expression(slope(d)),paste(" =", params[19:26,5],"\n"),
                      "movement rate = ", params[19:26,3], ", step length = ", params[19:26,6], ", R = ", params[19:26,7],
                      ", ini.config = ",params[19:26,9],sep="")  
titles.mov=paste("Simple Random Walk (COMB ", CODE[27:29],") \n", "movement rate = ", params[27:29,3],
                 ", step length = ", params[27:29,6], ", ini.config = ",params[27:29,9], sep="")
titles.skellam= paste("Reaction-Diffusion with Exponential Growth (COMB ", CODE[30:35],") \n",
                      expression(b0), paste(" =", params[30:35,1],", "), expression(d0), paste(" =", params[30:35,2],"\n"),
                      "movement rate = ", params[30:35,3], ", step length = ", params[30:35,6],
                      ", ini.config = ",params[30:35,9],sep="")
titles.fisher1= paste("Reaction-Diffusion with Logistic Growth (COMB ", CODE[36:43],") \n",
                      expression(b0), paste(" =", params[36:43,1],", "), expression(slope(b)),paste(" =", params[36:43,4],", "),
                      expression(d0), paste(" =", params[36:43,2],", "), expression(slope(d)),paste(" =", params[36:43,5],"\n"),
                      "movement rate = ", params[36:43,3], ", step length = ", params[36:43,6],
                      ", ini.config = ",params[36:43,9],sep="")
titles.fisher2= paste("Reaction-Diffusion with Logistic Growth (COMB ", CODE[44:51],") \n",
                      expression(b0), paste(" =", params[44:51,1],", "), expression(slope(b)),paste(" =", params[44:51,4],", "),
                      expression(d0), paste(" =", params[44:51,2],", "), expression(slope(d)),paste(" =", params[44:51,5],"\n"),
                      "movement rate = ", params[44:51,3], ", step length = ", params[44:51,6], ", R = ", params[44:51,7],
                      ", ini.config = ",params[44:51,9],sep="")
titles=c(titles.exp,titles.log1,titles.log2,titles.mov,titles.skellam,titles.fisher1,titles.fisher2)
titles                      
###########


###############
####### Plots
## First create a directory for the plots (./TWoResults here)
###############
paths=paste("./TWoTests/Comb-",CODE,sep="")
paths
ind1=c(1,4,2,5,3,6)

pdf("./TWoResults/Nt/Nt_Simple-Birth-Death.pdf",width = 8,height = 12)
par(mfrow=c(3,2))
for(i in ind1)
{
  teste=Nt.data(wdir = paths[i],tmax=params$tmax[i])
  setwd("./TWoResults/Nt")
  plot.Nt(teste,growth=params$b0[i]-params$d0[i],sum.incl = params$incl.birth[i]+params$incl.death[i],name = titles[i])
}
par(mfrow=c(1,1))
dev.off()

ind=c(19:23,25)

pdf("./TWoResults/Nt/Nt_General-Birth-Death-gd2.pdf",width = 8,height = 12)
par(mfrow=c(3,2))
for(i in ind)
{
  teste=Nt.data(wdir = paths[i],tmax=params$tmax[i])
  setwd("./TWoResults/Nt")
  plot.Nt(teste,growth=params$b0[i]-params$d0[i],sum.incl = params$incl.birth[i]+params$incl.death[i],name = titles[i])
}
par(mfrow=c(1,1))
dev.off()

ind=c(19:23,25)

pdf("./TWoResults/Nt/Nt_General-Birth-Death-ld.pdf",width = 8,height = 12)
par(mfrow=c(3,2))
for(i in ind)
{
  teste=Nt.data(wdir = paths[i],tmax=params$tmax[i])
  setwd("./TWoResults/Nt")
  plot.Nt(teste,growth=params$b0[i]-params$d0[i],sum.incl = params$incl.birth[i]+params$incl.death[i],name = titles[i])
}
par(mfrow=c(1,1))
dev.off()



#################
#################
ind=27:29
pdf("./TWoResults/Velocity/Vel_Simple-Random-Walk.pdf",width = 8,height = 8)
par(mfrow=c(2,2))
for(i in ind)
{
  teste=velocity(work.dir = paths[i],tmax = params$tmax[i],nrep = 50)
  setwd("./TWoResults/Velocity")
  plot.vel(teste,Dcoef = (params$step.length[i]^2)*params$movem.rate[i]/4, growth=params$b0[i]-params$d0[i],
           name = titles[i])
}
par(mfrow=c(1,1))
dev.off()
#######################
pdf("./TWoResults/Distribution/Dist_Simple-Random-Walk.pdf",width = 8,height = 12)
par(mfrow=c(length(ind),2))
for(i in ind)
{
  Dcoefs=(params$step.length[i]^2)*params$movem.rate[i]/4
  plot.sdXY(FILE = paste(paths[i],"/output-00001.txt",sep=""),Dcoef = Dcoefs,tmax = params$tmax[i],name = titles[i])
}
par(mfrow=c(1,1))
dev.off()

###############
###############
ind=c(30,33,31,34,32,35)
pdf("./TWoResults/Velocity/Vel_Skellam-b.pdf",width = 8,height = 12)
par(mfrow=c(3,2))
for(i in ind)
{
  Dcoefs=(params$step.length[i]^2)*params$movem.rate[i]/4
  growths=params$b0[i]-params$d0[i]
  teste=velocity(work.dir = paths[i],tmax = params$tmax[i],nrep = 50)
  setwd("./TWoResults/Velocity")
  plot.vel(teste,Dcoef = Dcoefs, b0=params$b0[i],d0=params$d0[i], yrange = c(-1,30),name = titles[i])
}
par(mfrow=c(1,1))
dev.off()
###
pdf("./TWoResults/Nt/Nt_Skellam.pdf",width = 8,height = 12)
par(mfrow=c(3,2))
for(i in ind)
{
  teste=Nt.data(wdir = paths[i],tmax=params$tmax[i])
  setwd("./TWoResults/Nt")
  plot.Nt(teste,growth=params$b0[i]-params$d0[i],sum.incl = params$incl.birth[i]+params$incl.death[i],name = titles[i])
}
par(mfrow=c(1,1))
dev.off()


ind=c(36:43)
pdf("./TWoResults/Velocity/Vel_Fisher-gd-b.pdf",width = 8,height = 16)
par(mfrow=c(4,2))
for(i in ind)
{
  Dcoefs=(params$step.length[i]^2)*params$movem.rate[i]/4
  growths=params$b0[i]-params$d0[i]
  teste=velocity(work.dir = paths[i],tmax = params$tmax[i],nrep = 20)
  setwd("./TWoResults/Velocity")
  plot.vel(teste,Dcoef = Dcoefs, b0=params$b0[i],d0=params$d0[i], yrange = c(-100,100),name = titles[i])    
}
par(mfrow=c(1,1))
dev.off()
pdf("./TWoResults/Nt/Nt_Fisher-gd.pdf",width = 8,height = 12)
par(mfrow=c(4,2))
for(i in ind)
{
  teste=Nt.data(wdir = paths[i],tmax=params$tmax[i])
  setwd("./TWoResults/Nt")
  plot.Nt(teste,growth=params$b0[i]-params$d0[i],sum.incl = params$incl.birth[i]+params$incl.death[i],name = titles[i])
}
par(mfrow=c(1,1))
dev.off()

ind=c(44:47)
pdf("./TWoResults/Velocity/Vel_Fisher-ld.pdf",width = 8,height = 8)
par(mfrow=c(2,2))
for(i in ind)
{
  teste=velocity(work.dir = paths[i],tmax = params$tmax[i],nrep = 20)
  setwd("./TWoResults/Velocity")
  plot.vel(teste,Dcoef = (params$step.length[i]^2)*params$movem.rate[i]/4, growth=params$b0[i]-params$d0[i],name = titles[i])
}
par(mfrow=c(1,1))
dev.off()
pdf("./TWoResults/Nt/Nt_Fisher-ld.pdf",width = 8,height = 8)
par(mfrow=c(2,2))
for(i in ind)
{
  teste=Nt.data(wdir = paths[i],tmax=params$tmax[i])
  setwd("./TWoResults/Nt")
  plot.Nt(teste,growth=params$b0[i]-params$d0[i],sum.incl = params$incl.birth[i]+params$incl.death[i],name = titles[i])
}
par(mfrow=c(1,1))
dev.off()

ind=c(48:51)
pdf("./TWoResults/Velocity/Vel_Fisher-ld-c2b.pdf",width = 8,height = 8)
par(mfrow=c(2,2))
for(i in ind)
{
  Dcoefs=(params$step.length[i]^2)*params$movem.rate[i]/4
  growths=params$b0[i]-params$d0[i]
  teste=velocity(work.dir = paths[i],tmax = params$tmax[i],nrep = 20)
  setwd("./TWoResults/Velocity")
  plot.vel(teste,Dcoef = Dcoefs, b0=params$b0[i],d0=params$d0[i], yrange = c(-150+2*sqrt(growths*Dcoefs),50+2*sqrt(growths*Dcoefs)),
           name = titles[i])
}
par(mfrow=c(1,1))
dev.off()
pdf("./TWoResults/Nt/Nt_Fisher-ld-c2.pdf",width = 8,height = 8)
par(mfrow=c(2,2)) 
for(i in ind)
{
  teste=Nt.data(wdir = paths[i],tmax=params$tmax[i])
  setwd("./TWoResults/Nt")
  plot.Nt(teste,growth=params$b0[i]-params$d0[i],sum.incl = params$incl.birth[i]+params$incl.death[i],
          name = titles[i])
}
par(mfrow=c(1,1))
dev.off()

#####################
#####################
testes=data.Acrit("./TWoTests/Critical_Patch/Comb-02")

pdf("./TWoResults/Acrit/Skellam_Comb-02.pdf",width = 8,height = 5)
plot.Acrit(teste[[1]])
dev.off()

