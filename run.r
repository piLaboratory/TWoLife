source("TWoLife.R")
source("Analytical-pckg.R")

params.Combs=read.table("small-LHS.txt",header=T) # Read the LHS table

cb1=paste("00",1:9,sep="")
cb2=paste("0",10:99,sep="")
cb3=paste(100)
CODE=c(cb1,cb2,cb3)

land=Landscape(cover=1,type="b",cell.size=30) # Generates the landscape.

combs=35
for (i in combs)
{
	print(system.time(multiRun(PATH=paste("~/Desktop/Teste-Vchalom/Chalom/Comb-",CODE[i],sep=""),
						 nrep=1,
						 b0=params.Combs$b0[i],
						 d0=params.Combs$d0[i],
						 m0=params.Combs$mov.rate[i],
						 inc.b=params.Combs$sum.slopes[i],
						 inc.d=0,
						 step=params.Combs$step[i],
						 radius=params.Combs$R[i],
						 dens.t=1,
						 config=0,
						 N0=params.Combs$N0[i],
						 lands=land,
						 tm=500)))
}

