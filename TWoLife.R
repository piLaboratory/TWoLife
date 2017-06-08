# Função "wrapper" para a chamada em C:
#
# Os passos abaixo foram adaptados de http://users.stat.umn.edu/~geyer/rc/

Sys.setenv("PKG_CPPFLAGS" = "-fopenmp -DPARALLEL") # liga biblioteca de paralelismo
system("rm TWoLife.so") #limpa sources velhos
system("rm TWoLife.o") #limpa sources velhos
system ("R CMD SHLIB TWoLife.cpp") ## compila no R
dyn.load("TWoLife.so") ## carrega os source resultantes como biblioteca dinamica no R

# Generates the landscape with specified conditions. 
# numb.cells represents both the lenght AND width of the landscape, so numb.cells=100 creates a 100x100 landscape
# Land.shape can be 0 = XXX or 1 = XXX.
# Bound.condition can be 0 = XXX or 1 = XXX. 
Landscape <- function (numb.cells = 100, cell.size = 1, land.shape = 1, type=c("random","blob"), bound.condition=0, cover=1) {
	type=match.arg(type)
	if(cover < 0 || cover > 1) {
		stop("Error creating landscape. Cover must be between 0 and 1")
	}
	# scape represents the actual landscape
	scape <- rep(1, numb.cells*numb.cells)
	if(cover < 1) {
		NtoRemove=round((1-cover)*numb.cells*numb.cells);
		if(type=="random") {
			while(NtoRemove>0)
			{
				i=round(runif(1,0,numb.cells-1));
				j=round(runif(1,0,numb.cells-1));
				# tests to see if this point has already been removed
				if(scape[1+numb.cells*j+i] == 1) {
					NtoRemove = NtoRemove - 1
					scape[1+numb.cells*j+i] = 0
				}
			}
		}
		if(type=="blob") {
			i=round(runif(1,0,numb.cells-1));
			j=round(runif(1,0,numb.cells-1));
			while(NtoRemove>0)
			{
				# tests to see if this point has already been removed
				if(scape[1+numb.cells*j+i] == 1) {
					NtoRemove = NtoRemove - 1
					scape[1+numb.cells*j+i] = 0
				}
				# Draft a new point to be removed (random walk!)
				if(sample(1:2,1) == 1) {
					i = i + (-1)**sample(1:2,1)
				} else {
					j = j + (-1)**sample(1:2,1)
				}
				if(i == -1) { i=numb.cells-1}
				if(i == numb.cells) { i=1}
				if(j == -1) { j=numb.cells-1}
				if(j == numb.cells) { j=1}
			}
		}
	}
	land <- list(numb.cells = numb.cells, cell.size=cell.size, land.shape=land.shape, type=type, bound.condition=bound.condition, cover=cover, scape=scape)
	class(land) <- "landscape"
	return(land)
}

TWoLife <- function (
					 raio=0.1, 
					 N=80, 
					 AngVis=360, 
					 passo=5, 
					 move=0.5, 
					 taxa.basal=0.6, 
					 taxa.morte=0.1, 
					 incl.birth=0.5/0.01, 
					 incl.death=0, 
					 density.type=0,
					 death.mat=7, 
					 landscape,
					 tempo=20,
           ini.config=0,
           out.code=1) 
{
	if(class(landscape) != "landscape") {
		stop("Error in function TWoLife: you must provide a valid landscape. See ?Landscape")
	}
  if(raio>landscape$numb.cells*landscape$cell.size/2)
  {stop("Error in function TWoLife: the radius must be lower than or equal to the half of landscape side (radius <= numb.cells*cell.size/2)")}
	
  saida.C <- .C("TWoLife", 
              as.double(raio),# 1
              as.integer(N),# 2
              as.double(AngVis),# 3
              as.double(passo),# 4
              as.double(move),# 5
              as.double(taxa.basal),# 6
              as.double(taxa.morte),# 7
              as.double(incl.birth),# 8
              as.double(incl.death),# 9
              as.integer(landscape$numb.cells),# 10
              as.double(landscape$cell.size),# 11
              as.integer(landscape$land.shape),# 12
              as.integer(density.type),# 13
              as.double(death.mat), # 14
              as.integer(ini.config), #15
              as.integer(landscape$bound.condition), #16
              as.integer(landscape$scape), #17
              as.double(tempo), #18
              as.integer(0), # 19
              as.double(rep(0, 5000)), # 20
              as.double(rep(0,5000)), # 21 
              as.integer(out.code)
              ## verificar se precisa definir o tamanho e se isto nao dará problemas (dois ultimos argumentos)
				  )
	n <- saida.C[[19]]
	x <- saida.C[[20]]
	y <- saida.C[[21]]
	x <- x[1:n]; y <- y[1:n]  
	return(data.frame(x=x,y=y))
}
npop<-10
 ## Um teste rapido
  land <- Landscape(cover=1,type="b",cell.size=100, bound.condition=1)#Nao usar bound.condition=0 pq essa parte da apply_boundary terá de ser adaptada
 # ## Uma rodada: coordenadas dos sobreviventes apos t=20
 teste <- TWoLife(raio=1560,
 				 N=npop,
 				 AngVis=360,
 				 passo=10,
 				 move=0,
 				 taxa.basal=0.2,
 				 taxa.morte=0, 
 				 incl.birth=1529.076,
 				 incl.death=0,
 				 density.type=1,
 				 death.mat=1,
 				 landscape=land,
 				 tempo=30,
          ini.config=1,
          out.code=234)


 TWoPlot <- function(pop, land, col1="gray20", col2="gray70") {
 	n = land$numb.cells
 	s <- seq(-n*land$cell.size/2, n*land$cell.size/2, length=n) # creates the x- and y- sequences for image
 	if (sum(land$scape) == n*n) { 
 		color = col1
 	} else {
 		color = c(col2, col1)
 	}
 	image(s, s, matrix(land$scape,ncol=n), col=color)
 	points(pop, pch=4, col=2)
 }
 TWoPlot(teste, land)
plot(teste, xlim=c(-100,100), ylim=c(-100,100))
print(dim(teste))

dados = file("output-00234.txt", "r")
dpaisagem = readLines(dados, n=9)
npatches = strtoi(unlist(strsplit(dpaisagem [4], " "))[3])


patches <-matrix( rep(0, 3*npop), nrow = npop, ncol = 3 )
colonizacao <- rep(0, npatches+1)
extincao <- rep(0, npatches+1)
patchPop <- rep(0, npatches+1)
nascimentos <- rep(0, npatches+1)
mortes <- rep(0, npatches+1)

for(i in 1:npop)
{
	lin = readLines(dados, n=1)
	lin<-strsplit(lin, " ")
	line <- unlist(lin)
	patches[i,1]<-strtoi(line[3]) + 1
	patches[i,2] <- 0
	patches[i,3] <- strtoi(line[2])
	patchPop[strtoi(line[3])+1] <- patchPop[strtoi(line[3])+1]+1 
}
maxID = npop

lin = readLines(dados, n=1)
while(lin != "EOF")
{
	lin<-strsplit(lin, " ")
	line <- unlist(lin)
	acao <-strtoi(line[2]) 
	p <- strtoi(line[4])+1 #lembrando que no R os indices não começam do 0
	ID <- strtoi(line[3])
	if(acao == 0)
	{
		patchPop[p] <- patchPop[p] - 1
		ind <- match(ID, patches[,3])
		if(length(patches[,3])>2)
			patches<-patches[-ind,]
		else
		{
				patches<-patches[-ind,]	
				patches<-matrix(patches, nrow=1, ncol=3)
		}

		if(patchPop[p] == 0)
		{
			extincao[p] <- extincao[p] + 1
			mortes[p] <- mortes[p] + 1
		}
	}
	
	if(acao == 1)
	{
		maxID<-maxID+1
		patchPop[p] <- patchPop[p] + 1
		nascimentos[p] <- nascimentos[p] + 1
		patches<-rbind(patches, c(p, 0, maxID))
	}
	if(acao == 2)
	{
		ind <- match(ID, patches[,3])
		if( p!=patches[ind, 1]  )
		{	
			patches[ind,2] <- patches[ind,1]
			patches[ind,1] <- p

			patchPop[patches[ind,1]] <- patchPop[patches[ind,1]] + 1
			patchPop[patches[ind,2]] <- patchPop[patches[ind,2]] - 1

			if(patchPop[patches[ind, 2]]==0)
				extincao[patches[ind, 2]] <- extincao[patches[ind, 2]] + 1
			if(patchPop[p]==1)
				colonizacao[p] <- colonizacao[p] + 1
		}
	}
	lin = readLines(dados, n=1)
}

## Tamanho de populacao apos t=6 de 100 repeticoes
#pop.size<- numeric()
#for (i in 1:20) 
#  {
#    pop.size[i] = 
#      nrow(
#          TWoLife(raio=0.1, N=80, AngVis=360, passo=5, move=0.1, taxa.basal=0.6,
#                    taxa.morte=0.1, 
#                    incl.birth=0.5/0.01, incl.death=0, numb.cells=200, cell.size=1, land.shape=1, 
#                    density.type=0, death.mat=7,bound.condition=0, cover=1, tempo=6))
#  }

## esperado: capacidade de suporte 
# Support <- function(taxa.basal=0.6, taxa.morte=0.1, incl.birth=0.5/0.01, 
# 					incl.death=0, numb.cells=200, cell.size=2) {
# 	densi.max = (taxa.basal-taxa.morte)/(incl.birth+incl.death)
# 	return ((numb.cells*cell.size)^2 * densi.max)
# }
## Media das simulacoes
#print(pop.size - Support())
#print(mean(pop.size - Support()))
