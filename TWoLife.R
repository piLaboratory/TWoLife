# Função "wrapper" para a chamada em C:
#
# Os passos abaixo foram adaptados de http://users.stat.umn.edu/~geyer/rc/

Sys.setenv("PKG_CPPFLAGS" = "-fopenmp") # liga biblioteca de paralelismo
system("rm TWoLife.{so,o}") #limpa sources velhos
system ("R CMD SHLIB TWoLife.cpp") ## compila no R
dyn.load("TWoLife.so") ## carrega os source resultantes como biblioteca dinamica no R

# Generates the landscape with specified conditions. 
# numb.cells represents both the lenght AND width of the landscape, so numb.cells=100 creates a 100x100 landscape
# Land.shape can be 0 = XXX or 1 = XXX.
# Bound.condition can be 0 = XXX or 1 = XXX. The only type supported right now is random. 
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
					 tempo=20) 
{
	if(class(landscape) != "landscape") {
		stop("Error in function TWoLife: you must provide a valid landscape. See ?Landscape")
	}
	saida.C <- .C("TWoLife", as.double(raio), as.integer(N),as.double(AngVis), as.double(passo),
				  as.double(move), as.double(taxa.basal), as.double(taxa.morte), 
				  as.double(incl.birth), as.double(incl.death), as.integer(landscape$numb.cells), 
				  as.double(landscape$cell.size), as.integer(landscape$land.shape), as.integer(density.type), 
				  as.double(death.mat), as.integer(landscape$bound.condition), as.integer(landscape$scape),
				  as.double(tempo), as.integer(0),
				  as.double(rep(0, 5000)), as.double(rep(0,5000)) ## verificar se precisa definir o tamanho e se isto nao dará problemas
				  )
	n <- saida.C[[18]]
	x <- saida.C[[19]]
	y <- saida.C[[20]]
	x <- x[1:n]; y <- y[1:n]  
	return(data.frame(x=x,y=y))
}

## Um teste rapido
land <- Landscape(cover=0.5)
## Uma rodada: coordenadas dos sobreviventes apos t=20
teste <- TWoLife(raio=0.1,
				 N=80,
				 AngVis=360,
				 passo=5,
				 move=3,
				 taxa.basal=0.6,
				 taxa.morte=0.1, 
				 incl.birth=0.5/0.01,
				 incl.death=0,
				 density.type=0,
				 death.mat=7,
				 landscape=land,
				 tempo=10)
TWoPlot <- function(pop, land) {
	n = land$numb.cells
	s <- seq(-n*land$cell.size/2, n*land$cell.size/2, length=n) # creates the x- and y- sequences for image
	if (sum(land$scape) == n*n) { 
		color = "gray20"
	} else {
		color = c("gray70", "gray20")
	}
	image(s, s, matrix(land$scape,ncol=n), col=color)
	points(pop, pch=4, col=2)
}
#plot(teste1, xlim=c(-100,100), ylim=c(-100,100))
#dim(teste1)
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
Support <- function(taxa.basal=0.6, taxa.morte=0.1, incl.birth=0.5/0.01, 
					incl.death=0, numb.cells=200, cell.size=2) {
	densi.max = (taxa.basal-taxa.morte)/(incl.birth+incl.death)
	return ((numb.cells*cell.size)^2 * densi.max)
}
## Media das simulacoes
#print(pop.size - Support())
#print(mean(pop.size - Support()))
