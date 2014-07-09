# Função "wrapper" para a chamada em C:
#
# Os passos abaixo foram adaptados de http://users.stat.umn.edu/~geyer/rc/

Sys.setenv("PKG_CPPFLAGS" = "-fopenmp") # liga biblioteca de paralelismo
system("rm TWoLife.{so,o}") #limpa sources velhos
system ("R CMD SHLIB TWoLife.cpp") ## compila no R
dyn.load("TWoLife.so") ## carrega os source resultantes como biblioteca dinamica no R

TWoLife <- function (
    raio, 
    N, 
    AngVis, 
    passo, 
    move, 
    taxa.basal, 
    taxa.morte, 
    incl.birth, 
    incl.death, 
    numb.cells, 
    cell.size, 
    land.shape, 
    density.type,
    death.mat, 
    bound.condition, 
    cover, 
    tempo) 
{
    saida.C <- .C("TWoLife", as.double(raio), as.integer(N),as.double(AngVis), as.double(passo),
                  as.double(move), as.double(taxa.basal), as.double(taxa.morte), 
                  as.double(incl.birth), as.double(incl.death), as.integer(numb.cells), 
                  as.double(cell.size), as.integer(land.shape), as.integer(density.type), 
                  as.double(death.mat), as.integer(bound.condition), as.double(cover),
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
## Uma rodada: coordenadas dos sobreviventes apos t=20
teste1 <- TWoLife(raio=0.1,
                    N=80,
                    AngVis=360,
                    passo=5,
                    move=3,
                    taxa.basal=0.6,
                    taxa.morte=0.1, 
                    incl.birth=0.5/0.01,
                    incl.death=0,
                    numb.cells=200,
                    cell.size=1,
                    land.shape=1, 
                    density.type=0,
                    death.mat=7,
                    bound.condition=0,
                    cover=1,
                    tempo=20)
plot(teste1, xlim=c(-100,100), ylim=c(-100,100))
dim(teste1)
## Tamanho de populacao apos t=6 de 100 repeticoes
pop.size<- numeric()
for (i in 1:100) 
  {
    pop.size[i] = 
      nrow(
          TWoLife(raio=0.1, N=80, AngVis=360, passo=5, move=0.1, taxa.basal=0.6,
                    taxa.morte=0.1, 
                    incl.birth=0.5/0.01, incl.death=0, numb.cells=200, cell.size=1, land.shape=1, 
                    density.type=0, death.mat=7,bound.condition=0, cover=1, tempo=6))
  }

## esperado: capacidade de suporte (numb.cells*cell.size)^2 * densi_max
## densi_max= (taxa_basal-taxa_morte)/(incl_b+incl_d)
200^2*(0.6-0.1)/(0.5/.01)
## Media das simulacoes
pop.size
mean(pop.size)
hist(pop.size)
