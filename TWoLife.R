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
readLandscape <- function (prop, label)
{
	land <- rep(-1, 512*512)
	nome <- paste("../paisagens-naturais/", prop, "_", label, "_class_1.txt", sep="")
	paisagem <- file(nome, "r")
	for(i in seq(0, 511))
	{
		lin = readLines(paisagem, n=1)
		lin <- strsplit(lin, " ")
		line <- unlist(lin)
		for(j in seq(0, 511))
		{
			land[1+512*j+i] <- strtoi(line[j+1])
		}
	}
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

projetoFelipe<-function()
{#Lembrar de citar Mandai em qualquer produto deste trabalho

analises<-function(arquivo, npop){

dados = file(arquivo, "r")
dpaisagem = readLines(dados, n=9)
npatches = strtoi(unlist(strsplit(dpaisagem [4], " "))[3])

IDlist<-rep(0, npop)
patches <-matrix( rep(0, 3*npop), nrow = npop, ncol = 3 )
colonizacao <- rep(0, npatches+1) #Conta sempre que um individuo chega num fragmento vazio
extincao <- rep(0, npatches+1)#Conta sempre que um fragmento fica vazio (Note que é de se esperar que as entradas de colonizacao e de extincao difiram apenas de +-1)
migracao1 <- rep(0, npatches+1) #Conta sempre que um individuo entra num fragmento, estando ele vazio ou nao. (Note que o individuo pode ficar entrando e saindo do fragmento e sempre será contado)
migracao2 <- rep(0, npatches+1)# Conta apenas quando o individuo vem de um fragmento diferente, ou seja, nao conta quando o individuo vai para a matriz e volta
emigracao <- rep(0, npatches+1)#Conta sempre que um individuo sai da paisagem

patchPop <- rep(0, npatches+1)
nascimentos <- rep(0, npatches+1)
mortes <- rep(0, npatches+1)

for(i in 1:npop)
{
	lin = readLines(dados, n=1)
	lin<-strsplit(lin, " ")
	line <- unlist(lin)
	patches[i,1]<-strtoi(line[3]) + 1
	patches[i,2] <- -1
	patches[i,3] <- -1
	IDlist[i] <- strtoi(line[2])
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
		mortes[p] <- mortes[p] + 1
		patchPop[p] <- patchPop[p] - 1
		ind <- match(ID, IDlist)
		if(length(IDlist)>2)
		{
			patches<-patches[-ind,]
			IDlist<-IDlist[-ind]
		}
		else
		{
				patches<-patches[-ind,]
				patches<-matrix(patches, nrow=1, ncol=3)
				IDlist<-IDlist[-ind]
		}

		if(patchPop[p] == 0)
			extincao[p] <- extincao[p] + 1

	}

	if(acao == 1)
	{
		maxID<-maxID+1
		patchPop[p] <- patchPop[p] + 1
		nascimentos[p] <- nascimentos[p] + 1
		patches<-rbind(patches, c(p, -1, -1))
		IDlist<- c(IDlist, maxID)
	}
	if(acao == 2)
	{
		ind <- match(ID, IDlist)
		if( p!=patches[ind, 1]  )
		{
			patches[ind,3]<- patches[ind,2]
			patches[ind,2] <- patches[ind,1]
			patches[ind,1] <- p

			patchPop[patches[ind,1]] <- patchPop[patches[ind,1]] + 1
			patchPop[patches[ind,2]] <- patchPop[patches[ind,2]] - 1

			if(patchPop[patches[ind, 2]]==0)
				extincao[patches[ind, 2]] <- extincao[patches[ind, 2]] + 1

			if(patchPop[p]==1)
				colonizacao[p] <- colonizacao[p] + 1

			migracao1[p] <- migracao1[p] + 1

			if(patches[ind,2]>1 || patches[ind,2]!=1 || patches[ind,3]!=patches[ind,1]  )
				migracao2[p] <- migracao2[p] + 1

		}
	}

	if(acao == 3)
	{
		emigracao[p] <- emigracao[p] + 1
		patchPop[p] <- patchPop[p] - 1
		ind <- match(ID, IDlist)
		if(length(IDlist)>2)
		{
			patches<-patches[-ind,]
			IDlist<-IDlist[-ind]
		}
		else
		{
				patches<-patches[-ind,]
				patches<-matrix(patches, nrow=1, ncol=3)
				IDlist<-IDlist[-ind]
		}

		if(patchPop[p] == 0)
			extincao[p] <- extincao[p] + 1

	}

	lin = readLines(dados, n=1)
}

arqm<-paste0("mortes-", arquivo)
arqn<-paste0("nascimentos-", arquivo)
arqe<-paste0("extincao-", arquivo)
arqc<-paste0("colonizacao-", arquivo)
arqm1<-paste0("migracao1-", arquivo)
arqm2<-paste0("migracao2-", arquivo)
arqem<-paste0("emigracao-", arquivo)

file.create(arqm)
file.create(arqn)
file.create(arqe)
file.create(arqc)
file.create(arqm1)
file.create(arqm2)
file.create(arqem)


for (i in 0:npatches)
{

	cat(paste(i,mortes[i+1]), file = arqm, append = TRUE, sep = "\n")
	cat(paste(i,nascimentos[i+1]), file = arqn, append = TRUE, sep = "\n")
	cat(paste(i,extincao[i+1]), file = arqe, append = TRUE, sep = "\n")
	cat(paste(i,colonizacao[i+1]), file = arqc, append = TRUE, sep = "\n")
	cat(paste(i,migracao1[i+1]), file = arqm1, append = TRUE, sep = "\n")
	cat(paste(i,migracao2[i+1]), file = arqm2, append = TRUE, sep = "\n")
	cat(paste(i,emigracao[i+1]), file = arqem, append = TRUE, sep = "\n")
}

cat("EOF", file = arqm, append = TRUE, sep = "\n")
cat("EOF", file = arqn, append = TRUE, sep = "\n")
cat("EOF", file = arqe, append = TRUE, sep = "\n")
cat("EOF", file = arqc, append = TRUE, sep = "\n")
cat("EOF", file = arqm1, append = TRUE, sep = "\n")
cat("EOF", file = arqm2, append = TRUE, sep = "\n")
cat("EOF", file = arqem, append = TRUE, sep = "\n")

}

library("pse")

#Diretório que receberá os outputs para esta combinação de parâmetros
#Sempre que for começar a rodar lembrar de colocar "run 1" em iteration.txt
ite <- file("iteration.txt", "r")
linha <- readLines(ite, n=1)
nite <- strtoi(unlist(strsplit(linha, " "))[2])
dir <- paste("LHS-", nite, sep = "")
dir.create(dir)
cat(paste("run ", nite+1, sep = ""), file = "iteration.txt", append = FALSE, sep = "\n")

#Leitura dos parâmetros de "Hipercubo.csv"
LHS <- file("Hipercubo.csv", "r")
cab <- readLines(LHS, n=nite)
param <- readLines(LHS, n=1)
param <- strsplit(param, ",")
param <- unlist(param)

setwd(dir)

#Arquivo metadata com os valores dos parâmetros
file.create("METADATA.txt")
cat(date(), file = "METADATA.txt", append = TRUE, sep = "\n")
cat("", file = "METADATA.txt", append = TRUE, sep = "\n")
cat(paste("logR:", param[2]), file = "METADATA.txt", append = TRUE, sep = "\n")
cat(paste("logs:", param[3]), file = "METADATA.txt", append = TRUE, sep = "\n")
cat(paste("logdelta:", param[4]), file = "METADATA.txt", append = TRUE, sep = "\n")
cat(paste("loglamb0:", param[5]), file = "METADATA.txt", append = TRUE, sep = "\n")
cat(paste("logmu0:", param[6]), file = "METADATA.txt", append = TRUE, sep = "\n")
cat(paste("logb:", param[7]), file = "METADATA.txt", append = TRUE, sep = "\n")
cat(paste("logm:", param[8]), file = "METADATA.txt", append = TRUE, sep = "\n")
cat(paste("dm:", param[9]), file = "METADATA.txt", append = TRUE, sep = "\n")
cat("EOF", file = "METADATA.txt", append = TRUE, sep = "\n")

R <- exp(as.double(param[2]))
s <- exp(as.double(param[3]))
delta <- exp(as.double(param[4]))
lamb0 <- exp(as.double(param[5]))
mu0 <- exp(as.double(param[6]))
b <- exp(as.double(param[7]))
m <- exp(as.double(param[8])) #Aqui caguei no pse.R, ai no Hipercubo.csv esta variável saiu com nome d, mas era pra ser m
dm <- as.double(param[9])

HabProp <- seq(5, 85, 5)
label <- 1:5
Replica <- 1:5
for(i in HabProp)
{
	for(j in label)
	{
		for(k in Replica)
		{
			scape <- readLandscape(i, j) #será implementado assim que eu tiver acesso às paisagens da Elaine
			land <- list(numb.cells = 512, cell.size = 0.03, bound.condition = 1, land.shape = 1, scape = scape )
			class(land) <- "landscape"
			o.c <- i*100 + j*10 + k #Identificador do output
			#Lembrar de ajustar o tempo pra chegar no equilibrio
			TWoLife(raio=R, N=100, AngVis=360, passo=s, move=delta, taxa.basal=lamb0, taxa.morte=mu0, incl.birth=b, incl.death=m, density.type=1, death.mat=dm, landscape = land, tempo=50, ini.config=1, out.code=o.c)

		}
	}
}
files <- list.files(pattern="out*", full.names=F, recursive=FALSE)
for(i in files)
{
	analises(i, 20)
}
setwd("..")

}

projetoFelipe()



#Acho que esta função não vai estar sendo necessária
#modelRun <- function (my.pars) {
#    return(mapply(oneRun, my.pars[,1], my.pars[,2], my.pars[,3], my.pars[,4]))
#}

#file.create("iteration.txt")
#cat("run 1", file = "iteration.txt", append = TRUE, sep = "\n")

#TWoLife <- function (raio=0.1, N=80, AngVis=360, passo=5, move=0.5, taxa.basal=0.6, taxa.morte=0.1, incl.birth=0.5/0.01, incl.death=0, density.type=0, death.mat=7, landscape, tempo=20, ini.config=0, out.code=1)
#Nao entraram
#N
#AngVis
#ini.config (lembrar de igualar a 1)
#density.type (lembrar de igualar a 1)
#landscape (how the hell eu vou enfiar isso no mapply) ideia: fazer vetor de paisagens
#out.code (como enfiar um contador no mapply?)
#nboot e maxIt não sei escolher

#analises("output-00234.txt", npop)

#TODO:
#ler paisagens
#colocar correlações OK
#arrumar analises OK
#definir valores das variáveis constantes (e opções tipo nboot e maxIt)
