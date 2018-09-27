library(doMC)
registerDoMC(4)
analises<-function(arquivo, npop){

dados = file(arquivo, "r")
dpaisagem = readLines(dados, n=9)
npatches = strtoi(unlist(strsplit(dpaisagem [4], " "))[3])

IDlist<-rep(0, npop)
patches <-matrix( rep(0, 3*npop), nrow = npop, ncol = 3 )
#colonizacao <- rep(0, npatches+1) #Conta sempre que um individuo chega num fragmento vazio (antes esse era o colonizacao1 e o próximo era o colonizacao2, porém ficou meio inútil e trabalhoso manter esta forma de contar colonizações, motivo pelo qual deixei de fora)
colonizacao <- rep(0, npatches+1) #Conta apenas quando um individuo chega num fragmento vazio vindo de um fragmento diferente, ou seja, nao conta quando o individuo vai para a matriz e volta
extincao <- rep(0, npatches+1) #Conta sempre que um fragmento fica vazio (Note que é de se esperar que as entradas de colonizacao e de extincao difiram apenas de +-1)
migracao1 <- rep(0, npatches+1) #Conta sempre que um individuo entra num fragmento, estando ele vazio ou nao. (Note que o individuo pode ficar entrando e saindo do fragmento e sempre será contado)
migracao2 <- rep(0, npatches+1) #Conta apenas quando o individuo vem de um fragmento diferente, ou seja, nao conta quando o individuo vai para a matriz e volta
emigracaoPaisagem <- rep(0, npatches+1) #Conta sempre que um individuo sai da paisagem
emigracao <- rep(0, npatches+1) #Conta sempre que um individuo sai do fragmento FALTA IMPLEMENTAR NO CODIGO
textinguiu <- rep(0, npatches+1) #Armazena o instante da última extinção
tcolonizou <- rep(0, npatches+1) #Armazena o instante da última colonização
tvazio <- list() #Armazena o tempo que o fragmento ficou vazio
tocupado <- list() #Armazena o tempo que o fragmento ficou cheio
tvazio[[npatches+2]] <- 0
tocupado[[npatches+2]] <- 0
tvazio[[npatches+2]] <- NULL
tocupado[[npatches+2]] <- NULL


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
#patchPopHistory <- c(0)
#colonizacaoHistory <- c(0)
#extincaoHistory <- c(0)
#migracao1History <- c(0)
#migracao2History <- c(0)
#emigracaoHistory <- c(0)
#emigracaoPaisagemHistory <- c(0)
#nascimentosHistory <- c(0)
#mortesHistory <- c(0)



tcolonizacao <- c(0)
textincao <- c(0)
tmigracao1 <- c(0)
tmigracao2 <- c(0)
temigracao <- c(0)
temigracaoPaisagem <- c(0)
tnascimentos <- c(0)
tmortes <- c(0)

tpatchPop <- c(0)

#Iniciar o patchPopHistory
lin = readLines(dados, n=1)
while(lin != "EOF")
{
	lin<-strsplit(lin, " ")
	line <- unlist(lin)
	acao <-strtoi(line[2])
	p <- strtoi(line[4])+1 #lembrando que no R os indices não começam do 0
	ID <- strtoi(line[3])
	t <- as.double(line[1])
	if(acao == 0)
	{
		mortes[p] <- mortes[p] + 1
		#mortesHistory <- c(mortesHistory, p)
		#tmortes <- c(tmortes, t)
		patchPop[p] <- patchPop[p] - 1
		#patchPopHistory <- c(patchPopHistory, p) 
		#tpatchPop <- c(tpatchPop, t)
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
		{
			extincao[p] <- extincao[p] + 1
			#extincaoHistory <- c(extincaoHistory, p)
			#textincao <- c(textincao, t)
			if(tcolonizou[p] > 0)
			{
				tocupado[[p]] <- c(tocupado[[p]], t - tcolonizou[p])
			}
			textinguiu[p] <- t
		}

	}

	if(acao == 1)
	{
		maxID<-maxID+1
		patchPop[p] <- patchPop[p] + 1
		#patchPopHistory <- c(patchPopHistory, p)
		#tpatchPop <- c(tpatchPop, t)
		nascimentos[p] <- nascimentos[p] + 1
		#nascimentosHistory <- c(nascimentosHistory, p)
		#tnascimentos <- c(tnascimentos, t)
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

			emigracao[patches[ind, 2]] <- emigracao[patches[ind, 2]] + 1
			#emigracaoHistory <- c(emigracaoHistory, patches[ind,2])
			#temigracao <- c(temigracao, t)

			patchPop[patches[ind,1]] <- patchPop[patches[ind,1]] + 1
			patchPop[patches[ind,2]] <- patchPop[patches[ind,2]] - 1
			#patchPopHistory <- c(patchPopHistory, patchPop)
			#tpatchPop <- c(tpatchPop, t)

			if(patchPop[patches[ind, 2]]==0)
			{
				extincao[patches[ind, 2]] <- extincao[patches[ind, 2]] + 1
				#extincaoHistory <- c(extincaoHistory, patches[ind, 2])
				#textincao <- c(textincao, t)
				if(tcolonizou[patches[ind, 2]] > 0)
				{
					tocupado[[patches[ind, 2]]] <- c(tocupado[[patches[ind, 2]]], t - tcolonizou[patches[ind, 2]])
				}
				textinguiu[patches[ind, 2]] <- t
			}

			

			migracao1[p] <- migracao1[p] + 1
			#migracao1History <- c(migracao1History, p)
			#tmigracao1 <- c(tmigracao1, t)
			
			if(patches[ind,2]>1 || patches[ind,2]!=1 || patches[ind,3]!=patches[ind,1]  )
			{
				migracao2[p] <- migracao2[p] + 1
				#migracao2History <- c(migracao2History, p)
				#tmigracao2 <- c(tmigracao2, t)
				
				if(patchPop[p]==1)
				{
					colonizacao[p] <- colonizacao[p] + 1
					#tcolonizacao <- c(tcolonizacao, t)
					#colonizacaoHistory <- c(colonizacaoHistory, p)
					if(textinguiu [p]> 0)
					{
						tvazio[[p]] <- c(tvazio[[p]], t - textinguiu[p])
					}
					tcolonizou[p] <- t
				}
				
			}
		}
	}

	if(acao == 3)
	{
		emigracaoPaisagem[p] <- emigracaoPaisagem[p] + 1
		#emigracaoPaisagemHistory <- c(emigracaoPaisagemHistory, p)
		#temigracaoPaisagem <- c(temigracaoPaisagem, t)
		patchPop[p] <- patchPop[p] - 1
		#patchPopHistory <- matrix( c(patchPopHistory, patchPop), nrow = npatches+1, ncol = ncol(patchPopHistory)+1 )
		#tpatchPop <- c(tpatchPop, t)
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
		{
			extincao[p] <- extincao[p] + 1
			#extincaoHistory <- c(extincaoHistory, p)
			#textincao <- c(textincao, t)
			if(tcolonizou[p] > 0)
			{
				tocupado[[p]] <- c(tocupado[[p]], t - tcolonizou[p])
			}
			textinguiu[p] <- t
		}

	}

	lin = readLines(dados, n=1)
}

arqm<-paste0("analises/mortes-", arquivo)
arqn<-paste0("analises/nascimentos-", arquivo)
arqe<-paste0("analises/extincao-", arquivo)
arqc<-paste0("analises/colonizacao-", arquivo)
arqm1<-paste0("analises/migracao1-", arquivo)
arqm2<-paste0("analises/migracao2-", arquivo)
arqem<-paste0("analises/emigracao-", arquivo)
arqemP<-paste0("analises/emigracaoPaisagem-", arquivo)
arqab<-paste0("analises/abundancia-", arquivo)
#arqtvazio<-paste0("analises/tvazio-", arquivo)
#arqtocupado<-paste0("analises/tocupado-", arquivo)


file.create(arqm)
file.create(arqn)
file.create(arqe)
file.create(arqc)
file.create(arqm1)
file.create(arqm2)
file.create(arqem)
file.create(arqemP)
file.create(arqab)
file.create(arqtvazio)
file.create(arqtocupado)


#for(j in length(tpatchPop):1)
#	cat(paste( tpatchPop[j], paste(patchPopHistory[,j], collapse = " "), sep = " " ), file = arqab, append = TRUE, sep = "\n")
#for(j in length(tmortes):1)
#	cat(paste( tmortes[j], paste(mortesHistory[,j], collapse = " "), sep = " " ), file = arqm, append = TRUE, sep = "\n")
#for(j in length(tnascimentos):1)
#	cat(paste( tnascimentos[j], paste(nascimentosHistory[,j], collapse = " ")), file = arqn, append = TRUE, sep = "\n")
#for(j in length(textincao):1)
#	cat(paste( textincao[j], paste(extincaoHistory[,j], collapse = " ")), file = arqe, append = TRUE, sep = "\n")
#for(j in length(tcolonizacao):1)
#	cat(paste( tcolonizacao[j], paste(colonizacaoHistory[,j], collapse = " ")), file = arqc, append = TRUE, sep = "\n")
#for(j in length(tmigracao1):1)
#	cat(paste( tmigracao1[j], paste(migracao1History[,j], collapse = " ")), file = arqm1, append = TRUE, sep = "\n")
#for(j in length(tmigracao2):1)
#	cat(paste( tmigracao2[j], paste(migracao2History[,j], collapse = " ")), file = arqm2, append = TRUE, sep = "\n")
#for(j in length(temigracao):1)
#	cat(paste( temigracao[j], paste(emigracaoHistory[,j], collapse = " ")), file = arqem, append = TRUE, sep = "\n")
#for(j in length(temigracaoPaisagem):1)
#	cat(paste( temigracaoPaisagem[j], paste(emigracaoPaisagemHistory[,j], collapse = " ")), file = arqemP, append = TRUE, sep = "\n")

for(i in 0:npatches)
{
	cat(paste(i,patchPop[i+1]), file = arqab, append = TRUE, sep = "\n")
	cat(paste(i,mortes[i+1]), file = arqm, append = TRUE, sep = "\n")
	cat(paste(i,nascimentos[i+1]), file = arqn, append = TRUE, sep = "\n")
	cat(paste(i,extincao[i+1]), file = arqe, append = TRUE, sep = "\n")
	cat(paste(i,colonizacao[i+1]), file = arqc, append = TRUE, sep = "\n")
	cat(paste(i,migracao1[i+1]), file = arqm1, append = TRUE, sep = "\n")
	cat(paste(i,migracao2[i+1]), file = arqm2, append = TRUE, sep = "\n")
	cat(paste(i,emigracao[i+1]), file = arqem, append = TRUE, sep = "\n")
	cat(paste(i,emigracaoPaisagem[i+1]), file = arqem, append = TRUE, sep = "\n")
	
#	if(length(tvazio[[i+1]]) > 0)
#	{
#		cat(i, file = arqtvazio, append = TRUE, sep = "")
#		for(j in tvazio[[i+1]])
#		{
#			cat(paste(" ", j, sep =""), file = arqtvazio, append = TRUE, sep = "")
#		}
#		cat("", file = arqtvazio, append = TRUE, sep = "\n")
#	}

#	if(length(tocupado[[i+1]]) > 0)
#	{
#		cat(i, file = arqtocupado, append = TRUE, sep = "")
#		for(j in tocupado[[i+1]])
#		{
#			cat(paste(" ", j, sep =""), file = arqtocupado, append = TRUE, sep = "")
#		}
#		cat("", file = arqtocupado, append = TRUE, sep = "\n")
#	}
}

cat("EOF", file = arqm, append = TRUE, sep = "\n")
cat("EOF", file = arqn, append = TRUE, sep = "\n")
cat("EOF", file = arqe, append = TRUE, sep = "\n")
cat("EOF", file = arqc, append = TRUE, sep = "\n")
cat("EOF", file = arqm1, append = TRUE, sep = "\n")
cat("EOF", file = arqm2, append = TRUE, sep = "\n")
cat("EOF", file = arqem, append = TRUE, sep = "\n")
cat("EOF", file = arqemP, append = TRUE, sep = "\n")
cat("EOF", file = arqab, append = TRUE, sep = "\n")
#cat("EOF", file = arqtvazio, append = TRUE, sep = "\n")
#cat("EOF", file = arqtocupado, append = TRUE, sep = "\n")

close(dados)
}

for(j in 1:10)
{
	setwd(paste("LHS-", j, sep =""))
	files <- list.files(pattern="out*", full.names=F, recursive=FALSE)
	dir.create("output")
	dir.create("analises")
	foreach(i = files) %dopar%
	{
		analises(i, 300)
		file.rename(i, paste("output/", i, sep = ""))
	}
	setwd("..")
}
