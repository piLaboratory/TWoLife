#include "paisagem.hpp"
#include "individuo.hpp"
#include <R.h>
#include <Rmath.h>

paisagem::paisagem(double raio, int N, double angulo_visada, double passo, double move, double taxa_basal, 
				   double taxa_morte, double incl_b, double incl_d, int numb_cells, double cell_size, int land_shape, 
				   int density_type, double death_mat, int bound_condition, int scape[]):
	tamanho(numb_cells*cell_size),
	N(N),
	tempo_do_mundo(0),
	numb_cells(numb_cells),
	cell_size(cell_size),
	landscape_shape(land_shape),
	boundary_condition(bound_condition)
{	
		for(int i=0;i<this->numb_cells;i++)
		{
			for(int j=0;j<this->numb_cells;j++)
			{
				// transforma o vetor scape recebido do construtor em uma matriz
				this->landscape[i][j]=scape[j*numb_cells+i];
			}
		}
		
	// Calculo do raio dependendo do tipo de densidade
	if(density_type==0)
	{
		raio = this->tamanho/sqrt(M_PI);
	}
	/* Coloca os indivíduos na paisagem por meio da função populating() */	
	this->populating(raio,N,angulo_visada,passo,move,taxa_basal,taxa_morte,incl_b,incl_d,death_mat,density_type);
	
}
		
void paisagem::populating(double raio, int N, double angulo_visada, double passo, double move, double taxa_basal,
						  double taxa_morte, double incl_b, double incl_d, double death_m,
						  int dens_type)
{
	// Considerar diferentes possibilidades de posições iniciais. TBI.
	for(int i=0; i<this->N; i++)
    {
        this->popIndividuos.push_back(new individuo(
                                                    0,//posicao x
                                                    0,//posicao y
                                                    0,//especie
                                                    taxa_morte,//taxa de morte
                                                    runif(0,360),// orientacao
                                                    angulo_visada,//angulo de visada
                                                    passo,//tamanho do passo
                                                    move, //taxa de movimentacao
                                                    raio,//tamanho do raio
                                                    taxa_basal,// taxa máxima de nascimento
                                                    99, // semente de numero aleatorio
													incl_b,
													incl_d,
                                                    death_m,
													dens_type));
        // como o popAgentes eh um ponteiro de vetores, ao adicionar enderecos das variaveis, usamos os new. Dessa forma fica mais rapido
        //pois podemos acessar apenas o endereco e nao ficar guardando todos os valores
    }
}

void paisagem::update()
{
    if(this->popIndividuos.size()>0)
    {    
		// Este for loop pode ser paralelizado, pois o que acontece com cada individuo eh independente
		#pragma omp parallel for
        for(unsigned int i=0; i<this->popIndividuos.size(); i++)
        {
            this->atualiza_vizinhos(this->popIndividuos[i]);//atualiza os vizinhos
            this->atualiza_habitat(this->popIndividuos[i]);//retorna o tipo de habitat
        }
		// Este loop não é parelelizado, APESAR de ser independente, para garantir que as funcoes
		// aleatorias sao chamadas sempre na mesma ordem (garante reprodutibilidade)
        for(unsigned int i=0; i<this->popIndividuos.size(); i++)
        {
            this->popIndividuos[i]->update();   //e atualiza o individuo i da populacao
        }

        this->realiza_acao();//escolhe tempo, indica e faz
    }
}


//para quando tiver especies, a definir...
//int paisagem::conta_especies()
//{
//    {
//        vector<int> especie_individuos;
//        int contador=0;
//        spIndividuos.assign(this->popIndividuos.size(),0);
//        int espe;
//        for(int i=0; i<this->popIndividuos.size(); i++)
//           {
//            espe = this->get_individuo(i)->get_sp();
//            especie_individuos[espe]++;
//        }
//        for(int z=0; z<this->popIndividuos.size();z++)
//        {
//            if(especie_individuos[z]!=0)
//            contador++;
//        }
//    return contador;
//    }
//}

void paisagem::realiza_acao() //TODO : criar matriz de distancias como atributo do mundo e atualiza-la apenas quanto ao individuos afetado nesta funcao)
{
    int menor=0;
    double menor_tempo = this->popIndividuos[0]->get_tempo();

    for(unsigned int i=1; i<this->popIndividuos.size(); i++)
    {
        if(this->popIndividuos[i]->get_tempo()<menor_tempo)
        {
            menor = i;
            menor_tempo = this->popIndividuos[i]->get_tempo();
        }
    }

    this->tempo_do_mundo = this->tempo_do_mundo+menor_tempo;
    int acao = this->popIndividuos[menor]->sorteia_acao();

    switch(acao) //0 eh morte, 1 eh nascer, 2 eh andar
    {
    case 0:
        delete this->popIndividuos[menor];
        this->popIndividuos.erase(this->popIndividuos.begin()+menor);
        break;

    case 1:
        individuo* chosen;
        //Novo metodo para fazer copia do individuo:
        chosen = new individuo(*this->popIndividuos[menor]);
        this->popIndividuos.push_back(chosen);
        break;

    case 2: 
        this->popIndividuos[menor]->anda();
		this->apply_boundary(popIndividuos[menor]);
		break;
    }
}

// metodo para condicao de contorno, argumento é um ponteiro para um individuo
//TODO: conferir se a combinacao x , y da condicao esta gerando o efeito desejado
//TBI: condicao periodica do codigo antigo feito com Garcia. Verificar se estah correta
// (veja p. ex. um unico individuo apenas se movimentando)
void paisagem::apply_boundary(individuo * const ind) //const
{
	double rad = (double)ind->get_raio();
	switch(this->boundary_condition)
	{
			
		case 0:
		if(this->landscape_shape==0)
		{
			if(rad*rad < (double)ind->get_x()*(double)ind->get_x()+(double)ind->get_y()*(double)ind->get_y())
			{
				for(unsigned int i=0; i<popIndividuos.size();i++)
				{
					if(this->popIndividuos[i]->get_id()==(int)ind->get_id())
					{
						delete this->popIndividuos[i];
						this->popIndividuos.erase(this->popIndividuos.begin()+i);
					}
				}
			}
		}
		
		if(this->landscape_shape==1)
		{
			if((double)ind->get_x()>=this->numb_cells*this->cell_size/2 || //>= porque na paisagem quadrado as bordas mais distantes de 0 iniciariam um proximo pixel que estaria fora da paisagem. Ou teriamos que assumir que esses pixels mais extremos tenha uma área maior, o que daria um trabalho adicional para implementar uma situação irreal.
			   (double)ind->get_x()<(this->numb_cells*this->cell_size/2)*(-1)||
			   (double)ind->get_y()>this->numb_cells*this->cell_size/2 ||
			   (double)ind->get_y()<=(this->numb_cells*this->cell_size/-2))
			{
				for(unsigned int i=0; i<popIndividuos.size();i++)
				{
					if(this->popIndividuos[i]->get_id()==(int)ind->get_id())
					{
						delete this->popIndividuos[i];
						this->popIndividuos.erase(this->popIndividuos.begin()+i);
					}
				}
			}			   
		}		
		break;
	
		case 1:
		if(ind->get_x()<0)
			ind->set_x(this->tamanho+ind->get_x());
		if(ind->get_x()>this->tamanho)
			ind->set_x(ind->get_x()-this->tamanho);
		if(ind->get_y()<0)
			ind->set_y(this->tamanho+ind->get_y());
		if(ind->get_y()>this->tamanho)
			ind->set_y(ind->get_y()-this->tamanho);
		break;
	}
	
	/* TBI
	case 2: reflexiva
	*/
}


double paisagem::calcDist(const individuo* a1, const individuo* a2) const //Virou método da paisagem pois não estava conseguindo usar as propriedades privadas da paisagem nessa função (p. ex. this->boundary_condition). Tive que tirar o primeiro const.
{
	switch(this->boundary_condition)
	{
		case 0:
			return sqrt(((double)a1->get_x()-(double)a2->get_x())*((double)a1->get_x()-(double)a2->get_x())+((double)a1->get_y()-(double)a2->get_y())*((double)a1->get_y()-(double)a2->get_y()));
			break;
	
		case 1:
			double x1=a1->get_x();
			double x2=a2->get_x();
			double y1=a1->get_y();
			double y2=a2->get_y();
			double dx= x1 > x2 ? x1 - x2 : x2 - x1; // escolhe o menor entre x1-x2 e x2-x1
			double dy= y1 > y2 ? y1 - y2 : y2 - y1;
			dx = dx > this->tamanho - dx ? this->tamanho - dx : dx; // escolhe o menor entre tam-dx e dx
			dy = dy > this->tamanho - dy ? this->tamanho - dy : dy;
			return sqrt(dx*dx + dy*dy);
			break;
	}
}

/*
  Sempre adicione const aos argumentos de métodos quando o método não
  deve alterá-los. Previne vários erros e pode otimizar compilação
*/

void paisagem::atualiza_vizinhos(individuo * const ag1) const //acessando os vizinhos dos agentes
{  
	vector <individuo*> listViz;
	if(ag1->get_densType()==0) //dens_type poderia voltar como propriedade da paisagem. Facilitariam as coisas. Como muitas propriedades e métodos deste código, elas podem ser interpretadas das duas formas (como do individuo ou como da paisagem). O que está dando confusão é que estamos fazendo um IBM, mas para algumas situações estamos querendo simular dinâmicas cujas variáveis de interesse são propriedades populacionais e não do indivíduo. Se aceito, limar o método get_densType() do individuo.h.
	{		
		for(unsigned int j=0; j<popIndividuos.size(); j++)
		{
			individuo* ag2=this->popIndividuos[j];
			if(ag1==ag2) continue;			   
			listViz.push_back(ag2);
		}
	}
	else
	{
		double rad = (double)ag1->get_raio();
		for(unsigned int j=0; j<popIndividuos.size(); j++)
		{
			individuo* ag2=this->popIndividuos[j];
			if(ag1==ag2) continue;   
			double d=this->calcDist(ag1,ag2);
			if(d<=rad) {listViz.push_back(ag2);}		
		}
	}
	ag1->set_vizinhos(listViz);
		
}

void paisagem::atualiza_habitat(individuo * const ind) const
{
	// Tinha um IF com landscape_shape que eu removi. Não entendi como a paisagem ser circular 
	// interfere em ser habitat ou não: isso deve interferir na apply_boundary apenas, certo?
	// Also: Tinha uma inversão do y que eu também não entendi e removi
	// A.C. 10.07.13
	int hx,hy;
	hx= (double)ind->get_x()/this->cell_size+this->numb_cells/2;
	hy= (double)ind->get_y()/this->cell_size+this->numb_cells/2;
	ind->set_habitat(this->landscape[hx][hy]);
}

