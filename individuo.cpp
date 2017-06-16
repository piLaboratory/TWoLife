#include "individuo.hpp"
#include <cmath>
#include "paisagem.hpp"
#include <cstdlib>
#include <iostream>
#include <R.h>
#include <Rmath.h>

using namespace std;

unsigned long individuo::MAXID = 0;

/** \brief Exceção no indivíduo.
 *
 * Esta classe implementa uma exceção ocorrida na classe indivíduo. Atualmente, é usada para detectar valores impossíveis
 * passados no construtor, mas vai ser expandida para compreender mais casos (\ref TBI). */
class individuo_exception: public exception
{
	virtual const char* what() const throw()
	{
		return "Aconteceu uma exceção na classe indivíduo!";
	}
} myex;

individuo::individuo(double x, double y, int especie, double taxa_morte,
                     double orientacao, double angulo_visada,
                     double passo, double move, double raio,
                     double taxa_basal, int semente, //retirar int semente 
					 double incl_b, double incl_d,
					 double death_mat, int dens_type):
	id(++MAXID), // pega o próximo ID livre
	x(x),
	y(y),
	especie(especie),
	taxa_morte(taxa_morte),
	move(move),
	passo(passo),
	orientacao(orientacao),
	ang_visada(angulo_visada),
	raio(raio),
	taxa_basal(taxa_basal),
	semente(semente),
	incl_birth(incl_b),
	incl_death(incl_d),
	const_d_matrix(death_mat),
	dens_type(dens_type)
{
	if (taxa_morte < 0) throw myex;
	if (passo < 0) throw myex;
	if (move < 0) throw myex;
	if (raio < 0) throw myex;
	if (taxa_basal < 0) throw myex;	
	
	if(incl_birth!=0 && incl_death!=0)
	{
		this->densi_max = (taxa_basal-taxa_morte)/(incl_b+incl_d);
		this->birth_death_eq = taxa_morte+incl_d*((taxa_basal-taxa_morte)/(incl_b+incl_d));
	}
}


/** \brief função para reprodução assexuada de um indivíduo
 * 
 * Faz uma cópia do indivíduo pai. Recebe um ponteiro dereferenciado (nomeado rhs) e usa isto para
 * Criar novo individuo com os mesmos valores de atributo do pai, exceto
 *  - id (veja \ref individuo::get_id)
 *  - vizinhos (veja \ref individuo::set_vizinhos)
 *  - tempo para evento (veja \ref individuo::get_tempo) 
 * Usa notacao :atributo(valor) ao inves de atribuicão. 
 * Funcao chamada por paisagem::realiza_acao() quando a ação é um nascimentto
 * @param rhs ponteiro dereferenciado para o pai */
individuo::individuo(const individuo& rhs)
    :id(++MAXID),
      x(rhs.x),
      y(rhs.y),
      especie(rhs.especie),
      taxa_morte(rhs.taxa_morte),
      move(rhs.move),
      passo(rhs.passo),
      ang_visada(rhs.ang_visada),
      raio(rhs.raio),
      taxa_basal(rhs.taxa_basal),
      tipo_habitat(rhs.tipo_habitat),
      semente(rhs.semente),
	  incl_birth(rhs.incl_birth),
	  incl_death(rhs.incl_death),
	  const_d_matrix(rhs.const_d_matrix),
	  dens_type(rhs.dens_type),
	  birth_death_eq(rhs.birth_death_eq)
	  
{ //precisamos dessa chave e da que fecha ela?
	
}

/** \brief Método de atualização dos indivíduos 
 * Esta função é a camada de atualização dos indivíduos
 * A cada execução deste método:
 * - É atualizada a taxa de nascimento de acordo com a densidade de vizinhos no raio de vizinhanca obtido com paisagem::atualiza_vizinhos
 * - Sorteia o tempo de acordo com as novas taxas
 * As taxas de morte e movimentação no momento fixas. Mas tambem serão funções da densidade de vizinhos (\ref TBI).
 */
void individuo::update(double dens)
{
  double densi = dens; // densidade inclui n de vizinhos + o individuo
  if(this->tipo_habitat==0) 
	{
		this->birth = 0;
		// ToDo: Implementar aqui modelo mais geral para mortalidade na matriz. Aqui a denso dependencia é igual à do habitat, só muda a mortalidade basal que é maior que no habitat.
		this->death = this->const_d_matrix*this->taxa_morte+this->incl_death*densi; 
	}
  else 
	{
		this->birth = this->taxa_basal-this->incl_birth*densi;
		this->death = this->taxa_morte+this->incl_death*densi;
		if(this->birth<0){this->birth=0;}
	}
	
  this->sorteiaTempo();
}


void individuo::sorteiaTempo()
{
  this->tempo_evento = rexp(1.0/(this->move+this->birth+this->death));
}


int individuo::sorteia_acao()
{
    vector<double> probscum;
    double total = this->death+this->birth+this->move;
    probscum.push_back(this->death/total);
    probscum.push_back(probscum[0]+(this->birth/total)); 
    probscum.push_back(probscum[1]+(this->move/total));  
    double evento;
    evento = runif(0.0,1.0); 

    int decisao;
	for(unsigned int i=0; i<probscum.size()-1; i++)
    {
        if(probscum[i]>evento)//encontrar o primeiro maior valor
        {
            decisao = i;
            return decisao;//retorna a decisao tomada para o mundo (0 = morte, 1 = nascer, 2 = andar)
			break;
		}
	}
	return probscum.size()-1;
}

void individuo::anda(bool aleatorio)
{
    //a cada movimentacao ele muda de direcao; random walk dentro de um angulo de visao
	if (aleatorio) {
		this->orientacao = runif(-180.0,180.0);//random way point para sortear uma direcao qualquer
	} else {
		this->orientacao+= runif(-ang_visada/2.0, ang_visada/2.0);//random way point
	}
    double oriRad=this->orientacao*M_PI/180.0;//tranforma em radianos para poder calcular as distancias das posicoes x e y
    double dx= cos(oriRad)*this->passo;
    double dy= sin(oriRad)*this->passo;
    this->x+=dx;
    this->y+=dy;
}
