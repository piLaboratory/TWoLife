#ifndef PAISAGEM_H
#define PAISAGEM_H
#include "individuo.hpp"
#include<vector>
#include<algorithm>
#include<cstdlib>
#include<cmath>
///** Dimensão maxima do mapa, usada para rotinas de fragmentação TBI */
#define dim 10000 /*Em algum momento pode ser alterado para um argumento do
		    construtor. No momento não é prioritário. */
//Isto aqui não está aparecendo no doxygen
using namespace std;

/** \brief A classe paisagem implementa o mundo onde os agentes estão.
 *
 * Esta classe é responsável pela criação dos indivíduos, por manter o relógio do mundo, e por estabelecer a comunicação entre os indivíduos \sa \ref paisagem::update */
class paisagem
{
private:

    //propriedades privadas
    /** Tamanho do lado da paisagem */
    const double tamanho;
    /** Número de indivíduos no começo da simulação */
	const unsigned long N;
	/** Vetor dos indivíduos da populacao */
	vector <individuo*> popIndividuos;
	/** Número de pixels do lado da paisagem */
    const int numb_cells;
    /** Tamanho do lada de cada pixel da imagem usada para a classificação da paisagem */
    const double cell_size;
    /** Forma da paisagem  (0 = circular, 1 = quadrada)*/
	const int landscape_shape;
	/** Tipo de condição de contorno (0 = absortiva, 1 = periódica, 2 = reflexiva) */
	const int boundary_condition;
	/** Matriz de pixels da paisagem, assume os valores 0 quando o pixel está na matriz e 1 quando está no habitat */
	int landscape[dim][dim];//[linha][coluna] temporariamente substituido or valor fixo
	/** Matriz com a determinação do fragmento a que cada pixel pertence (0 para matriz; 1, 2, 3, ... para fragmentos */
	int patches[dim][dim];
	/** Número de fragmentos da paisagem, desconsiderando-se a matriz. */
	int numb_patches;
	/** Ponteiro para o vetor contendo a área de cada fragmento */
	double* patch_area;
	/** Posição dos indivíduos no início da simulação (0 = origem; 1 = aleatória com distribuição uniforme na paisagem; 2 = aleatória com distribuição normal na paisagem)*/
	const int initialPos;

	//metodos privados
	void populating(
					/** Raio no qual os indivíduos percebem vizinhos. Deve ser menor que o tamanho do mundo */
					const double raio,
					/** Número de indivíduos no começo da simulação */
					const int N,
					/** Ângulo de visada dos indivíduos da população */
					const double angulo_visada,
					/** Passo de caminhada dos indivíduos da população */
					const double passo,
					/** Taxa de movimentação dos individuos da população */
					const double move,
					/** Taxa de nascimento de um indivíduo no habitat favorável e sem vizinhos */
					const double taxa_basal,
					/** Taxa de mortalidade dos indivíduos */
					const double taxa_morte,
					/** Inclinação da curva de denso-depedência para natalidade */
					const double incl_b,
					/** Inclinação da curva de denso-depedência para mortalidade */
					const double incl_d,
					/** Constante que indica quantas vezes a mortalidade basal na matriz é maior que no habitat */
					const double death_m,
					/** Tipo de densidade (0 = GLOBAL, 1 = LOCAL) */
					const int dens_type
					);

	/** Atualiza a lista de vizinhos de um indivíduo */
    void atualiza_vizinhos(individuo * const ind) const;
    /** Informa o indivíduo o tipo de habitat (matriz ou habitat) correspondente à sua atual posição, atualizando \ref individuo::tipo_habitat */
    void atualiza_habitat(individuo * const ind) const;
    /** Informa o indivíduo o fragmento correspondente à sua atual posição, atualizando \ref individuo::patch_label */
    void atualiza_patch(individuo * const ind) const;
    /** Aplica a condição de contorno após a movimentação */
	bool apply_boundary(individuo * const ind); //const; // metodo para aplicação da condicao de contorno


public:
	//vector <individuo*> popIndividuos;

	/** Contador de quanto tempo já transcorreu no mundo simulado */
    double tempo_do_mundo;

	//metodos publicos
	/** Construtor da classe paisagem */
    paisagem(
			/** Raio no qual os indivíduos percebem vizinhos. Deve ser menor que o tamanho do mundo */
			const double raio,
			/** Número de indivíduos no começo da simulação */
			const int N,
			/** Ângulo de visada dos indivíduos */
			const double angulo_visada,
            /** Passo de caminhada dos indivíduos */
			const double passo,
			/** Taxa de movimentação dos indivíduos */
			const double move,
			/** Taxa de nascimento de um indivíduo no habitat favorável e sem vizinhos */
			const double taxa_basal,
			/** Taxa de morte dos indivíduos */
			const double taxa_morte,
			/** Inclinação da curva de denso-depedência para natalidade */
			const double incl_b,
			/** Inclinação da curva de denso-depedência para mortalidade */
			const double incl_d,
			/** Número de pixels do lado da paisagem */
			const int numb_cells,
			/** Tamanho do lada de cada pixel da imagem usada para a classificação da paisagem */
			const double cell_size,
			/** Forma da paisagem  (0 = circular, 1 = quadrada)*/
			const int land_shape,
			/** Tipo de densidade ("g" = GLOBAL, "l" = LOCAL) */
			const int density_type,
			/** Constante que indica quantas vezes a mortalidade basal na matriz é maior que no habitat */
			const double death_mat,
			/** Posição dos indivíduos no início da simulação (0 = origem; 1 = aleatória com distribuição uniforme na paisagem; 2 = aleatória com distribuição normal na paisagem)*/
			const int inipos,
			/** Condição de contorno (0 = absortiva, 1 = periódica, 2 = reflexiva)*/
			const int bound_condition,
			/** Vetor de cobertura de habitat na paisagem */
			int scape[]
			); //construtor

	/** Atualiza as  */
    void update();//atualizador
    /** Seleciona o individuo com menor tempo até o proximo evento para realizar uma ação */
	int sorteia_individuo();
	/** Após a seleção do indivíduo que realizará a ação, sorteia uma das três ações possíveis de acordo com suas respectivas taxas e retorna a ação sorteada para a paisagem
	 * (0 = morte, 1 = nascimento, 2 = movimentação) */
	int sorteia_acao(const int lower){return this->popIndividuos[lower]->sorteia_acao();}
	/** Realiza a ação sorteado do indivíduo selecionado. Além disso informa se o indivíduo sai da paisagem por emigração */
	bool realiza_acao(int acao, int lower);//vai pegar os tempos de cada individuo e informa qual foi o escolhido e manda ele fazer */
	/** Atualiza o contador de tempo, somando o tempo para o evento do indivíduo selecionado */
	void atualiza_tempo(const int lower){this->tempo_do_mundo = this->tempo_do_mundo + this->popIndividuos[lower]->get_tempo();}
	/** Retorna o número total de indivíduos na paisagem */
    const int conta_individuos() const{return popIndividuos.size();}
	/** Retorna um vetor contendo todos os indivíduos na paisagem */
    individuo* get_individuos(int i) const {return popIndividuos[i];}
	/** Retorna o número de espécies existentes na paisagem \ref TBI */
    const int conta_especies() const;
	/** Retorna o tamanho da paisagem (definido no construtor) */
    const double get_tamanho() const {return this->tamanho;}
	/** Calcula a distância entre dois indivíduos */
	double calcDist(const individuo* a1, const individuo* a2) const;
	/** Calcula a densidade de vizinhos de um indivíduo */
	double calcDensity(const individuo* ind1) const;

    /** Retorna false se o indivíduo estava no ambiente no passo de tempo 0, e true se ele nasceu durante a simulação.
	 * Usado para pintar os indivíduos nascidos de um cor diferente dos individuos originais */
    const bool nascido(individuo * const ind) const {return ind->get_id() > this->N;}
	/** Função recursiva utilizada para encontrar os freagmentos da paisagem */
    void find_patches(int x, int y, int current_label);
	/** Retorna o número de fragmentos da paisagem */
    int get_numb_patches(){return numb_patches;}
	/** Retorna a área de um dado fragmento da paisagem */
    double get_patch_area(int i) const {return this->patch_area[i];}

};

#endif // PAISAGEM_H
