#ifndef INDIVIDUO_H
#define INDIVIDUO_H
#include <vector>

using namespace std;

/** \brief A classe individuo representa um agente da simulação.
 *
 *  Esta classe contém informações referentes a cada indivíduo, incluindo sua localização e estado, uma lista contendo ponteiros para os vizinhos próximos, tamanho do passo de movimentação, etc. Esta classe NÃO contém métodos que sejam de responsabilidade do ambiente, como o método de atualizar os vizinhos */
class individuo
{
private:
    //propriedades privadas
	/** Identificador único de cada indivíduo */
    const unsigned long id;
	/** Identificador máximo já usado */
    static unsigned long MAXID;
	/** Posição X da localização do indivíduo */
    double x;
	/** Posição Y da localização do indivíduo */
	double y;
	/** Espécie a qual o indivíduo pertence. Cada espécie será representada por um número único, TBI*/
    const int especie;
    double birth;//taxa de nascimento
    const double taxa_morte;//taxa intrínseca de morte quando no habitat e na ausencia de outros individuos na vizinhanca local
    const double move; //taxa de movimentacao
    const double passo;//tamanho do passo
    double orientacao;//orientacao
    const double ang_visada;//angulo de visada
    double tempo_evento;//tempo sorteado para a proxima acao
    double densi_max;//densidade maxima suportada pelo individuo
    double raio;//raio de percepcao da densidade
    vector<individuo*> lisViz;//vetor de vizinhanca    ----//CAMADA PERCEPTIVA
    const double taxa_basal;//taxa intrínseca de nascimento quando no habitat e na ausencia de outros individuos na vizinhanca local
    int tipo_habitat;//CAMADA PERCEPTIVA
    const int semente;//semente para gerar os numeros aleatorios
	
	double death;// taxa de morte
	double birth_death_eq; // taxa de natalidade e mortalidade quando N = K ou dens = dens_max
	const double incl_birth;//inclinação da curva de denso-depedência para natalidade
	const double incl_death;//inclinação da curva de denso-depedência para mortalidade
	const double const_d_matrix;// constante que indica quantas vezes a mortalidade basal na matriz eh maior que no habitat. Pensar em como considerar diferentes mortalidades (constantes) em diferentes tipos de matriz
	const int dens_type; //forma como a densidade √© calculada (0 = GLOBAL, 1 = LOCAL)
	//metodos privados
    void sorteiaTempo(); //tempo para ocorrer um evento
    void sampleTimeEXP();
    void sampleTimeRW();
	
public:
	/** Construtor da classe individuo. Deve ser chamado pela paisagem para posicionar os 
	 * individuos que já estão na paisagem no início da simulação. */
    individuo(
			/** Posição X do indivíduo */
			double x, 
			/** Posição Y do indivíduo */
			double y, 
			/** Espécie a qual o indivíduo pertence. Cada espécie é representada por um número único. \ref TBI */
			const int especie, 
			/** Taxa de crescimento intrínseco quando no habitat e na ausência de outros indivíduos na vizinhanca local */
			const double taxa_morte,
			/** Ângulo de orientação, veja \ref individuo::anda */
			double orientacao, 
			/** Ângulo de visada, veja \ref individuo::anda */
			const double angulo_visada,
			/** Tamanho do passo dado a cada movimentação, veja \ref individuo::anda */
            const double passo,
			/** Taxa de movimentação */
			const double move,
			/** Raio de percepção da densidade */
			const double raio, 
			/** Taxa de crescimento intrínseco quando no habitat e na ausência de outros indivíduos na vizinhanca local */
            const double taxa_basal,
			/** Semente para gera√ß√£o de n√∫meros aleat√≥rios */
			const int semente,// ??????????????
			const double incl_birth,
			const double incl_death,
			const double death_mat,
			const int dens_type);
	/** Construtor de c√≥pia, usado para gerar novos indiv√≠duos por reprodu√ß√£o ASSEXUADA.
	 *  Todos os dados do individuo pai ser√£o copiados, com a exce√ß√£o de:
	 *  - id (veja \ref individuo::get_id)
	 *  - vizinhos (veja \ref individuo::set_vizinhos)
	 *  - tempo para evento (veja \ref individuo::get_tempo) 
	 *  */
    individuo(/** Individuo pai */ const individuo& rhs); 

    /** Reinicia o contador de indivíduos **/
    static void reset_id() {MAXID = 0;}
	/** Retorna o identificador único deste indivíduo */
	const unsigned long get_id() const {return this->id;}
	/** Atualiza a lista de vizinhos deste indivíduo. Deve ser chamada a cada passo de tempo pela \ref paisagem */
    void set_vizinhos (/** Lista dos vizinhos */ const vector<individuo*> lis){this->lisViz=lis;}
    /**        **/
    const vector<individuo*> get_NBHood() const {return this->lisViz;}
    /**        **/
    void include_Neighbour(individuo * const agent){
	    for (unsigned int i = 0; i<this->lisViz.size(); i++)
		    if (agent == this->lisViz[i])
			    return; // already included
		    this->lisViz.push_back(agent);
    }
    void drop_Neighbour(individuo * const agent){
	    for (unsigned int i = 0; i<this->lisViz.size(); i++)
		    if (agent == this->lisViz[i])
			    this->lisViz.erase(this->lisViz.begin() + i);
    }
    /** Atualiza o tipo de hábitat no qual o indivíduo está. Deve ser chamada a cada passo de tempo pela \ref paisagem. */
    void set_habitat (const int tipo){this->tipo_habitat=tipo;}
	/** Atualiza a posi√ß√£o X do invid√≠duo */
    void set_x(/** Nova posi√ß√£o */double i){this->x =i;}
	/** Atualiza a posi√ß√£o Y do invid√≠duo */
    void set_y(/** Nova posi√ß√£o */double i){this->y =i;}
	/** Retorna a posi√ß√£o X do indiv√≠duo */
    inline const double get_x() const {return this->x;}
	/** Retorna a posição Y do indivíduo */
    inline const double get_y() const {return this->y;}
	/** Retorna o raio de percep√ßa√£o do indiv√≠duo */
    const double get_raio() const {return this->raio;}
	/** Retorna o tipo de densidade que afeta o indiv√≠duo */
    const int get_densType() const {return this->dens_type;}
	/** Returns the number of individuals inside the neighbourhood of the individual (it includes the focal individual) */
    const int NBHood_size() const {return this->lisViz.size()+1;}

    // outros metodos publicos
	/** Retorna o tempo sorteado para o próximo evento acontecer com este indivíduo.
	 * \sa 
	 * \ref individuo::update
	 * \ref paisagem::update */
    const double get_tempo(){return this->tempo_evento;}
	/** Sorteia uma das tr√™s a√ß√µes poss√≠veis e retorna a a√ß√£o sorteada para a paisagem
	 * (0 = morte, 1 = nascimento, 2 = movimento), baseado nas pr√≥prias taxas */
    int sorteia_acao();
	/** Atualiza a taxa de nascimento e/ou de morte baseado na densidade de indiv√≠duos dentro do raio de percep√ß√£o e sorteia o tempo
	 * do pr√≥ximo evento baseado nas taxas de nascimento, morte e movimenta√ß√£o 
	 * \sa \ref individuo::get_tempo */
    void updateEXPi();
    void updateLOGi(double dens);
    void updateRWi();
    void updateSKEXPi();
    void updateSKLOGi(double dens);
	/** Faz com que o indiv√≠duo ande um passo, do tamanho passo. A orienta√ß√£o na qual o indiv√≠duo vai andar √© a orienta√ß√£o atual
	 * (definida no construtor como orientacao) mais um √¢ngulo aleat√≥rio dentro do √¢ngulo de visada (angulo_visada). A defini√ß√£o de um 
	 * √¢ngulo de visada de 360 graus equivale a uma caminhada aleat√≥ria. */
    void anda(/** Passe aleatorio = true para for√ßar uma caminhada aleat√≥ria */ bool aleatorio = true);
};


#endif // INDIVIDUO_H
