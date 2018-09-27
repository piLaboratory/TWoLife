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
    /** Taxa de natalidade final, após a consideração de todos os efeitos */
    double birth;
    /** Taxa intrínseca de mortalidade quando no habitat e na ausência de outros indivíduos na vizinhança local */
    const double taxa_morte;
    /** Taxa de movimentação */
    const double move; 
    /** Tamanho do passo dado a cada movimentação */
    const double passo;
    /** Ângulo de orientação sorteado para a movimentação */
    double orientacao;
    /** Ângulo de visada */
    const double ang_visada;
    /** Tempo sorteado para o próximo evento acontecer com este indivíduo */
    double tempo_evento;
    /** Densidade máxima suportada pelo indivíduo*/
    double densi_max;
    /** Raio de percepção da densidade */
    double raio;
    /** Vetor de vizinhanca */
    vector<individuo*> lisViz;//vetor de vizinhanca    ----//CAMADA PERCEPTIVA
    /** Taxa intrínseca de natalidade quando no habitat e na ausência de outros indivíduos na vizinhança local */
    const double taxa_basal;
    /** Identificdor do tipo de habitat do pixel correspondente à atual posição do indivíduo (0 = matriz; 1 = habitat) */
    int tipo_habitat;//CAMADA PERCEPTIVA
    /** Identificdor do fragmento correspondente à atual posição do indivíduo (0 = matriz; 1, 2, ... = fragmentos de habitat) */
    int patch_label;
    /** Semente para gerar os numeros aleatorios */
    const int semente;
	/** Taxa de mortalidade final, após a consideração de todos os efeitos */
	double death;
	/** Taxa de natalidade e mortalidade quando N = K ou dens = dens_max */
	double birth_death_eq;
	/** Inclinação da curva de denso-depedência para natalidade */
	const double incl_birth;
	/** Inclinação da curva de denso-depedência para mortalidade */
	const double incl_death;
	/** Constante que indica quantas vezes a mortalidade basal na matriz eh maior que no habitat */
	const double const_d_matrix;// Pensar em como considerar diferentes mortalidades (constantes) em diferentes tipos de matriz 
	/** Forma como a densidade é calculada (0 = GLOBAL, 1 = LOCAL)*/
	const int dens_type;
	//metodos privados
	/** Gera um número aleatório, segundo distribuição exponencial, correspondente ao tempo para a ocorrência do próximo evento, levando em consideração as taxas de mortalidade, natalidade e movimentação */
    void sorteiaTempo();
	
public:
	/** Construtor da classe individuo. Deve ser chamado pela paisagem para posicionar os 
	 * indivíduos que já estão na paisagem no início da simulação. */
    individuo(
			/** Posição X da localização do indivíduo */
			double x, 
			/** Posição Y da localização do indivíduo */
			double y, 
			/** Espécie a qual o indivíduo pertence. Cada espécie é representada por um número único. \ref TBI */
			const int especie, 
			/** Taxa intrínseca de mortalidade quando no habitat e na ausência de outros indivíduos na vizinhança local */
			const double taxa_morte,
			/** Ângulo de orientação sorteado para a movimentação, veja \ref individuo::anda */
			double orientacao, 
			/** Ângulo de visada, veja \ref individuo::anda */
			const double angulo_visada,
			/** Tamanho do passo dado a cada movimentação, veja \ref individuo::anda */
            const double passo,
			/** Taxa de movimentação */
			const double move,
			/** Raio de percepção da densidade */
			const double raio, 
			/** Taxa intrínseca de natalidade quando no habitat e na ausência de outros indivíduos na vizinhança local */
            const double taxa_basal,
			/** Semente para gerar os numeros aleatorios */
			const int semente,
			/** Inclinação da curva de denso-depedência para natalidade */
			const double incl_birth,
			/** Inclinação da curva de denso-depedência para mortalidade */
			const double incl_death,
			/** Constante que indica quantas vezes a mortalidade basal na matriz é maior que no habitat */
			const double death_mat,
			/** Forma como a densidade é calculada (0 = GLOBAL, 1 = LOCAL)*/
			const int dens_type);
	/** Construtor de cópia, usado para gerar novos indivíduos por reprodução assexuada.
	 *  Todos os dados do individuo pai serão copiados, com a exceção de:
	 *  - id (veja \ref individuo::get_id)
	 *  - vizinhos (veja \ref individuo::set_vizinhos)
	 *  - tempo para evento (veja \ref individuo::get_tempo) 
	 *  */
    individuo(/** Indivíduo pai */ const individuo& rhs); 

    /** Reinicia o contador de indivíduos */
    static void reset_id() {MAXID = 0;}
	/** Retorna o identificador único deste indivíduo */
	const unsigned long get_id() const {return this->id;}
	/** Atualiza a lista de vizinhos deste indivíduo. Deve ser chamada a cada passo de tempo pela \ref paisagem */
    void set_vizinhos (/** Lista dos vizinhos */ const vector<individuo*> lis){this->lisViz=lis;}
	/** Atualiza o tipo de hábitat no qual o indivíduo está. Deve ser chamada a cada passo de tempo pela \ref paisagem. */
    void set_habitat (const int tipo){this->tipo_habitat=tipo;}
    /** Atualiza o fragmento do indivíduo */
    void set_patch (const int label){this->patch_label=label;}
	/** Atualiza a posição X do indivíduo */
    void set_x(/** Nova posição */double i){this->x =i;}
	/** Atualiza a posição Y do indivíduo */
    void set_y(/** Nova posição */double i){this->y =i;}
	/** Retorna a posição X do indivíduo */
    inline const double get_x() const {return this->x;}
	/** Retorna a posição Y do indivíduo */
    inline const double get_y() const {return this->y;}
	/** Retorna o raio de percepção do indivíduo */
    const double get_raio() const {return this->raio;}
	/** Retorna o tipo de densidade que afeta o indivíduo */
    const int get_densType() const {return this->dens_type;}
	/** Returns the number of individuals inside the neighbourhood of the individual (it includes the focal individual) */
    const int NBHood_size() const {return this->lisViz.size()+1;}
	/** Retorna o fragmento no qual o indivíduo está */
    const int get_patch() const {return this->patch_label;}

    // outros metodos publicos
	/** Retorna o tempo sorteado para o próximo evento acontecer com este indivíduo.
	 * \sa 
	 * \ref individuo::update
	 * \ref paisagem::update */
    const double get_tempo(){return this->tempo_evento;}
	/** Sorteia uma das três ações possíveis de acordo com suas respectivas taxas e retorna a ação sorteada para a paisagem
	 * (0 = morte, 1 = nascimento, 2 = movimentação) */
    int sorteia_acao();
	/** Atualiza a taxa de natalidade e/ou de mortalidade baseado na densidade de indivíduos dentro do raio de percepção e sorteia o tempo
	 * do próximo evento baseado nas taxas de natalidade, de mortalidade e de movimentação 
	 * \sa \ref individuo::get_tempo */
    void update(double dens);
	/** Faz com que o indivíduo ande um passo, de tamanho passo. A orientação na qual o indivíduo vai andar é a orientação atual
	 * (definida no construtor como orientacao) mais um ângulo aleatório dentro do ângulo de visada (angulo_visada). A definição de um 
	 * ângulo de visada de 360 graus equivale a uma caminhada aleatória */
    void anda(/** Passe aleatorio = true para forçar uma caminhada aleatória */ bool aleatorio = true);
};
//Falta mexer no doxygen dos construtores e da update

#endif // INDIVIDUO_H
