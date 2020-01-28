#ifndef INDIVIDUO_H
#define INDIVIDUO_H
#include <vector>

using namespace std;
/*
 The individual class represents an agent of the simulation
 
 This class conatins informations pertinent to an individual, including it's location, state, a vector containing pointers to close neighbours, lenght distance of a dispersal event, rate of dispersal, etc. This class DOES NOT contain methods of landscape/paisagem responsability, mutch like the update of the neighbour vector of each individual.
 */
class individuo
{
private:
    
    // Private properties

	// Unique individual identifier
    const unsigned long id;
	// Current maximum individual identifier
    static unsigned long MAXID;
	// The X coordinate of the individual
    double x;
	// The Y coordinate of the individual
	double y;
    // The species identifier number of the individual (TBI)
    const int especie;
    // Final/current birth rate of the individual (all effects considered)
    double birth;
    // Final/current death rate of the individual (all effects considered)
    const double taxa_morte;
    // Dispersal rate of the individual
    const double move; 
    //The Lenght distance of a dispersal event
    const double passo;
    // Angle that an individual is curently facing
    double orientacao;
    // Angle used for orientation when dispersing
    const double ang_visada;
    // Drafted time require for the individual to execute an action
    double tempo_evento;
    // Maximum density the individual is capable of bearing
    double densi_max;
    //Radius of density dependant influence
    double raio;
    // Vector for storing the neightbours of the individual
    vector<individuo*> lisViz;//vetor de vizinhanca    ----//CAMADA PERCEPTIVA
    // The basal birth rate (The rate at which the individuals give birth on a habitat patch without neigbours)
    const double taxa_basal;
    // The basal death rate (The rate at which the individuals die on a habitat patch without neigbours)
    const double taxa_morte;
    // Identfier of the Type of habitat the individual is currently on (0 = matriz; 1 = habitat)
    double tipo_habitat;//CAMADA PERCEPTIVA
    // Identfier of the patch of habitat the individual is currently on  (0= matrix, 1= patch1, 2= patch2... n= patchn, n+1=pathn+1)
    int patch_label;
    // Seed for generating random numers for the individual
    const int semente;
    // Final/current death rate of the individual (all effects considered)
	double death;
	// Birth and death rates when they are at equilibrium (populational K)
	double birth_death_eq;
	// The slope of the birth density dependance function
	const double incl_birth;
	// The slope of the death density dependance function
	const double incl_death;
	// Constant that indicates how many times higher the death rate should be on non-habitat pixels
    const double const_d_matrix;// Pensar em como considerar diferentes mortalidades (constantes) em diferentes tipos de matriz
	// Density type (0 = global, 1 = local/within a individual radius)
	const int dens_type;
    
    // The genetical optimum environmental value for a individual (where its mortality rate is the lowest)
    const double env_optimum;
    // The phenotipical optimum environmental value for a individual (where its mortality rate is the lowest)
    const double phenotype_mean
    // The standard deviation of environmental usage by a individual, how generalist it is
    const double width_sd;
    // The randomly selected displacement value (This is used to dislocate the fitness distribution of the individual)
    const double rdn_noise;
    
    
	// Private Methods
    // Fucntion that generates random numbers, following a exponential distribution, corresponding to the time needed to execute the next action, taking in account the birth, death and dispersal rates.
    void sorteiaTempo();
	
public:
    /*
     Constructor of the individual/individuo class
     Must be called by the landscape/paisagem class to position individuals at the star of the simualtion
    */
    
    individuo(
			// The X coordinate of the individual
			double x,
			// The y coordinate of the individual
			double y, 
			// The species identifier number of the individual (TBI)
			const int especie, 
			// The basal death rate (The rate at which the individuals die on a habitat patch without neigbours)
			const double taxa_morte,
            // Angle that an individual is curently facing
			double orientacao, 
			// Angle used for orientation when dispersing
			const double angulo_visada,
			//The Lenght distance of a dispersal event
            const double passo,
			// Dispersal rate of the individuals
			const double move,
			// Radius of density dependant influence
			const double raio, 
			// The basal birth rate (The rate at which the individuals give birth on a habitat patch without neigbours)
            const double taxa_basal,
			// Seed for generating random numers for the individual
			const int semente,
			// The slope of the birth density dependance function
			const double incl_birth,
			// The slope of the death density dependance function
			const double incl_death,
			// Constant that indicates how many times higher the death rate should be on non-habitat pixels
			const double death_mat,
			// Density type (0 = global, 1 = local/within a individual radius)
			const int dens_type,
            double phenotype_mean,
            double width_sd);
    
    /* Copy constructor, used for generating new individuals by assexual reproduction
     All the characteristics of the parent individual will be copyed, exept:
        - id (veja \ref individuo::get_id/ will be set to new value)
        -list of Neighbours (veja \ref individuo::set_vizinhos/ will be updated)
        - time until nex event (veja \ref individuo::get_tempo/ will be drafted)
     Paran const individuo& rhs - Parent individual
    */
    individuo(/** Indivíduo pai */ const individuo& rhs); 

    // Function that sets the maximum id to 0
    static void reset_id();
    
    // Function that returns the ID of an individual
    const unsigned long get_id() const;
    
    /*
     Function that passes the imputed vector of neighbours of an individual (individuals within a radius distance) to the individuals.
     Must be called at each time step by the landscape
     Paran const vector<individuo*> lis - The list of individuals within aradius dinstace of the focal individual
    */
    void set_vizinhos (/** Lista dos vizinhos */ const vector<individuo*> lis);
    
    /*
     Function that updates the habitat type of the pixeld the individual is currently on
     Mut be called at each time step by the landscape (0= matrix 1= habitat)
     Paran const int tipo - Pixel adress on the landscape matrix
     */
    void set_habitat (const int tipo);
    
    /*
     Function that updates the fragment identfier of the pixel the individual is currently on
     Paran const int label - Pixel adress on the patch id matrix
     
     */
    void set_patch (const int label);
    
    /*
     Function that updates the X coordinate of the individual
     Paran double i - The new x Coordinate
    */
    void set_x(/** Nova posição */double i);
    
	/*
     Function that updates the y coordinate of the individual
     Paran double i - The new y Coordinate
    */
    void set_y(/** Nova posição */double i);
    
    // Function that returns the x coordinate of the individual
    inline const double get_x() const;
    
	// Function that returns the y coordinate of the individual
    inline const double get_y() const;
    
	// Function that returns the density dependance radius of the individual
    const double get_raio() const;
    
    // Function that returns the density type afecting the individual (0= global, 1=local)
    const int get_densType() const;
    
    // Function that returns the number of individuals within the readius of the focal individuals for density calculations (it also includes the focal individual)
    const int NBHood_size() const;
    
    // Function that retuns the fragment identfier of the pixel the individual is currently on 0= matrix, 1= patch1, 2= patch2... n= patchn, n+1=pathn+1)
    const int get_patch() const;

    // Other Public Methods
    
	/** Retorna o tempo sorteado para o próximo evento acontecer com este indivíduo.
	 * \sa 
	 * \ref individuo::update
	 * \ref paisagem::update */
    // Function that Returns the drafted time for executing an action by the individual based on its birth, death and dispersal rates
    const double get_tempo();
    
    // Function that drafts one of the three possible actions (0 = death, 1 = birth, 2 = dispersal) (acounting for their respective taxes) and returns the drafted action to the landscape
    int sorteia_acao();
    
      //\sa \ref individuo::get_tempo
    /*
     Function that updates both death and birth rates of individuals based on the number of individuals within the area of their density dependance radius
     It them calls a function to draft the time needed for execution of that individual based on its birth, death and dispersal rates
     Paran double dens - The number of individuals pondered by the area of the density dependance radius
     */
    void update(double dens);
    
    /* Function that dislocates the XY corrdinates of an individual by a "passo" lenght distance and oriented by the angle the individuals is currently facing added of a random component (angulo_vizada) (if angulo_vizada== 360 -> Random Walk)
     Paran bool aleatorio = true -
    */
    void anda(/** Passe aleatorio = true para forçar uma caminhada aleatória */ bool aleatorio = true);
};

/* Function that returns the value of the probability density function for the normal distribution
   Paran double x - Quantile of interest
   Paran double mean - Mean value of the distribution
   Paran double sd - Standad deviation of the distribution
  */
  double dnorm(double x ,double mean=0, double sd=1);

/* Function that dislocates the XY corrdinates of an individual due to the selection of an sampled location within its "passo" radius distance.
   Paran double possibilitities[][] - array containing the sampled possible migration locations (x: first collum, y: second collum, fitness score:third collum)
*/
  void habitat_selection(double &possibilitities[][]);

//Falta mexer no doxygen dos construtores e da update

#endif // INDIVIDUO_H
