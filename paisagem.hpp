#ifndef PAISAGEM_H
#define PAISAGEM_H
#include "individuo.hpp"
#include<vector>
#include<algorithm>
#include<cstdlib>
#include<cmath>

//Maximum landscape side lenght, used in fragmentation  routines TBI
#define dim 10000

using namespace std;

/*
 The paisagem/landscape class is responsible for generating the domain were the agents(individuals) reside
 This class is responsable for creating individuals (populating()), updating the world clock (atualiza_tempo()), and stablishes comunication between individuals(update()) */
class paisagem
{
private:

    // Private properties
    
    // Lenght of a square landscape side (if landscape_shape==1)
    const double tamanho;
    //Number of individuals at the start of the simulation
	const unsigned long N;
	//Vector for storing the individuals of the population
	vector <individuo*> popIndividuos;
	//Number of pixels in a square landscape side (if land_shape == 1)
    const int numb_cells;
    // Resolution Lenght of a square pixel side
    const double cell_size;
    //Shape of the landscape (0= circle, 1= square)
	const int landscape_shape;
	//The boundary condition type affects how individuals interact with the edges of the landscape (0= absortive, 1= periodical (pacman),2= reflexive)
	const int boundary_condition;
    // Matrix containing the evironmental values of the landscape pixels (0= matrix, 1= habitat)
	
    // Landscape moved to Public
    
	int patches[dim][dim];
	// number of non_matrix patches
	int numb_patches;
	// pointer for a vector storing the area of each patch
	double* patch_area;
	//The initial postion of individuals (0 = origin, 1 = random, 2 = normaly distributed with mean on origin)
	const int initialPos;

	//Private methods
    
    // Function responsible for calling the constructor of the individual/individuo class to create N individuals at the start of the simulation
	void populating(
					//Radius of density dependant influence
					const double raio,
					// Number of individuals at the start of the simulation
					const int N,
					// Angle used for orientation when dispersing
					const double angulo_visada,
					// The Lenght distance of a dispersal event
					const double passo,
					// The rate at which the individuals disperse
					const double move,
					  // The basal birth rate (The rate at which the individuals give birth on a habitat patch without neigbours)
					const double taxa_basal,
					// The basal death rate (The rate at which the individuals die on a habitat patch without neigbours)
					const double taxa_morte,
					// The slope of the birth density dependance function
					const double incl_b,
					// The slope of the death density dependance function
					const double incl_d,
					// Constant that indicates how many times higher the death rate should be on non-habitat pixels
					const double death_m,
					// Density type (0 = global, 1 = local/within a individual radius)
					const int dens_type,
                    double phenotype_mean,
                    double width_sd
					);

    /*
     Updates an individuals neighbourhood list (each of the individuals within a radius distance of the focal individual)
     Paran individuo * const ind - an object of the individual/individuo class
    */
    void atualiza_vizinhos(individuo * const ind) const;
    
    /*
     Updates the individual the environmental value of the pixel corresponent to its current coordinate (0= matrix, 1= habitat) \ref individuo::tipo_habitat
     Paran individuo * const ind - an object of the individual/individuo class
     */
    void atualiza_habitat(individuo * const ind) const;
    
     /*
     Updates the individual the patch identification value of the pixel corresponent to its current coordinate (0= matrix, 1= patch1, 2= patch2... n= patchn, n+1=pathn+1) ndividuo::patch_label
     Paran individuo * const ind - an object of the individual/individuo class
    */
    void atualiza_patch(individuo * const ind) const;
    
    /*
     Aplies the boundary condition after a dispersal event
     Paran individuo * const ind - an object of the individual/individuo class
    */
	bool apply_boundary(individuo * const ind);

public:

    // The world counter used for storing how much time has already passed
    double tempo_do_mundo;
    
    int landscape[dim][dim];//[row][col] temporarely substitute for a fixed "dim"
    // Matrix determining the pixels of a patch (0= matrix, 1= patch1, 2= patch2... n= patchn, n+1=pathn+1) Previously on Private

	//Public methods
	// Constructor of the landscape/paisagem class
    paisagem(
			//Radius of density dependant influence
			const double raio,
			// Number of individuals at the start of the simulation
			const int N,
			// Angle used for orientation when dispersing
			const double angulo_visada,
            // The Lenght distance of a dispersal event
			const double passo,
			// The rate at which the individuals disperse
			const double move,
             // The basal birth rate (The rate at which the individuals give birth on a habitat patch without neigbours)
			const double taxa_basal,
			// The basal death rate (The rate at which the individuals die on a habitat patch without neigbours)
			const double taxa_morte,
			// The slope of the birth density dependance function
			const double incl_b,
			// The slope of the death density dependance function
			const double incl_d,
            //Number of pixels in a square landscape side (if land_shape == 1)
             const int numb_cells,
			// Resolution Lenght of a square pixel side
			const double cell_size,
			//Shape of the landscape (0= circle, 1= square)
			const int land_shape,
			// Density type (0 = global, 1 = local/within a individual radius)
			const int density_type,
			// Constant that indicates how many times higher the death rate should be on non-habitat pixels
			const double death_mat,
			//The initial postion of individuals (0 = origin, 1 = random, 2 = normaly distributed with mean on origin)
			const int inipos,
			//The boundary condition type affects how individuals interact with the edges of the landscape (0= absortive, 1= periodical (pacman),2= reflexive)
			const int bound_condition,
			// Vector containing the evironmental values of the landscape pixels (0= matrix, 1= habitat)
			double scape[]
			);

	// Function that calls other functions of the individual/individuo class to update the vector of individuals of the landscapepaisagem object (atualiza_vizinhos,  atualiza_habitat, atualiza_patch, update())
    void update();
    
    // Function that uses the get_tempo function of the individual/individuo class to draft times for each individual within th landscape and returns the individual with the least amount of time required for the execution of the next action
	int sorteia_individuo();
    
    /*
     Function that recieves the selected individual and calls a function member of the individual/individuao class to randomly select one of the three possible actions for that individual to perform (0= death, 2= birth, 3= dispersal)
     Paran const int lower - The postion of the individual with the lowest drafted time
    */
    int sorteia_acao(const int lower);
    
	/*
     Function that executes the selected action. It also returns a positive boolean value if the individual disperses out of the landscape boundary
     Paran int acao - A value representing the selected action (0= death, 2= birth, 3= dispersal)
     Paran int lower - The postion of the individual with the lowest drafted time
    */
	bool realiza_acao(int acao, int lower);
    
    /*
     Function that updates the world clock, adding time required for the action execution by the selection individual
     const int lower - The postion of the individual with the lowest drafted time
     */
	void atualiza_tempo(const int lower);
    
    
     // Function that calls a function of the individual/individuo class to return the total number of individuals currently in the landscape
    const int conta_individuos() const;
    
    /*
     Function that returs an individual from within the landscape vector of individuals
     Paran int i - The position(i vector) of the individual to be returned
     */
    individuo* get_individuos(int i) const;
    
    // Function that returns the curent number of species in the landscape
    const int conta_especies() const;
    
    // Function that returns the lenght of square landscape side
    const double get_tamanho() const ;
    
    /* Function that computes and returns the distance between a pair of individuals
     Paran const individuo* a1 - First individual
     Paran const individuo* a2 - Second individual
	*/
	double calcDist(const individuo* a1, const individuo* a2) const;
    
    /* Function that computes and returns the density of individuals within a radius distance area from a focal individual
     Paran const individuo* a1 - Focal individual
    */
	double calcDensity(const individuo* ind1) const;

    /*
     Function that returs a negative boolean value if the inserted individual was present at the start of the simulation and a positive value is the individuals was born during the simulation run time (used for colouring original individuals diferently)
     Paran individuo * const ind - Individual to be tested
     */
    const bool nascido(individuo * const ind);
    
    /* Recursive function used to find and identify the witch pixels belong to each of the landscape patches
     Paran int x - The initial value of the x coordinate
     Paran int y - The initial value of the y coordinate
     Paran int current_label - The initial path identification number
    */
    void find_patches(int x, int y, int current_label);
    
    // Function that returns the number of fragments on thw landscape
    int get_numb_patches();
    
    /* Function that computes and returns the area of a fragment in the landscape
    Paran int i - The patch number
     */
    double get_patch_area(int i) const ;
    
    // Function that decides upon, and calls other functions to perform, the desired method of dispersion (Random walk or habitat selection). In the habitat selection case this function also samples x Points within the individuals "passo" radius and the indivuduals relative fitness on that location.
    // const int lower - The postion of the individual with the lowest drafted time
    double paisagem::walk(int lower);
};

#endif // PAISAGEM_H
