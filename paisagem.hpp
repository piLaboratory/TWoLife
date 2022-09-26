#ifndef PAISAGEM_H
#define PAISAGEM_H
#include "individuo.hpp"
#include<vector>
#include<algorithm>
#include<cstdlib>
#include<cmath>
///** Maximum landscape side lenght, used in fragmentation  routines TBI */
#define dim 10000 /*Em algum momento pode ser alterado para um argumento do
construtor. No momento não é prioritário. */
//Isto aqui não está aparecendo no doxygen
using namespace std;

/** \brief A classe paisagem implementa o mundo onde os agentes estão.
 *
 * The paisagem/landscape class is responsible for generating the domain were the agents(individuals) reside
 * This class is responsable for creating individuals, updating the world clock, and stablishes comunication between individuals
 * \sa \ref paisagem::update */

class paisagem
{
private:
    
    // Private properties
    /** Lenght of a square landscape side */
    const double tamanho;
    /** Number of individuals at the start of the simulation */
    const unsigned long N;
    /** Vector for storing the individuals of the population */
    vector <individuo*> popIndividuos;
    /** Number of pixels in a square landscape side */
    const int numb_cells;
    /** Resolution Lenght of a square pixel side */
    const double cell_size;
    /** Shape of the landscape (0= circle, 1= square)*/
    const int landscape_shape;
    /** The boundary condition type affects how individuals interact with the edges of the landscape (0= absortive, 1= periodical (pacman),2= reflexive) */
    const int boundary_condition;
    /** Matrix containing the evironmental values of the landscape pixels (0= matrix, 1= habitat)*/
    double landscape[dim][dim];
    /** Matriz coatining the fragment ID that each pixel resides (0 for  matriz; 1, 2, 3, ... for fragments) */
    int patches[dim][dim];
    /** number of non-matrix patches */
    int numb_patches;
    /** Pointer for a vector storing the area of each patch */
    double* patch_area;
    /** The initial postion of individuals (0 = origin, 1 = random, 2 = normaly distributed with mean on origin,  3 = inputed initial coordinates)*/
    const int initialPos;
    
    //Private methods
    
    // Function responsible for calling the constructor of the individual/individuo class to create N individuals at the start of the simulation
    
    void populating(
                    /** Radius of density dependant influence */
                    const double raio,
                    /** Number of individuals at the start of the simulation */
                    const int N,
                    /** Angle used for orientation when dispersing */
                    const double angulo_visada,
                    /** The Lenght distance of a dispersal event (random walk) or maximum dispersal radius (Selection step) */
                    const double passo,
                    /** The rate at which the individuals disperse */
                    const double taxa_move,
                    /** The basal birth rate (The rate at which the individuals give birth on a habitat patch without neigbours) */
                    const double taxa_basal,
                    /** The basal death rate (The rate at which the individuals die on a habitat patch without neigbours) */
                    const double taxa_morte,
                    /** The slope of the birth density dependance function */
                    const double incl_b,
                    /** The slope of the death density dependance function*/
                    const double incl_d,
                    /** Constant that indicates how many times higher the death rate should be on non-habitat pixels */
                    const double death_m,
                    /** Constant that indicates how many times lower the movement rate should be on non-habitat pixels */
                    const double move_m,
                    /** Density type (0 = global, 1 = local/within a individual radius) */
                    const int dens_type,
                    /** Vector containing the x coordinates initial individuals */
                    double initialPosX[],
                    /** Vector containing the y coordinates initial individuals */
                    double initialPosY[],
                    /** Vector containing the genotypical trait means of the initial individuals */
                    double genotype_means[],
                    /** Vector containing the standard deviations of the initial individuals */
                    double width_sds[],
                    /** Vector containing the  numer of points initial individuals sample when selecting habitats */
                    int points[],
                    /**  Booling value switching simuluations to neutral state (all individuals acting as an average individual) */
                    bool Null);
    
    /** Creates an individuals neighbourhood list (each of the individuals within a radius distance of the focal individual), or all of them (global)
     Paran individuo * const ind - an object of the individual/individuo class */
    void set_vizinhos(individuo * const ind) const;
    
    /** Updates an individuals neighbourhood list (each of the individuals within a radius distance of the focal individual), or all of them (global)
     Paran individuo * const ind - an object of the individual/individuo class
     Paran int acao  - an variable containing the  action performed by the individual  */
    void atualiza_vizinhos(int acao, int ind) const;
    
    /** Updates in the individual the environmental value of the pixel corresponent to its current coordinate (if binary: 0= matrix, 1= habitat)
     Paran individuo * const ind - an object of the individual/individuo class
     \ref individuo::tipo_habitat */
    void atualiza_habitat(individuo * const ind) const;
    
    /** Updates in the individual the patch identification value of the pixel corresponent to its current coordinate
     Paran individuo * const ind - an object of the individual/individuo class
     \ref individuo::patch_label */
    void atualiza_patch(individuo * const ind) const;
    
    /** Aplies the boundary condition after a dispersal event
     Paran individuo * const ind - an object of the individual/individuo class */
    bool apply_boundary(individuo * const ind);
    
    
public:
    
    /** The world counter used for storing how much time has already passed */
    double tempo_do_mundo;
    
    //Public methods
    /** Constructor of the landscape/paisagem class */
    paisagem(
             /** Radius of density dependant influence */
             const double raio,
             /** Number of individuals at the start of the simulation */
             const int N,
             /** Angle used for orientation when dispersing */
             const double angulo_visada,
             /** The Lenght distance of a dispersal event (random walk) or maximum dispersal radius (Selection step) */
             const double passo,
             /** The rate at which the individuals disperse */
             const double taxa_move,
             /** The basal birth rate (The rate at which the individuals give birth on a habitat patch without neigbours) */
             const double taxa_basal,
             /** The basal death rate (The rate at which the individuals die on a habitat patch without neigbours) */
             const double taxa_morte,
             /** The slope of the birth density dependance function */
             const double incl_b,
             /** The slope of the death density dependance functione */
             const double incl_d,
             /** Number of pixels in a square landscape side  */
             const int numb_cells,
             /** Resolution Lenght of a square pixel side*/
             const double cell_size,
             /** Shape of the landscape (0= circle, 1= square)*/
             const int land_shape,
             /** Density type (0 = global, 1 = local/within a individual radius) */
             const int density_type,
             /** Constant that indicates how many times higher the death rate should be on non-habitat pixels*/
             const double death_mat,
             /** Constant that indicates how many times lower the movement rate should be on non-habitat pixels */
             const double move_mat,
             /** The initial postion of individuals (0 = origin, 1 = random, 2 = normaly distributed with mean on origin)*/
             const int inipos,
             /** The boundary condition type affects how individuals interact with the edges of the landscape (0= absortive, 1= periodical (pacman),2= reflexive)*/
             const int bound_condition,
             /** Vector containing the evironmental values of the landscape pixels (0= matrix, 1= habitat) */
             double scape[],
             /** Vector containing the x coordinates initial individuals */
             double initialPosX[],
             /** Vector containing the y coordinates initial individuals */
             double initialPosY[],
             /** Vector containing the genotypical trait means of the initial individuals */
             double genotype_means[],
             /** Vector containing the standard deviations of the initial individuals */
             double width_sds[],
             /** Vector containing the  numer of points initial individuals sample when selecting habitats */
             int points[],
             /** Vector containing the  numer of points initial individuals sample when selecting habitats */
             bool Null
             );
    
    /**  Function that calls other functions of the individual/individuo class to update the vector of individuals of the landscape/paisagem object (atualiza_vizinhos,  atualiza_habitat, atualiza_patch, update()) 
     Paran int acao - A value representing the selected action (0= death, 1= birth, 2= dispersal)
     Paranint  ind - Index of the Individual performing the action
     */
    void update(int acao, int ind);
    
    /** Function that uses the get_tempo function of the individual/individuo class to draft times for each individual within the landscape and returns the individual with the least amount of time required for the execution of the next action */
    int sorteia_individuo();
    
    /** Function that recieves the selected individual and calls a function member of the individual/individuao class to randomly select one of the three possible actions for that individual to perform (0= death, 1= birth, 2= dispersal)
     Paran const int lower - The postion of the individual with the lowest drafted time */
    int sorteia_acao(const int lower){return this->popIndividuos[lower]->sorteia_acao();}
    
    /** Function that executes the selected action. It also returns a positive boolean value if the individual disperses out of the landscape boundary
     Paran int acao - A value representing the selected action (0= death, 1= birth, 2= dispersal)
     Paran int lower - The postion of the individual with the lowest drafted time */
    bool realiza_acao(int acao, int lower);
    
    /** Function that updates the world clock, adding time required for the action execution by the selection individual
     const int lower - The postion of the individual with the lowest drafted time */
    void atualiza_tempo(const int lower){this->tempo_do_mundo = this->tempo_do_mundo + this->popIndividuos[lower]->get_tempo();}
    
    /** Function that calls a function of the individual/individuo class to return the total number of individuals currently in the landscape */
    const int conta_individuos() const{return popIndividuos.size();}
    
    /** Function that returs an individual from within the landscape vector of individuals
     Paran int i - The position(i vector) of the individual to be returned */
    individuo* get_individuos(int i) const {return popIndividuos[i];}
    
    /** Function that returns the curent number of species in the landscape
     \ref TBI */
    const int conta_especies() const;
    
    /** Function that returns the lenght of square landscape side */
    const double get_tamanho() const {return this->tamanho;}
    
    /** Function that computes and returns the distance between a pair of individuals
     Paran const individuo* a1 - First individual
     Paran const individuo* a2 - Second individual */
    double calcDist(const individuo* a1, const individuo* a2) const;
    
    /** Function that computes and returns the density of individuals within a radius distance area from a focal individual
     Paran const individuo* a1 - Focal individual */
    double calcDensity(const individuo* ind1) const;
    
    /** Function that returs a negative boolean value if the inserted individual was present at the start of the simulation and a positive value is the individuals was born during the simulation run time (used for colouring original individuals diferently)
     Paran individuo * const ind - Individual to be tested */
    const bool nascido(individuo * const ind) const {return ind->get_id() > this->N;}
    
    /** Recursive function used to find and identify the witch pixels belong to each of the landscape patches
     Paran int x - The initial value of the x coordinate
     Paran int y - The initial value of the y coordinate
     Paran int current_label - The initial path identification number */
    void find_patches(int x, int y, int current_label);
    
    /** Function that returns the number of fragments on the landscape */
    int get_numb_patches(){return numb_patches;}
    
    /** Function that computes and returns the area of a fragment in the landscape
     Paran int i - The patch number */
    double get_patch_area(int i) const {return this->patch_area[i];}
    
    /** Function that decides upon, and calls other functions to perform, the desired method of dispersion (Random walk or habitat selection).
     In the habitat selection case this function also samples x Points within the individuals "passo" radius and the indivuduals relative fitness on that location.
     const int lower - The index of the individual with the lowest drafted time */
    bool walk(int lower);
    
};

#endif // PAISAGEM_H

