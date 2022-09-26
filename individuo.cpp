#include "individuo.hpp"
#include <cmath>
#include "paisagem.hpp"
#include <cstdlib>
#include <iostream>
#include <R.h>
#include <Rmath.h>

using namespace std;

// Initializes max ID
unsigned long individuo::MAXID = 0;

/** \brief Individual exception
 Displays warning messages for impossible values
 (\ref TBI). */
class individuo_exception: public exception
{
    virtual const char* what() const throw()
    {
        return "Aconteceu uma exceção na classe indivíduo!";
    }
} myex;

// Class Constructor
individuo::individuo(double x,
                     double y,
                     int especie,
                     double taxa_morte,
                     double orientacao,
                     double angulo_visada,
                     double passo,
                     double taxa_move,
                     double raio,
                     double taxa_basal,
                     int semente, 
                     double incl_b,
                     double incl_d,
                     double death_mat,
                     double move_mat,
                     int dens_type,
                     vector <double> genotype,
                     vector <double> width,
                     int points):
// This section is similar to regular atribution and is run before brackets
id(++MAXID), // pega o próximo ID livre
x(x),
y(y),
especie(especie),
taxa_morte(taxa_morte),
orientacao(orientacao),
ang_visada(angulo_visada),
passo(passo),
taxa_move(taxa_move),
move(taxa_move),
raio(raio),
taxa_basal(taxa_basal),
semente(semente),
incl_birth(incl_b),
incl_death(incl_d),
const_d_matrix(death_mat),
const_m_matrix(move_mat),
dens_type(dens_type),
points(points)
{
    // Checks for impossible values
    if (taxa_morte < 0) throw myex;
    if (passo < 0) throw myex;
    if (taxa_move < 0) throw myex;
    if (raio < 0) throw myex;
    if (taxa_basal < 0) throw myex;
    
    if(incl_birth!=0 && incl_death!=0) // Checks if the birth and death rates are density dependant
    {
        this->densi_max = (taxa_basal-taxa_morte)/(incl_b+incl_d); // Sets the maximum density
        this->birth_death_eq = taxa_morte+incl_d*((taxa_basal-taxa_morte)/(incl_b+incl_d)); // Sets the point where rates are equal
    }
    
    this->genotype_mean = genotype; //Sets the individuals genotype value (or values)
    
    this->width_sd = width; //Sets the individuals width value (or values)
    
    this->rdn_noise = 0; //Sets random noises (mutation) to zero (for future implementation)
    
    for (int i=0; i<this->genotype_mean.size(); i++) { // passes by the genotype means of the individual
        
        this->env_optimum.push_back(this->rdn_noise+this->genotype_mean[i]); // Generates phenotype by adding a random value to the genotype value
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

/** Copy constructor, used for generating new individuals by assexual reproduction
 All the characteristics of the parent individual will be copyed, exept:
    - id (veja \ref individuo::get_id/ will be set to new value)
    -list of Neighbours (veja \ref individuo::set_vizinhos/ will be updated)
    - time until nex event (veja \ref individuo::get_tempo/ will be drafted)
 Paran const individuo& rhs - Parent individual */
individuo::individuo(const individuo& rhs):
id(++MAXID),
x(rhs.x),
y(rhs.y),
especie(rhs.especie),
taxa_morte(rhs.taxa_morte),
ang_visada(rhs.ang_visada),
passo(rhs.passo),
taxa_move(rhs.taxa_move),
move(rhs.taxa_move),
raio(rhs.raio),
taxa_basal(rhs.taxa_basal),
semente(rhs.semente),
incl_birth(rhs.incl_birth),
incl_death(rhs.incl_death),
const_d_matrix(rhs.const_d_matrix),
const_m_matrix(rhs.const_m_matrix),
dens_type(rhs.dens_type),
tipo_habitat(rhs.tipo_habitat),
birth_death_eq(rhs.birth_death_eq)

{
    this->orientacao = runif(0,360); // Sets an ramdom orientation
    
    if (rhs.genotype_mean.size()==1) // Checks if this is not a null simulation
    {
        this->rdn_noise= 0; //Sets random noises (recombination/ mutation) to zero (for future implementation)
        //this->genotype_mean.push_back(rhs.genotype_mean[0] + runif(-1.0,1.0));
        this->genotype_mean.push_back(rhs.genotype_mean[0]+rdn_noise); //Sets the individuals genotype value
        this->width_sd= rhs.width_sd; //Sets the individuals width value (or values)
        //this->rdn_noise= 0; //Sets random noises (recombination/ mutation) to zero (for future implementation)
        this->env_optimum.push_back(this->rdn_noise+this->genotype_mean[0]); // Generates the genotype of the offspring by adding a random value to its parents genotype value
        
    }
    else // Checks if this is a null simulation
    {
        
        this->rdn_noise= 0;
        
        for (int i=0; i<rhs.genotype_mean.size(); i++) { // passes by the genotype means of the individual
            
            this->genotype_mean.push_back(this->rdn_noise+rhs.genotype_mean[i]); // Generates the genotype of the offspring by adding a random value to its parents genotype value
        }
        
        //this->genotype_mean= rhs.genotype_mean;
        this->width_sd = rhs.width_sd; //Sets the individuals width values
        this->rdn_noise = 0;  //Sets random noises (recombination/ mutation) to zero (for future implementation)
        
        for (int i=0; i<this->genotype_mean.size(); i++) { // passes by the genotype means of the individual
            
            this->env_optimum.push_back(this->rdn_noise+this->genotype_mean[i]); // Generates phenotype by ading a random value to the genotype values
        }
        
    }
    this->points=rhs.points; //Sets the points equal to the parents
}

/** \brief Método de atualização dos indivíduos
 * Esta função é a camada de atualização dos indivíduos
 * A cada execução deste método:
 * - É atualizada a taxa de nascimento de acordo com a densidade de vizinhos no raio de vizinhanca obtido com paisagem::atualiza_vizinhos
 * - Sorteia o tempo de acordo com as novas taxas
 * As taxas de morte e movimentação no momento fixas. Mas tambem serão funções da densidade de vizinhos (\ref TBI).
 */

/** Function that updates both death and birth rates of individuals based on the number of individuals within the area of their density dependance radius
It them calls a function to draft the time needed for execution of that individual based on its birth, death and dispersal rates*/
void individuo::update(double dens)
{
    double densi = dens; // Density of  neighbour individuals (includes focal one)
    
    
    if (this->width_sd[0] == 0) { // Checks if individuals are completely specialists
        if(this->tipo_habitat==0) // Checks if the individual is currently on the matrix
        {
            this->birth = 0; // Sets birth rate to 0
            // ToDo: Implementar aqui modelo mais geral para mortalidade na matriz. Aqui a denso dependencia é igual à do habitat, só muda a mortalidade basal que é maior que no habitat.
            this->death = this->const_d_matrix*this->taxa_morte+this->incl_death*densi; // Sets the higer mortality rate atributed to nonhabita patches
            this->move = this->const_m_matrix*this->taxa_move; // Sets the higher/lower movement rates atributed to nonhabita patches
        }
        else
        {
            this->birth = this->taxa_basal-this->incl_birth*densi; // Computes the actual birth rate on habitat patch (that is influenced by the density of neighbours)
            this->death = this->taxa_morte+this->incl_death*densi; // Computes the actual death rate on habitat patch (that is influenced by the density of neighbours)
            this->move = this->taxa_move; // Sets the standard movement rates
        }
    }
    else{
        
        //this->move = this->taxa_move;
        
        this->birth = this->taxa_basal-this->incl_birth*densi; // Computes the actual birth rate on habitat patch (that is influenced by the density of neighbours)
        this->death = ((this->const_d_matrix*this->taxa_morte)-((dnorm_sum(this->tipo_habitat, this->env_optimum, this->width_sd)/dnorm(this->env_optimum[0],this->env_optimum[0],this-> width_sd[0]))*((this->const_d_matrix*this->taxa_morte)-this->taxa_morte))); // Computes the actual death rate on habitat patch (that is influenced by the suitability of its current habitat)
        
    }
    
    if(this->birth<0) //Checks if the birth rate is lower than possible
    {
        this->birth=0; // Sets to the lowest possible value to Birth
        
    }
    
    this->sorteiaTempo(); // Calls the function to draft the time needed to execute the next action
}

// Function that drafts one of the three possible actions (acounting for their respective taxes) and returns the drafted action to the landscape
void individuo::sorteiaTempo()
{
    this->tempo_evento = rexp(1.0/(this->move+this->birth+this->death)); //Sets the time to the individual
}

// Function that drafts one of the three possible actions (acounting for their respective taxes) and returns the drafted action to the landscape
int individuo::sorteia_acao()
{
    vector<double> probscum; //creates a vector for storing the probabilities
    double total = this->death+this->birth+this->move; // Sums all rates
    probscum.push_back(this->death/total); // stores the death event probability at the first position
    probscum.push_back(probscum[0]+(this->birth/total)); // stores the birth event probability at the secon position
    probscum.push_back(probscum[1]+(this->move/total)); // stores the dispersal event probability at the third position
    double evento; // temporary variable for storing  a random value
    evento = runif(0.0,1.0); // Samples between 0 and 1
    
    int decisao; // temporary variable for storing dracted event
    for(unsigned int i=0; i<probscum.size()-1; i++) // Goes through all events
    {
        if(probscum[i]>evento)// Finds the first action probability higher than the drafted number
        {
            decisao = i; //Stores the drated action
            return decisao; // Returns the drafted action to the landscape (0 = death, 1 = birth, 2 = dispersal)
            break;
        }
    }
    return probscum.size()-1; // returns sorted action
}

//Function that dislocates the XY corrdinates of an individual by a "passo" lenght distance and oriented by the angle the individuals is currently facing added of a random component (angulo_vizada) (if angulo_vizada== 360 -> Random Walk)
void individuo::anda()
{
    this->orientacao+= runif(-this->ang_visada/2.0, this->ang_visada/2.0);//random way point to draft any direction
    if(this->orientacao < 0) // Checks for impossible (negative)values
        this->orientacao +=360; // Transforms in absolute dgrees
    else if(this->orientacao >= 360) // Checks for impossible (big)values
        this->orientacao -= 360; // Fixes spilage
    double oriRad=this->orientacao*M_PI/180.0;// Tranforms to radians to calculate the new XY coordinates
    
    double dx= cos(oriRad)*this->passo; // Calculates the new x coordinate
    double dy= sin(oriRad)*this->passo; // Calculates the new y coordinate
    this->x+=dx; // Sets new x coordinate
    this->y+=dy; // Sets new y coordinate
}

//Function that returns the value of the probability density function for the normal distribution
double individuo::dnorm(double x ,double mean, double sd){
    
    return (1/(sd*sqrt(2*M_PI))*(exp(-1*(pow(x-mean, 2)/(2*pow(sd, 2))))));
    
}

//Function that returns the value  the sum of several probability density function for the normal distribution
double individuo::dnorm_sum( double x ,vector <double> mean, vector<double> sd){
    
    double probcumsum=0; // Temporary variable for storing the sum
    
    for (int i=0; i<mean.size(); i++) { // Passes by each normal distribution
        probcumsum= probcumsum+ dnorm(x, mean[i],sd[i]); // Adds it to the counter
        
    }
    return probcumsum/mean.size(); // Returns the mean probability density
    
}

//Function that selects from within a given set of coordinates based on the dispersing individuals preference ranks via a softmax function
void individuo::habitat_selection(double possibilitities[][3])
{
    double scores[this->points]; // Temporary vector for storing the scores
    double cumsum=0, choice=0, score=0; // Temporary variable for storing the drafting components
    int final_pos= 0; // Temporary variable for storing the chosen index
    
    for (int i=0; i<this->points; i++) { // Passes by the points
        
        scores[i] = exp(dnorm_sum(possibilitities[i][2], this->genotype_mean, this->width_sd)); // Atributes a score to the coordinate based on the probability density function of the migrating individual
        cumsum += scores[i]; // Sums that score to the total
    }
    
    choice= runif(0.0,1.0); // Drafts between 0 and 1
    
    for (int i=0; i<this->points; i++) { // Passes by the points
        
        score += (scores[i]/cumsum); // Computes the comparative value of the ranks
        
        if (score > choice) { // Checks the comparative value of the ranks has exceeded the drafted value
            
            final_pos= i; // Choses the index of the selected coordinate
            break;
        }
    }
    
    this->x = possibilitities[final_pos][0]; // Sets new x coordinate
    this->y = possibilitities[final_pos][1]; // Sets new y coordinate
}

