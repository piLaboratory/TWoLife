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
					 double death_mat, int dens_type, vector <double> genotype_mean, vector <double> width_sd):
// This is run before the function
	id(++MAXID), // Updates the MAXID and uses ut to construct the individual
	x(x), // Sets the desire/drafted/inputed x Coordinate
	y(y), // Sets the desire/drafted/inputed y Coordinate
	especie(especie), // Sets the desire/drafted/inputed species ID
	taxa_morte(taxa_morte), // Sets the imputed basal death rate
	move(move), // Sets the inputed dispersial rate
	passo(passo), // Sets the inputed dispersial lenght
	orientacao(orientacao), // Sets the desire/drafted/inputed angle that an individual is curently facing
	ang_visada(angulo_visada), // Sets the inputed angle used for orientation when dispersing
	raio(raio),//  Sets the inputed radius of density dependant influence
	taxa_basal(taxa_basal), // Sets the inputed basal birth rate
	semente(semente), //Sets the inputed random seed
	incl_birth(incl_b), //Sets the inputed birth density dependant inclination
	incl_death(incl_d), //Sets the inputed death density dependant inclination
	const_d_matrix(death_mat), //Sets the inputed Constant that indicates how many times higher the death rate should be on non-habitat pixels
	dens_type(dens_type), //Sets the inputed density type

{
    // Checks if there is an impossible parametervalue
	if (taxa_morte < 0) throw myex;
	if (passo < 0) throw myex;
	if (move < 0) throw myex;
	if (raio < 0) throw myex;
	if (taxa_basal < 0) throw myex;	
	
    
	if(incl_birth!=0 && incl_death!=0) // Checks if the birth and death rates are density dependant
	{
		this->densi_max = (taxa_basal-taxa_morte)/(incl_b+incl_d); // Sets the maximum density
		this->birth_death_eq = taxa_morte+incl_d*((taxa_basal-taxa_morte)/(incl_b+incl_d)); // Sets the point where rates are equal
	}
    
    this->genotype_mean= genotype_mean;
    this->width_sd= width_sd;
    this->rdn_noise= runif(-1.0,1.0);
    
    for (int i=0; i<this->genotype_mean.size(); i++) {
        
        this->env_optimum.push_back(this->rdn_noise+this->genotype_mean[i]);
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

/* Copy constructor, used for generating new individuals by assexual reproduction
 All the characteristics of the parent individual will be copyed, exept:
    - id (veja \ref individuo::get_id/ will be set to new value)
    -list of Neighbours (veja \ref individuo::set_vizinhos/ will be updated)
    - time until nex event (veja \ref individuo::get_tempo/ will be drafted)
 Paran const individuo& rhs - Parent individual
*/
individuo::individuo(const individuo& rhs)
    :id(++MAXID), // Updates the MAXID and uses ut to construct the individual
      x(rhs.x), // Sets x coordinate equal to the parent
      y(rhs.y), // Sets y coordinate equal to the parent
      especie(rhs.especie), // Sets species ID equal to the parent
      taxa_morte(rhs.taxa_morte), // Sets basal death equal to the parent
      move(rhs.move), // Sets dispersal rate equal to the parent
      passo(rhs.passo), // Sets dispersal distance lenght equal to the parent
      ang_visada(rhs.ang_visada), // Sets the angle used for orientation when dispersing equal to the parent
      raio(rhs.raio), // Sets the radius of density dependant influence equal to the parent
      taxa_basal(rhs.taxa_basal), // Sets basal birth rate equal to the parent
      tipo_habitat(rhs.tipo_habitat), // Sets the current habitat type equal to the parent
      semente(rhs.semente), // Sets  the random seed equal to the parent
	  incl_birth(rhs.incl_birth), // Sets the birth density dependant inclination equal to the parent
	  incl_death(rhs.incl_death), // Sets death density dependant inclination equal to the parent
	  const_d_matrix(rhs.const_d_matrix), // Sets the inputed Constant that indicates how many times higher the death rate should be on non-habitat pixelsequal to the parent
	  dens_type(rhs.dens_type), // Sets the inputed density type equal to the parent
	  birth_death_eq(rhs.birth_death_eq) // Sets  birth death equilibrium point equal to the parent

{
    if (rhs.genotype_mean.size()==1)
    {
        
        this->genotype_mean.push_back(rhs.genotype_mean[0] + runif(-1.0,1.0));
        this->width_sd= rhs.width_sd;
        this->rdn_noise= runif(-1.0,1.0);
        this->env_optimum.push_back(this->rdn_noise+this->genotype_mean[0]);
        
    }
    else
    {
        
        this->rdn_noise= runif(-1.0,1.0);
        
        for (int i=0; i<genotype_mean.size(); i++) {
            
             this->genotype_mean.push_back(this->rdn_noise+genotype_mean[i]);
        }
        
        //this->genotype_mean= rhs.genotype_mean;
        this->width_sd= rhs.widt_sd;
        this->rdn_noise= runif(-1.0,1.0);
        
        for (int i=0; i<this->genotype_mean.size(); i++) {
            
            this->env_optimum.push_back(this->rdn_noise+this->genotype_mean[i]);
        }
        
    }
    
    
	
}

/** \brief Método de atualização dos indivíduos 
 * Esta função é a camada de atualização dos indivíduos
 * A cada execução deste método:
 * - É atualizada a taxa de nascimento de acordo com a densidade de vizinhos no raio de vizinhanca obtido com paisagem::atualiza_vizinhos
 * - Sorteia o tempo de acordo com as novas taxas
 * As taxas de morte e movimentação no momento fixas. Mas tambem serão funções da densidade de vizinhos (\ref TBI).
 */

/*
 Function that updates both death and birth rates of individuals based on the number of individuals within the area of their density dependance radius
It them calls a function to draft the time needed for execution of that individual based on its birth, death and dispersal rates
 */
void individuo::update(double dens)
{
  double densi = dens; // Density of  neighbour individuals (includes focal one)
  
    if (HEG == false) {
        if(this->tipo_habitat==0) // Checks if the individual is currently on the matrix
          {
              this->birth = 0; // Sets birth rate to 0
              
              // ToDo: Implementar aqui modelo mais geral para mortalidade na matriz. Aqui a denso dependencia é igual à do habitat, só muda a mortalidade basal que é maior que no habitat.
              this->death = this->const_d_matrix*this->taxa_morte+this->incl_death*densi; // Sets the higer mortality rate atributed to nonhabita patches
          }
        else
          {
              this->birth = this->taxa_basal-this->incl_birth*densi; // Computes the actual birth rate on habitat patch (that is influenced by the density of neighbours)
              this->death = this->taxa_morte+this->incl_death*densi; // Computes the actual death rate on habitat patch (that is influenced by the density of neighbours)
          }
    }
    else{
        
        this->birth = this->taxa_basal-this->incl_birth*densi; // Computes the actual birth rate on habitat patch (that is influenced by the density of neighbours)
        this->death = this->const_d_matrix-((dnorm_sum(this->tipo_habitat, this->genotype_mean, this->width_sd)/dnorm_sum(this->genotype_mean,this->genotype_mean,this-> width_sd))*(this->const_d_matrix-this->taxa_morte)); // Computes the actual death rate on habitat patch (that is influenced by the suitability of its current habitat)
        
    }
    
    if(this->birth<0) //Checks if the birth rate is lower than possible
    {
        this->birth=0; // Sets to the lowest possible value to Birth
        
    }
    
    
    
  this->sorteiaTempo(); // Calls the function to draft the time needed to execute the next action
}

// Generates random numbers, following a exponential distribution, corresponding to the time needed to execute the next action, taking in account the birth, death and dispersal rates.
void individuo::sorteiaTempo()
{
  this->tempo_evento = rexp(1.0/(this->move+this->birth+this->death)); //Sets the time to the individual
}

// Function that drafts one of the three possible actions (acounting for their respective taxes) and returns the drafted action to the landscape
int individuo::sorteia_acao()
{
    vector<double> probscum; //creates a vector for storing the probabilities
    double total = this->death+this->birth+this->move; // Sums all rates
    probscum.push_back(this->death/total); // stores the death probability at the first position
    probscum.push_back(probscum[0]+(this->birth/total)); // stores the birth probability at the secon position
    probscum.push_back(probscum[1]+(this->move/total));  // stores the dispersal probability at the third position
    
    double evento;
    evento = runif(0.0,1.0); // drafts a double between 0 and 1 to be used as treshhold value for selection the action

    int decisao;
	for(unsigned int i=0; i<probscum.size()-1; i++) // Goes through all events
    {
        if(probscum[i]>evento)// Finds the first action probability higher than the drafted number
        {
            decisao = i; // Stores the drated action
            return decisao;// Returns the drafted action to the landscape (0 = death, 1 = birth, 2 = dispersal)
			break;
		}
	}
	return probscum.size()-1; // returns 2 (dispersal) (DUVIDA)
}

//Function that dislocates the XY corrdinates of an individual by a "passo" lenght distance and oriented by the angle the individuals is currently facing added of a random component (angulo_vizada) (if angulo_vizada== 360 -> Random Walk)
void individuo::anda(bool aleatorio)
{
    
    //if (Habitat_Selection == False ) {
        
        if (aleatorio) // checks if random walk (defalt) is selected
        {
            this->orientacao = runif(-180.0,180.0);//random way point to draft any direction
        } else {
            this->orientacao+= runif(-ang_visada/2.0, ang_visada/2.0);//random way point within specified angle distance
        }
        double oriRad=this->orientacao*M_PI/180.0; // Tranforms to radians to calculate the new XY coordinates

        double dx= cos(oriRad)*this->passo; // Calculates the new x coordinate
        double dy= sin(oriRad)*this->passo; // Calculates the new y coordinate
        this->x+=dx; // Sets new x coordinate
        this->y+=dy; // Sets new y coordinate

    /*}
    else
    {
        
        double possibilitities[this->points][2];
        double scores[this->points];
        double cumsum=0, choice=0, dist;
        
        for (int i=0; i<this->points; i++) {
            
            choice=runif(0.0,360.0);
            dist=runif(0.0,this->passo);
            
            possibilitities[i][0]=this->x+cos(choice)*dist;
            possibilitities[i][1]=this->y+sin(choice)*dist;
            
            scores[i]->dnorm(landscape[possibilitities[i][0]][possibilitities[i][1]], this->genotype_mean, this->width_sd);
            cumsum+=exp(scores[i]);
        }
        
        choice= runif(0.0,1.0);
        
        for (int i=0; i<this->points; i++) {
            
            scores[i]= exp(scores[i])/cumsum;
            if (scores[i]>choice) {
                
                choice= i;
                break;
            }
            
        }
        
        this->x=possibilitities[i][0]; // Sets new x coordinate
        this->y=possibilitities[i][1]; // Sets new y coordinate
    }*/
}

void individuo::habitat_selection(double &possibilitities[][])
{
  double scores[this->points];
  double cumsum=0, choice=0, score=0;
  
  for (int i=0; i<this->points; i++) {
      
      scores[i] = exp(dnorm_sum(possibilitities[i][2], this->genotype_mean, this->width_sd));
      cumsum += scores[i];
  }
  
  choice= runif(0.0,1.0);
  
  for (int i=0; i<this->points; i++) {
      
      score += scores[i]/cumsum;
      
      if (score > choice) {
          
          choice= i;
          break;
      }
      
  }
  
  this->x=possibilitities[choice][0]; // Sets new x coordinate
  this->y=possibilitities[choice][1]; // Sets new y coordinate
}

// Function that resets the maximum ID
 static void individuo::reset_id()
{
    MAXID = 0; // Sets the MAXID to 0
    
}

// Function that returns the ID of an individual
const unsigned long individuo::get_id() const
{
    return this->id; // returns the ID of an individual
    
}

//Function that passes the imputed vector of neighbours of an individual (individuals within a radius distance) to the individuals
void individuo::set_vizinhos (const vector<individuo*> lis)
{
    this->lisViz=lis; //passes the imputed vector of neighbours of an individual (individuals within a radius distance) to the individuals
    
}

// Function that updates the habitat type of the pixeld the individual is currently on
void individuo::set_habitat (const double tipo)
{
    this->tipo_habitat=tipo; //Sets the habitat type of the individual to the current one
    
}

//Function that updates the fragment identfier of the pixel and individual is currently on
void individuo::set_patch (const int label)
{
    this->patch_label=label; //Sets the patch label of the individual to the current one
    
}

//Function that updates the X coordinate od the individual
void individuo::set_x(/** Nova posição */double i)
{
    this->x =i; //Sets the x position after dispersal
    
}

//Function that updates the y coordinate od the individual
 void individuo::set_y(/** Nova posição */double i)
{
    this->y =i; //Sets the y position after dispersal
    
}

// Function that returns the x coordinate of the individual
inline const double individuo::get_x() const
{
    return this->x; // Returns the x coordinate of the individual
    
}

// Function that returns the y coordinate of the individual
inline const double individuo::get_y() const
{
    return this->y; // Returns the y coordinate of the individual
    
}

// Function that returns the density dependance radius of the individual
const double individuo::get_raio() const
{
    return this->raio; // Returns the density dependance radius of the individual
    
}

// Function that returns the density type afecting the individual
const int individuo::get_densType() const
{
    return this->dens_type;// returns the density type afecting the individual
    
}

// Function that returns the number of individuals within the readius of the focal individuals for density calculations (it also includes the focal individual)
const int individuo::NBHood_size() const
{
    return this->lisViz.size()+1; // Returns the number of individuals within the readius of the focal individuals for density calculations (it also includes the focal individual)
    
}

// Function that retuns the fragment identfier of the pixel the individual is currently on
   const int individuo::get_patch() const
{
    return this->patch_label;// Retuns the fragment identfier of the pixel the individual is currently on
    
}

// Function that returns the drafted time for executing an action by the individual
const double individuo::get_tempo()
{
    return this->tempo_evento; // Returns the drafted time for executing an action by the individual
    
}

//Function that returns the value of the probability density function for the normal distribution
double individuo::dnorm(double x ,double mean=0, double sd=1){
    
    return (1/(sd*sqrt(2*M_PI))*(exp(-1*(pow(x-mean, 2)/(2*pow(sd, 2))))));
                
}

double individuo::dnorm_sum( double x ,vector <double> mean, vector<double> sd){
    
    double probcumsum=0;
    
    for (int i=0; i<mean.size(); i++) {
        probcumsum= probcumsum+ dnorm(x, mean[i],sd[i]);
        
    }
    return probcumsum/mean.size();
    
}
