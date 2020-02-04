#include "paisagem.hpp"
#include "individuo.hpp"
#include <R.h>
#include <Rmath.h>

// Class Constructor
paisagem::paisagem(double raio, int N, double angulo_visada, double passo, double move, double taxa_basal, double taxa_morte, double incl_b, double incl_d, int numb_cells, double cell_size, int land_shape, int density_type, double death_mat, int inipos, int bound_condition, double scape[], double genotype_mean[], double width_sd[], bool Null):
    // This section is similar to regular atribution and is run before brackets
	tamanho(numb_cells*cell_size),
	N(N),
	tempo_do_mundo(0),
	numb_cells(numb_cells),
	cell_size(cell_size),
	landscape_shape(land_shape),
	boundary_condition(bound_condition),
	initialPos(inipos)
{
		for(int i=0;i<this->numb_cells;i++) // Passes through the rows
		{
			for(int j=0;j<this->numb_cells;j++) // Passes through the columns
			{
				
				this->landscape[i][j]=scape[j*numb_cells+i]; // Re-shapes the scape vector into a matrix
			}
		}

	// Atributes values to the landscape pixels, determining the patch they belong to
	int component = 0;// Counter variable used for identyfing the patches
	for (int k = 0; k < this->numb_cells; ++k) // Passes through the rows
		for (int l = 0; l < this->numb_cells; ++l) // Passes through the columns
			if (!this->patches[k][l] && this->landscape[k][l]) find_patches(k, l, ++component); // Checks if each landscape pixel has already been assigned a patch identfication value AND if it is a non-matrix patch
	this->numb_patches = component; // Stores the number of habitat patches in the landscape

	this->patch_area = new double[this->numb_patches+1]; //Points to a new vector object used for storing the area of all habitats and the matrix

	for (unsigned int j = 0; j<numb_patches+1; j++) //Passes through al numb_patches
		this->patch_area[j] = 0;// Initializes the vector giving it 0 as a null value

	for(int i = 0; i < this->numb_cells; i++) //Passes through the rows
		for(int j = 0; j < this->numb_cells; j++) //Passes through the columns
			this->patch_area[patches[i][j]] += 1;// Indexes the patch_area vector by the patch identfication numbers of each pixel in the landscape so that it adds the number o pixels of each patch
	for (unsigned int j = 0; j<numb_patches+1; j++) //Passes through the patches
		this->patch_area[j] = this->patch_area[j]*this->cell_size*this->cell_size;// multiplies the number of pixels of each patch by the area of a pixel

	// Modification of the radius for the global density case
    if(density_type==0)// Checks if the density type is set to zero (global)
    {
        raio = this->tamanho/sqrt(M_PI);// Changes the perception radius to the equivalent of a radius of a circle with the same area of a square with the lenght of the landscape (REFERENCE - Acho que esta errado, mesmo no caso do circulo deveria ser o diametro para contemplar a maior distancia possivel)
    }
    
    // Inserts N individuals in the landscape through the populating() function
    this->populating(raio,N,angulo_visada,passo,move,taxa_basal,taxa_morte,incl_b,incl_d,death_mat,density_ty, genotype_mean,
                     width_sd, Null);
	
	for(unsigned int i=0; i<this->popIndividuos.size(); i++)// Passes through each individuals
	{
		this->atualiza_vizinhos(this->popIndividuos[i]); // Calls a function of the individual/individuo class to update the neighbours of a individual (the individuals within a radius distance of the focal individual)
		this->atualiza_habitat(this->popIndividuos[i]); // Calls a function of the individual/individuo class to update the habitat type the individual is currently on
		this->atualiza_patch(this->popIndividuos[i]); //Calls a function of the individual/individuo class to update the patch identfier of the patch he is currently on
        }


	for(unsigned int i=0; i<this->popIndividuos.size(); i++) // Passes through each individual
	{
		double dsty=this->calcDensity(popIndividuos[i]); // Calls a function of the individual/individuo class to determine the density of individuals within a radius distance area from a focal individual
		this->popIndividuos[i]->update(dsty);   // Calls a function of the individual/individuo class to update the characteristics of each individual
	}

}

// Function responsible for calling the constructor of the individual/individuo class to create N individuals at the start of the simulation
void paisagem::populating(double raio, int N, double angulo_visada, double passo, double move, double taxa_basal, double taxa_morte, double incl_b, double incl_d, double death_m, int dens_type, double genotype_mean[], double width_sd[], bool Null)
{
	individuo::reset_id(); // Restars the id counter of the individuals
    
    if (Null==FALSE) {
        vector<double> genotype(1), width(1);
    }
    else{
        
        for (int i=0; i<N; i++) {
            genotype.push_back(genotype_mean[i]);
            width.push_back(width_sd[][i]);
        }
        
    }
    

	if(this->initialPos==0) // Checks if the initialPos is set to 0 (origin)
	{
		for(int i=0; i<this->N; i++) // Goess through the amount of initial individuals selected
		{
            if (Null==FALSE) {
                genotype= genotype_mean[i];	
                width= width_sd[][i];
            }

            
			this->popIndividuos.push_back(new individuo(//As popAgents is a pointer of vectors, when adding variable addresses we use "new". This way should be faster, as we can acsess only the adress instead of keep storing the values
														0,// X Coordinate
														0,// Y Coordinate
														0,// Species ID
														taxa_morte,// The basal death rate (The rate at which the individuals die on a habitat patch without neigbours)
														runif(0,360),// Angle that an individual is curently facing
														angulo_visada,// Angle used for orientation when dispersing
														passo,// The Lenght distance of a dispersal event
														move,// The rate at which the individuals disperse
														raio,// Density dependance radius
														taxa_basal,//The basal birth rate (The rate at which the individuals give birth on a habitat patch without neigbours)
														99, // Randomn seed
														incl_b,//The slope of the birth density dependece function
														incl_d,//The slope of the death density dependece function
														death_m,// Constant that indicates how many times higher the death rate should be on non-habitat pixels
														dens_type,//Density type (0 = global, 1 = local/within a individual radius)
                                                        genotype,
                                                        width
                                                        ));
		}
	}
	if(this->initialPos==1) // Checks if the initialPos is set to 1 (Random initial positions)
	{
		for(int i=0; i<this->N; i++) // Goess through the amount of initial individuals selected
		{
            if (Null==FALSE) {
                genotype= genotype_mean[i];
                width= width_sd[][i];
            }
            
			this->popIndividuos.push_back(new individuo(
														runif(this->tamanho/(-2),this->tamanho/2), // X Coordinate
														runif(this->tamanho/(-2),this->tamanho/2), // Y Coordinate
														0,// Species ID
														taxa_morte,// The basal death rate (The rate at which the individuals die on a habitat patch without neigbours)
														runif(0,360),// Angle that an individual is curently facing
														angulo_visada,// Angle used for orientation when dispersing
														passo,// The Lenght distance of a dispersal event
														move,// The rate at which the individuals disperse
														raio,// Density dependance radius
														taxa_basal,//The basal birth rate (The rate at which the individuals give birth on a habitat patch without neigbours)
														99, // Randomn seed
														incl_b,//The slope of the birth density dependece function
														incl_d,//The slope of the death density dependece function
														death_m,// Constant that indicates how many times higher the death rate should be on non-habitat pixels
														dens_type,//Density type (0 = global, 1 = local/within a individual radius)
                                                        genotype_mean,
                                                        width_sd
                                                        ));
		}
    }
	if(this->initialPos==2) // Checks if the initialPos is set to 2 (Random initial positions with normal distribution) TBI: tornar os parametros da rnorm livres
	{
		for(int i=0; i<this->N; i++) // Goess through the amount of initial individuals selected
		{
            if (Null==FALSE) {
                genotype= genotype_mean[i];
                width= width_sd[][i];
            }
            
			this->popIndividuos.push_back(new individuo(
														rnorm(0,sqrt(move)*passo),// X Coordinate
														rnorm(0,sqrt(move)*passo),// Y Coordinate
														0,// Species ID
														taxa_morte,// The basal death rate (The rate at which the individuals die on a habitat patch without neigbours)
														runif(0,360),// Angle that an individual is curently facing
														angulo_visada,// Angle used for orientation when dispersing
														passo,// The Lenght distance of a dispersal event
														move,// The rate at which the individuals disperse
														raio,// Density dependance radius
														taxa_basal,//The basal birth rate (The rate at which the individuals give birth on a habitat patch without neigbours)
														99, // Randomn seed
														incl_b,//The slope of the birth density dependece function
														incl_d,//The slope of the death density dependece function
														death_m,// Constant that indicates how many times higher the death rate should be on non-habitat pixels
														dens_type,//Density type (0 = global, 1 = local/within a individual radius)
                                                        genotype_mean,
                                                        width_sd
                                                        ));
		}
	}
}

// Updater function
void paisagem::update()
{
    if(this->popIndividuos.size()>0) // checks if there ate any currently alive individuals
    {
	// This loop can be parallelized, as what happens with each individual is independent
	#ifdef PARALLEL
	#pragma omp parallel for
	#endif
        
        for(unsigned int i=0; i<this->popIndividuos.size(); i++) //Goess through all currently alive individuals
        {
            this->atualiza_vizinhos(this->popIndividuos[i]);// updates an individuals neighbourhood list (each of the individuals within a radius distance of the focal individual)
            this->atualiza_habitat(this->popIndividuos[i]);// updates the individual the environmental value of the pixel corresponent to its current coordinate (0= matrix, 1= habitat)
            this->atualiza_patch(this->popIndividuos[i]);// updates the individual the patch identification value of the pixel corresponent to its current coordinate (0= matrix, 1= patch1, 2= patch2... n= patchn, n+1=pathn+1)
        }
        // This loop is not parallelized, in spite of being independent, to garantee that function with random components are always called in the sma order (reproductbility)
        
        for(unsigned int i=0; i<this->popIndividuos.size(); i++)//Goess through all currently alive individuals
        {
			double dsty=this->calcDensity(popIndividuos[i]); // returns the dentity of individuals within a radius distance area from each individual
            this->popIndividuos[i]->update(dsty); // Updates the characteristics of each individual
        }

	}
}
// Function that drafts an individual
int paisagem::sorteia_individuo()
{
	// time for next event and simulation time update
	int menor=0; // Stores the first individual of the vector as the one with curently the least required time.
	double menor_tempo = this->popIndividuos[0]->get_tempo(); // Drafts a time from a exponential equation used as the needed time to execute a action

	for(unsigned int i=1; i<this->popIndividuos.size(); i++) //Goess through all individuals
	{
		if(this->popIndividuos[i]->get_tempo()<menor_tempo) // Drafts a time for each individual and checks if it is lower than the curent lower one
		{
			menor = i; // Stores the current smaller individual
			menor_tempo = this->popIndividuos[i]->get_tempo(); // Stores the current lower time
		}
	}
	return menor; // Returns the individual with the lowest time
}

// Function to dract an action for the drafted individual
int paisagem::sorteia_acao(const int lower)
{
    return this->popIndividuos[lower]->sorteia_acao(); // Returns the individual with the lowest time
}

bool paisagem::realiza_acao(int acao, int lower) //TODO : criar matriz de distancias como atributo do mundo e atualiza-la apenas quanto ao individuos afetado nesta funcao)
{
	bool emigrou = false; // Initializes a boolean variable used to record if an individual has left the boundaries of the landscape as false
    
    switch(acao) //(0= death, 2= birth, 3= dispersal)
    {
    case 0:
        delete this->popIndividuos[lower]; // Deletes the object in the chosen position
        this->popIndividuos.erase(this->popIndividuos.begin()+lower); // Erases the position that object occupied in the vector of individuals
        break;

    case 1:
            
        individuo* chosen; // Declares a new object of the individual class
        chosen = new individuo(*this->popIndividuos[lower]);// Calls a function to create a new individuals with similar characteristics to the parent one
        this->popIndividuos.push_back(chosen); // Stores the new individual at the postition of the landscape individuals vector
        break;

    case 2:
            
        this->walk(lower);
            
	emigrou = this->apply_boundary(popIndividuos[lower]); //Apply the boundary condition if the individual disperses out of the landscape boundary Changes the boolean value to pisitive if a individual went ou of the landscape boundaries
		break;
    }
	return emigrou; // Returns a positive boolean value if a individual went ou of the landscape boundaries
}

// Time updater function
void paisagem::atualiza_tempo(const int lower)
{
    this->tempo_do_mundo = this->tempo_do_mundo + this->popIndividuos[lower]->get_tempo(); // Adds the time required for the execution of the seleced action by the drafted individual to the world clock
    
}

// Individual counter function
 const int conta_individuos() const
{
    return popIndividuos.size();
}

// Funtion that returns a given individual
individuo* paisagem::get_individuos(int i) const
{
    return popIndividuos[i]; // Returns the selected individual
}

// Function that returns the lenght of a square landscape side
const double paisagem::get_tamanho() const
{
    return this->tamanho; // Returns the lenght of a square landscape side
}

// metodo para condicao de contorno, argumento é um ponteiro para um individuo
//TODO: conferir se a combinacao x , y da condicao esta gerando o efeito desejado
//TBI: condicao periodica do codigo antigo feito com Garcia. Verificar se estah correta FE:Verifiquei e tinha um problema
// (veja p. ex. um unico individuo apenas se movimentando)

// Funtion to apply the desired boundary conditions when  a dispersal event leaves an individual out of the landscape boundarys
bool paisagem::apply_boundary(individuo * const ind) //const
{
	bool emigrou = false;// Initializes a boolean variable used to record if an individual has left the boundaries of the landscape as false
	double rad = (double)ind->get_raio();// Obtains the radius of the individual that dispersed
    
	switch(this->boundary_condition) //(0= absortive, 1= periodical (pacman),2= reflexive)
	{

		case 0:
		if(this->landscape_shape==0) // Checks if the landscape shape is set to 0 (circle)
		{
			if(rad*rad < (double)ind->get_x()*(double)ind->get_x()+(double)ind->get_y()*(double)ind->get_y()) // Nao entendi o que esta sendo checado
			{
				for(unsigned int i=0; i<popIndividuos.size();i++) //Goess through all individuals
				{
					if(this->popIndividuos[i]->get_id()==(int)ind->get_id()) // Checks if the individual to be killed (erased) is the right one
					{
						delete this->popIndividuos[i]; // Deletes the object in the chosen position
						this->popIndividuos.erase(this->popIndividuos.begin()+i); // Erases the position that object occupied in the vector of individuals
					}
				}
				emigrou = true; // Sets the boolean value to positive simbolizing the individual emigrated
			}
		}
            

		if(this->landscape_shape==1) // Checks if the landscape shape is set to 1 (square)
		{
			if((double)ind->get_x()>=(this->numb_cells*this->cell_size/2) || //Checks if the x coordinate is higher than the maximum value OR
			   (double)ind->get_x()<((this->numb_cells*this->cell_size/2)*(-1))|| // if the x coordinate is lower than the minimum value OR
			   (double)ind->get_y()>(this->numb_cells*this->cell_size/2) || // if the y coordinate is higher than the maximum value Or
			   (double)ind->get_y()<=((this->numb_cells*this->cell_size/2)*(-1)) // if the x coordinate is lower than the minimum value
               
               //>= porque na paisagem quadrada as bordas mais distantes de 0 iniciariam um proximo pixel que estaria fora da paisagem. Ou teriamos que assumir que esses pixels mais extremos tenha uma área maior, o que daria um trabalho adicional para implementar uma situação irreal.)
               // Porque isso nao foi implementado em todos os casos?
			{
				for(unsigned int i=0; i<popIndividuos.size();i++) //Goess through all individuals
				{
					if(this->popIndividuos[i]->get_id()==(int)ind->get_id()) // Checks if the individual to be killed (erased) is the right one
                        //DUVIDA: porque tem int?
					{
						delete this->popIndividuos[i]; // Deletes the object in the chosen position
						this->popIndividuos.erase(this->popIndividuos.begin()+i); // Erases the position that object occupied in the vector of individuals
					}
				}
				emigrou = true; // Sets the boolean value to positive simbolizing the individual emigrated
			}
		}
		break;

		case 1:
		if(ind->get_x()<(this->numb_cells*this->cell_size/2)*(-1)) //Checks if the x coordinate is higher than the maximum value
			ind->set_x(this->tamanho+ind->get_x()); // Changes the x coordinate to represent the periodical condition (left corner to right corner)
               
		if(ind->get_x()>=this->numb_cells*this->cell_size/2) //Checks if the x coordinate is lower than the minimum value
			ind->set_x(ind->get_x()-this->tamanho); // Changes the x coordinate to represent the periodical condition (lright corner to left corner)
               
		if(ind->get_y()<(this->numb_cells*this->cell_size/2)*(-1)) //Checks if the y coordinate is higher than the maximum value
			ind->set_y(this->tamanho+ind->get_y());// Changes the y coordinate to represent the periodical condition (bottom corner to top corner)
               
		if(ind->get_y()>=this->numb_cells*this->cell_size/2 ) //Checks if the y coordinate is lower than the maximum value
			ind->set_y(ind->get_y()-this->tamanho);// Changes the y coordinate to represent the periodical condition (up corner to bottom corner)
		break;
	}
	return emigrou; // Returns a positive boolean value if the individual emigrated
	/* TBI
	case 2: reflexiva
	*/
}

// Function that calculate the distance between two individuals
double paisagem::calcDist(const individuo* a1, const individuo* a2) const
//Virou método da paisagem pois não estava conseguindo usar as propriedades privadas da paisagem nessa função (p. ex. this->boundary_condition). Tive que tirar o primeiro const.
{
	switch(this->boundary_condition) //(0= absortive, 1= periodical (pacman),2= reflexive)
	{
		case 0: // Euclidian Distance
			return sqrt(((double)a1->get_x()-(double)a2->get_x())*((double)a1->get_x()-(double)a2->get_x())+((double)a1->get_y()-(double)a2->get_y())*((double)a1->get_y()-(double)a2->get_y())); //Returns the euclidian distance between the two individuals
			break;

		case 1: // The Euclidian distance acounting for the periodical effect (individuals on diferent edges may be close to one another)
			double x1=a1->get_x(); // Gets the x coordinate of the first individuals
			double x2=a2->get_x(); // Gets the x coordinate of the second individuals
			double y1=a1->get_y(); // Gets the y coordinate of the first individuals
			double y2=a2->get_y(); // Gets the y coordinate of the second individuals
			double dx= x1 > x2 ? x1 - x2 : x2 - x1; // choses the lower between x1-x2 e x2-x1
			double dy= y1 > y2 ? y1 - y2 : y2 - y1;  //choses the lower between y1-y2 e y2-y1
			dx = dx > this->tamanho - dx ? this->tamanho - dx : dx; // choses the lower between tam-dx e dx
			dy = dy > this->tamanho - dy ? this->tamanho - dy : dy; //choses the lower between tam-dy e dy
			return sqrt(dx*dx + dy*dy); // Computes and returns the distance
			break;
	}
}

// A function to calculate de density of individuals according to density type (global or local) and considering landscape boundary effects in the calculation.
double paisagem::calcDensity(const individuo* ind1) const
{
	// Functions for local density calculation

	/* 1. Circular area defining a region in which denso-dependence occurs: landscape boundary effects.
	 In this case, density is the number of individuals inside the circle divided by circle area.
	 This is the same calculation as for global density, except by the cases in which landscape boundary affects
	 the area of the circle.
	 */
    
    double density; // Creates a temporary variable
    density=ind1->NBHood_size()/(M_PI*(ind1->get_raio()*ind1->get_raio())); // Computes the standard case of density (when an individual radius area is completely within the landscape)

	// Condition giving the boundary effects cases
	if(ind1->get_densType()==1) // checks if the density type is set to 1 (local/ within a radius distance)
	{
		if(ind1->get_x()*ind1->get_x()>((this->tamanho/2)-ind1->get_raio())*((this->tamanho/2)-ind1->get_raio()) || // Checks if the x coordinate is within a dadius distance from the edge OR
		   ind1->get_y()*ind1->get_y()>((this->tamanho/2)-ind1->get_raio())*((this->tamanho/2)-ind1->get_raio())) // if the y coordinate is within a dadius distance from the edge
		{
			// OBS: The absolute values retain their distance from the edge, so the compuatations from all quadrants can be done with a single set of equations
			double modIx = fabs(ind1->get_x()); // Stores the the absolute value of x
			double modIy = fabs(ind1->get_y()); // Stores the the absolute value of y
			double XYmax = this->tamanho/2; // Store the lenght of quadrant side
			vector<double>secX;
			vector<double>secY;

			// Functions for adjusted local density calculation, according to the specific case
			// 1)
			if(modIx>XYmax-ind1->get_raio() && modIy<=XYmax-ind1->get_raio()) // Checks if the x coordinate is within a dadius distance from the edge AND if the y coordinate is not
			{
				secX.push_back(XYmax);
				secX.push_back(XYmax);
				secY.push_back(modIy+sqrt(ind1->get_raio()*ind1->get_raio()-((XYmax-modIx)*(XYmax-modIx))));
				secY.push_back(modIy-sqrt(ind1->get_raio()*ind1->get_raio()-((XYmax-modIx)*(XYmax-modIx))));

				double distSecs = secY[0]-secY[1];
				double height = XYmax - modIx;
				double theta = acos(1-(distSecs*distSecs/(2*ind1->get_raio()*ind1->get_raio()))); // angle in radians
				double adjArea = M_PI*ind1->get_raio()*ind1->get_raio() - (theta*ind1->get_raio()*ind1->get_raio()/2 - (distSecs*height/2));

				density = ind1->NBHood_size()/adjArea;

			}

			// 2)
			if(modIx<=XYmax-ind1->get_raio() && modIy>XYmax-ind1->get_raio()) // Checks if the y coordinate is within a dadius distance from the edge AND if the x coordinate is not
			{
				secY.push_back(XYmax);
				secY.push_back(XYmax);
				secX.push_back(modIx+sqrt(ind1->get_raio()*ind1->get_raio()-((XYmax-modIy)*(XYmax-modIy))));
				secX.push_back(modIx-sqrt(ind1->get_raio()*ind1->get_raio()-((XYmax-modIy)*(XYmax-modIy))));

				double distSecs = secX[0]-secX[1];
				double height = XYmax - modIy;
				double theta = acos(1-(distSecs*distSecs/(2*ind1->get_raio()*ind1->get_raio()))); // angle in radians
				double adjArea = M_PI*ind1->get_raio()*ind1->get_raio() - (theta*ind1->get_raio()*ind1->get_raio()/2 - (distSecs*height/2));

				density = ind1->NBHood_size()/adjArea;
			}


			if(modIx>XYmax-ind1->get_raio() && modIy>XYmax-ind1->get_raio())
			{

				// 3)
				if((modIx-XYmax)*(modIx-XYmax)+(modIy-XYmax)*(modIy-XYmax)>ind1->get_raio()*ind1->get_raio())
				{
					secX.push_back(modIx+sqrt(ind1->get_raio()*ind1->get_raio()-((XYmax-modIy)*(XYmax-modIy))));
					secY.push_back(XYmax);
					secX.push_back(modIx-sqrt(ind1->get_raio()*ind1->get_raio()-((XYmax-modIy)*(XYmax-modIy))));
					secY.push_back(XYmax);
					secX.push_back(XYmax);
					secY.push_back(modIy+sqrt(ind1->get_raio()*ind1->get_raio()-((XYmax-modIx)*(XYmax-modIx))));
					secX.push_back(XYmax);
					secY.push_back(modIy-sqrt(ind1->get_raio()*ind1->get_raio()-((XYmax-modIx)*(XYmax-modIx))));

					double distSecs = sqrt((secX[3]-secX[1])*(secX[3]-secX[1])+(secY[1]-secY[3])*(secY[1]-secY[3]));
					double distSecs2 = sqrt((secX[2]-secX[0])*(secX[2]-secX[0])+(secY[0]-secY[2])*(secY[0]-secY[2]));
					double theta = acos(1-(distSecs*distSecs/(2*ind1->get_raio()*ind1->get_raio()))); // angle in radians
					double phi = acos(1-(distSecs2*distSecs2/(2*ind1->get_raio()*ind1->get_raio()))); // angle in radians
					double adjArea = (2*M_PI-theta)*ind1->get_raio()*ind1->get_raio()/2 + phi*ind1->get_raio()*ind1->get_raio()/2 + (secX[0]-secX[1])*(XYmax-modIy)/2 + (secY[2]-secY[3])*(XYmax-modIx)/2;

					density = ind1->NBHood_size()/adjArea;

				}
				// 4)
				else
				{
					secX.push_back(modIx-sqrt(ind1->get_raio()*ind1->get_raio()-((XYmax-modIy)*(XYmax-modIy))));
					secY.push_back(XYmax);
					secX.push_back(XYmax);
					secY.push_back(modIy-sqrt(ind1->get_raio()*ind1->get_raio()-((XYmax-modIx)*(XYmax-modIx))));

					double distSecs = sqrt((secX[1]-secX[0])*(secX[1]-secX[0])+(secY[0]-secY[1])*(secY[0]-secY[1]));
					double theta = acos(1-(distSecs*distSecs/(2*ind1->get_raio()*ind1->get_raio()))); // angle in radians
					double adjArea = theta*ind1->get_raio()*ind1->get_raio()/2 + (XYmax-secX[0])*(XYmax-modIy) + (XYmax-secY[1])*(XYmax-modIx);

					density = ind1->NBHood_size()/adjArea;
				}

			}
		}
	}

	/* 2. Density kernel (TBI).

	 if(ind1->get_densType()==2) {}

	 */
	return density;
}


/*
  Sempre adicione const aos argumentos de métodos quando o método não
  deve alterá-los. Previne vários erros e pode otimizar compilação
*/

// Function that updates an individuals neighbourhood list
void paisagem::atualiza_vizinhos(individuo * const ag1) const //acessando os vizinhos dos agentes
{
	vector <individuo*> listViz; // Creates a vector of individuals to store the neighbours
    
	if(ag1->get_densType()==0) // Checks if the density type is set to 0 (global)
        
        //dens_type poderia voltar como propriedade da paisagem. Facilitariam as coisas. Como muitas propriedades e métodos deste código, elas podem ser interpretadas das duas formas (como do individuo ou como da paisagem). O que está dando confusão é que estamos fazendo um IBM, mas para algumas situações estamos querendo simular dinâmicas cujas variáveis de interesse são propriedades populacionais e não do indivíduo. Se aceito, limar o método get_densType() do individuo.h.
	{
		for(unsigned int j=0; j<popIndividuos.size(); j++) //Goes through all the individuals
		{
			individuo* ag2=this->popIndividuos[j]; //Stores a copy of each individual
			if(ag1==ag2) continue; //Except if it is the focal individual itself
			listViz.push_back(ag2);// Stores the copy of each other individual in a neighbours vector
		}
	}
	else //Checks if the density type is not set to 0 (1= local)
	{
		double rad = (double)ag1->get_raio(); // Gets the radius of the focal individual
		for(unsigned int j=0; j<popIndividuos.size(); j++) //Goes through all currently alive individuals
		{
			individuo* ag2=this->popIndividuos[j]; //Stores a copy of each individual
			if(ag1==ag2) continue; //Except if it is the focal individual itself
			double d=this->calcDist(ag1,ag2); // Checks if the copyed individual is within a radius distance
			if(d<=rad) {listViz.push_back(ag2); // Stores the copy of each other individual within a radius distance into a neighbours vector
                
            }
		}
	}
	ag1->set_vizinhos(listViz);

}

// Function that updates if an individual is currently on habitat or matrix
void paisagem::atualiza_habitat(individuo * const ind) const
{
	// Tinha um IF com landscape_shape que eu removi. Não entendi como a paisagem ser circular
	// interfere em ser habitat ou não: isso deve interferir na apply_boundary apenas, certo?
	// Also: Tinha uma inversão do y que eu também não entendi e removi
	// A.C. 10.07.13

	// Um termo (-1) foi removido erroneamente por A.C.. Para o hy, o sentido em que o número de células aumenta é o
	//contrário do sentido em que as coordenadas aumentam. Portanto a multiplicação por - 1 é necessária.
	// M.A. 12.09.14
	int hx,hy; //Initializes variables for storing pixels coordinats
	hx= (double)ind->get_x()/this->cell_size+this->numb_cells/2;// Finds the pixel X coodenate
	hy= ((double)ind->get_y()/this->cell_size)*(-1)+this->numb_cells/2; //Finds the pixel y coodenate
	ind->set_habitat(this->landscape[hx][hy]);// Sets the patch label of the individual to the label of the patch it is currently on
}

// Function to update the patch the individual is currently in
void paisagem::atualiza_patch(individuo * const ind) const
{
	int hx,hy; //Initializes variables for storing pixels coordinats
	hx= (double)ind->get_x()/this->cell_size+this->numb_cells/2; // Finds the pixel X coodenate
	hy= ((double)ind->get_y()/this->cell_size)*(-1)+this->numb_cells/2; //Finds the pixel y coodenate
	ind->set_patch(this->patches[hx][hy]); // Sets the patch label of the individual to the label of the patch it is currently on
}

// Retursive function that finds and identfy each landscape path
void paisagem::find_patches(int x, int y, int current_label)
{
  if (x < 0 || x == this->numb_cells) return; // Checks if the x coordinate is out of bounds
  if (y < 0 || y == this->numb_cells) return; // Checks if the y coordinate is out of bounds
  if (this->patches[x][y] || !this->landscape[x][y]) return; // Checks if the curret pixel is already labeled AND if is a habitat pixel

  // mark the current cell with its respective patch ID value
  this->patches[x][y] = current_label;

  // recursively mark the neighbors
  find_patches(x + 1, y, current_label);
  find_patches(x, y + 1, current_label);
  find_patches(x - 1, y, current_label);
  find_patches(x, y - 1, current_label);
}

//Function that returns if an individual was born during simulation runtime
const bool paisagem::nascido(individuo * const ind) const
{
    return ind->get_id() > this->N;// Returns positive boolean if the individual was born during the simulation
}

//funtion that returns the number of patches on the landscape
int paisagem::get_numb_patches()
{
    return numb_patches;//Returns the number of patches on the landscape
    
}
 double paisagem::get_patch_area(int i) const
{
    return this->patch_area[i];
    
}
               
double paisagem::walk(int lower){
                
    if(selection=TRUE)
    {
        double possibilitities[this->popIndividuos[lower]->points][3];
                    
        for (int i=0; i<this->popIndividuos[lower]->points; i++)
        {
                        
            choice=runif(0.0,360.0);
            dist=runif(0.0,this->popIndividuos[lower]->passo);
                        
            possibilitities[i][0]=this->popIndividuos[lower]->x+cos(choice)*dist;
            possibilitities[i][1]=this->popIndividuos[lower]->y+sin(choice)*dist;
                        
            possibilitities[i][2]<-dnorm(landscape[possibilitities[i][0]][possibilitities[i][1]], this->popIndividuos[lower]->genotype_mean, this->popIndividuos[lower]->width_sd);
        }
            this->popIndividuos[lower]->habitat_selection(possibilitities);
                    
                    
    }
                
    else
    {
        this->popIndividuos[lower]->anda(); // Calls a function that changes the X and Y coordinates of the individuals
    }
}
                

                


//TBI
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
