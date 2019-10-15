/*
 TWoLife.cpp file used for integrating R/C
 This file defines the TWoLife function in C, to be called by R
 This function receives the simulation parameters and executes,
 Based on http://users.stat.umn.edu/~geyer/rc/
 
 Note that functions in C++ need to be indicated with the " extern "C" "
 So that it can be easily accessed by R
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <R.h>
#include <Rmath.h>
#include "paisagem.cpp" 
#include "individuo.cpp" 
#include <sstream>
#include <string>
using namespace std;

// Function Declarations
extern "C" void TWoLife (double * raio, int * N, double * angulo_visada, double * passo, double * move, double * taxa_basal, double * taxa_morte, double * incl_b, double * incl_d, int * numb_cells, double * cell_size, int * land_shape, int * density_type, double * death_mat, int * inipos, int * bound_condition, int * scape,double * tempo, int * nPop, double * x, double * y, int * outCode);

/*
 Function to be called by R
 All the parameters are Pointers!
 
 Paran double* raio - Radius of density dependant influence
 Paran int* N - Number of individuals at the start of the simulation
 Paran double* angulo_visada - Angle used for orientation when dispersing
 Paran double* passo - The Lenght distance of a dispersal event
 Paran double* move - The rate at which the individuals disperse
 Paran double* taxa_basal - The basal birth rate (The rate at which the individuals give birth on a habitat patch without neigbours)
 Paran double * taxa_morte - The basal death rate (The rate at which the individuals die on a habitat patch without neigbours)
 Paran double * incl_b - The slope of the birth density dependece function
 Paran double * incl_d - The slope of the death density dependence function
 Paran int * numb_cells - Number of pixels in a square landscape side (if land_shape == 1)
 Paran double * cell_size - Resolution Lenght of a square pixel side
 Paran int * land_shape - Shape of the landscape (0= circle, 1= square)
 Paran int * density_type - Density type (0 = global, 1 = local/within a individual radius)
 Paran double * death_mat - Constant that indicates how many times higher the death rate should be on non-habitat pixels
 Paran int * inipos - The initial postion of individuals (0 = origin, 1 = random, 2 = normaly distributed with mean on origin)
 Paran int * bound_condition - The boundary condition type affects how individuals interact with the edges of the landscape (0= absortive, 1= periodical (pacman),2= reflexive)
 Paran int * scape - Vector containing the evironmental values of the landscape pixels (0= matrix, 1= habitat)
 Paran double * tempo - The maximum number of events to happen until the end of the simulation
 Paran int * nPop - Variable for storing the number of individuals at the end of the simulation
 Paran double * x - Vector to store the final x coodinate of the living individuals at the end of the simulation
 double * y - Vector to store the final y coodinate of the living individuals at the end of the simulation
 Paran int * outCode -  A value used for identifying the simulatio run
 */


// Function Definition to be called by R
extern "C" void TWoLife (double * raio, int * N, double * angulo_visada, double * passo, double * move, double * taxa_basal, double * taxa_morte, double * incl_b, double * incl_d, int * numb_cells, double * cell_size, int * land_shape, int * density_type, double * death_mat, int * inipos, int * bound_condition, int * scape,double * tempo, int * nPop, double * x, double * y, int * outCode)
{
    
	GetRNGstate(); //Sets the random seed based on the R (For stand-alone one could use "srand(time(0));")
    
    // Creates an object of the paisagem/landscape class using the R inputed values as parameters
    // This consequently creates a vector o objects of the individuo/individual class as an atribute of the paisagem/landscape
    // The conditions are completely set for the first time step of the simulation
	paisagem* floresta = new paisagem(raio[0], N[0], angulo_visada[0], passo[0], move[0], taxa_basal[0], taxa_morte[0], incl_b[0], incl_d[0], numb_cells[0], cell_size[0], land_shape[0], density_type[0], death_mat[0], inipos[0], bound_condition[0], scape);
    
    
    // This sequence creates an attribute containing the output file name. The template is output-00000.txt.
    string fileNAME = "output-00000.txt"; // Creates a string variable and inserts the template name in it
    stringstream tmps; //Creates an object of the strinstream class
    tmps<<outCode[0]; // Stores the simulation number imputed
    string addToName = tmps.str(); // Copies the simulation run number into another object
    int fnSize = fileNAME.size(); // Stores the number o digits of the template name on an object
    int tmpsSize = addToName.size(); // Stores the number o digits of the simulation run number on an object
    fileNAME.erase(fileNAME.begin()+fnSize-4-tmpsSize,fileNAME.begin()+fnSize-4);// Erases the space needed for storing the simulation run number
    fileNAME.insert(fnSize-4-tmpsSize,addToName);// Inserts the simulation ru number
	
	ofstream outputSIM; // ofstream for the output file
	outputSIM.open(fileNAME.c_str()); // Names the output file

	outputSIM << "Lado: " << floresta->get_tamanho() << endl; // Outputs the lenght of a side of the landscape (if land_shape == 1)
	outputSIM << "Tamanho pixel: " << cell_size[0] << endl; // Outputs the resolution Lenght of a square pixel side
	outputSIM << "Pixels lado: " << numb_cells[0] << endl; // Outputs the Number of pixels in a side of a square landscape (if land_shape == 1)
	outputSIM << "N fragmentos: " << floresta->get_numb_patches() << endl; // Outputs the number of patches in the landscape
	outputSIM << "Area fragmentos: "; // Outputs the area of each fragment
	double sum = 0; // creates a variable for storing the pach area sum
	//outputSIM << floresta->get_patch_area(0) << " "; // Turn on if you want to also output the area of the matrix
	for (int i = 1; i <= floresta->get_numb_patches(); i++)  //Goes through each patch
	{
		sum += floresta->get_patch_area(i); // Sums each patch area
		outputSIM << floresta->get_patch_area(i) << " ";// Outputs each patch area
	}
	outputSIM << "\nProporção de habitat: " << sum/(floresta->get_tamanho()*floresta->get_tamanho())<< "\n\n\n\n";n //Outputs the habitat/matrix proportion

	for(unsigned int i=0; i<floresta->conta_individuos();i++) //Goes through each individual
	{
		outputSIM << floresta->tempo_do_mundo << " " << floresta->get_individuos(i)->get_id() << " "  << floresta->get_individuos(i)->get_patch() << " " << floresta->get_individuos(i)->get_x() << " " << floresta->get_individuos(i)->get_y() << endl; // Outputs the starting conditions of each individual
	}
    
    // Loop that constitutes the simulation time steps
    while (floresta->tempo_do_mundo < tempo[0] && floresta->conta_individuos() > 0) // goes through all the desired time steps
	{
		int ind_neo = floresta->sorteia_individuo(); // Calls a function(Landscape/paisagem) to return the individual with the lowest randowly drafted time
		int acao = floresta->sorteia_acao(ind_neo); // Calls a function(Landscape/paisagem) to call a function(individual/individuao) to randomly determine and return an action to be performed by the before drafted individual
		floresta->atualiza_tempo(ind_neo); // Calls a function(Landscape/paisagem) to update the simulation clock by adding the time drafted by the selected individual

		int ID_neo = floresta->get_individuos(ind_neo)->get_id(); // Calls a function(Landscape/paisagem) that returns the unique ID of the drafted individual
		double x_neo = floresta->get_individuos(ind_neo)->get_x(); // Calls a function(Landscape/paisagem) that returns the current X coordinate of the drafted individual
		double y_neo = floresta->get_individuos(ind_neo)->get_y(); // Calls a function(Landscape/paisagem) that returns the current y coordinate of the drafted individual
		int patch_neo = floresta->get_individuos(ind_neo)->get_patch(); // Calls a function(Landscape/paisagem) that returns the path the drafted individual is currently on
		
		bool emigrou = floresta->realiza_acao(acao, ind_neo); // Calls a function(Landscape/paisagem) that returns performs the action (0=death, 1=birth, 2=dispersal) and returns a true boolean value if an individual reaches the edge of the spatial domain
		floresta->update(); // Calls a function(Landscape/paisagem) to update the landscape object and consequently also updates each object of the individual class

		if(acao==2 && !emigrou)// Checks if there was a dispersal event AND if the final coordinates of that individual is withn the landscape boundary
		{
			x_neo = floresta->get_individuos(ind_neo)->get_x();// Stores the curent x coordinate of the drafted individuals
			y_neo = floresta->get_individuos(ind_neo)->get_y();// Stores the curent y coordinate of the drafted individuals
			patch_neo = floresta->get_individuos(ind_neo)->get_patch(); //Stores the path the drafted individual is currently on
		}
		else if(acao==2 && emigrou) // Checks if there was a dispersal event AND if the final coordinates of that individual is outside of the landscape boundary
			acao = 3;

		outputSIM << floresta->tempo_do_mundo << " " << acao << " " << ID_neo << " "  << patch_neo << " " << x_neo << " " << y_neo << endl; // Outputs which individuals did what at each time step
	}
    
	outputSIM<< "EOF\n";
	outputSIM.close(); //ends the output file
	
	// This section stores the coordenates of each living individuals at the end of the simulation (to "export" it to R)
	*nPop = floresta->conta_individuos(); //Discovers and stores the number of living individuals
	for (int i =0; i < *nPop; i ++) {// Goes through all the individuals
		x[i] = floresta->get_individuos(i)->get_x(); //Stores the X coordinate
		y[i] = floresta->get_individuos(i)->get_y(); //Stores the y coordinate
	} //DUVIDA: porque x[i] e y[i] nao tem asterisco antes?
    
    
	delete floresta; // Deletes the object of the landscape/paisagem class
    
	PutRNGstate(); //Ends the random seed based on the R (For stand-alone one could use "srand(time(0));")
}


