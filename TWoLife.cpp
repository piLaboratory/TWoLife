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

/** \file TWoLife.cpp \brief File used for R/C++ interaction
 *
 * This file defines the TWoLife function in C, to be be called by R */

/** The function recieves the inputed parameters and executes the simulation
 // Based on the text avaiable on http://users.stat.umn.edu/~geyer/rc/
 //
 // Observe that funtions in C++ need to be called as extern "C"
 // as to be easely accessd by the R interface
 //
 // See the description off the \ref paisagem class for parameter meenings
 // */
 extern "C" void TWoLife (double * raio,
                         int * N,
                         double * angulo_visada,
                         double * passo,
                         double * taxa_move,
                         double * taxa_basal,
                         double * taxa_morte,
                         double * incl_b,
                         double * incl_d,
                         int * numb_cells,
                         double * cell_size,
                         int * land_shape,
                         int * density_type,
                         double * death_mat,
                         double * move_mat,
                         int * inipos,
                         int * bound_condition,
                         double * scape,
                         double * tempo,
                         double * initialPosX,
                         double * initialPosY,
                         double * genotype_means,
                         double  * width_sds,
                         int * points,
                         bool * Null,// Output "inputs"
                         int * nPop,
                         double * x,
                         double * y,
                         double * m,
                         char ** outCode)
 {
   // Names the output file name
   stringstream tmps; // Creates an object of the stringstream to recieve the output name
   tmps<<outCode[0]; // Stores the inputed name in the stringstream object
   string fileNAME = tmps.str(); // Stores the desired file name as Filename
   
   GetRNGstate(); /* This comands calls R's random number generator */

/** Generates an object of the paisagem class using the inputed simulation parameters  */
paisagem* floresta = new paisagem(raio[0],
                                  N[0],
                                   angulo_visada[0],
                                                passo[0],
                                                     taxa_move[0],
                                                              taxa_basal[0],
                                                                        taxa_morte[0],
                                                                                  incl_b[0],
                                                                                        incl_d[0],
                                                                                              numb_cells[0],
                                                                                                        cell_size[0],
                                                                                                                 land_shape[0],
                                                                                                                           density_type[0],
                                                                                                                                       death_mat[0],
                                                                                                                                                move_mat[0],
                                                                                                                                                        inipos[0],
                                                                                                                                                              bound_condition[0],
                                                                                                                                                                             scape,
                                                                                                                                                                             initialPosX,
                                                                                                                                                                             initialPosY,
                                                                                                                                                                             genotype_means,
                                                                                                                                                                             width_sds,
                                                                                                                                                                             points,
                                                                                                                                                                             Null[0]);

ofstream outputSIM; // Creates ofstream for the output file
outputSIM.open(fileNAME.c_str()); // Opens the ofstream for the output file

// Outputs the simulation starting conditions
outputSIM << "Lado: " << floresta->get_tamanho() << endl; // Outputs the lenght of the landscape side
outputSIM << "Tamanho pixel: " << cell_size[0] << endl; // Outputs the lenght of a pixel side
outputSIM << "Pixels lado: " << numb_cells[0] << endl; // Outputs the number of pixels on the landscapeside
outputSIM << "N fragmentos: " << floresta->get_numb_patches() << endl; // Outputs the number of fragments on the landscape
outputSIM << "Area fragmentos: "; // Start Outputing the area of each of the fragments


double sum = 0;// Creates a accumulator variable
//outputSIM << floresta->get_patch_area(0) << " "; //Imprimir area matriz?
for (int i = 1; i <= floresta->get_numb_patches(); i++) // Passes through each patch
{
  sum += floresta->get_patch_area(i); // Sums a patch area to the acumulator
  outputSIM << floresta->get_patch_area(i) << " "; // Outputs the area of each fragment
}

outputSIM << "\nProporção de habitat: " << sum/(floresta->get_tamanho()*floresta->get_tamanho())<< "\n\n\n\n"; // Outputs the proportion of the landscape that is habitable




for(unsigned int i=0; i<floresta->conta_individuos();i++) // Passes through each starting individual
{
  // Outputs all the starting individuals and their informations
  outputSIM << floresta->tempo_do_mundo << " " << floresta->get_individuos(i)->get_id() << " "  << floresta->get_individuos(i)->get_patch() << " " << floresta->get_individuos(i)->get_x() << " " << floresta->get_individuos(i)->get_y() << " " << floresta->get_individuos(i)->get_genotype_mean() << endl;
}

/** Performs the simulation  */
while (floresta->tempo_do_mundo < tempo[0] && floresta->conta_individuos() > 0 && floresta->conta_individuos() < 2000) // Performs the simulation until the maximum desired time is reached
{
  int ind_neo = floresta->sorteia_individuo(); // Drafts an individual to perform an action
  int acao = floresta->sorteia_acao(ind_neo); // Drafts an action for the individual to perform
  floresta->atualiza_tempo(ind_neo);  // Updstes the world Clock
  
  int ID_neo = floresta->get_individuos(ind_neo)->get_id(); // Obtains the ID of the drafted individual
  double x_neo = floresta->get_individuos(ind_neo)->get_x(); // Obtains the X coordinate of the drafted individual
  double y_neo = floresta->get_individuos(ind_neo)->get_y(); // Obtains the Y coordinate of the drafted individual
  int patch_neo = floresta->get_individuos(ind_neo)->get_patch(); // Obtains the curent patch of the drafted individual
  double m_neo = floresta->get_individuos(ind_neo)->get_genotype_mean();
  
  bool emigrou = floresta->realiza_acao(acao, ind_neo); // Makes the drafted individual perform the chosen action
  
  if(acao==2 && emigrou) //Checks if the action performed was a migration event AND the drafted individual left the landscape boundary
    floresta->update(0, ind_neo); // Updates the drafted individual migration event as a death event
  else // Comtemplates all other cases
    floresta->update(acao, ind_neo); // Updates the individual according to the performed action
  
  if(acao==2 && !emigrou) //Checks if the action performed was a migration event AND the drafted individuals new location is within the landscape boundary
  {
    x_neo = floresta->get_individuos(ind_neo)->get_x(); // Obtains the individuals current x coordinate
    y_neo = floresta->get_individuos(ind_neo)->get_y(); // Obtains the individuals current y coordinate
    patch_neo = floresta->get_individuos(ind_neo)->get_patch(); // Obtains the individuals current patch
  }
  else if(acao==2 && emigrou) //Checks if the action performed was a migration event AND the drafted individual left the landscape boundary
    acao = 3; // Sets the output entry of the individual as an emigration event
  
  // Outputs the actions performed by the drafted individual
  outputSIM << floresta->tempo_do_mundo << " " << acao << " " << ID_neo << " "  << patch_neo << " " << x_neo << " " << y_neo << " " << m_neo << endl;
}


outputSIM<< "EOF\n";// Outbuts a warning stating the simulation was run to completion
outputSIM.close(); //Closes the output file

//Outputs the remaining individuals
*nPop = floresta->conta_individuos(); // Checks how many individuals are left
for (int i =0; i < *nPop; i ++) { // Passes through each remaining individual
  x[i] = floresta->get_individuos(i)->get_x(); // Outputs its x coordinate
  y[i] = floresta->get_individuos(i)->get_y(); // Outputs its y coordinate
  m[i] = floresta->get_individuos(i)->get_genotype_mean();
}

delete floresta;

PutRNGstate();
 }
 
 