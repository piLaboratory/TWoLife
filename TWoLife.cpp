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

/** \file TWoLife.cpp \brief Arquivo usado na integração R/C
 * 
 * Este arquivo define a função TWoLife em C, que será chamada a partir do R */


/** Funcao que recebe os parâmetros da simulação e a executa,
 // baseada no texto em http://users.stat.umn.edu/~geyer/rc/
 //
 // NOTE que funções em C++ precisam ser indicadas como extern "C"
 // para poderem ser acessadas facilmente pela interface do R!
 //
 // Veja a descrição da classe \ref paisagem para o significado dos parâmetros
 // */
extern "C" void TWoLife (double * raio, int * N, double * angulo_visada, double * passo, double * move,
						 double * taxa_basal, double * taxa_morte, double * incl_b, double * incl_d,
						 int * numb_cells, double * cell_size, int * land_shape, int * density_type, 
						 double * death_mat, int * inipos, int * bound_condition, int * scape, double * tempo, int * nPop, double * x, double * y, int * outCode)
{
	// This sequence creates an attribute containing the output file name. The template is output-00000.txt.
	string fileNAME = "output-00000.txt";
	stringstream tmps;
	tmps<<outCode[0];
	string addToName = tmps.str();
	int fnSize = fileNAME.size();
    int tmpsSize = addToName.size();
	fileNAME.erase(fileNAME.begin()+fnSize-4-tmpsSize,fileNAME.begin()+fnSize-4);
	fileNAME.insert(fnSize-4-tmpsSize,addToName);
	
	GetRNGstate(); /* (mudar para doxygen):  este comando chama o engine de numeros aleatorios do R
					Por causa dela nossa biblioteca nao eh standalone */
	
	paisagem* floresta = new paisagem(raio[0], N[0], angulo_visada[0], passo[0], 
									  move[0], taxa_basal[0], taxa_morte[0], incl_b[0], 
									  incl_d[0], numb_cells[0], cell_size[0], land_shape[0],
									  density_type[0], death_mat[0], inipos[0], bound_condition[0], 
									  scape);
	
	ofstream outputSIM; // ofstream for the output file
	outputSIM.open(fileNAME.c_str());

	outputSIM << "Lado: " << floresta->get_tamanho() << endl;
	outputSIM << "Tamanho pixel: " << cell_size[0] << endl;
	outputSIM << "Pixels lado: " << numb_cells[0] << endl;
	outputSIM << "N fragmentos: " << floresta->get_numb_patches() << endl;
	outputSIM << "Area fragmentos: ";
	double sum = 0;
	//outputSIM << floresta->get_patch_area(0) << " "; //Imprimir area matriz?
	for (int i = 1; i <= floresta->get_numb_patches(); i++)
	{
		sum += floresta->get_patch_area(i);
		outputSIM << floresta->get_patch_area(i) << " ";		
	}
	outputSIM << "\nProporção de habitat: " << sum/(floresta->get_tamanho()*floresta->get_tamanho())<< "\n\n\n\n";

	

	for(unsigned int i=0; i<floresta->conta_individuos();i++)
	{
		outputSIM << floresta->tempo_do_mundo << " " << floresta->get_individuos(i)->get_id() << " "  << floresta->get_individuos(i)->get_patch() << " " << floresta->get_individuos(i)->get_x() << " " << floresta->get_individuos(i)->get_y() << endl;
	}
	while (floresta->tempo_do_mundo < tempo[0] && floresta->conta_individuos() > 0)
	{
		int ind_neo = floresta->sorteia_individuo();
		int acao = floresta->sorteia_acao(ind_neo);
		floresta->atualiza_tempo(ind_neo);

		int ID_neo = floresta->get_individuos(ind_neo)->get_id();
		double x_neo = floresta->get_individuos(ind_neo)->get_x();
		double y_neo = floresta->get_individuos(ind_neo)->get_y();
		int patch_neo = floresta->get_individuos(ind_neo)->get_patch();
		
		bool emigrou = floresta->realiza_acao(acao, ind_neo);
		floresta->update();

		if(acao==2 && !emigrou)
		{
			x_neo = floresta->get_individuos(ind_neo)->get_x();
			y_neo = floresta->get_individuos(ind_neo)->get_y();
			patch_neo = floresta->get_individuos(ind_neo)->get_patch();
		}
		else if(acao==2 && emigrou)
			acao = 3;

		outputSIM << floresta->tempo_do_mundo << " " << acao << " " << ID_neo << " "  << patch_neo << " " << x_neo << " " << y_neo << endl;
	}
	outputSIM<< "EOF\n";
	outputSIM.close(); //end of output file
	
	
	*nPop = floresta->conta_individuos();
	for (int i =0; i < *nPop; i ++) {
		x[i] = floresta->get_individuos(i)->get_x();
		y[i] = floresta->get_individuos(i)->get_y();
	} //DUVIDA: porque x[i] e y[i] nao tem asterisco antes?
	delete floresta;
	PutRNGstate();
}


