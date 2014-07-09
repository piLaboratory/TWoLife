#include <iostream>
#include <vector>
#include <cstdlib>
#include <R.h>
#include <Rmath.h>
#include "paisagem.cpp" 
#include "individuo.cpp" 

using namespace std;

/** \file TWoLife.cpp \brief Arquivo usado na na integração C++/R
 * 
 * Este arquivo define a função TwoLife em C, que será usada a partir do R. */


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
						   double * death_mat, int * bound_condition, double * cover, double * tempo, int * nPop, double * x, double * y)
{
  GetRNGstate(); /* (mudar para doxygen):  este comando chama o engine de numeros aleatorios do R
		    Por causa dela  nossa biblioteca nao eh standalone */
	paisagem* floresta = new paisagem(raio[0], N[0], angulo_visada[0], passo[0], 
					  move[0], taxa_basal[0], taxa_morte[0], incl_b[0], 
					  incl_d[0], numb_cells[0], cell_size[0], land_shape[0],
					  density_type[0], death_mat[0], bound_condition[0], 
					  cover[0]);
	while (floresta->tempo_do_mundo < tempo[0] && floresta->conta_individuos() > 0) {
			floresta->update();
	}
	*nPop = floresta->conta_individuos();
	for (int i =0; i < *nPop; i ++) {
		x[i] = floresta->get_individuos(i)->get_x();
		y[i] = floresta->get_individuos(i)->get_y();
	} //DUVIDA: porque x[i] e y[i] nao tem asterisco antes?
	delete floresta;
	PutRNGstate();
}
