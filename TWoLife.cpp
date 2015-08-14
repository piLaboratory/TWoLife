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
                         double * death_mat, int * inipos, int * bound_condition, int * scape, double * tempo, int * outCode)
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
    /*
     for(unsigned int i=0; i<floresta->conta_individuos();i++)
     {
     outputSIM << floresta->tempo_do_mundo << " " << floresta->get_individuos(i)->get_id() << " " << floresta->get_individuos(i)->get_x() << " " << floresta->get_individuos(i)->get_y() << endl;
     }
     */
    
    /////// Defining the simulation CASE (i.e. pure exponential, pure logistic, random walk or Skellam) /////////
    int SIM_Type = 99;
    if(incl_b[0]+incl_d[0]==0 && move[0]==0)
    {
        SIM_Type = 1;
    }
    else
    {
        if(incl_b[0]+incl_d[0]>0 && move[0]==0)
        {
            if(density_type[0]==1)
            {
                SIM_Type = 3;
            }
            else {SIM_Type = 2;}
        }
        else
        {
            if(taxa_basal[0]+taxa_morte[0]==0 && move[0]>0)
            {
                SIM_Type = 4;
            }
            else
            {
                if(incl_b[0]+incl_d[0]==0 && move[0]>0)
                {
                    SIM_Type = 5;
                }
                else
                {
                    if(incl_b[0]+incl_d[0]>0 && move[0]>0 && density_type[0]==0)
                    {
                        SIM_Type = 6;
                    }
                    else
                    {
                        //if(incl_b[0]+incl_d[0]>0 && move[0]>0 && density_type[0]==1)
                        //{
                        SIM_Type = 7;
                        //}
                    }
                }
            }
        }
    }
    
    //////////////////////////////////////////////////////////////////////////////////////////
    switch (SIM_Type)
    {
        case 1:
            for(unsigned int i=0; i<floresta->conta_individuos();i++)
            {
                outputSIM << floresta->tempo_do_mundo << " " << floresta->get_individuos(i)->get_id() << " " << endl;
            }
            
            while (floresta->tempo_do_mundo < tempo[0] && floresta->conta_individuos() > 0)
            {
                double t_ant = floresta->tempo_do_mundo;
                int lowerInd=floresta->updateEXP();
                if(t_ant < (int)floresta->tempo_do_mundo)
                {
                    for(unsigned int i=0; i<floresta->conta_individuos();i++)
                    {
                        outputSIM << (int)floresta->tempo_do_mundo << " " << floresta->get_individuos(i)->get_id() << " " << endl;
                    }
                }
                floresta->realiza_acao(lowerInd);
            }
            break;
            
        case 2:
            for(unsigned int i=0; i<floresta->conta_individuos();i++)
            {
                outputSIM << floresta->tempo_do_mundo << " " << floresta->get_individuos(i)->get_id() << " " << endl;
            }
            
            while (floresta->tempo_do_mundo < tempo[0] && floresta->conta_individuos() > 0)
            {
                double t_ant = floresta->tempo_do_mundo;
                int lowerInd=floresta->updateLOG();
                if(t_ant < (int)floresta->tempo_do_mundo)
                {
                    for(unsigned int i=0; i<floresta->conta_individuos();i++)
                    {
                        outputSIM << (int)floresta->tempo_do_mundo << " " << floresta->get_individuos(i)->get_id() << " " << endl;
                    }
                }
                floresta->realiza_acao(lowerInd);
            }
            break;
            
        case 3:
            for(unsigned int i=0; i<floresta->conta_individuos();i++)
            {
                outputSIM << floresta->tempo_do_mundo << " " << floresta->get_individuos(i)->get_id() << " " << endl;
            }
            
            while (floresta->tempo_do_mundo < tempo[0] && floresta->conta_individuos() > 0)
            {
                double t_ant = floresta->tempo_do_mundo;
                int lowerInd=floresta->updateLOGL();
                if(t_ant < (int)floresta->tempo_do_mundo)
                {
                    for(unsigned int i=0; i<floresta->conta_individuos();i++)
                    {
                        outputSIM << (int)floresta->tempo_do_mundo << " " << floresta->get_individuos(i)->get_id() << " " << endl;
                    }
                }
                floresta->realiza_acao(lowerInd);
            }
            break;
            
        case 4:
            for(unsigned int i=0; i<floresta->conta_individuos();i++)
            {
                outputSIM << floresta->tempo_do_mundo << " " << floresta->get_individuos(i)->get_id() << " " << floresta->get_individuos(i)->get_x() << " " << floresta->get_individuos(i)->get_y() << endl;
            }
            
            while (floresta->tempo_do_mundo < tempo[0] && floresta->conta_individuos() > 0)
            {
                double t_ant = floresta->tempo_do_mundo;
                int lowerInd=floresta->updateRW();
                if(t_ant < (int)floresta->tempo_do_mundo)
                {
                    for(unsigned int i=0; i<floresta->conta_individuos();i++)
                    {
                        outputSIM << (int)floresta->tempo_do_mundo << " " << floresta->get_individuos(i)->get_id() << " " << floresta->get_individuos(i)->get_x() << " " << floresta->get_individuos(i)->get_y() << endl;
                    }
                }
                floresta->doActionRW(lowerInd);
            }
            break;
            
        case 5:
            for(unsigned int i=0; i<floresta->conta_individuos();i++)
            {
                outputSIM << floresta->tempo_do_mundo << " " << floresta->get_individuos(i)->get_id() << " " << floresta->get_individuos(i)->get_x() << " " << floresta->get_individuos(i)->get_y() << endl;
            }
            
            while (floresta->tempo_do_mundo < tempo[0] && floresta->conta_individuos() > 0)
            {
                double t_ant = floresta->tempo_do_mundo;
                int lowerInd=floresta->updateSKEXP();
                if(t_ant < (int)floresta->tempo_do_mundo)
                {
                    for(unsigned int i=0; i<floresta->conta_individuos();i++)
                    {
                        outputSIM << (int)floresta->tempo_do_mundo << " " << floresta->get_individuos(i)->get_id() << " " << floresta->get_individuos(i)->get_x() << " " << floresta->get_individuos(i)->get_y() << endl;
                    }
                }
                floresta->realiza_acao(lowerInd);
            }
            break;
            
        case 6:
            for(unsigned int i=0; i<floresta->conta_individuos();i++)
            {
                outputSIM << floresta->tempo_do_mundo << " " << floresta->get_individuos(i)->get_id() << " " << floresta->get_individuos(i)->get_x() << " " << floresta->get_individuos(i)->get_y() << endl;
            }
            
            while (floresta->tempo_do_mundo < tempo[0] && floresta->conta_individuos() > 0)
            {
                double t_ant = floresta->tempo_do_mundo;
                int lowerInd=floresta->updateSKLOGG();
                if(t_ant < (int)floresta->tempo_do_mundo)
                {
                    for(unsigned int i=0; i<floresta->conta_individuos();i++)
                    {
                        outputSIM << (int)floresta->tempo_do_mundo << " " << floresta->get_individuos(i)->get_id() << " " << floresta->get_individuos(i)->get_x() << " " << floresta->get_individuos(i)->get_y() << endl;
                    }
                }
                floresta->realiza_acao(lowerInd);
            }
            break;
            
        case 7:
            for(unsigned int i=0; i<floresta->conta_individuos();i++)
            {
                outputSIM << floresta->tempo_do_mundo << " " << floresta->get_individuos(i)->get_id() << " " << floresta->get_individuos(i)->get_x() << " " << floresta->get_individuos(i)->get_y() << " " << floresta->get_individuos(i)->NBHood_size() << endl;
            }
            
            while (floresta->tempo_do_mundo < tempo[0] && floresta->conta_individuos() > 0)
            {
                double t_ant = floresta->tempo_do_mundo;
                int lowerInd=floresta->updateSKLOGL();
                if(t_ant < (int)floresta->tempo_do_mundo)
                {
                    for(unsigned int i=0; i<floresta->conta_individuos();i++)
                    {
                        outputSIM << (int)floresta->tempo_do_mundo << " " << floresta->get_individuos(i)->get_id() << " " << floresta->get_individuos(i)->get_x() << " " << floresta->get_individuos(i)->get_y() << " " << floresta->get_individuos(i)->NBHood_size() << endl;
                    }
                }
                floresta->doActionSKLOGL(lowerInd);
                //cout << floresta->get_individuos(lowerInd)->get_id() << endl;
            }
            break;
    }
    
    if(floresta->conta_individuos()==0)
    {
        switch (SIM_Type)
        {
            case 1:
                outputSIM << floresta->tempo_do_mundo << " " << "NA" << endl;
                break;
            
            case 2:
                outputSIM << floresta->tempo_do_mundo << " " << "NA" << endl;
            break;
        
            case 3:
                outputSIM << floresta->tempo_do_mundo << " " << "NA" << endl;
                break;
                
            case 4:
                outputSIM << floresta->tempo_do_mundo << " " << "NA" << " " << "NA" << " " << "NA" << endl;
                break;
                
            case 5:
                outputSIM << floresta->tempo_do_mundo << " " << "NA" << " " << "NA" << " " << "NA" << endl;
                break;
                
            case 6:
                outputSIM << floresta->tempo_do_mundo << " " << "NA" << " " << "NA" << " " << "NA" << endl;
                break;
                
            case 7:
                outputSIM << floresta->tempo_do_mundo << " " << "NA" << " " << "NA" << " " << "NA" << " " << "NA" << endl;
                break;
        }
    }
    outputSIM.close(); //end of output file
                
    delete floresta;
    PutRNGstate();
}
        
        
