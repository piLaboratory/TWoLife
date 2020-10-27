#include "paisagem.hpp"
#include "individuo.hpp"
#include <R.h>
#include <Rmath.h>

paisagem::paisagem(double raio, int N, double angulo_visada, double passo, double taxa_move, double taxa_basal,
                   double taxa_morte, double incl_b, double incl_d, int numb_cells, double cell_size, int land_shape,
                   int density_type, double death_mat, double move_mat, int inipos, int bound_condition, double scape[],double initialPosX[],
                                                                                                                                          
                                                                                                                                          double initialPosY[],
                                                                                                                                                            double genotype_means[],
                                                                                                                                                                                 double width_sds[], int points [],
                                                                                                                                                                                                 bool Null):
  tamanho(numb_cells*cell_size),
  N(N),
  tempo_do_mundo(0),
  numb_cells(numb_cells),
  cell_size(cell_size),
  landscape_shape(land_shape),
  boundary_condition(bound_condition),
  initialPos(inipos)
{
  for(int i=0;i<this->numb_cells;i++)
  {
    for(int j=0;j<this->numb_cells;j++)
    {
      // transforma o vetor scape recebido do construtor em uma matriz
      this->landscape[i][j]=scape[j*numb_cells+i];
      this->patches[i][j]=0;
    }
  }
  
  //Atribui valores à matriz patches, que determina o fragmento a que cada pixel pertence
  int component = 0;
  for (int k = 0; k < this->numb_cells; ++k)
    for (int l = 0; l < this->numb_cells; ++l)
      if (!this->patches[k][l] && this->landscape[k][l]) find_patches(k, l, ++component);
      this->numb_patches = component;
      
      this->patch_area = new double[this->numb_patches+1];
      
      for (unsigned int j = 0; j<numb_patches+1; j++)
        this->patch_area[j] = 0;
      
      for(int i = 0; i < this->numb_cells; i++)
        for(int j = 0; j < this->numb_cells; j++)
          this->patch_area[patches[i][j]] += 1;
      for (unsigned int j = 0; j<numb_patches+1; j++)
        this->patch_area[j] = this->patch_area[j]*this->cell_size*this->cell_size;
      
      // Calculo do raio dependendo do tipo de densidade. 0 = global, 1 = local (raio), 2 = kernel.
      if(density_type==0)
      {
        raio = this->tamanho/sqrt(M_PI);
      }
      /* Coloca os indivíduos na paisagem por meio da função populating() */
      this->populating(raio,N,angulo_visada,passo,taxa_move,taxa_basal,taxa_morte,incl_b,incl_d,death_mat,move_mat,density_type,initialPosX, initialPosY,genotype_means,width_sds,points, Null);
      
      for(unsigned int i=0; i<this->popIndividuos.size(); i++)
      {
        this->set_vizinhos(this->popIndividuos[i]);//atualiza os vizinhos
        this->atualiza_habitat(this->popIndividuos[i]);//retorna o tipo de habitat
        this->atualiza_patch(this->popIndividuos[i]);
      }
      
      
      for(unsigned int i=0; i<this->popIndividuos.size(); i++)
      {
        double dsty=this->calcDensity(popIndividuos[i]);
        this->popIndividuos[i]->update(dsty);   //e atualiza o individuo i da populacao
      }
      
}

void paisagem::populating(double raio,
                          int N,
                          double angulo_visada,
                          double passo,
                          double taxa_move,
                          double taxa_basal,
                          double taxa_morte,
                          double incl_b,
                          double incl_d,
                          double death_m,
                          double move_m,
                          int dens_type,
                          double initialPosX[],
                          double initialPosY[],
                          double genotype_means[],
                          double width_sds[],
                          int points[],                
                          bool Null)
{
  individuo::reset_id(); // reinicia o contador de id dos individuos
  // Considerar diferentes possibilidades de posições iniciais. TBI.
  
  vector<double> genotype, width;
  
  if (Null==FALSE) {
    
    genotype.resize(1),
    width.resize(1);
  }
  else{
    
    for (int i=0; i<N; i++) {
      genotype.push_back(genotype_means[i]);
      width.push_back(width_sds[i]);
    }
  }
  
  if(this->initialPos==0)
  {
    for(int i=0; i<this->N; i++)
    {
      if (Null==FALSE) {
        genotype[0] = genotype_means[i];
        width[0] = width_sds[i];
      }
      
      this->popIndividuos.push_back(new individuo(
          0,//posicao x
          0,//posicao y
          0,//especie
          taxa_morte,//taxa de morte
          runif(0,360),// orientacao
          angulo_visada,//angulo de visada
          passo,//tamanho do passo
          taxa_move, //taxa de movimentacao
          raio,//tamanho do raio
          taxa_basal,// taxa máxima de nascimento
          99, // semente de numero aleatorio
          incl_b,
          incl_d,
          death_m,
          move_m,
          dens_type,
          genotype,
          width,
          points[i]));
      // como o popAgentes eh um ponteiro de vetores, ao adicionar enderecos das variaveis, usamos os new. Dessa forma fica mais rapido
      //pois podemos acessar apenas o endereco e nao ficar guardando todos os valores
    }
  }
  if(this->initialPos==1) // Random initial positions (initialPos==1)
  {
    for(int i=0; i<this->N; i++)
    {
      if (Null==FALSE) {
        genotype[0] = genotype_means[i];
        width[0] = width_sds[i];
      }
      
      this->popIndividuos.push_back(new individuo(
          runif(this->tamanho/(-2),this->tamanho/2), //posicao x
          runif(this->tamanho/(-2),this->tamanho/2), //posicao y
          0,//especie
          taxa_morte,//taxa de morte
          runif(0,360),// orientacao
          angulo_visada,//angulo de visada
          passo,//tamanho do passo
          taxa_move, //taxa de movimentacao
          raio,//tamanho do raio
          taxa_basal,// taxa máxima de nascimento
          99, // semente de numero aleatorio
          incl_b,
          incl_d,
          death_m,
          move_m,
          dens_type,
          genotype,
          width,
          points[i]));
    }
  }
  if(this->initialPos==2) // Random initial positions with normal distribution. TBI: tornar os parametros da rnorm livres
  {
    for(int i=0; i<this->N; i++)
    {
      if (Null==FALSE) {
        genotype[0] = genotype_means[i];
        width[0] = width_sds[i];
      }
      
      this->popIndividuos.push_back(new individuo(
          rnorm(0,sqrt(taxa_move)*passo),//posicao x
          rnorm(0,sqrt(taxa_move)*passo),//posicao y
          0,//especie
          taxa_morte,//taxa de morte
          runif(0,360),// orientacao
          angulo_visada,//angulo de visada
          passo,//tamanho do passo
          taxa_move, //taxa de movimentacao
          raio,//tamanho do raio
          taxa_basal,// taxa máxima de nascimento
          99, // semente de numero aleatorio
          incl_b,
          incl_d,
          death_m,
          move_m,
          dens_type,
          genotype,
          width,
          points[i]));
    }
  }
  
  if(this->initialPos==3) // Sets given initial positions to individuals
  {
    
    for(int i=0; i<this->N; i++)
    {
      if (Null==FALSE) {
        genotype[0] = genotype_means[i];
        width[0] = width_sds[i];
      }
      
      
      this->popIndividuos.push_back(new individuo(
          initialPosX[i],//posicao x
                     initialPosY[i],//posicao y
                                0,//especie
                                taxa_morte,//taxa de morte
                                runif(0,360),// orientacao
                                angulo_visada,//angulo de visada
                                passo,//tamanho do passo
                                taxa_move, //taxa de movimentacao
                                raio,//tamanho do raio
                                taxa_basal,// taxa máxima de nascimento
                                99, // semente de numero aleatorio
                                incl_b,
                                incl_d,
                                death_m,
                                move_m,
                                dens_type,
                                genotype,
                                width,
                                points[i]));
    }
  }
}

void paisagem::update(int acao, int ind)
{
  if(this->popIndividuos.size()>0)
  {
    switch (acao)
    {
    case 0:
      break;
    case 1:
      this->atualiza_habitat(this->popIndividuos[this->popIndividuos.size()-1]);
      this->atualiza_patch(this->popIndividuos[this->popIndividuos.size()-1]);
      break;
    case 2:
      this->atualiza_habitat(this->popIndividuos[ind]);
      this->atualiza_patch(this->popIndividuos[ind]);
      break;
    }
    
    for(unsigned int i=0; i<this->popIndividuos.size(); i++)
    {
      double dsty=this->calcDensity(popIndividuos[i]);
      this->popIndividuos[i]->update(dsty);   //e atualiza o individuo i da populacao
    }
  }
}

int paisagem::sorteia_individuo()
{
  // time for next event and simulation time update
  int menor=0;
  double menor_tempo = this->popIndividuos[0]->get_tempo();
  
  for(unsigned int i=1; i<this->popIndividuos.size(); i++)
  {
    if(this->popIndividuos[i]->get_tempo()<menor_tempo)
    {
      menor = i;
      menor_tempo = this->popIndividuos[i]->get_tempo();
    }
  }
  return menor;
}

bool paisagem::realiza_acao(int acao, int lower) //TODO : criar matriz de distancias como atributo do mundo e atualiza-la apenas quanto ao individuos afetado nesta funcao)
{
  bool emigrou = false;
  switch(acao) //0 eh morte, 1 eh nascer, 2 eh andar
  {
  case 0:
    this->atualiza_vizinhos(acao, lower);
    delete this->popIndividuos[lower];
    this->popIndividuos.erase(this->popIndividuos.begin()+lower);
    break;
    
  case 1:
    individuo* chosen;
    //Novo metodo para fazer copia do individuo:
    chosen = new individuo(*this->popIndividuos[lower]);
    this->popIndividuos.push_back(chosen);
    this->atualiza_vizinhos(acao, lower);
    break;
    
  case 2:
    
    emigrou = this->walk(lower);
    
    if(emigrou)
      this->realiza_acao(0, lower);
    else if(!emigrou)
      this->atualiza_vizinhos(acao, lower);
    break;
  }
  return emigrou;
}

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


// metodo para condicao de contorno, argumento é um ponteiro para um individuo
//TODO: conferir se a combinacao x , y da condicao esta gerando o efeito desejado
//TBI: condicao periodica do codigo antigo feito com Garcia. Verificar se estah correta FE:Verifiquei e tinha um problema
// (veja p. ex. um unico individuo apenas se movimentando)
bool paisagem::apply_boundary(individuo * const ind) //const
{
  bool emigrou = false;
  double rad = (double)ind->get_raio();
  switch(this->boundary_condition)
  {
    
  case 0:
    if(this->landscape_shape==0)
    {
      if(rad*rad < (double)ind->get_x()*(double)ind->get_x()+(double)ind->get_y()*(double)ind->get_y())
      {
        emigrou = true;
      }
    }
    
    if(this->landscape_shape==1)
    {
      if((double)ind->get_x()>=this->numb_cells*this->cell_size/2 || //>= porque na paisagem quadrado as bordas mais distantes de 0 iniciariam um proximo pixel que estaria fora da paisagem. Ou teriamos que assumir que esses pixels mais extremos tenha uma área maior, o que daria um trabalho adicional para implementar uma situação irreal.
         (double)ind->get_x()<(this->numb_cells*this->cell_size/2)*(-1)||
         (double)ind->get_y()>this->numb_cells*this->cell_size/2 ||
         (double)ind->get_y()<=(this->numb_cells*this->cell_size/-2))
      {
        emigrou = true;
      }
    }
    break;
    
  case 1:
    if(ind->get_x()<(this->numb_cells*this->cell_size/2)*(-1))
      ind->set_x(this->tamanho+ind->get_x());
    if(ind->get_x()>=this->numb_cells*this->cell_size/2)
      ind->set_x(ind->get_x()-this->tamanho);
    if(ind->get_y()<(this->numb_cells*this->cell_size/2)*(-1))
      ind->set_y(this->tamanho+ind->get_y());
    if(ind->get_y()>=this->numb_cells*this->cell_size/2 )
      ind->set_y(ind->get_y()-this->tamanho);
    break;
    
  case 2:
    if(ind->get_x() < -this->tamanho/2)
      ind->set_x( -this->tamanho/2 + abs(this->tamanho/2 - ind->get_x()) );
    if(ind->get_x() >= this->tamanho/2)
      ind->set_x( this->tamanho/2 - abs(this->tamanho/2 - ind->get_x()) );
    if(ind->get_y() < -this->tamanho/2)
      ind->set_y( -this->tamanho/2 + abs(this->tamanho/2 - ind->get_y()) );
    if(ind->get_y() >= this->tamanho/2)
      ind->set_x( this->tamanho/2 - abs(this->tamanho/2 - ind->get_y()) );
    break;
  }
  return emigrou;
  
}


double paisagem::calcDist(const individuo* a1, const individuo* a2) const //Virou método da paisagem pois não estava conseguindo usar as propriedades privadas da paisagem nessa função (p. ex. this->boundary_condition). Tive que tirar o primeiro const.
{
  switch(this->boundary_condition)
  {
  case 0:
    return sqrt(((double)a1->get_x()-(double)a2->get_x())*((double)a1->get_x()-(double)a2->get_x())+((double)a1->get_y()-(double)a2->get_y())*((double)a1->get_y()-(double)a2->get_y()));
    break;
    
  case 1:
    double x1=a1->get_x();
    double x2=a2->get_x();
    double y1=a1->get_y();
    double y2=a2->get_y();
    double dx= x1 > x2 ? x1 - x2 : x2 - x1; // escolhe o menor entre x1-x2 e x2-x1
    double dy= y1 > y2 ? y1 - y2 : y2 - y1;
    dx = dx > this->tamanho - dx ? this->tamanho - dx : dx; // escolhe o menor entre tam-dx e dx
    dy = dy > this->tamanho - dy ? this->tamanho - dy : dy;
    return sqrt(dx*dx + dy*dy);
    break;
  }
}

// A function to calculate de density of individuals according to density type (global or local) and considering landscape boundary effects in the calculation.
double paisagem::calcDensity(const individuo* ind1) const
{
  double density;
  
  if((M_PI*(ind1->get_raio()*ind1->get_raio()))>(this->tamanho*this->tamanho)){
    
    density=ind1->NBHood_size()/(this->tamanho*this->tamanho);
  }
  
  else{
    density=ind1->NBHood_size()/(M_PI*(ind1->get_raio()*ind1->get_raio()));
  }
  
  // Functions for local density calculation
  
  /* 1. Circular area defining a region in which denso-dependence occurs: landscape boundary effects.
   In this case, density is the number of individuals inside the circle divided by circle area.
   This is the same calculation as for global density, except by the cases in which landscape boundary affects
   the area of the circle.
   */
  
  // Condition giving the boundary effects cases
  if(ind1->get_densType()==1)
  {
    if(ind1->get_x()*ind1->get_x()>((this->tamanho/2)-ind1->get_raio())*((this->tamanho/2)-ind1->get_raio()) ||
       ind1->get_y()*ind1->get_y()>((this->tamanho/2)-ind1->get_raio())*((this->tamanho/2)-ind1->get_raio()))
    {
      // temporary objects
      double modIx = fabs(ind1->get_x());
      double modIy = fabs(ind1->get_y());
      double XYmax = this->tamanho/2;
      vector<double>secX;
      vector<double>secY;
      
      // Functions for adjusted local density calculation, according to the specific case
      // 1)
      if(modIx>XYmax-ind1->get_raio() && modIy<=XYmax-ind1->get_raio())
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
      if(modIx<=XYmax-ind1->get_raio() && modIy>XYmax-ind1->get_raio())
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

void paisagem::set_vizinhos(individuo * const ag1) const //acessando os vizinhos dos agentes
{
  vector <individuo*> listViz;
  if(ag1->get_densType()==0) //dens_type poderia voltar como propriedade da paisagem. Facilitariam as coisas. Como muitas propriedades e métodos deste código, elas podem ser interpretadas das duas formas (como do individuo ou como da paisagem). O que está dando confusão é que estamos fazendo um IBM, mas para algumas situações estamos querendo simular dinâmicas cujas variáveis de interesse são propriedades populacionais e não do indivíduo. Se aceito, limar o método get_densType() do individuo.h.
  {
    for(unsigned int j=0; j<popIndividuos.size(); j++)
    {
      individuo* ag2=this->popIndividuos[j];
      if(ag1==ag2) continue;
      listViz.push_back(ag2);
    }
  }
  else
  {
    double rad = (double)ag1->get_raio();
    for(unsigned int j=0; j<popIndividuos.size(); j++)
    {
      individuo* ag2=this->popIndividuos[j];
      if(ag1==ag2) continue;
      double d=this->calcDist(ag1,ag2);
      if(d<=rad) {listViz.push_back(ag2);}
    }
  }
  ag1->set_vizinhos(listViz);
  
}

void paisagem::atualiza_vizinhos(int acao, int ind) const //acessando os vizinhos dos agentes
{
  vector <individuo*> listViz;
  switch(acao)
  {
  case 0: // Performs the neighbourhood update in case of a death event
    for(unsigned int i=0; i<this->popIndividuos[ind]->NBHood_size()-1; i++) // Passes through each of the drafted individuals neighbours
    {
      for(unsigned int j=0; j<this->popIndividuos[ind]->lisViz[i]->NBHood_size()-1; j++)// Passes through each of the drafted individuals neighbours neighbours
      {
        if(this->popIndividuos[ind]->lisViz[i]->lisViz[j]->get_id() == this->popIndividuos[ind]->get_id()) // Checks if the current neighbours neighbour has the id of the drafted individual
        {
          this->popIndividuos[ind]->lisViz[i]->lisViz.erase(this->popIndividuos[ind]->lisViz[i]->lisViz.begin()+j); // Erases the drafted individual from the neighbourhood of its neighbours
        }
      }
    }
    break; // Breaks switch
    
  case 1:// Performs the neighbourhood update in case of a birth event
    for(unsigned int i=0; i<this->popIndividuos[ind]->NBHood_size()-1; i++)// Passes through each of the drafted individuals neighbours
    {
      this->popIndividuos[ind]->lisViz[i]->lisViz.push_back(this->popIndividuos[this->popIndividuos.size()-1]);// Adds the last born individual to its neighbours neighbourhood
      listViz.push_back(this->popIndividuos[ind]->lisViz[i]);// Stores the drafted individuals neighbours
    }
    
    this->popIndividuos[ind]->lisViz.push_back(this->popIndividuos[this->popIndividuos.size()-1]);// Adds the last born individual to its parent neighbourhood
    listViz.push_back(this->popIndividuos[ind]);// Adds the last born individual to its own neighbourhood
    this->popIndividuos[this->popIndividuos.size()-1]->set_vizinhos(listViz);// Calls a funtion that stores the last born individuals neighbours in its neighbourhood
    break; // Breaks switch
    
  case 2:// Performs the neighbourhood update in case of a migration event
    
    if(this->popIndividuos[ind]->get_densType()!=0){// Checks if the density type is not set to global
      
      for(unsigned int i=0; i<this->popIndividuos[ind]->NBHood_size()-1; i++) // Passes through each of the drafted individuals neighbours
      {
        for(unsigned int j=0; j<this->popIndividuos[ind]->lisViz[i]->NBHood_size()-1; j++)// Passes through each of the drafted individuals neighbours neighbours
        {
          if(this->popIndividuos[ind]->lisViz[i]->lisViz[j]->get_id() == this->popIndividuos[ind]->get_id()) // Checks if the current neighbours neighbour has the id of the drafted individual
          {
            this->popIndividuos[ind]->lisViz[i]->lisViz.erase(this->popIndividuos[ind]->lisViz[i]->lisViz.begin()+j); // Erases the drafted individual from the neighbourhood of its neighbours
          }
        }
      }
      
      double rad = (double) this->popIndividuos[ind]->get_raio(); // Obtains the drafted individuals radius
      for(unsigned int k=0; k<this->popIndividuos.size(); k++) // Passes through each living individual
      {
        if(k==ind) // Checks if the compared individual isnt itself
          continue;
        if(this->calcDist(this->popIndividuos[ind], this->popIndividuos[k]) <= rad) // Checks if the compared individual is within a radius distance
        {
          listViz.push_back(this->popIndividuos[k]); // Stores the drafted individuals neighbours
          this->popIndividuos[k]->lisViz.push_back(this->popIndividuos[ind]); // Stores the drafted individual on its neighbours neighbourhood
        }
      }
      this->popIndividuos[ind]->set_vizinhos(listViz); // Calls a funtion that stores the neighbours in the individuals neighbourhood
      break; // Breaks switch
    
    }
  }
}

void paisagem::atualiza_habitat(individuo * const ind) const
{
  // Tinha um IF com landscape_shape que eu removi. Não entendi como a paisagem ser circular
  // interfere em ser habitat ou não: isso deve interferir na apply_boundary apenas, certo?
  // Also: Tinha uma inversão do y que eu também não entendi e removi
  // A.C. 10.07.13
  
  // Um termo (-1) foi removido erroneamente por A.C.. Para o hy, o sentido em que o número de células aumenta é o
  //contrário do sentido em que as coordenadas aumentam. Portanto a multiplicação por - 1 é necessária.
  // M.A. 12.09.14
  int hx,hy;
  hx= (double)ind->get_x()/this->cell_size+this->numb_cells/2;
  hy= ((double)ind->get_y()/this->cell_size)*(-1)+this->numb_cells/2;
  ind->set_habitat(this->landscape[hx][hy]);
}

void paisagem::atualiza_patch(individuo * const ind) const
{
  int hx,hy;
  hx= (double)ind->get_x()/this->cell_size+this->numb_cells/2;
  hy= ((double)ind->get_y()/this->cell_size)*(-1)+this->numb_cells/2;
  ind->set_patch(this->patches[hx][hy]);
}

void paisagem::find_patches(int x, int y, int current_label)
{
  if (x < 0 || x == this->numb_cells) return; // out of bounds
  if (y < 0 || y == this->numb_cells) return; // out of bounds
  if (this->patches[x][y] || !this->landscape[x][y]) return; // already labeled or not marked with 1 in m
  
  // mark the current cell
  this->patches[x][y] = current_label;
  
  // recursively mark the neighbors
  find_patches(x + 1, y, current_label);
  find_patches(x, y + 1, current_label);
  find_patches(x - 1, y, current_label);
  find_patches(x, y - 1, current_label);
  
  // 8-rule
  find_patches(x + 1, y + 1, current_label);
  find_patches(x - 1, y + 1, current_label);
  find_patches(x - 1, y - 1, current_label);
  find_patches(x + 1, y - 1, current_label);
}

bool paisagem::walk(int lower){
  
  bool emigrou =FALSE;
  
  if(this->popIndividuos[lower]->get_points()>0)
  {
    double possibilitities[this->popIndividuos[lower]->get_points()][3];
    
    double choice, dist;
    int hx,hy;
    
    for (int i=0; i<this->popIndividuos[lower]->get_points(); i++)
    {
      
      choice=runif(0.0,360.0);
      dist=runif(0.0,this->popIndividuos[lower]->get_passo());
      
      possibilitities[i][0]=this->popIndividuos[lower]->get_x()+cos(choice)*dist;
      possibilitities[i][1]=this->popIndividuos[lower]->get_y()+sin(choice)*dist;
      
      
      
      switch(this->boundary_condition){
      
      case 0:
        
        if(possibilitities[i][0]<(this->numb_cells*this->cell_size/2)*(-1) || possibilitities[i][0]>=this->numb_cells*this->cell_size/2 || possibilitities[i][1]<(this->numb_cells*this->cell_size/2)*(-1) || possibilitities[i][1]>=this->numb_cells*this->cell_size/2)
        {
          
          possibilitities[i][2] = -10000;
        }
        
      case 1:
        if(possibilitities[i][0]<(this->numb_cells*this->cell_size/2)*(-1))
          possibilitities[i][0]=(this->tamanho+possibilitities[i][0]);
        if(possibilitities[i][0]>=this->numb_cells*this->cell_size/2)
          possibilitities[i][0]=(possibilitities[i][0]-this->tamanho);
        if(possibilitities[i][1]<(this->numb_cells*this->cell_size/2)*(-1))
          possibilitities[i][1]=(this->tamanho+possibilitities[i][1]);
        if(possibilitities[i][1]>=this->numb_cells*this->cell_size/2 )
          possibilitities[i][1]=(possibilitities[i][1]-this->tamanho);
        
        hx= (double)possibilitities[i][0]/this->cell_size+this->numb_cells/2;
        hy= ((double)possibilitities[i][1]/this->cell_size)*(-1)+this->numb_cells/2;
        
        possibilitities[i][2] = this->landscape[hx][hy];
        
        break;
        
      case 2:
        
        while (possibilitities[i][0]<(this->numb_cells*this->cell_size/2)*(-1) || possibilitities[i][0]>=this->numb_cells*this->cell_size/2 || possibilitities[i][1]<(this->numb_cells*this->cell_size/2)*(-1) || possibilitities[i][1]>=this->numb_cells*this->cell_size/2) {
          
          choice=runif(0.0,360.0);
          dist=runif(0.0,this->popIndividuos[lower]->get_passo());
          
          possibilitities[i][0]=this->popIndividuos[lower]->get_x()+cos(choice)*dist;
          possibilitities[i][1]=this->popIndividuos[lower]->get_y()+sin(choice)*dist;
        }
        hx= (double)possibilitities[i][0]/this->cell_size+this->numb_cells/2;
        hy= ((double)possibilitities[i][1]/this->cell_size)*(-1)+this->numb_cells/2;
        
        possibilitities[i][2] = this->landscape[hx][hy];
        
        break;
      }
      
    }
    
    this->popIndividuos[lower]->habitat_selection(possibilitities);
    
    if (this->boundary_condition==0) {
      emigrou = this->apply_boundary(popIndividuos[lower]);
    }
    
  }
  
  else
  {
    this->popIndividuos[lower]->anda(); // Calls a function that changes the X and Y coordinates of the individuals
    emigrou = this->apply_boundary(popIndividuos[lower]);
  }
  return emigrou;
}


