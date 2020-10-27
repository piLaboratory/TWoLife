#include "paisagem.hpp"
#include "individuo.hpp"
#include <R.h>
#include <Rmath.h>

// Class Constructor
paisagem::paisagem(double raio,
                   int N,
                   double angulo_visada,
                   double passo,
                   double taxa_move,
                   double taxa_basal,
                   double taxa_morte,
                   double incl_b,
                   double incl_d,
                   int numb_cells,
                   double cell_size,
                   int land_shape,
                   int density_type,
                   double death_mat,
                   double move_mat,
                   int inipos,
                   int bound_condition,
                   double scape[],
                   double initialPosX[],
                   double initialPosY[],
                   double genotype_means[],
                   double width_sds[],
                   int points [],
                   bool Null):
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
            this->patches[i][j]=0; // Initializas the patches vector giving it 0 as a null value
        }
    }
    
    // Atributes values to the landscape pixels, determining the patch they belong to
    int component = 0; // Counter variable used for identyfing the patches
    for (int k = 0; k < this->numb_cells; ++k) // Passes through the rows
    for (int l = 0; l < this->numb_cells; ++l) // Passes through the columns
    if (!this->patches[k][l] && this->landscape[k][l]) find_patches(k, l, ++component); // Checks if each landscape pixel has already been assigned a patch identfication value AND if it is a non-matrix patch
    this->numb_patches = component; // Stores the number of habitat patches in the landscape
    
    this->patch_area = new double[this->numb_patches+1]; //Points to a new vector object used for storing the area of all habitats and the matrix
    
    
    for (unsigned int j = 0; j<numb_patches+1; j++) //Passes through al numb_patches
    this->patch_area[j] = 0; // Initializes the vector giving it 0 as a null value
    
    for(int i = 0; i < this->numb_cells; i++) //Passes through the rows
    for(int j = 0; j < this->numb_cells; j++) //Passes through the columns
    this->patch_area[patches[i][j]] += 1; // Indexes the patch_area vector by the patch identfication numbers of each pixel in the landscape so that it adds the number o pixels of each patch
    for (unsigned int j = 0; j<numb_patches+1; j++) //Passes through the patches
    this->patch_area[j] = this->patch_area[j]*this->cell_size*this->cell_size; // multiplies the number of pixels of each patch by the area of a pixel
    
    
    // Modification of the radius for the global density case
    if(density_type==0) // Checks if the density type is set to zero (global)
    {
        raio = this->tamanho/sqrt(M_PI); // Changes the perception radius to the equivalent of a radius of a circle with the same area of a square with the lenght of the landscape
    }
    
    // Inserts N individuals in the landscape through the populating() function
    this->populating(raio,N,angulo_visada,passo,taxa_move,taxa_basal,taxa_morte,incl_b,incl_d,death_mat,move_mat,density_type,initialPosX, initialPosY,genotype_means,width_sds,points, Null);
    
    for(unsigned int i=0; i<this->popIndividuos.size(); i++) // Passes through each individuals
    {
        this->set_vizinhos(this->popIndividuos[i]); // Calls a function of the individual/individuo class to update the neighbours of a individual (the individuals within a radius distance of the focal individual)
        this->atualiza_habitat(this->popIndividuos[i]); // Calls a function of the individual/individuo class to update the habitat type the individual is currently on
        this->atualiza_patch(this->popIndividuos[i]); // Calls a function of the individual/individuo class to update the patch identfier of the patch he is currently on
    }
    
    
    for(unsigned int i=0; i<this->popIndividuos.size(); i++) // Passes through each individual
    {
        double dsty=this->calcDensity(popIndividuos[i]); // Calls a function of the individual/individuo class to determine the density of individuals within a radius distance area from a focal individual
        this->popIndividuos[i]->update(dsty);   // Calls a function of the individual/individuo class to update the characteristics of each individual
    }
    
}

// Function responsible for calling the constructor of the individual/individuo class to create N individuals at the start of the simulation
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
    individuo::reset_id(); // Restars the id counter of the individuals
    
    vector<double> genotype, width; // Creates vectors for storing either a genotype/wisght value (normal case), or all of them (null case)
    
    if (Null==FALSE) {// Checks if this is not a Null simulation run
        
        genotype.resize(1), // Resizes de genotype vectors for containing only a value
        width.resize(1); // Resizes de width vectors for containing only a value
    }
    else{ // Checks if this  a Null simulation run
        
        for (int i=0; i<N; i++) { // Pass trough each individual
            genotype.push_back(genotype_means[i]); // Stores each genotype value of the population in the vector
            width.push_back(width_sds[i]); // Stores each width value of the population in the vector
        }
    }
    
    // Considerations of diferent types of starting positions
    if(this->initialPos==0) // Checks if the initialPos is set to 0 (origin)
    {
        for(int i=0; i<this->N; i++) // Goess through the amount of initial individuals selected
        {
            if (Null==FALSE) { // Checks if this is not a Null simulation run
                genotype[0] = genotype_means[i]; // Stores the genotype value of the current individual
                width[0] = width_sds[i]; // Stores the width value of the current individual
            }
            
            //OBS: As popAgents is a pointer of vectors, when adding variable addresses we use "new". This way should be faster, as we can acsess only the adress instead of keep storing the values
            
            // Creates an individual in the origin
            this->popIndividuos.push_back(new individuo(
                                                        0,//X Coordinate
                                                        0,//Y Coordinate
                                                        0,//Species ID
                                                        taxa_morte,//The basal death rate (The rate at which the individuals die on a habitat patch without neigbours)
                                                        runif(0,360),// Angle that an individual is curently facing
                                                        angulo_visada,//Angle used for orientation when dispersing
                                                        passo,//The Lenght distance of a dispersal event
                                                        taxa_move, // The rate at which the individuals disperse
                                                        raio,// Density dependance radius
                                                        taxa_basal,// The basal birth rate (The rate at which the individuals give birth on a habitat patch without neigbours)
                                                        99, // Randomn seed
                                                        incl_b, //The slope of the birth density dependece function
                                                        incl_d, //The slope of the death density dependece function
                                                        death_m, //Constant that indicates how many times higher the death rate should be on non-habitat pixels
                                                        move_m, // Constant that indicates how many times lower the movement rate should be on non-habitat pixels
                                                        dens_type, //Density type (0 = global, 1 = local/within a individual radius)
                                                        genotype, // Individuals genotype trait mean,
                                                        width, // Individuals genotype trait width,
                                                        points[i])); // Number of sampling points
            
        }
    }
    if(this->initialPos==1) // Checks if the initialPos is set to 1 (Random initial positions)
    {
        for(int i=0; i<this->N; i++) // Goess through the amount of initial individuals selected
        {
            if (Null==FALSE) { // Checks if this is not a Null simulation run
                genotype[0] = genotype_means[i]; // Stores the genotype value of the current individual
                width[0] = width_sds[i]; // Stores the width value of the current individual
            }
            
            // Creates an individual at random locations
            this->popIndividuos.push_back(new individuo(
                                                        runif(this->tamanho/(-2),this->tamanho/2), //X Coordinate
                                                        runif(this->tamanho/(-2),this->tamanho/2), //Y Coordinate
                                                        0,//Species ID
                                                        taxa_morte,//The basal death rate (The rate at which the individuals die on a habitat patch without neigbours)
                                                        runif(0,360),// Angle that an individual is curently facing
                                                        angulo_visada, //Angle used for orientation when dispersing,
                                                        passo,//The Lenght distance of a dispersal event
                                                        taxa_move, // The rate at which the individuals disperse
                                                        raio, // Density dependance radius
                                                        taxa_basal, // The basal birth rate (The rate at which the individuals give birth on a habitat patch without neigbours)
                                                        99, // Randomn seed
                                                        incl_b, //The slope of the birth density dependece function
                                                        incl_d, //The slope of the death density dependece function
                                                        death_m, //Constant that indicates how many times higher the death rate should be on non-habitat pixels
                                                        move_m, // Constant that indicates how many times lower the movement rate should be on non-habitat pixels
                                                        dens_type, //Density type (0 = global, 1 = local/within a individual radius)
                                                        genotype, // Individuals genotype trait mean,
                                                        width, // Individuals genotype trait width,
                                                        points[i])); // Number of sampling points
        }
    }
    if(this->initialPos==2) // Checks if the initialPos is set to 2, Random initial positions with normal distribution.
    {
        for(int i=0; i<this->N; i++) // Goess through the amount of initial individuals selected
        {
            if (Null==FALSE) { // Checks if this is not a Null simulation run
                genotype[0] = genotype_means[i]; // Stores the genotype value of the current individual
                width[0] = width_sds[i]; // Stores the width value of the current individual
            }
            
            // Creates an individual at Random initial positions with normal distribution.
            this->popIndividuos.push_back(new individuo(
                                                        rnorm(0,sqrt(taxa_move)*passo),//X Coordinate
                                                        rnorm(0,sqrt(taxa_move)*passo),//Y Coordinate
                                                        0,//Species ID
                                                        taxa_morte,//The basal death rate (The rate at which the individuals die on a habitat patch without neigbours)
                                                        runif(0,360),// Angle that an individual is curently facing
                                                        angulo_visada,//Angle used for orientation when dispersing
                                                        passo,//The Lenght distance of a dispersal event
                                                        taxa_move, // The rate at which the individuals disperse
                                                        raio,// Density dependance radius
                                                        taxa_basal,// The basal birth rate (The rate at which the individuals give birth on a habitat patch without neigbours)
                                                        99, // Randomn seed
                                                        incl_b, //The slope of the birth density dependece function
                                                        incl_d, //The slope of the death density dependece function
                                                        death_m, //Constant that indicates how many times higher the death rate should be on non-habitat pixels
                                                        move_m, // Constant that indicates how many times lower the movement rate should be on non-habitat pixels
                                                        dens_type, //Density type (0 = global, 1 = local/within a individual radius)
                                                        genotype, // Individuals genotype trait mean,
                                                        width, // Individuals genotype trait width,
                                                        points[i])); // Number of sampling points
        }
    }
    
    if(this->initialPos==3) // Checks if the initialPos is set to 3, Sets given initial positions to individuals
    {
        
        for(int i=0; i<this->N; i++) // Goess through the amount of initial individuals selected
        {
            if (Null==FALSE) { // Checks if this is not a Null simulation run
                genotype[0] = genotype_means[i]; // Stores the genotype value of the current individual
                width[0] = width_sds[i]; // Stores the width value of the current individual
            }
            
            // Creates an individual at given coordinate.
            this->popIndividuos.push_back(new individuo(
                                                        initialPosX[i],//X Coordinate
                                                        initialPosY[i],//Y Coordinatey
                                                        0,//Species ID
                                                        taxa_morte,//The basal death rate (The rate at which the individuals die on a habitat patch without neigbours)
                                                        runif(0,360), // Angle that an individual is curently facing
                                                        angulo_visada,//Angle used for orientation when dispersing
                                                        passo,//The Lenght distance of a dispersal event
                                                        taxa_move, // The rate at which the individuals disperse
                                                        raio,// Density dependance radius
                                                        taxa_basal,// The basal birth rate (The rate at which the individuals give birth on a habitat patch without neigbours)
                                                        99, // Randomn seed
                                                        incl_b, //The slope of the birth density dependece function
                                                        incl_d, //The slope of the death density dependece function
                                                        death_m, //Constant that indicates how many times higher the death rate should be on non-habitat pixels
                                                        move_m, // Constant that indicates how many times lower the movement rate should be on non-habitat pixels
                                                        dens_type, //Density type (0 = global, 1 = local/within a individual radius)
                                                        genotype,  // Individuals genotype trait mean,
                                                        width, // Individuals genotype trait width,
                                                        points[i])); // Number of sampling points
        }
    }
}

// Updater function
void paisagem::update(int acao, int ind)
{
    if(this->popIndividuos.size()>0) // checks if there are any currently alive individuals
    {
        switch (acao) // Switches between possible actions (0= death, 1= birth, 2= dispersal)
        {
            case 0: //0= death
                break;
            case 1: //1= birth
                this->atualiza_habitat(this->popIndividuos[this->popIndividuos.size()-1]); // updates in the individual the environmental value of the pixel corresponent to its current coordinate
                this->atualiza_patch(this->popIndividuos[this->popIndividuos.size()-1]); // updates the individual the patch identification value of the pixel corresponent to its current coordinate
                break;
            case 2: // 2= dispersal
                this->atualiza_habitat(this->popIndividuos[ind]); // updates the individual the environmental value of the pixel corresponent to its current coordinate
                this->atualiza_patch(this->popIndividuos[ind]); // updates the individual the patch identification value of the pixel corresponent to its current coordinate
                break;
        }
        
        for(unsigned int i=0; i<this->popIndividuos.size(); i++) //Goes through all currently alive individuals
        {
            double dsty=this->calcDensity(popIndividuos[i]); // returns the dentity of individuals within a radius distance area from each individual
            this->popIndividuos[i]->update(dsty);   // Updates the rates of each individual
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

//Function that executes the selected action. It also returns a positive boolean value if the individual disperses out of the landscape boundary.
bool paisagem::realiza_acao(int acao, int lower) //TODO : criar matriz de distancias como atributo do mundo e atualiza-la apenas quanto ao individuos afetado nesta funcao)
{
    bool emigrou = false; // Initializes a boolean variable used to record if an individual has left the boundaries of the landscape as false
    switch(acao) //(0= death, 1= birth, 2= dispersal)
    {
        case 0: //0= death
            this->atualiza_vizinhos(acao, lower);
            delete this->popIndividuos[lower]; // Deletes the object in the chosen position
            this->popIndividuos.erase(this->popIndividuos.begin()+lower); // Erases the position that object occupied in the vector of individuals
            break;
            
        case 1: // 1= birth
            individuo* chosen; // Declares a new object of the individual class
            //Novo metodo para fazer copia do individuo:
            chosen = new individuo(*this->popIndividuos[lower]); // Calls a function to create a new individuals with similar characteristics to the parent one
            this->popIndividuos.push_back(chosen); // Stores the new individual at the postition of the landscape individuals vector
            this->atualiza_vizinhos(acao, lower);
            break;
            
        case 2: // 2= birth
            
            emigrou = this->walk(lower); // Performs movement action and returs if the invidual left the landscape boundaries
            
            if(emigrou) //checks if an invidual left the landscape boundaries
                this->realiza_acao(0, lower); //Treats an emigrating individual as if it died
            else if(!emigrou) //checks if an invidual left the landscape boundaries
                this->atualiza_vizinhos(acao, lower); //Apply the boundary condition if the individual disperses out of the landscape boundary Changes the boolean value to pisitive if a individual went ou of the landscape boundaries
            break;
    }
    return emigrou; // Returs true if an invidual left the landscape boundaries
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


// Funtion to apply the desired boundary conditions when  a dispersal event leaves an individual out of the landscape boundarys
bool paisagem::apply_boundary(individuo * const ind) //const
{
    bool emigrou = false; // Initializes a boolean variable used to record if an individual has left the boundaries of the landscape as false
    double rad = (double)ind->get_raio(); // Obtains the radius of the individual that dispersed
    
    switch(this->boundary_condition) // Switches between conditions (0= absortive, 1= periodical (pacman),2= reflexive)
    {
            
        case 0: //0= absortive
            if(this->landscape_shape==0) // Checks if the landscape shape is set to 0 (circle)
            {
                if(rad*rad < (double)ind->get_x()*(double)ind->get_x()+(double)ind->get_y()*(double)ind->get_y())
                    // Checks if the individuals is out of bounds
                {
                    emigrou = true; // Sets emigration to true
                }
            }
            
            if(this->landscape_shape==1) // Checks if the landscape shape is set to 1 (square)
            {
                if((double)ind->get_x()>=this->numb_cells*this->cell_size/2 || //Checks if the x coordinate is higher than the maximum value OR
                   // obs: >= porque na paisagem quadrado as bordas mais distantes de 0 iniciariam um proximo pixel que estaria fora da paisagem. Ou teriamos que assumir que esses pixels mais extremos tenha uma área maior, o que daria um trabalho adicional para implementar uma situação irreal.
                   (double)ind->get_x()<(this->numb_cells*this->cell_size/2)*(-1)|| // if the x coordinate is lower than the minimum value OR
                   (double)ind->get_y()>this->numb_cells*this->cell_size/2 || // if the y coordinate is higher than the maximum value Or
                   (double)ind->get_y()<=(this->numb_cells*this->cell_size/-2)) // if the x coordinate is lower than the minimum value
                {
                    emigrou = true; // Sets emigration to true
                }
            }
            break;
            
        case 1: // 1= periodical
            if(ind->get_x()<(this->numb_cells*this->cell_size/2)*(-1)) //Checks if the x coordinate is higher than the maximum value
                ind->set_x(this->tamanho+ind->get_x()); // Changes the x coordinate to represent the periodical condition (left corner to right corner)
            if(ind->get_x()>=this->numb_cells*this->cell_size/2) //Checks if the x coordinate is lower than the minimum value
                ind->set_x(ind->get_x()-this->tamanho); // Changes the x coordinate to represent the periodical condition (lright corner to left corner)
            if(ind->get_y()<(this->numb_cells*this->cell_size/2)*(-1)) //Checks if the y coordinate is higher than the maximum value
                ind->set_y(this->tamanho+ind->get_y()); // Changes the y coordinate to represent the periodical condition (bottom corner to top corner)
            if(ind->get_y()>=this->numb_cells*this->cell_size/2 ) //Checks if the y coordinate is lower than the maximum value
                ind->set_y(ind->get_y()-this->tamanho); // Changes the y coordinate to represent the periodical condition (up corner to bottom corner)
            break;
            
        case 2: // 2= reflexive
            if(ind->get_x() < -this->tamanho/2) //Checks if the x coordinate is lower than the minimum value
                ind->set_x( -this->tamanho/2 + abs(this->tamanho/2 - ind->get_x()) ); // Sets the x coordinate to the boundary minus its excess lengh
            if(ind->get_x() >= this->tamanho/2) //Checks if the x coordinate is higher than the maximum value
                ind->set_x( this->tamanho/2 - abs(this->tamanho/2 - ind->get_x()) );// Sets the x coordinate to the boundary minus its excess lengh
            if(ind->get_y() < -this->tamanho/2) //Checks if the y coordinate is lower than the minimum value
                ind->set_y( -this->tamanho/2 + abs(this->tamanho/2 - ind->get_y()) );// Sets the y coordinate to the boundary minus its excess lengh
            if(ind->get_y() >= this->tamanho/2) //Checks if the y coordinate is higher than the maximum value
                ind->set_x( this->tamanho/2 - abs(this->tamanho/2 - ind->get_y()) );// Sets the y coordinate to the boundary minus its excess lengh
            break;
    }
    return emigrou; // Returns a positive boolean value if the individual emigrated
    
}

// Function that calculate the distance between two individuals
double paisagem::calcDist(const individuo* a1, const individuo* a2) const //Virou método da paisagem pois não estava conseguindo usar as propriedades privadas da paisagem nessa função (p. ex. this->boundary_condition). Tive que tirar o primeiro const.
{
    switch(this->boundary_condition) //switches between boundary conditoion cases (0= absortive, 1= periodical (pacman),2= reflexive)
    {
        case 0: // Euclidian Distance
            return sqrt(((double)a1->get_x()-(double)a2->get_x())*((double)a1->get_x()-(double)a2->get_x())+((double)a1->get_y()-(double)a2->get_y())*((double)a1->get_y()-(double)a2->get_y()));
            break;
            
        case 1: // The Euclidian distance acounting for the periodical effect (individuals on diferent edges may be close to one another)
            double x1=a1->get_x(); // Gets the x coordinate of the first individuals
            double x2=a2->get_x(); // Gets the x coordinate of the second individuals
            double y1=a1->get_y(); // Gets the y coordinate of the first individuals
            double y2=a2->get_y(); // Gets the y coordinate of the second individuals
            double dx= x1 > x2 ? x1 - x2 : x2 - x1; // choses the lower between x1-x2 e x2-x1
            double dy= y1 > y2 ? y1 - y2 : y2 - y1; //choses the lower between y1-y2 e y2-y1
            dx = dx > this->tamanho - dx ? this->tamanho - dx : dx; // choses the lower between tam-dx e dx
            dy = dy > this->tamanho - dy ? this->tamanho - dy : dy; //choses the lower between tam-dy e dy
            return sqrt(dx*dx + dy*dy); // Computes and returns the distance
            break;
    }
}

// A function to calculate de density of individuals according to density type (global or local) and considering landscape boundary effects in the calculation.
double paisagem::calcDensity(const individuo* ind1) const
{
    double density; // Creates a temporary variable for storing density
    
    if((M_PI*(ind1->get_raio()*ind1->get_raio()))>(this->tamanho*this->tamanho)){ // check if the individual radius would be higher than the landscape
        
        density=ind1->NBHood_size()/(this->tamanho*this->tamanho); // Computes the standard case of density using the total area od the world
    }
    
    else{
        density=ind1->NBHood_size()/(M_PI*(ind1->get_raio()*ind1->get_raio())); // Computes the standard case of density (when an individual radius area is completely within the landscape)
    }
    
    // Functions for local density calculation
    
    /* 1. Circular area defining a region in which denso-dependence occurs: landscape boundary effects.
     In this case, density is the number of individuals inside the circle divided by circle area.
     This is the same calculation as for global density, except by the cases in which landscape boundary affects
     the area of the circle.
     */
    
    // Condition giving the boundary effects cases
    if(ind1->get_densType()==1) // checks if the density type is set to 1 (local/ within a radius distance)
    {
        if(ind1->get_x()*ind1->get_x()>((this->tamanho/2)-ind1->get_raio())*((this->tamanho/2)-ind1->get_raio()) || // Checks if the x coordinate is within a dadius distance from the edge OR
           ind1->get_y()*ind1->get_y()>((this->tamanho/2)-ind1->get_raio())*((this->tamanho/2)-ind1->get_raio())) // if the y coordinate is within a dadius distance from the edge
        { // OBS: The absolute values retain their distance from the edge, so the compuatations from all quadrants can be done with a single set of equations
            
            // temporary objects
            double modIx = fabs(ind1->get_x()); // Stores the the absolute value of x
            double modIy = fabs(ind1->get_y()); // Stores the the absolute value of y
            double XYmax = this->tamanho/2; // Stores the lenght of quadrant side
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
    return density; // Returns the individuals curent density
}


/*
 Sempre adicione const aos argumentos de métodos quando o método não
 deve alterá-los. Previne vários erros e pode otimizar compilação
 */

// Function that updates an individuals neighbourhood list
void paisagem::set_vizinhos(individuo * const ag1) const //acessando os vizinhos dos agentes
{
    vector <individuo*> listViz; // Creates a vector of individuals to store the neighbours
    if(ag1->get_densType()==0) // Checks if the density type is set to 0 (global)
        //dens_type poderia voltar como propriedade da paisagem. Facilitariam as coisas. Como muitas propriedades e métodos deste código, elas podem ser interpretadas das duas formas (como do individuo ou como da paisagem). O que está dando confusão é que estamos fazendo um IBM, mas para algumas situações estamos querendo simular dinâmicas cujas variáveis de interesse são propriedades populacionais e não do indivíduo. Se aceito, limar o método get_densType() do individuo.h.
    {
        for(unsigned int j=0; j<popIndividuos.size(); j++) //Goes through all the individuals
        {
            individuo* ag2=this->popIndividuos[j]; //Stores a copy of each individual
            if(ag1==ag2) continue; //Except if it is the focal individual itself
            listViz.push_back(ag2); // Stores the copy of each other individual in a neighbours vector
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
            if(d<=rad) {listViz.push_back(ag2);} // Stores the copy of each other individual within a radius distance into a neighbours vector
        }
    }
    ag1->set_vizinhos(listViz); // Passes to an individual its neighbourhood list
    
}

// Function for updating individuals neighbours depending on their last action
void paisagem::atualiza_vizinhos(int acao, int ind) const //acessando os vizinhos dos agentes
{
    vector <individuo*> listViz; // Temporary object for string neighbours
    switch(acao) // Switches between actions (0= death, 1= birth, 2= dispersal)
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
    hx= (double)ind->get_x()/this->cell_size+this->numb_cells/2; // Finds the pixel X coodenate
    hy= ((double)ind->get_y()/this->cell_size)*(-1)+this->numb_cells/2; //Finds the pixel y coodenate
    ind->set_habitat(this->landscape[hx][hy]); // Sets the patch label of the individual to the label of the patch it is currently on
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
    
    // marks the current cell with its respective patch ID value
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

// Function that decides upon, and calls other functions to perform, the desired method of dispersion (Random walk or habitat selection). In the habitat selection case this function also samples x Points within the individuals "passo" radius and the indivuduals relative fitness on that location.
bool paisagem::walk(int lower){
    
    bool emigrou =FALSE; // Initilizes boolean variable that mars if the indiviual emigrated
    
    if(this->popIndividuos[lower]->get_points()>0) // Checks if the the individuals sample their step radius
    {
        double possibilitities[this->popIndividuos[lower]->get_points()][3]; // Creates a vector for storing the possible migration locations
        
        double choice, dist; // Temporary variables for storing angle and lenght
        int hx,hy; // Temporary variables for storing the possible coordinates
        
        for (int i=0; i<this->popIndividuos[lower]->get_points(); i++) // Passes by each possible location
        {
            
            choice=runif(0.0,360.0); // Samples an angle
            dist=runif(0.0,this->popIndividuos[lower]->get_passo()); // Samples a distance
            
            possibilitities[i][0]=this->popIndividuos[lower]->get_x()+cos(choice)*dist; // Stores the possible x
            possibilitities[i][1]=this->popIndividuos[lower]->get_y()+sin(choice)*dist; // Stores the possible y
            
            
            
            switch(this->boundary_condition){ // switches between baoundary conditions (0= absortive, 1= periodical (pacman),2= reflexive)
                    
                case 0: // 0= absortive
                    
                    if(possibilitities[i][0]<(this->numb_cells*this->cell_size/2)*(-1) || //Checks if the x coordinate is higher than the maximum value OR
                       possibilitities[i][0]>=this->numb_cells*this->cell_size/2 || // if the x coordinate is lower than the minimum value OR
                       possibilitities[i][1]<(this->numb_cells*this->cell_size/2)*(-1) || // if the y coordinate is higher than the maximum value Or
                       possibilitities[i][1]>=this->numb_cells*this->cell_size/2)// if the x coordinate is lower than the minimum value
                    {
                        
                        possibilitities[i][2] = -10000; // Sets the out of boundaries locations rank to a very small number
                    }
                    
                case 1: // 1= periodical (pacman)
                    if(possibilitities[i][0]<(this->numb_cells*this->cell_size/2)*(-1)) //Checks if the x coordinate is higher than the maximum value OR
                        possibilitities[i][0]=(this->tamanho+possibilitities[i][0]); // Transtports coordinate to the other edge
                    if(possibilitities[i][0]>=this->numb_cells*this->cell_size/2) // if the x coordinate is lower than the minimum value OR
                        possibilitities[i][0]=(possibilitities[i][0]-this->tamanho); // Transtports coordinate to the other edge
                    if(possibilitities[i][1]<(this->numb_cells*this->cell_size/2)*(-1))// if the y coordinate is higher than the maximum value Or
                        possibilitities[i][1]=(this->tamanho+possibilitities[i][1]); // Transtports coordinate to the other edge
                    if(possibilitities[i][1]>=this->numb_cells*this->cell_size/2 ) // if the x coordinate is lower than the minimum value
                        possibilitities[i][1]=(possibilitities[i][1]-this->tamanho); // Transtports coordinate to the other edge
                    
                    hx= (double)possibilitities[i][0]/this->cell_size+this->numb_cells/2; // Discovers the x pixel index of the coordinate
                    hy= ((double)possibilitities[i][1]/this->cell_size)*(-1)+this->numb_cells/2; // Discovers the y pixel index of the coordinate
                    
                    possibilitities[i][2] = this->landscape[hx][hy]; // Stores the habitat value of the current pixel
                    
                    break;
                    
                case 2:// 2= reflexive
                    // Prevents out of boundary corrdinates from being sampled
                    while (possibilitities[i][0]<(this->numb_cells*this->cell_size/2)*(-1) || //Checks if the x coordinate is higher than the maximum value OR
                           possibilitities[i][0]>=this->numb_cells*this->cell_size/2 || // if the x coordinate is lower than the minimum value OR
                           possibilitities[i][1]<(this->numb_cells*this->cell_size/2)*(-1) || // if the y coordinate is higher than the maximum value Or
                           possibilitities[i][1]>=this->numb_cells*this->cell_size/2) { // if the x coordinate is lower than the minimum value
                        
                        choice=runif(0.0,360.0); // re-Samples an angle
                        dist=runif(0.0,this->popIndividuos[lower]->get_passo()); // re-Samples a distance
                        
                        possibilitities[i][0]=this->popIndividuos[lower]->get_x()+cos(choice)*dist;// Stores the possible x
                        possibilitities[i][1]=this->popIndividuos[lower]->get_y()+sin(choice)*dist;// Stores the possible y
                    }
                    hx= (double)possibilitities[i][0]/this->cell_size+this->numb_cells/2; // Discovers the x pixel index of the coordinate
                    hy= ((double)possibilitities[i][1]/this->cell_size)*(-1)+this->numb_cells/2; // Discovers the y pixel index of the coordinate
                    
                    possibilitities[i][2] = this->landscape[hx][hy]; // Stores the habitat value of the current pixel
                    
                    break;
            }
            
        }
        
        this->popIndividuos[lower]->habitat_selection(possibilitities);// Passes the possible locations for the individual to chose from
        
        if (this->boundary_condition==0) {
            emigrou = this->apply_boundary(popIndividuos[lower]);
        }
        
    }
    
    else  // Checks if the the individuals dont sample their step radius (ramdom)
    {
        this->popIndividuos[lower]->anda(); // Calls a function that changes the X and Y coordinates of the individuals via randomwalk
        emigrou = this->apply_boundary(popIndividuos[lower]); // Applyes boundary condition
    }
    return emigrou;
}


