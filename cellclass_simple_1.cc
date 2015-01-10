#include "cellclass.h"
#include <cmath>
#include <cfloat>
#include <fstream>
//#include "/Users/dallon/Desktop/Files/Work/Research/cadherin/randomc/randomc.h"
#include "randomc.h"
#include "stocc.h" 

//double degreeG=40;//this changes the wedge where the cadherins can go
//double degreeB=40;//this changes the wedge where the cadherins can go

// Create and initialize the random number generator
CRandomMersenne rg(12335);
StochasticLib2 sto(14992);
//TRandomMersenne rg(1233);
int index_square1[9][2]={0,0,-1,-1,-1,0,-1,1,0,-1,0,1,1,-1,1,1,1,0};


// definition of private functions

void Cell::calculate_cadherin_length(int cadherin_index){
  calculate_cadherin_vector(cadherin_index);
	cadherin[cadherin_index].length = sqrt(cadherin[cadherin_index].vector[0]*cadherin[cadherin_index].vector[0] + cadherin[cadherin_index].vector[1]*cadherin[cadherin_index].vector[1]);
	//cout<<cadherin[cadherin_index].length<<" "<<cadherin[cadherin_index].vector[0]<<" "<<cadherin[cadherin_index].vector[1];
}

void Cell::calculate_cadherin_force(int cadherin_index){
}

void Cell::calculate_cadherin_vector(int cadherin_index){
  cadherin[cadherin_index].vector[0] = cadherin[cadherin_index].location[0]-center[0];
  cadherin[cadherin_index].vector[1] = cadherin[cadherin_index].location[1]-center[1];
}

void Cell::calculate_cadherin_front_indicator(int cadherin_index){
	double tmp;
	tmp = cadherin[cadherin_index].vector[0]*front[0]+ cadherin[cadherin_index].vector[1]*front[1];//front is a unit vector, [1,0]
	cadherin[cadherin_index].front_indicator = 0;// front	
	if( tmp < 0 ){
		cadherin[cadherin_index].front_indicator = 1;// back		
		}		
}

void Cell::calculate_integrin_length(int integrin_index){
  integrin[integrin_index].length = sqrt(integrin[integrin_index].vector[0]*integrin[integrin_index].vector[0] + integrin[integrin_index].vector[1]*integrin[integrin_index].vector[1]);
}

void Cell::calculate_integrin_force(int integrin_index){
}

void Cell::calculate_integrin_vector(int integrin_index){
  integrin[integrin_index].vector[0] = integrin[integrin_index].location[0]-center[0];
  integrin[integrin_index].vector[1] = integrin[integrin_index].location[1]-center[1];
}

void Cell::calculate_integrin_front_indicator(int integrin_index){
  double tmp;
  tmp = integrin[integrin_index].vector[0]*front[0]+ integrin[integrin_index].vector[1]*front[1];
  integrin[integrin_index].front_indicator = 0;// front
  if( tmp < 0 )
      integrin[integrin_index].front_indicator = 1;// back
}


// definition of public functions

// Initialize the cell variables since I am a lame c++ programer
void Cell::initialize( double xc, double yc, Cadherin init_cad[], Integrin init_int[],double frontx, double fronty, int type){
  center[0] = xc;
  center[1] = yc;
  front[0] = frontx;
  front[1] = fronty;
  celltype = type;
for( int i = 0; i < ncad; i++){  
  cadherin[i] = init_cad[i];
  calculate_cadherin_vector(i);
  calculate_cadherin_front_indicator(i);
 }
for( int i = 0; i < nint; i++){  
  integrin[i] = init_int[i];
 }
}

// Constructor
Cell::Cell(){
}
Cell::Cell(double x, double y, Cadherin initial_cadherins[], Integrin initial_integrins[], double front1, double front2, int type1){
  center[0] = x;
  center[1] = y;
  front[0] = front1;
  front[1] = front2;
  celltype = type1;
for( int i = 0; i < ncad; i++){
  cadherin[i] = initial_cadherins[i];
  calculate_cadherin_vector(i);
  calculate_cadherin_front_indicator(i);
 }
for( int i = 0; i < nint; i++){  
  integrin[i] = initial_integrins[i];
 }
}
// Destructor
Cell::~Cell(){
}

//gets

void Cell::get_center(double value[], int lim){
 for (int i=0; i < lim; i++){
    value[i] = center[i];
  }
}
void Cell::get_front(double value[], int lim){
 for (int i=0; i < lim; i++){
    value[i] = front[i];
  }
}
int Cell::get_type(int dummy){
  int value;
  value = celltype; 
  return value; 
}

// get cadherin stuff
void Cell::get_cadherin_location(int cadherin_index, double value[], int lim){
  for (int i=0; i < lim; i++){
    value[i] = cadherin[cadherin_index].location[i];
  }
}
double Cell::get_cadherin_length(int cadherin_index){
  double value; 	
  calculate_cadherin_length(cadherin_index);
  value = cadherin[cadherin_index].length;
  //cout << value<<"\n";
  return value;
	
}
int Cell::get_ncadherin(int dummy){
  int value;
  value = ncadherin;
  return value;
}
double Cell::get_cadherin_force(int cadherin_index){
  double value;
  value = cadherin[cadherin_index].force;
  return value;
}
double Cell::get_cadherin_kspring(int cadherin_index, int attach){
  double value;
  value = cadherin[cadherin_index].kspring[attach];
  return value;
}
double Cell::get_cadherin_mu(int cadherin_index){
  double value;
  value = cadherin[cadherin_index].mu;
  return value;
}
int Cell::get_cadherin_attach(int cadherin_index){
  int value;
  value = cadherin[cadherin_index].attach;
  return value;
}
double Cell::get_cadherin_time(int cadherin_index){
  double value;
  value = cadherin[cadherin_index].time;
  return value;
}

void Cell::get_cadherin_vector(int cadherin_index, double value[], int lim){
  for (int i=0; i < lim; i++){
    value[i] = cadherin[cadherin_index].vector[i];
  }
}
int Cell::get_cadherin_front_indicator(int cadherin_index){
  int value;
calculate_cadherin_front_indicator(cadherin_index);
  value = cadherin[cadherin_index].front_indicator;
  return value;
}
//get integrin stuff
void Cell::get_integrin_location(int integrin_index, double value[], int lim){
  for (int i=0; i < lim; i++){
    value[i] = integrin[integrin_index].location[i];
  }
}
double Cell::get_integrin_length(int integrin_index){
  double value;
  value = integrin[integrin_index].length;
  return value;
}
int Cell::get_nintegrin(int dummy){
  int value;
  value = nintegrin;
  return value;
}
double Cell::get_integrin_force(int integrin_index){
  double value;
  value = integrin[integrin_index].force;
  return value;
}
double Cell::get_integrin_kspring(int integrin_index){
  double value;
  value = integrin[integrin_index].kspring;
  return value;
}
double Cell::get_integrin_mu(int integrin_index){
  double value;
  value = integrin[integrin_index].mu;
  return value;
}
int Cell::get_integrin_attach(int integrin_index){
  int value;
  value = integrin[integrin_index].attach;
  return value;
}
double Cell::get_integrin_time(int integrin_index){
  double value;
  value = integrin[integrin_index].time;
  return value;
}
void Cell::get_integrin_vector(int integrin_index, double value[], int lim){
  for (int i=0; i < lim; i++){
    value[i] = integrin[integrin_index].vector[i];
  }
}
int Cell::get_integrin_front_indicator(int integrin_index){
  int value;
  value = integrin[integrin_index].front_indicator;
  return value;
}
//sets
void Cell::set_center(double value[], int lim){
  for( int i=0; i<lim; i++){
    center[i]=value[i];
  }
}
void Cell::set_front(double value[], int lim){
  for( int i=0; i<lim; i++){
    front[i]=value[i];
  }
}

void Cell::set_type(int type){
	celltype=type; 
}
//set cadherins
void Cell::set_cadherin_attach(int cadherin_index, int value){
  cadherin[cadherin_index].attach = value; //-1 means free, cellnumber*ncad+cadherin_index.
}
void Cell::set_ncadherin(int value){
  ncadherin = value;
}
void Cell::set_cadherin_location(int cadherin_index, double value[], int lim){
  for (int i=0; i < lim; i++){
    cadherin[cadherin_index].location[i] = value[i];
  }
}
void Cell::set_cadherin_time(int cadherin_index, double value){
  cadherin[cadherin_index].time = value; // gives the detach time
}
//set integrins
void Cell::set_integrin_attach(int integrin_index, int value){
  integrin[integrin_index].attach = value; //-1 means free,0 means attached to substrate, cellnumber*nint+integrin_index.
}
void Cell::set_nintegrin(int value){
  nintegrin = value;
}
void Cell::set_integrin_location(int integrin_index, double value[], int lim){
  for (int i=0; i < lim; i++){
    integrin[integrin_index].location[i] = value[i];
  }
}
void Cell::set_integrin_time(int integrin_index, double value){
  integrin[integrin_index].time = value; // gives the attach or detach time
}
//update cadherins

void Cell::update_cadherin_attachment_initial(int cadherin_index, double dt, double average_spore, double average_stalk[], double average_different, int restart, double average_unbound_cad_spore, double average_unbound_cad_stalk, int cell_index, Cell *cell){
  Cadherin& cur_int = cadherin[cadherin_index];
  double newx =   cur_int.location[0];
  double newy =   cur_int.location[1];
  //Set a certain portion to stay unattached.
  if(rg.Random()<prob_attach_initial){// decide cadherin binds
    //attach the cadherin
    int i = static_cast<int> (newx/(lng_scale));
    int j = static_cast<int> (newy/(lng_scale));
    // Check the first grid then randomly loop through the 8 grid
    // squares until a cell is located
    // Randomly loop through the values of k
    int kindex[lgrd];
    int jindex[9];
    // initialize  kindex
    for (int ii = 0; ii<lgrd; ii++)
	    kindex[ii]= ii;
	  // randomize  kindex	
    for( int ii = lgrd-1; ii>=1; ii--){
    	int jj = rg.IRandom(0,ii);
    	int tmp = kindex[ii];
    	kindex[ii] = kindex[jj];
    	kindex[jj] = tmp;
    } 
    // initialize  jindex
    for (int ii = 0; ii<9; ii++)
	    jindex[ii]= ii;
    // randomize  jindex	
    for( int ii = 8; ii>=2; ii--){
    	int jj = rg.IRandom(1,ii) ;
    	int tmp = jindex[ii];
    	jindex[ii] = jindex[jj];
    	jindex[jj] = tmp;
    } 

      
    int tmpi;
    int tmpj;
    for( int jj = 0; jj<9; jj++){
	    tmpi = i+index_square1[jindex[jj]][0];

  	  if (tmpi > lnx-1)
  	    tmpi = lnx -1;
  	  if (tmpi < 0)
  	    tmpi = 0;
  	  tmpj = j+index_square1[jindex[jj]][1];
  	  if (tmpj > lny-1)
  	    tmpj = lny -1;
  	  if (tmpj <0)
    	  tmpj = 0;
    	for( int k = 0; k<lgrd-1;k++){
    	  if(localize_node_grid[tmpi][tmpj][kindex[k]] != -1){
    	    int itmp = localize_node_grid[tmpi][tmpj][kindex[k]];
    	    int icell = itmp;
    	    //test if the node is on the same cell
    	    if(icell != cell_index){//can try to bind  
    	      // test to see if node has been bound this time step
    	      for( int iin=0; iin<ncad; iin++){
      	      if( cell[icell].get_cadherin_attach(iin) == -1){
            		// cadherin should be detached
            		cur_int.attach =itmp*ncad+iin;
            		cell[icell].set_cadherin_attach(iin,ncad*cell_index+cadherin_index);
            		int icelltype; 
            		int dummy;
            		icelltype=cell[icell].get_type(dummy); //0 is spore cell, 1 is stalk 
            		
                int cellindextype;
            		cellindextype=cell[cell_index].get_type(dummy);

                  if (icelltype == cellindextype){
                    if (icelltype ==0){ // green to green attachment
                      cur_int.time = sto.Poisson(average_spore);// find new attach time average_time  is in seconds, cell types are different
                      cur_int.time = cur_int.time/3600.;// convert to hours
                    }else{ //blue to blue attach
                      cur_int.time = sto.Poisson(average_stalk[restart]);// find new attach time average_time  is in seconds, cell types are different
                      cur_int.time = cur_int.time/3600.;// convert to hours
                    }
                  }else{ // green to blue attachment
                		cur_int.time = sto.Poisson(average_different);// find new attach time average_time  is in seconds, cell types are different
                		cur_int.time = cur_int.time/3600.;// convert to hours
                  }
            			  
              		cell[icell].set_cadherin_time(iin,cur_int.time);
              		double vector[2];
              		for (int i=0; i<2; i++)
              		  vector[i]=cur_int.location[i];			  
              		cell[icell].set_cadherin_location(iin,vector,2);
              		return;
                
          	  }
	          }//end loop through cadherins 
	        }
	      } 
	    }// end k loop
    }// end jj loop  
  }// end if to bind
  else{// do not bind
		int icelltype; 
		int dummy;
		icelltype=cell[cell_index].get_type(dummy); //0 is spore cell, 1 is stalk 
		if (icelltype == 0) {//if cell is green, spore
			cur_int.time = sto.Poisson(average_unbound_cad_spore);// find new attach time average_time  is in seconds 
		}else {
			cur_int.time = sto.Poisson(average_unbound_cad_stalk);// find new attach time average_time  is in seconds 
		}
      cur_int.time =cur_int.time/3600;  //convert to hours
      cur_int.attach =-1;
    }
}


void Cell::update_cadherin_attachment_new(int cadherin_index, double dt,double t, double average_spore,double average_stalk[],double average_different, int restart, double average_unbound_cad_spore, double average_unbound_cad_stalk, int cell_index, Cell *cell){
  Cadherin& cur_int = cadherin[cadherin_index];
	
	int iicelltype; 
	int idummy;
	iicelltype=cell[cell_index].get_type(idummy); //0 is spore cell, 1 is stalk
	
  // generate a random number from 0 to 1 from a uniform distribution
  double temprand;
  temprand=rg.Random();  // random number between 0 and 1  
		
	double degreeG=35;//this changes the wedge where the cadherins can go
	double degreeB=35;//this changes the wedge where the cadherins can go
	temprand=(2*temprand-1)*degreeG*pi/180;
  
	double new_dir[2];

	//add temprand to the direction the cell is moving
  double new_length;
  new_length = rg.Random()*basic_length_cad +5 +7; // make length random adding onto basic length of 5 microns
	double change_dir = 1;

  double value[2];//I think this is a vector with x and y coordiantes of the front of the cell
  get_front(value, 2);
  double frontdeg;
	
	
  if( value[0] != 0){
		frontdeg = atan(value[1]/value[0]);//right no this is 0
  }else{
		frontdeg= pi/2.;
		if(value[1]<0){
			frontdeg = 3*pi/2;
		}
  }
	
  temprand = temprand +frontdeg;

  //THIS IS WHERE WE DETERMINE HOW MANY OF OUR BLUE CADHERINS ATTACH Back vs Front
  if(iicelltype==1){//if cell is blue
    if(new_dir[0]<0){//if the cadherin is in the back
      if (rg.Random()<cad_prob_change_back_to_front){// change direction of pseudopod
        temprand = temprand - pi;
      }
    }
  }else {//if cell is green
    if(new_dir[0]<0){//if the cadherin is in the back
      if (rg.Random()<prob_change_back_to_front){// change direction of pseudopod
        temprand = temprand-pi;
      }
    }
  }

  new_dir[0] = cos(temprand)*new_length;
  new_dir[1] = sin(temprand)*new_length;
  
	/*
	//THIS IS WHERE WE DETERMINE HOW MANY OF OUR BLUE CADHERINS ATTACH Back vs Front
	if(iicelltype==1){//if cell is blue
		if(new_dir[0]<0){//if the cadherin is in the back
			if (rg.Random()<cad_prob_change_back_to_front){// change direction of pseudopod
				new_dir[0] = new_dir[0]*-1;//
			}
		}
	}else {//if cell is green
		if(new_dir[0]<0){//if the cadherin is in the back
			if (rg.Random()<prob_change_back_to_front){// change direction of pseudopod
				new_dir[0] = new_dir[0]*-1;//
			}
		}
	}
*/

  double newx =  center[0] + new_dir[0];
  double newy  = center[1] +new_dir[1];
	calculate_cadherin_vector(cadherin_index);
	calculate_cadherin_front_indicator(cadherin_index);

  cur_int.location[0] = newx;
  cur_int.location[1] = newy;
  {// binds with another cadherin
    int i = static_cast<int> (newx/(lng_scale));
    int j = static_cast<int> (newy/(lng_scale));
      
    // Check the first grid then randomly loop through the 8 grid
    // squares until a node is located
    // Randomly loop through the values of k
    int kindex[lgrd];
    int jindex[9];
    // initialize  kindex
    for (int ii = 0; ii<lgrd; ii++)
	    kindex[ii]= ii;
    // randomize  kindex	
    for( int ii = lgrd-1; ii>=1; ii--){
	    int jj = rg.IRandom(0,ii);
	    int tmp = kindex[ii];
	    kindex[ii] = kindex[jj];
	    kindex[jj] = tmp;
    } 
    // initialize  jindex
    for (int ii = 0; ii<9; ii++)
	    jindex[ii]= ii;
    // randomize  jindex	
    for( int ii = 8; ii>=2; ii--){
	    int jj = rg.IRandom(1,ii) ;
	    int tmp = jindex[ii];
	    jindex[ii] = jindex[jj];
	    jindex[jj] = tmp;
    } 
      
    int tmpi;
    int tmpj;
    for( int jj = 0; jj<9; jj++){
	    tmpi = i+index_square1[jindex[jj]][0];
	    if (tmpi > lnx-1)
	      tmpi = lnx -1;
	    if (tmpi < 0)
	      tmpi = 0;
	    tmpj = j+index_square1[jindex[jj]][1];
	    if (tmpj > lny-1)
	      tmpj = lny -1;
	    if (tmpj <0)
	      tmpj = 0;
	    for( int k = 0; k<lgrd-1;k++){
	      if(localize_node_grid[tmpi][tmpj][kindex[k]] != -1){
	        // test to see if node has been bound this time step
	        int itmp = localize_node_grid[tmpi][tmpj][kindex[k]];
	        int icell = itmp;
	        for(int inode =0; inode<ncad; inode++){
	          //test if the node is on the same cell
	          if(icell != cell_index){//can try to bind
			        int icelltype; 
			        int dummy;
			        icelltype=cell[icell].get_type(dummy);       //0 is spore cell, 1 is stalk 
			        int cellindextype;
			        cellindextype=cell[cell_index].get_type(dummy);
	            if( cell[icell].get_cadherin_attach(inode) == -1 ||  cell[icell].get_cadherin_time(inode) < t){//start here
		            // cadherin should be detached
			  
			          cur_int.attach =itmp*ncad+inode;
			          int iattach = ncad*cell_index+cadherin_index;
			          cell[icell].set_cadherin_attach(inode,iattach);


                if (icelltype == cellindextype){
                  if (icelltype ==0){ // green to green attachment
                    cur_int.time = sto.Poisson(average_spore);// find new attach time average_time  is in seconds, cell types are different
                    cur_int.time = cur_int.time/3600.+t;// convert to hours and add t
                  }else{ //blue to blue attach
                    cur_int.time = sto.Poisson(average_stalk[restart]);// find new attach time average_time  is in seconds, cell types are different
                    cur_int.time = cur_int.time/3600.+t;// convert to hours and add t
                  }
                }else{ // green to blue attachment
                  cur_int.time = sto.Poisson(average_different);// find new attach time average_time  is in seconds, cell types are different
                  cur_int.time = cur_int.time/3600.+t;// convert to hours and add t
                }


		            cell[icell].set_cadherin_time(inode,cur_int.time);

		            double vector[2];
		            // Set position of cadherin which was free to the location of the reaching cadherin 
		            for (int i=0; i<2; i++)
		              vector[i]=cur_int.location[i];
		            cell[icell].set_cadherin_location(inode,vector,2);
		            return;
	            }
	          }
	        }//end inode loop
	      }
	    }// end k loop
    }// end jj loop
  }// this ends scoping of bind with another cadherin
  //Nothing to bind to so remain unbound
	int icelltype; 
	int dummy;
	icelltype=cell[cell_index].get_type(dummy); //0 is spore cell, 1 is stalk 
	if (icelltype == 0) {//if cell is green, spore
		cur_int.time = sto.Poisson(average_unbound_cad_spore);// find new attach time average_time  is in seconds 
	}else {
		cur_int.time = sto.Poisson(average_unbound_cad_stalk);// find new attach time average_time  is in seconds 
	}
  cur_int.time =cur_int.time/3600;  //convert to hours
  cur_int.attach =-1;
   
}


// Current code is written so all cadherin are attached immediately after detached
void Cell::update_cadherin_attachment(int cadherin_index, double dt){
  Cadherin& cur_int = cadherin[cadherin_index];
  // generate a random number from 0 to 1 from a uniform distribution
  double temprand;
  if(cur_int.front_indicator == 1){// cadherin is in the back 
    //    if(cur_int.time<t){// detach the cadherin and reattach it
        temprand=rg.Random();  // random number between 0 and 1
	double degree = 10;
	temprand = 2*(temprand -.5)*degree*pi/180; // gives an angle from -degree to degree in radians
	double new_dir[2];
	double new_length = 30;
	new_dir[0] = cos(temprand)*new_length;
	new_dir[1] = sin(temprand)*new_length;
	cur_int.location[0] = center[0] + new_dir[0];
	cur_int.location[1] = center[1] + new_dir[1];

	cur_int.vector[0] = new_dir[0];
	cur_int.vector[1] = new_dir[1];
	cur_int.front_indicator = 0;  // set the indicator to make the cadherin at the front of the cell
  }
}
void Cell::update_cadherins(double dt,double t, double average_spore,double average_stalk[],double average_different,int restart, double average_unbound_cad_spore, double average_unbound_cad_stalk, int ni, int index_cell,Cell *cell){
    calculate_cadherin_vector(ni);
    calculate_cadherin_front_indicator(ni);
    update_cadherin_attachment_new(ni,dt,t,average_spore,average_stalk, average_different,restart, average_unbound_cad_spore, average_unbound_cad_stalk,index_cell,cell);
    calculate_cadherin_vector(ni);
    calculate_cadherin_front_indicator(ni);
}

void Cell::update_cadherins_initial(double dt, double average_spore,double average_stalk[],double average_different,int restart, double average_unbound_cad_spore, double average_unbound_cad_stalk, int ni, int index_cell,Cell *cell){
    calculate_cadherin_vector(ni);
    calculate_cadherin_front_indicator(ni);
    update_cadherin_attachment_initial(ni,dt,average_spore,average_stalk,average_different,restart, average_unbound_cad_spore, average_unbound_cad_stalk,index_cell,cell);
    calculate_cadherin_vector(ni);
    calculate_cadherin_front_indicator(ni);
}

//update integrins
void Cell::update_integrin_attachment_initial(int integrin_index, double dt, double average_substrate_spore, double average_substrate_stalk[],int restart, double average_unbound_int_spore,double average_unbound_int_stalk, int cell_index, Cell *cell){
  Integrin& cur_int = integrin[integrin_index];
  double newx =   cur_int.location[0];
  double newy =   cur_int.location[1];
  int cell_indextype;  
  int dummy;
  cell_indextype=cell[cell_index].get_type(dummy);
  //Set a certain portion to stay unattached.  
  if(rg.Random()<prob_attach_initial_int){// decide integrin binds
    // If the cell is too far from substrate this routine is not called.
    // find distance to substrate.  The substrate is at y=0
    cur_int.attach = 0;// integrin is attached to substrate 
    cur_int.location[1]=0.0;
    if (cell_indextype==0){	
      cur_int.time = sto.Poisson(average_substrate_spore);// find new attach time average_time  is in seconds
    }else{
      cur_int.time = sto.Poisson(average_substrate_stalk[restart]);// find new attach time average_time  is in seconds
    }
    cur_int.time = cur_int.time/3600.;// convert to hours
    calculate_integrin_vector(integrin_index);
    calculate_integrin_front_indicator(integrin_index);
  }	  // end if for binding to substrate
  else{// do not bind
    if(cell_indextype==1){ //cell is stalk cell
      cur_int.time = sto.Poisson(average_unbound_int_stalk);// find new attach time average_time  is in seconds
      cur_int.time =cur_int.time/3600;  //convert to hours
    }else{
      cur_int.time = sto.Poisson(average_unbound_int_spore);// find new attach time average_time  is in seconds for spore cell
      cur_int.time =cur_int.time/3600;  //convert to hours
    }
    cur_int.attach =-1;
  }
}

void Cell::update_integrin_attachment_new(int integrin_index, double dt,double t, double average_substrate_spore,double average_substrate_stalk[],int restart, int cell_index, Cell *cell){
  Integrin& cur_int = integrin[integrin_index];
  int cell_indextype;
  int dummy;
  cell_indextype=cell[cell_index].get_type(dummy);
  double temprand;
  temprand=rg.Random();  // random number between 0 and 1
  double degree = 35; //this changes the wedge where the integrins can go
    //    temprand = 2*(temprand -.5)*degree*pi/180; // gives an angle from -degree to degree in radians
	temprand = (temprand)*degree*pi/180; // gives an angle from 0 to degree in radians
	//temprand = -pi/2.;
  double new_dir[2];
  double prob=0;

	//add temprand to the direction the cell is moving
  double new_length;
  new_length = rg.Random()*basic_length_int; // make length random adding onto basic length of 5 microns
  // double change_dir = 1;
  // prob = prob_change_from_front_to_back;
  //      change_dir = -1;
  // 	 //prob_change =.8;
  //   prob = prob_change_back_to_front;
  //   //      prob=2/dt;
  //   new_length = new_length*basic_length_back_factor;
  // //prob=.5; //comment out 
  // if (rg.Random()<prob){// change direction of pseudopod
  //      change_dir = -1*change_dir;
  // }
/*
  double value[2];
  double frontdeg;
  get_front(value, 2);
*/
	
	new_dir[0] = cos(temprand)*new_length;
  new_dir[1] = sin(temprand)*new_length;
  double newx =  center[0] + new_dir[0];
  double newy  = center[1] +new_dir[1];

  cur_int.location[0] = newx;
  cur_int.location[1] = newy;
  // not called if cell is not close enough to substrate. The substrate is at y=0
	cur_int.attach = 0;// integrin is attached to substrate
	cur_int.location[1]=0.; 
	if (cell_indextype==0){	
	  cur_int.time = sto.Poisson(average_substrate_spore);// find new attach time average_time  is in seconds
	}else{
	  cur_int.time = sto.Poisson(average_substrate_stalk[restart]);// find new attach time average_time  is in seconds
	}
	cur_int.time = cur_int.time/3600.+t;// convert to hours and add t.
}


void Cell::update_integrins(double dt,double t, double average_substrate_spore, double average_substrate_stalk[], int restart, int ni, int index_cell,Cell *cell){
    update_integrin_attachment_new(ni,dt,t,average_substrate_spore, average_substrate_stalk, restart, index_cell,cell);
    calculate_integrin_vector(ni);
    calculate_integrin_front_indicator(ni);
}

void Cell::update_integrins_initial(double dt, double average_substrate_spore, double average_substrate_stalk[],int restart,double average_unbound_int_spore,double average_unbound_int_stalk, int ni, int index_cell,Cell *cell){
  update_integrin_attachment_initial(ni,dt,average_substrate_spore,average_substrate_stalk,restart,average_unbound_int_spore,average_unbound_int_stalk,index_cell,cell);
}
//not used

void Cell::move_center(double dt){
  // This moves one cell on a fixed substrate with springs attaching the cadherins.  
  // The rest length of the springs is assumed to be zero.
  double tmp[2];
  double sum_k = 0;
  double drag_coefficient_cell = 10;
  tmp[0] = 0;
  tmp[1] = 0;
  for (int i=0; i<ncad; i++){
    if(cadherin[i].attach == 1){
      sum_k = sum_k + cadherin[i].kspring[0];
      tmp[0] = tmp[0] + cadherin[i].vector[0]*cadherin[i].kspring[0];
      tmp[1] = tmp[1] + cadherin[i].vector[1]*cadherin[i].kspring[0];
    }
  }
  for (int i=0; i<2; i++){
    tmp[i] = tmp[i]/sum_k;
    //This line is for the cell moving in the fluid with the drag_coefficient_cell
    //   center[i] = center[i] + (center[i] - tmp[i])*exp(-sum_k/drag_coefficient_cell*dt);
    //This line assumes cells moves with respect to lattice for the drag.  It uses a 
    //  midpoint type approximation for the integral.
    //   center[i] = tmp[i] + (center[i] - tmp[i])*exp(-sum_k/drag_coefficient_cell*dt)+exp(-sum_k/drag_coefficient_cell*dt*.5)*(L[i] -Lold[i]);
        center[i] = tmp[i]*(1.-exp(-sum_k/drag_coefficient_cell*dt)) + center[i];//correct same as fluid line

  }
}

void Cell::print_output(ofstream *fout, int n, double dt){
  // Print the output suitable for plot_cell.m
  // Print the time
  *fout << n*dt << "\n";
  // Print the cell center
  *fout << center[0] << " " << center[1] << "\n";
  // Print the cell type
  *fout << celltype << "\n"; 
  // Print the cadherin information
  for (int i=0; i<ncad; i++){
	  calculate_cadherin_front_indicator(i);
      double vector[2];
      get_cadherin_location(i, vector,2) ;    
    *fout << vector[0] << " " << vector[1] << " " << cadherin[i].attach << " "
	<< cadherin[i].front_indicator << "\n";
    //    cout<<cadherin[i].attach<<" att "<<i<<" i "<<"\n";
  }
  // Print the integrin information
  for (int i=0; i<nint; i++){
      double vector[2];
      get_integrin_location(i, vector,2) ;    
    *fout << vector[0] << " " << vector[1] << " " << integrin[i].attach << " "
	<< integrin[i].front_indicator << "\n";
    //    cout<<integrin[i].attach<<" att "<<i<<" i "<<"\n";
  }
}

void Cell::print_output_force(ofstream *fout, int n, double dt){
  // Print the output suitable for matlab forces on substrate
  for (int i=0; i<nint; i++){
    Integrin& cur_int = integrin[i];
    // Check if the integrin is in contact with the substrate
    if(cur_int.attach==0){
     // Print the time
      //*fout << n*dt << "\n";
      //*fout << cur_int.location[0] << "\n";
      double tmplength = sqrt(cur_int.vector[0]*cur_int.vector[0]+cur_int.vector[1]*cur_int.vector[1]);
      double stretch = tmplength-5;
      //      if(stretch<0){
	    //	stretch=0;
      //      }
      if(tmplength<1e-8){
	      tmplength=1.;
      }
      *fout<< cur_int.vector[0]/tmplength*stretch*cur_int.kspring<<"  ";
    }
  }
}
