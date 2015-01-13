/*************************************************************************************
 This code was initially written by Dr. John Dallon, Brigham Young University 
 
 
 This code was further modified by Nathan Wohlgemuth, Brigham Young University
 *************************************************************************************/

/****************************************Notes****************************************
 1. Cell speeds of 30 microns per hour and shape of 100 microns stretched on a substrate 
 (from Journal of Cell Science 111, 2423-2432 1998.  "Epidermal growth factor alters 
 fibroblast migration speed and directional persistence reciprocally and  in a 
 matrix-dependent manner" by Magaret F. Ware, Alan Wells, and Douglas A Lauffenburger.)
 2. Rho kinase mediates serum-induced contraction in fibroblast fibers independent of 
 myosin LC20 phosphorylation (Am J Physiol Cell Physiol 284: C599-C606, 2003) 
 3. Hiromi Nobe1, Koji Nobe1, Fabeha Fazal2, Primal de Lanerolle2, and Richard J. Paul1
 measured forces of 5.3e-7 Newtons for fibroblasts .003 Dynes from Balaban et al.,2001 
 and 1 Dyne per 10 square micron patch Munevar et al., 2001 
 4. For Dd Uchida and Yumura 2004 report that contact sites have
 duration of 19+-8 sec, more contacts in front of cell as back
 retracts, 5 contact gave speeds from 10-1microns/min.  15-20 contact
 give speeds of 0-4.
 5. Units are microns and hours and grams. 
 *************************************************************************************/

// Changing the reach length of the cadherins from 15 to 5 makes a more slender slug
// If the integrin length is 5 the back of the slug cannot stay attached 10 is better.

#include <cmath>
#include <fstream> 
#include <sstream>
#include <iostream>

//#include "util.h" // I cannot find this file anywhere (Jan 10, 2015)
#include "cellclass.h"
#include "randomc.h"
#include "stocc.h"
#include "/usr/local/include/cvode/cvode.h"
#include "/usr/local/include/cvode/cvode_spgmr.h"
#include "/usr/local/include/nvector/nvector_serial.h" 

namespace patch { // Fixes a known g++ compiler bug (gcc 4.8.0 and higher have this bug fixed); Added Jan 10, 2015
	template <typename T> std::string to_string(const T& n) {
		std::ostringstream stm;
        stm << n;
        return stm.str();
    }
}

using std::cout;

/* function declarations */
void initialize_cells(int ncx,int ncy, Cell *cell, ofstream& foutc);
void localize_node_grid_update_all (Cell *cell);
void function_node(double u[], double savf[],int nn,Cell *cell);
int check_flag(void* flagvalue, char *funcname, int opt);
int move_nodes(double t, N_Vector cc, N_Vector fval, void *f_data);	// routine defining the nonlinear function
void update_cadherin_cell_center_location(void *cvode_mem, int ncadherins, int cadherins_update[],int nintegrins, int integrins_update[], double dt, Cell * cell);
int Psolve(realtype t, N_Vector cc, N_Vector fval, N_Vector r, N_Vector z, realtype gamma, realtype delta, int lr, void *f_data, N_Vector tmp);

// TRandomMersenne rg(1283);  
// TRandomMersenne rg1(49495);
// int seed = (unsigned int) time(NULL);
int seed = 3894; 
CRandomMersenne rg1(seed);
StochasticLib2 sto1(13948); 
// TRandomMersenne rg1(seed); 

int  localize_node_grid[lnx][lny][lgrd];
double basic_length_cad; // in micrometer
double basic_length_int; // in micrometer
double basic_length_back_factor;
int initialrestart = 1; // 1 makes A cells move, 0 makes A cells stop

// cell values
double average_spore = 25;// average attach time in sec from another cadherin
double average_stalk[2] = {25,25};// average attach time in sec from another cadherin
double average_different = 25;// average attach time in sec from another cadherin
double average_unbound_cad_spore = 5;// average time spore cadherin is unbound in sec*
double average_unbound_cad_stalk = 5;// average time stalk cadherin is unbound in sec*

double average_unbound_int_spore = 10;// average time integrin is unbound in sec
double average_unbound_int_stalk = 10;// average time integrin is unbound in sec
double average_substrate_spore = 20; // 30 average detach in sec time from substrate or average time bound
double average_substrate_stalk[2] = {20,20}; // 300000 average detach in sec time from substrate for stopped A cell and 20 for moving A cell or average time bound

double prob_attach_initial = 0.5; // probability that the initial state of the cadherin is bound
double prob_attach_initial_int= 0.5; // probability that the initial state of the integrin is bound
//double prob_factor_back; // multiplied to prob_attach_substrate if in back l_2 for back
double init[2]; 
double velocity[2]; 
double cadviscosityfactor; //multiplied to viscosity for cadherins

int cadherins_update[totalcad];
int number_cadherins_update;

int integrins_update[totalint];
int number_integrins_update; 

double cad_averand = 0;
double cad_avedir[2] = {0,0};
double cad_numang = 0; 
double int_averand = 0;
double int_avedir[2] = {0,0};
double int_numang = 0; 
double avejindex[9] = {0,0,0,0,0,0,0,0,0};

double cell_forces[2*cx*cy]; // the forces on each cell

int main() {
	double tstart = 0;
	int ntime = int(tend/dt);
	int flag;
	int jindex[ncells];
	int nnmax = 2*(ncells+totalcad);
	
	ifstream fin; 
	fin.open("cellmotion.dat"); // open file to read
	if(fin.good()) { // if file opening reported no errors
		fin >> basic_length_int;
		fin >> basic_length_cad;
		//fin >> prob_factor_back; was 5 not used anywhere in code
		fin>>seed;
	}
	else { // if file opening reported errors
		cout << "cellmotion.dat not found" << endl;
	}
    fin.close();
	
    //prob_factor_back = prob_factor_back*.1; 
	
	fin.clear(); // reset the fin objects error flags
    fin.open("params.dat"); // if file opening reported no errors
	if(fin.good()) {
		fin >> average_spore;
		fin >> average_unbound_cad_spore;
	}
	else { // if file opening reported errors
		cout << "params.dat not found" << endl;
	}
    fin.close();
	
	ofstream foutm;
	foutm.open("mass_centers.dat",ios::app);
	foutm << average_spore << " " << average_unbound_cad_spore << " ";
	
	cout << "Running " << average_spore << "/" << average_unbound_cad_spore << endl;
	
	ofstream foutr; // attach ratios
	foutr.open(("attach_ratios_" + patch::to_string(average_spore) + "_" + patch::to_string(average_unbound_cad_spore) + ".dat").c_str());
	
	ofstream foutf; // force output
	foutf.open(("force_output_" + patch::to_string(average_spore) + "_" + patch::to_string(average_unbound_cad_spore) + ".dat").c_str());
	
	ofstream foutc; // cell output
	foutc.open(("cell_output_" + patch::to_string(average_spore) + "_" + patch::to_string(average_unbound_cad_spore) + ".dat").c_str());
	
	ofstream foutcf; // cell forces
	foutcf.open(("cell_forces_" + patch::to_string(average_spore) + "_" + patch::to_string(average_unbound_cad_spore) + ".dat").c_str());
	
	
	ifstream pin; 
	pin.open("cadviscosityfactor.txt"); // open file to read
	if(pin.good()) { // if file opening reported no errors
		pin>>cadviscosityfactor; //cadherin viscosity factor
	}
	else { // if file opening reported errors
		cout << "cadviscosityfactor.txt not found" << endl;
	}
	pin.close();
	
	Cell * cell = new Cell[ncells]; // get memory for array of cells
	initialize_cells(cx,cy,cell,foutm); //was foutc
	
	
	//TRandomMersenne rg1(seed);
	/* SETUP for the nonlinear solver CVODE */ 
	void *cvode_mem;
	
	cvode_mem = NULL;//is this where it needs to be set to the largest size of data we can get from caherins; check this if getting segmentation fault
	N_Vector sc;
	sc = N_VNew_Serial(nnmax);
	N_VConst_Serial(1.,sc); 
	/* Call KINCreate / KINMalloc to initialize KINSOL :
	 nvSpec is the nvSpec pointer used in the serial version
	 A pointer to KINSOL problem memory is returned and stored in cvode_mem . */
	
	cvode_mem = CVodeCreate(CV_BDF,CV_NEWTON);  // for stiff problems
	// cvode_mem= CVodeCreate(CV_ADAMS,CV_FUNCTIONAL); // for non stiff problems
	if ( check_flag (( void *) cvode_mem, " CVodeCreate ", 0)) return (1);
	flag = CVodeInit(cvode_mem, move_nodes,tstart,sc);   
	if ( check_flag (( void *) cvode_mem , " CVodeInit ", 0)) return (1);
	int maxstep = 20000;
	flag = CVodeSetMaxNumSteps(cvode_mem, maxstep);
	if ( check_flag ( &flag , "CVodeSetMaxNumSteps", 1)) return(1) ;
	
	localize_node_grid_update_all(cell);
	for(int j = 0; j < ncells; j++) {
		int restart = cell[j].restart;
		for(int i = 0; i<ncad; i++) {
			if(cell[j].get_cadherin_attach(i) == -1) {
				cell[j].update_cadherins_initial(dt,average_spore,average_stalk,average_different,restart, average_unbound_cad_spore,average_unbound_cad_stalk,i,j,cell);
			}
			int lim = 2;		  
			double vector[lim]; // for cadherin vector
			cell[j].get_cadherin_vector(i, vector, lim);
		}
	}
	// initialize the integrins by deciding if they should bind to substrate
	for(int j = 0; j < ncells; j++) {
		// check to see if the cell is close enough to the substrate
		double center[2];			       
		cell[j].get_center(center,2);
		int restart = cell[j].restart;
		if(center[1] < 7) { // 5 is the unstretched integrin length the substrate is at y=0 so if the cell center is 7+5 it is close enough to attach; 
			for(int i = 0; i < nint; i++) {
				cell[j].update_integrins_initial(dt, average_substrate_spore, average_substrate_stalk, restart, average_unbound_int_spore, average_unbound_int_stalk, i, j, cell);
			}
		}
	}
	
	foutc << ncells << "\n"; // Output the number of cells
	foutc << ncad << "\n"; // Output the number of cadherins
	foutc << nint << "\n"; // Output the number of integrins
	// output at time t = 0
	for(int j = 0; j < ncells; j++) {
		cell[j].print_output(foutc,0,dt);
	}
	cout << "Time = 0" << endl;
	
	// move the cell
	double prob;
	for(int it = 0; it <= ntime; it++) {
		int cadattach = 0;
		int intattach = 0;
		double t = (it+1)*dt;
		number_cadherins_update = 0;
		number_integrins_update = 0;
		
		for(int j = 0; j<ncells; j++) {    
			// initialize array
			for(int i = 0; i<ncad; i++) {  
				cadherins_update[j*ncad+i] = 0;
			}
			for(int i = 0; i<nint; i++) {  
				integrins_update[j*nint+i] = 0;
			} 
		}
		
		//loop through cells and determine if cadherins and integrins should release or not		
		for(int ii = 0; ii < ncells; ii++) { // start cadherin
			int iicelltype = cell[ii].get_type(); 
			for(int jj = 0; jj < ncad; jj++) {      
				double detach = cell[ii].get_cadherin_time(jj);
				int tempattach = cell[ii].get_cadherin_attach(jj);
				
				if(tempattach != -1) {			
					double length = cell[ii].get_cadherin_length(jj);  
					if(maxlength < length) {
						detach = 0;		  
					}
					cadattach++; //count how many caherins are attached
				}
				
				if(tempattach != -1) {
					int ajj = tempattach%ncad;			  
					int aii = static_cast<int>((tempattach - ajj)/ncad); //ncad (number of cadherins)=10
					
					int atempattach = cell[aii].get_cadherin_attach(ajj);	
					
					
					if(detach < t) {//time to change, detaches jj cadherin		  
						cell[ii].set_cadherin_attach(jj,-1);
						double tmp;
						if (iicelltype == 0) {//if cell is green, spore
							tmp = sto1.Poisson(average_unbound_cad_spore);//determine detach time//this is just to make a change
							tmp = tmp/3600 + t;
						}
						else {
							tmp = sto1.Poisson(average_unbound_cad_stalk);//determine detach time
							tmp = tmp/3600 + t;
						}
						
						cell[ii].set_cadherin_time(jj,tmp);//set detach time		
						
						if(atempattach != -1) {//detaches other cell's cadherin, ajj  
							cell[aii].set_cadherin_attach(ajj,-1);
							double atmp;
							if (iicelltype == 0) {//if cell is green, spore
								atmp = sto1.Poisson(average_unbound_cad_spore);//determine detach time
								atmp = atmp/3600 + t;
							}
							else {
								atmp = sto1.Poisson(average_unbound_cad_stalk);//determine detach time
								atmp = atmp/3600 + t;
							}
							cell[aii].set_cadherin_time(ajj,atmp);//set detach time
						}	  
					}
				}
			} //end cadherin
			
			for(int jj = 0; jj < nint; jj++) {
				double detach = cell[ii].get_integrin_time(jj);
				double length = cell[ii].get_integrin_length(jj);
				int tempattach = cell[ii].get_integrin_attach(jj);
				int num = cell[ii].get_type();
				int celltype = cell[ii].get_type();
				if(tempattach == 0) { //count how many integrins are attached
					intattach++;
				}
				if(detach < t) { //time to change				 
					if(cell[ii].get_integrin_attach(jj) == 0) { //detach   
						double tmp;
						cell[ii].set_integrin_attach(jj,-1);
						int iitype = cell[ii].get_type();				
						if(iitype == 1) { //if cell is celltype A, set tmp
							tmp = sto1.Poisson(average_unbound_int_stalk);
							tmp = tmp/3600 +t;
						}
						else { //if cell is not celltype A, set tmp
							tmp = sto1.Poisson(average_unbound_int_spore);
							tmp = tmp/3600 +t;
						}
						cell[ii].set_integrin_time(jj,tmp);   
					} // All cells reset integrin times
				} // end if
			} //end integrin  
		} // end cell loop
		
		// loop through the cells and attach the cadherins and integrins 		  
		for(int ii = 0; ii < ncells; ii++) jindex[ii] = ii;
		// randomize  jindex	
		for(int ii = ncells-1; ii >= 1; ii--) {
			int jj = rg1.IRandom(0,ii) ; 
			int tmp = jindex[ii];
			jindex[ii] = jindex[jj];
			jindex[jj] = tmp;
		}
		
		for(int j = 0; j < ncells; j++) {
			int restart = cell[jindex[j]].restart;
			
			for(int i = 0; i<ncad; i++) {
				if(cell[jindex[j]].get_cadherin_time(i) < t) {
					if(restart==0) { // if stopped it must be blue, so we won't check for blue celltype
					}
					else {
						cell[jindex[j]].update_cadherins(dt,t,average_spore,average_stalk,average_different,restart, average_unbound_cad_spore,average_unbound_cad_stalk,i,jindex[j],cell);//updating all cells in random order
					}	    
				}  		  
			}
			
			for(int i = 0; i < nint; i++) {      
				if(cell[jindex[j]].get_integrin_time(i)<t) {//if it's time for the integrin to attach
					//check to see if the cell is close enough to the substrate
					double center[2];			       
					cell[jindex[j]].get_center(center,2);
					if(center[1] < 7) { // 5 is the unstretched integrin length the substrate is at y=0 so if the cell center is 7+5 it is close enough to attach;	  
						int restart = cell[jindex[j]].restart;  
						if(cell[jindex[j]].get_type() == 0) { //if cell is green
							cell[jindex[j]].update_integrins(dt,t,average_substrate_spore, average_substrate_stalk, restart,i,jindex[j],cell);//update integrins
						}
						else { //if cell is blue
							if(restart == 1) { //and going
								cell[jindex[j]].update_integrins(dt,t,average_substrate_spore, average_substrate_stalk, restart,i,jindex[j],cell);//update integrins
							}
							else { //if cell is stopped
								//do nothing; don't attach
							}
						}
					}
				}
			}
		}
		
		//count the cadherins and integrins which need to be moved
		for(int ii = 0; ii < ncells; ii++) {
			for(int jj = 0; jj < ncad; jj++) { 
				if(cell[ii].get_cadherin_attach(jj) != -1) {
					cadherins_update[number_cadherins_update] = ii*ncad+jj;
					number_cadherins_update++;
				}
			}
			for(int jj = 0; jj < nint; jj++) {
				if(cell[ii].get_integrin_attach(jj)!=-1) {
					integrins_update[number_integrins_update] = ii*nint+jj;
					number_integrins_update++;
				}
			}    
		}
		
		if((number_cadherins_update+number_integrins_update) != 0) {
			update_cadherin_cell_center_location(cvode_mem, number_cadherins_update,cadherins_update, number_integrins_update,integrins_update,dt,cell);// THIS MOVES THE CENTERS AND THE NODES
		}
		localize_node_grid_update_all(cell);
		
		// output data 
		if(it%20 == 0) {  
			foutf << dt*it << "\n";
			for(int j = 0; j<ncells; j++) {
				cell[j].print_output(foutc,it+1,dt);
				cell[j].print_output_force(foutf,it,dt);
			}
			foutf << "\n";
		}
		else if(it%20 == 1) {
			for(int j=0; j<ncells; j++) {
				foutcf << cell_forces[2*j] << " " << cell_forces[2*j+1] << "\n";
			}
		}
		foutr << it*dt << " " << cadattach << " " << intattach << "\n";
	}// end to time loop
	
	cout << cad_averand/cad_numang << "\n";
	cout << cad_avedir[0]/cad_numang << "  " << cad_avedir[1]/cad_numang << "\n";
	for(int ii = 0; ii < 9; ii++) {
		cout << avejindex[ii]/cad_numang << " ";
	}
	cout << "\n";
	cout << int_averand/int_numang << "\n";
	cout << int_avedir[0]/int_numang << "  " << int_avedir[1]/int_numang << "\n";
	
	int totx = 0;
	for(int i = 0; i < ncells; i++) {
		double value[1];
		cell[i].get_center(value,1);
		totx += value[0];
	}
	foutm << 1.0*totx/ncells << "\n";
	
	foutc.close();
	foutf.close();
	foutr.close();
	foutcf.close();
	foutm.close();
	
	CVodeFree(&cvode_mem);
	return 0;
}

// end of main

/*************************************************************************************
 auxillary function definitions
 *************************************************************************************/



/************************************** move_the_cell ********************************
 *************************************************************************************/

void update_cadherin_cell_center_location(void *cvode_mem, int ncadherins, int cadherins_update[], int nintegrins, int integrins_update[], double dt,Cell *cell) {
	N_Vector cc;
	int flag;
	int dim = 2;
	int nn = dim*ncadherins;
	int nnint = dim*nintegrins;
	double tstart = 0.0;
	double tend;
	double rtol = 1e-6;
	double atol = 1e-6;
	int nnnmax = 2*(ncells+ncadherins);
	double tmpx[ncells*(ncad+1)*dim]; // holds the location of the cell center and cadherins
	double su[ncells*(ncad+1)*dim]; // scaling arra
	
	for(int m = 0; m <= ncells*(ncad+1)*dim-1; m++) {
		su[m] = 1.0;
		tmpx[m] = 0.0;
	}
	
	//    load in the cell center locations and ncadherin
	for( int i = 0; i < ncells; i++) {
		int ii = dim*i;
		double center[dim];
		cell[i].set_ncadherin(nn);
		cell[i].set_nintegrin(nnint);
		cell[i].get_center(center,dim);
		tmpx[ii] = center[0];
		tmpx[ii+1] = center[1];
	}
	//load in the cadherins locations
	for(int i = 0; i < ncadherins; i++) {
		int ii = dim*(i+ncells);
		double vector[dim];
		int itmp = cadherins_update[i];
		int inode = itmp%ncad;
		int icell = static_cast<int>((itmp - inode)/ncad);
		cell[icell].get_cadherin_location(inode,vector,dim);
		tmpx[ii] = vector[0]; // the x coordinate
		tmpx[ii+1] = vector[1];	//the y coordinate
	}
    cc = N_VMake_Serial((ncells+totalcad)*dim,tmpx);
	
	flag = CVodeReInit(cvode_mem, tstart, cc);   
	if(check_flag(&flag, "CVodeReInit", 1)) {
		return;
	}
	
	flag = CVodeSetUserData(cvode_mem, cell);
	if(check_flag(&flag, "CVodeSetUserData", 1)) {
		return;
	}
	
	flag = CVodeSStolerances(cvode_mem, rtol, atol);
	if(check_flag(&flag, "CVodeSStolerances", 1)) {
		return;
	}
	
	const int mmax = 15;      //dimension of krylov space
	
	flag = CVSpgmr(cvode_mem, PREC_NONE, mmax);
	if(check_flag(&flag, "CVSpgmr", 1)) {
		return;
	}
	double tret;
	tend = tstart + dt;
	
	flag = CVode(cvode_mem, tend,cc,&tret,CV_NORMAL);
	
	/* Call KINSol and print output concentration profile */
	if ( check_flag ( &flag , "CVode", 1)) {
		return;
	}
	
	for(int i = 0; i < ncells*dim + nn; i++) {
		tmpx[i] = NV_Ith_S(cc,i);
	}
	
	// Free memory
	N_VDestroy_Serial(cc);
	
	// Update the cell locations
	for(int i = 0; i < ncells; i++) {
		int ii = dim*(i);
		double center[2];
		center[0] = tmpx[ii];
		center[1] = tmpx[ii+1];
		cell[i].set_center(center,dim);
	}
	
    //update the cadherin locations
    for(int i = 0; i < ncadherins; i++) {
		int ii = dim*(i+ncells);
		double vector[dim];
		int itmp = cadherins_update[i];
		int inode = itmp%ncad;
		int icell = static_cast<int>((itmp - inode)/ncad);
		int attach = cell[icell].get_cadherin_attach(inode);
		vector[0] = tmpx[ii]; // the x coordinate
		vector[1] = tmpx[ii+1];	//the y coordinate
		cell[icell].set_cadherin_location(inode,vector,dim);
		inode = attach%ncad;
		icell = static_cast<int>((attach - inode)/ncad);
		cell[icell].set_cadherin_location(inode,vector,dim);
	}
	return;
} // end update_cadherin_cell_center_location()


/************************************************************************************
 Initialize cells
 ************************************************************************************/

void initialize_cells(int ncx,int ncy, Cell *cell, ofstream& fout) {
	double dx = 5; // 10 microns apart for each cell
	double dy = 5; // 10 microns apart for each cell
	
	int totx = 0;
	for(int k = 0; k < ncx; k++) {
		for(int j = 0; j < ncy; j++) {
			int ii = k + j*ncx;
			double xc = (k+1)*dx+1000;
			double yc = (j+1)*dy;
			totx += xc;
			double fx = 2.0*(rg1.Random() - 0.5);
			double fy = 2.0*(rg1.Random() - 0.5);
			fx = 1;
			fy = 0;
			double tmp = sqrt(fx*fx +fy*fy);
			fx /= tmp;
			fy /= tmp;
			int type = 0;
			type = 0;
			if(type == 0) {
				cell[ii].restart=initialrestart;
			}
			else {
				cell[ii].restart=0;
			}
			Cadherin cadherin_initial[ncad];
			for(int kk = 0; kk < ncad; kk++) {
				int i = kk;
				cadherin_initial[i].location[0]=xc+6*cos(2*pi*kk/ncad);
				cadherin_initial[i].location[1]=yc+6*sin(2*pi*kk/ncad);
				cadherin_initial[i].length = 5.0; // microns
				cadherin_initial[i].force = 0;
				cadherin_initial[i].mu = set_mu;
				cadherin_initial[i].kspring[0] = set_spring[type]*cad_spring_factor[type];
				cadherin_initial[i].kspring[1] = diff_spring;
				cadherin_initial[i].attach = -1;    
				cadherin_initial[i].time = 0;
			}
			// Initialize integrins
			Integrin integrin_initial[nint];
			for(int kk = 0; kk < nint; kk++) {
				int i = kk;
				integrin_initial[i].location[0]=xc;// set them at the cell center
				integrin_initial[i].location[1]=yc;
				integrin_initial[i].length = 5.; // microns
				integrin_initial[i].force = 0;
				integrin_initial[i].mu = set_mu;
				integrin_initial[i].kspring = set_spring[type];
				integrin_initial[i].attach = -1;    
				integrin_initial[i].time = 0;
			}
			cell[ii].initialize(xc,yc,cadherin_initial,integrin_initial,fx,fy,type);
		}    
	}
	fout << 1.0*totx/ncells << " ";
	return;
} // end initialize_cells

/*************************************************************************************
 localize_node_grid_update_all
 This function assigns each square grid all the cells in it. 
 *************************************************************************************/

void localize_node_grid_update_all(Cell *cell) {
	int l, m, i, j;
	
	/*  Set localize_node_grid to -1 */
	for(l = 0; l < lnx; l++) {
		for(m = 0; m < lny; m++) {
			for(int n = 0; n < lgrd; n++) {
				localize_node_grid[l][m][n] = -1;
			}
		}
	}
	for(int icell = 0; icell < ncells; icell++) {
		double vector[2];
		cell[icell].get_center(vector,2);
		/* determine which grid the cell is currently in */
		i = static_cast<int>((vector[0] + shift)/lng_scale);
		j = static_cast<int>((vector[1] + shift)/lng_scale);
		if(i > lnx-1) {
			i = lnx -1;
		}
		if(i < 0) {
			i = 0;
		}
		if(j > lny-1) {
			j = lny -1;
		}
		if(j < 0) {
			j = 0;
		}
		/* add the node to the grid it belongs to */
		int k = 0;
		while(localize_node_grid[i][j][k] != -1 && k<lgrd) {
			k++;
		}
		localize_node_grid[i][j][k] = icell;
		if(k == lgrd-1) {
			cout << k << " lngrid is too small " << endl;
		}
	} // end of icell loop
} // end localize_node_grid_update_all()


/*************************************************************************************
 Nonlinear ODEs to solve
 This function defines the RHS of the ODE system. 
 dv/dt = F_spring - drag = F_spring - v * mu
 By neglecting inertia, one get dx/dt = 1/mu * F_spring
 *************************************************************************************/
int move_nodes(double t, N_Vector cc, N_Vector fval, void *f_data) {
	int nn = NV_LENGTH_S(cc);
	
	double u[nn];
	double savf[nn];
	Cell *cell;
	
	for(int i = 0; i < nn; i++) { // initialize the arrays since we only use part of them
		u[i]=0;
		savf[i]=0;
	}
	
	cell = (Cell*)f_data;
	
	for(int i = 0; i < nn; i++) {
		u[i]=NV_Ith_S(cc,i);
	}
	
	function_node(u, savf, nn, cell);
	// Save the forces on each cell to global variable. Ideally, we'd only do the last time of the iterative solver.
	for(int i = 0; i < 2*cx*cy; i++) {
		cell_forces[i] = savf[i];
	} 
	for(int i = 0; i < nn; i++) {
		NV_Ith_S(fval, i) = savf[i];
    }
	return 0;
} // end move_nodes()

void function_node(double u[], double savf[], int nn, Cell *cell) {
	// zero out the forces
	// double vector[2];
	double vector2[2];
	int dim = 2;
	double tmpx;
	double tmpy;
	double dist;
	double mur;
	mur = 1.0/set_mu;
	for(int i = 0; i < nn; i++) {
		savf[i] = 0.0;
	}
	for(int i = 0; i < ncells; i++) {
		int ii = 2*i;
		// push away from the substrate
		if(u[ii+1] < 5) { // substrate is at zero and add 5 to push away.  For no substrate move to -500
			savf[ii+1] = savf[ii+1] + mur*10*cell[i].get_cadherin_kspring(1,0)*(exp(5-u[ii+1])-1); 
		}
		// add a body spring to keep the cell centers at least 5 microns away from each other
		for(int j = 0; j < ncells; j++) {
			if(j != i) {
				tmpx = u[ii] - u[2*j];
				tmpy = u[ii+1] - u[2*j+1];
				dist = sqrt(tmpx*tmpx+tmpy*tmpy);
				if(dist < 6) { // push cells away from each other
					savf[ii] = savf[ii] + mur*10*cell[i].get_cadherin_kspring(1,0)*(exp(6-dist)-1)*tmpx/dist;
					savf[ii+1] = savf[ii+1] +mur*10*cell[i].get_cadherin_kspring(1,0)*(exp(6-dist)-1)*tmpy/dist;
				}
			} // endif not same cell
		} // end for loop; j through jcells
	} // end for loop; i through ncells
	
	// Do the integrin forces
	int ntmpint = cell[0].get_nintegrin();
	for(int i = 0; i < ntmpint; i += dim) {
		int itmp = integrins_update[(i)/2]; 
		int inode = itmp%nint;
		int icell = static_cast<int>((itmp - inode)/nint); 
		double vector[2];
		cell[icell].get_integrin_location(inode,vector,2);
		// add spring constant terms
		double tmpx = vector[0]-u[dim*(icell)];
		double tmpy = vector[1]-u[dim*(icell)+1];
		double tmplength = tmpx*tmpx+tmpy*tmpy;
		tmplength = sqrt(tmplength);
		double stretch = tmplength -5; // assumes normal integrin distance is 5 microns
		if(tmplength > 1e-8) {
			tmplength = 1/tmplength;
		}
		else {
			tmplength = 0; 
		}
		int type = cell[icell].get_type();
		double tmptype = int_spring_factor[type];
		
		double tmpfx = cell[icell].get_integrin_kspring(inode)*stretch*tmpx*tmplength*tmptype; // hookes law
		double tmpfy = cell[icell].get_integrin_kspring(inode)*stretch*tmpy*tmplength*tmptype; // hookes law
		
		savf[2*(icell)] = savf[dim*(icell)] + mur*tmpfx;
		savf[2*(icell)+1] = savf[dim*(icell)+1] + mur*tmpfy;
	} //end integrin forces
	// Do the cadherin forces
	int ntmp = cell[0].get_ncadherin(); // gets cell's general cadherin # 
	for(int i = dim*ncells; i < ntmp + dim*ncells; i += dim) { // dim = 2
		int itmp = cadherins_update[(i-dim*ncells)/2]; // identifies cadherin to be updated ???
		int inode = itmp%ncad; // local cadherin #
		int icell = static_cast<int>((itmp - inode)/ncad); // cell #
		double tmpx = u[i] - u[dim*(icell)];
		double tmpy = u[i+1] - u[dim*(icell) + 1];
		double tmplength = tmpx*tmpx+tmpy*tmpy;
		tmplength = sqrt(tmplength);
		double stretch = tmplength - 5; // assumes normal cadherin distance is 5 microns
		if(tmplength > 1e-8) {
			tmplength = 1/tmplength;
		}
		else {
			tmplength = 0;
		}
		int attach = cell[icell].get_cadherin_attach(inode); // returns general cadherin # of other cell's cadherin it's attached to
		int ainode = attach%ncad; 
		int aicell = static_cast<int>((attach - ainode)/ncad);
		int tmptype = 0;
		int icelltype = cell[icell].get_type(); // 0 is spore cell, 1 is stalk 
		int aicelltype = cell[aicell].get_type();
		if(icelltype != aicelltype) {
			tmptype = 1;
		}
		
		double tmpfx = cell[icell].get_cadherin_kspring(inode, tmptype)*stretch*tmpx*tmplength; // hookes law
		double tmpfy = cell[icell].get_cadherin_kspring(inode, tmptype)*stretch*tmpy*tmplength; // hookes law
		
		savf[2*(icell)] = savf[dim*(icell)] + mur*tmpfx;
		savf[2*(icell)+1] = savf[dim*(icell)+1] + mur*tmpfy;
		
		//move cadherins
		savf[i] = savf[i] - mur*tmpfx/cadviscosityfactor;
		savf[i+1] = savf[i+1] - mur*tmpfy/cadviscosityfactor;
		
		// add spring constant terms
		tmpx = u[i] - u[dim*aicell];
		tmpy = u[i+1] - u[dim*aicell + 1];
		tmplength = tmpx*tmpx + tmpy*tmpy;
		tmplength = sqrt(tmplength);
		stretch = tmplength - 5; // assumes normal cadherin distance is 5 microns
		if(tmplength > 1e-8) {
			tmplength = 1/tmplength;
		}
		else {
			tmplength = 0;
		}
		tmpfx = cell[aicell].get_cadherin_kspring(ainode, tmptype)*stretch*tmpx*tmplength; // hookes law
		tmpfy = cell[aicell].get_cadherin_kspring(ainode, tmptype)*stretch*tmpy*tmplength; // hookes law
		
		savf[i] = savf[i] - mur*tmpfx/cadviscosityfactor;
		savf[i+1] = savf[i+1] - mur*tmpfy/cadviscosityfactor;
		// add to the cell center
		savf[2*(aicell)] = savf[dim*(aicell)] + mur*tmpfx;
		savf[2*(aicell)+1] = savf[dim*(aicell)+1]  + mur*tmpfy;
	} // end cadherin forces	
	return;
} // end function_node()


/*
 * Check function return value ...
 * opt == 0 means SUNDIALS function allocates memory so check if
 * returned NULL pointer
 * opt == 1 means SUNDIALS function returns a flag so check if
 * flag >= 0
 * opt == 2 means function allocates memory so check if returned
 * NULL pointer
 */
int check_flag(void* flagvalue, char* funcname, int opt) {
	int* errflag;
	
	/* Check if SUNDIALS function returned NULL pointer - no memory allocated */
	if(opt == 0 && flagvalue == NULL) {
		fprintf(stderr, "\n SUNDIALS_ERROR : %s()  failed  -  returned   NULL   pointer \n\n", funcname);
		return 1;
	}
	else if(opt == 1) { /* Check if flag < 0 */
		errflag = (int *) flagvalue;
		if (* errflag < 0) {
			fprintf (stderr, "\n SUNDIALS_ERROR : %s()  failed   with   flag  = %d \n\n", funcname, *errflag);
			return 1;
		}
	}
	else if(opt == 2 && flagvalue == NULL) { /* Check if function returned NULL pointer - no memory allocated */
		fprintf(stderr, "\n MEMORY_ERROR : %s() failed  -  returned   NULL pointer \n\n", funcname);
		return 1;
	}
	return 0;
} // end check_flag()
