/*************************************************************************************
      This code was initially written by Dr. John Dallon, Brigham Young University 
      
 
	  This code was further modified by Nathan Wohlgemuth, Brigham Young University
*************************************************************************************/

/*************************Notes*******************************************************
1. cell speeds of 30 microns per hour and shape of 100 microns stretched on a substrate 
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

5. Units are microns and hours 
and grams. 
*************************************************************************************/

//NOTES  Changing the reach length of the cadherins from 15 to 5 makes a more slender slug
// If the integrin length is 5 the back of the slug cannot stay attached 10 is better.


#include <cmath>
#include <fstream> 
#include <sstream>

//#include "util.h" // I cannot find this file anywhere (Jan 10, 2015)
#include "cellclass.h"
#include "randomc.h"
#include "stocc.h"
#include "/usr/local/include/cvode/cvode.h"
#include "/usr/local/include/cvode/cvode_spgmr.h"
#include "/usr/local/include/nvector/nvector_serial.h" 
 
namespace patch {
    template <typename T> std::string to_string(const T& n) {
        std::ostringstream stm;
        stm << n;
        return stm.str();
    }
}

/* function declarations */
void initialize_cells(int ncx,int ncy, Cell *cell, ofstream& foutc);
void localize_node_grid_update_all (Cell *cell);
void   function_node(double u[], double savf[],int nn,Cell *cell); 
int check_flag(void* flagvalue, char *funcname, int opt);
int  move_nodes (double t, N_Vector cc, N_Vector fval, void *f_data);	// routine defining the nonlinear function
void update_cadherin_cell_center_location(void *cvode_mem, int ncadherins, int cadherins_update[],int nintegrins, int integrins_update[], double dt, Cell * cell);
int Psolve(realtype t, N_Vector cc, N_Vector fval, N_Vector r, N_Vector z, realtype gamma, realtype delta, int lr, void *f_data, N_Vector tmp);

//double signal( double time);
//using namespace std;     
// assume domain is lnx*lng_scale 
 
   // TRandomMersenne rg(1283);  
// TRandomMersenne rg1(49495);
//int seed = (unsigned int) time(NULL);
int seed = 3894;//1283; 
CRandomMersenne rg1(seed);
StochasticLib2 sto1(13948); 
// TRandomMersenne rg1(seed); 
 
int  localize_node_grid[lnx][lny][lgrd];
double basic_length_cad; // in micrometer
double basic_length_int; // in micrometer
double basic_length_back_factor;
int initialrestart=1; //1 makes A cells move, 0 makes A cells stop

 
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

double prob_attach_initial=.5; // probability that the initial state of the cadherin is bound
double prob_attach_initial_int=.5; // probability that the initial state of the integrin is bound
//double prob_factor_back;// multiplied to prob_attach_substrate if in back l_2 for back
double init[2]; 
double velocity[2]; 
double cadviscosityfactor;//multiplied to viscosity for cadherins

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
	int ntime = int (tend/dt);
	int flag;
	int jindex[ncells];
	int nnmax = 2*(ncells+totalcad);
	
	ifstream fin; 
	fin.open("cellmotion.dat"); // open file to read
	if(fin.good()) {
		fin >> basic_length_int;
		fin >> basic_length_cad;
		//fin >> prob_factor_back; was 5 not used anywhere in code
		fin>>seed;
	}
	else {
		cout << "cellmotion.dat not found" << endl;
	}
    fin.close();

    //prob_factor_back = prob_factor_back*.1; 
  
    fin.open("params.dat");
	if(fin.good()) {
		fin >> average_spore;
		fin >> average_unbound_cad_spore;
	}
	else {
		cout << "params.dat not found" << endl;
	}
    fin.close();
  
	ofstream foutm;
	foutm.open("mass_centers.dat",ios::app);
	foutm << average_spore << " " << average_unbound_cad_spore << " ";

	cout << "Running " << average_spore << "/" << average_unbound_cad_spore << "\n";

	ofstream foutr;
	foutr.open(("attach_ratios_" + patch::to_string(average_spore) + "_" + patch::to_string(average_unbound_cad_spore) + ".dat").c_str());
		
	ofstream foutf;
	foutf.open(("force_output_" + patch::to_string(average_spore) + "_" + patch::to_string(average_unbound_cad_spore) + ".dat").c_str());
  
	ofstream foutc;
	foutc.open(("cell_output_" + patch::to_string(average_spore) + "_" + patch::to_string(average_unbound_cad_spore) + ".dat").c_str());

	ofstream foutcf;
	foutcf.open(("cell_forces_" + patch::to_string(average_spore) + "_" + patch::to_string(average_unbound_cad_spore) + ".dat").c_str());


	ifstream pin; 
	pin.open("cadviscosityfactor.txt"); // open file to read
	if(pin.good()) {
		pin>>cadviscosityfactor; //cadherin viscosity factor
	}
	else {
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
	for(int j = 0; j<ncells; j++){
		int restart = cell[j].restart;
		for(int i = 0; i<ncad; i++){
			if(cell[j].get_cadherin_attach(i)==-1) {
				cell[j].update_cadherins_initial(dt,average_spore,average_stalk,average_different,restart, average_unbound_cad_spore,average_unbound_cad_stalk,i,j,cell);
			}
			int lim=2;		  
			double vector[lim];//for cadherin vector
			cell[j].get_cadherin_vector(i, vector, lim);
	
			if (abs(vector[0])<1) {
			}
		}
	}
	//initialize the integrins by deciding if they should bind to substrate
	for(int j = 0; j<ncells; j++){
		//check to see if the cell is close enough to the substrate
		double center[2];			       
		cell[j].get_center(center,2);
		int restart = cell[j].restart;
		if(center[1] < 7){// 5 is the unstretched integrin length the substrate is at y=0 so if the cell center is 7+5 it is close enough to attach; 
			for(int i = 0; i<nint; i++){
				cell[j].update_integrins_initial(dt,average_substrate_spore,average_substrate_stalk, restart,average_unbound_int_spore, average_unbound_int_stalk,i,j,cell);
			}
		}
	}
  
	// Output some preliminary stuff (just the number of cells)
	foutc << ncells << "\n";
	foutc << ncad << "\n";
	foutc << nint << "\n";
	// output at time t = 0
	for(int j = 0; j<ncells; j++) {
		cell[j].print_output(foutc,0,dt);
	}
	cout << "Time = " << 0*dt << "\n";


	// move the cell
	double prob;
	for(int it = 0; it <= ntime; it++) {
		int cadattach = 0;
		int intattach =0;
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
			int dummy;
			int iicelltype = cell[ii].get_type(dummy); 
			for(int jj = 0; jj < ncad; jj++){      
				double detach = cell[ii].get_cadherin_time(jj);
				int tempattach = cell[ii].get_cadherin_attach(jj);
		  
				if(tempattach != -1) {			
					double length = cell[ii].get_cadherin_length(jj);  
					if(maxlength < length){
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

			for(int jj=0; jj < nint; jj++) {
				double detach = cell[ii].get_integrin_time(jj);
				double length = cell[ii].get_integrin_length(jj);
				int celltype;
				int dummy;
				int tempattach = cell[ii].get_integrin_attach(jj);
			  
				int num;
				num = cell[ii].get_type(dummy);

				celltype = cell[ii].get_type(dummy);

				if(tempattach == 0) { //count how many integrins are attached
					intattach++;
				}

				if(detach < t) { //time to change
				 
					if(cell[ii].get_integrin_attach(jj) == 0) { //detach   
						double tmp;
						cell[ii].set_integrin_attach(jj,-1);
						int iitype;
						cell[ii].get_type(iitype);				
						if(iitype == 1) { //if cell is celltype A, set tmp
							tmp = sto1.Poisson(average_unbound_int_stalk);
							tmp = tmp/3600 +t;
						}
						else { //if cell is not celltype A, set tmp
							tmp = sto1.Poisson(average_unbound_int_spore);
							tmp = tmp/3600 +t;
						}
						cell[ii].set_integrin_time(jj,tmp);   
					} //All cells reset integrin times
				}
			} //end integrin  
		} // end cell loop


		  
		// loop through the cells and attach the cadherins and integrins 
		  
		for(int ii = 0; ii<ncells; ii++) jindex[ii]= ii;
		// randomize  jindex	
		for(int ii = ncells-1; ii>=1; ii--) {
			int jj = rg1.IRandom(0,ii) ; 
			int tmp = jindex[ii];
			jindex[ii] = jindex[jj];
			jindex[jj] = tmp;
		}

		for(int j = 0; j < ncells; j++){
			int restart = cell[jindex[j]].restart;
			
			for(int i = 0; i<ncad; i++){
				if(cell[jindex[j]].get_cadherin_time(i) < t){
					if(restart==0) { // if stopped it must be blue, so we won't check for blue celltype
					}
					else {
						cell[jindex[j]].update_cadherins(dt,t,average_spore,average_stalk,average_different,restart, average_unbound_cad_spore,average_unbound_cad_stalk,i,jindex[j],cell);//updating all cells in random order
					}	    
				}  		  
			}
			
			for(int i = 0; i < nint; i++) {      
				if(cell[jindex[j]].get_integrin_time(i)<t){//if it's time for the integrin to attach
					//check to see if the cell is close enough to the substrate
					double center[2];			       
					cell[jindex[j]].get_center(center,2);
					if(center[1] < 7) { // 5 is the unstretched integrin length the substrate is at y=0 so if the cell center is 7+5 it is close enough to attach;	  
						int restart = cell[jindex[j]].restart;
						int dummy;  
						if(cell[jindex[j]].get_type(dummy) == 0) { //if cell is green
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
				if(cell[ii].get_integrin_attach(jj)!=-1){
					integrins_update[number_integrins_update] = ii*nint+jj;
					number_integrins_update++;
				}
			}    
		}

		if((number_cadherins_update+number_integrins_update) != 0){
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
			for (int j=0; j<ncells; j++) {
				foutcf << cell_forces[2*j] << " " << cell_forces[2*j+1] << "\n";
			}
		}
		//cout<<"Time = "<<(it+1)*dt<<"\n";


		foutr << it*dt << " " << cadattach << " " << intattach << "\n";
		
	}// end to time loop

	cout << cad_averand/cad_numang << "\n";
	cout << cad_avedir[0]/cad_numang << "  " << cad_avedir[1]/cad_numang << "\n";
	for(int ii = 0; ii <9; ii++){
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

void update_cadherin_cell_center_location(void *cvode_mem, int ncadherins, int cadherins_update[], int nintegrins, int integrins_update[], double dt,Cell *cell){
  N_Vector cc;
  int flag;
  int dim =2;
  int nn = dim*(ncadherins);
  int nnint = dim*(nintegrins);
  double tstart =0.;
  double tend;
  double rtol=1e-6;
  double atol=1e-6;
  int nnnmax = 2*(ncells+ncadherins);
  //cout <<" nn "<<nn<<" "<<ncadherins<<" "<<ncells<<"\n";
/* Allocate memory , and set problem data , initial values , tolerances */
  //  cc = N_VNew_Serial ( nn );
  //  sc = N_VNew_Serial ( nn );
  // Solve the nonlinear system using nksol F(x) = 0  F(x) is defined in Fspring. 
  // Define several working arrays and variables
  {
    //    const int mmax = 20;      //dimension of krylov space
    //    const int lrw = 8*nn + (nn+1)*mmax +2*mmax*mmax+3*mmax+40; // set the size of rwork 
    //    const int liw = 20 +mmax+20;	// set the size of iwork
    //    double ftol = dt ;	// set the stopping criteria 
    //    double stptol = 1.e-20;	// set the smallest step size
    //    int iopt = 1;        // iopt=0 no optional inputs
    //    int mf = 3;	         // mf=1 dogleg, mf=2 arnoldi, mf=3 gmres
    //    int mdif = 0;	// mdif=0 no user defined jac
    //    int ipflag = 0;	// ipflag=0 no preconditioning
    
 
    //    int iterm;			//  output flag

    //    double rwork[lrw];
    //    int iwork[liw];
    //    double tmpx[nn];           // holds the location of the cell center and cadherins
            double tmpx[ncells*(ncad+1)*dim];           // holds the location of the cell center and cadherins
    //    double tmpx[ncells*dim+nn+1];           // holds the location of the cell center and cadherins
	//    double tmpx[nnnmax];           // holds the location of the cell center and cadherins
    //    double tmpx[nnmax];           // holds the location of the cell center and cadherins
    //    double of[nn];             //the residual of the nonlinear function
    //    double su[nn];             // scaling array
            	  double su[ncells*(ncad+1)*dim];        // scaling arra
    //    double su[ncells*dim+nn+1];        // scaling arra
	//	  double su[nnnmax];             // scaling array
    //    for(int m=0; m<=nn-1; m++)
	//	  for(int m=0; m<=nnnmax-1; m++)

	  for(int m=0; m<=ncells*(ncad+1)*dim-1; m++)
	  //    for(int m=0; m<=ncells*dim+nn; m++)
      {
	su[m]=1.0;
	tmpx[m]=0.;
      }
    //  sc = N_VMake_Serial(nn,su);
    //      sc = N_VMake_Serial(ncells*22,su);
      //  sc = N_VMake_Serial(nnmax,su);
    // for (int i = 0; i < 8; i++)
    //      {
    //	iwork[i] = 0;
    //      }
    //    iwork[7]=1000;
    //    rwork[0]=100.*dt; // for fast moving cells may need to increase this
      // rwork[0]=1.; // for fast moving cells may need to increase this
      //     rwork[1]=dt;
    //    rwork[2]=0.;
    //    iwork[0]=mmax;
    //    iwork[6]=mmax;
    


    
    //    load in the cell center locations and ncadherin
    for( int i = 0; i<ncells; i++){
      int ii = dim*(i);
      double center[dim];
      cell[i].set_ncadherin(nn);
      cell[i].set_nintegrin(nnint);
      cell[i].get_center(center,dim);
      tmpx[ii] = center[0];
      tmpx[ii+1] = center[1]; 
      //    cout<< " center "<<center[0]<<" "<<center[1]<<" dim "<<dim<<" i "<<i<<"\n";
    }
    //load in the cadherins locations
    for( int i = 0; i<ncadherins; i++){
	int ii = dim*(i+ncells);
	double vector[dim];
	int itmp = cadherins_update[i];
	int inode = itmp%ncad;
	int icell = static_cast<int>((itmp - inode)/ncad);
	cell[icell].get_cadherin_location(inode,vector,dim);
	//	cout<<"icell "<<icell<<"inode "<<inode<<"vect "<<vector[0]<<" "<<vector[1]<<" "<<ncadherins<<"\n";
	tmpx[ii] = vector[0];       // the x coordinate
	tmpx[ii+1] = vector[1];	//the y coordinate
    } 
    //    cc = N_VMake_Serial(nn,tmpx);
    //  cc = N_VMake_Serial(ncells*(ncad+1)*dim,tmpx);
    cc = N_VMake_Serial((ncells+totalcad)*dim,tmpx);
       // cc = N_VMake_Serial(ncells*dim+nn+1,tmpx);
    //      cc = N_VMake_Serial(nnnmax,tmpx);


     flag = CVodeReInit(cvode_mem, tstart,cc);   
     if ( check_flag ( &flag , "CVodeReInit", 1)) return ;
     flag = CVodeSetUserData (cvode_mem , cell);
     if ( check_flag ( &flag , "CVodeSetUserData", 1)) return ;
     flag = CVodeSStolerances(cvode_mem,rtol,atol);
     if ( check_flag ( &flag , "CVodeSStolerances", 1)) return ;

     // double hmax=dt/4.;
     // if ( check_flag ( &flag , " KINMalloc ", 1)) return(1);
     const int mmax = 15;      //dimension of krylov space
     flag = CVSpgmr (cvode_mem ,PREC_NONE, mmax );
     //    flag = KINSpbcg (cvode_mem , mmax );
     // flag = KINSptfqmr(cvode_mem , mmax );
     if ( check_flag ( &flag , "CVSpgmr", 1)) return;
     double tret;
     tend=tstart+dt;
 
     //flag = CVode(cvode_mem, tend,cc,&tret,CV_ONE_STEP);

     flag = CVode(cvode_mem, tend,cc,&tret,CV_NORMAL);
  
     /* Call KINSol and print output concentration profile */
     if ( check_flag ( &flag , "CVode", 1)) return;



      /* fortran nonlinear solver  */
//    nksol_ (&nn,&nn,tmpx,of, &move_nodes_, &jac_,su,su, &ftol, &stptol,rwork,&lrw,iwork,&liw,&iopt, &iterm,&pset_,&psol_,&mf,&mdif,&ipflag);
//      {
//	printf ("+++++++++++++++++++++iterm %d of %f %f\n", iterm,of[0],of[1]);
//	printf ("+++++++++++++++++++++iterm %f of %f %f\n", tmpx[0],tmpx[1],tmpx[2]);
//      }
      /*  count how many errors from nksol */
//    if(iterm!=1)
//      cerror=cerror+1;
// for( int i = 0; i<ncells*(ncad+1)*dim; i++){
 for( int i = 0; i<ncells*dim+nn; i++){
   tmpx[i] = NV_Ith_S(cc,i);
   //   cout<<"tmmmm "<<tmpx[i]<<"\n";
 }
 //free memory
  N_VDestroy_Serial(cc);
  // update the cell locations
  for( int i = 0; i<ncells; i++){
      int ii = dim*(i);
      double center[2];
      center[0] = tmpx[ii];
      center[1] = tmpx[ii+1];
      // cout << "t "<<tmpx[ii]<<" "<<tmpx[ii+1]<<" c "<<center[0]<<" "<<center[1]<<"\n";
      //	       cout<<" "<<NV_Ith_S(cc,ii)<<"\n";
      cell[i].set_center(center,dim);
      //              cout<< " ***center "<<center[0]<<" "<<center[1]<<" dim "<<dim<<" i "<<i<<"\n";
    }
    //update the cadherin locations
    for( int i = 0; i<ncadherins; i++){
      //              cout<< " i*** "<<i<<"\n";
	int ii = dim*(i+ncells);
	double vector[dim];
	int itmp = cadherins_update[i];
	int inode = itmp%ncad;
	int icell = static_cast<int>((itmp - inode)/ncad);
	int attach = cell[icell].get_cadherin_attach(inode);
	vector[0] = tmpx[ii];       // the x coordinate
	vector[1] = tmpx[ii+1];	//the y coordinate
	cell[icell].set_cadherin_location(inode,vector,dim);
	inode = attach%ncad;
	icell = static_cast<int>((attach - inode)/ncad);
	  //	  cout<< " here ***"<<i<<" " <<icell<<" "<<inode<<" "<<attach<<"\n";
	cell[icell].set_cadherin_location(inode,vector,dim);	
	  //              cout<< " i "<<i<<"\n";
    } 
  // update the localize_node_grid
    //  localize_node_grid_update_all(cell);
    // no need to update the integrin locations - they do not move.
  }
}


/************************************************************************************
                            Initialize cells 
************************************************************************************/

void initialize_cells(int ncx,int ncy, Cell *cell, ofstream& fout){
  double dx = 5; // 10 microns apart for each cell
  double dy = 5; // 10 microns apart for each cell

  
  //  double set_spring = 4.4e7;// give a force of about .005 dynes at 15 microns
  //  double set_mu = 1.76e6;
  // double set_spring = 4.;// give a force of about .00530 dynes at 15 microns
  //  double set_mu = .176;

  //october
  int totx = 0;
  for (int k =0; k<ncx; k++){
    for (int j =0; j<ncy; j++){
      int ii = k+j*ncx;
      double xc = (k+1)*dx+1000;
      double yc = (j+1)*dy;
      totx += xc;
      double fx = 2.*(rg1.Random()-.5);
      double fy = 2*(rg1.Random()-.5);
      //      fx = 1000-xc;
      //      fy = 1000-yc;
      fx =1;
      fy=0;
      double tmp = sqrt(fx*fx +fy*fy);
      fx = fx/tmp;
      fy = fy/tmp;
      int type;
/*      
		int temptype=((k+j)%3);  //this is to make the cell types dispersed through the slug
      if(temptype<1){
	type=1;
      }else{
	type=0;
      }
 */     

      //this splits the cell into two types
      //1 is stalk, 0 is spore; 
     //cout << ncx<<"\n";
 //   double split = ncx / 2.;
		//cout << "split"<<split<<"\n";
//		if(k>split && j<3){
//	    type =1;
//    }else{
	    type =0;
//    }


      if (type == 0) {
		  cell[ii].restart=initialrestart;
		  //cout << "greenrestart = 1; first timestep?";
	  }else{
		  cell[ii].restart=0;
		  //cout << "bluerestart = 0; first timestep?";
	  }
      
      // direct cells
    
      //      int direction = 1;
      //           if(fy <0) 
      //	     direction = -1;
      //
	 //                 fx = direction;
	 //               fy = 0;
      //      if(fx<0){
      //	fx =0;
      //	fy =direction;
      //	}
      // end direct cells
      Cadherin cadherin_initial[ncad];
      //      for (int i=0; i<ncad; i++){
 	for (int kk=0; kk<ncad; kk++){
	  int i = kk;
	  cadherin_initial[i].location[0]=xc+6*cos(2*pi*kk/ncad);
	  cadherin_initial[i].location[1]=yc+6*sin(2*pi*kk/ncad);
	  cadherin_initial[i].length = 5.; // microns
	  cadherin_initial[i].force = 0;
	  cadherin_initial[i].mu = set_mu;
		//cout << "setmu initial"<< set_mu<< "\n";
	  cadherin_initial[i].kspring[0] = set_spring[type]*cad_spring_factor[type];
	  cadherin_initial[i].kspring[1] = diff_spring;

	  cadherin_initial[i].attach = -1;    
	  cadherin_initial[i].time = 0;
	  //cadherin_initial[i].front_indicator = 0;
	}
	//     }

	// Initialize integrins
       Integrin integrin_initial[nint];
      //      for (int i=0; i<ncad; i++){
 	for (int kk=0; kk<nint; kk++){
	  int i = kk;
	  integrin_initial[i].location[0]=xc;// set them at the cell center
	  integrin_initial[i].location[1]=yc;
	  integrin_initial[i].length = 5.; // microns
	  integrin_initial[i].force = 0;
	  integrin_initial[i].mu = set_mu;
	  integrin_initial[i].kspring = set_spring[type];

	  integrin_initial[i].attach = -1;    
	  integrin_initial[i].time = 0;
	  //integrin_initial[i].front_indicator = 0;
	}

	cell[ii].initialize(xc,yc,cadherin_initial,integrin_initial,fx,fy,type);

	  //	  cell[ii].print_output(fout,0,0,col);
	  //      for(int jk=0; jk<10; jk++){
	  //	  int nni =0;
	  //	  int nnj= 0;
	  //	  cell[ii].get_cadherin_node(jk,nni,nnj);
	  //	}
    }    
  }
  //october
  fout << 1.0*totx/ncells << " ";

}

/*************************************************************************************
                         localize_node_grid_update_all
           This function assigns each square grid all the cells in it. 
*************************************************************************************/

void localize_node_grid_update_all(Cell *cell)
{
  int l, m, i, j;

  /*  Set localize_node_grid to -1 */

  for (l = 0; l < lnx; l++)
    for (m = 0; m < lny; m++){
      for (int n = 0; n < lgrd; n++)
	  localize_node_grid[l][m][n] = -1;
    }

  for (int icell = 0; icell < ncells; icell++){
      double vector[2];
      cell[icell].get_center(vector,2);

      /* determine which grid the cell is currently in */
      i =static_cast<int>((vector[0]+shift) / lng_scale) ;
      j = static_cast<int>((vector[1]+shift) / lng_scale);

  	if (i > lnx-1)
  	  i = lnx -1;
  	if (i < 0)
  	  i = 0;
  	if (j > lny-1)
  	  j = lny -1;
  	if (j <0)
  	  j = 0;

      /* add the node to the grid it belongs to */
      int k = 0;
      //cout <<localize_node_grid[i][j][k]<<"\n";
      while (localize_node_grid[i][j][k] != -1 && k<lgrd)
  	{
  	  //   cout <<"bb  lgrd "<<lgrd<<" "<<k<<"\n";
  	  k = k+1;
  	}
      localize_node_grid[i][j][k] = icell;
      if(k == lgrd-1){
		 cout << k<<" lngrid is too small \n";
      }
     


  }          // end of icell loop

}            // end of the function definition



/*************************************************************************************
                              Nonlinear ODEs to solve

This function defines the RHS of the ODE system. 
dv/dt = F_spring - drag = F_spring - v * mu
By neglecting inertia, one get dx/dt = 1/mu * F_spring

*************************************************************************************/


int  move_nodes (double t, N_Vector cc, N_Vector fval, void *f_data){
  int nn=NV_LENGTH_S(cc);
  
  double u[nn];
  double savf[nn];
  Cell *cell;

  
  for (int i=0; i<nn; i++){// initialize the arrays since we only use part of them
    u[i]=0;
    savf[i]=0;
  }

  cell = (Cell*)f_data;

  for (int i=0; i<nn; i++){
    u[i]=NV_Ith_S(cc,i) ;
    //   cout<<" nn "<<nn<<" u "<<u[i]<<" "<<NV_Ith_S(cc,i)<<"\n";
  }

  function_node(u,savf,nn,cell);
//September save the forces on each cell to global variable. Ideally, we'd only do the last time of the iterative solver.
  for (int i = 0; i < 2*cx*cy; i++){
  	cell_forces[i] = savf[i];
  } 
  for (int i=0; i<nn; i++){
    NV_Ith_S(fval,i) = savf[i];
    //     cout<<" oososo "<<savf[i]<<" "<<NV_Ith_S(fval,i)<<"\n";
    }
  return(0);
}
void   function_node(double u[], double savf[],int nn,Cell *cell){
  // zero out the forces
  //   double vector[2];
      double vector2[2];
  int dim =2;
    double tmpx;
    double tmpy;
    double dist;
    double mur ;
    mur = 1./set_mu;
  for (int i = 0; i<nn; i++)
    savf[i]= 0.;
  for (int i =0; i<ncells; i++){
    int ii = 2*i;
    // push away from the substrate
    if(u[ii+1]<5){ //substrate is at zero and add 5 to push away.  For no substrate move to -500
      savf[ii+1] = savf[ii+1] + mur*10*cell[i].get_cadherin_kspring(1,0)*(exp(5-u[ii+1])-1); 
    }
    //savf[ii+1] = (savf[ii+1] - mur*.00262989*40); //Force of gravity on each cell //commented out sept 
    //  cout<<"sav "<<savf[ii]<<" "<<savf[ii+1]<<" ii "<<i<<" cell"<<u[ii+1]<<"\n";
    //  cout<<" vec "<<vector[0]<<" "<<vector[1]<<" "<<u[ii]<<" "<<dt<<" "<<set_mu<<"\n";
    //  cout <<"  d"<<set_mu*(u[ii]-vector[0])/dt<<"\n";
    // add a body spring to keep the cell centers at least 5 microns away from each other
    for (int j =0; j<ncells; j++){
      if(j != i){
	//	cell[j].get_center(vector2,2);
	//	tmpx = u[ii]-vector2[0];
	//	tmpy = u[ii+1]-vector2[1];
	tmpx = u[ii]-u[2*j];
	tmpy = u[ii+1]-u[2*j+1];
	dist = sqrt(tmpx*tmpx+tmpy*tmpy);
	if(dist<6){// push cells away from each other
	  savf[ii] = savf[ii] + mur*10*cell[i].get_cadherin_kspring(1,0)*(exp(6-dist)-1)*tmpx/dist;
	  savf[ii+1] = savf[ii+1] +mur*10*cell[i].get_cadherin_kspring(1,0)*(exp(6-dist)-1)*tmpy/dist;
	  //              savf[ii] = savf[ii] + (exp(10-dist)-1)*tmpx/dist;
	  //	  	  savf[ii+1] = savf[ii+1] + (exp(10-dist)-1)*tmpy/dist;
	  //	  	  savf[ii] = savf[ii] + 10*cell[i].get_cadherin_kspring(1)*(10-dist)*tmpx/dist;
	  //	  	  	  savf[ii+1] = savf[ii+1] + 10*cell[i].get_cadherin_kspring(1)*(10-dist)*tmpy/dist;
	}
      }//endif not same cell
    }// end loop for j through other cells
  }

  // Do the integrin forces *********************
  int ntmpint = cell[0].get_nintegrin(2);
   for (int i =0; i<ntmpint; i=i+dim){
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
    if(tmplength>1e-8){
      tmplength = 1/tmplength;
      }
    else{
      tmplength = 0; 
    }
    int dummy;
    int type = cell[icell].get_type(dummy);
    double tmptype;
    tmptype = int_spring_factor[type];


    double tmpfx = cell[icell].get_integrin_kspring(inode)*stretch*tmpx*tmplength*tmptype; // hookes law
    double tmpfy = cell[icell].get_integrin_kspring(inode)*stretch*tmpy*tmplength*tmptype; // hookes law
    


    //if(stretch<0){
    //  tmpfx=0.;
    //  tmpfy=0.;
    //}

    //  cout<<"tmpfx "<<tmpfx<<" icell "<<icell<<"itmp"<<itmp<<"inode"<<inode<<"\n";
    //  cout<<"tmpfy "<<tmpfy<<" icell "<<icell<<"itmp"<<itmp<<"inode"<<inode<<"\n";
    //  cout<<" wwl "<<(i-2*ncells)/2<<"i "<<i<<" ncells"<<ncells<<"\n";
    // add to the cell center
    // check if cadherin is attached to surface or other cadherin
    //    if(attach==-2){// if cadherin spring is compressed set the forces to zero
      //      if(stretch<0){
	//        tmpfx=0.;
	//	tmpfy=0.;
      //      }
      
    //    }

    savf[2*(icell)] = savf[dim*(icell)] + mur*tmpfx;
    savf[2*(icell)+1] = savf[dim*(icell)+1] + mur*tmpfy;
	//    savf[2*(icell)] = savf[2*icell]+100;
	//        savf[2*(icell)+1] = savf[2*icell+1]+0;
   }//end integrin forces
  // Do the cadherin forces ************************************
  int ntmp = cell[0].get_ncadherin(2);//gets cell's general cadherin # 
   for (int i = dim*ncells; i<ntmp+dim*ncells; i=i+dim){//dim = 2
    int itmp = cadherins_update[(i-dim*ncells)/2];// identifies cadherin to be updated ???
    int inode = itmp%ncad;//local cadherin #
    int icell = static_cast<int>((itmp - inode)/ncad);//cell #
    //    double vector[2];
    //    cell[icell].get_cadherin_location(inode,vector,2);
    // add spring constant terms
    double tmpx = u[i]-u[dim*(icell)];
    double tmpy = u[i+1]-u[dim*(icell)+1];
    double tmplength = tmpx*tmpx+tmpy*tmpy;
    tmplength = sqrt(tmplength);
    double stretch = tmplength -5; // assumes normal cadherin distance is 5 microns
    if(tmplength>1e-8){
      tmplength = 1/tmplength;
      }
    else{
      tmplength = 0;
    }
    
    int attach = cell[icell].get_cadherin_attach(inode);//returns general cadherin # of other cell's cadherin it's attached to
    int ainode = attach%ncad; 
    int aicell = static_cast<int>((attach - ainode)/ncad);
	   //int aattach = cell[aicell].get_cadherin_attach(ainode);
	   //cout << "attach:"<<attach<<" "<< "Aattach:"<<aattach<<"\n";
    int tmptype = 0;
    int dummy;
    int icelltype = cell[icell].get_type(dummy);       //0 is spore cell, 1 is stalk 
    int aicelltype = cell[aicell].get_type(dummy);
    if(icelltype!=aicelltype){
      tmptype = 1;
	}


    double tmpfx = cell[icell].get_cadherin_kspring(inode, tmptype)*stretch*tmpx*tmplength; // hookes law
    double tmpfy = cell[icell].get_cadherin_kspring(inode, tmptype)*stretch*tmpy*tmplength; // hookes law
 
    //if(stretch<0){
    //  tmpfx=0.;
    //  tmpfy=0.;
    //}
    //  cout<<"tmpfx "<<tmpfx<<" icell "<<icell<<"itmp"<<itmp<<"inode"<<inode<<"\n";

    //  cout<<"tmpfy "<<tmpfy<<" icell "<<icell<<"itmp"<<itmp<<"inode"<<inode<<"\n";
    //  cout<<" wwl "<<(i-2*ncells)/2<<"i "<<i<<" ncells"<<ncells<<"\n";
    // add to the cell center
    // check if cadherin is attached to surface or other cadherin
    
    //    if(attach==-2){// if cadherin spring is compressed set the forces to zero
      //      if(stretch<0){
	//        tmpfx=0.;
	//	tmpfy=0.;
      //      }
      
    //    }

    savf[2*(icell)] = savf[dim*(icell)] + mur*tmpfx;
    savf[2*(icell)+1] = savf[dim*(icell)+1] + mur*tmpfy;
	//    savf[2*(icell)] = savf[2*icell]+100;
	//        savf[2*(icell)+1] = savf[2*icell+1]+0;
	  //move cadherins
	  savf[i] = savf[i] - mur*tmpfx/cadviscosityfactor;
	  savf[i+1] = savf[i+1] - mur*tmpfy/cadviscosityfactor;


	  // add spring constant terms
	  tmpx = u[i]-u[dim*(aicell)];
	  tmpy = u[i+1]-u[dim*(aicell)+1];
	  tmplength = tmpx*tmpx+tmpy*tmpy;
	  tmplength = sqrt(tmplength);
	  stretch = tmplength -5; // assumes normal cadherin distance is 5 microns
	  if(tmplength>1e-8){
	    tmplength = 1/tmplength;}
	  else{
	    tmplength = 0;
	  }
	  tmpfx = cell[aicell].get_cadherin_kspring(ainode, tmptype)*stretch*tmpx*tmplength; // hookes law
	  tmpfy = cell[aicell].get_cadherin_kspring(ainode, tmptype)*stretch*tmpy*tmplength; // hookes law
	  
	  savf[i] = savf[i] - mur*tmpfx/cadviscosityfactor;
	  savf[i+1] = savf[i+1] - mur*tmpfy/cadviscosityfactor;
	  // add to the cell center
	  savf[2*(aicell)] = savf[dim*(aicell)] + mur*tmpfx;
	  savf[2*(aicell)+1] = savf[dim*(aicell)+1]  + mur*tmpfy;
	  //  cout<<savf[2*(icell)+1]<<"  "<<icell<<" "<<mur*tmpfy/.001<<"\n";	
	  //  cout<<tmpfx<<"_ "<<inode<<"_ "<<icell<<"\n";	  	
   } // end cadherin forces	

}
   

/*
 * Check function return value ...
 * opt == 0 means SUNDIALS function allocates memory so check if
 * returned NULL pointer
 * opt == 1 means SUNDIALS function returns a flag so check if
 * flag >= 0
 * opt == 2 means function allocates memory so check if returned
 * NULL pointer
 */
 int check_flag (void * flagvalue , char * funcname , int opt )
{
 int * errflag ;
/* Check if SUNDIALS function returned NULL pointer - no memory allocated */
 if (opt == 0 && flagvalue == NULL ) {
   fprintf ( stderr ,
	     "\n SUNDIALS_ERROR : %s()  failed  -  returned   NULL   pointer \n\n",
	     funcname );
   return (1);
 }
 /* Check if flag < 0 */
 else if ( opt == 1) {
   errflag = ( int *) flagvalue ;
   if (* errflag < 0) {
     fprintf (stderr ,
	      "\n SUNDIALS_ERROR : %s()  failed   with   flag  = %d \n\n",
	      funcname , * errflag );
     return (1);
   }
 }
 /* Check if function returned NULL pointer - no memory allocated */
 else if ( opt == 2 && flagvalue == NULL ) {
   fprintf ( stderr ,
	     "\n MEMORY_ERROR : %s() failed  -  returned   NULL pointer \n\n",funcname );
   return (1);
 }
 return (0);
}
/*
//Define the signal function durations is 3 minutes
double signal(double time)
{
  double value;
  if(time>maxsignalduration/2.)
    value=2.*(maxsignalduration-time)/(maxsignalduration)+.01;
  else
    value=2.*time/maxsignalduration;
  return(value);
}
*/
