#ifndef _DECVAR_H_
#define _DECVAR_H_

#define frontflg front_indicator 
#define updateflg time
#define attachflg attach
#define totalcad upsize
#define totalint upsizeint
 
// cell param
static const int nint = 10;
static const int ncad = 15;
//double basic_length_cad = 30; // in micrometer
//double basic_length_int = 30; // in micrometer 

static const double Spore_spring = 2.6e3*1.25; // spring constant SPORE
static const double Stalk_spring= 2.6e3*1.25; // spring constant STALK
static const double diff_spring= 2.6e3*1.25; // spring constant when they are different types
static const double set_spring[2]= {Spore_spring, Stalk_spring};

static const double cad_spring_spore = 1;
static const double cad_spring_stalk = 1;
static const double cad_spring_factor[2] = {cad_spring_spore, cad_spring_stalk};

static const double int_spring_spore = 1; //Make these 0 to eliminate integrin forces
static const double int_spring_stalk = 1;
static const double int_spring_factor[2] = {int_spring_spore, int_spring_stalk};

static const double set_mu = .2714; // will be .101788 it was .2714; // 70.e6 cell center leaves body 2.1 is drag coefficient for 1.5 (was 4) micron radius sphere in water. (?)

//double prob_attach_substrate//= .9; // .3  prob that the cadherin binds to substrate  l_2
static const double prob_change = .5; // .8 prob that cadherins bind to cadherin l4
//static const double cad_prob_change_from_front_to_back = .5;// .2 change direction from front to back; 
/*
static const double cad_prob_change_from_front_to_back_blue = .20;// .2 change direction from front to back;
static const double cad_prob_change_back_to_front=.80;// .6 change direction from back to front;WE ARE USING THIS ONE NOW FOR BLUE
static const double prob_change_from_front_to_back = .2;// .2 change direction from front to back; is the code really using this? -no, all ints are infront I think
static const double prob_change_back_to_front=.70;// .6 change direction from back to front: WE ARE USING THIS NOW TOO FOR GREEN
*/
static const double maxlength = 20;


//static const int dalcff=2;  // front indicator and factor to alterforce cannot be 1



/* cell layout */ 
static const int cx = 20;
static const int cy = 8; 
static const int ncells = cx*cy;
static const int totalint = ncells*nint; // total number of integrins
static const int totalcad = ncells*ncad; // total number of cadherins

/* grid layout
   lng_scale is the grid size, 
   lnx, lny are the number of grid points in the x, y dimensions
   lgrd is the maximum number of integrins a square grid can accomodate 
   lgrd is the maximum number of integrins a square grid can accomodate */
static const double lng_scale = 12; // in microns
static const int lnx = 10000;
static const int lny = 100;
static const int lgrd = 100;

//static const int lgrdh = 100; 

/* time integration */
static const double dt = .0166/60;  // in hours 
static const double tend = 1;// in hours

/* other constants */
static const double pi = 3.14159265358979323846264338327950288419716939937510;


/* externs */
extern int localize_node_grid[lnx][lny][lgrd];
//extern int localize_hand_grid[lnx][lny][lgrdh];
extern double basic_length_cad; // in micrometer
extern double basic_length_int; // in micrometer
extern double basic_length_back_factor;
extern double prob_attach_initial_int;
extern double prob_attach_initial;
//extern double prob_factor_back;// multiplied to prob_attach_substrate if in back l_2 for back
static const double shift = 10.;// shifts the localize_node_grid away from the edges since substrate is at 0

extern double cad_averand;
extern double cad_avedir[2];
extern double cad_numang ;

extern double int_averand;
extern double int_avedir[2];
extern double int_numang ; 

extern double avejindex[9];


#endif
