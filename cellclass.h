// Cell class define
#ifndef _CELL_H_
#define _CELL_H_
//#include <iostream>
//#include <fstream>
using namespace std;
#include <string>
#include <iostream>

using std::string;

#include "decvar.h"

struct Cadherin
{
  double location[2]; // the coordinates of the cadherin
  double vector[2];   // the vector starting at the cell center and ending at the cadherin, i.e. location-center
  double length;      // the length of vector
  double force;       // the magnitude of the force on the element
  double kspring[2];     // the spring constant
  double mu;          // the dashpot coefficient
  double time;         // time for the cadherin to dettach
  int attach;         // -1 means it is free and 0 and up means attached to cadherin (indicates cellnumber*ncad+cadherin_index to identify cadherin to which it is attached)
 // int front_indicator; // dalcff means it is in the front and 1 means it is in the back
	
};
struct Integrin 
{
  double location[2]; // the coordinates of the integrin
  double vector[2];   // the vector starting at the cell center and ending at the integrin, i.e. location-center
  double length;      // the length of vector
  double force;       // the magnitude of the force on the element
  double kspring;     // the spring constant
  double mu;          // the dashpot coefficient
  double time;         // time for the integrin to dettach
  int attach;         // -1 means it is free and 0 means attached to substrate
//  int front_indicator; // 0 means it is in the front and 1 means it is in the back
};
class Cell
{
 private:
  double center[2];         // the coordinates of the center of the cell
  Cadherin cadherin[ncad];  // the cadherins
  Integrin integrin[nint];  // the integrins
 // double front[2];         // the vector indicating the front of the cell
  int ncadherin; // a variable the same for all cells which is used in move_nodes
  int nintegrin; // number of integrins attached to substrate
  int celltype; //0 is spore cell, 1 is stalk 

  void calculate_cadherin_length(int cadherin_index); 
  void calculate_cadherin_force(int cadherin_index);    
  void calculate_cadherin_vector(int cadherin_index);    

  void calculate_integrin_length(int integrin_index); 
  void calculate_integrin_force(int integrin_index);    
  void calculate_integrin_vector(int integrin_index);    
  //void calculate_integrin_front_indicator(int integrin_index);

  void update_cadherin_attachment(int cadherin_index, double dt);
  void update_cadherin_attachment_new(int cadherin_index, double dt,double t, double average_spore,double average_stalk[], double average_different, int restart, double average_unbound_cad_spore, double average_unbound_cad_stalk,int cell_index, Cell *cell);
  void update_cadherin_attachment_initial(int cadherin_index, double dt,double average_spore,double average_stalk[], double average_different,int restart, double average_unbound_cad_spore, double average_unbound_cad_stalk,int cell_index, Cell *cell);

  void update_integrin_attachment(int integrin_index, double dt);
  void update_integrin_attachment_new(int integrin_index, double dt,double t, double average_substrate_spore,double average_substrate_stalk[], int restart,int cell_index, Cell *cell);
  void update_integrin_attachment_initial(int integrin_index, double dt,double average_substrate_spore,double average_substrate_stalk[],int restart, double average_unbound_int_spore, double average_unbound_int_stalk,int cell_index, Cell *cell);

 public:
  void initialize( double xc, double yc, Cadherin init_cad[], Integrin init_int[], double frontx, double fronty, int celltype);
  Cell(double x, double y, Cadherin initial_cad[],Integrin init_int[], double front1, double front2, int celltype1);
  Cell();
  ~Cell();
  int restart;

 // void calculate_cadherin_front_indicator(int cadherin_index);


  void get_center(double value[], int lim);
 // void get_front(double value[], int lim);

  double get_cadherin_length(int cadherin_index);
  double get_cadherin_force(int cadherin_index);
  double get_cadherin_kspring(int cadherin_index, int tmptype);
  double get_cadherin_mu(int cadherin_index);
  int get_ncadherin(int dummy);
  int get_cadherin_attach(int cadherin_index);
  double get_cadherin_time(int cadherin_index);
  void get_cadherin_location(int cadherin_index, double value[], int lim);
  void get_cadherin_vector(int cadherin_index, double value[], int lim);
 // int get_cadherin_front_indicator(int cadherin_index);
  double get_integrin_length(int integrin_index);
  double get_integrin_force(int integrin_index);
  double get_integrin_kspring(int integrin_index);
  double get_integrin_mu(int integrin_index);
  int get_nintegrin(int dummy);
  int get_integrin_attach(int integrin_index);
  double get_integrin_time(int integrin_index);
  void get_integrin_location(int integrin_index, double value[], int lim);
  void get_integrin_vector(int integrin_index, double value[], int lim);
//  int get_integrin_front_indicator(int integrin_index);
  int get_type(int type);

  void set_center(double value[], int lim);
  //void set_front(double value[], int lim);
  void set_type(int value);

  void set_cadherin_attach(int cadherin_index, int value);
  void set_cadherin_time(int cadherin_index, double value);
  void set_ncadherin(int value);
  void set_cadherin_location(int cadherin_index, double value[], int lim);

  void set_integrin_attach(int integrin_index, int value);
  void set_integrin_time(int integrin_index, double value);
  void set_nintegrin(int value);
  void set_integrin_location(int integrin_index, double value[], int lim);


  void move_center(double dt);
  void print_output(ofstream *fout, int n, double dt);
  void print_output_force(ofstream *fout,int n, double dt);
  void update_cadherins(double dt, double t,double average_spore,double average_stalk[], double average_different,int restart,double average_unbound_cad_spore, double average_unbound_cad_stalk, int inode,int cell_index,Cell *cell);
  void update_cadherins_initial(double dt, double average_spore,double average_stalk[], double average_different,int restart, double average_unbound_cad_spore, double average_unbound_cad_stalk,int ni, int index_cell,Cell *cell);
  void update_integrins(double dt, double t,double average_substrate_spore,double average_substrate_stalk[],int restart, int inode,int cell_index,Cell *cell);
  void update_integrins_initial(double dt, double average_substrate_spore,double average_substrate_stalk[],int restart,double average_unbound_int_spore, double average_unbound_int_stalk,int ni, int index_cell,Cell *cell);
};
#endif
