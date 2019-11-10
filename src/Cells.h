#include <getopt.h>
#include <stdlib.h> /*free, exit, calloc, malloc, atoi, strtof*/
#include <math.h>   /*isfinite, sqrt, log*/
#include <stdio.h>  
#include <string.h> /*strerr, strlen, strcat */
#include <errno.h>  /*errno*/
#include <limits.h> /*ULONG_MAX INT_MAX*/
#include "Queue.h"
#include "mersenne.h"
#include "cartesian.h"

#define MAX_LINE_LEN 100
#define NO_IDX ULONG_MAX
#define N_CELL_TYPES 4
enum nodeType_t   { NORMAL, CANCER, INFECT_INT, INFECT_EXT, VACANT, N_NODE_TYPES };
enum cellAction_t { KILL, GROW, N_ACTIONS };
enum visited_t    { NOTSEEN, SEEN, VISITED };
enum infection_t  { NONE, RANDOM, CENTER, PERIMETER, MULTINODE };

/*******************/
/*** Structures  ***/
/*******************/
/* Event */
struct event_t {
   enum cellAction_t cell_Action;
   enum nodeType_t node_Type;
};

/* Cell Lattice */
struct lattice_t {
   unsigned long int *n_Target_Nodes[N_CELL_TYPES]; 
                               //2nd dim size: n_Nbrs_Max+1
   struct node_t *** target_Nodes[N_CELL_TYPES]; 
                               //size:[4][n_Nbrs_Max+1][n_Nodes]
                               //[SOURCE_Type][#target neighbors][cell addy]

   unsigned long int  n_Cells[N_CELL_TYPES]; // # of each TYPE of Cell
   struct node_t ** cells[N_CELL_TYPES]; //size: [n_Nodes];

   unsigned long int    n_Nodes;
   struct node_t ** nodes;

   double  time;
   unsigned long int count;  /* # of events since lattice initialization */
   unsigned long int center_Idx;
   unsigned int    n_Nbrs_Max;
   unsigned int    n_Nbrs_Avg;
   char*  coord_Filename;
   int    n_Dimensions;
};

/* Cell */
struct node_t {
   struct node_t ** neighbors; //Array size: n_Nbrs_Max
   int n_Neighbors;
   Point *pos;

   enum nodeType_t type;
   int interior; //indicates Infected cells position in tumor
   int n_Nbrs[N_CELL_TYPES]; //Number of neighbors of type

   unsigned long int targets_Idx[N_CELL_TYPES]; 
               // Locations(s) in target_Nodes (NO_IDX == not present)
   
   unsigned long int list_Idx;  
               // Cell's location in cells[this.type] array NO_IDX if empty

   unsigned long int Idx;       
               // Cells' location in lttc->nodes array

   /* variables used for BFS */
   enum visited_t visited;          
   unsigned int distance;
};


/*=================*/
/*=== Functions ===*/
/*=================*/
float  argtof( const char *nptr );
enum   nodeType_t TargetOf( const enum nodeType_t type );
void LinReg( float const * const X, 
             float const * const Y, 
	     long unsigned int const N,
	     float * slope,
	     float * intercept );

int  SeedCancer( struct node_t* center, unsigned long int cancer_Count, 
                                        unsigned int cancer_Radius );
int  InjectVirus( );        
int  VirusBallAroundNode( struct node_t* node, unsigned long int n_Infect );
int  Simulate( );
int  OutputStats();
int  OutputSummary();
struct node_t* DetermineEventRecipient( const struct event_t event);
struct event_t DetermineEventType( );
void   ChangeNodeType( struct node_t* node, enum nodeType_t type );
int  VerifyState( );
void OpenFiles();
void ParameterFile(float time_Extinct);
void CloseFiles();

/*****************************/
/* Network Related Functions */
/*****************************/
struct lattice_t* NewLattice( unsigned long int n_Nodes, 
                unsigned int n_Dimensions,
                unsigned int n_Nbrs_Max, 
                unsigned long int center_Idx,
                enum nodeType_t new_Type ) ;
Point* CenterOfMass( struct node_t** nodes, unsigned long int N );
struct node_t* FindNodeNear( Point* point, struct node_t ** nodelist, 
                                           unsigned long int N);
int    BFS (struct node_t* node_u, 
            enum nodeType_t network_Type,
	    int max_Count,
	    enum nodeType_t new_Type );
void   RemoveCell (unsigned long int* const del_Idx, 
                   struct node_t** array, 
		   unsigned long int* const n_cells, 
		   unsigned long int* const last_Idx);
void   AddCell( struct node_t* const cell, 
                struct node_t** array, 
		unsigned long int* const n_array, 
		unsigned long int* const cell_Idx );
unsigned long int*   GetNetworkDetails( struct node_t ** nodes, 
                                        unsigned long int N );
struct node_t* FindEdge( struct node_t** nodes, unsigned long int N, 
                         struct node_t* start, double dir );
void   DeleteLattice( struct lattice_t* lattice );
int    ReadDataFile( const char* const data, 
                  unsigned long int * n_Nodes, 
                  int* n_Dimensions,
                  unsigned int * n_Nbrs_Max,
                  unsigned long int * center_Idx ) ;
/* Uses Global lttc */
void   OutputGrid();
void   OutputLatticeState( char* filename ); 
int ReadNetworkFile( const char* const network, 
                enum nodeType_t new_Type ) ;
int ReadCoordinateFile(  const char* const coords );
int ReadLattice(const char* const data, 
                const char* const network, 
		const char* const coords );

/***********************/
/* MergeSort Functions */
/***********************/
void MergeSort ( double** list, unsigned long int N );
void MergeSplit( double** list, unsigned long int b_idx, 
                                unsigned long int e_idx, double** temp );
void Merge     ( double** list, unsigned long int b_idx, 
                                unsigned long int m_idx, 
		                unsigned long int e_idx, double** temp );

/****************/
/***Parameters***/
/****************/
/* FIXME Should these all be global? */
extern struct lattice_t * lttc;
extern FILE*  vlpt;
extern FILE*  vpts;
extern FILE *fpp, *fpl;

/* Pass from Main to Simulate */
extern float output_Interval; /* frequency of lattice output */
extern int   print_Stats;
extern int   print_Summary;
extern int   print_Grid;
extern int   print_Tumor_Sizes;
extern int   print_Watch;
extern char* print_File;
extern int   print_Files;
extern int   verify_State;
extern int   rng_Seed;
extern int   t_Max_Simulation; 
extern int   t_Infect; 
extern int   do_Equalization;


/* Pass from Main to Simulate() to InjectVirus() */
extern double percent_Infected;
extern double percent_Interior;
extern enum infection_t infect_Type;

/* Pass from Main to Simulate() to DetermineEventType() */
extern float  rates[N_CELL_TYPES][N_ACTIONS];
