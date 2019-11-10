#include "Cells.h"

struct lattice_t * lttc = NULL;
char *program_name = NULL;
int print_Stats = 0;
int print_Summary = 0;
int print_Grid = 0;
int print_Tumor_Sizes = 0;
int print_Watch = 0;
char *print_File;
int print_Files = 0;
int verify_State = 0;
double percent_Infected = 0.0;
double percent_Interior = 100;
int   do_Equalization = 0;
enum infection_t infect_Type = NONE;
FILE *fpp, *fpl;

int    t_Max_Simulation = 75;
int    t_Infect         =  7;
int    rng_Seed = 11;
float  output_Interval=.5; /* frequency of lattice output */
float  rates[N_CELL_TYPES][N_ACTIONS] = { { 0, 0 },
                                          { 0, 0 },
/*  Virus Interior */                     { 0, 0 },
/*  Virus Exterior */                     { 0, 0 } };

/* 
 * CenterOfMass() accepts a pointer to an array of N "node_t nodes". 
 * Returns NULL if dimensions don't match or are != 2 || 3 
 * finds the average of each node's position (nodes[i]->pos) coordinates.
 * Stores the coordinates in a Point
 * Returns the address of this point. 
 */
Point* CenterOfMass( struct node_t** nodes, const unsigned long N )
{
   unsigned long nodeIt;
   int n_dims = nodes[0]->pos->N;
   int dimIt;
   double *temp = calloc( n_dims, sizeof(double) );

   if( n_dims != 2 && n_dims != 3 )
      return NULL;

   for( nodeIt=0; nodeIt < N; nodeIt++ ) {
      for( dimIt=0; dimIt < n_dims; dimIt++) {
         temp[ dimIt ] += nodes[nodeIt]->pos->data[dimIt];
      }
   }
   for( dimIt=0; dimIt < n_dims; dimIt++) {
      temp[dimIt] /= N;
   }

   Point* center = NULL;
   switch( n_dims ) {
   case 2:
      center = NewPoint( n_dims, temp[0], temp[1] );
      break;
   case 3:
      center = NewPoint( n_dims, temp[0], temp[1], temp[2] );
      break;
   default:
      break;
   }

   free( temp );

   return center;
} /* End of CenterOfMass() */


struct node_t* FindNodeNear(Point *point,  
                            struct node_t ** node_list, 
                          unsigned long N )
{
/* Return address of node in node_list with shortest euclidian distance 
   to point */
   unsigned long nodeIt;
   unsigned long min_Idx = 0;
   float dx=0;
   float dy=0;
   float dz=0;
   float *distance = malloc( N * sizeof( *distance) );

   for( nodeIt = 0; nodeIt < N; nodeIt++){
      dx = point->data[0] - node_list[nodeIt]->pos->data[0];
      dy = point->data[1] - node_list[nodeIt]->pos->data[1];
      if( point->N == 3 ){
         dz = point->data[2] - node_list[nodeIt]->pos->data[2];
      }
      distance[nodeIt] = (dx*dx + dy*dy + dz*dz);
   }

   for( nodeIt = 1; nodeIt < N; nodeIt++){
      if( distance[nodeIt] < distance[min_Idx] )
         min_Idx = nodeIt;
   }

   free(distance);

   return node_list[min_Idx];
} /* End of FindNodeNear() */

/* Perform a Breadth First Search of all nodes connected to node_u by 
      network_Type neighbors.  
   Assume connected nodes have visited==NOTSEEN && distance==INT_MAX 
   Sets node->visited to VISITED and node->distance to shortest number
      of links between node_u and node. 
   Returns the number of nodes visited */
int BFS (struct node_t* node_u,        /* Node to begin search from */
         enum nodeType_t network_Type, /* Only visit type nbrs */
         int max_Count,                /* Stop after visiting max_Count */
         enum nodeType_t new_Type )    /* Change nodes to new_Type */ {
   int count=0;
   struct node_t* node_v;
   struct Queue * const Q = NewQueue( 8 );
   node_u->visited = SEEN;
   node_u->distance = 0;
   Enqueue( Q, node_u );

   while( (node_u = Dequeue(Q)) ){
      if( count >= max_Count ){
         continue;
      }
      count++;
      if( node_u->type != new_Type )
         ChangeNodeType( node_u, new_Type );
      node_u->visited = VISITED;
      int nbrIt;
      for( nbrIt = 0; nbrIt < node_u->n_Neighbors; nbrIt++ ) {
         node_v = node_u->neighbors[nbrIt];
         if( node_v->visited == NOTSEEN && node_v->type == network_Type ) {
            node_v->visited = SEEN;
            node_v->distance = node_u->distance+1;
            if( Enqueue( Q, node_v ) == 0 )
               printf(" Problem growing Queue\n" );
         }
      }
   }
   DeleteQueue(Q);
   return count;
}

/* NextNetwork() returns first of "N" nodes in "nodes" that has 
     node->visited set to NOTSEEN or NULL if none are found.
   Used along with BFS() to find disconnected portions of a network */
struct node_t* NextNetwork( struct node_t** nodes, unsigned long N ){
   int nodeIt;
   for( nodeIt = 0; nodeIt < N; nodeIt++ ) {
      if( nodes[nodeIt]->visited == NOTSEEN ){
         return nodes[nodeIt];
      }
   }
   return NULL;
}

/* GetNetworkDetails() ****************************************************
   Resets node->visited & node->distance for all cancer nodes in lttc.
   Returns pointer to newly allocated array that needs to be free()d.
   The first element of returned array is the number of tumors N.
   The second to Nth+1 element indicate the sizes of these tumors.
   Elements aren't sorted except by chance
**************************************************************************/
unsigned long* GetNetworkDetails( struct node_t ** nodes, unsigned long N )
{
   struct node_t *node_u;
   int nodeIt;
   int array_Size = 4;
   unsigned long* sizes = calloc( array_Size, sizeof(*sizes) );
   unsigned long* temp;
   if( sizes == NULL ) {
      fprintf(stderr, "Memory error while calculating tumor sizes.\n");
      return sizes;
   }
   if( N == 0 ){
      return sizes;
   }

   for( nodeIt = 0; nodeIt < N; nodeIt++ ){
      nodes[nodeIt]->visited = NOTSEEN;
      nodes[nodeIt]->distance = INT_MAX;
   }

   node_u = nodes[0];
   while( node_u ) {
      sizes[0] = sizes[0] + 1;
      if( sizes[0] >= array_Size-1 ){
         temp = realloc(sizes, array_Size*2 * sizeof(*temp) );
         if( temp == NULL ){
            fprintf(stderr, "Memory error calculating tumor sizes.\n");
            return sizes;
         }
         sizes = temp;
         array_Size *= 2;
      }
      sizes[ sizes[0] ] = BFS(node_u, nodes[0]->type, 
	                     INT_MAX, nodes[0]->type);
      node_u = NextNetwork( nodes, N );
   }

   return sizes;
}

float argtof( const char *nptr ) {
/* Checks for invalid values from strtof() */
   float num;
   num = strtof( nptr, NULL );

   if( !isfinite(num) || num != num ) {
      fprintf( stderr, "Argument is inf or nan. Using default of 0\n" );
      num = 0;
   }

   return num;
}/* End of argtof() */

/** LinReg() ***********************************************************/
/*  Calculate slope of the N points (Xi, Yi) using Linear Regression.  */
/*  If intercept isn't NULL, that is calculated as well.               */
/*  Formulas used are those derived in "Mathematical Statistics with   */
/*      Applications 6th Edition" by Wackerly on page 538.             */
/***********************************************************************/
void LinReg( float const * const X, 
             float const * const Y, 
	     long unsigned int const N,
	     float * slope,
	     float * intercept ) {
  int It;
  float Sxy=0;
  float Sx =0;
  float Sy =0;
  float Sxx=0;

  for(It = 0; It < N; It++) {
    Sxy += X[It]*Y[It];
    Sx  += X[It];
    Sy  +=       Y[It];
    Sxx += X[It]*X[It];
  }

  *slope = (Sxy - Sx*Sy/N)/(Sxx-Sx*Sx/N);
  if( intercept != NULL )
    *intercept = Sy/N - *slope * Sx/N;

  return;
}

/** simulate() ***********************************************************
*  Runs one step of the simulation
*  Expects global struct lattice_t lttc to be in a consistant state
*  Outputs data if indicated by global print_* variables
*************************************************************************/
int Simulate( ) {
   static float output_Timer = 0;
   static int count = 0;
   if( lttc->time == 0 ) {
      output_Timer = output_Interval;
   }

   const struct event_t e_t = DetermineEventType( );
   if( e_t.node_Type == VACANT )
   {
      fprintf(stderr, "Event Type is Vacant?\n");
      return -1;
   }
   struct node_t * const event_Cell = DetermineEventRecipient( e_t );

   if( event_Cell == NULL ) {
      fprintf(stderr, "Failed to select a cell to Grow or Kill\n");
      return -1;
   }

   switch ( e_t.cell_Action ) {
   case KILL:
      ChangeNodeType( event_Cell, VACANT);
      break;
   case GROW:
      ChangeNodeType( event_Cell, e_t.node_Type);
      break;
   default:
      fprintf(stderr, "e_t.cell_Action isn't KILL or GROW?\n");
      fprintf(stderr, "  Continuing w/o taking action or updating time.\n");
      return -1;
   }

   
   /* Output Statistics & Info */
   if( lttc->time > output_Timer ) {
      int newlines=0;
      Point* center = NULL;
      output_Timer += output_Interval;

      if( print_Files ){
         /*
            Files being printed contain rng_Seed & the number of cancer cells at 
            injection followed by a list of space separated quantities indicating 
            the number of Cancer & Infected cells at the intervals specified below.
              vpts.dat (output_Interval beginning at 1 ending at 126) && 
              vlpt.dat (Day 0 + 5th output_Interval beginning at 1)
         */
	 unsigned long int n_Tumor = lttc->n_Cells[CANCER]     +
				     lttc->n_Cells[INFECT_INT] +
				     lttc->n_Cells[INFECT_EXT] ;
         if( count == 5 ){
            count = 0;
            fprintf(fpl," %lu ", n_Tumor );
         }
         if( lttc->time < 126 ) {
            fprintf(fpp," %lu ", n_Tumor );
         }
         count++;
      }


      if( print_Tumor_Sizes && lttc->n_Cells[CANCER]){
         if( !newlines ) {
            printf("\n\nTime: %5.2f; ", lttc->time);
            newlines = 1;
         }
         unsigned long* sizes;
         int It;
         sizes = GetNetworkDetails(lttc->cells[CANCER], 
	                           lttc->n_Cells[CANCER]);
         printf("Number of Tumors: %2lu; Sizes:", sizes[0] );
         for( It = 1; It < sizes[0]+1; It++ ) {
            printf(" %4lu ", sizes[It] );
         }
         center = CenterOfMass( lttc->cells[CANCER], 
	                       lttc->n_Cells[CANCER] );
         printf("\nCenter of Mass: (");
         for( It = 0; It < lttc->n_Dimensions; It++ )
            printf("%f,", center->data[It]);
         printf(")\n");

         DeletePoint( center );

         free(sizes);
      }
      if( print_Grid ) {
         if( !newlines ) {
            printf("\n\nTime: %5.2f;\n", lttc->time);
            newlines = 1;
         }
         OutputGrid( lttc );
      }
      if( print_Stats ) {
         OutputStats( lttc );
      }
      if( print_Watch ){
         printf("Press Return.\n");
         char c = getchar();
         if(c == 'q')
            exit(-1);
      }
      if( print_Summary             &&
          lttc->n_Cells[CANCER] > 0 &&
          newlines == 0             &&
          print_Watch == 0  ) {
         OutputSummary( lttc );
      }
      if( print_File ){
         char filename[MAX_LINE_LEN];
         snprintf(filename, MAX_LINE_LEN, "%s%.1f.csv", 
	                                print_File, lttc->time);
         OutputLatticeState( filename );
      }

   }
   return 1;
}/* End of Simulate() */



/* Read parameters from data file formatted as: */
/* n_Nodes n_Dimensions n_Nbrs_Max [center_Idx] */
int ReadDataFile( const char* const data,
                  unsigned long * n_Nodes,
                  int* n_Dimensions,
                  unsigned * n_Nbrs_Max,
                  unsigned long * center_Idx ) {
   FILE* fh;
   int return_Val;
   fh = fopen( data, "r" );
   if( fh == NULL ) {
        fprintf(stderr, "Failed to open data file %s.\n", data);
        return( -1 );
   }
   return_Val = fscanf(fh, "%lu %d %u %lu",
                       n_Nodes,
                       n_Dimensions,
                       n_Nbrs_Max,
                       center_Idx);
   switch( return_Val ) {
   case 3:
      *center_Idx = 0;
      break;
   case 4:
      if ( *center_Idx > *n_Nodes) {
         fprintf(stderr, "Value for center > number of nodes.\n");
         fprintf(stderr, "Using 1 instead. \n");
         *center_Idx = 1;
      }
      /* decrement since file indexes are 1 based & code indices are 0 based */
      (*center_Idx)--;
      break;
   default:
      fprintf(stderr, "Problem reading from datafile %s\n", data);
      return( -1 );
   }
   fclose( fh );/* Finished reading datafile */

   return 0;
} /* End ReadDataFile() */


/* Allocate Space for Network and Lattice members */
struct lattice_t* NewLattice( unsigned long n_Nodes,
                unsigned n_Dimensions,
                unsigned n_Nbrs_Max,
                unsigned long center_Idx,
                enum nodeType_t new_Type ) {

   int nbrIt;
   int It;
   struct lattice_t *L;
   char *error_Msg = "Failed to Allocate memory for lattice\n";
   if( n_Dimensions != 2 && n_Dimensions != 3 ){
      fprintf(stderr, "Invalind number of dimensions\n");
      return(NULL);
   }

   L = calloc (1, sizeof *L);
   if( L == NULL ) {
      fprintf(stderr, "%s", error_Msg);
      return(NULL);
   }
   L->nodes = calloc( n_Nodes , sizeof( struct node_t * ) );
   if( L->nodes == NULL ) {
      fprintf(stderr, "%s", error_Msg);
      return(NULL);
   }
   for( It=0; It < N_CELL_TYPES; It++ ) {
     L->n_Cells[ It ] = 0;
     L->target_Nodes[It] = calloc( n_Nbrs_Max+1, sizeof(struct node_t**));
     if( L->target_Nodes[It] == NULL ) {
        fprintf(stderr, "%s", error_Msg);
        return(NULL);
     }
     L->n_Target_Nodes[It]  = calloc( n_Nbrs_Max+1, sizeof(unsigned long) );
     L->cells[It]        = calloc( n_Nodes, sizeof(struct node_t *) );
     if( L->n_Target_Nodes[It]   == NULL ||
         L->cells[It] == NULL ) {
        fprintf(stderr, "%s", error_Msg);
        return(NULL);
     }
     for( nbrIt=0; nbrIt <= n_Nbrs_Max; nbrIt++ ) {
        L->target_Nodes[It][nbrIt] 
	                     = calloc( n_Nodes, sizeof(struct node_t*));
        if( L->target_Nodes[It][nbrIt] == NULL ){
           fprintf(stderr, "%s", error_Msg);
           return(NULL);
        }
     }
   }

   /* Initialization of lattice */
   L->time              = 0;
   L->count             = 0;
   L->n_Nbrs_Avg        = 0;
   L->n_Dimensions      = n_Dimensions;
   L->center_Idx        = center_Idx;
   L->n_Nodes           = n_Nodes;
   L->n_Nbrs_Max        = n_Nbrs_Max;
   L->n_Cells[new_Type] = n_Nodes;


   return L;
}



/* Each line of "network" contains two fields
   The first is the number of neighbors that node has
   The second is a comma separated list of indices for each neighbor 
   The node's index is it's line number-1 
   The neighbor's index is also decremented by 1 to have 0-based indices */
int ReadNetworkFile( const char* const network, enum nodeType_t new_Type ) {
   FILE* fh;
   unsigned long nodeIt;
   int return_Val;
   unsigned n_Neighbors, nbrIt;

   fh = fopen( network, "r" );
   if( fh == NULL ) {
        fprintf(stderr, "Failed to open network file %s. %s\n",
              network, strerror(errno));
        return( -1 );
   }

   for( nodeIt = 0; nodeIt < lttc->n_Nodes; nodeIt++ ) {
      struct node_t *node = lttc->nodes[ nodeIt ]; 
      /* Get number of neighbors for node[ nodeIt ] */
      return_Val = fscanf( fh, "%u", &n_Neighbors );
      if( return_Val != 1 ) {
         fprintf(stderr, "Problem reading number of neighbors from network file %s\n", network);
         fprintf(stderr, "Line Number: %ld of %ld\n", nodeIt+1, lttc->n_Nodes);
         return( -1 );
      }
      lttc->n_Nbrs_Avg += n_Neighbors;
      node->n_Neighbors = n_Neighbors;
      node->n_Nbrs[new_Type] = n_Neighbors;
      node->neighbors = malloc( sizeof(struct node_t*) * n_Neighbors);
      if( node->neighbors == NULL ){
         fprintf(stderr, "Memory error while creating neighbors.\n");
         return(-1);
      }

      /* Populate neighbor lists */
      for( nbrIt = 0; nbrIt < n_Neighbors; nbrIt++ ) {
         unsigned long nbrIdx;
         return_Val = fscanf( fh, "%lu", &nbrIdx);
         if( return_Val != 1 ) {
            fprintf(stderr, "Problem reading network file %s\n", network);
            fprintf(stderr, "Populating neighbor list #%d\n", nbrIt+1);
            fprintf(stderr, "Node: %ld\n", nodeIt+1);
            return( -1 );
         }
	 /* decrement since file indexes are 1 based & code indices are 0 based */
         nbrIdx--;
         node->neighbors[ nbrIt ] = lttc->nodes[ nbrIdx ];
      }
   }
   fclose( fh );
   lttc->n_Nbrs_Avg /= lttc->n_Nodes;

   return 0;
}


/* Each line of file coords gives that node's cartesian coordinates */
/* The optional 3rd (or 4th) field gives that node's type */
int ReadCoordinateFile(  const char* const coords ){
   unsigned long nodeIt;
   struct node_t* cell;
   FILE* coords_fh = fopen( coords, "r" );

   if( coords_fh == NULL ) {
        fprintf(stderr, "Failed to open coordinate file %s. %s\n",
              coords, strerror(errno));
        return( -1 );
   }

   for( nodeIt = 0; nodeIt < lttc->n_Nodes; nodeIt++ ) {
      char line[MAX_LINE_LEN];
      char type = '0';
      double x, y, z;
      fgets(line, MAX_LINE_LEN, coords_fh);
      int r_val;

      cell = lttc->nodes[nodeIt];
      switch( lttc->n_Dimensions ) {
      case 3:
         r_val = sscanf( line, "%lf %lf %lf %c", &x, &y, &z, &type);
         cell->pos = NewPoint( lttc->n_Dimensions, x, y, z);
         break;
      case 2:
         r_val = sscanf( line, "%lf %lf %c", &x, &y, &type);
         cell->pos = NewPoint( lttc->n_Dimensions, x, y);
         break;
      default:
         return -1;
      }

      if( r_val == lttc->n_Dimensions ){
	 type = 1; 
      }


      switch( type ){
      case '0':
	 ChangeNodeType( cell, VACANT );
	 break;
      case '1':
         ChangeNodeType( cell, NORMAL );
         break;
      case '2':
         ChangeNodeType( cell, CANCER );
         break;
      case '3':
         ChangeNodeType( cell, INFECT_EXT );
         break;
      default:
         break;
      }
   }
   fclose( coords_fh );
   return 0;
}


struct node_t * NewNode ( unsigned long nodeIt,
                enum nodeType_t new_Type ) {
   int typeIt;
   struct node_t * node;
   node = calloc( 1, sizeof(struct node_t));
   if( node == NULL ){
      fprintf(stderr, "Failed to allocate memory for node %lu\n", nodeIt);
      return(NULL);
   }
   node->type     = new_Type;
   node->list_Idx = nodeIt;
   node->Idx      = nodeIt;
   node->visited  = NOTSEEN;
   node->distance = INT_MAX;
   for(typeIt = 0; typeIt < N_CELL_TYPES; typeIt++) {
      node->targets_Idx[typeIt] = NO_IDX;
      node->n_Nbrs[typeIt] = 0;
   }

   return node;
}

int ReadLattice(const char* const data, 
                const char* const network, 
                const char* const coords) {
   unsigned long nodeIt;
   const int new_Type = NORMAL;
   int return_Val;

   unsigned long center_Idx;
   int n_Dimensions;
   unsigned long n_Nodes;
   unsigned n_Nbrs_Max;

   /* Read Data File */
   return_Val = ReadDataFile( 
	          data, &n_Nodes, &n_Dimensions, &n_Nbrs_Max, &center_Idx);
   if( return_Val != 0 ) {
      fprintf(stderr, "Problem reading from data file %s\n", data);
      return( -1 );
   }

   /* Create Lattice */
   lttc = 
      NewLattice( n_Nodes, n_Dimensions, n_Nbrs_Max, center_Idx, new_Type );
   if( lttc == NULL ) {
      fprintf(stderr, "Problem creating lattice.\n");
      return( -1 );
   }

   /* Allocate Space for & Initialize nodes. */
   for( nodeIt = 0; nodeIt < n_Nodes; nodeIt++ ) {
      lttc->nodes[nodeIt] = NewNode( nodeIt, new_Type );
      if( lttc == NULL ) {
         fprintf(stderr, "Problem creating node %lu.\n", nodeIt);
         return( -1 );
      }
      lttc->cells[new_Type][nodeIt] = lttc->nodes[nodeIt];
   }

   /* Read Network File */
   return_Val = ReadNetworkFile( network, new_Type );
   if( return_Val != 0 ) {
      fprintf(stderr, "Problem reading from Network file %s\n", network);
      return( -1 );
   }

   /* Read Coordinate File */
   if( coords != NULL )
   {
      return_Val = ReadCoordinateFile( coords );
      if( return_Val != 0 ) {
	 fprintf(stderr, "Problem reading from Coordinate file %s\n", coords);
	 return( -1 );
      }
   }

   return 0;
}/* End of ReadLattice() */



/* AddNeighborsToTargetNodes() *************************************
   Add node's neighbor's to correct "lttc->target_Nodes" arrays
*******************************************************************/
void AddNeighborsToTargetNodes( struct node_t* const node ) {
   enum nodeType_t type = node->type;
   int nbrIt;
   if( ( type == INFECT_INT       ) && 
       ( node->n_Nbrs[CANCER] > 0 ) ){
      for( nbrIt=0; nbrIt < node->n_Neighbors; nbrIt++){
	 struct node_t *nbr = node->neighbors[nbrIt];
	 if( nbr->type == CANCER ){
	    if( nbr->targets_Idx[INFECT_INT] != NO_IDX ) {
	       struct node_t **array =
		     lttc->target_Nodes[INFECT_INT][ nbr->n_Nbrs[INFECT_INT]-1 ];
	       unsigned long *last_Idx =
		  &lttc->n_Target_Nodes[INFECT_INT][ nbr->n_Nbrs[INFECT_INT]-1 ];

	       RemoveCell(&nbr->targets_Idx[ INFECT_INT ], array, last_Idx,
			  &array[*last_Idx-1]->targets_Idx[INFECT_INT] );
	    }

	    AddCell(nbr, lttc->target_Nodes[INFECT_INT][nbr->n_Nbrs[INFECT_INT]],
		  &lttc->n_Target_Nodes[INFECT_INT][nbr->n_Nbrs[INFECT_INT]],
		  &nbr->targets_Idx[INFECT_INT] );
	 }
      }
   }
   if( ( type == INFECT_EXT       ) && 
       ( node->n_Nbrs[CANCER] > 0 ) ){
      for( nbrIt=0; nbrIt < node->n_Neighbors; nbrIt++){
	 struct node_t *nbr = node->neighbors[nbrIt];
	 if( nbr->type == CANCER ){
	    if( nbr->targets_Idx[INFECT_EXT] != NO_IDX ) {
	       struct node_t **array =
		     lttc->target_Nodes[INFECT_EXT][ nbr->n_Nbrs[INFECT_EXT]-1 ];
	       unsigned long *last_Idx =
		  &lttc->n_Target_Nodes[INFECT_EXT][ nbr->n_Nbrs[INFECT_EXT]-1 ];

	       RemoveCell(&nbr->targets_Idx[ INFECT_EXT ], array, last_Idx,
			  &array[*last_Idx-1]->targets_Idx[INFECT_EXT] );
	    }

	    AddCell(nbr, lttc->target_Nodes[INFECT_EXT][nbr->n_Nbrs[INFECT_EXT]],
		  &lttc->n_Target_Nodes[INFECT_EXT][nbr->n_Nbrs[INFECT_EXT]],
		  &nbr->targets_Idx[INFECT_EXT] );
	 }
      }
   }
   if( type == NORMAL || type == CANCER ) {
      for( nbrIt=0; nbrIt < node->n_Neighbors; nbrIt++){
         struct node_t *nbr = node->neighbors[nbrIt];
         if( nbr->type == VACANT ){
            if( nbr->targets_Idx[type] != NO_IDX )
            {
               struct node_t **array =
                     lttc->target_Nodes[type][ nbr->n_Nbrs[type]-1 ];
               unsigned long *last_Idx =
                  &lttc->n_Target_Nodes[type][ nbr->n_Nbrs[type]-1 ];

               RemoveCell(&nbr->targets_Idx[ type ], array, last_Idx,
                          &array[*last_Idx-1]->targets_Idx[type] );
            }

            AddCell(nbr, lttc->target_Nodes[type][nbr->n_Nbrs[type]],
                  &lttc->n_Target_Nodes[type][nbr->n_Nbrs[type]],
                  &nbr->targets_Idx[type] );
         }
      }
   }

   return;
}
/* RemoveFromTargetNodes() *****************************************
   Removes node's index from all "lttc->target_Nodes" arrays
   Does not set node->targets_Idx values to NO_IDX
*******************************************************************/
void RemoveFromTargetNodes( struct node_t* const node ) {
   enum nodeType_t typeIt;

   for( typeIt = 0; typeIt < N_CELL_TYPES; typeIt++ ) {
      if( node->targets_Idx[typeIt] != NO_IDX ) {
         struct node_t **array = 
               lttc->target_Nodes[typeIt][ node->n_Nbrs[typeIt] ];
         unsigned long *last_Idx =  
            &lttc->n_Target_Nodes[typeIt][ node->n_Nbrs[typeIt] ];

         RemoveCell(&node->targets_Idx[ typeIt ], array, last_Idx,
                    &array[*last_Idx-1]->targets_Idx[typeIt] );
      }
   }
   return;
}
/* AddToTargetNodes() **********************************************
   Add node to correct "lttc->target_Nodes" arrays
   Updates node->targets_Idx
*******************************************************************/
void AddToTargetNodes( struct node_t* const node ) {
   
   enum nodeType_t type = node->type;

   if( type == CANCER )
   {
      if( node->n_Nbrs[INFECT_INT] > 0 )
         AddCell(node, lttc->target_Nodes[INFECT_INT][node->n_Nbrs[INFECT_INT]],
               &lttc->n_Target_Nodes[INFECT_INT][node->n_Nbrs[INFECT_INT]],
               &node->targets_Idx[INFECT_INT] );
      if( node->n_Nbrs[INFECT_EXT] > 0 )
         AddCell(node, lttc->target_Nodes[INFECT_EXT][node->n_Nbrs[INFECT_EXT]],
               &lttc->n_Target_Nodes[INFECT_EXT][node->n_Nbrs[INFECT_EXT]],
               &node->targets_Idx[INFECT_EXT] );
   }
   if( type == VACANT )
   {
      if( node->n_Nbrs[NORMAL] > 0 )
         AddCell(node, lttc->target_Nodes[NORMAL][node->n_Nbrs[NORMAL]],
               &lttc->n_Target_Nodes[NORMAL][node->n_Nbrs[NORMAL]],
               &node->targets_Idx[NORMAL] );
      if( node->n_Nbrs[CANCER] > 0 )
         AddCell(node, lttc->target_Nodes[CANCER][node->n_Nbrs[CANCER]],
               &lttc->n_Target_Nodes[CANCER][node->n_Nbrs[CANCER]],
               &node->targets_Idx[CANCER] );
   }
   return;
}
/* RemoveNeighborsFromTargetNodes() */
void RemoveNeighborsFromTargetNodes( struct node_t* const node ) {
   
   if (node->type == VACANT)
      return;


   int nbrIt;
   for (nbrIt = 0; nbrIt < node->n_Neighbors; nbrIt++) {
      struct node_t * nbr = node->neighbors[nbrIt];
      /* Adjust n_Nbrs of each neighbor */
      int n_Nbrs = nbr->n_Nbrs[node->type]--;
      /* If neighbor is a target, move to proper list */
      if( nbr->targets_Idx[node->type] != NO_IDX ) {
	 unsigned long *last_Node_Idx = 
	    &lttc->n_Target_Nodes[node->type][n_Nbrs];
	 struct node_t* last_Node =
	    lttc->target_Nodes[node->type][n_Nbrs][(*last_Node_Idx)-1];

	 RemoveCell( &nbr->targets_Idx[node->type],
		      lttc->target_Nodes[node->type][n_Nbrs],
		      last_Node_Idx,
		     &last_Node->targets_Idx[node->type] );
	 if( n_Nbrs-1 > 0 ) {
	    AddCell( nbr,
		     lttc->target_Nodes[node->type][n_Nbrs-1],
		    &lttc->n_Target_Nodes[node->type][n_Nbrs-1],
		    &nbr->targets_Idx[node->type] );
	 }
      }
   }

   return;
}

/** ChangeNodeType( node, type ) ************************************
  This routine will change "node"'s type to "type" (CANCER or INFECT)
********************************************************************/
void ChangeNodeType( struct node_t* node, enum nodeType_t type )
{
   int nbrIt;
   static int recurse = 0;

   if( type == INFECT_EXT || type == INFECT_INT ){
      if( 100.0*node->n_Nbrs[NORMAL]/node->n_Neighbors <= percent_Interior )
	 type = INFECT_INT;
      else
	 type = INFECT_EXT;
   }

   if( node->type == type )
      return;


   RemoveNeighborsFromTargetNodes( node );
   RemoveFromTargetNodes( node );
   /* Remove from Cell List */
   if( node->type != VACANT )
      RemoveCell( 
          &node->list_Idx,
           lttc->cells[node->type],
          &lttc->n_Cells[node->type],
          &lttc->cells[node->type][lttc->n_Cells[node->type]-1]->list_Idx );

   node->type = type;
   /* Check neighbor's virus type */
   if( recurse == 0 ){
      for (nbrIt = 0; nbrIt < node->n_Neighbors; nbrIt++) {
	 recurse = 1;
	 struct node_t*  nbr = node->neighbors[nbrIt];
	 enum nodeType_t nbr_Type = nbr->type;
	 if( nbr_Type == INFECT_EXT || nbr_Type == INFECT_INT ){
	    if( nbr->n_Nbrs[CANCER]/nbr->n_Neighbors > percent_Interior/100 )
	       ChangeNodeType( nbr, INFECT_INT );
	    else
	       ChangeNodeType( nbr, INFECT_EXT );
	 }
	 recurse = 0;
      }
   }

   /* Add to list of Cells */
   if( type != VACANT ) {
      AddCell(
            node,
            lttc->cells[type],
            &lttc->n_Cells[type],
            &node->list_Idx );
      for (nbrIt = 0; nbrIt < node->n_Neighbors; nbrIt++)
         node->neighbors[nbrIt]->n_Nbrs[node->type]++;
   }
   AddToTargetNodes( node );
   AddNeighborsToTargetNodes( node );

   return;
}/* End of ChangeNodeType */

/** SeedCancer( ) **********************************************************
  Places initial cancer cells in lattice.
  Replaces lttc->center and cancer_Count of it's nearest neighbors w/
          or until distance of cancer_Radius is reached (which ever is 
	  "greater") with cancer cells.
***************************************************************************/
int SeedCancer( struct node_t* center, 
                unsigned long cancer_Count, 
                unsigned cancer_Radius ) {
   struct node_t *node_u;
   unsigned long nodeIt;
   for( nodeIt=0; nodeIt < lttc->n_Nodes; nodeIt++ ) {
      lttc->nodes[nodeIt]->visited = NOTSEEN;
      lttc->nodes[nodeIt]->distance = INT_MAX;
   }

   node_u = center;

   int count = 0;
   struct node_t *node_v;
   struct Queue * const Q = NewQueue( 8 );
   node_u->distance = 0;
   node_u->visited = SEEN;
   Enqueue( Q, node_u );
   while( (node_u = Dequeue(Q) ) ){
      if( ( node_u->distance > cancer_Radius && node_u->distance < INT_MAX )
            && count >= cancer_Count ){
         continue;
      }
      count++;
      ChangeNodeType( node_u, CANCER );

      node_u->visited = VISITED;
      int nbrIt;
      for( nbrIt = 0; nbrIt < node_u->n_Neighbors; nbrIt++ ) {
         node_v = node_u->neighbors[nbrIt];
         if( node_v->visited == NOTSEEN ) {
            node_v->visited = SEEN;
            node_v->distance = node_u->distance + 1;
            if( Enqueue( Q, node_v ) == 0 )
               printf(" Problem growing Queue\n" );
         }
      }
   }
   DeleteQueue( Q );
   return 0;
}/* End of SeedCancer() */


/** InjectVirus( ) *********************************************
  Changes roughly percent_Infected CANCER cells to INFECT cells
***************************************************************/
int InjectVirus( ) {
   unsigned long n_Infect = lttc->n_Cells[CANCER] * percent_Infected/100;
   struct node_t *center;
   Point* c_point;
   switch( infect_Type ){
   case RANDOM:
      while( lttc->n_Cells[INFECT_INT]+lttc->n_Cells[INFECT_EXT]< n_Infect ){
         unsigned long random;
         do {
            random = lttc->n_Cells[CANCER] * genrand();
         } while( random >= lttc->n_Cells[CANCER] );
         ChangeNodeType( lttc->cells[CANCER][random], INFECT_EXT );
      }
      break;
   case CENTER: {
      c_point = CenterOfMass( lttc->cells[CANCER], 
                              lttc->n_Cells[CANCER] );
      center = FindNodeNear( c_point, lttc->cells[CANCER], 
                                      lttc->n_Cells[CANCER] );
      VirusBallAroundNode( center, n_Infect );
      DeletePoint( c_point );
      break; }
   case MULTINODE: {
      c_point = CenterOfMass( lttc->cells[CANCER], 
                              lttc->n_Cells[CANCER] );
      center = FindNodeNear( c_point, lttc->cells[CANCER], 
                                      lttc->n_Cells[CANCER] );

      int It;
      for( It = 0; It < 3; It++ ) {
         double random = genrand() * 2 * M_PI;
         struct node_t* node = FindEdge( lttc->cells[CANCER], 
                                         lttc->n_Cells[CANCER], 
                                         center, random );
         VirusBallAroundNode( node, n_Infect/3 );
      }
      DeletePoint( c_point );
      break; }
   case PERIMETER: {
      c_point = CenterOfMass( lttc->cells[CANCER], lttc->n_Cells[CANCER]  );
      double* distance[2];
      distance[0] = malloc( lttc->n_Cells[CANCER] * sizeof(double));
      distance[1] = malloc( lttc->n_Cells[CANCER] * sizeof(double));
      unsigned long nodeIt;
      for( nodeIt=0; nodeIt < lttc->n_Cells[CANCER]; nodeIt++){
         struct node_t *cell = lttc->cells[CANCER][nodeIt];
         distance[0][nodeIt] = cell->Idx;
         distance[1][nodeIt] = Point2PointDistance( c_point, cell->pos);
      }
      MergeSort( distance, lttc->n_Cells[CANCER]);
      for( nodeIt=0; nodeIt < n_Infect; nodeIt++)
         ChangeNodeType( lttc->nodes[ (int)distance[0][nodeIt] ], INFECT_EXT );
      free( distance[0] );
      free( distance[1] );
      DeletePoint( c_point );
      break; }
   default:
      return -1;
   }
   return 0;

}/* End of InjectVirus() */

/* MergeSort() ********************************************************
   Performs a merge sort on the N data pairs in list.
   list[1] contains the values to sort on
   list[0][i] contains a value that remains associated with list[1][i]
**********************************************************************/
void MergeSort( double** list, unsigned long N ){
   double* temp[2];
   temp[0] = malloc( N * sizeof( double ));
   temp[1] = malloc( N * sizeof( double ));
   MergeSplit( list, 0, N, temp );
   free( temp[0] );
   free( temp[1] );
   return;
}

/* MergeSplit() *******************************************************
   Helper function for MergeSort().  
   Should not be called outside MergeSort()
**********************************************************************/
void MergeSplit( double** list,
      unsigned long b_idx,
      unsigned long e_idx,
                 double** temp ){
   if( e_idx - 1 <= b_idx )
      return; //Length 1 list is sorted

   unsigned long m_idx = ( e_idx + b_idx ) / 2;
   MergeSplit( list, b_idx, m_idx, temp );
   MergeSplit( list, m_idx, e_idx, temp );
   Merge( list, b_idx, m_idx, e_idx, temp);
   unsigned long It;
   for( It = b_idx; It < e_idx; It++ ) {
      list[0][It] = temp[0][It];
      list[1][It] = temp[1][It];
   }
}

/* Merge() *******************************************************
   Helper function for MergeSort().  
   Should not be called outside MergeSort()
**********************************************************************/
void Merge( double** list, unsigned long b_idx,
                          unsigned long m_idx,
                          unsigned long e_idx,
            double** temp ){
   unsigned long It;
   unsigned long LH = b_idx;
   unsigned long RH = m_idx;

   for(It = b_idx; It < e_idx; It++){
      if( (LH < m_idx && RH >= e_idx) ||  
          (LH < m_idx && list[1][LH] > list[1][RH])){
         temp[0][It] = list[0][LH];
         temp[1][It] = list[1][LH];
         LH++;
      } else {
         temp[0][It] = list[0][RH];
         temp[1][It] = list[1][RH];
         RH++;
      }
   }
   return;
}


/* FindEdge() *********************************************************
   construct line L with direction "dir" that intersects start

   set node to "start"
   while( node's neighbors are all cancer )
      mark node
      find neighbor "nbr" of node closest to L
      set node to "nbr"

   returns address of node along L that is at the "edge" of the tumor
**********************************************************************/
struct node_t* FindEdge( struct node_t** nodes, unsigned long N, 
                         struct node_t*  start, double dir ){
   unsigned long nodeIt = 0;
   int n_Dims = nodes[0]->pos->N;
   for( nodeIt=0; nodeIt < N; nodeIt++ )
      nodes[nodeIt]->visited = NOTSEEN;

   struct node_t* node = start;
   node->visited = VISITED;
   Vector *v = NULL;
   switch( n_Dims ) {
   case 2:
      v = NewVector(n_Dims, cos( dir ), sin( dir ));
      break;
   case 3:
      v = NewVector(n_Dims, cos( dir ), sin( dir ), start->pos->data[2]);
      break;
   default:
      break;
   }
   Line *l = NewLine(start->pos, v );

   while( node->n_Nbrs[CANCER] == node->n_Neighbors ){
      int nbrIt = 0;
      double s_dist, dist;
      struct node_t* temp = NULL;
      for( nbrIt = 0; nbrIt < node->n_Neighbors && temp == NULL; nbrIt++) {
         if( node->neighbors[nbrIt]->visited == NOTSEEN ){
            s_dist = Point2LineDistance( node->neighbors[nbrIt]->pos, l );
            temp = node->neighbors[nbrIt];
         }
      }
      for( ; nbrIt < node->n_Neighbors; nbrIt++) {
         if( node->neighbors[nbrIt]->visited == VISITED )
            continue;

         dist = Point2LineDistance( node->neighbors[nbrIt]->pos, l );
         if( dist < s_dist ){
            temp = node->neighbors[nbrIt];
            s_dist = dist;
         }
      }
      node = temp;
      node->visited = VISITED;
   }
   DeleteVector( v );
   DeleteLine( l );

   return node;
}
/** VirusBallAroundNode() ********************************************
  Creates a densly packed ball of n_Infect infected cells centered
  around node.
*********************************************************************/
int VirusBallAroundNode( struct node_t * node, unsigned long n_Infect )
{
   unsigned long nodeIt;
   int count = 0;
   Point *point = node->pos;

   for( nodeIt = 0; nodeIt < lttc->n_Cells[CANCER]; nodeIt++ ) {
      lttc->cells[CANCER][nodeIt]->visited = NOTSEEN;
      lttc->cells[CANCER][nodeIt]->distance = UINT_MAX;
   }
   while( lttc->n_Cells[CANCER] > 0 && count < n_Infect )
   {
      if( node->type != INFECT_INT || node->type != INFECT_EXT )
         count += BFS(node, CANCER, n_Infect-count, INFECT_EXT );
      node = FindNodeNear(point, lttc->cells[CANCER], 
	                        lttc->n_Cells[CANCER]);
   }

   return count;
}
/* DetermineEventRecipient() ****************************************
  selects a node for action.
   if killing type T, randomly select node of type T.

   if growing type T, randomly select from targetof(T) nodes weighted
        by the number of neighboring source nodes.
*********************************************************************/
struct node_t * DetermineEventRecipient(const struct event_t e_t ) {
   struct node_t* recp = NULL;
   double random;
   enum nodeType_t type = e_t.node_Type;

   do {
      random = genrand();
   } while( random == 1 || random == 0 );

   switch (e_t.cell_Action) {
   case KILL:
      random *= lttc->n_Cells[type];
      recp = lttc->cells[type][ (unsigned long int) random ];
      break;
   case GROW:
      {
         //Select a cell to replicate to replicates onto.
         int n_nbrIt;
         unsigned long int n_Links=0;
         for(n_nbrIt = 1; n_nbrIt <= lttc->n_Nbrs_Max; n_nbrIt++)
            n_Links += lttc->n_Target_Nodes[ type ][n_nbrIt] * n_nbrIt;
         random *= n_Links;
         n_nbrIt = 1;
         while( lttc->n_Target_Nodes[ type ][n_nbrIt]*n_nbrIt <= random ) {
            random  -= lttc->n_Target_Nodes[ type ][n_nbrIt] * n_nbrIt;
            n_nbrIt++;
         }
         recp = lttc->target_Nodes[type][n_nbrIt][(unsigned long int)(random/n_nbrIt)];
      }
      break;
   default:
      fprintf(stderr, "Action was not KILL or GROW\n");
      return NULL;
   }
   if( recp == NULL ) {
      fprintf(stderr, "\nSelected cell is NULL.\n");
      fprintf(stderr, "lttc->time = %f\n", lttc->time );
      return NULL;
   }
   return recp;
}/* End of DetermineEventRecipient() */ 


/**************************************************************************
  Print Percentage of Node Types
**************************************************************************/
int OutputStats() {
   /* This is coded so that the output when only Normal cells exist is the
      same as the output created with Chetan's code */
   printf("  time: %6.2f   Density: %14.6f",
         lttc->time, (float)lttc->n_Cells[NORMAL]/lttc->n_Nodes);
   if( lttc->n_Cells[CANCER] > 0 ){
      printf("   Cancer Density: %14.6f",
                     (float)lttc->n_Cells[CANCER]/lttc->n_Nodes);
   }
   long unsigned int n_Virus = lttc->n_Cells[INFECT_INT] +
                               lttc->n_Cells[INFECT_EXT];
   if( n_Virus > 0 ) {
      if( lttc->n_Cells[CANCER] == 0 )
         printf("                                 " );
      printf("   Virus Density: %14.6f",
                     (float)n_Virus/lttc->n_Nodes);
   }
   printf("\n");


   return 0;
}/* End of OutputStats() */

/**************************************************************************
  Print Percentage of Cancer Nodes
**************************************************************************/
int OutputSummary() {
   printf(" %6.2f, ", (float)lttc->n_Cells[CANCER]/lttc->n_Nodes);

   return 0;
}/* End of OutputSummary() */

/** OutputGrid( lattice ) *****************
* Print 2d lattice to stdout as ASCII art
******************************************/
void OutputGrid() {
   const struct node_t* cell = lttc->nodes[0];
   const float width = sqrt( lttc->n_Nodes );
   int count = 0;

   if( width != (unsigned long) width ) {
      fprintf( stderr, "Lattice isn't square. Refusing to print grid.\n" );
      return;
   }

   unsigned long nodeIt;
   for( nodeIt=0; nodeIt < lttc->n_Nodes; nodeIt++ ) {
      cell = lttc->nodes[nodeIt];
      switch (cell->type) {
      case NORMAL:
         printf("\x1b[32m .\x1b[0m ");
         break;
      case CANCER:
         printf("\x1b[31m c\x1b[0m ");
         break;
      case INFECT_INT:
         printf(" i ");
         break;
      case INFECT_EXT:
         printf(" I ");
         break;
      default:
         printf(" _ ");
         break;
      }
      count++;
      if( count%(int)width%10 == 0 )
         printf(" | ");
      if( count % (int)width == 0 ){
         printf("\n");
	 if(count % ((int)width*10) == 0 ){
	    int i;
	    for( i=0; i<width + width/10; i++ )
	       printf(" = ");
	    putchar('\n');
	 }
      }
   }
   return;
}/* End of OutputGrid() */ 

/** OutputLatticeState( lattice ) ******************************************
* Print current lattice state to file for analysis
***************************************************************************/
void OutputLatticeState(char* filename) {
   FILE* fh;
   const struct node_t* cell = lttc->nodes[0];

   fh = fopen( filename, "w" );
   if( fh == NULL ) {
        fprintf(stderr, "Failed to open %s for output. %s\n",
              filename, strerror(errno));
        exit(-1);
   }

   unsigned long nodeIt;
   for( nodeIt=0; nodeIt < lttc->n_Nodes; nodeIt++ ) {
      cell = lttc->nodes[nodeIt];

      int dimIt;
      for( dimIt=0; dimIt < lttc->n_Dimensions; dimIt++)
         fprintf(fh, "%f ", cell->pos->data[dimIt]);

      switch (cell->type) {
      case NORMAL:
         fprintf(fh, " 1 \n");
         break;
      case CANCER:
         fprintf(fh, " 2 \n");
         break;
      case INFECT_INT:
      case INFECT_EXT:
         fprintf(fh, " 3 \n");
         break;
      default: //NORMAL
         fprintf(fh, " 0 \n");
      }
   }
   fclose( fh );

   return;
}/* End of OutputLatticeState() */


/** DetermineEventType( lattice ) ***********************************
* Calculate the probabilities of event types and randomly select one
* 6 bins:  
*   3 KILL are sized rate_i * n_Cells_i;
*   3 GROW are sized rate_i * (n_links from i to TargetOf(i))/n_Nbrs_Avg
* Adjust lattice time
***************************************************************************/
struct event_t DetermineEventType( ){

   double ev_Bin[ N_CELL_TYPES ][ N_ACTIONS ] = {{0,0},{0,0},{0,0},{0,0}};
   double random;
   struct event_t ev;
   int actionIt;
   double summation = 0;

   int typeIt, n_NbrsIt;
   for( typeIt = NORMAL; typeIt < N_CELL_TYPES; typeIt++ ) {
      int n_Links=0;
      ev_Bin[typeIt][KILL] = rates[typeIt][KILL] * lttc->n_Cells[typeIt];

      for( n_NbrsIt=0; n_NbrsIt <= lttc->n_Nbrs_Max; n_NbrsIt++ ) {
         n_Links += lttc->n_Target_Nodes[ typeIt ][ n_NbrsIt ] * n_NbrsIt;
      }
      ev_Bin[ typeIt ][ GROW ] =
               rates[ typeIt ][GROW] * (double) n_Links / lttc->n_Nbrs_Avg;
   }

   for( typeIt = NORMAL; typeIt < N_CELL_TYPES; typeIt++ ) {
      for(actionIt=KILL; actionIt < N_ACTIONS; actionIt++) {
         summation += ev_Bin[ typeIt ][actionIt];
         ev_Bin[ typeIt ][ actionIt ] = summation;
      }
   }

   do {
      random = genrand();
   } while (random == 0 || random == 1 );

   lttc->time += -log( random ) / ev_Bin[INFECT_EXT][KILL];
   lttc->count++;

   do {
      random = genrand();
   } while (random == 0 || random == 1 );
   random = random * ev_Bin[N_CELL_TYPES-1][N_ACTIONS-1];

   for( ev.node_Type=NORMAL; ev.node_Type < N_CELL_TYPES; ev.node_Type++ )
   {
      for(ev.cell_Action=KILL; ev.cell_Action < N_ACTIONS; ev.cell_Action++)
      {
         if( random < ev_Bin[ev.node_Type][ev.cell_Action] ) {
            return ev;
         }
      }
   }
   fprintf(stderr, "\n\n\n\nNo event selected.\n");
   fprintf(stderr, "Number of NORMAL cells: %lu\n", lttc->n_Cells[NORMAL]);
   fprintf(stderr, "Number of CANCER cells: %lu\n", lttc->n_Cells[CANCER]);
   fprintf(stderr, "Number of INFECT_INT cells: %lu\n", lttc->n_Cells[INFECT_INT]);
   fprintf(stderr, "Number of INFECT_EXT cells: %lu\n", lttc->n_Cells[INFECT_EXT]);
   ev.node_Type=VACANT;

   return( ev );
}/* End of DetermineEventType() */



/** RemoveCell( cell_Idx, array, n_cells, index ) ************************
* Remove cell at "del_idx" from "array" by moving array[*n_cells] to it's
* location. Update "n_cells" in array, "del_idx" and index of former last
* node "last_idx".
*
* Doesn't perform sanity checks on arguments.
* If del_idx, n_cells, and last_idx aren't valid for array bad things
*     will happen.
************************************************************************/
void RemoveCell (unsigned long* const del_Idx,
                 struct node_t** array,
                 unsigned long * const n_cells,
                 unsigned long * const last_Idx) {
   const int idx_a = *del_Idx;
   const int idx_b = *last_Idx;

   if( *del_Idx < *last_Idx ) {
      array[ idx_a ] = array[ idx_b ];
      *last_Idx = idx_a;
   }
   (*n_cells)--;
   array[ *n_cells ] = NULL;
   *del_Idx = NO_IDX;

   return;
}/* End of RemoveCell() */

/** AddCell( cell, array, n_Array, index ) *********************************
* Add cell to array with n_Array nodes.
* Set cell's index into array.
***************************************************************************/
void AddCell( struct node_t* const cell,
              struct node_t** array,
              unsigned long*const  n_Array,
              unsigned long*const  cell_Idx ){
   array[*n_Array] = cell;
   *cell_Idx = *n_Array;
   (*n_Array)++;
   return;
}/* End of AddCell() */


/** TargetOf( nodeType )****************************************************
* Function takes a cell type as an argument and returns its target.
* If type is invalid or has no target, program quits
***************************************************************************/
enum nodeType_t TargetOf( const enum nodeType_t type ) {
      switch( type ){
      case INFECT_INT:
      case INFECT_EXT:
         return CANCER;
      case NORMAL:  //Treat the same as CANCER.
      case CANCER:
         return VACANT;
      case VACANT:
         fprintf(stderr, "Asked for target type of VACANT node\n");
         exit(-1);
      default:
         fprintf(stderr, "Asked for target type of undefined node type\n");
         exit(-1);
      }

}/* End of TargetOf() */



void DeleteLattice( struct lattice_t* lattice ){
   unsigned long int nodeIt = 0;
   enum nodeType_t typeIt = 0;
   int nbrIt=0;
   for( typeIt = 0; typeIt < N_CELL_TYPES; typeIt++ ) {
      for( nbrIt=0; nbrIt <= lattice->n_Nbrs_Max; nbrIt++ ) {
	 free( lattice->target_Nodes[typeIt][nbrIt] );
	 lattice->target_Nodes[typeIt][nbrIt] = NULL ;
      }
   }
   for( typeIt = 0; typeIt < N_CELL_TYPES; typeIt++ ) {
      free( lattice->n_Target_Nodes[typeIt] );
      lattice->n_Target_Nodes[typeIt] = NULL;
      free( lattice->target_Nodes[typeIt] );
      lattice->target_Nodes[typeIt] = NULL;
      free( lattice->cells[typeIt] );
      lattice->cells[typeIt] = NULL;

   }
   for( nodeIt = 0; nodeIt < lattice->n_Nodes; nodeIt++ ) {
      DeletePoint( lattice->nodes[nodeIt]->pos );
      lattice->nodes[nodeIt]->pos = NULL ;
      free( lattice->nodes[nodeIt]->neighbors );
      lattice->nodes[nodeIt]->neighbors = NULL ;
      free( lattice->nodes[nodeIt] );
      lattice->nodes[nodeIt] = NULL ;
   }

   free( lattice->nodes );
   lattice->nodes = NULL;
   free( lattice );
   lattice = NULL;
   return;
}

void OpenFiles(){
      fpp=fopen("vpts.dat", "a+");
      fpl=fopen("vlpt.dat", "a+");
      unsigned long int n_Tumor = lttc->n_Cells[CANCER]     +
	                          lttc->n_Cells[INFECT_INT] +
                        	  lttc->n_Cells[INFECT_EXT] ;
      fprintf(fpp," %d ",rng_Seed);
      fprintf(fpl," %d ",rng_Seed);
      fprintf(fpp," %lu ", n_Tumor);
      fprintf(fpl," %lu ", n_Tumor);
}

void CloseFiles(){
   int time_Check     = (int) floor(lttc->time);
   int time_Check_fpl = (int) floor(lttc->time);
   unsigned long int n_Tumor = lttc->n_Cells[CANCER]     +
			       lttc->n_Cells[INFECT_INT] +
			       lttc->n_Cells[INFECT_EXT] ;
   while( time_Check_fpl % 5 != 0 )
      time_Check_fpl--;

   while( time_Check < 126) {
      time_Check = time_Check+1;
      fprintf(fpp," %lu ", n_Tumor);
   }
   fprintf(fpp,"\n");
   fclose(fpp);


   if( time_Check_fpl < t_Max_Simulation ) {
      time_Check_fpl+=5;
      fprintf(fpl," %lu ", n_Tumor);
   }
   fprintf(fpl,"\n");
   fclose(fpl);
}

void ParameterFile(float time_Extinct){
   FILE *fp;
   fp=fopen("vout.dat", "a+");
   fprintf(fp," %d ",rng_Seed);
   fprintf(fp," %f ",rates[NORMAL][GROW]);
   fprintf(fp," %f ",rates[NORMAL][KILL]);
   fprintf(fp," %f ",rates[CANCER][GROW]);
   fprintf(fp," %f ",rates[CANCER][KILL]);
   fprintf(fp," %f ",rates[INFECT_INT][GROW]);
   fprintf(fp," %f ",rates[INFECT_INT][KILL]);
   fprintf(fp," %f ",rates[INFECT_EXT][GROW]);
   fprintf(fp," %f ",rates[INFECT_EXT][KILL]);
   fprintf(fp," %f ",percent_Infected);
   if( time_Extinct > 0 ) { 
     fprintf(fp," %f ",time_Extinct); 
   }
   else {
     fprintf(fp," %f ",lttc->time);
   }
   fprintf(fp," %d ",(int)lttc->n_Cells[NORMAL]);
   fprintf(fp," %d ",(int)lttc->n_Cells[CANCER]);
   fprintf(fp," %d ",(int)lttc->n_Cells[INFECT_INT]);
   fprintf(fp," %d ",(int)lttc->n_Cells[INFECT_EXT]);
   unsigned long int n_Virus = lttc->n_Cells[INFECT_INT] +
			       lttc->n_Cells[INFECT_EXT] ;
   
   if( lttc->n_Cells[CANCER] == 0 ) { fprintf(fp," 1\n"); }
   if( lttc->n_Cells[CANCER]  > 0 && 
       n_Virus               == 0 ) { fprintf(fp," 2\n"); }
   if( lttc->n_Cells[CANCER]  > 0 && 
       n_Virus                > 0 ) { fprintf(fp," 3\n"); }
   fclose(fp);
   return;
}
