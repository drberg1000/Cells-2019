#include "cartesian.h"

/**************************************/
/**************************************/
/*                                    */
/*          Public Routines           */
/*                                    */
/**************************************/
/**************************************/

/***********************************************************/
/*                         POINTS                          */
/***********************************************************/
/* Creates a new Point;  */
/* Arguments 2->N+1 are expected to be cartesian coordinates */
Point* NewPoint( short unsigned N, ... ) {
   va_list arguments;                     
   va_start ( arguments, N );           
   Point *P = NewPoint2(N, arguments);
   va_end ( arguments );

   return P;                      
}
void DeletePoint( Point* P ) {
   if( P != NULL )
      free( P->data );
   free( P );
   P = NULL;
   return;
}


/***********************************************************/
/*                         VECTORS                         */
/***********************************************************/
/* Creates a new Vector;  Be sure to pass floats/doubles */
/* Arguments 2->N+1 are expected to be cartesian coordinates */
Vector* NewVector( short unsigned N, ... ) {
   va_list arguments;
   va_start ( arguments, N );           
   Vector *V = (Vector *) NewPoint2(N, arguments);
   va_end ( arguments );

   return  V;
}
void DeleteVector( Vector* V ) {
   DeletePoint( (Point*) V );
   return;
}

/***********************************************************/
/*                         LINES                           */
/***********************************************************/

/* Create a new Line with slope S that intersects point P */
Line*  NewLine( Point* P, Vector* S ){
   Line *L = malloc( sizeof( *L ));
   if( L == NULL || P->N != S->N ) {
      free( L );
      L = NULL;
      return NULL;
   }
   
   L->N = P->N;
   L->P = P;
   L->D = S;
   return L;
}
void DeleteLine( Line* L ){
   free(L);
   L = NULL;
}

/***********************************************************/
/*                    UTILITY FUNCTIONS                    */
/***********************************************************/
double Point2PointDistance( Point* P1, Point* P2 ){
   double distance = 0;
   unsigned short It;
   for(It = 0; It < P1->N; It++){
      double temp = (P1->data[It]-P2->data[It]);
      distance +=  temp * temp;
   }

   return sqrt(distance);
}

double Point2LineDistance( Point *P, Line *L ){
   double slope=0;
   double temp1=0;
   double temp2=0;
   double distance = 0;
   double vector[P->N];
   unsigned short It;

   if( P->N != L->N )
      return 0;
   
   /* Calculate Slope */
   for(It = 0; It < P->N; It++){
      temp1 += L->D->data[It] * (L->P->data[It] - P->data[It]);
      temp2 += L->D->data[It] * L->D->data[It];
   }
   slope = -1 * temp1 / temp2;

   for(It = 0; It < P->N; It++)
      vector[It] = L->P->data[It] + L->D->data[It] * slope - P->data[It];


   for(It = 0; It < P->N; It++)
      distance += vector[It] * vector[It];

   return sqrt(distance);
}

/**************************************/
/**************************************/
/*                                    */
/*         Private Routines           */
/*                                    */
/**************************************/
/**************************************/

/* Creates a new point; */
Point* NewPoint2( short unsigned N, va_list arguments){
   /* Create Space */
   Point *P = malloc( sizeof( *P ) );
   if( P == NULL )
      return P;
   P->data = malloc( N * sizeof( *P->data ));
   if( P->data == NULL ) {
      free( P );
      P = NULL;
      return P;
   }

   /* Initialize */
   P->N = N;
   int It;
   for( It = 0; It < N; It++ )        {
      P->data[It] =  va_arg( arguments, double );
   }

   return P;
}
