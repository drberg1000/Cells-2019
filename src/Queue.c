#include  "Queue.h"

/* NewQueue() *************************************************************/
/*  Size indicates initial size of array used to store Queue              */
/*  Queue will double in size as needed to make room for new objects.     */
/**************************************************************************/
struct Queue* const NewQueue( int size ) {
   struct Queue* const Q = malloc( sizeof(struct Queue) );
   if( Q == NULL ){
      return NULL;
   }
   Q->q = calloc( size, sizeof(void*) );
   if( Q->q == NULL ){
      return NULL;
   }
   Q->tail = 0;
   Q->head = 0;
   Q->length = size;

   return Q;
}
/* DeleteQueue() **********************************************************/
/* Frees memory associated with Q                                         */
/**************************************************************************/
void DeleteQueue( struct Queue* const Q ) {
   free( Q->q );
   free( Q );
   return;
}

/* Enqueue() **************************************************************/
/*  adds object x to Queue Q                                              */
/*  Q will double in size as needed to make room for new objects          */
/**************************************************************************/
int Enqueue( struct Queue* const Q, void* x ) {
   if( Q->head == (Q->tail + 1) ) {
      if( GrowQueue( Q ) == 0 ) 
         return 0;
   }
   Q->q[Q->tail] = x;
   if(Q->tail == Q->length-1)
      Q->tail = 0;
   else
      Q->tail += 1;
   return 1;
}

/* Enqueue() **************************************************************/
/*  adds object x to Queue Q                                              */
/*  Q will not shrink regardless of how many elements are removed         */
/**************************************************************************/
void* Dequeue( struct Queue* const Q ) {
   void* x = Q->q[Q->head];
   Q->q[Q->head] = NULL;

   if( Q->head == Q->tail )
      return NULL;
   if( Q->head == Q->length-1 )
      Q->head = 0;
   else
      Q->head += 1;

   return x;
}


/* GrowQueue() ************************************************************/
/*   Doubles the size of Queue and moves elements as needed to maintain   */
/*   the appropriate structure                                            */
/**************************************************************************/
int GrowQueue( struct Queue* const Q ) {
   void** temp;
   temp = realloc(Q->q, Q->length*2 * sizeof(void *));
   if( temp == NULL ){
      return 0;
   }
   Q->q = temp;
   if( Q->tail < Q->head ){
      int It;
      for(It = 0; It < Q->tail; It++)
         Q->q[Q->length+It] = Q->q[It];

      Q->tail += Q->length;
   }
   Q->length *= 2;

   return 1;
}

/* Display() **************************************************************/
/* Outputs contents of Q and indicates where the head and tail are.       */
/*    Helpful for debugging or understanding how the Queue works.         */
/*                                                                        */
/* Assumes that Q's elements are char*.                                   */
/* Add argument providing a function to print actual type.                */
/**************************************************************************/
void Display( struct Queue* const Q ) {
   int It = 0;
   for(It = 0; It < Q->length; It++) {
      if( (It <  Q->tail &&     It >= Q->head) ||
          (It <  Q->tail && Q->tail < Q->head) ||
          (It >= Q->head && Q->tail < Q->head) ) {
         char c = *(char*)Q->q[It];
         printf("%c", c );
      }
      else printf("-");
   }
   puts("");

   for(It = 0; It < Q->length; It++ ){
      if(It == Q->head)
         putchar('h');
      else if(It == Q->tail)
         putchar('t');
      else if( (It <  Q->tail &&     It >= Q->head) ||
          (It <  Q->tail && Q->tail < Q->head) ||
          (It >= Q->head && Q->tail < Q->head) )
         putchar('-');
      else
         putchar(' ');
   }
   puts("");
}
