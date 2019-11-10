#ifndef DB_QUEUE_HEADER_FILE
#define DB_QUEUE_HEADER_FILE

#include <stdio.h>
#include <stdlib.h>

/******************/
/* Public Methods */
/******************/

/* Returns the address of a new Queue with space for "size" void pointers */
struct Queue* const   NewQueue( int size );

/* Adds a new void pointer (x) to Queue Q.  If Q is full, it's size doubles.
   Returns 0 if Q tried to grow but failed, otherwise returns 1           */
int              Enqueue( struct Queue* const Q, void* x  );

/* Returns and removes a void pointer from Queue Q. */
void*            Dequeue( struct Queue* const Q );

/* Free's Queue Q and associated memory */
void           DeleteQueue( struct Queue* const Q );

/* Displays value at each address in Q as characters separated by spaces 
   on one line and markers indicating head and tail on the second.  Useful
   primarily as a way to visualize how the Queue works.  */
/* Future version will allow specification of an argument to print 
   non-character types */
void             Display( struct Queue* const Q );



/*******************/
/* Private Methods */
/*******************/

struct Queue {
   int tail;
   int head;
   int length;
   void** q;
};

/* Routine to double the size of a full Queue and 
   move elements to the new location.               */
int            GrowQueue( struct Queue* const Q );

#endif
