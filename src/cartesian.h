#ifndef DB_CARTESIAN_H
#define DB_CARTESIAN_H

#include <math.h>
#include <stdarg.h>
#include <stdlib.h>

typedef struct Point Point;
typedef struct Vector Vector;
typedef struct Line Line;

/* Constructors allow creation of 1D & 2D objects in N-dimensional space. 
   N is an unsigned short integer which is at least 2^16-1 */
Point* NewPoint( unsigned short N, ... /*doubles*/ ) ;
Vector* NewVector( unsigned short N, ... /*doubles*/ ) ;
Line* NewLine( Point* P, Vector* Dir ); /* P is on line, 
    				           Dir is vector parallel to line */

/* If P & L have different dimensions Point2LineDistance returns 0;
   Otherwise it returns the shortest distance from P to L */
double Point2LineDistance( Point *P, Line *L );
double Point2PointDistance( Point* P1, Point* P2);

/* Destructors */
void DeleteVector( Vector* V );
void DeletePoint( Point* P );
void DeleteLine( Line* L );

/* Private Members */
Point* NewPoint2( unsigned short N, va_list arguments);

struct Point {
   double  * data;
   unsigned short N;
};
struct Vector {
   double  * data;
   unsigned short N;
};
struct Line {
  Point* P;  //Point on line
  Vector* D; //Direction of Line
  unsigned short N;
};

#endif
