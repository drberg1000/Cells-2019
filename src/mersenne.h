#ifndef DB_MERSENNE_HEADER_FILE
#define DB_MERSENNE_HEADER_FILE
/*   genrand() generates one pseudorandom real number (double) */
/* which is uniformly distributed on [0,1]-interval, for each  */
/* call. sgenrand(seed) set initial values to the working area */
/* of 624 words. Before genrand(), sgenrand(seed) must be      */
/* called once. (seed is any 32-bit integer except for 0).     */

/* initializing the array with a NONZERO seed */
void sgenrand(unsigned long seed);

/* generating reals from [0,1] */
double genrand();

#endif
