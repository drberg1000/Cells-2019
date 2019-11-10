#include <stdlib.h> 
#include <getopt.h>
#include <stdio.h>  
#include "Cells.h"

/* Arguments */
static struct option const long_options[] = {
   /* The third value in this structure is not the default setting.
      To determine defaults, please see Parameters in Cells.h */
  /* Long Name      | Argument type    | ? | Short Name */
   {"Equalize",        no_argument,       0, 'E' },
   {"NormalGrowth",    required_argument, 0, 'N' },
   {"NormalDeath",     required_argument, 0, 'n' },
   {"CancerGrowth",    required_argument, 0, 'C' },
   {"CancerDeath",     required_argument, 0, 'c' },
   {"ViralIntGrowth",  required_argument, 0, 'I' },
   {"ViralIntDeath",   required_argument, 0, 'i' },
   {"ViralExtGrowth",  required_argument, 0, 'J' },
   {"ViralExtDeath",   required_argument, 0, 'j' },
   {"PercentInfected", required_argument, 0, 'V' },
   {"InfectionType",   required_argument, 0, 'T' },
   {"PercentInterior", required_argument, 0, 'P' },
   {"CancerRadius",    required_argument, 0, 'r' },
   {"CancerCount",     required_argument, 0, 'q' },
   {"RandomSeed",      required_argument, 0, 'R' },
   {"TimeMax",         required_argument, 0, 'm' },
   {"TimeInfect",      required_argument, 0, 'v' },

   {"OutputInterval",  required_argument, 0, 'o' },
   {"Grid",            no_argument,       0, 'g' },
   {"Stats",           no_argument,       0, 's' },
   {"Summary",         no_argument,       0, 'S' },
   {"TumorSizes",      no_argument,       0, 't' },
   {"Watch",           no_argument,       0, 'W' },
   {"File",            optional_argument, 0, 'f' },
   {"Files",           no_argument,       0, 'F' },
   {"Help",            no_argument,       0, 'H' },
   {"help",            no_argument,       0, 'h' },
   {0,                0,                 0,  0  }
   /* Add new short options to the while condition */
   /* Update Usage() to reflect changes */
   /* Update Switch() statement to reflect changes */
};

void NormalEquilibrium();
void Usage();

/** main() ************************************************************
*   Processes commandline arguments
*   Determines when Normal equilibrium is reached
*   Calls Cancer & Virus Seeding at appropriate times
***********************************************************************/
int main (int argc, char* argv[] ) {
   char *data;
   char *network;
   char *coords;
   int optc;
   int long_index = 0;

   int cancer_Radius = 1;
   int cancer_Count = 1;

   /* Process Arguments */
   while( (optc = getopt_long( argc, argv, 
                               "hHEFN:n:C:c:I:i:J:j:m:v:gsfo:P:V:r:SR:tWq:T:", 
			       long_options, &long_index )) != EOF ) {
      switch( optc ) {
         case 'H':
         case 'h':
            Usage();
            break;
         case 'r':
            cancer_Radius = argtof( optarg );
            break;
         case 'q':
            cancer_Count = argtof( optarg );
            break;
	 case 'E':
	    do_Equalization = 1;
	    break;
         case 'N':
            rates[NORMAL][GROW] = argtof( optarg );
            break;
         case 'n':
            rates[NORMAL][KILL] = argtof( optarg );
            break;
         case 'C':
            rates[CANCER][GROW] = argtof( optarg );
            break;
         case 'c':
            rates[CANCER][KILL] = argtof( optarg );
            break;
         case 'I':
            rates[INFECT_INT][GROW] = argtof( optarg );
            break;
         case 'i':
            rates[INFECT_INT][KILL] = argtof( optarg );
            break;
         case 'J':
            rates[INFECT_EXT][GROW] = argtof( optarg );
            break;
         case 'j':
            rates[INFECT_EXT][KILL] = argtof( optarg );
            break;
         case 'g':
            print_Grid = 1;
            break;
         case 'S':
            print_Summary = 1;
            break;
         case 'R':
            rng_Seed = (int) argtof( optarg ); 
            break;
         case 's':
            print_Stats = 1;
            break;
         case 't':
            print_Tumor_Sizes = 1;
            break;
         case 'W':
            print_Watch = 1;
            break;
	 case 'P':
	    percent_Interior = argtof( optarg );
	    break;
         case 'o':
            output_Interval = argtof( optarg );
            break;
	 case 'm':
	    t_Max_Simulation = argtof( optarg );
	    break;
	 case 'v':
	    t_Infect = argtof( optarg );
	    break;
         case 'V':
            percent_Infected = argtof( optarg );
	    if( percent_Infected > 100 || percent_Infected < 0 ){
	       fprintf(stderr, "PercentInfected must be (0, 100)\n");
	       return -1;
	    }
	    if( percent_Infected < 1 )
	    {
	       fprintf(stderr, "PercentInfected < 1.\n");
	       fprintf(stderr, "Infecting <1%% of Cancer cells\n");
	    }
            break;
	 case 'T':
	    if     ( strcmp(optarg, "RANDOM")    == 0 ) 
	       infect_Type = RANDOM;
	    else if( strcmp(optarg, "CENTER")    == 0) 
	       infect_Type = CENTER;
	    else if( strcmp(optarg, "PERIMETER") == 0 ) 
	       infect_Type = PERIMETER;
	    else if( strcmp(optarg, "MULTI")     == 0 ) 
	       infect_Type = MULTINODE;
	    else 
	    {
	       fprintf(stderr, "Invalid Infection Type\n");
	       return 0;
	    }
	    break;
         case 'f':
            if( optarg )
               print_File = optarg;
            else
               print_File = "Lattice_State_";
            break;
         case 'F':
            print_Files = 1;
            break;
         case '?':
            printf("Unknown Option.\n");
            printf("Usage: %s [OPTIONS] DAT_FILE NET_FILE XYZ_FILE\n", 
		  argv[0]);
            printf("Try '%s --help' for more information.\n", argv[0]);
            return 0;
         default:
            return 0;

      }
   }
   if((argc - optind) != 3 && (argc-optind) != 4 ) {
      printf("Usage: %s [OPTIONS] DAT_FILE NET_FILE XYZ_FILE\n", argv[0]);
      printf("Try '%s --help' for more information.\n", argv[0]);
      return 0;
   }
   data    = argv[optind  ];
   network = argv[optind+1];
   coords  = argv[optind+2];
   /* Arguments Processed */ 

   /* Setup Simulation */
   if( print_Summary ) {
      printf(" %d, %f, %f, %f, %f, ", rng_Seed, 
                	    rates[NORMAL][GROW], rates[NORMAL][KILL], 
	                    rates[CANCER][GROW], rates[CANCER][KILL] );
      fflush(stdout);
   }
   sgenrand(rng_Seed);
   int return_Val = 0;
   return_Val += ReadLattice( data, network, coords );
   if( return_Val != 0 ) {
      fprintf( stderr, "Problem reading in Lattice\n");
      return( -1 );
   }

   /* Run Simulation until Normal Equilibrium is Reached */
   if( do_Equalization == 1 ){
      NormalEquilibrium();
      if( print_Stats ) {
	 printf(" \n");
	 printf("  Normal Cell Equilibrium Stopped!!\n");
	 printf("  Equilibrium Value: %20.17f     \n", 
		  (double) lttc->n_Cells[NORMAL]/lttc->n_Nodes);
	 printf("  Restating Time and Seeding Cancer\n");
	 printf(" \n");
      }
   }

   /* Time (re)starts when Cancer is seeded */
   lttc->time = 0;
   if( cancer_Count > 0 )
      SeedCancer( lttc->nodes[lttc->center_Idx], cancer_Count, cancer_Radius);
   
   float time_Extinct=0;
   /*
   int time_Check = 1;
   int time_Check_fpl = 1;
   */
   if( print_Files == 1 && t_Infect != 0 ) {
      OpenFiles();
   }
   if( print_File && t_Infect != 0 ) {
      char filename[MAX_LINE_LEN];
      snprintf(filename, MAX_LINE_LEN, "%s%.1f.csv",
		     print_File, lttc->time);
      OutputLatticeState( filename );
   }
   

   /* Run for t_Infect Days (default 7) before injecting virus */
   int success;
   while( lttc->time < t_Infect ) {
      if( lttc->n_Cells[CANCER] == 0 ) {
         printf("\n");
         fprintf(stderr, "Cancer Cells died out before injecting virus\n");
         exit(-1);
      }
      if( lttc->n_Cells[NORMAL] == 0 ) {
         printf("\n");
         fprintf(stderr, "Normal Cells died out before injecting virus\n");
         exit(-1);
      }
      success = Simulate( lttc );
      if ( success == -1 ) {
         fprintf(stderr, "Simulate returned -1.\n");
      }
   }
   if( infect_Type != NONE ){
      InjectVirus( lttc );
      if( print_Grid ) {
	OutputGrid();
      }
   }

   if( print_Files == 1 && t_Infect == 0 ) {
      OpenFiles();
   }
   if( print_File && t_Infect == 0 ) {
      char filename[MAX_LINE_LEN];
      snprintf(filename, MAX_LINE_LEN, "%s%.1f.csv",
		     print_File, lttc->time);
      OutputLatticeState( filename );
   }

   /* Run all but last days of Simulation */
   int n_LinReg = 10;
   while( lttc->time <= t_Max_Simulation - n_LinReg*output_Interval ) {
      unsigned long int n_Infect = lttc->n_Cells[INFECT_INT] + 
	                           lttc->n_Cells[INFECT_EXT] ;
      if( print_Files == 1 ) {
         /* Stop if Cancer or Virus Dies Out */
         if( lttc->n_Cells[CANCER] == 0 && 
                      time_Extinct == 0 ) { 
            time_Extinct = lttc->time; 
         }
         if(      percent_Infected  > 0 && 
                          n_Infect == 0 && 
                      time_Extinct == 0   ) {
            time_Extinct = lttc->time;
         }
      }
      if( ( ( lttc->n_Cells[CANCER] == 0 && n_Infect == 0 ) || 
	    ( percent_Infected > 0 && n_Infect == 0       )    ) 
	    && lttc->time > 126 ) {
          if( print_Stats ) {
             OutputStats( lttc );
          }
	  if( print_Grid ) {
             OutputGrid();
	  }


         if( print_Files == 1 ) {
            ParameterFile( time_Extinct );
         }
	 break;
      }


      success = Simulate( lttc );
      if ( success == -1 ) {
         fprintf(stderr, "Simulate returned -1.\n");
      }
   }

   /* Run last days of Simulation */
   int timeIt = 0;
   int typeIt;
   float output_Timer = t_Max_Simulation - n_LinReg*output_Interval ;
   float *X_Vals = calloc( n_LinReg, sizeof(float) );
   float *Y_Vals[N_CELL_TYPES];
   for( typeIt = NORMAL; typeIt < N_CELL_TYPES; typeIt++ ) {
      Y_Vals[typeIt] = calloc( n_LinReg, sizeof(float) ) ; 
   }
   while( lttc->time < t_Max_Simulation ) {
      /* Save Time & n_Cells every output_Interval */
      if( lttc->time > output_Timer ) {
         output_Timer += output_Interval;
	 X_Vals[timeIt] = (float) lttc->time;
	 for( typeIt = NORMAL; typeIt <= INFECT_EXT; typeIt++ ) {
	    Y_Vals[typeIt][timeIt] = (float) lttc->n_Cells[typeIt];
	 }
	 timeIt++;
      }
      /* Stop if Cancer or Virus Dies Out */
      unsigned long int n_Infect = lttc->n_Cells[INFECT_INT] + 
	                           lttc->n_Cells[INFECT_EXT] ;
      if( ( lttc->n_Cells[CANCER] == 0 && n_Infect == 0 ) || 
	  (       percent_Infected > 0 && n_Infect == 0 )    ) {
	 free(X_Vals);
	 X_Vals = NULL;
	 for( typeIt = NORMAL; typeIt <= INFECT_EXT; typeIt++ ) {
	    free(Y_Vals[typeIt]);
	    Y_Vals[typeIt] = NULL;
	 }

         if( print_Files == 1 ) {
            CloseFiles();
         }

	 return 0;
      }
      success = Simulate( lttc );
      if ( success == -1 ) {
         fprintf(stderr, "Simulate returned -1.\n");
      }
   }

   /* Simulation Complete -- Wrap-Up */
   if( print_Files == 1 ) {
      ParameterFile( 0.0 );
      CloseFiles();
   }
   if( print_Summary )
   {
      float slope;
      printf(" Slopes: ");
      for( typeIt = NORMAL; typeIt <= N_CELL_TYPES; typeIt++ )
      {
         LinReg( X_Vals, Y_Vals[typeIt], n_LinReg, &slope, NULL );
         printf("%f ", slope);
      }
   
      printf(" Counts: "); 
      for( typeIt = NORMAL; typeIt <= N_CELL_TYPES; typeIt++ )
	 printf("%lu ", lttc->n_Cells[typeIt]); 
      printf("\n");
   }
   free(X_Vals);
   X_Vals = NULL;
   for( typeIt = NORMAL; typeIt < N_CELL_TYPES; typeIt++ ) {
      free(Y_Vals[typeIt]);
      Y_Vals[typeIt] = NULL;
   }

   DeleteLattice( lttc );
   return 0;

}/* End of main() */

void NormalEquilibrium(){
   int success;
   float convergance_Interval = .5;
   float output_Timer = convergance_Interval;
   int timer_Reset = 0;
   int n_Avg = 5;
   int n_Cnt = 30;
   int *converge = calloc ( n_Avg*2, sizeof(int) );
   int It, conv_It;
   for( It = 0; 1; It++ ) {
      if( lttc->n_Cells[NORMAL] == 0 ) {
         printf("\n");
         fprintf(stderr, "Normal Cells died out prior to equilibirum.\n");
         free( converge );
	 exit(-1);
      }
      success = Simulate( lttc );
      if ( success == -1  ) {
         printf("\n");
         fprintf(stderr, "Simulate returned -1.\n");
	 continue;
      }
      timer_Reset=0;

      /* Convergance Check.  If found: continue 30 steps & return */
      if( lttc->time > output_Timer ) {
         output_Timer += convergance_Interval;
         timer_Reset = 1;
      }
      if( timer_Reset == 1 ) {
         converge[conv_It % (2*n_Avg)] = lttc->n_Cells[NORMAL];
         if( conv_It  == 9 ) {
            int i, sum_A = 0, sum_B = 0;
            sum_A = 0;
            sum_B = 0;
            for( i = 0; i < n_Avg; i++ ) {
               sum_A += converge[ i ];
               sum_B += converge[ i + n_Avg ];
            }

            if( sum_A < sum_B ) {
               if( print_Stats ) {
                  printf(" \n  Normal Cell Equilibrium Flagged!! \
                        %20.17f      %20.17f     \n \n",
                        (double) sum_A/(n_Avg*lttc->n_Nodes), 
                        (double) sum_B/(n_Avg*lttc->n_Nodes));
               }
               int j=0;
               do{
                  success = Simulate( lttc );
                  if ( success == -1 ) {
		     fprintf(stderr, "Simulate returned -1.\n");
                     continue;
                  }
                  if( lttc->time > output_Timer ) {
                     output_Timer += convergance_Interval;
                     timer_Reset = 1;
                  }
                  if( timer_Reset == 1 )
                     j++;
                  timer_Reset = 0;
               } while ( j <= n_Cnt );
               break;
            }
         }
         conv_It = (conv_It+1) % (2*n_Avg); 
      }
   }

   free( converge );
   converge = NULL;
   return;
}

/** Usage()****************************************************
* Print usage information
**************************************************************/
void Usage () {
  fprintf(stderr, "Usage: Cells [option(s)] datFILE, netFILE, coordFILE\n");
  fprintf(stderr, " The simulator options are:\n\
  -N --NormalGrowth=NUM     Growth rate of normal cells (default 0)\n\
  -n --NormalDeath=NUM      Death  rate of normal cells (default 0)\n\
  -C --CancerGrowth=NUM     Growth rate of cancer cells (default 0)\n\
  -c --CancerDeath=NUM      Death  rate of cancer cells (default 0)\n\
  -I --ViralIntGrowth=NUM   Growth rate of Interior viral  cells (default 0)\n\
  -i --ViralIntDeath=NUM    Death  rate of Interior viral  cells (default 0)\n\
  -J --ViralExtGrowth=NUM   Growth rate of Exterior viral  cells (default 0)\n\
  -j --ViralExtDeath=NUM    Death  rate of Exterior viral  cells (default 0)\n\
  -R --RandomSeed=NUM       Seed for random number generator (default 11)\n\
  -E --Equalize             Don't start time until normal population is stable\n\
  -r --CancerRadius=NUM     Cancer ball radius in neighbors (default 1)\n\
  -q --CancerCount=NUM      Number of cancer nodes to seed (default 1) \n\
  -v --TimeInfect=INTEGER   Time between cancer and virus seeding (default 7)\n\
  -V --PercentInfected=NUM  Percentage of cancer nodes that become viral\n\
  -T --InfectionType=TYPE   TYPE=[RANDOM|CENTER|PERIMETER|MULTINODE]\n\
  -P --PercentInterior=NUM  Percent of Cancer neighbors for interior\n\
  -m --TimeMax=INTEGER      Maximum time the simulator should run (default 75)\n\
The output options are:\n\
  -o --OutputInterval=NUM   Set output frequency (default .5)\n\
  -g --Grid                 Print square 2d grid to stdout.\n\
  -s --Stats                Print stats to stdout.\n\
  -S --Summary              Print summary to stdout.\n\
  -t --TumorSizes           Print number of tumors and their size\n\
  -W --Watch                Wait for return key after every output\n\
  -F --Files                Output to vpts.dat, vlpt.dat, and vout.dat\n\
  -f --File=BASENAME        Save lattice states to BASENAME_time.csv\n");
  exit(-1);
}/* End of Usage() */

