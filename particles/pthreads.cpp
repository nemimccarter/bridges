#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <pthread.h>
#include "common.h"
#include <vector>

// constants for bin solution
#define density 0.0005
#define mass 0.01
#define cutoff 0.01
#define min_r (cutoff / 100)
#define dt 0.0005

//
//  global variables
//
int n, n_threads,no_output=0;
particle_t *particles;
FILE *fsave,*fsum;
pthread_barrier_t barrier;
pthread_mutex_t mutex=PTHREAD_MUTEX_INITIALIZER;
double gabsmin=1.0,gabsavg=0.0;
std::vector<particle_t*> *bins;
int bpr, numbins;
pthread_mutex_t *binsMutex;

//
//  check that pthreads routine call was successful
//
#define P( condition ) {if( (condition) != 0 ) { printf( "\n FAILURE in %s, line %d\n", __FILE__, __LINE__ );exit( 1 );}}

//calculate particle's bin number
int binNum(particle_t &p, int bpr)
{
    return (floor(p.x / cutoff) + bpr * floor(p.y / cutoff));
}

//
//  This is where the action happens
//
void *thread_routine( void *pthread_id )
{
    int navg,nabsavg=0;
    double dmin,absmin=1.0,davg,absavg=0.0;
    int thread_id = *(int*)pthread_id;

    int particles_per_thread = (numbins + n_threads - 1) / n_threads;
    int first = min(thread_id * particles_per_thread, numbins);
    int last = min((thread_id + 1) * particles_per_thread, numbins);
 
//    int particles_per_thread = (n + n_threads - 1) / n_threads;
//    int first = min(  thread_id    * particles_per_thread, n );
//    int last  = min( (thread_id+1) * particles_per_thread, n );
    
    //
    //  simulate a number of time steps
    //
    for( int step = 0; step < NSTEPS; step++ )
    {
        dmin = 1.0;
        navg = 0;
        davg = 0.0;
        //
        //  compute forces
        //
        for( int particle_index = first; particle_index < last; particle_index++ )
        {
            particles[particle_index].ax = particles[particle_index].ay = 0;
            
            // find current particle's bin, handle boundaries
            int cbin = binNum(particles[particle_index], bpr);
            int lowi = -1, highi = 1, lowj = -1, highj = 1;

            if (cbin < bpr)
                lowj = 0;
            if (cbin % bpr == 0)
                lowi = 0;
            if (cbin % bpr == (bpr - 1))
                highi = 0;
            if (cbin >= bpr * (bpr - 1))
                highj = 0;

            for (int i = lowi; i <= highi; i++ )
                for (int j = lowj; j <= highj; j++)
                {
                    int nbin = cbin + i + bpr * j;
                    for (int k = 0; k < bins[nbin].size(); k++)
                        apply_force(particles[particle_index], *bins[nbin][k], &dmin, &davg, &navg);
                }
            }
        
        pthread_barrier_wait( &barrier );
        
        if( no_output == 0 )
        {
          //
          // Computing statistical data
          // 
          if (navg) {
            absavg +=  davg/navg;
            nabsavg++;
          }
          if (dmin < absmin) absmin = dmin;
	}

        //
        //  move particles
        //
        for( int i = first; i < last; i++ ) 
            move( particles[i] );
        
        pthread_barrier_wait( &barrier );
        
        //
        //  save if necessary
        //
        if (no_output == 0) 
          if( thread_id == 0 && fsave && (step%SAVEFREQ) == 0 )
            save( fsave, n, particles );
        
    }
     
    if (no_output == 0 )
    {
      absavg /= nabsavg; 	
      //printf("Thread %d has absmin = %lf and absavg = %lf\n",thread_id,absmin,absavg);
      pthread_mutex_lock(&mutex);
      gabsavg += absavg;
      if (absmin < gabsmin) gabsmin = absmin;
      pthread_mutex_unlock(&mutex);    
    }

    return NULL;
}

//
//  benchmarking program
//
int main( int argc, char **argv )
{    
    //
    //  process command line
    //
    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set the number of particles\n" );
        printf( "-p <int> to set the number of threads\n" );
        printf( "-o <filename> to specify the output file name\n" );
        printf( "-s <filename> to specify a summary file name\n" );
        printf( "-no turns off all correctness checks and particle output\n");        
        return 0;
    }
    
    n = read_int( argc, argv, "-n", 1000 );
    n_threads = read_int( argc, argv, "-p", 2 );
    char *savename = read_string( argc, argv, "-o", NULL );
    char *sumname = read_string( argc, argv, "-s", NULL );

    fsave = savename ? fopen( savename, "w" ) : NULL;
    fsum = sumname ? fopen ( sumname, "a" ) : NULL;

    if( find_option( argc, argv, "-no" ) != -1 )
      no_output = 1;

    // init bins and mutex
    bins = new std::vector<particle_t*>[numbins];
    binsMutex = new pthread_mutex_t[numbins];
    for (int m = 0; m < numbins; m++)
        pthread_mutex_init(&binsMutex[m], NULL);

    //
    //  allocate resources
    //
    particles = (particle_t*) malloc( n * sizeof(particle_t) );
    set_size( n );
    init_particles( n, particles );

    pthread_attr_t attr;
    P( pthread_attr_init( &attr ) );
    P( pthread_barrier_init( &barrier, NULL, n_threads ) );

    int *thread_ids = (int *) malloc( n_threads * sizeof( int ) );
    for( int i = 0; i < n_threads; i++ ) 
        thread_ids[i] = i;

    pthread_t *threads = (pthread_t *) malloc( n_threads * sizeof( pthread_t ) );

    //
    //  do the parallel work
    //
    double simulation_time = read_timer( );
    for( int i = 1; i < n_threads; i++ ) 
        P( pthread_create( &threads[i], &attr, thread_routine, &thread_ids[i] ) );
    
    thread_routine( &thread_ids[0] );
    
    for( int i = 1; i < n_threads; i++ ) 
        P( pthread_join( threads[i], NULL ) );
    simulation_time = read_timer( ) - simulation_time;
   
    printf( "n = %d, simulation time = %g seconds", n, simulation_time);

    if( find_option( argc, argv, "-no" ) == -1 )
    {
      gabsavg /= (n_threads*1.0);
      // 
      //  -the minimum distance absmin between 2 particles during the run of the simulation
      //  -A Correct simulation will have particles stay at greater than 0.4 (of cutoff) with typical values between .7-.8
      //  -A simulation were particles don't interact correctly will be less than 0.4 (of cutoff) with typical values between .01-.05
      //
      //  -The average distance absavg is ~.95 when most particles are interacting correctly and ~.66 when no particles are interacting
      //
      printf( ", absmin = %lf, absavg = %lf", gabsmin, gabsavg);
      if (gabsmin < 0.4) printf ("\nThe minimum distance is below 0.4 meaning that some particle is not interacting ");
      if (gabsavg < 0.8) printf ("\nThe average distance is below 0.8 meaning that most particles are not interacting ");
    }
    printf("\n");

    //
    // Printing summary data
    //
    if( fsum)
        fprintf(fsum,"%d %d %g\n",n,n_threads,simulation_time); 
    
    //
    //  release resources
    //
    P( pthread_barrier_destroy( &barrier ) );
    P( pthread_attr_destroy( &attr ) );
    free( thread_ids );
    free( threads );
    free( particles );
    if( fsave )
        fclose( fsave );
    if( fsum )
        fclose ( fsum );
    
    return 0;
}
