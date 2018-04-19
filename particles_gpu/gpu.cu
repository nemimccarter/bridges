#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <cuda.h>
#include "common.h"

#define NUM_THREADS 256
#define b(x, y, p, n) (bins[(x)+(y)*(n)+(p)*(n)*(n)])
#define c(x, y, n) (counts[(x)+(y)*(n)]) // check if there are particles in a bin

int bs; // bin size
int np; // number of max particles per bin
int* bins;
int* counts;

int* d_bins;
int* d_counts;
/*
cudaMalloc((void**) &d_bins, bs * bs * np * sizeof(int));
cudaMalloc((void**) &d_counts, bs * bs * sizeof(int));
cudaMemset(d_counts, 0, bs * bs * sizeof(int));
*/
extern double size;

int get_binsize()
{
  return (int)(size / cutoff) + 1;
}

//
//  benchmarking program
//
__global__ void assign_bins_gpu(particle_t * particles, int * bins, int * counts,
int n, int bs)
{
 int tid = threadIdx.x + blockIdx.x * blockDim.x;
 if(tid >= n) return;
 int x = (int)(particles[tid].x/cutoff);
 int y = (int)(particles[tid].y/cutoff);
 int new_c = atomicAdd(&c(x, y, bs), 1);
 b(x, y, new_c, bs) = tid;
}

__device__ void apply_force_gpu(particle_t &particle, particle_t &neighbor)
{
  double dx = neighbor.x - particle.x;
  double dy = neighbor.y - particle.y;
  double r2 = dx * dx + dy * dy;
  if( r2 > cutoff*cutoff )
      return;
  //r2 = fmax( r2, min_r*min_r );
  r2 = (r2 > min_r*min_r) ? r2 : min_r*min_r;
  double r = sqrt( r2 );

  //
  //  very simple short-range repulsive force
  //
  double coef = ( 1 - cutoff / r ) / r2 / mass;
  particle.ax += coef * dx;
  particle.ay += coef * dy;

}

__global__ void compute_forces_gpu(particle_t * particles, int n)
{

  int id = threadIdx.x + blockIdx.x * blockDim.x;
  
  if(id >= bs*bs) 
    return; // return if out of bounds
  
  int bx = id % bs; // get the x bin position
  int by = id / bs; // get the y bin position
  int nx, ny;
  
  // Get thread (particle) ID
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if(tid >= n) return;

  particles[tid].ax = particles[tid].ay = 0;
  for(int j = 0 ; j < n ; j++)
    apply_force_gpu(particles[tid], particles[j]);

}

__global__ void compute_forces_bin_gpu(particle_t * particles, int bs, int * bins,
int * counts)
{
 int id = threadIdx.x + blockIdx.x * blockDim.x;
 if(id >= bs*bs) return; // return if out of bounds
 int bx = id % bs; // get the x bin position
 int by = id / bs; // get the y bin position
 int nx, ny;
 // … work goes here
}

__global__ void move_gpu (particle_t * particles, int n, double size)
{

  // Get thread (particle) ID
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if(tid >= n) return;

  particle_t * p = &particles[tid];
    //
    //  slightly simplified Velocity Verlet integration
    //  conserves energy better than explicit Euler method
    //
    p->vx += p->ax * dt;
    p->vy += p->ay * dt;
    p->x  += p->vx * dt;
    p->y  += p->vy * dt;

    //
    //  bounce from walls
    //
    while( p->x < 0 || p->x > size )
    {
        p->x  = p->x < 0 ? -(p->x) : 2*size-p->x;
        p->vx = -(p->vx);
    }
    while( p->y < 0 || p->y > size )
    {
        p->y  = p->y < 0 ? -(p->y) : 2*size-p->y;
        p->vy = -(p->vy);
    }

}



int main( int argc, char **argv )
{   
    bs = get_binsize();
    np = bs / 10;
 
    // This takes a few seconds to initialize the runtime
    cudaThreadSynchronize(); 

    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set the number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        return 0;
    }
    
    int n = read_int( argc, argv, "-n", 1000 );

    char *savename = read_string( argc, argv, "-o", NULL );
    
    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );

    // GPU particle data structure
    particle_t * d_particles;
    cudaMalloc((void **) &d_particles, n * sizeof(particle_t));

    set_size( n );

    init_particles( n, particles );

    cudaThreadSynchronize();
    double copy_time = read_timer( );

    // Copy the particles to the GPU
    cudaMemcpy(d_particles, particles, n * sizeof(particle_t), cudaMemcpyHostToDevice);

    cudaThreadSynchronize();
    copy_time = read_timer( ) - copy_time;
    
    //
    //  simulate a number of time steps
    //
    cudaThreadSynchronize();
    double simulation_time = read_timer( );

    // assign bins    

    for( int step = 0; step < NSTEPS; step++ )
    {
	// assign bins
    	assign_bins_gpu <<< blks, NUM_THREADS >>> (d_particles, d_bins, d_counts, n, bs);
        
	clear_accel_gpu <<< blks, NUM_THREADS >>> (d_particles, n);
	//
        //  compute forces
        //

	int bin_blks = (bs * bs + NUM_THREADS - 1) / NUM_THREADS;
	int blks = (n + NUM_THREADS - 1) / NUM_THREADS;
	compute_forces_gpu <<< blks, NUM_THREADS >>> (d_particles, n);
        
        //
        //  move particles
        //
	move_gpu <<< blks, NUM_THREADS >>> (d_particles, n, size);
        int bin_blks = (bs * bs + NUM_THREADS - 1) / NUM_THREADS;
	compute_forces_bin_gpu <<< bin_blks, NUM_THREADS >>> (d_particles, bs, d_bins, d_counts);

	move_gpu <<< blks, NUM_THREADS >>> (d_particles, n, size); 
        //
        //  save if necessary
        //
        if( fsave && (step%SAVEFREQ) == 0 ) {
	    // Copy the particles back to the CPU
            cudaMemcpy(particles, d_particles, n * sizeof(particle_t), cudaMemcpyDeviceToHost);
            save( fsave, n, particles);
	}
    }
    cudaThreadSynchronize();
    simulation_time = read_timer( ) - simulation_time;
    
    printf( "CPU-GPU copy time = %g seconds\n", copy_time);
    printf( "n = %d, simulation time = %g seconds\n", n, simulation_time );
    
    free( particles );
    cudaFree(d_particles);
    if( fsave )
        fclose( fsave );
    
    return 0;
}