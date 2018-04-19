#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include "common.h"


#define MAX_ENTRIES 100 

//
//  benchmarking program
//
int main( int argc, char **argv )
{
    int n[MAX_ENTRIES],i,count=0,num,p[MAX_ENTRIES];
    double t[MAX_ENTRIES],slope[MAX_ENTRIES-1],ss[MAX_ENTRIES],sse[MAX_ENTRIES],ws[MAX_ENTRIES],grade,ssgrade,wsgrade,sse_avg,ws_avg;
    double lt[MAX_ENTRIES],ln[MAX_ENTRIES],b2,sx=0.0,sx2=0.0,sxy=0.0,sy=0.0;

    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help \n" );
        printf( "-s <filename> to specify name of summary file \n");
        printf( "-v to specify what to autograde (serial,pthreads,openmp,mpi) \n" );
        return 0;
    }
    
    char *savename = read_string( argc, argv, "-s", NULL );
    FILE *fread = savename ? fopen( savename, "r" ) : NULL;

    char *autoname = read_string( argc, argv, "-v", NULL );
     
    if (strcmp(autoname,"serial")==0){
      if(fread) 
        while( fscanf (fread,"%d %lf",&n[count],&t[count]) != EOF ) 
          count++;
     
      for (i=0; i<count-1;i++) {
        slope [i] = ( log(t[i+1]) - log(t[i]) ) / ( log(n[i+1]*1.0) - log(n[i]*1.0) );
      }
      for (i=0; i<count; i++) {
        lt[i] = log(t[i]);
        ln[i] = log(n[i]*1.0);
      }
      for (i=0; i<count; i++) {
        sx += ln[i];
        sy += lt[i];
        sxy += ln[i]*lt[i];
        sx2 += ln[i]*ln[i];
      }

      b2 = (sxy - (sx * sy)/ (count * 1.0) ) / (sx2 - (sx * sx)/ (count * 1.0));
 
      printf("\nSerial code is O(N^slope)");
      printf("\nSlope estimates are :");
      for (i=0; i<count-1; i++){
        printf(" %lf",slope[i]);
      }
      printf("\nSlope estimate for line fit is: %lf\n", b2);
	  
	  if (b2 < 1.3) grade = 100.00;
	  else if (b2 < 1.5) grade = 75.00 + (1.5-b2)/(0.2) * 25.00;
	  else if (b2 < 2) grade = (2-b2)/0.5 * 75.00;
	  else grade =0.0;
	  
	  printf("Serial Grade = %7.2lf",grade);
      printf("\n\n");
    }

    if (strcmp(autoname,"pthreads")==0 || strcmp(autoname,"openmp")==0 || strcmp(autoname,"mpi")==0){
      if(fread){
        fscanf (fread,"%d %lf",&n[count],&t[count]);
        count++; p[0]=1;
        while( fscanf (fread,"%d %d %lf",&n[count],&p[count],&t[count]) != EOF )
          count++;
      }
 
      num = count/2;
 
      ss[0] = sse[0] = ws[0] = t[0]/t[1];
      for (i=2; i<=num;i++) {
         ss[i-1] = t[0]/t[i];
		 sse[i-1] = ss[i-1]/p[i];
         ws[i-1] = t[0]/t[i+num-1];
      }
 
      printf("\nStrong scaling estimates are :\n");
      for (i=0; i<num; i++){
        printf(" %7.2lf",ss[i]);
      }
      printf(" (speedup)\n");
      for (i=0; i<num; i++){
        printf(" %7.2lf",sse[i]);
      }
      printf(" (efficiency)    for\n");
      for (i=0; i<num; i++){
        printf(" %7d",p[i+1]);
      }
      printf(" threads/processors\n\n");

	  sse_avg=0.0;
	  for (i=0; i<num; i++){
        sse_avg+=sse[i];
      }
	  sse_avg/=num;
	  
	  printf("Average strong scaling efficiency: %7.2lf \n\n",sse_avg);
	  
      printf("Weak scaling estimates are :\n");
      for (i=0; i<num; i++){
        printf(" %7.2lf",ws[i]);
      }
      printf(" (efficiency)    for\n");
      for (i=0; i<num; i++){
        printf(" %7d",p[i+1]);
      }
      printf(" threads/processors\n\n");

	  ws_avg=0.0;
	  for (i=0; i<num; i++){
        ws_avg+=ws[i];
      }
	  ws_avg/=num;
	  
	  printf("Average weak scaling efficiency: %7.2lf \n\n",ws_avg);

	  if (sse_avg > 0.8) ssgrade=100.00;
	  else if (sse_avg > 0.5) ssgrade=75.00 + 25.00 * (sse_avg-0.5)/0.3; 
	  else ssgrade = sse_avg/0.5 * 75.00;
	  
	  if (ws_avg > 0.8) wsgrade=100.00;
	  else if (ws_avg > 0.5) wsgrade=75.00 + 25.00 * (ws_avg-0.5)/0.3; 
	  else wsgrade = ws_avg/0.5 * 75.00;
	  
	  grade= 0.5 * ssgrade + 0.5 * wsgrade;
	  
	  printf("\n%s Grade = %7.2f\n\n",autoname,grade);
    }

    fclose( fread );
    
    return 0;
}
