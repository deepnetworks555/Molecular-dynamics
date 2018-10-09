//simulation of a number of particles (N) which have  Lennard-jones interaction with each other
//g++ -Wall -o moldynamics moldynamics.cpp
//Fall 2012
#include <iostream>
#include <stdio.h>		//for writting to file
#include <math.h>
#include <cstdlib>
#include <cstdio>
#include "gasdev.h"
#include "ran1.h"
using namespace std;


//+++++++++++Parameters++++++++++++++
const int N=300;			//number of particles 300
const double Rho=0.83;			//density
const double dt = 0.001;    		//time step
const double tempref=1.5;		//temperature
//++++++++++++++++++++++++++++++++++++	

double length;		                  	  //length of each side  of the cubic box----sigma is unit of length
double rcut = 2.5;           	 	          //cut-off radius 	
double position[N][3],position0[N][3];  	  //positions
double velocity[N][3],velocity0[N][3];		   //velocities
double acceleration[N][3];	  		  //accelerations
int stepnumber=0;
double temperature();
double velocityscale();
double PE=0;
const double pi=3.14159265358979;
const int nBins = 300;
const int nhis=300;
int vBins[nBins];
double dv;
double g[nhis];
long idum=(-24);

 		
void initialize(double rho)
{
	cout<< " Initializing...\n" ;
        length = pow((N / rho) , 1.0/3);
        cout <<"\nbox length = "<< length<<"\n";
	int n3 = 1;

/*for simple cube
	while ( n3*n3*n3 < N ) n3++; 
	//initialize positions
	int ix=0,iy=0,iz=0;
	for (int i=0; i < N; i++){
		position[i][0] = (ix + 0.5) * length/n3;
		position[i][1] = (iy + 0.5) * length/n3;
		position[i][2] = (iz + 0.5) * length/n3;


		ix++;
		if (ix == n3){
		    ix=0;
		    iy++;
		    if (iy == n3){
		      iy=0;
		      iz++;
		    }
		}
	}
*/

	//**********fcc latice***********************
	while (4 * n3 * n3 * n3 < N) { ++n3; }
	double a = length / n3;
	//assign positions of the atoms in basic cell
	double xBasic[4]={0.25 ,0.75 ,0.75 ,0.25 };
	double yBasic[4]={0.25 ,0.75 ,0.25 ,0.75 };
	double zBasic[4]={0.25 ,0.25 ,0.75 ,0.75 };
	

	int n=0;
	for(int x =0;x < n3; x++)
	    for(int y =0;y < n3; y++)
		for(int z =0;z < n3; z++)
		    for(int k=0;k < 4; k++)
			if(n < N){
		    		position[n][0]=(x + xBasic[k]) * a;
		    		position[n][1]=(y + yBasic[k]) * a;
		    		position[n][2]=(z + zBasic[k]) * a;
				++n;
		    	}

	//**********************************initializing velocities**************************	 

	for (int i = 0; i < N; i++)
	     for (int k = 0; k < 3; k++)
		  velocity[i][k] = gasdev(&idum);		//using box-Muller for gaussian

	//***************changing velocities so that the velocity of Center Of Mass becomes zero*************
	double cmv[3] = { 0 , 0 , 0};
	for (int i=0; i < N; i++)
	    for ( int k=0; k < 3; k++)
	         cmv[k] += velocity[i][k];
	cmv[0] /=N;
	cmv[1] /=N;
	cmv[2] /=N;
	for (int i = 0; i < N; i++)
	     for (int k = 0; k < 3; k++)
		 velocity[i][k] -= cmv[k];

	printf("velocity of CM before: (%1.1e,%1.1e,%1.1e)\n",cmv[0],cmv[1],cmv[2]);
	for (int i=0; i < N; i++)
	  		  for ( int k=0; k < 3; k++)
	      		   cmv[k] += velocity[i][k];
			cmv[0] /=N;
			cmv[1] /=N;
			cmv[2] /=N;

	printf("velocity of CM after setting it to zero: (%1.1e,%1.1e,%1.1e)\n",cmv[0],cmv[1],cmv[2]);
	//******************************************************************
	
	//scale velocities for the desired temperature
	velocityscale();	 


}

void computeAccelerations()
{
     PE=0;
     for (int i=0; i < nhis;i++) {g[i]=0;}
     for (int i=0; i < N; i++)			//set all accelerations to zero
	 for (int k=0; k < 3; k++)
	     acceleration[i][k] = 0;    

     
     for (int i=0; i < N-1; i++)
	  for (int j = i+1; j < N; j++) {
	       double f = 0;
 	       double rij[3];			//dictance between two particles
    	       double rSqr = 0;			//square of dictance between two particles
	       for (int k = 0; k < 3; k++){
		    rij[k] = position[i][k] - position[j][k];
		    //Minimum image convention
		    if ( rij[k] > (0.5 * length))  {rij[k] -= length;}
		    if ( rij[k] < (-0.5 * length)) {rij[k] += length;}
		    
		    rSqr += (rij[k] * rij[k]);

	       }
	       //Calculating potential energy
	       double r2 = 1.0 / rSqr;
	       double ri6 = (rSqr*rSqr*rSqr);
	       double r6 = 1.0 / ri6;
	       
	       //pair distribution
	       double deltaR=(length)/nhis;
	       if (rSqr < length*length){
		   double r=sqrt(rSqr);
		   int ig = r/deltaR;
		   g[ig]+=2;
		   g[ig]/=(4 * pi * rSqr*deltaR);
		   g[ig] /=Rho*N;
	       }
	
	       if ( rSqr < rcut * rcut ){				//cosidering cut-off distance
		   f = 48 * ( r6 * r6 * r2 -0.5 * r6* r2) ;  
		   PE += 4 * (r6 * r6 - r6);
	       	   for (int k = 0; k < 3; k++){
		   	 acceleration[i][k] += f * rij[k];
		  	 acceleration[j][k] -= f * rij[k];
		   }
 	      }
           }

}

double VelocityVerlet(double dt){			//Velocity Verlet Algorithm
	computeAccelerations();
	
	for (int i=0; i < N; i++){
	    for (int k=0; k < 3; k++){
		position[i][k] += velocity[i][k] * dt + 0.5 * acceleration[i][k] * dt * dt;	//calculating new positions
		velocity[i][k] += 0.5 * acceleration[i][k] * dt;
		
		
		//Apply periodic boundary conditions:
		if (position[i][k] < 0) {position[i][k] += length; }
		if (position[i][k] > length) {position[i][k] -= length; }	
	    }
	}
	computeAccelerations();
	for (int i=0; i < N; i++)
	    for (int k=0; k < 3; k++)
		velocity[i][k] += 0.5 * acceleration[i][k] * dt;   
	    


	//Calculating kinetic energy
	double KE = 0;
	for (int i=0; i < N; i++)
	     KE += velocity[i][0]*velocity[i][0]+velocity[i][1]*velocity[i][1]+velocity[i][2]*velocity[i][2];
	KE *= 0.5;

return (KE);

}

double temperature()
{	
	//cout<< " calcualting the temperature...\n" ;
	double sum = 0;
	for (int i=0; i < N; i++)
	     sum += velocity[i][0]*velocity[i][0]+velocity[i][1]*velocity[i][1]+velocity[i][2]*velocity[i][2];
	sum /= 3*(N-1);
	//cout<< "\n  the temperature is "<<sum ;
	
return (sum);
}

double velocityscale()
{

	//***********scaling velocities for desired temperature**************
	//cout<< "scaling velocities\n " ;
	double te=temperature();
	double B = sqrt(tempref/te);
	for (int i = 0; i < N; i++)
	    for(int k = 0; k < 3; k++)
			velocity[i][k] *= B;


}

double MSD()
{   	//Mean squared displacement
	 double MSD=0;
         for (int i=0; i < N; i++){
		double x2 = (position[i][0] - position0[i][0]) * (position[i][0] - position0[i][0]);
		double y2 = (position[i][1] - position0[i][1]) * (position[i][1] - position0[i][1]);
		double z2 = (position[i][2] - position0[i][2]) * (position[i][2] - position0[i][2]);
		MSD += x2+y2+z2;
	 }
	 MSD /=N;
return (MSD);
}	   

void velocitydistribution()
{	
	double vMax=-10000;
	
	double v[N];
     	for(int i=0; i < N ;i++){
		v[i] = sqrt(velocity[i][0]*velocity[i][0] + velocity[i][1]*velocity[i][1] + velocity[i][2]*velocity[i][2]);
		//choose maximum velocity
		if (v[i] > vMax) { vMax = v[i];}
	}
	dv = vMax / nBins;
	for(int i=0; i < N ;i++){
		int Bin = v[i]/dv;
		vBins[Bin]++;
	}
	
}



int main()
{
	double cmv[3] = { 0 , 0 , 0};
	FILE * pFile1 ,* pFile2,* pFile3,* pFile4,* pFile5,* pFile6;
	pFile1 = fopen ("energy.txt","w");
	pFile4 = fopen ("graph.txt","w");
	pFile5 = fopen ("MSD.txt","w");
	pFile6 = fopen ("pairdis.txt","w");
	double KE,totE,temp;
	initialize(Rho);

	for (int i=0; i < 10000; i++){
	    KE=VelocityVerlet(dt);
	       if(i<1000)
		if ((i % 50) == 0)
			velocityscale();

	       if (i==1000)
		   for (int i=0; i < N; i++){

		position0[i][0] = position[i][0];
		position0[i][1] = position[i][1];
		position0[i][2] = position[i][2];
		
		velocity0[i][0] = velocity[i][0];
		velocity0[i][1] = velocity[i][1];
		velocity0[i][2] = velocity[i][2];
	}

	  	  	
	    totE = KE + PE;
	    temp = temperature();

		if (i==300){
			
			for (int i=0; i < N; i++)
	  		  for ( int k=0; k < 3; k++)
	      		   cmv[k] += velocity[i][k];
			cmv[0] /=N;
			cmv[1] /=N;
			cmv[2] /=N;
		}
          
	    if ((i % 50) == 0){
	        fprintf (pFile1, "%5.5d  %f  %f  %f  %f\n",i,totE,PE,KE,temp);
                if(i>1000)
                  fprintf (pFile5, "%5.5d  %f \n",i,MSD());
	    }
	   
	     if (i == 2000){
		velocitydistribution();
		for(int ii=0;ii < nBins ;ii++)
			fprintf (pFile4, "%.2f  %i \n",ii*dv,vBins[ii]);
	      }

	}

	printf("velocity of CM at step 300: (%1.1e,%1.1e,%1.1e)\n",cmv[0],cmv[1],cmv[2]);
	for(int i=0;i < nhis ;i++)
	   fprintf (pFile6, "%f  %f \n",i*length/(nhis),g[i]);
	
	
	fclose (pFile1);
	fclose (pFile4);
	fclose (pFile5);
	fclose (pFile6);




	
return 0;
}



		
