//============================================================================
// Name        : Rydberg_Wavefunction_Numerov_Integration.c
// Author      : Nikola Sibalic // Edited by Gerard Higgins 24.01.2017 for Sr+
// Version     : 0.9
// Copyright   : BSD 3 clause
// Description : Rydberg_Wavefunction_Numerov_Integration in C
//============================================================================
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define DEBUG_OUTPUT

double innerLimit;
double outerLimit;
double step;
double init1;
double init2;

int l;
double s,j;
double stateEnergy;
double alphaC;
double alpha;
int Z;
double a1,a2,a3,a4,rc;  // depends on l - determined by Python in advance

inline double EffectiveCharge(double r){
	// returns effective charge of the core - Eq. (18b) of Marinescu PRA 49 982 (1994) and Eq. (3.26) of Aymar Rev. Mod. Phys. 68 4 (1996)
	return 2.0+((double)Z-2.)*exp(-a1*r)-r*(a3+a4*r)*exp(-a2*r);
}

inline double CorePotential(double r){
    // l dependent core potential (angular momentum of e) at radius r
    // returns Core potential - Eq. (18a) of Marinescu PRA 49 982 (1994) and Eq. (3.26) of Aymar Rev. Mod. Phys. 68 4 (1996)
    return -EffectiveCharge(r)/r-alphaC/(2.*pow(r,4))*(1.-exp(-pow(r/rc,6)));
}

double commonTerm1;


inline double Potenital(double r, double step){
	// l<4
	//return CorePotential(r)+pow(alpha,2)/(2.0*pow(r,3))*commonTerm1; // here is the potential plus the spin-orbit potential alpha*l.s/(2r**3)
	double V1 = CorePotential(r);
	double V2 = CorePotential(r+step);
	double dVdr = (V2-V1)/step;
	double Vso = pow(alpha,2)*commonTerm1*dVdr/(2.0*r*pow((1.0-pow(alpha,2)*V1/2.0),2));
	return V1 + Vso;
}

inline double Potenital2(double r){
	// l>=4
    // act as if it is a Hydrogen atom
    return -2./r;
}

double commonTerm2;

inline double kfun(double r, double step){
	// with potential for l<4
	return 2.*(stateEnergy-Potenital(r, step))-commonTerm2/pow(r,2); // E = -0.5*d2/dr2 + l(l+1)/2r2 + Vcore + Vso => -d2/dr2 = 2(E-Vcore-Vso) - l(l+1)/r2
}

inline double kfun2(double r){
	// with potential for l>4
	return 2.*(stateEnergy-Potenital2(r))-commonTerm2/pow(r,2);
}


int main(int argc, char** argv) {
	// Numerov arguments: innerLimit,outerLimit,kfun,step,init1,init2
	double innerLimit,outerLimit,step,init1,init2;
	char* outputFile;

	if (argc == 19) {
		//innerLimit,outerLimit,kfun,step,init1,init2
		innerLimit = atof(argv[1]); // atof converts a string to a float/double
		outerLimit = atof(argv[2]);
		step = atof(argv[3]);
		init1 = atof(argv[4]);
		init2 = atof(argv[5]);
		outputFile = argv[6];

		// parameters for potential
		l = atoi(argv[7]); // atoi converts a string to an integer
		s = atof(argv[8]);
		j = atof(argv[9]);
		stateEnergy = atof(argv[10]);
		alphaC = atof(argv[11]);
		alpha = atof(argv[12]);
		Z = atoi(argv[13]);
		a1 = atof(argv[14]);
		a2 =  atof(argv[15]);
		a3 =  atof(argv[16]);
		a4 =  atof(argv[17]);
		rc =  atof(argv[18]);

#ifdef DEBUG_OUTPUT
		printf("innerLimit\t=\t%.3f\nouterLimit\t=\t%.3f\nstep\t=\t%.3f\ninit1\t=\t%.3f\ninit2\t=\t%.3f\n",innerLimit,outerLimit,step,init1,init2);
		printf("Outputfile:\t%s\n",outputFile);
		printf("l\t=\t%i\ns\t=\t%.1f\nj\t=\t%.1f\n",l,s,j);

		printf("stateEnergy\t=\t%.7f\nalphaC\t=\t%.3f\nalpha\t=\t%.3f\nZ\t=\t%i\n",stateEnergy,alphaC,alpha,Z);

		printf("a1\t\t%.4f\na2\t\t%.4f\na3\t\t%.4f\na4\t\t%.4f\nrc\t\t%.4f\n",a1,a2,a3,a4,rc);
#endif
	}
	else{
		printf("Wrong argument number!\n");
		return 1;
	}

	// let's speed up calculation by calculating some common terms beforehand
	commonTerm1 = (j*(j+1.0)-((double)l)*(l+1.0)-s*(s+1.))/2.0; // this is L.S
	commonTerm2 = ((double)l)*(l+1.); // this is <L.L>

	int divergencePoint;

	int totalLength =  (int)((outerLimit-innerLimit)/step);

#ifdef DEBUG_OUTPUT
	printf("Index = %i\n",totalLength);
	printf("Index should be about = %.2f\n",((outerLimit-innerLimit)/step));
#endif

	int br = totalLength;
	double* sol = (double*) malloc(br*sizeof(double)); // sol is the solutions, malloc allocates a memory block
	double* rad = (double*) malloc(br*sizeof(double)); // rad is the radial distances

	if (!sol || !rad){
		printf("Memory allocaiton failed! Aborting.");
		return 1;
	}

    // for l<4

	if (l<4){

		br = br-1; // c uses zero-numbering. We solve first at the outerlimit.
	    double r = outerLimit;
	    double step2 = step*step;
	    sol[br] = (2*(1-5.0/12.0*step2*kfun(r, step))*init1-(1+1/12.0*step2*kfun(r+step, step))*init2)/(1+1/12.0*step2*kfun(r-step, step));
		rad[br] = r;

		r = r-step;
		br = br-1;

		sol[br] = (2*(1-5.0/12.0*step2*kfun(r, step))*sol[br+1]-(1+1/12.0*step2*kfun(r+step, step))*init1)/(1+1/12.0*step2*kfun(r-step, step));
		rad[br] = r;

		double maxValue = 0;

	    double checkPoint = 0;
	    double fromLastMax = 0;

	    while (br>checkPoint){
	        br = br-1;
	        r = r-step;
	        sol[br] = (2*(1-5.0/12.0*step2*kfun(r, step))*sol[br+1]-(1+1/12.0*step2*kfun(r+step, step))*sol[br+2])/(1+1/12.0*step2*kfun(r-step, step));
	        rad[br] = r;
	        if (fabs(sol[br])>maxValue){ // fabs give the absolute value as a float (actually as a double)
	            maxValue = fabs(sol[br]);
	        }
	        else{
	            fromLastMax += 1;
	            if (fromLastMax>50){ // if it has been 50 steps are taken without having the maximum as the solution, then the while loop is stopped.
	                checkPoint = br;
	            }
	        }
	    }

	    divergencePoint = 0;
	    while ((br>0)&&(divergencePoint == 0)){
	        br = br-1;
	        r = r-step;
	        sol[br] = (2*(1-5.0/12.0*step2*kfun(r, step))*sol[br+1]-(1+1/12.0*step2*kfun(r+step, step))*sol[br+2])/(1+1/12.0*step2*kfun(r-step, step));
	        rad[br] = r;
	        if ((divergencePoint==0)&&(fabs(sol[br])>maxValue)){ // if the solution is bigger than the max value, then the solution is diverging
	            divergencePoint = br;
	            while ((fabs(sol[divergencePoint])>fabs(sol[divergencePoint+1])) && (divergencePoint<checkPoint)){
	                divergencePoint +=1;
	            }
	            if (divergencePoint>checkPoint){
	                printf("ERROR: Numerov error\n");
	                return 1;
	            }
	        }
	    }


	} // end of if l<4
	else{ //l>=4

		br = br-1;
	    double r = outerLimit;
	    double step2 = step*step;
	    sol[br] = (2*(1-5.0/12.0*step2*kfun2(r))*init1-(1+1/12.0*step2*kfun2(r+step))*init2)/(1+1/12.0*step2*kfun2(r-step));
		rad[br] = r;

		r = r-step;
		br = br-1;

		sol[br] = (2*(1-5.0/12.0*step2*kfun2(r))*sol[br+1]-(1+1/12.0*step2*kfun2(r+step))*init1)/(1+1/12.0*step2*kfun2(r-step));
		rad[br] = r;

		double maxValue = 0;

	    double checkPoint = 0;
	    double fromLastMax = 0;

	    while (br>checkPoint){
	        br = br-1;
	        r = r-step;
	        sol[br] = (2*(1-5.0/12.0*step2*kfun2(r))*sol[br+1]-(1+1/12.0*step2*kfun2(r+step))*sol[br+2])/(1+1/12.0*step2*kfun2(r-step));
	        rad[br] = r;
	        if (fabs(sol[br])>maxValue){
	            maxValue = fabs(sol[br]);
	        }
	        else{
	            fromLastMax += 1;
	            if (fromLastMax>50){
	                checkPoint = br;
	            }
	        }
	    }

	    divergencePoint = 0;
	    printf("better?\n");
	    while ((br>0)&&(divergencePoint == 0)){
	        br = br-1;
	        r = r-step;
	        sol[br] = (2*(1-5.0/12.0*step2*kfun2(r))*sol[br+1]-(1+1/12.0*step2*kfun2(r+step))*sol[br+2])/(1+1/12.0*step2*kfun2(r-step));
	        rad[br] = r;
	        if ((divergencePoint==0)&&(fabs(sol[br])>maxValue)){
	            divergencePoint = br;
	            while ((fabs(sol[divergencePoint])>fabs(sol[divergencePoint+1])) && (divergencePoint<checkPoint)){
	                divergencePoint +=1;
	            }
	            if (divergencePoint>checkPoint){
	                printf("ERROR: Numerov error\n");
	                return 0;
	            }
	        }
	    }

	}



    //return rad,sol,divergencePoint
    char filename[500];
	sprintf(filename, "%s_rad.dat",outputFile);
	FILE *frad = fopen(filename,"wb");
	if (frad!=NULL)
	{
		fwrite(rad,sizeof(double),totalLength,frad);
		fclose(frad);
	}
	else{
		printf("Error occured when trying to open file (frad) to save wavefunction!\n");
		return 1;
	}


	sprintf(filename, "%s_sol.dat",outputFile);
	FILE *fsol = fopen(filename,"wb");
	if (fsol!=NULL)
	{
		fwrite(sol,sizeof(double),totalLength,fsol);
		fclose(fsol);
	}
	else{
		printf("Error occured when trying to open file (fsol) to save wavefunction!\n");
		return 1;
	}


	sprintf(filename, "%s_divergnece.dat",outputFile);
	FILE *fdiv = fopen(filename,"wb");
	if (fdiv!=NULL)
	{
		fwrite(&divergencePoint,sizeof(int),1,fdiv);
		fclose(fdiv);
	}
	else{
		printf("Error occured when trying to open file (fdiv) to save wavefunction!\n");
		return 1;
	}

    free(rad); // deallocate memory block
    free(sol);
    return 0;
}
