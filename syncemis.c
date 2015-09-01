#include <stdio.h>
#include <math.h>

double getp(double B, double omega, double n);

int main(void)

{ 
  //put in the file name to write to
  char filename[20] = "data.txt" ;
  

 //let's make some arrays with varying values of B, n, and omega

  double B_array[500]; //declaring the arrays
  double n_array[500]; 
  double omega_array[500]; 
  unsigned long arr_length =sizeof(B_array)/sizeof(double) ;
  double omega_init = 10e8;

  //let's calculate some values for reasonable estimates of the B and n
  
  double alpha = .1; //viscosity parameter
  double m6 = 4.0; // mass of sgr A* in units of 10^6 solar masses
  double m3 = 1.0; // accretion rate in 10^-3 Eddington
  double r = 10.0; // radius in Rschwarz

  //now let's get some rough estimates for values of density and B field based on these numbers
  
  double ncent =2.0e10 * pow(alpha, -1) * pow(m6, -1) * m3 * pow(r, -1.5); // g/cm^3
  ncent = ncent * (1.0/1000.0)* pow(10,6); // converting to kg/m^3 

  //but in our equation we want a number density of electrons, so we need to divide out the mass
  ncent = ncent / 9.11e-31; //ncent now a number density
    
  double Bcent= 2.07e4 * pow(alpha, -.5) * pow(m6, -.5) * pow(m3, .5) * pow(r, -5.0/4.0); //Gauss
  Bcent = Bcent /10.0e4; // converting to Tesla

  printf("the value of ncent is: %le nelecs/m^3\n", ncent);
  printf("the value of Bcent is: %f Tesla\n", Bcent);
 

  int i;
  for (i=0; i< arr_length;i++)
    {
      B_array[i]=Bcent - Bcent*i/500.0;
      n_array[i]=ncent - ncent*i/500.0;

      omega_array[i] = omega_init + omega_init*pow(2,i/10.0); //

 
    }

  
  //printf("the value of ptot is, %le\n",getp(Bcent,omega_array[0],ncent)); 
 
  //writing data to file

   FILE *fp;
   fp = fopen(filename,"w");

   for (i=0;i<arr_length;i++)
          {
	    double val =  getp(Bcent, omega_array[i], ncent);

	    //if (val==0)
	    // {
	    //	printf("ok well that's wrong");
	    //	  }

	      
	    	fprintf(fp, "%le , %le\n", omega_array[i], val);

	  }

   fclose(fp);




  return 0;
}

double getp(double B, double omega, double n)  
{
 
  //defining some constants

  double const q = 1.602e-19; // Coulombs
  double const m = 9.109e-31; // kg
  double const c= 2.998e8; // m/s 

  double p = 2.0; //power law index of electron distribution
  double gamma_factor = 1.038 * 2.127; // not general, specific to p=2

  double ptot = ((sqrt(3.0)*pow(q,3)*B*n)/(2*M_PI*m*pow(c,2)*(p+1)))*gamma_factor*pow(m*c*omega/(3*q*B),-(p-1)/2.0);
    
  return ptot;

 }
