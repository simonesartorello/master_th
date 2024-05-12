////////////////////////////////////////////////////////////
// sampler.cpp
////////////////////////////////////////////////////////////

#include <iostream> 
#include <cfloat>   // constants related to 
                    // floating-point types (DBL_MAX)
#include <cmath>    
#include <math.h>    
#include <iomanip>  
#include <ctime>     
#include <fstream>  
#include <random>   


////////////////////////////////////////////////////////////
// SETTING PROBLEM PARAMETERS
////////////////////////////////////////////////////////////
// set number of test particles to sample
#define SAMPLES   10000
////////////////////////////////////////////////////////////


using namespace std;

////////////////////////////////////////////////////////////
// MAIN
////////////////////////////////////////////////////////////
int main(){

  double x,y,z;                     
  double vx,vy,vz;                  
  double theta, phi, phi1, theta1;  
  double radius, velocity;     
  // variables to store random numbers generated
  double s1,s2,s3,s6,s7;  
  double E0, rr;                   
  int n;   
  int samples = SAMPLES;            


// Random number generation setup
////////////////////////////////////////////////////////////
  //random_device rd;
  //mt19937 gen(rd());

  // Mersenne Twister random number generator
  mt19937 gen(7);      

  uniform_real_distribution<double> dis1(0.0, nextafter(1.0, 
                                          DBL_MAX));
  // uniform_real_distribution does [start, stop) but we 
  // want [start, stop] nextafter passes the next largest 
  // value instead
  uniform_real_distribution<double> dis2(0.0,1.0);
  uniform_real_distribution<double> dis3(0.0,1.0);
  uniform_real_distribution<double> dis4(0.0,1.0);
  uniform_real_distribution<double> dis5(0.0,0.1);
  uniform_real_distribution<double> dis6(0.0, nextafter(1.0, 
                                          DBL_MAX));
  uniform_real_distribution<double> dis7(0.0,1.0);


////////////////////////////////////////////////////////////
//  sampling random test particles 
////////////////////////////////////////////////////////////
  
  ofstream data1;
  data1.open("samples.txt");
  cout << setiosflags(ios::scientific);
  data1 << setiosflags(ios::scientific); 

  n = 0;

  while (n < samples + 1){

// populating position space
////////////////////////////////////////////////////////////
    s1 = dis1(gen);             // generate random s1 
    s2 = dis2(gen); 
    s3 = dis3(gen);
    // randomizing angular coordinates 
    phi = s2*2.0*M_PI;          // [0,2PI)
    theta = acos(2.0*s1-1);     // [-1,+1]
    // radius computed using the  inverse cumulative 
    // distribution function of the Plummer model with a = 1
    radius = 1.0/sqrt((pow(s3,-2.0/3.0)-1.0));   

    // cartesian coordinates
    x = radius * cos(phi) * sin(theta);
    y = radius * sin(phi) * sin(theta);
    z = radius * cos(theta);
////////////////////////////////////////////////////////////



// rejection sampling to populate velocity space 
////////////////////////////////////////////////////////////
    double q = 0.0;
    double g = 0.1;          // upper limit target dist
    while (g > q*q*pow((1.0-q*q),3.5)){
      q = dis4(gen);
      g = dis5(gen);
    }
////////////////////////////////////////////////////////////
    // computing velocity magnitude
    velocity = q * sqrt(2.0) * pow((1.0 + radius*radius),
                                    (-0.25));
    s6 = dis6(gen);    
    s7 = dis7(gen);
    // randomizing angular coordinates 
    phi1 = s7*2.0*M_PI;         // [0,2PI)
    theta1 = acos(2.0*s6-1);    // [-1,+1]

    // cartesian coordinates
    vx = velocity * cos(phi1) * sin(theta1);
    vy = velocity * sin(phi1) * sin(theta1);
    vz = velocity * cos(theta1);
////////////////////////////////////////////////////////////



    // computing radius and energy for each sampled 
    // test particle
    rr = sqrt(x*x + y*y + z*z + 1.0);
    E0 = 0.5*(vx*vx + vy*vy + vz*vz) - 1.0/rr;

    // check energy and radius values and saving them 
    // in a file: we aim for a constrained energy range and 
    // impose a maximum allowable radius to ensure more 
    // controlled initial conditions
    if (E0 >= -1.0 && E0 < 0.0 && rr < 13){
      n += 1;           // increment accepted values counter
      data1  << x << "   "  << y << "   "  << z << "   "  
            << vx << "   " << vy << "   " << vz << endl;
    } 

  }
  
  data1.close();
}