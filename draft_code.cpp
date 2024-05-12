////////////////////////////////////////////////////////////
// tesi.cpp
////////////////////////////////////////////////////////////
/* Code structure
  1 SETTING PROBLEM PARAMETERS
  2 CHOOSING WHICH MODULE TO RUN
  3 MAIN
  4 Random number generation setup
  5 loading test particles sampled file 
  6 NUMERICAL METHOD
  6.1 starting time iterations
  7 PERTURBATION 
  8 NOISE
  8.1 fractional velocity volume function 
////////////////////////////////////////////////////////////
  NOTE:
  when i had to run a high number 
  of simulations i modify the code so that i could give as 
  input different arguments to automate it as mush as 
  possible with a bash code, so that i could speed up the 
  process and reduce possible errors. 

  int main(int argc, char* argv[]){

    if (argc < number of arguments) {

    // for the case of oscillations around eq or 
    // contraction (remember to change argc < 3)
    // cerr << "Usage: ./your_program <perturbation_factor> 
                            <perturbation_freq>" << endl;

    // for the computation of lyapunov exponents during 
    // damped oscillations (remember to change argc < 6)
    // factor radius to choose the radius "length"
    // seed  to choose seed value
    // trac to choose the test particle to simulate  
    //          for k times with different noise values
    // shift to choose epsilon per lyapunov
    // cerr << "Usage: ./your_program <factor_radius> 
                  <seed> <trac> <numb> <shift>" << endl;
    return 1;
  }

  // stod converts cline argument to a double
  // double perturbation_factor = stod(argv[1]); 
  // double perturbation_freq   = stod(argv[2]);


  // when running code "for lyapunov"
  // stod converts first com. line argument to a double
  // double factor_radius = stod(argv[1]); 
  // stoi converts first command-line argument to a int
  // int seed   = stoi(argv[2]);
  // int trac = stoi(argv[3]); 
  // int numb   = stoi(argv[4]);
  // double shift   = stod(argv[5]);
*/
////////////////////////////////////////////////////////////


#include <iostream>   // input/output operations
#include <cfloat>     // for DBL_MAX
#include <cmath>
#include <math.h>
#include <iomanip>
#include <ctime>
#include <fstream>    // file handling
#include <random>


////////////////////////////////////////////////////////////
// 1 SETTING PROBLEM PARAMETERS
////////////////////////////////////////////////////////////
#define SAMPLES   10000            // N_test particles
#define PARTICLES 1000000000000    // N_tot system 
#define INCREMENT 0.01             // dt 
#define MAXSTEPS  1000             // integration steps 

#define SHIFT     0.0
#define FACTOR_RADIUS 1.0

#define CUTOFF_RADIUS 5.0          // cutoff radius
////////////////////////////////////////////////////////////

#define PERTURBATION_FACTOR 0.05   // perturb factor
#define PERTURBATION_FREQ   3.0    // perturb frequency
#define CUTOFF_RADIUS       5.0    // cutoff radius


////////////////////////////////////////////////////////////
// 2 CHOOSING WHICH MODULE TO RUN
////////////////////////////////////////////////////////////
#define NUMERICALMETHOD 1   // 0 = verlet position
                            // 1 = mannella SLO 

// PERTURBATION = choose if oscillation (which one) or not
#define PERTURBATION    2   // 0 = NO
                            // 1 = oscillations around eq
                            // 2 = exp damped oscillations 

// NOISE = choose of DYNAMICAL FRICTION+NOISE effects or not
#define NOISE           1   // IF NUMERICALMETHOD 1 
                            // noised no = 0 or yes = 1

// fractional velocity volume function PSI
#define PSI             0     
// PSI = 0 if psi = erf(v/sigma*sqrt(2)) -
// - 2v*exp(-v^2 / 2sigma^2) / sigma*sqrt(2PI);
// PSI = 1 if psi = 1;
////////////////////////////////////////////////////////////



using namespace std;

////////////////////////////////////////////////////////////
// 3 MAIN
////////////////////////////////////////////////////////////
int main(){

  // neq = number of equations
  int neq = 3;          
  int nsteps, n; 
  // variable to label test particles 
  int k;  


  double shift = SHIFT; 
  double factor_radius = FACTOR_RADIUS;

  double t0 = 0.0;
  //set increment
  double t, dt = INCREMENT;              
  //set number of steps
  int maxsteps = MAXSTEPS;               
  // number of samples
  int samples = SAMPLES;
  double perturbation_factor = PERTURBATION_FACTOR;
  double perturbation_freq   = PERTURBATION_FREQ;
  double cutoff_radius       = CUTOFF_RADIUS;
  double E, E0, E_mean, Ekin, Epot, Etot, E2, Ekin2, Epot2;
  double xx[samples][neq],vv[samples][neq];
  double particles = PARTICLES;



  double force[samples][neq];
  double rr, x,y,z, rr2, vx ,vy, vz, vel;
  //scale lenght a, A = a*a 
  double a, A;                     
  double squared_x, squared_y, squared_z; 
  double rmsx, rmsy, rmsz; 
  double squared_vx, squared_vy, squared_vz; 
  double rmsvx, rmsvy, rmsvz; 
  double squared_r, squared_vr; 
  double rmsr, rmsvr, radial_velocity;
  double tau, c1,c2,d1;
  double eta = 0.0;
  double G = 1.0, M = 1.0; 
  double sigma3, rho, coulomb, numberdensity;
  // velocity dispersion sigma
  double sigma;                    
  // fractional velocity volume function psi
  double psi;                    
  // zeta = K + U of i particle 
  double zeta;                     
  // sqrt(a1*a1 + a2*a2 + a3*a3)
  double ddd;                      


  //for emittance computation
  double sum_x, sum_y, sum_z, sum_vx, sum_vy, sum_vz;
  double x_mean, y_mean, z_mean, 
          vx_mean, vy_mean, vz_mean;
  double r_mean, vr_mean;
  double difference_x, difference_y, difference_z, 
          difference_vx, difference_vy, difference_vz;
  double sumSquaredDiff_x, sumSquaredDiff_y, 
          sumSquaredDiff_z, sumSquaredDiff_vx, 
          sumSquaredDiff_vy, sumSquaredDiff_vz;
  double variance_x, variance_y, variance_z, 
          variance_vx, variance_vy, variance_vz;
  double sumCov_x_vx, sumCov_y_vy, sumCov_z_vz;
  double covar_x_vx, covar_y_vy, covar_z_vz;
  double emit_x, emit_y,emit_z,emit;

  double displacement;
  double r_t0;
  double displacement2;
  double displacement2SUM;
  double MSD;


  double Efin[samples];
  // test particles lost (if E>=0 and r>=r_cutoff) counter
  double LOSTcounter;                     
  double freqLOST;

  //variables used when loading initial conditions file
  int ROWS = SAMPLES;      // number of rows in the file
  int COLS = 6;            // number of columns in the file 


////////////////////////////////////////////////////////////
// 4 Random number generation setup
////////////////////////////////////////////////////////////
  // Mersenne Twister random number generator
  mt19937 gen(1); 
  // when computing lyapunov exp, change the 
  // seed for each run    
  // mt19937 gen(seed);



  double a1, a2, a3;
  normal_distribution<double> normal_a1(0.0,1.0);
  normal_distribution<double> normal_a2(0.0,1.0);
  normal_distribution<double> normal_a3(0.0,1.0);

  double b1, b2, phi_b1, theta_b2;
  uniform_real_distribution<double> random_b1(0.0, 
                                  nextafter(1.0, DBL_MAX));
  uniform_real_distribution<double> random_b2(0.0,1.0);

  

////////////////////////////////////////////////////////////
// 5 loading test particles sampled file 
////////////////////////////////////////////////////////////
  cout << setiosflags(ios::scientific);
  // array to store loaded data
  double data[ROWS][COLS];  
  
  // input sampler file
  std::ifstream inputFile("samples.txt"); 
  // check for errors
  if (!inputFile) {
      std::cout << "failed to open file" 
                << std::endl;
      return 1; 
  }
  // iterate over each element of the array and read 
  // corresponding data from the file using the >> operator
  for (int i = 0; i < ROWS; i++) {
      for (int j = 0; j < COLS; j++) {
          inputFile >> setiosflags(ios::scientific);
          inputFile >> data[i][j];
      }
  }
  inputFile.close();

  // saving contents of the data array in xx and vv 
  // by iterating over each element 
  for (int i = 0; i < ROWS; i++) {
    for (int j = 0, k = 0; j < 3 && k < 3; j++, k++) {
      xx[i][k] = factor_radius * data[i][j] + shift;
      // + when not computing lyap, to have a cluster 
      // set factor_radius = 1.0e-4 (example) 
      // + shift = 0.00001 for lyap computation 
      // (it's the difference in the init.conditions)
      // + data[i][j] with fixed i = tracer and choose it 
      // every time to run it more times with different 
      // noise values, changing random generator
      // seed everytime in this case
      }
  }
  for (int i = 0; i < ROWS; i++) {
    for (int j = 3, k = 0; j < 6 && k < 3; j++, k++) {
      // vv[i][k] = data[i][j]; 
      
      vv[i][k] = 0.0; // if contraction case with or 
      // without lyap computation
      }
  }



////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////
// 6 NUMERICAL METHOD
////////////////////////////////////////////////////////////

  ofstream fraction;
  ofstream coordinates;
  ofstream energy0;
  ofstream energy1;
  ofstream energy2;
  ofstream rootmeansq;

  // fraction of escapers
  fraction.open("frac_esc.txt");
  // initial energy distribution
  energy0.open("initial_dist.txt");
  // final energy distribution
  energy1.open("final_dist.txt");
  // energy distribution at intermediate steps
  energy2.open("energyFRACTIONtime.txt");
  // coordinates xx, vv during evolution
  coordinates.open("coordinates.txt");

  // file where i save t, positions, velocities, emittance
  // and all the other useful quantities;
  rootmeansq.open("rootmeansq.txt");

  cout << setiosflags(ios::scientific);
  fraction << setiosflags(ios::scientific);
  energy0 << setiosflags(ios::scientific);
  energy1 << setiosflags(ios::scientific);
  energy2 << setiosflags(ios::scientific);
  coordinates << setiosflags(ios::scientific);
  rootmeansq << setiosflags(ios::scientific);



// cycle to save initial values of energy and radius of 
// each test particle at t0
////////////////////////////////////////////////////////////
  for (int k = 0; k < samples; k++){
    //E0 = 0.0;
    rr2 = sqrt(xx[k][0]*xx[k][0] + xx[k][1]*xx[k][1] + 
                xx[k][2]*xx[k][2] + 1.0);
    E0 = 0.5*(vv[k][0]*vv[k][0] + vv[k][1]*vv[k][1] + 
                vv[k][2]*vv[k][2]) - 1.0/rr2;
    // Etot += E0;
    energy0 << k << ' ' << E0 << ' ' << rr2 << endl;
  }
  // Etot/samples to check mean energy of the test 
  // particles cluster I am using for the simulation 
  //cout <<  Etot/samples << endl;

// flag to check the first time number of 
// escapers >= 50% total samples, initially set to true
bool firstTimeHalfReached = true;


// 6.1 starting time iterations
////////////////////////////////////////////////////////////
  for (int i = 1; i < maxsteps+1; i++){
    // cycle over timesteps

    //initializing variables to zero
    squared_x  = 0.0;
    squared_y  = 0.0;
    squared_z  = 0.0;
    squared_vx = 0.0;
    squared_vy = 0.0;
    squared_vz = 0.0;
    squared_r  = 0.0;
    squared_vr = 0.0;
    sum_x      = 0.0;
    sum_y      = 0.0;
    sum_z      = 0.0;
    sum_vx     = 0.0;
    sum_vy     = 0.0;
    sum_vz     = 0.0;
    rr2        = 0.0;

    Ekin = 0.0;
    Epot = 0.0;

    LOSTcounter = 0.0;


    for (k = 0; k < samples; k++){
      //cycle over each test particle
      t = t0 + dt*i;
      r_t0 = sqrt(xx[k][0]*xx[k][0] + xx[k][1]*xx[k][1] 
                + xx[k][2]*xx[k][2] + 1.0);

////////////////////////////////////////////////////////////
    #if NUMERICALMETHOD == 0     //Verlet position
////////////////////////////////////////////////////////////
      for (n=0; n < neq; n++){
        xx[k][n] = xx[k][n]+0.5*dt*vv[k][n];
      }
    
////////////////////////////////////////////////////////////
// PERTURBATION 
////////////////////////////////////////////////////////////
      #if PERTURBATION == 0
      a = 1.0;
      A = a*a;
      #elif PERTURBATION == 1
      a = 1.0 + perturbation_factor*
                sin(perturbation_freq*t);
      A = a*a;
      #endif
      //cout << a << endl;
      x = xx[k][0];
      y = xx[k][1];
      z = xx[k][2];

      rr = sqrt(x*x + y*y + z*z + A);
////////////////////////////////////////////////////////////

      // computing force along x, y, z axis for each tracer
      force[k][0] = -G*M*x/(rr*rr*rr);
      force[k][1] = -G*M*y/(rr*rr*rr);
      force[k][2] = -G*M*z/(rr*rr*rr);
    }



    for (k=0; k < samples; k++){
      for (n=0; n<neq; n++){
        //drift
        vv[k][n] = vv[k][n] +dt*force[k][n];
      }

      for (n=0; n<neq; n++){
        //kick
        xx[k][n] = xx[k][n] + 0.5*dt*vv[k][n];
      }
////////////////////////////////////////////////////////////
    #elif NUMERICALMETHOD == 1 
    // Mannella Symplectic Low Order (SLO)
////////////////////////////////////////////////////////////
      for (n=0; n < neq; n++){     
        //evolve position by half a step (drift)
        xx[k][n] = xx[k][n] + 0.5*dt*vv[k][n]; 
      }
      

////////////////////////////////////////////////////////////
// 7 PERTURBATION 
////////////////////////////////////////////////////////////
      #if PERTURBATION == 0
      a = 1.0;
      A = a*a;
      #elif PERTURBATION == 1
      a = 1.0 + perturbation_factor*
                sin(perturbation_freq*t);
      A = a*a;
      #elif PERTURBATION == 2
      a = 1.0 + exp(-t*t/2)*cos(perturbation_freq*t);
      A = a*a;  
      #endif
      
      x = xx[k][0];
      y = xx[k][1];
      z = xx[k][2];

      // r^2 = x*x + y*y + z*z
      // a^2 = A = a*a
      // rr = sqrt(r^2 + a^2)
      rr  = sqrt(x*x + y*y + z*z + A); 
      vel = sqrt(vv[k][0]*vv[k][0] + vv[k][1]*vv[k][1] 
            + vv[k][2]*vv[k][2]);


////////////////////////////////////////////////////////////
// 8 NOISE 
////////////////////////////////////////////////////////////
      #if NOISE == 1

      rho = ((3.0/(4.0*M_PI)) * M * A ) / (pow(rr*rr,2.5));
      coulomb = 10.0; // coulomb logarithm
      numberdensity = particles/((4.0/3.0)*M_PI*rr*rr*rr);
      // tau = mean interparticle distance
      tau = pow(numberdensity,(-1.0/3.0)); 

////////////////////////////////////////////////////////////
// 8.1 fractional velocity volume function form choice
////////////////////////////////////////////////////////////
      // fractional velocity volume function psi(vel) = 1.0 
      #if PSI == 1
      // cube of sqrt of velocity dispersion  
      sigma3 = pow(G*M/(6.0*rr),1.5);     
      // dynamical friction coefficient 
      eta = 4.0*M_PI*G*G*(2.0 * 1.0/particles)
                          *rho*coulomb/sigma3;        

      // fractional velocity volume function psi(vel) 
      // with the "right" functional form
      #elif PSI == 0
      // sqrt of velocity dispersion 
      sigma = sqrt(G*M/(6.0*rr)); 
      // fractional velocity volume function psi(vel)
        if (vel == 0.0){
          psi = erf(vel/(sigma*sqrt(2))) - 
                2*vel*exp(-vel*vel/(2*sigma*sigma))/
                (sigma*sqrt(2*M_PI));
          eta = 0.0;
        }
        else{
          psi = erf(vel/(sigma*sqrt(2))) - 
          2*vel*exp(-vel*vel/(2*sigma*sigma))/
          (sigma*sqrt(2*M_PI));
          // dynamical friction coefficient 
          eta = 4.0*M_PI*G*G*(2.0 * 1.0/particles)
                            *rho*coulomb*psi/(vel*vel*vel);  
        }
      #endif
////////////////////////////////////////////////////////////


      // sampling 3 different varaibles from a gaussian 
      // with mean = 0.0, variance = 1.0
      a1 = normal_a1(gen);
      a2 = normal_a2(gen);
      a3 = normal_a3(gen);
      ddd = sqrt(a1*a1 +a2*a2 +a3*a3);

      // energy of each particle used for the computation
      // of noise  
      zeta = abs(0.5 * (vv[k][0]*vv[k][0] 
                      + vv[k][1]*vv[k][1] 
                      + vv[k][2]*vv[k][2]) - 1.0/rr);
      zeta = zeta * ddd;

      // mannella SLO coefficients
      c1 = 1.0 - eta*dt*0.5;
      c2 = 1.0/(1.0 + eta*dt*0.5);
      d1 = sqrt(2.0*zeta*eta*dt);

      // projection of noise component randomizing 
      // angular coordinates 
      b1 = random_b1(gen); // generate b1 in range [0,1)
      b2 = random_b2(gen); 
      phi_b1   = b1*2.0*M_PI;
      theta_b2 = acos(2.0*b2-1.0);

////////////////////////////////////////////////////////////
      #elif NOISE == 0
      c1 = 1;
      c2 = 1;
      d1 = 0;
      #endif
////////////////////////////////////////////////////////////

      force[k][0] = -G*M*x/(rr*rr*rr);
      force[k][1] = -G*M*y/(rr*rr*rr);
      force[k][2] = -G*M*z/(rr*rr*rr);
    }

    for (k = 0; k < samples; k++){
      //evolve velocities by a full step (kick)
      vv[k][0] = c2 * (c1 * vv[k][0] + dt * force[k][0] 
                    + d1 * cos(phi_b1) * sin(theta_b2));
      vv[k][1] = c2 * (c1* vv[k][1] + dt*force[k][1]
                    + d1 * sin(phi_b1) * sin(theta_b2));
      vv[k][2] = c2 * (c1 * vv[k][2] + dt * force[k][2]
                    + d1*cos(theta_b2));


      for (n=0; n<neq; n++){   
        //evolve position by half a step (drift)
        xx[k][n] = xx[k][n] + 0.5*dt*vv[k][n];
      }

////////////////////////////////////////////////////////////
    #endif
////////////////////////////////////////////////////////////

      // saving position and velocity coordinates 
      if (i%100 == 0 ){
        coordinates << xx[k][0] << "   " << xx[k][1] 
           << "   " << xx[k][2] << "   " << vv[k][0] 
           << "   " << vv[k][1] << "   " << vv[k][2] 
           << endl;
      }

      // for the computation of the MSD
      displacement  = rr2 - r_t0;
      displacement2 = displacement * displacement;

      displacement2SUM += displacement2;
      
      // sums of squared values of position and velocity
      squared_x  += xx[k][0]*xx[k][0];
      squared_y  += xx[k][1]*xx[k][1];
      squared_z  += xx[k][2]*xx[k][2];
      squared_vx += vv[k][0]*vv[k][0];
      squared_vy += vv[k][1]*vv[k][1];
      squared_vz += vv[k][2]*vv[k][2];

      // sums values of position and velocity
      sum_x  += xx[k][0];
      sum_y  += xx[k][1];
      sum_z  += xx[k][2];
      sum_vx += vv[k][0];
      sum_vy += vv[k][1];
      sum_vz += vv[k][2];

      // radial velocity computation
      radial_velocity = (xx[k][0]*vv[k][0] +
                         xx[k][1]*vv[k][1] + 
                         xx[k][2]*vv[k][2]) / 
                         (sqrt(xx[k][0]*xx[k][0] + 
                               xx[k][1]*xx[k][1] + 
                               xx[k][2]*xx[k][2]));

      // sums of squared radius and radial velocity
      squared_r  += (xx[k][0]*xx[k][0] + xx[k][1]*xx[k][1] 
                    + xx[k][2]*xx[k][2]);
      squared_vr += radial_velocity*radial_velocity;

      rr2 = sqrt(xx[k][0]*xx[k][0] + xx[k][1]*xx[k][1] 
               + xx[k][2]*xx[k][2] + A);

      // energy computation
      Ekin = 0.5*(vv[k][0]*vv[k][0] + vv[k][1]*vv[k][1] 
                + vv[k][2]*vv[k][2]);
      Epot = - 1.0/(rr2);

      Ekin2 += 0.5*(vv[k][0]*vv[k][0] + vv[k][1]*vv[k][1] 
                + vv[k][2]*vv[k][2]);
      Epot2 += - 1.0/(rr2);


      E  = Ekin + Epot;
      E2 = Ekin2 + Epot2;

      // saving energy value and radius for each particle
      if (i%100 == 0){
        energy2 << k << ' ' << E << ' ' << rr2 << endl;
      }


////////////////////////////////////////////////////////////
      // check energy to count test particles lost 
      if ( E >= 0.0 && rr2 >= cutoff_radius){
        LOSTcounter += 1 ;
      }
////////////////////////////////////////////////////////////      
    }

    // counting N(escapers)/N_tot vs time
    if (i%10 == 0){
      fraction << t << "  " << LOSTcounter/samples <<  endl; 
    }


    if (i == 3000000){
      // open file append mode to save all this values
      // in the same file
      ofstream esc_fixed_time("Nesc_t30000.txt",
                              ios::app);

      if (!esc_fixed_time.is_open()) {
        cerr << "error opening file" << endl;
        return 1;
      }

      esc_fixed_time << perturbation_freq << "  " 
                     << LOSTcounter/samples << "  " 
                     << t << "  " << perturbation_factor 
                     << endl;
    }

    if (LOSTcounter >= SAMPLES/2 && firstTimeHalfReached){
      ofstream half_escapers("t_half_esc.txt", ios::app);

      if (!half_escapers.is_open()) {
        cerr << "error opening file" << endl;
        return 1;
      }

      half_escapers << perturbation_freq << "  " 
                    << LOSTcounter/samples << "  " 
                    << t << "  " << perturbation_factor 
                    << endl;

      // set the flag to false so that this block won't 
      // execute again after half
      // number of initial test particles has gone
      firstTimeHalfReached = false;
    }

    // computing rms
    rmsx  = sqrt(squared_x/samples);
    rmsy  = sqrt(squared_y/samples);
    rmsz  = sqrt(squared_z/samples); 
    rmsvx = sqrt(squared_vx/samples);
    rmsvy = sqrt(squared_vy/samples);
    rmsvz = sqrt(squared_vz/samples); 
    rmsr  = sqrt(squared_r/samples);
    rmsvr = sqrt(squared_vr/samples);

    // computing mean energy
    E_mean = E2/samples;

    // computing mean values of x,y,z,vx,vy,vz
    x_mean  = sum_x/samples;
    y_mean  = sum_y/samples;
    z_mean  = sum_z/samples;
    vx_mean = sum_vx/samples;
    vy_mean = sum_vy/samples;
    vz_mean = sum_vz/samples;

    r_mean = x_mean*x_mean + y_mean*y_mean 
            + z_mean*z_mean;
    vr_mean = x_mean*vx_mean + y_mean*vy_mean 
            + z_mean*vz_mean / r_mean;

    // reinitializing variables to zero
    sumSquaredDiff_x  = 0.0;
    sumSquaredDiff_y  = 0.0;
    sumSquaredDiff_z  = 0.0;
    sumSquaredDiff_vx = 0.0;
    sumSquaredDiff_vy = 0.0; 
    sumSquaredDiff_vz = 0.0;
    sumCov_x_vx       = 0.0;
    sumCov_y_vy       = 0.0;
    sumCov_z_vz       = 0.0;


    for (int k = 0; k < samples; k++){
      // computing k-variable - mean value of the set 
      // k-variable belongs to
      difference_x  = xx[k][0] - x_mean;
      difference_y  = xx[k][1] - y_mean;
      difference_z  = xx[k][2] - z_mean;
      difference_vx = vv[k][0] - vx_mean;
      difference_vy = vv[k][1] - vy_mean;
      difference_vz = vv[k][2] - vz_mean;

      // computing sum of squared differences of step 
      // above (for variance)
      sumSquaredDiff_x  += difference_x * difference_x;
      sumSquaredDiff_y  += difference_y * difference_y;
      sumSquaredDiff_z  += difference_z * difference_z;
      sumSquaredDiff_vx += difference_vx * difference_vx;
      sumSquaredDiff_vy += difference_vy * difference_vy;
      sumSquaredDiff_vz += difference_vz * difference_vz;

      // computing sum of squared differences 
      // (for covariance)
      sumCov_x_vx += difference_x * difference_vx;
      sumCov_y_vy += difference_y * difference_vy;
      sumCov_z_vz += difference_z * difference_vz;
    } 

    // computing variance of x,y,z,vx,vy,vz values
    variance_x  = sumSquaredDiff_x/samples;
    variance_y  = sumSquaredDiff_y/samples;
    variance_z  = sumSquaredDiff_z/samples;
    variance_vx = sumSquaredDiff_vx/samples;
    variance_vy = sumSquaredDiff_vy/samples;
    variance_vz = sumSquaredDiff_vz/samples;

    // computing covariance of (x,vx) (y,vy) (z,vz) values
    covar_x_vx = sumCov_x_vx/samples;
    covar_y_vy = sumCov_y_vy/samples;
    covar_z_vz = sumCov_z_vz/samples;


    // computing emittance 
    emit_x = sqrt(variance_x * variance_vx - covar_x_vx 
                                           * covar_x_vx);
    emit_y = sqrt(variance_y * variance_vy - covar_y_vy 
                                           * covar_y_vy);
    emit_z = sqrt(variance_z * variance_vz - covar_z_vz 
                                           * covar_z_vz);
    emit = pow((emit_x*emit_y*emit_z),1.0/3.0);


    // mean square displacement
    MSD = displacement2SUM/samples;


    if (i%100 == 0){
      rootmeansq << t << "  " << rmsx << "  " << rmsy 
                 << "  " << rmsz << "  " << rmsvx << "  " 
                 << rmsvy << "  " << rmsvz << "  " 
                 << E_mean << "  " 
                 << rr2 << "  " << a << "  " << rmsr 
                 << "  " << rmsvr << "  " << emit << "  "
                 << MSD << endl;
      
      // to have mean values only when necessary (i.e. for 
      // the mean trajectory for lyapunov 
      // exp computation) and avoid long output
      //rootmeansq << t << "  " << x_mean <<"  " << y_mean 
      //           << "  " << z_mean << "  " << vx_mean << 
      //           "  " << vy_mean << "  " << vz_mean << 
      //           "  " << E_mean << "  "
      //           << rr2 << "  " << a << "  " << rmsr 
      //           << "  " << rmsvr << "  " << emit << " " 
      //           << r_mean << "  " << vr_mean << endl;
    }

    
  }

  for (int k = 0; k < samples; k++){
    //E0 = 0.0;
    rr2 = sqrt(xx[k][0]*xx[k][0] + xx[k][1]*xx[k][1] 
                                 + xx[k][2]*xx[k][2] + A);
    E0 = 0.5 * (vv[k][0]*vv[k][0] + vv[k][1]*vv[k][1] 
              + vv[k][2]*vv[k][2]) - 1.0/rr2;
    // saving final values of energy and radius 
    // for each particle
    energy1 << k << ' ' << E0 << ' ' << rr2 << endl;
  }

  fraction.close();
  energy0.close();
  energy1.close();
  energy2.close();
  coordinates.close();
  rootmeansq.close();
  
  return 0;
}

// end main
////////////////////////////////////////////////////////////
