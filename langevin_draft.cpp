////////////////////////////////////////////////////////////////////////////////
/* Code structure ()

*/
////////////////////////////////////////////////////////////////////////////////


#include <iostream> //for input/output operations
#include <cfloat>   // for DBL_MAX
#include <cmath>
#include <math.h>
#include <iomanip>
#include <ctime>
#include <fstream>  //for file handling
#include <random>

double dot_product(const double* v1, const double* v2, int size);
void normalize_vector(double* v, int size);


////////////////////////////////////////////////////////////////////////////////
// SETTING PROBLEM PARAMETERS
////////////////////////////////////////////////////////////////////////////////
#define SAMPLES   1               // number of times I "sample" the same tracer with different noise values: example, each tracer 5 times, then take mean etc...
#define TRACERS_CLUSTER 100       // size of the tracers cluster used to study the system
//PARTICLES 1000000000000         // number of particles in the system 
#define INCREMENT 0.02            // dt 
#define MAXSTEPS  50000           // number of step to perform integration
#define SHIFT 0.0                 // set = 1.0 to shift from the center and remember modify epsilons to have a point in phase space
                                  // or in the case of shifting initial condition

// perturbation a = PERTURBATION_FACTOR*sin(PERTURBATION_FREQ*t)
#define PERTURBATION_FACTOR 0.05      // perturbation factor
#define PERTURBATION_FREQ   3.0       // perturbation frequency
#define CUTOFF_RADIUS       5.0       // cutoff radius


////////////////////////////////////////////////////////////////////////////////
// CHOOSING WHICH PART TO RUN
////////////////////////////////////////////////////////////////////////////////
#define NUMERICALMETHOD 1     // 0 = verlet position, 1 = mann SLO 
#define PERTURBATION    1     // 0 = NO, 1 = YES
#define NOISE           1     // if NUMERICALMETHOD 1 --> choose between noised no = 0 or yes = 1

// PSI = fractional velocity volume function. 
#define PSI             0     // PSI = 0 if psi = erf(v/sigma*sqrt(2)) - 2v*exp(-v^2 / 2sigma^2) / sigma*sqrt(2PI), PSI = 1 if psi =1


////////////////////////////////////////////////////////////////////////////////



using namespace std;

///////////////////////////////////////////////////////////////////////////////
// main
///////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[]){

  if (argc < 6) {
    // error in case perturbation factor and frequency not provided
    cerr << "Usage: ./your_program <particles> <seed> <trac> <numb> <shift>" << endl;
    return 1;
  }

  // particles = number of background particles
  long long particles = stoll(argv[1]); // stod converts the first command-line argument to a double
  // seed = random number generator seed
  int seed   = stoi(argv[2]);
  // trac is the tracer number i'm "following"
  int trac = stoi(argv[3]); // stod converts the first command-line argument to a int
  // numb = just to name the output file associated to each tracer
  int numb   = stoi(argv[4]);
  // shiftlyap = value of the shift btw the two tracers to let me compute the lyap exp
  double shiftlyap   = stod(argv[5]);


  int neq = 3, nsteps, n, k;              // neq number of equations

  // Initialize tangent vectors and orthogonal matrix
  double Q[neq][neq]; // Assuming neq is the dimension of your system
  double v0[neq];

  double lambda[neq];
  for (int i = 0; i < neq; ++i) {
    lambda[i] = 0.0;
  }

  //variables used in sampling part
  double theta, phi, phi1, theta1;
  double radius, velocity, s1,s2,s3,s6,s7, x1;
  //double shiftwrtorigin = SHIFT; //position shift wrt origin
  double epsilon1 = 1.0, epsilon2 = 1.0;

  double t0 = 0.0;
  double t, dt = INCREMENT;               //set increment
  int maxsteps = MAXSTEPS;                //set number of steps
  int samples = SAMPLES;                  // number of samples
  int tracers_cluster = TRACERS_CLUSTER;
  double Y[neq];
  double E, E0, E0mean, E_mean, Ekin, Epot, Ew[tracers_cluster];
  double xx[tracers_cluster][neq],vv[tracers_cluster][neq];
  //double particles = PARTICLES;



  double force[tracers_cluster][neq];
  double rr, x,y,z, rr2, norm1, norm2, vx ,vy, vz, vel;
  double a, A;                            //scale lenght a, A = a*a 
  double squared_x, squared_y, squared_z, rmsx, rmsy, rmsz; 
  double squared_vx, squared_vy, squared_vz, rmsvx, rmsvy, rmsvz; 
  double squared_r, squared_vr, rmsr, rmsvr, radial_velocity;
  double tau, noise, c1,c2,d1;
  double gamma = 0.0; 
  double T = 0.0;                         //solo per mann SHO
  double eta = 0.0;
  double eta1 = 0.0, eta2 = 0.0;          //mann SHO
  double G = 1.0, M = 1.0, sigma3, rho, coulomb, numberdensity;
  double sigma;                           // velocity dispersion
  double psi;                             // fractional velocity volume function 
  double zeta;                            // K + U of i particle 
  double ddd;                             // sqrt(a1*a1 + a2*a2 + a3*a3)


  //for emittance computation
  double sum_x, sum_y, sum_z, sum_vx, sum_vy, sum_vz;
  double x_mean, y_mean, z_mean, vx_mean, vy_mean, vz_mean;
  double r_mean, vr_mean;
  double difference_x, difference_y, difference_z, 
          difference_vx, difference_vy, difference_vz;
  double sumSquaredDiff_x, sumSquaredDiff_y, sumSquaredDiff_z, 
          sumSquaredDiff_vx, sumSquaredDiff_vy, sumSquaredDiff_vz;
  double variance_x, variance_y, variance_z, 
          variance_vx, variance_vy, variance_vz;
  double sumCov_x_vx, sumCov_y_vy, sumCov_z_vz;
  double covar_x_vx, covar_y_vy, covar_z_vz;
  double emit_x, emit_y,emit_z,emit;


  double Efin[tracers_cluster];
  double LOSTcounter;                     // number of samples lost (E>=0 and r>=r_cutoff)
  double freqLOST;


  double x_mean0,y_mean0,z_mean0;
  double vx_mean0, vy_mean0, vz_mean0;
  double x_sum0,y_sum0,z_sum0;
  double vx_sum0,vy_sum0,vz_sum0;


  //variables used when loading initial conditions file
  int ROWS = TRACERS_CLUSTER;                // number of rows in the file
  int COLS = 6;                           // number of columns in the file 


// Random number generation setup
//////////////////////////////////////////////////////////////////
  //random_device rd;
  //mt19937 gen(rd());
  mt19937 gen(seed);                         // Mersenne Twister random number generator

  uniform_real_distribution<double> dis1(0.0, nextafter(1.0, DBL_MAX));
  // uniform_real_distribution does [start, stop), but we want [start, stop]
  //  nextafter passes the next largest value instead
  uniform_real_distribution<double> dis2(0.0,1.0);
  uniform_real_distribution<double> dis3(0.0,1.0);
  uniform_real_distribution<double> dis4(0.0,1.0);
  uniform_real_distribution<double> dis5(0.0,0.1);
  uniform_real_distribution<double> dis6(0.0, nextafter(1.0, DBL_MAX));
  uniform_real_distribution<double> dis7(0.0,1.0);
  normal_distribution<double> gaussian(0.0,1.0);


  double a1, a2, a3;
  normal_distribution<double> normal_a1(0.0,1.0);
  normal_distribution<double> normal_a2(0.0,1.0);
  normal_distribution<double> normal_a3(0.0,1.0);

  double b1, b2, phi_b1, theta_b2;
  uniform_real_distribution<double> random_b1(0.0, nextafter(1.0, DBL_MAX));
  uniform_real_distribution<double> random_b2(0.0,1.0);

  

///////////////////////////////////////////////////////////////////////////////
// loading initial_cond.dat file if already have one
///////////////////////////////////////////////////////////////////////////////


  cout << setiosflags(ios::scientific);
  double data[ROWS][COLS];  // array to store the data from the .dat file

  std::ifstream inputFile("/Users/amministratore/Desktop/paper_lyap/100cdn/initial_cond_E0.dat"); //create an object of the ifstream class and open the .dat file using the file name
  if (!inputFile) {
      std::cout << "Failed to open file with initial conditions. error in edoc17.cpp" << std::endl;
      return 1;             // exit the program if the file cannot be opened
  }
  // iterate over each element of the array and read the corresponding data from the file using the >> operator
  for (int i = 0; i < ROWS; i++) {
      for (int j = 0; j < COLS; j++) {
          inputFile >> setiosflags(ios::scientific);
          inputFile >> data[i][j];
      }
  }
  inputFile.close();

  // saving contents of the data array in xx and vv by iterating over each element 
  // here i create # samples of the same tracer chosen with trac input
  for (int i = 0; i < tracers_cluster; i++) {
    for (int j = 0, k = 0; j < 3 && k < 3; j++, k++) {
      xx[i][k] = data[i][j] + shiftlyap;    // change data[i] to data[tracer] to work with a single tracer only
      }
  }
  for (int i = 0; i < tracers_cluster; i++) {
    for (int j = 3, k = 0; j < 6 && k < 3; j++, k++) {
      vv[i][k] = 0.0; //data[0][j];
      }
  }



  ofstream data2;
  data2.open("samples_vnull.txt");
  cout << setiosflags(ios::scientific);
  data2 << setiosflags(ios::scientific); 


  //for (k=0; k<ROWS; k++){
  //    data2 << xx[k][0] << "   " << xx[k][1] << "   " << xx[k][2] << "   "   <<vv[k][0] << "   " << vv[k][1] << "   " << vv[k][2] << endl;
  //} 
  //for (k=0; k<samples; k++){
  //    cout << trac << ' ' << xx[k][0] << "   " << xx[k][1] << "   " << xx[k][2] << "   "   <<vv[k][0] << "   " << vv[k][1] << "   " << vv[k][2] << endl;
  //}

  data2.close();

///////////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////////
// NUMERICAL METHODS
///////////////////////////////////////////////////////////////////////////////

  //ofstream fraction;
  //ofstream coordinates;
  //ofstream energy0;
  //ofstream energy1;
  //ofstream energy2;
  ofstream rootmeansq;

  //fraction.open("frac_esc_contr1.txt");
  //energy0.open("initial_dist_contr1.txt");
  //energy1.open("final_dist_contr1.txt");
  //energy2.open("energyFRACTIONtime_contr1.txt");
  //coordinates.open("coordinates2.txt");
  rootmeansq.open("rootmeansq" + to_string(static_cast<int>(numb)) + ".txt");

  cout << setiosflags(ios::scientific);
  //fraction << setiosflags(ios::scientific);
  //energy0 << setiosflags(ios::scientific);
  //energy1 << setiosflags(ios::scientific);
  //energy2 << setiosflags(ios::scientific);
  //coordinates << setiosflags(ios::scientific);
  rootmeansq << setiosflags(ios::scientific);



// cycle to save initial values
///////////////////////////////////////////////////////////////////////////////
  for (int k = 0; k < tracers_cluster; k++){
    //E0 = 0.0;
    rr2 = sqrt(xx[k][0]*xx[k][0] + xx[k][1]*xx[k][1] + xx[k][2]*xx[k][2] + 1.0);
    //E0 += 0.5* (vv[k][0]*vv[k][0] + vv[k][1]*vv[k][1] + vv[k][2]*vv[k][2]) - (1.0/(rr2))*0.5;
    E0 = 0.5*(vv[k][0]*vv[k][0] + vv[k][1]*vv[k][1] + vv[k][2]*vv[k][2]) - 1.0/rr2;
    //energy0 << k << ' ' << E0 << ' ' << rr2 << endl;
    x_sum0 += xx[k][0];
    y_sum0 += xx[k][1];
    z_sum0 += xx[k][2];
    vx_sum0 += vv[k][0];
    vy_sum0 += vv[k][1];
    vz_sum0 += vv[k][2];
  }

  x_mean0 = x_sum0/tracers_cluster;
  y_mean0 = y_sum0/tracers_cluster;
  z_mean0 = z_sum0/tracers_cluster;
  vx_mean0 = vx_sum0/tracers_cluster;
  vy_mean0 = vy_sum0/tracers_cluster;
  vz_mean0 = vz_sum0/tracers_cluster;
  //E0mean= E0/samples;
//cout << E0mean << endl;

// define a flag to check the first time number of escapers >= 50% total samples
bool firstTimeHalfReached = true;


rootmeansq << t0 << "  " << x_mean0 <<"  " << y_mean0 << "  " << z_mean0 << "  " << vx_mean0 << "  " << vy_mean0 << "  " << vz_mean0 << "  "<< 0 << "  "<< 0 <<
      " " << 0 << " " << 0 << " " << 0 << " " << 0 << " " << 0 << " " << 0 << " " << 0 << endl;

// method
///////////////////////////////////////////////////////////////////////////////
  for (int i = 1; i < maxsteps+1; i++){

    //reinitializing variables to zero
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
    //E = 0.0;

    for (k = 0; k < tracers_cluster; k++){
      t = t0 + dt*i;
      LOSTcounter = 0.0;


      //cout << t0 << "  " << xx[k][0] << "  " << xx[k][1] << "  " << xx[k][2] 
      //   << "  " << vv[k][0] << "  " << vv[k][1] << "  " << vv[k][2]  <<endl;

///////////////////////////////////////////////////////////////////////////////
    #if NUMERICALMETHOD == 0     //Verlet position
///////////////////////////////////////////////////////////////////////////////
      for (n=0; n < neq; n++){
        xx[k][n] = xx[k][n]+0.5*dt*vv[k][n];
      }
    
///////////////////////////////////////////////////////////////////////////////
//choose if PERTURBATION or not
///////////////////////////////////////////////////////////////////////////////
      #if PERTURBATION == 0
      a = 1.0;
      A = a*a;
      #elif PERTURBATION == 1
      a = 1.0 + PERTURBATION_FACTOR*sin(PERTURBATION_FREQ*t);
      A = a*a; //Plummer potential has a sort of softening because it prevents radius goes to zero 
      #endif
      //cout << a << endl;
      x = xx[k][0];
      y = xx[k][1];
      z = xx[k][2];

      rr =sqrt(x*x + y*y + z*z + A);
///////////////////////////////////////////////////////////////////////////////

      force[k][0] = -G*M*x/(rr*rr*rr);
      force[k][1] = -G*M*y/(rr*rr*rr);
      force[k][2] = -G*M*z/(rr*rr*rr);
    }



    for (k=0; k < samples; k++){
      for (n=0; n<neq; n++){
        vv[k][n] = vv[k][n] +dt*force[k][n];
      }

      for (n=0; n<neq; n++){
        xx[k][n] = xx[k][n] + 0.5*dt*vv[k][n];
      }
///////////////////////////////////////////////////////////////////////////////
    #elif NUMERICALMETHOD == 1 // Mannella Symplectic Low Order (SLO)
///////////////////////////////////////////////////////////////////////////////

      //evolve position by half a step (drift)
      xx[k][0] = xx[k][0] + 0.5*dt*vv[k][0];
      xx[k][1] = xx[k][1] + 0.5*dt*vv[k][1];
      xx[k][2] = xx[k][2] + 0.5*dt*vv[k][2];
      
      

///////////////////////////////////////////////////////////////////////////////
//choose if PERTURBATION or not
///////////////////////////////////////////////////////////////////////////////
      #if PERTURBATION == 0
      a = 1.0;
      A = a*a;
      #elif PERTURBATION == 1
      a = 1.0 + exp(-t*t/2)*cos(PERTURBATION_FREQ*t);
      A = a*a; //Plummer potential has a sort of softening because it prevents radius goes to zero 
      #endif
      
      x = xx[k][0];
      y = xx[k][1];
      z = xx[k][2];

      // r^2 = x*x + y*y + z*z , a^2 = A = a*a
      rr = sqrt(x*x + y*y + z*z + A); // rr = sqrt(r^2 + a^2)
      vel = sqrt(vv[k][0]*vv[k][0] + vv[k][1]*vv[k][1] + vv[k][2]*vv[k][2]);
      //cout << vel << endl;


///////////////////////////////////////////////////////////////////////////////
//choose if NOISE or not
///////////////////////////////////////////////////////////////////////////////
      #if NOISE == 1

      rho = ((3.0/(4.0*M_PI)) * M * A ) / (pow(rr*rr,2.5));
      coulomb = 10.0; // coulomb logarithm
      noise = gaussian(gen);
      numberdensity = particles/((4.0/3.0)*M_PI*rr*rr*rr);
      tau = pow(numberdensity,(-1.0/3.0));  //mean interparticle distance
      //tau = (tgamma(1.0/3.0)/3.0)*pow(3.0/(4.0*M_PI*numberdensity),1.0/3.0);

      #if PSI == 1
      sigma3 = pow(G*M/(6.0*rr),1.5);       //cube of sqrt of velocity dispersion  
      // fractional velocity volume function psi(vel) = 1.0 
      eta = 4.0*M_PI*G*G*(2.0 * 1.0/particles)*rho*coulomb/sigma3;        // dynamical friction coefficient 

      #elif PSI == 0
      sigma = sqrt(G*M/(6.0*rr));           //sqrt of velocity dispersion 
      // fractional velocity volume function psi(vel)
        if (vel == 0.0){
          psi = erf(vel/(sigma*sqrt(2))) - 2*vel*exp(-vel*vel/(2*sigma*sigma))/(sigma*sqrt(2*M_PI));
          eta = 0.0;
        }
        else{
          psi = erf(vel/(sigma*sqrt(2))) - 2*vel*exp(-vel*vel/(2*sigma*sigma))/(sigma*sqrt(2*M_PI));
          eta = 4.0*M_PI*G*G*(2.0 * 1.0/particles)*rho*coulomb*psi/(vel*vel*vel);   // dynamical friction coefficient 
        }
      //cout << psi << ' ' << eta << endl;
      #endif


      // estracting 3 different number from a gaussian with mean = 0.0, variance =1.0
      a1 = normal_a1(gen);
      a2 = normal_a2(gen);
      a3 = normal_a3(gen);
      ddd = sqrt(a1*a1 +a2*a2 +a3*a3);

      zeta=  abs(0.5*(vv[k][0]*vv[k][0] + vv[k][1]*vv[k][1] + vv[k][2]*vv[k][2]) - 1.0/rr);
      zeta = zeta * ddd;

      //if (i == 10000){
      //  cout << a1 << "  " << a2 <<"  " << a3 << "  " << zeta << endl;
      //}


      c1 = 1.0 - eta*dt*0.5;
      c2 = 1.0/(1.0 + eta*dt*0.5);
      d1 = sqrt(2.0*zeta*eta*dt);


      b1 = random_b1(gen); // generate b1 in the range [0,1)
      b2 = random_b2(gen); // range [0,1)
      phi_b1 = b1*2.0*M_PI;
      theta_b2 = acos(2.0*b2-1.0);

      //cout << phi_b1 << ' ' << theta_b2 << ' ' << d1 << endl;

      //if (i == 10000){
      //  cout << c1 << "  " << c2 <<"  " << d1 << "  " << phi_b1 << ' ' << theta_b2 << ' ' << ' ' << zeta << endl;
      //}


      #elif NOISE == 0
      c1 = 1;
      c2 = 1;
      d1 =0;
      #endif
///////////////////////////////////////////////////////////////////////////////

      force[k][0] = -G*M*x/(rr*rr*rr);
      force[k][1] = -G*M*y/(rr*rr*rr);
      force[k][2] = -G*M*z/(rr*rr*rr);
    

      vv[k][0] = c2*(c1* vv[k][0] + dt*force[k][0] + d1*cos(phi_b1) * sin(theta_b2));
      vv[k][1] = c2*(c1* vv[k][1] + dt*force[k][1] + d1*sin(phi_b1) * sin(theta_b2));
      vv[k][2] = c2*(c1* vv[k][2] + dt*force[k][2] + d1*cos(theta_b2));

      //if (i == 1){
      //  cout << c1 << "  " << c2 <<"  " << d1 << "  " << phi_b1 << "  " << theta_b2 << "  " << zeta << endl;
      //}

      //for (n=0; n<neq; n++){        //evolve velocities by a full step (kick)
      //  vv[k][n] = c2*(c1* vv[k][n] + dt*force[k][n] + d1*noise);
      //}

      //evolve position by half a step (drift)

      //if (i == 10){
        //cout << vv[k][0] << "  " << vv[k][1] << "  " << vv[k][2] << endl;
      //}

      xx[k][0] = xx[k][0] + 0.5*dt*vv[k][0];
      xx[k][1] = xx[k][1] + 0.5*dt*vv[k][1];
      xx[k][2] = xx[k][2] + 0.5*dt*vv[k][2];
      //}

      //if (i == 5){
      //  cout << t << "  " << xx[k][0] <<"  " << xx[k][1] << "  " << xx[k][2] << "  " << vv[k][0] << "  " << vv[k][1] << "  " << vv[k][2] << endl;
      //}
      

/*
    // Perturb the state variables x, y, and z
    v0[0] = 1.0e-4 + xx[7][0]; // Perturb x
    v0[1] = 1.0e-4 + xx[7][1]; // Perturb y
    v0[2] = 1.0e-4 + xx[7][2]; // Perturb z


    // orthogonalization (Gram-Schmidt)
    for (int i = 0; i < neq; ++i) {
        for (int j = 0; j < neq; ++j) {
            v0[j] -= Q[i][j] * dot_product(v0, Q[i], neq);
        //cout << v0[j] << endl;
        }
    }
    
double magnitude = 0.0;

    //  magnitude of the vector
    for (int i = 0; i < neq; ++i) {
        magnitude += v0[i] * v0[i];
    }

    // avoid division by zero
    if (magnitude > 0.0) {
        magnitude = sqrt(magnitude);

        // normalize vector
        for (int i = 0; i < neq; ++i) {
            v0[i] /= magnitude;
        }
    }

    // Lyapunov exponent computation
    for (int i = 0; i < neq; ++i) {
        lambda[i] += log(fabs(dot_product(v0, Q[i], neq))) / dt;
    }*/
///////////////////////////////////////////////////////////////////////////////
    #endif
///////////////////////////////////////////////////////////////////////////////


      //if (i%100 == 0 ){
        //coordinates << xx[k][0] << "   " << xx[k][1] << "   " << xx[k][2] << "   "   <<vv[k][0] << "   " << vv[k][1] << "   " << vv[k][2] << endl;
        //coordinates << sqrt(xx[k][0]*xx[k][0] + xx[k][1]*xx[k][1] + xx[k][2]*xx[k][2]) << "   " << sqrt(vv[k][0]*vv[k][0] + vv[k][1]*vv[k][1] + vv[k][2]*vv[k][2]) << endl;
      //}
      //coordinates << xx[k][0] << "   " << xx[k][1] << "   " << xx[k][2] << "   "   <<vv[k][0] << "   " << vv[k][1] << "   " << vv[k][2] << endl;

      squared_x  += xx[k][0]*xx[k][0];
      squared_y  += xx[k][1]*xx[k][1];
      squared_z  += xx[k][2]*xx[k][2];
      squared_vx += vv[k][0]*vv[k][0];
      squared_vy += vv[k][1]*vv[k][1];
      squared_vz += vv[k][2]*vv[k][2];

      sum_x  += xx[k][0];
      sum_y  += xx[k][1];
      sum_z  += xx[k][2];
      sum_vx += vv[k][0];
      sum_vy += vv[k][1];
      sum_vz += vv[k][2];


      radial_velocity = (xx[k][0]*vv[k][0] + xx[k][1]*vv[k][1] + xx[k][2]*vv[k][2]) / (sqrt(xx[k][0]*xx[k][0] + xx[k][1]*xx[k][1] + xx[k][2]*xx[k][2]));

      squared_r += (xx[k][0]*xx[k][0] + xx[k][1]*xx[k][1] + xx[k][2]*xx[k][2]);
      squared_vr += radial_velocity*radial_velocity;

      rr2 = sqrt(xx[k][0]*xx[k][0] + xx[k][1]*xx[k][1] + xx[k][2]*xx[k][2] + A);

      //E = 0.0;
      Ekin = 0.5*(vv[k][0]*vv[k][0] + vv[k][1]*vv[k][1] + vv[k][2]*vv[k][2]);
      Epot = - 1.0/(rr2);
      E = Ekin + Epot;

      /*
      if (i%100 == 0){
        //energy2 << k << ' ' << E << ' ' << rr2 << endl;
      }
      //energy << E << endl;
      */


      //data1 << t << "  " << sqrt(xx[k][0]*xx[k][0] + xx[k][1]*xx[k][1] + xx[k][2]*xx[k][2]) << "  " << xx[k][0] << "   " << xx[k][1] << "  " << xx[k][2] << "  " 
      //      << E << "  " << endl;

///////////////////////////////////////////////////////////////////////////////
// check Energy to count number of tracers lost 
      //if ( E >= 0.0 && rr2 >= CUTOFF_RADIUS){
      //  LOSTcounter += 1 ;
      //}
      //if ( E < 0 ){
        //energy <<  E <<  endl;

      //}
///////////////////////////////////////////////////////////////////////////////

    //energy << t << " " << E <<  endl;

      
    }


    /*
    //energy  << Efin[k] <<  endl;
    if (i%10 == 0){
      //fraction << t << "  " << LOSTcounter/samples <<  endl; //to count N(escapers)/N_tot vs time
    } 
    */

/*
    if (i == 3000000){
      // open the file in append mode
      ofstream esc_fixed_time("Nesc_t30000VSomega_UNPERT.txt", ios::app);

      if (!esc_fixed_time.is_open()) {
        cerr << "Error opening file for appending!" << endl;
        return 1;
      }

      // write output to the file
      esc_fixed_time << PERTURBATION_FREQ << "  " << LOSTcounter/samples << "  " << t << endl;

      // file automatically close when outputFile goes out of scope
    }

    if (LOSTcounter >= SAMPLES/2 && firstTimeHalfReached){
      // open file in append mode
      ofstream half_escapers("t_half_escVSomega_UNPERT.txt", ios::app);

      if (!half_escapers.is_open()) {
        cerr << "Error opening file for appending!" << endl;
        return 1;
      }

      // write output to the file
      half_escapers << PERTURBATION_FREQ << "  " << t << "  " << LOSTcounter/samples << endl;

      // set flag to false so that this block won't execute again
      firstTimeHalfReached = false;
    }
*/

    //computing rms
    rmsx  = sqrt(squared_x/tracers_cluster);
    rmsy  = sqrt(squared_y/tracers_cluster);
    rmsz  = sqrt(squared_z/tracers_cluster); 
    rmsvx = sqrt(squared_vx/tracers_cluster);
    rmsvy = sqrt(squared_vy/tracers_cluster);
    rmsvz = sqrt(squared_vz/tracers_cluster); 
    rmsr  = sqrt(squared_r/tracers_cluster);
    rmsvr = sqrt(squared_vr/tracers_cluster);

    //computing mean energy
    E_mean = E/tracers_cluster;

    //computing mean values of x,y,z,vx,vy,vz
    x_mean  = sum_x/tracers_cluster;
    y_mean  = sum_y/tracers_cluster;
    z_mean  = sum_z/tracers_cluster;
    vx_mean = sum_vx/tracers_cluster;
    vy_mean = sum_vy/tracers_cluster;
    vz_mean = sum_vz/tracers_cluster;

    r_mean = x_mean*x_mean + y_mean*y_mean + z_mean*z_mean;
    vr_mean = x_mean*vx_mean + y_mean*vy_mean + z_mean*vz_mean / r_mean;

    //cout << sum_x << " " << sum_y << " " << sum_z << " " << sum_vx << " " << sum_vy << " " << sum_vz  << endl; 

    //reinitializing variables to zero
    sumSquaredDiff_x  = 0.0;
    sumSquaredDiff_y  = 0.0;
    sumSquaredDiff_z  = 0.0;
    sumSquaredDiff_vx = 0.0;
    sumSquaredDiff_vy = 0.0; 
    sumSquaredDiff_vz = 0.0;
    sumCov_x_vx       = 0.0;
    sumCov_y_vy       = 0.0;
    sumCov_z_vz       = 0.0;


    for (int k = 0; k < tracers_cluster; k++){
      //computing k-variable - mean value of the set k-variable belongs to
      difference_x  = xx[k][0] - x_mean;
      difference_y  = xx[k][1] - y_mean;
      difference_z  = xx[k][2] - z_mean;
      difference_vx = vv[k][0] - vx_mean;
      difference_vy = vv[k][1] - vy_mean;
      difference_vz = vv[k][2] - vz_mean;

      //cout << difference_x << " " << difference_y << " " << difference_z << " " << difference_vx << " " << difference_vy << " " << difference_vz  << endl; 

      //computing sum of squared differences of step above (for variance)
      sumSquaredDiff_x  += difference_x * difference_x;
      sumSquaredDiff_y  += difference_y * difference_y;
      sumSquaredDiff_z  += difference_z * difference_z;
      sumSquaredDiff_vx += difference_vx * difference_vx;
      sumSquaredDiff_vy += difference_vy * difference_vy;
      sumSquaredDiff_vz += difference_vz * difference_vz;

      //computing sum of squared differences (for covariance)
      sumCov_x_vx += difference_x * difference_vx;
      sumCov_y_vy += difference_y * difference_vy;
      sumCov_z_vz += difference_z * difference_vz;
    } 

    //computing variance of x,y,z,vx,vy,vz values
    variance_x  = sumSquaredDiff_x/tracers_cluster;
    variance_y  = sumSquaredDiff_y/tracers_cluster;
    variance_z  = sumSquaredDiff_z/tracers_cluster;
    variance_vx = sumSquaredDiff_vx/tracers_cluster;
    variance_vy = sumSquaredDiff_vy/tracers_cluster;
    variance_vz = sumSquaredDiff_vz/tracers_cluster;

    //computing covariance of (x,vx) (y,vy) (z,vz) values
    covar_x_vx = sumCov_x_vx/tracers_cluster;
    covar_y_vy = sumCov_y_vy/tracers_cluster;
    covar_z_vz = sumCov_z_vz/tracers_cluster;


    //computing emittance 
    emit_x = sqrt(variance_x * variance_vx - covar_x_vx * covar_x_vx);
    emit_y = sqrt(variance_y * variance_vy - covar_y_vy * covar_y_vy);
    emit_z = sqrt(variance_z * variance_vz - covar_z_vz * covar_z_vz);
    //cout << emit_x << " " << emit_y << " " << emit_z << endl;
    emit = pow((emit_x*emit_y*emit_z),1.0/3.0);


    
    //if (i%100 == 0){
      rootmeansq << t << "  " << x_mean <<"  " << y_mean << "  " << z_mean << "  " << vx_mean << "  " << vy_mean << "  " << vz_mean << "  "<< E_mean << "  "<< 2*(Ekin/tracers_cluster)/fabs(Epot/tracers_cluster)<<
      " " << rr2 << " " << a << " " << rmsr << " " << rmsvr << " " << emit << " " << r_mean << " " << vr_mean << endl;
    //}

    
  }

        // average lyapunov exponents over time
//for (int i = 0; i < neq; ++i) {
//    lambda[i] /= maxsteps;
//}
//cout << lambda[0] << " " << lambda[1] << " " << lambda[2] << endl;


  for (int k = 0; k < tracers_cluster; k++){
    //E0 = 0.0;
    rr2 = sqrt(xx[k][0]*xx[k][0] + xx[k][1]*xx[k][1] + xx[k][2]*xx[k][2] + A);
    //E0 += 0.5* (vv[k][0]*vv[k][0] + vv[k][1]*vv[k][1] + vv[k][2]*vv[k][2]) - (1.0/(rr2))*0.5;
    E0 = 0.5*(vv[k][0]*vv[k][0] + vv[k][1]*vv[k][1] + vv[k][2]*vv[k][2]) - 1.0/rr2;
    //energy1 << k << ' ' << E0 << ' ' << rr2 << endl;
  }


  //fraction.close();
  //energy0.close();
  //energy1.close();
  //energy2.close();
  //coordinates.close();
  rootmeansq.close();
  
  

  return 0;
}

// end main
///////////////////////////////////////////////////////////////////////////////



double dot_product(const double* v1, const double* v2, int size) {
    double result = 0.0;
    for (int i = 0; i < size; ++i) {
        result += v1[i] * v2[i];
    }
    return result;
}


// normalize vector
void normalize_vector(double* v, int size) {
    double magnitude = 0.0;

    // vector magnitude computation
    for (int i = 0; i < size; ++i) {
        magnitude += v[i] * v[i];
    }

    // no division by zero
    if (magnitude > 0.0) {
        magnitude = sqrt(magnitude);

        // normalize vector
        for (int i = 0; i < size; ++i) {
            v[i] /= magnitude;
        }
    }
}
