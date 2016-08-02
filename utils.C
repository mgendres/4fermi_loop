#include <math.h>
#include <stdlib.h>
using namespace std;

#include "params.h"
#include "utils.h"

double trap_size(double m, double k) { return pow( m * k , -0.25 ); }
double trap_freq(double m, double k) { return pow( k/m , 0.5); }

double v(double k, int x ) { return k*x*x/2.0; }

int min(int a, int b) {
  if (a<=b) { return a; } else { return b; }
}

int max(int a, int b) {
  if (a>=b) { return a; } else { return b; }
}



void print_params(std::ostream& out) {

  out << endl;
  out << "Lattice shape: (2L+1) x T with L = " << L << " and T = " << T << endl;
  out << "Number of species: 4" << endl;
  out << "Number of fermions: ";
  out << "( " <<  N[0];
  out << ", " <<  N[1];
  out << ", " <<  N[2];
  out << ", " <<  N[3];
  out << " )" <<  endl;

  out << "Coupling g: ";
  out << g << endl;

  out << "Masses: ";
  out << "( " << m[0];
  out << ", " << m[1];
  out << ", " << m[2];
  out << ", " << m[3];
  out << " )" << endl;

  out << "Spring constants: ";
  out << "( " << k[0];
  out << ", " << k[1];
  out << ", " << k[2];
  out << ", " << k[3];
  out << " )" << endl;

  out << "Trap size:";
  out << "( " << trap_size(m[0],k[0]);
  out << ", " << trap_size(m[1],k[1]);
  out << ", " << trap_size(m[2],k[2]);
  out << ", " << trap_size(m[3],k[3]);
  out << " )" << endl;

  out << "Trap frequency:";
  out << "( " << trap_freq(m[0],k[0]);
  out << ", " << trap_freq(m[1],k[1]);
  out << ", " << trap_freq(m[2],k[2]);
  out << ", " << trap_freq(m[3],k[3]);
  out << " )" << endl;

  out << "Random generator (ranlxd) seed: ";
  out << rseed << endl;

  out << "Random generator (ranlxd) luxlev: ";
  out << luxlev << endl;

  out << "Measurement (blocking) frequency: ";
  out << meas_freq << " (";
  out << block_freq << ")" << endl;
  out << endl;
  out << endl;

}


int to_idx(int n, int t) {
  return t+n*T;
}

int to_bond_idx(int x, int t) {
  return t+(L+x)*T;
}

double z(int dx, double m) {
  double ans = 1.0;
  ans /= pow(2.0, dx+1);
  ans /= sqrt(1.0 + 2.0/m);
  ans *= pow(1.0 + 1.0/m + sqrt(1.0 + 2.0/m), dx + 1) \
       - pow(1.0 + 1.0/m - sqrt(1.0 + 2.0/m), dx + 1);
  return ans;
}

double d_log_z_d_log_m_inverse(int dx, double m) {
  double ans = 0.0;
  ans +=  pow(1.0 + 1.0/m + sqrt(1.0 + 2.0/m), dx) * ( 1.0 + 1.0/sqrt(1+2.0/m) );
  ans -=  pow(1.0 + 1.0/m - sqrt(1.0 + 2.0/m), dx) * ( 1.0 - 1.0/sqrt(1+2.0/m) );
  ans *= dx+1;
  ans /= ( pow(1.0 + 1.0/m + sqrt(1.0 + 2.0/m), dx + 1) - pow(1.0 + 1.0/m - sqrt(1.0 + 2.0/m), dx + 1) );
  ans -= 1.0 / ( 1 + 2.0/m );
  ans /= m;
  return ans;
}


