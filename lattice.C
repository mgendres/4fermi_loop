#include <iostream>
#include <math.h>
using namespace std;
#include "lattice.h"
#include "params.h"
#include "utils.h"


lattice::lattice(int N_) : cutoff(100) {

  N = N_;
  x = new int [ (N+2)*T ];

}

lattice::~lattice() {
  delete [] x;
}

void lattice::check_config() {
  cout << "Checking configs..." << endl;
  for (int t=0; t<T; ++t) {
    if ( x[to_idx(0,t)] != -L-1 ) {
      cout << "Error found: (t, n) = " << t << " " << 0 << " " << endl;
    };
    if ( x[to_idx(N+1,t)] != L+1 ) {
      cout << "Error found: (t, n) = " << t << " " << N+1 << " " << endl;
    };
    for (int n=0; n<=N; ++n) {
      if (x[to_idx(n,t)] >= x[to_idx(n+1,t)]) {
        cout << "Error found: (t, n) = " << t << " " << n << " " << endl;
      }
    }
  }
}


void lattice::init_config(int dx) {
  int idx;
  for (int t=0; t<T; ++t) {
    idx = to_idx(0,t);
    x[idx] = -L-1;
    idx = to_idx(N+1,t);
    x[idx] = L+1;
    for (int n=1; n<=N; ++n) {
      idx = to_idx(n,t);
      x[idx] = -N/2 + n - 1;
      x[idx] *= dx+1;
    }
  }
}

void lattice::print_config() {
  int idx;
  for (int t=0; t<T; ++t) {
    for (int n=0; n<=N+1; ++n) {
      idx = to_idx(n,t);
      cout << x[idx] << " ";
    }
    cout << endl;
  }
  cout << endl;
}

double lattice::kinetic_energy(double m) {

  int dx;
  int tm;
  int x_n_t;
  int x_n_tm;
  int x_nm_t;
  int x_nm_tm;

  double f0 = (1.0 + m - m*sqrt( (2.0 + m ) / m )) / (2.0 + m);
  double f1 = f0 + 1.0/(2.0+m);

  double kin = 0.0; 
  if (2*L+1 < cutoff) {
    for (int t=0; t<T; ++t ) {
      tm = (t-1+T)%T;
      for (int n=1; n<=N+1; ++n ) {
        x_n_t   = x[to_idx(n, t  )];
        x_n_tm  = x[to_idx(n, tm )];
        x_nm_t  = x[to_idx(n-1, t )];
        x_nm_tm = x[to_idx(n-1, tm )];
        dx  = min( x_n_t, x_n_tm);
        dx -= max( x_nm_t, x_nm_tm );
        kin += d_log_z_d_log_m_inverse( dx-1, m );
        kin += abs( x_n_t - x_n_tm);
      }
      kin -= d_log_z_d_log_m_inverse( 2*L+1, m );
      //kin -= f0 + f1*(2*L+1);
    }
  } else {
    for (int t=0; t<T; ++t ) {
      tm = (t-1+T)%T;
      for (int n=1; n<=N+1; ++n ) {
        x_n_t   = x[to_idx(n, t  )];
        x_n_tm  = x[to_idx(n, tm )];
        x_nm_t  = x[to_idx(n-1, t )];
        x_nm_tm = x[to_idx(n-1, tm )];
        dx  = min( x_n_t, x_n_tm);
        dx -= max( x_nm_t, x_nm_tm );
        if (dx-1 < cutoff) {
          kin += d_log_z_d_log_m_inverse( dx-1, m ) - (f0 + f1*(dx-1) );
        }
        kin += (1-f1)*abs( x_n_t - x_n_tm);
      }
      kin -= N/(2.0+m);
    }
  }
  kin /= T;

  return -kin;
}

double lattice::potential_energy(double k) {
  int idx;
  double pot = 0.0; 
  for (int t=0; t<T; ++t ) {
    for (int n=1; n<=N; ++n ) {
      idx = to_idx(n,t);
      pot += v( k, x[idx]  );
    }
  }
  return pot/T;
}

int* lattice::get_handle() {
  return x;
}
