#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <ctime>

// Global parameters
#include "params.h"
#include "utils.h"
#include "mt19937.h"
#include "lattice.h"

// Z_n(m) can be computed to machine precision only up to n = 584.
// Therefore, 2*L+1 < 584 unless an approximation to Z_n(m) is used
// "cutoff" is the value of n at which Z_n(m) is approximated 
const int cutoff = 100;

double mean_occupancy(int *occupancy) {
  int ans = 0;
  for (int i=0; i< (2*L+1)*T; ++i) {
    ans += occupancy[i]/4;
  }
  return (double) ans/ (double) T;
}

void print_occupancy(int *occupancy) {
  for (int t=0; t<T; ++t) {
    for (int x=-L; x <= L; ++x) {
      cout << occupancy[ to_bond_idx(x,t) ]/4 << " ";
    }
    cout << endl;
  }
  cout << endl;
}

int main(void)
{

  cout << setprecision(15);
  print_params(cout);

  // Lattice confis 
  lattice *lat[4];
  for (int s=0; s<4; ++s) {
    lat[s] = new lattice( N[s] );
   }

  // Initialize random number generator
  mt19937 rng(0);
  if (rseed < 0) {
    rng.init_genrand( time(0) );
  } 
  if (rseed > 0) {
    rng.init_genrand(rseed);
  }
  if (rseed==0) {
    rng.read_state("rand.bin");
  }// else {
    // Warm up the generator
    for (int i=0; i<20; i++) { rng.genrand_real2(); }
//  }


  if (!restart) {
    for (int s=0; s<4; ++s) { lat[s]->init_config(0); }
  } else {
    // Import configuration
    ifstream conf_file ("conf.bin", ios::in| ios::binary);
    for (int s=0; s<4; ++s) {
      conf_file.read( (char*) lat[s]->get_handle() ,  (N[s] + 2)*T * sizeof(int)  );
    }
    conf_file.close();
  }

//  for (int s=0; s<4; ++s) { print_config(x, s); }

  // Check configurations
  for (int s=0; s<4; ++s) { lat[s]->check_config(); }

  // Tabulate the z values
  double r[4];
  double r0[4];
  double **z_tab = new double* [4];
  for (int s=0; s<4; ++s) {
    r[s] =  0.5* ( 1.0 + 1.0/m[s] + sqrt(1.0 + 2.0/m[s]) );
    r0[s] =  pow(2.0*m[s],2);
    z_tab[s] = new double[2*L+1];
    for (int dx=0; dx< 2*L+1; ++dx) {
      z_tab[s][dx] = z(dx, m[s]);
//      cout << z_tab[s][dx] << endl;
    }
  }

  // Tablate potential
  double **v_tab = new double* [4];
  for (int s=0; s<4; ++s) {
    v_tab[s] = new double [2*L+1];
    for (int x=-L; x<=L; ++x) {
      v_tab[s][L+x] = exp( -v(k[s], x) );
    }
  }
//  for (int x=-L; x<=L; ++x) {
//    cout << x << " " << v_tab[2][L+x] << endl;
//  }
  

  // Occupancy table 
  int *x;
  int *occupancy= new int [(2*L+1)*T];
  for (int i=0; i< (2*L+1)*T; ++i) { occupancy[i] = 0; }
  for (int s=0; s<4; ++s) {
    x = lat[s]->get_handle();
    for (int n=1; n<=N[s]; ++n) {
      for (int t=0; t<T; ++t) {
        int idx = to_idx(n,t);
        occupancy[ to_bond_idx( x[idx], t) ] += 1;
      }
    }
  }

//for (int s=0; s<4; ++s) {
//  cout << mean_pot( x, s) << "   ";
//}
//cout << endl;

//  for (int s=0; s<4; ++s) {
//    cout << mean_kin( x, s) << "   ";
//  }
//  cout << endl;

  int x_trial;
  int crossQ;
  double w;
  int idx_n_t;
  int idx_n_tp;
  int idx_n_tm;
  int idx_np_t;
  int idx_np_tp;
  int idx_np_tm;
  int idx_nm_t;
  int idx_nm_tp;
  int idx_nm_tm;
  int x_;
  int x_n_tp;
  int x_n_tm;
  int x_np_t;
  int x_nm_t;
  int x_np_tp;
  int x_np_tm;
  int x_nm_tp;
  int x_nm_tm;
  int tp, tm;
  int y1, y2, y3;
  double mean_occupancy_block = 0.0;
  double mean_pot_block[4] = { 0.0, 0.0, 0.0, 0.0 };
  double mean_kin_block[4] = { 0.0, 0.0, 0.0, 0.0 };


  time_t start;
  time_t end;

  start = time(0);

  ofstream out_file;
  if (restart==0) {
     out_file.open( "out" );
     print_params(out_file);
   } else {
     out_file.open( "out", ios::app );
  } 
  out_file << setprecision(15);

  for (long int i=1; i<=n_conf; ++i) {

    if (0) {
      for (int s=0; s<4; ++s) { lat[s]->print_config(); }
    }


    // Compute observable on every meas_freq config
   if (i%meas_freq==0) {
      mean_occupancy_block += mean_occupancy(occupancy);
      for (int s=0; s<4; ++s) {
        mean_pot_block[s] += lat[s]->potential_energy( k[s] );
        mean_kin_block[s] += lat[s]->kinetic_energy( m[s] );
      }
    }

    if (i%(meas_freq*block_freq)==0) {
      // output some results
      out_file << mean_occupancy_block / block_freq << " ";
      mean_occupancy_block = 0.0;
      for (int s=0; s<4; ++s) {
        out_file << mean_pot_block[s] / block_freq << " ";
        mean_pot_block[s] = 0.0;
      }
      for (int s=0; s<4; ++s) {
        out_file << mean_kin_block[s] / block_freq << " ";
        mean_kin_block[s] = 0.0;
      }
      out_file << endl;
    }

    for (int s=0; s<4; ++s) {
      x = lat[s]->get_handle();
      for (int t=0; t<T; ++t) {
        tp = (t+1+T)%T;
        tm = (t-1+T)%T;


        idx_n_t   = to_idx(0,t);
        idx_n_tp  = to_idx(0,tp);
        idx_n_tm  = to_idx(0,tm);
        idx_np_t  = idx_n_t  + T;
        idx_np_tp = idx_n_tp + T;
        idx_np_tm = idx_n_tm + T;
        idx_nm_t  = idx_n_t  - T;
        idx_nm_tp = idx_n_tp - T;
        idx_nm_tm = idx_n_tm - T;

        for (int n=1; n<=N[s]; ++n) {

          idx_n_t   += T;
          idx_n_tp  += T;
          idx_n_tm  += T;
          idx_np_t  += T;
          idx_np_tp += T;
          idx_np_tm += T;
          idx_nm_t  += T;
          idx_nm_tp += T;
          idx_nm_tm += T;

          x_      = x[idx_n_t];
          x_n_tp  = x[idx_n_tp];
          x_n_tm  = x[idx_n_tm];
          x_np_t  = x[idx_np_t];
          x_nm_t  = x[idx_nm_t];
          x_np_tp = x[idx_np_tp];
          x_np_tm = x[idx_np_tm];
          x_nm_tp = x[idx_nm_tp];
          x_nm_tm = x[idx_nm_tm];

          if ( rng.genrand_real2()  > 0.5 ) {
            x_trial = x_ + 1;
            crossQ = (x_trial != x_np_t  ) &&
                     (x_trial != x_np_tp ) &&
                     (x_trial != x_np_tm );
          } else {
            x_trial = x_ - 1;
            crossQ = (x_trial != x_nm_t  ) &&
                     (x_trial != x_nm_tp ) &&
                     (x_trial != x_nm_tm );
          }

          if ( crossQ ) {

            w  = v_tab[s][ L + x_trial];
            w /= v_tab[s][ L + x_ ]; 

            y1 = min( x_n_tm , x_trial );
            y2 = min( x_n_tm , x_ );
            if (y1!=y2){
              y3 = max( x_nm_tm, x_nm_t ) + 1;
              if (y1-y3 < cutoff) {
                w *= z_tab[s][ y1 - y3 ]; 
                w /= z_tab[s][ y2 - y3 ];
              } else {
                if ( y1 > y2) { w *= r[s]; } else { w /= r[s]; } }
            }

            y1 = max( x_n_tm , x_trial );
            y2 = max( x_n_tm , x_ );
            if (y1!=y2){
              y3 = min( x_np_tm, x_np_t ) - 1;
              if (y3-y1 < cutoff) {
                w *= z_tab[s][ y3 - y1 ];
                w /= z_tab[s][ y3 - y2 ];
              } else {
                if ( y1 < y2) { w *= r[s]; } else { w /= r[s]; }
              }
            }

            y1 = min( x_trial, x_n_tp);
            y2 = min( x_, x_n_tp);
            if (y1!=y2){
              y3 = max( x_nm_t , x_nm_tp ) + 1;
              if (y1-y3 < cutoff) {
                w *= z_tab[s][ y1 - y3 ];
                w /= z_tab[s][ y2 - y3 ];
              } else {
                if ( y1 > y2) { w *= r[s]; } else { w /= r[s]; }
              }
            }

            y1 = max( x_trial   , x_n_tp );
            y2 = max( x_, x_n_tp );
            if (y1!=y2){
              y3 = min( x_np_t, x_np_tp ) - 1;
              if (y3-y1 < cutoff) {
                w *= z_tab[s][ y3 - y1 ];
                w /= z_tab[s][ y3 - y2 ];
              } else {
                if ( y1 < y2) { w *= r[s]; } else { w /= r[s]; }
              }
            }

            y1 =  abs( x_trial   - x_n_tm );
            y1 += abs( x_trial   - x_n_tp );
            y2 =  abs( x_ - x_n_tm );
            y2 += abs( x_ - x_n_tp );
            if (y1 != y2) {
              if ( y1 < y2 ) { w *= r0[s]; } else { w /= r0[s]; }
            }

            y1 =  ( occupancy[ to_bond_idx( x_trial ,t) ] + 1 ) / 4;
            y1 += ( occupancy[ to_bond_idx( x_, t) ] - 1 ) / 4;
            y2 =   occupancy[ to_bond_idx( x_trial ,t) ] / 4;
            y2 +=   occupancy[ to_bond_idx( x_, t) ]  / 4;
            if (y1 != y2) {
              if (y1 < y2) { w /= 1.0 + g; } else { w *= 1.0 + g;  };
            }

            if ( w > rng.genrand_real2() ) {
              occupancy[ to_bond_idx(x_     ,t) ] -= 1;
              occupancy[ to_bond_idx(x_trial,t) ] += 1;
              x[idx_n_t] = x_trial;
            }

          }


        }
      }
    }
  }

  out_file.close();

  {
    // Store configs in a file
    ofstream conf_file ("conf.bin", ios::out | ios::binary);
    for (int s=0; s<4; ++s) {
      conf_file.write( (char*) lat[s]->get_handle(),  (N[s] + 2)*T * sizeof(int) );
    }
    conf_file.close();
    // Store random state in a file
    rng.write_state("rand.bin");
  }

//  for (int s=0; s<4; ++s) { print_config(x, s); }
//  print_occupancy(occupancy);


//  // Compute z ratios... check numerical precision for large n and all m
//  for (int dx=0; dx<1000; ++dx) {
//    double m=1.6;
//    cout <<  dx << " ";
//    cout <<  z(dx+1, m) / z(dx, m) << " ";
//    cout <<  z(dx-1, m) / z(dx, m) << " ";
//    cout << endl;
//  }

//  // Compute z ratios... check numerical precision for large n and all m
//  for (int dx=0; dx<201; ++dx) {
//    cout <<  dx << " ";
//    cout <<  d_log_z_d_log_m_inverse(dx, m[1]) << " ";
//    cout << endl;
//  }


  for (int s=0; s<4; ++s) {
    delete [] z_tab[s];
    delete [] v_tab[s];
    delete lat[s];
  }
  delete [] z_tab;
  delete [] v_tab;
  delete [] occupancy;


  end = time(0);
  cout << "Elapsed time: " << end - start << " seconds" << endl;
  cout << "Configs/hour: " << n_conf*3600.0 / (double) (end - start) << endl;


  return(0);

}

