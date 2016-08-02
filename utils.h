#include <iostream>

#ifndef INCLUDED_UTILS
#define INCLUDED_UTILS

double trap_size(double, double);
double trap_freq(double, double);

double v(double, int);

int min(int, int);
int max(int, int);

void print_params(std::ostream& );

int to_idx(int, int);
int to_bond_idx(int, int);

double z(int, double);
double d_log_z_d_log_m_inverse(int, double);

#endif
