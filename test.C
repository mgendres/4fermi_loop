#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <ctime>



#include "params.h"
#include "lattice.h"


int main() {

  lattice *x[4];
  x[0] = new lattice(4);
  x[1] = new lattice(4);
  x[2] = new lattice(4);
  x[3] = new lattice(4);

  cout << x[0]->operator[](3) << endl;
  cout << *x[0][3] << endl;
//  x[3] = 3;
//  cout << x[3] << endl;

  return 0;
}
