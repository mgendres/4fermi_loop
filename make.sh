#!/bin/bash
#g++ -g -O main.C ranlux-3.3/ranlxd.c ranlux-3.3/ranlxs.c 
g++ -O3 mt19937.C utils.C lattice.C  main.C 
