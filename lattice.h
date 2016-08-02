#ifndef INCLUDED_LATTICE
#define INCLUDED_LATTICE

class lattice
{
  private:

    int *x;
    int N;
    const int cutoff;

    lattice& operator=(const lattice&);
    lattice(const lattice&);

  public:
    explicit lattice(int);
    ~lattice();

    int* get_handle();

    void check_config(void); 
    void print_config(void); 
    void init_config(int);
    void write_config(char*);
    void read_config(char*);

    double potential_energy(double);
    double kinetic_energy(double);

};

#endif
