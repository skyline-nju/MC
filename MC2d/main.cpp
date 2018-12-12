#include <iostream>
#include <ctime>

#include "output.h"

using namespace std;

int main(int argc, char argv[]) {
  double L = 20;
  double phi = 0.3;
  double sigma = 1;
  int N = cal_particle_number_2(phi, L, L, sigma);
  PBC_2 pbc(L, L);
  Ran myran(1);
  int mc_steps = 20000;
  double delta_trans = 0.1;

  const auto t_beg = chrono::system_clock::now();

  cout << "total particles: " << N << endl;
  cout << "domain length: " << L << endl;
  cout << "packing fraction: " << phi << endl;
  cout << "total time steps: " << mc_steps << endl;
  cout << "start simulation" << endl;
  cout << "--------" << endl;

  HardSphere hs(sigma);
  hs.create_rand(N, &myran, pbc);
  ini_snap_output(phi, L, L);
  for (int i = 0; i < mc_steps; i++) {
    hs.mc_move(delta_trans, &myran, pbc);
    if (i % 100 == 0) {
      if (i % 5000 == 0)
        cout << "t = " << i << endl;
      output_snap(hs, i);
    }
  }
  const auto t_end = chrono::system_clock::now();
  chrono::duration<double> elapsed_sec = t_end - t_beg;
  cout << "--------" << endl;
  cout << "finish simulation, elapsed time: " << elapsed_sec.count() << "s" << endl;
}