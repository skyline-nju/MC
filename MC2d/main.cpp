#include <iostream>
#include "ioData.h"

using namespace std;

int main(int argc, char argv[]) {
  double L = 20;
  double phi = 0.4;
  double sigma = 1;
  int N = cal_particle_number_2(phi, L, L, sigma);
  PBC_2 pbc(L, L);
  Ran myran(1);
  int mc_steps = 20000;
  double delta_trans = 0.1;
  double delta_ori = 0.1;

  //HardSphere hs(sigma);
  //hs.create_rand(N, &myran, pbc);
  //set_snap(phi, L, L);
  //for (int i = 0; i < mc_steps; i++) {
  //  hs.mc_move(delta_trans, &myran, pbc);
  //  if (i % 100 == 0)
  //    output_snap(hs, i);
  //}

  //double theta_max = PI / 2;
  //double eps1 = -5;
  //double eps2 = -5;
  //HardJanus js(sigma, 0.1, theta_max);
  //js.create_rand(N, &myran, pbc);
  //set_Janus_snap(phi, L, L, theta_max, eps1, eps2);

  //auto func_enegy = [&js, eps1, eps2](double &u, Vec_2f &r_ij, Vec_2f &ori_i, Vec_2f &ori_j) {
  //  return js.cal_energy(u, r_ij, ori_i, ori_j, eps1, eps2);
  //};
  //for (int i = 0; i < mc_steps; i++) {
  //  int trivial_move_count = js.mc_move(
  //    delta_trans, delta_ori, &myran, pbc, func_enegy);
  //  if (i % 10 == 0) {
  //    output_snap(js, i, 2);
  //    std::cout << "i = " << i << "\tacc rate = "
  //      << double(js.get_n()) / trivial_move_count << std::endl;
  //  }
  //}

  double theta_max = 60.0 / 180 * PI;
  double uAA = -4;
  double uBB = 2;
  double uAB = 1;
  HardABA aba(sigma, 0.05, theta_max);
  aba.create_rand(N, &myran, pbc);
  set_ABA_snap(phi, L, L, theta_max, uAA, uBB);

  auto func_enegy = [&aba, uAA, uBB, uAB](double &u, Vec_2f &r_ij, Vec_2f &ori_i, Vec_2f &ori_j) {
    return aba.cal_energy(u, r_ij, ori_i, ori_j, uAA, uBB, uAB);
  };

  for (int i = 0; i < mc_steps; i++) {
    int trivial_move_count = aba.mc_move(
      delta_trans, delta_ori, &myran, pbc, func_enegy);
    if (i % 20 == 0) {
      output_snap(aba, i);
      std::cout << "i = " << i << "\tacc rate = "
        << double(aba.get_n()) / trivial_move_count << std::endl;
    }
  }

}