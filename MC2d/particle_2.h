/**
 * @file particle_2.h
 * @author Yu Duan (duanyu100@yeah.net)
 * @brief 
 * @version 0.1
 * @date 2018-12-12
 * 
 * Class for 2d particle. Only the HardSphere class is implemented.
 * 
 * @copyright Copyright (c) 2018
 * 
 */
#ifndef PARTICLE_2H
#define PARTICLE_2H
#include <vector>
#include <cmath>
#include <iostream>
#include "vect.h"
#include "boundary.h"
#include "rand.h"

const double PI = 3.14159265358979;

typedef Vec_2<double> Vec_2f;

/**
 * @brief class for hard sphere in 2D. The two particles cannot be overlaped.
 * The minimum distance bwteen them is the diameter.
 * 
 */
class HardSphere {
public:
  HardSphere(double diameter);
  bool check_overlap(int i, int j, const PBC_2& pbc) const;
  bool check_overlap(int i, const Vec_2f &a, const PBC_2 &pbc) const;
  void create_rand_pos(int nPar, Ran *myran, const PBC_2& pbc);
  void create_rand(int nPar, Ran *myran, const PBC_2& pbc);
  void trans_move(const Vec_2f &old_pos, Vec_2f &new_pos,
                  double delta_trans, Ran* myran, const PBC_2& pbc) const;
  void trans_move(int i, Vec_2f &new_pos, double delta_trans,
                  Ran* myran, const PBC_2 &pbc) const;
  bool trivial_move(double delta_trans, Ran *myran, const PBC_2& pbc);
  void mc_move(double delta_trans, Ran *myran, const PBC_2& pbc);
  int get_n() const { return pos.size();}
  double get_sigma() const { return sigma; }

  std::vector<Vec_2f> pos;  // array of particle position
protected:
  double sigma; // diameter of one particle
  double sigma_square;
};

inline bool HardSphere::check_overlap(int i, int j, const PBC_2& pbc) const {
  return pbc.nearest_dis_square(pos[i], pos[j]) <= sigma_square;
}

inline bool HardSphere::check_overlap(int i, const Vec_2f &a, const PBC_2 &pbc) const {
  return pbc.nearest_dis_square(pos[i], a) <= sigma_square;
}

/**
 * @brief A attempled translational move.
 * 
 * @param old_pos Old position before move
 * @param new_pos New position after move
 * @param delta_trans Magnitude of translational move
 * @param myran Random number generator
 * @param pbc Class for perodic boundary condition
 */
inline void HardSphere::trans_move(const Vec_2f& old_pos, Vec_2f &new_pos,
                                   double delta_trans, Ran* myran,
                                   const PBC_2& pbc) const {
  new_pos.x = old_pos.x + (myran->doub() - 0.5) * delta_trans;
  new_pos.y = old_pos.y + (myran->doub() - 0.5) * delta_trans;
  pbc.wrap(new_pos);
}

/**
 * @brief A attempled translational move.
 * 
 * @param i The index of particle to move
 * @param new_pos New position after move
 * @param delta_trans Magnitude of translational move
 * @param myran Random number generator
 * @param pbc Class for perodic boundary condition
 */
inline void HardSphere::trans_move(int i, Vec_2f &new_pos, double delta_trans,
                                   Ran *myran, const PBC_2& pbc) const {
  new_pos.x = pos[i].x + (myran->doub() - 0.5) * delta_trans;
  new_pos.y = pos[i].y + (myran->doub() - 0.5) * delta_trans;
  pbc.wrap(new_pos);
}

// calculate particle number from the packing fraction phi
template <typename T>
int cal_particle_number_2(double phi, T Lx, T Ly, T sigma) {
  double phi_max = PI / (2 * std::sqrt(3.0));
  if (phi > phi_max) {
    std::cout << "Input packing fraction phi = " << phi
      << " is larger than phi_max = " << phi_max << std::endl;
    exit(1);
  } else {
    double a = sigma * 0.5;
    return int(std::round(phi * Lx * Ly / (PI * a * a)));
  }
}
#endif
