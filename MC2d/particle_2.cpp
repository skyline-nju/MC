#include "particle_2.h"

HardSphere::HardSphere(double diameter) {
  sigma = diameter;
  sigma_square = sigma * sigma;
}

/**
 * @brief Create particles with random position
 * 
 * @param nPar Number of particles to create
 * @param myran Random number generator
 * @param pbc Class for Periodic boundary condition
 */
void HardSphere::create_rand_pos(int nPar, Ran *myran, const PBC_2& pbc) {
  pos.reserve(nPar);
  while (pos.size() < nPar) {
    Vec_2f xy(myran->doub() * pbc.Lx(), myran->doub() * pbc.Ly());
    bool is_overlapped = false;
    for (int i = 0; i < pos.size(); i++) {
      if (check_overlap(i, xy, pbc)) {
        is_overlapped = true;
        break;
      }
    }
    if (!is_overlapped) {
      pos.push_back(xy);
    }
  }
}

void HardSphere::create_rand(int nPar, Ran * myran, const PBC_2 & pbc) {
  create_rand_pos(nPar, myran, pbc);
}

/**
 * @brief Attempt to move a randomly selected particle
 * 
 * @param delta_trans Magnitude of translational move
 * @param myran Random number generator
 * @param pbc Periodic boundary condition
 * @return true The trivial move is successful
 * @return false The trivial move is unsuccessful
 */
bool HardSphere::trivial_move(double delta_trans, Ran *myran, const PBC_2&pbc) {
  int k = int(myran->doub() * pos.size());
  Vec_2f new_pos;
  trans_move(pos[k], new_pos, delta_trans, myran, pbc);
  bool is_overlapped = false;
  for (int i = 0; i < pos.size(); i++) {
    if (i != k && check_overlap(i, new_pos, pbc)) {
      is_overlapped = true;
      break;
    }
  }
  if (is_overlapped) {
  } else {
    pos[k] = new_pos;
  }
  return !is_overlapped;
}

/**
 * @brief Perform one MC move
 * 
 * @param delta_trans Magnitude of translational move
 * @param myran Random number generator
 * @param pbc Periodic boundary condition
 */
void HardSphere::mc_move(double delta_trans, Ran *myran, const PBC_2 & pbc) {
  int valid_move = 0;
  while (valid_move < pos.size()) {
    if (trivial_move(delta_trans, myran, pbc)) {
      valid_move++;
    }
  }
}
