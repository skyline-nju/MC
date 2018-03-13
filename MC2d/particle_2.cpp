#include "particle_2.h"

HardSphere::HardSphere(double diameter) {
  sigma = diameter;
  sigma_square = sigma * sigma;
}

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

bool HardSphere::trivial_move(double delta_trans, Ran *myran,
  const PBC_2&pbc) {
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

void HardSphere::mc_move(double delta_trans, Ran *myran, const PBC_2 & pbc) {
  int valid_move = 0;
  while (valid_move < pos.size()) {
    if (trivial_move(delta_trans, myran, pbc)) {
      valid_move++;
    }
  }
}

HardJanus::HardJanus(double diameter, double lambda, double theta_m):
  HardSphere(diameter), theta_max(theta_m) {
  double sigma_ext = diameter * (1 + lambda);
  sigma_ext_square = sigma_ext * sigma_ext;
  cos_theta_max = std::cos(theta_max);
}

void HardJanus::create_rand(int nPar, Ran * myran, const PBC_2 & pbc) {
  create_rand_pos(nPar, myran, pbc);
  create_rand_ori(nPar, myran);
}

void HardJanus::create_rand_ori(int nPar, Ran *myran) {
  ori.reserve(nPar);
  for (int i = 0; i < nPar; i++) {
    double theta = myran->doub() * PI * 2;
    ori.emplace_back(std::cos(theta), std::sin(theta));
  }
}

bool HardJanus::trivial_move(double delta_trans, double delta_ori,
  double eps1, double eps2, Ran * myran, const PBC_2 & pbc) {
  int j = int(myran->doub() * pos.size());
  double u_cur = 0;
  for (int i = 0; i < pos.size(); i++) {
    if (i != j) {
      Vec_2f r_ij(pos[j] - pos[i]);
      pbc.nearest_dis(r_ij);
      cal_energy(u_cur, r_ij, ori[i], ori[j], eps1, eps2);
    }
  }

  Vec_2f new_pos;
  trans_move(pos[j], new_pos, delta_trans, myran, pbc);
  Vec_2f new_ori;
  ori_move(j, new_ori, delta_ori, myran);

  double u_next = 0;
  bool flag_overlap = false;
  for (int i = 0; i < pos.size(); i++) {
    if (i != j) {
      Vec_2f r_ij(new_pos - pos[i]);
      pbc.nearest_dis(r_ij);
      if (!cal_energy(u_next, r_ij, ori[i], new_ori, eps1, eps2)) {
        flag_overlap = true;
        break;
      }
    }
  }
  bool flag_accept = false;
  if (!flag_overlap) {
    double du = u_next - u_cur;
    if (du <= 0 || myran->doub() < std::exp(-du)) {
      flag_accept = true;
      pos[j] = new_pos;
      ori[j] = new_ori;
    }
  }
  return flag_accept;
}

bool HardJanus::cal_energy(double &u, const Vec_2f& r_ij,
const Vec_2f& ori_i, const Vec_2f& ori_j, double uAA, double uBB) const {
  double r_ij_square = r_ij.square();
  if (r_ij_square < sigma_square) {
    return false;
  } else if (r_ij_square < sigma_ext_square) {
    double r_cos_theta_max = cos_theta_max * std::sqrt(r_ij_square);
    double ori_i_dot_r_ij = ori_i.dot(r_ij);
    if (ori_i_dot_r_ij > r_cos_theta_max) {
      if (- (ori_j.dot(r_ij)) > r_cos_theta_max) {
        u += uAA;
      }
    } else if (-ori_i_dot_r_ij > r_cos_theta_max) {
      if (ori_j.dot(r_ij) > r_cos_theta_max) {
        u += uBB;
      }
    }
  }
  return true;
}

bool HardJanus::cal_energy(double &u, const Vec_2f& r_ij,
  const Vec_2f& ori_i, const Vec_2f& ori_j, double uAA) const {
  double r_ij_square = r_ij.square();
  if (r_ij_square < sigma_square) {
    return false;
  } else if (r_ij_square < sigma_ext_square) {
    double r_cos_theta_max = cos_theta_max * std::sqrt(r_ij_square);
    if (ori_i.dot(r_ij) > r_cos_theta_max && -(ori_j.dot(r_ij)) > r_cos_theta_max) {
      u += uAA;
    }
  }
  return true;
}

HardABA::HardABA(double diameter, double lambda, double theta_m):
  HardJanus(diameter, lambda, theta_m) {
  
}

bool HardABA::cal_energy(double & u, const Vec_2f & r_ij,
  const Vec_2f & ori_i, const Vec_2f & ori_j, double eps, bool flagAA) const {
  double r_ij_square = r_ij.square();
  if (r_ij_square < sigma_square) {
    return false;
  } else if (r_ij_square < sigma_ext_square) {
    double r_cos_theta_max = cos_theta_max * std::sqrt(r_ij_square);
    double prod = ori_i.dot(r_ij);
    if (flagAA) {
      if (prod > r_cos_theta_max || prod < -r_cos_theta_max) {
        prod = ori_j.dot(r_ij);
        if (prod > r_cos_theta_max || prod < -r_cos_theta_max) {
          u += eps;
        }
      }
    } else {
      if (prod < r_cos_theta_max && prod > -r_cos_theta_max) {
        prod = ori_j.dot(r_ij);
        if (prod < r_cos_theta_max && prod > -r_cos_theta_max) {
          u += eps;
        }
      }
    }
  }
  return true;
}

bool HardABA::cal_energy(double & u, const Vec_2f & r_ij, const Vec_2f & ori_i, const Vec_2f & ori_j,
  double uAA, double uBB, double uAB) const {
  double r_ij_square = r_ij.square();
  if (r_ij_square < sigma_square) {
    return false;
  } else if (r_ij_square < sigma_ext_square) {
    double r_cos_theta_max = cos_theta_max * std::sqrt(r_ij_square);
    double prod = ori_i.dot(r_ij);
    if (prod > r_cos_theta_max || prod < -r_cos_theta_max) {
      prod = ori_j.dot(r_ij);
      if (prod > r_cos_theta_max || prod < -r_cos_theta_max) {
        u += uAA;
      } else {
        u += uAB;
      }
    } else if (prod < r_cos_theta_max && prod > -r_cos_theta_max) {
      prod = ori_j.dot(r_ij);
      if (prod < r_cos_theta_max && prod > -r_cos_theta_max) {
        u += uBB;
      } else {
        u += uAB;
      }
    } else {
      u += uAB;
    }
  }
  return true;
}

HardMix::HardMix(double sigma_L, double sigma_S,
  double theta_m_L, double theta_m_S, double lambda) {
  large_sphere = new HardJanus(sigma_L, lambda, theta_m_L);
  small_sphere = new HardABA(sigma_S, lambda, theta_m_S);
  double dis = sigma_L + sigma_S;
  double dis_ext = dis * lambda;
  inter_dis_square = dis * dis;
  inter_force_range_square = dis_ext * dis_ext;
}

HardMix::~HardMix() {
  delete large_sphere;
  delete small_sphere;
}

void HardMix::create_rand(double nL, double nS, Ran * myran, const PBC_2 & pbc) {
  HardSphere tmp(large_sphere->get_sigma());
  ntot = nL + nS;
  n_large = nL;
  n_small = nS;
  tmp.create_rand(ntot, myran, pbc);
  for (int i = 0; i < ntot; i++) {
    if (i < nL) {
      large_sphere->pos[i] = tmp.pos[i];
    } else {
      small_sphere->pos[i - nL] = tmp.pos[i];
    }
  }
  large_sphere->create_rand_ori(nL, myran);
  small_sphere->create_rand_ori(nS, myran);
}

bool HardMix::cal_mix_energy(double & u, const Vec_2f & r_ij, const Vec_2f & ori_i,
  const Vec_2f & ori_j, double u_SA_LA) {
  double r_ij_square = r_ij.square();
  if (r_ij_square < inter_dis_square) {
    return false;
  } else if (r_ij_square < inter_force_range_square) {
    double r = std::sqrt(r_ij_square);
    if (ori_i.dot(r_ij) > small_sphere->get_theta_m() * r &&
        ori_j.dot(r_ij) < -large_sphere->get_theta_m() * r) {
      u += u_SA_LA;
    }
  }
  return true;
}
