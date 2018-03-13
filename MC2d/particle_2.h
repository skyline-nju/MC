#ifndef PARTICLE_2H
#define PARTICLE_2H
#include <vector>
#include <cmath>
#include "vect.h"
#include "boundary.h"
#include "rand.h"
#include "comn.h"

typedef Vec_2<double> Vec_2f;

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
  void ori_move(const Vec_2f &old_ori, Vec_2f &new_ori,
    double delta_ori, Ran* myran) const;
  bool trivial_move(double delta_trans, Ran *myran, const PBC_2& pbc);
  void mc_move(double delta_trans, Ran *myran, const PBC_2& pbc);
  int get_n() const { return pos.size();}
  double get_sigma() const { return sigma; }

  std::vector<Vec_2f> pos;
protected:
  double sigma;
  double sigma_square;
};

inline bool HardSphere::check_overlap(int i, int j, const PBC_2& pbc) const {
  return pbc.nearest_dis_square(pos[i], pos[j]) <= sigma_square;
}

inline bool HardSphere::check_overlap(int i, const Vec_2f &a, const PBC_2 &pbc) const {
  return pbc.nearest_dis_square(pos[i], a) <= sigma_square;
}

inline void HardSphere::trans_move(const Vec_2f& old_pos, Vec_2f &new_pos,
  double delta_trans, Ran* myran, const PBC_2& pbc) const {
  new_pos.x = old_pos.x + (myran->doub() - 0.5) * delta_trans;
  new_pos.y = old_pos.y + (myran->doub() - 0.5) * delta_trans;
  pbc.wrap(new_pos);
}

inline void HardSphere::trans_move(int i, Vec_2f &new_pos, double delta_trans,
  Ran *myran, const PBC_2& pbc) const {
  new_pos.x = pos[i].x + (myran->doub() - 0.5) * delta_trans;
  new_pos.y = pos[i].y + (myran->doub() - 0.5) * delta_trans;
  pbc.wrap(new_pos);

}

inline void HardSphere::ori_move(const Vec_2f &old_ori, Vec_2f &new_ori,
  double delta_ori, Ran *myran) const {
  circle_point_picking(new_ori.x, new_ori.y, *myran);
  new_ori = old_ori + new_ori * delta_ori;
  new_ori.normalize();
}

class HardJanus: public HardSphere {
public:
  HardJanus(double diameter, double lambda, double theta_m);

  void create_rand(int nPar, Ran *myran, const PBC_2 & pbc);

  void create_rand_ori(int nPar, Ran *myran);

  void ori_move(int i, Vec_2f &new_ori, double delta_ori, Ran *myran) const;

  bool trivial_move(double delta_trans, double delta_ori,
    double eps1, double eps2, Ran *myran, const PBC_2& pbc);

  template <class Func>
  bool trivial_move(double delta_trans, double delta_ori,
    Ran* myran, const PBC_2& pbc, Func func_enegy);

  template <class Func>
  int mc_move(double delta_trans, double delta_ori, Ran *myran,
    const PBC_2 &pbc, Func func_enegy);

  bool cal_energy(double &u, const Vec_2f &r_ij, const Vec_2f &ori_i, const Vec_2f &ori_j,
    double uAA, double uBB) const;

  bool cal_energy(double &u, const Vec_2f &r_ij,
    const Vec_2f &ori_i, const Vec_2f &ori_j, double uAA) const;

  double get_theta_m() const { return theta_max; }

  std::vector<Vec_2f> ori;

protected:
  double sigma_ext_square;
  double cos_theta_max;
  double theta_max;
};

inline void HardJanus::ori_move(int i, Vec_2f &new_ori, double delta_ori,
  Ran *myran) const {
  circle_point_picking(new_ori.x, new_ori.y, *myran);
  new_ori = ori[i] + new_ori * delta_ori;
  new_ori.normalize();

}
template <class Func>
bool HardJanus::trivial_move(double delta_trans, double delta_ori,
  Ran *myran, const PBC_2& pbc, Func func_enegy) {
  int j = int(myran->doub() * pos.size());
  double u_cur = 0;
  for (int i = 0; i < pos.size(); i++) {
    if (i != j) {
      Vec_2f r_ij(pos[j] - pos[i]);
      pbc.nearest_dis(r_ij);
      func_enegy(u_cur, r_ij, ori[i], ori[j]);
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
      if (!func_enegy(u_next, r_ij, ori[i], new_ori)) {
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

template <class Func>
int HardJanus::mc_move(double delta_trans, double delta_ori,
  Ran *myran, const PBC_2 &pbc, Func func_enegy) {
  int valid_move_count = 0;
  int trivial_move_count = 0;
  while (valid_move_count < pos.size()) {
    if (trivial_move(delta_trans, delta_ori, myran, pbc, func_enegy)) {
      valid_move_count++;
    }
    trivial_move_count++;
  }
  return trivial_move_count;
}

class HardABA : public HardJanus {
public:
  HardABA(double diameter, double lambda, double theta_m);

  bool cal_energy(double &u, const Vec_2f& r_ij,
    const Vec_2f& ori_i, const Vec_2f& ori_j, double eps, bool flagAA=true) const;

  bool cal_energy(double &u, const Vec_2f& r_ij,
    const Vec_2f& ori_i, const Vec_2f& ori_j, double uAA, double uBB, double uAB) const;
};

class HardMix {
public:
  HardMix(double sigma_L, double sigma_S,
    double theta_m_L, double theta_m_S, double lambda);
  ~HardMix();

  void create_rand(double nL, double nS, Ran *myran, const PBC_2 & pbc);

  template <class Func>
  bool trivial_move(double delta_trans, double delta_ori,
    Ran *myran, const PBC_2& pbc, Func f1, Func f2, Func f3);

  bool cal_mix_energy(double &u, const Vec_2f &r_ij, const Vec_2f &ori_i,
    const Vec_2f &ori_j, double u_SA_LA);

  HardJanus * large_sphere;
  HardABA * small_sphere;

private:
  int ntot;
  int n_large;
  int n_small;
  double inter_dis_square;
  double inter_force_range_square;
};

template<class Func>
bool HardMix::trivial_move(double delta_trans, double delta_ori, Ran * myran,
  const PBC_2 & pbc, Func f11, Func f22, Func f12) {
  int j = int(myran->doub() * ntot);
  double u_cur = 0;
  Vec_2f new_pos;
  Vec_2f new_ori;
  double u_next = 0;
  bool flag_overlap = false;
  if (j < n_large()) {
    for (int i = 0; i < n_large; i++) {
      if (i != j) {
        Vec_2f r_ij(large_sphere->pos[j] - large_sphere->pos[i]);
        pbc.nearest_dis(r_ij);
        f11(u_cur, r_ij, large_sphere->ori[i], large_sphere->ori[j]);
      }
    }
    for (int i = 0; i < n_small; i++) {
      Vec_2f r_ij(large_sphere->pos[j] - samll_sphere->pos[i]);
      pbc.nearest_dis(r_ij);
      f12(u_cur, r_ij, small_sphere->ori[i], large_sphere->ori[j]);
    }
    large_shpere->trans_move(j, new_pos, delta_trans, myran, pbc);
    large_shpere->ori_move(j, new_ori, delta_ori, myran);
    for (int i = 0; i < n_large; i++) {
      if (i != j) {
        Vec_2f r_ij(new_pos - large_sphere->pos[i]);
        pbc.nearest_dis(r_ij);
        if (!f11(u_next, r_ij, large_sphere->ori[i], new_ori)) {
          flag_overlap = true;
          break;
        }
      }
    }
    if (!flag_overlap) {
      for (int i = 0; i < n_small; i++) {
        Vec_2f r_ij(new_pos - samll_sphere->pos[i]);
        pbc.nearest_dis(r_ij);
        if (!f12(u_next, r_ij, small_sphere->ori[i], new_ori)) {
          flag_overlap = true;
          break;
        }
      }
    }

  } else {
    int k = j - n_large;
    for (int i = 0; i < n_large; i++) {
      Vec_2f r_ij(small_sphere->pos[k] - large_sphere->pos[i]);
      pbc.nearest_dis(r_ij);
      f12(u_cur, -r_ij, small_sphere->ori[k], large_sphere->ori[j]);
    }
    for (int i = 0; i < n_small; i++) {
      if (i != k) {
        Vec_2f r_ij(small_sphere->pos[k] - small_sphere->pos[i]);
        pbc.nearest_dis(r_ij);
        f22(u_cur, r_ij, small_sphere->ori[i], small_sphere->ori[k]);
      }
    }
    small_shpere->trans_move(k, new_pos, delta_trans, myran, pbc);
    small_shpere->ori_move(k, new_ori, delta_ori, myran);
    for (int i = 0; i < n_large; i++) {
      Vec_2f r_ij(new_pos - large_sphere->pos[i]);
      pbc.nearest_dis(r_ij);
      if (!f12(u_next, -r_ij, new_ori, large_sphere->ori[j])) {
        flag_overlap = true;
        break;
      }
    }
    if (!flag_overlap) {
      for (int i = 0; i < n_small; i++) {
        if (i != k) {
          Vec_2f r_ij(new_pos - small_sphere->pos[i]);
          pbc.nearest_dis(r_ij);
          if (!f22(u_next, r_ij, small_sphere->ori[i], new_ori)) {
            flag_overlap = true;
            break;
          }
        }
      }
    }
  }

  bool flag_accept = false;
  if (!flag_overlap) {
    double du = u_next - u_cur;
    if (du <= 0 || myran->doub() < std::exp(-du)) {
      flag_accept = true;
      if (j >= n_large) {
        small_shpere->pos[j - n_large] = new_pos;
        small_shpere->ori[j - n_large] = new_ori;
      } else {
        large_shere->pos[j] = new_pos;
        large_shere->ori[j] = new_ori;
      }
    }
  }
  return flag_accept;
}


#endif
