/**
 * @file boundary.h
 * @author Yu Duan (duanyu100@yeah.net)
 * @brief 
 * @version 0.1
 * @date 2018-12-12
 * 
 * Boundary condition for the simulation domain with length (Lx, Ly).
 * Only the periodic boundary condition is implemented.
 * 
 * @copyright Copyright (c) 2018
 * 
 */
#ifndef BOUNDARY_H
#define BOUNDARY_H
#include "vect.h"

/**
 * class for periodic boundary condition in 2D
 */
class PBC_2 {
public:
  /**
   * @brief Construct a new PBC_2 object
   * 
   */
  PBC_2() : L(), half_L() {}

  /**
   * @brief Construct a new PBC_2 object
   * 
   * @param Lx Domain length in x direction
   * @param Ly Domain length in y direction
   */
  PBC_2(double Lx, double Ly) :
    L(Lx, Ly), half_L(0.5 * Lx, 0.5 * Ly) {}
  
  /**
   * @brief Construct a new PBC_2 object
   * 
   * @param L0  A 2-d vector (Lx, Ly). 
   */
  PBC_2(const Vec_2<double> &L0) :
    L(L0), half_L(0.5 * L0.x, 0.5 * L0.y) {}

  /**
   * @brief Initialize the class
   * 
   * @param Lx  Domain length in x direction
   * @param Ly  Domain length in x direction
   */
  void ini(double Lx, double Ly);

  /**
   * @brief Calculate the effective length of input vector under
   * periodic boundary condition.
   * 
   * @param dis Input vector, also as output vector
   */
  void nearest_dis(Vec_2<double> &dis) const;

  /**
   * @brief Calculate the nearest distance between two particle under periodic
   * boundary condition.
   * 
   * @tparam Par1 Template class for particle 1
   * @tparam Par2 Template class for particle 2
   * @param dis   The distance of two particles
   * @param p1    Particle 1
   * @param p2    Particle 2
   */
  template<class Par1, class Par2>
  void nearest_dis(Vec_2<double> &dis, const Par1 &p1, const Par2 &p2) const;

  /**
   * @brief Calculate the square of nearest distance between two particle under
   * periodic boundary condition
   * 
   * @tparam Par1 
   * @tparam Par2 
   * @param p1 
   * @param p2 
   * @return double 
   */
  template<class Par1, class Par2>
  double nearest_dis_square(const Par1 &p1, const Par2 &p2) const;

  /**
   * @brief Wrap the postion vector of input particle when it's outside of the
   * simulation domain.
   * 
   * @tparam Par  Template class for intput particle 
   * @param R   Input particle
   */
  template <class Par>
  void wrap(Par &R) const;

  // return Lx
  double Lx() const { return L.x; }

  // return Ly
  double Ly() const { return L.y; }

protected:
  Vec_2<double> L;
  Vec_2<double> half_L;
};

inline void PBC_2::ini(double Lx, double Ly) {
  L.x = Lx;
  L.y = Ly;
  half_L.x = 0.5 * Lx;
  half_L.y = 0.5 * Ly;
}

inline void PBC_2::nearest_dis(Vec_2<double> &dis) const {
  if (dis.x < -half_L.x) {
    dis.x += L.x;
  } else if (dis.x > half_L.x) {
    dis.x -= L.x;
  }
  if (dis.y < -half_L.y) {
    dis.y += L.y;
  } else if (dis.y > half_L.y) {
    dis.y -= L.y;
  }
}

template<class Par1, class Par2>
void PBC_2::nearest_dis(Vec_2<double>& dis, const Par1 & p1, const Par2 & p2) const {
  dis.x = p1.x - p2.x;
  dis.y = p1.y - p2.y;
  nearest_dis(dis);
}

template<class Par1, class Par2>
double PBC_2::nearest_dis_square(const Par1 &p1, const Par2 &p2) const {
  Vec_2<double> dR(p1.x - p2.x, p1.y - p2.y);
  nearest_dis(dR);
  return dR.square();
}

template <class Par>
void PBC_2::wrap(Par &R) const {
  if (R.x < 0) {
    R.x += L.x;
  } else if (R.x >= L.x) {
    R.x -= L.x;
  }
  if (R.y < 0) {
    R.y += L.y;
  } else if (R.y >= L.y) {
    R.y -= L.y;
  }
}
#endif
