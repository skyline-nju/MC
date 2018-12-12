#ifndef VECT_H
#define VECT_H
#include <cmath>

/**
 * @brief Class for 2-d vector
 * 
 * @tparam T 
 */
template <typename T>
struct Vec_2 {
  T x, y;

  Vec_2() : x(0), y(0) {}
  //Vec_2(T a) : x(a), y(a) {}
  Vec_2(const Vec_2<T> &a) : x(a.x), y(a.y) {}
  Vec_2(T x0, T y0) : x(x0), y(y0) {}

  void operator = (const Vec_2<T> &a) { x = a.x; y = a.y; }
  void operator += (const Vec_2<T> &a) { x += a.x; y += a.y; }
  void operator -= (const Vec_2<T> &a) { x -= a.x; y -= a.y; }
  void operator *= (double a) { x *= a; y *= a;}

  Vec_2 operator +(double a) const { return Vec_2<T>(x + a, y + a); }
  Vec_2 operator +(const Vec_2<T> &a) const { return Vec_2<T>(x + a.x, y + a.y); }
  friend inline Vec_2 operator +(double a, const Vec_2<T> &b) {
    return Vec_2<T>(a + b.x, a + b.y);
  }

  Vec_2 operator -(double a) const { return Vec_2<T>(x - a, y - a); }
  Vec_2 operator -(const Vec_2<T> &a) const { return Vec_2<T>(x - a.x, y - a.y); }
  friend inline Vec_2 operator -(double a, const Vec_2<T> &b) {
    return Vec_2<T>(a - b.x, a - b.y);
  }

  Vec_2 operator *(double a) const { return Vec_2<T>(x *a, y * a); }
  friend inline Vec_2 operator *(double a, const Vec_2<T> &b) {
    return Vec_2<T>(a * b.x, a * b.y);
  }

  Vec_2 operator /(double a) const { return Vec_2<T>(x / a, y / a); }

  double square() const { return x * x + y * y; }
  double dot(const Vec_2<T> & a) const { return x * a.x + y * a.y; }
  double cross(const Vec_2<T> &a) const { return x * a.y - y * a.x; }
  void normalize() {
    double one_over_r = 1 / std::sqrt(square());
    x *= one_over_r;
    y *= one_over_r;
  }

  void rotate(double theta) {
    double c = std::cos(theta);
    double s = std::sin(theta);
    double x_new = x * c - y * s;
    double y_new = x * s + y * c;
    x = x_new;
    y = y_new;
  }
};

#endif

