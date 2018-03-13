#include "ioData.h"

using namespace std;

ofstream fout_snap;
char snap_file[100];
Vec_2<double> system_size;

void set_snap(double phi, double Lx, double Ly) {
  mkdir("snap");
  snprintf(snap_file, 100, "snap%s%g_%.2f.extxyz", delimiter.c_str(), Lx, phi);
  fout_snap.open(snap_file);
  system_size.x = Lx;
  system_size.y = Ly;
}

void set_Janus_snap(double phi, double Lx, double Ly, double theta_max,
  double uAA, double uBB) {
  mkdir("snap_Janus");
  snprintf(snap_file, 100, "snap_Janus%s%g_%.2f_%g_%.2f_%.2f.extxyz",
    delimiter.c_str(), Lx, phi, theta_max / PI * 180, uAA, uBB);
  fout_snap.open(snap_file);
  system_size.x = Lx;
  system_size.y = Ly;
}

void set_ABA_snap(double phi, double Lx, double Ly, double theta_max,
  double uAA, double uBB) {
  mkdir("snap_ABA");
  snprintf(snap_file, 100, "snap_ABA%s%g_%.2f_%g_%.2f_%.2f.extxyz",
    delimiter.c_str(), Lx, phi, theta_max / PI * 180, uAA, uBB);
  fout_snap.open(snap_file);
  system_size.x = Lx;
  system_size.y = Ly;
}

void set_mix_snap(double phi, double Lx, double Ly, double theta_m1,
  double theta_m2, double uSA_LA, double uSB_SB, double uLB_LB) {
  mkdir("snap_mix");
  snprintf(snap_file, 100, "snap_ABA%s%g_%.2f_%g_%g_%g_%g_%g.extxyz",
    delimiter.c_str(), Lx, phi, theta_m1 / PI * 180, theta_m2 / PI * 180,
    uSA_LA, uSB_SB, uLB_LB);
  fout_snap.open(snap_file);
  system_size.x = Lx;
  system_size.y = Ly;

}

void output_snap(const HardSphere &hs, int i_step) {
  fout_snap << hs.get_n() << "\n";
  // comment line
  fout_snap << "Lattice=\"" << system_size.x << " 0 0 0 " <<
    system_size.y << " 0 0 0 1\" "
    << "Properties=species:S:1:pos:R:2 "
    << "Time=" << i_step;
  // particles
  for (auto it = hs.pos.cbegin(); it != hs.pos.cend(); ++it) {
    fout_snap << "\n" << std::fixed << std::setprecision(3) << "N\t"
      << (*it).x << "\t" << (*it).y;
  }
  fout_snap << endl;
}

void output_snap(const HardJanus &hj, int i_step, int mode) {
  int n = hj.get_n();
  if (mode == 0) {
    fout_snap << n << "\n";
    // comment line
    fout_snap << "Lattice=\"" << system_size.x << " 0 0 0 " <<
      system_size.y << " 0 0 0 1\" "
      << "Properties=species:S:1:pos:R:2:ori:R:2 "
      << "Time=" << i_step;
    // particles
    for (int i = 0; i < n; i++) {
      //double theta = std::atan2(hj.ori[i].y, hj.ori[i].x);
      fout_snap << "\n" << std::fixed << std::setprecision(3) << "A\t"
        << hj.pos[i].x << "\t" << hj.pos[i].y << "\t"
        << hj.ori[i].x << "\t" << hj.ori[i].y;
    }
  } else if (mode == 1) {
    fout_snap << n << "\n";
    // comment line
    fout_snap << "Lattice=\"" << system_size.x << " 0 0 0 " <<
      system_size.y << " 0 0 0 1\" "
      << "Properties=species:S:1:pos:R:2:ori:R:1 "
      << "Time=" << i_step;
    // particles
    for (int i = 0; i < n; i++) {
      double theta = std::atan2(hj.ori[i].y, hj.ori[i].x);
      fout_snap << "\n" << std::fixed << std::setprecision(3) << "A\t"
        << hj.pos[i].x << "\t" << hj.pos[i].y << "\t"
        << theta;
    }
  } else if (mode == 2) {
    fout_snap << n * 2 << "\n";
    // comment line
    fout_snap << "Lattice=\"" << system_size.x << " 0 0 0 " <<
      system_size.y << " 0 0 0 1\" "
      << "Properties=species:S:1:pos:R:2 "
      << "Time=" << i_step;
    // particles
    for (int i = 0; i < n; i++) {
      double x = hj.pos[i].x + 0.02 * hj.ori[i].x;
      double y = hj.pos[i].y + 0.02 * hj.ori[i].y;
      fout_snap << "\n" << std::fixed << std::setprecision(3) << "A\t"
        << hj.pos[i].x << "\t" << hj.pos[i].y;
      fout_snap << "\nB\t" << std::fixed << std::setprecision(3) << x << "\t" << y;
    }
  }
  
  fout_snap << endl;
}

void output_snap(const HardABA & aba, int i_step, int mode) {
  int n = aba.get_n();
  if (mode == 0) {
    fout_snap << n * 3 << "\n";
    // comment line
    fout_snap << "Lattice=\"" << system_size.x << " 0 0 0 " <<
      system_size.y << " 0 0 0 1\" "
      << "Properties=species:S:1:pos:R:2 "
      << "Time=" << i_step;
    // particles
    double sigma = aba.get_sigma();
    double a = 0.5 * sigma;
    double b = 0.98 * a;
    double theta_m = aba.get_theta_m();
    double sin_theta = std::sin(theta_m);
    double dis = a * std::cos(theta_m)
      - std::sqrt(b * b - a * a * sin_theta * sin_theta);

    for (int i = 0; i < n; i++) {
      double x1 = aba.pos[i].x + dis * aba.ori[i].x;
      double y1 = aba.pos[i].y + dis * aba.ori[i].y;
      double x2 = aba.pos[i].x - dis * aba.ori[i].x;
      double y2 = aba.pos[i].y - dis * aba.ori[i].y;
      fout_snap << "\n" << std::fixed << std::setprecision(3) << "B\t"
        << aba.pos[i].x << "\t" << aba.pos[i].y << "\nA\t" <<
        x1 << "\t" << y1 << "\nA\t" << x2 << "\t" << y2;
    }
  }
  fout_snap << endl;
}
