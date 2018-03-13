#ifndef IODATA_H
#define IODATA_H
#include <fstream>
#include <chrono>
#include <ctime>
#include "comn.h"
#include "particle_2.h"

void set_snap(double phi, double Lx, double Ly);

void set_Janus_snap(double phi, double Lx, double Ly, double theta_max,
  double uAA, double uBB);

void set_ABA_snap(double phi, double Lx, double Ly, double theta_max,
  double uAA, double uBB);

void set_mix_snap(double phi, double Lx, double Ly, double theta_max1,
  double theta_max2, double uSA_LA, double uSB_SB, double uLB_LB);

void output_snap(const HardSphere &hs, int i_step);

void output_snap(const HardJanus &js, int i_step, int mode = 0);

void output_snap(const HardABA &aba, int i_step, int mode = 0);

#endif
