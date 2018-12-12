#ifndef OUTPUT_H
#define OUTPUT_H
#include <fstream>
#include <chrono>
#include <ctime>
#include "particle_2.h"

#ifdef _MSC_VER
const std::string delimiter("\\");
#else
const std::string delimiter("/");
#endif

// create folder
void mkdir(const char *folder);

// initializing snapshot outputing
void ini_snap_output(double phi, double Lx, double Ly);

// output snapshot
void output_snap(const HardSphere &hs, int i_step);

#endif
