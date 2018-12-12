#include "output.h"
#include <iomanip>
#ifdef _MSC_VER
#include <io.h>
#else
#include <unistd.h>
#endif

using namespace std;

ofstream fout_snap;
char snap_file[100];
Vec_2<double> system_size;

void mkdir(const char *folder) {
#ifdef _MSC_VER
  if (_access(folder, 0) != 0)
#else
  if (access(folder, 0) != 0)
#endif
  {
    char command[100];
    snprintf(command, 100, "mkdir %s", folder);
    if (system(command))
      cout << "create folder: " << folder << " successfully" << endl;
  } else
    cout << "folder: " << folder << " already exists" << endl;
}

void ini_snap_output(double phi, double Lx, double Ly) {
  mkdir("snap");
  snprintf(snap_file, 100, "snap%s%g_%.2f.extxyz", delimiter.c_str(), Lx, phi);
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
