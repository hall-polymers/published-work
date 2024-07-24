/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef ANGLE_CLASS

AngleStyle(orient,AngleOrient)

#else

#ifndef LMP_ANGLE_ORIENT_H
#define LMP_ANGLE_ORIENT_H

#include "stdio.h"
#include "angle.h"

namespace LAMMPS_NS {

class AngleOrient : public Angle {
 public:
  AngleOrient(class LAMMPS *);
  virtual ~AngleOrient();
  virtual void compute(int, int);
  void coeff(int, char **);
  double equilibrium_angle(int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  double single(int, int, int, int);

 protected:
  double *kspr0, *kspr1, *kspr2,*theta1_0, *theta2_0, *phi0;
  int *style;

  void allocate();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Incorrect args for angle coefficients

Self-explanatory.  Check the input script or data file.

*/
