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

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "fix_ext_alignment.h"
#include "atom.h"
#include "atom_masks.h"
#include "update.h"
#include "modify.h"
#include "domain.h"
#include "region.h"
#include "respa.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"
#include "force.h"
#include "atom_vec_ellipsoid.h"
#include "neighbor.h"
#include "comm.h"
#include "math_const.h"
#include "math_vector.h"
#include <iostream>
#include <typeinfo>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;
using namespace std;

enum{NONE,CONSTANT,EQUAL,ATOM,ANGLE_F,ANGLE_V,ANGLE_U};
#define SMALL 0.001
/* ---------------------------------------------------------------------- */

FixExtAlignment::FixExtAlignment(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  xstr(NULL), ystr(NULL), zstr(NULL), estr(NULL), idregion(NULL), sforce(NULL)


{
  if (narg < 8) error->all(FLERR,"Illegal fix extalignment command");

  dynamic_group_allow = 1;
  scalar_flag = 1;
  vector_flag = 1;
  size_vector = 3;
  global_freq = 1;
  extscalar = 1;
  extvector = 1;
  respa_level_support = 1;
  ilevel_respa = 0;
  virial_flag = 1;

  // cosphi = force->cosphi;
  xstr = ystr = zstr = NULL;

  // int style_one;
  if (strcmp(arg[3],"angle_f") == 0)
      style_one = force->numeric(FLERR,"4");
  else if (strcmp(arg[3],"angle_v") == 0)
      style_one = force->numeric(FLERR,"5");
  else if (strcmp(arg[3],"angle_u") == 0)
      style_one = force->numeric(FLERR,"6");
  else
    error->all(FLERR,"Illegal extalignment style");

  // std::cout << "diego";
  // std::cout << style_one;
  // std::cout << "diego";

  if (1 == 1) {
    xvalue = force->numeric(FLERR,arg[4]);
    xstyle = CONSTANT;
  }
  if (1 == 1) {
    yvalue = force->numeric(FLERR,arg[5]);
    ystyle = CONSTANT;
  }
  if (1 == 1) {
    zvalue = force->numeric(FLERR,arg[6]);
    zstyle = CONSTANT;
  }
  if (1 == 1) {
    kalign = force->numeric(FLERR,arg[7]);
  }

  // std::cout << "diego";
  // std::cout << style_one;
  // std::cout << "diego";

  // optional args

  nevery = 1;
  iregion = -1;

  int iarg = 6;

  force_flag = 0;
  foriginal[0] = foriginal[1] = foriginal[2] = foriginal[3] = 0.0;

  maxatom = 1;
  memory->create(sforce,maxatom,4,"addforce:sforce");
}

/* ---------------------------------------------------------------------- */

FixExtAlignment::~FixExtAlignment()
{
  delete [] xstr;
  delete [] ystr;
  delete [] zstr;
  delete [] estr;
  delete [] idregion;
  memory->destroy(sforce);
  // memory->destroy(xvalue);
  // memory->destroy(yvalue);
  // memory->destroy(zvalue);
  // memory->destroy(style_one);
  // memory->destroy(kalign);
}

/* ---------------------------------------------------------------------- */

int FixExtAlignment::setmask()
{
  datamask_read = datamask_modify = 0;

  int mask = 0;
  mask |= POST_FORCE;
  mask |= THERMO_ENERGY;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixExtAlignment::init()
{
  // check variables

  // set index and check validity of region

  if (xstyle == CONSTANT || ystyle == CONSTANT || zstyle == CONSTANT)
    varflag = CONSTANT;

  if (strstr(update->integrate_style,"respa")) {
    ilevel_respa = ((Respa *) update->integrate)->nlevels-1;
    if (respa_level >= 0) ilevel_respa = MIN(respa_level,ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixExtAlignment::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet"))
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(ilevel_respa);
    post_force_respa(vflag,ilevel_respa,0);
    ((Respa *) update->integrate)->copy_f_flevel(ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixExtAlignment::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixExtAlignment::post_force(int vflag)
{
  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  imageint *image = atom->image;
  double v[6];
  int nlocal = atom->nlocal;

  double **torque = atom->torque;
  int *ellipsoid = atom->ellipsoid;
  AtomVecEllipsoid *avec = (AtomVecEllipsoid *) atom->style_match("ellipsoid");
  AtomVecEllipsoid::Bonus *bonus = avec->bonus;
  int newton_bond = force->newton_bond;

  double u0[3],uI[3];
  double values[3], valuesnorm,directionfield[3],uIcrossdf[3],uIdotdf,f1[3];
  double fterm1[3],fterm2[3];
  double *quatI;
  double quatItmp[4];
  double U,theta1, theta2, phi,dtheta1, dtheta2, dphi;
  double inverse;
  double a0,a1,a2,normuI,normdf, dUda0, dUda1, dUda2,normdirectionfield;
  double dUdfI[3], dUdfJ[3];
  double torq[3], torqI[3], torqJ[3];
  double GI[3];
  //int    style_one;

  // std::cout << "diego";
  // std::cout << style_one;
  // std::cout << "diego";

  // values[0] = xvalue; values[1] = yvalue; values[2] = zvalue;
  // valuesnorm = sqrt(xvalue*xvalue+yvalue*yvalue+zvalue*zvalue);


  if (update->ntimestep % nevery) return;

  // energy and virial setup

  if (vflag) v_setup(vflag);
  else evflag = 0;

  if (lmp->kokkos)
    atom->sync_modify(Host, (unsigned int) (F_MASK | MASK_MASK),
                      (unsigned int) F_MASK);

  // update region if necessary

  Region *region = NULL;
  if (iregion >= 0) {
    region = domain->regions[iregion];
    region->prematch();
  }

  // foriginal[0] = "potential energy" for added force
  // foriginal[123] = force on atoms before extra force added

  foriginal[0] = foriginal[1] = foriginal[2] = foriginal[3] = 0.0;
  force_flag = 0;

  // constant force
  // potential energy = - x dot f in unwrapped coords

  if (varflag == CONSTANT) {
    double unwrap[3];
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        if (region && !region->match(x[i][0],x[i][1],x[i][2])) continue;
        domain->unmap(x[i],image[i],unwrap);

        if (style_one == ANGLE_F){
          u0[0] = 1; u0[1] = 0; u0[2] = 0; //f
        }
        else if (style_one == ANGLE_V){
          u0[0] = 0; u0[1] = 1; u0[2] = 0; //v
        }
        else if (style_one == ANGLE_U){
          u0[0] = 0; u0[1] = 0; u0[2] = 1; //u
        }


        quatI = bonus[ellipsoid[i]].quat;
        quatItmp[0] = quatI[0]; quatItmp[1] = quatI[1]; quatItmp[2] = quatI[2]; quatItmp[3] = quatI[3];

        LAMMPS_NS::quat_vec_rot(uI,u0,quatItmp);
        LAMMPS_NS::vec_norm(uI);

        uIdotdf = xvalue*uI[0] + yvalue*uI[1] + zvalue*uI[2];

        if (uIdotdf > 0){
            directionfield[0] = xvalue;
            directionfield[1] = yvalue;
            directionfield[2] = zvalue;
            LAMMPS_NS::vec_norm(directionfield);
        }

        else{
            directionfield[0] = -1*xvalue;
            directionfield[1] = -1*yvalue;
            directionfield[2] = -1*zvalue;
            LAMMPS_NS::vec_norm(directionfield);
        }

        uIcrossdf[0] = uI[1]*directionfield[2] - uI[2]*directionfield[1];
        uIcrossdf[1] = uI[2]*directionfield[0] - uI[0]*directionfield[2];
        uIcrossdf[2] = uI[0]*directionfield[1] - uI[1]*directionfield[0];

        torque[i][0] += kalign*uIcrossdf[0];
        torque[i][1] += kalign*uIcrossdf[1];
        torque[i][2] += kalign*uIcrossdf[2];

        // uIdotdf = uI[0]*directionfield[0] + uI[1]*directionfield[1] + uI[2]*directionfield[2];
        //
        // if (uIdotdf > 1.0) uIdotdf = 1.0;
        // if (uIdotdf < -1.0) uIdotdf = -1.0;
        //
        // GI[0] = ( uIdotdf*directionfield[0] - uI[0] );
        // GI[1] = ( uIdotdf*directionfield[1] - uI[1] );
        // GI[2] = ( uIdotdf*directionfield[2] - uI[2] );
        //
        // f1[0] = zvalue * GI[0] ;
        // f1[1] = zvalue * GI[1] ;
        // f1[2] = zvalue * GI[2] ;

        // if (evflag) {
        //   v[0] = xvalue * unwrap[0]*cosphi*cosphi;
        //   v[1] = yvalue * unwrap[1]*cosphi*cosphi;
        //   v[2] = zvalue * unwrap[2]*cosphi*cosphi;
        //   v[3] = xvalue * unwrap[1]*cosphi*cosphi;
        //   v[4] = xvalue * unwrap[2]*cosphi*cosphi;
        //   v[5] = yvalue * unwrap[2]*cosphi*cosphi;
        //   v_tally(i,v);
        // }
      }

  }
}

/* ---------------------------------------------------------------------- */

void FixExtAlignment::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == ilevel_respa) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixExtAlignment::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   potential energy of added force
------------------------------------------------------------------------- */

double FixExtAlignment::compute_scalar()
{
  // only sum across procs one time

  if (force_flag == 0) {
    MPI_Allreduce(foriginal,foriginal_all,4,MPI_DOUBLE,MPI_SUM,world);
    force_flag = 1;
  }
  return foriginal_all[0];
}

/* ----------------------------------------------------------------------
   return components of total force on fix group before force was changed
------------------------------------------------------------------------- */

double FixExtAlignment::compute_vector(int n)
{
  // only sum across procs one time

  if (force_flag == 0) {
    MPI_Allreduce(foriginal,foriginal_all,4,MPI_DOUBLE,MPI_SUM,world);
    force_flag = 1;
  }
  return foriginal_all[n+1];
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */
double FixExtAlignment::memory_usage()
{
  double bytes = 0.0;
  if (varflag == ATOM) bytes = maxatom*4 * sizeof(double);
  return bytes;
}
