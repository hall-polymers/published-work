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

/* ----------------------------------------------------------------------
   Contributing author: Chris A Brackley (U Edinburgh) 
   Based on the simulation scheme described in 
       G Chirico and J Langowski, Biopolymers 34 p415-433 (1994)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "angle_orient.h"
#include "atom.h"
#include "neighbor.h"
#include "domain.h"
#include "comm.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "math_vector.h"

#include "atom_vec_ellipsoid.h"

using namespace LAMMPS_NS;
using namespace MathConst;

enum{ANGLE_F,ANGLE_V,ANGLE_U};
#define SMALL 0.001

/* ---------------------------------------------------------------------- */

AngleOrient::AngleOrient(LAMMPS *lmp) : Angle(lmp) {}

/* ---------------------------------------------------------------------- */

AngleOrient::~AngleOrient()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(kspr1);
    memory->destroy(kspr2);
    memory->destroy(kspr0);
    memory->destroy(theta1_0);
    memory->destroy(theta2_0);
    memory->destroy(phi0);
    memory->destroy(style);
  }
}

/* ---------------------------------------------------------------------- */

void AngleOrient::compute(int eflag, int vflag)
{
  int i,j,iDummy,n,m;
  double eangle,force1[3],force3[3],dx,dy,dz;


  eangle = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = 0;

  double **x = atom->x; // position vector
  double **f = atom->f; // force vector
  double **torque = atom->torque;
  int *ellipsoid = atom->ellipsoid;
  AtomVecEllipsoid *avec = (AtomVecEllipsoid *) atom->style_match("ellipsoid");
  AtomVecEllipsoid::Bonus *bonus = avec->bonus;
  int **anglelist = neighbor->anglelist;
  int nanglelist = neighbor->nanglelist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  double f0[3],fI[3],fJ[3],r[3], rn;
  double fterm1[3],fterm2[3];
  double *quatI, *quatJ;
  double quatItmp[4], quatJtmp[4];
  double U,theta1, theta2, phi,dtheta1, dtheta2, dphi;
  double inverse;
  double a0,a1,a2, dUda0, dUda1, dUda2;
  double dUdfI[3], dUdfJ[3];
  double torq[3], torqI[3], torqJ[3];
  
  // I think I fixed this to work with newton_bond = false
  //if (!newton_bond)
  //  error->all(FLERR,"'newton' flag for bonded interactions must be 'on'");

  for (n = 0; n < nanglelist; n++) {

    i = anglelist[n][0]; 
    j = anglelist[n][1]; 
    iDummy = anglelist[n][2]; 
    m = anglelist[n][3];

    if (style[m] == ANGLE_F){
      f0[0] = 1; f0[1] = 0; f0[2] = 0; //f
    }
    else if (style[m] == ANGLE_V){
      f0[0] = 0; f0[1] = 1; f0[2] = 0; //v
    }
    else if (style[m] == ANGLE_U){
      f0[0] = 0; f0[1] = 0; f0[2] = 1; //u
    }
    


    //get rIJ
    dx = x[i][0] - x[j][0]; // get vector between i and j
    dy = x[i][1] - x[j][1];
    dz = x[i][2] - x[j][2];
    r[0] = dx;
    r[1] = dy;
    r[2] = dz;
    rn = sqrt(dx*dx+dy*dy+dz*dz);


    quatI = bonus[ellipsoid[i]].quat;
    quatItmp[0] = quatI[0]; quatItmp[1] = quatI[1]; quatItmp[2] = quatI[2]; quatItmp[3] = quatI[3];
    quatJ = bonus[ellipsoid[j]].quat;
    quatJtmp[0] = quatJ[0]; quatJtmp[1] = quatJ[1]; quatJtmp[2] = quatJ[2]; quatJtmp[3] = quatJ[3];

    LAMMPS_NS::quat_vec_rot(fI,f0,quatItmp);
    LAMMPS_NS::quat_vec_rot(fJ,f0,quatJtmp);



    // calc (and condition) a0, a1, a2
    a0 = LAMMPS_NS::vec_dot(fI,fJ);
    a1 = LAMMPS_NS::vec_dot(fI,r)/rn;
    a2 = LAMMPS_NS::vec_dot(fJ,r)/rn;
    if (a0 > 1) a0 = 1;
    if (a0 < -1) a0 = -1;
    if (a1 > 1) a1 = 1;
    if (a1 < -1) a1 = -1;
    if (a2 > 1) a2 = 1;
    if (a2 < -1) a2 = -1;

    // calculate thetas and phi 
    phi = acos(a0);
    dphi = phi - phi0[m];
    
    // note that thetas are pi minus the input value
    // this is since rij is defined as ri - rj (instead of the typical rj - ri)
    theta1 = acos(a1);
    dtheta1 = theta1 - theta1_0[m];

    theta2 = acos(a2);
    dtheta2 = theta2 - theta2_0[m];

    // calc energy
    U = 0.5*(kspr1[m]*dtheta1*dtheta1 + kspr2[m]*dtheta2*dtheta2 + kspr0[m]*dphi*dphi);
    
    // calc partial derivatives
    inverse = sqrt(1-a0*a0);
    if (inverse < SMALL) inverse = SMALL;
    dUda0 = - kspr0[m]*dphi / inverse;

    inverse = sqrt(1-a1*a1);
    if (inverse < SMALL) inverse = SMALL;
    dUda1 = - kspr1[m]*dtheta1 / inverse;

    inverse = sqrt(1-a2*a2);
    if (inverse < SMALL) inverse = SMALL;
    dUda2 = - kspr2[m]*dtheta2 / inverse;

    // forces 
    fterm1[0] = -(fI[0]/rn - r[0]*a1/rn/rn) * dUda1;
    fterm1[1] = -(fI[1]/rn - r[1]*a1/rn/rn) * dUda1;
    fterm1[2] = -(fI[2]/rn - r[2]*a1/rn/rn) * dUda1;

    fterm2[0] = -(fJ[0]/rn - r[0]*a2/rn/rn) * dUda2;
    fterm2[1] = -(fJ[1]/rn - r[1]*a2/rn/rn) * dUda2;
    fterm2[2] = -(fJ[2]/rn - r[2]*a2/rn/rn) * dUda2;

    force1[0] = fterm1[0] + fterm2[0];
    force1[1] = fterm1[1] + fterm2[1];
    force1[2] = fterm1[2] + fterm2[2];

    if (newton_bond || i < nlocal) {
      f[i][0] += force1[0];
      f[i][1] += force1[1];
      f[i][2] += force1[2];
    }
    
    if (newton_bond || j < nlocal) {
      f[j][0] -= force1[0];
      f[j][1] -= force1[1];
      f[j][2] -= force1[2];
    }

    //torque (see Allen and Tildsley Apendix C)
    if (newton_bond || i < nlocal) {
      dUdfI[0] = r[0] / rn * dUda1 + fJ[0]*dUda0;
      dUdfI[1] = r[1] / rn * dUda1 + fJ[1]*dUda0;
      dUdfI[2] = r[2] / rn * dUda1 + fJ[2]*dUda0;

      torqI[0] = -(fI[1]*dUdfI[2] - fI[2]*dUdfI[1]);
      torqI[1] = -(fI[2]*dUdfI[0] - fI[0]*dUdfI[2]);
      torqI[2] = -(fI[0]*dUdfI[1] - fI[1]*dUdfI[0]);

      torque[i][0] += torqI[0];
      torque[i][1] += torqI[1];
      torque[i][2] += torqI[2];
    }

    if (newton_bond || j < nlocal) {
      dUdfJ[0] = r[0] / rn * dUda2 + fI[0]*dUda0;
      dUdfJ[1] = r[1] / rn * dUda2 + fI[1]*dUda0;
      dUdfJ[2] = r[2] / rn * dUda2 + fI[2]*dUda0;

      torqJ[0] = -(fJ[1]*dUdfJ[2] - fJ[2]*dUdfJ[1]);
      torqJ[1] = -(fJ[2]*dUdfJ[0] - fJ[0]*dUdfJ[2]);
      torqJ[2] = -(fJ[0]*dUdfJ[1] - fJ[1]*dUdfJ[0]);

      torque[j][0] += torqJ[0];
      torque[j][1] += torqJ[1];
      torque[j][2] += torqJ[2];
    }

    
    force3[0] = force3[1] = force3[2]  = 0.0;  // for virial calculation

    if (eflag) eangle = U;

    if (evflag) // tally energy (virial=0 because force=0)
      ev_tally(i,j,iDummy,nlocal,newton_bond,eangle,force1,force3,
               dx,dy,dz,0.0,0.0,0.0);

  }
}

/* ---------------------------------------------------------------------- */

void AngleOrient::allocate()
{
  allocated = 1;
  int n = atom->nangletypes;

  memory->create(kspr0,n+1,"angle:kspr");
  memory->create(kspr1,n+1,"angle:kspr");
  memory->create(kspr2,n+1,"angle:kspr");
  memory->create(theta1_0,n+1,"angle:theta1_0");
  memory->create(theta2_0,n+1,"angle:theta2_0");
  memory->create(phi0,n+1,"angle:phi0");
  memory->create(style,n+1,"angle:style");

  memory->create(setflag,n+1,"angle:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

void AngleOrient::coeff(int narg, char **arg)
{
  if (narg != 8) error->all(FLERR,"Incorrect args for angle coefficients");
  if (!allocated) allocate();

  int ilo,ihi;
  force->bounds(FLERR,arg[0],atom->nangletypes,ilo,ihi);

  int    style_one;
  if (strcmp(arg[1],"angle_f") == 0)
      style_one = ANGLE_F;
  else if (strcmp(arg[1],"angle_v") == 0)
      style_one = ANGLE_V;
  else if (strcmp(arg[1],"angle_u") == 0)
      style_one = ANGLE_U;
  else
    error->all(FLERR,"Illegal angle orient style");

  double kspr1_one = force->numeric(FLERR,arg[2]);
  double kspr2_one = force->numeric(FLERR,arg[3]);
  double kspr0_one = force->numeric(FLERR,arg[4]);
  //subtract from pi, this is since rij is defined from i to j, 
  double theta1_0_one = 180 - force->numeric(FLERR,arg[5]);
  double theta2_0_one = 180 - force->numeric(FLERR,arg[6]);
  double phi0_one = force->numeric(FLERR,arg[7]);


  theta1_0_one *= MY_PI / 180.0; //convert to radians
  theta2_0_one *= MY_PI / 180.0;
  phi0_one *= MY_PI / 180.0;


  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    kspr1[i] = kspr1_one;
    kspr2[i] = kspr2_one;
    kspr0[i] = kspr0_one;
    theta1_0[i] = theta1_0_one;
    theta2_0[i] = theta2_0_one;
    phi0[i] = phi0_one;
    style[i] = style_one;
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR,"Incorrect args for angle coefficients");
}

/* ----------------------------------------------------------------------
   used by SHAKE
------------------------------------------------------------------------- */

double AngleOrient::equilibrium_angle(int i)
{
  return 0.0;
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void AngleOrient::write_restart(FILE *fp)
{
  fwrite(&style[1],sizeof(int),atom->nangletypes,fp);
  fwrite(&kspr1[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&kspr2[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&kspr0[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&theta1_0[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&theta2_0[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&phi0[1],sizeof(double),atom->nangletypes,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void AngleOrient::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    fread(&style[1],sizeof(int),atom->nangletypes,fp);
    fread(&kspr1[1],sizeof(double),atom->nangletypes,fp);
    fread(&kspr2[1],sizeof(double),atom->nangletypes,fp);
    fread(&kspr0[1],sizeof(double),atom->nangletypes,fp);
    fread(&theta1_0[1],sizeof(double),atom->nangletypes,fp);
    fread(&theta2_0[1],sizeof(double),atom->nangletypes,fp);
    fread(&phi0[1],sizeof(double),atom->nangletypes,fp);
  }
  MPI_Bcast(&style[1],atom->nangletypes,MPI_INT,0,world);
  MPI_Bcast(&kspr1[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&kspr2[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&kspr0[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&theta1_0[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&theta2_0[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&phi0[1],atom->nangletypes,MPI_DOUBLE,0,world);

  for (int i = 1; i <= atom->nangletypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   used by ComputeAngleLocal
------------------------------------------------------------------------- */

double AngleOrient::single(int type, int i, int j, int iDummy)
{
  
  error->all(FLERR,"not implemented! should just copy from 'compute' above");
  //double eangle,f1[3],f3[3],dx,dy,dz;

  //int *ellipsoid = atom->ellipsoid;
  //AtomVecEllipsoid *avec = (AtomVecEllipsoid *) atom->style_match("ellipsoid");
  //AtomVecEllipsoid::Bonus *bonus = avec->bonus;

  //double u0[3],uI[3],uJ[3];
  //double *quatI, *quatJ;
  //double quatItmp[4], quatJtmp[4];
  //double U,theta,dtheta, dot;
  //double dUdtheta;
  //double invdarccosddot, darccosddot;

  //int m = type;

  //if (style[m] == ANGLE_F){
  //  u0[0] = 1; u0[1] = 0; u0[2] = 0; //f
  //}
  //else if (style[m] == ANGLE_V){
  //  u0[0] = 0; u0[1] = 1; u0[2] = 0; //v
  //}
  //else if (style[m] == ANGLE_U){
  //  u0[0] = 0; u0[1] = 0; u0[2] = 1; //u
  //}

  //quatI = bonus[ellipsoid[i]].quat;
  //quatItmp[0] = quatI[0]; quatItmp[1] = quatI[1]; quatItmp[2] = quatI[2]; quatItmp[3] = quatI[3];
  //quatJ = bonus[ellipsoid[j]].quat;
  //quatJtmp[0] = quatJ[0]; quatJtmp[1] = quatJ[1]; quatJtmp[2] = quatJ[2]; quatJtmp[3] = quatJ[3];

  //LAMMPS_NS::quat_vec_rot(uI,u0,quatItmp);
  //LAMMPS_NS::quat_vec_rot(uJ,u0,quatJtmp);


  //dot = LAMMPS_NS::vec_dot(uI,uJ);
  //if (dot > 1) dot = 1;
  //if (dot < -1) dot = -1;
  //theta = acos(dot);
  //dtheta = theta - theta0[m];
  //U = 0.5*kspr[m]*dtheta*dtheta;

  //printf("ERROR!! NOT IMPLEMENTED!!!\n");
  //exit(1);

  //dont need torques or forces
  //dUdtheta = kspr[m]*dtheta;
  //invdarccosddot = sqrt(1-dot*dot);
  //if (invdarccosddot < SMALL) invdarccosddot = SMALL;
  //darccosddot = -1.0 / invdarccosddot;

  //// torque I
  ////calc derivatives, its just chain rule
  //double dUduI[3], crossI[3];
  //
  //dUduI[0] = dUdtheta * darccosddot * uJ[0]; //uJ is correct, uI term disappears
  //dUduI[1] = dUdtheta * darccosddot * uJ[1];
  //dUduI[2] = dUdtheta * darccosddot * uJ[2];

  //crossI[0] = uI[1]*dUduI[2] - uI[2]*dUduI[1];
  //crossI[1] = uI[2]*dUduI[0] - uI[0]*dUduI[2];
  //crossI[2] = uI[0]*dUduI[1] - uI[1]*dUduI[0];

  //torque[i][0] += -crossI[0];  
  //torque[i][1] += -crossI[1]; 
  //torque[i][2] += -crossI[2];

  //// torque J
  //double dUduJ[3], crossJ[3];
  //dUduJ[0] = dUdtheta * darccosddot * uI[0]; //uI is correct, uJ term disappears 
  //dUduJ[1] = dUdtheta * darccosddot * uI[1];
  //dUduJ[2] = dUdtheta * darccosddot * uI[2];

  //crossJ[0] = uJ[1]*dUduJ[2] - uJ[2]*dUduJ[1];
  //crossJ[1] = uJ[2]*dUduJ[0] - uJ[0]*dUduJ[2];
  //crossJ[2] = uJ[0]*dUduJ[1] - uJ[1]*dUduJ[0];

  //torque[j][0] += -crossJ[0];  
  //torque[j][1] += -crossJ[1]; 
  //torque[j][2] += -crossJ[2];

  //no forces 
  //dx = x[j][0] - x[i][0]; // get vector between i and j
  //dy = x[j][1] - x[i][1];
  //dz = x[j][2] - x[i][2];

  //f1[0] = f1[1] = f1[2]  = 0.0;  // for virial calculation
  //f3[0] = f3[1] = f3[2]  = 0.0;  // for virial calculation

  //return U; // energy
}
