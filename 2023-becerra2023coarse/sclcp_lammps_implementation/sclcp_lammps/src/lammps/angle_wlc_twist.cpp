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
#include "angle_wlc_twist.h"
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

/* ---------------------------------------------------------------------- */

AngleWLCTwist::AngleWLCTwist(LAMMPS *lmp) : Angle(lmp) {}

/* ---------------------------------------------------------------------- */

AngleWLCTwist::~AngleWLCTwist()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(kalign);
    memory->destroy(ktwist);
    memory->destroy(twist0);
  }
}

/* ---------------------------------------------------------------------- */

void AngleWLCTwist::compute(int eflag, int vflag)
{
  int i,iP1,iDummy,n,type;
  double eangle,f1[3],f3[3],dx,dy,dz;
  double uI[3],uP1[3],  // axis of the atoms
    fI[3],fP1[3],
    vI[3],vP1[3];

  double *quatI, *quatP1;
  double quatItmp[4], quatP1tmp[4];

  double uIdotuiP1plus1, inv_uIdotuiP1plus1, cosalphaiplusgammai;
  double fIcrossfP1[3], vIcrossvP1[3], uIcrossuP1[3], HI[3], kHI[3];
  double fIcrossvP1[3], vIcrossfP1[3];
  double tI[3], inv_bI, bI2;
  double uIcrosstI[3];
  double f0[3], v0[3], u0[3];
  

  double GI[3],uIdottI;

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

  double PI, costwist0, sintwist0, twist0rad, alphaiplusgammai, cosalphaiplusgammaieq; 
  PI = 3.14159265359;

  if (!newton_bond)
    error->all(FLERR,"'newton' flag for bonded interactions must be 'on'");

  for (n = 0; n < nanglelist; n++) {
    i = anglelist[n][0]; 
    iP1 = anglelist[n][1]; 
    iDummy = anglelist[n][2]; 
    type = anglelist[n][3];

    /* ----------------------------------------------
       Get u,f,v axes of beads i and iP1
       ---------------------------------------------- */
   
    f0[0] = 1; f0[1] = 0; f0[2] = 0;   
    v0[0] = 0; v0[1] = 1; v0[2] = 0;  
    u0[0] = 0; u0[1] = 0; u0[2] = 1; 

    quatI = bonus[ellipsoid[i]].quat;
    quatItmp[0] = quatI[0]; quatItmp[1] = quatI[1]; quatItmp[2] = quatI[2]; quatItmp[3] = quatI[3];
    //fvu for site 1
    LAMMPS_NS::quat_vec_rot(fI,f0,quatItmp);
    LAMMPS_NS::quat_vec_rot(vI,v0,quatItmp);
    LAMMPS_NS::quat_vec_rot(uI,u0,quatItmp);
    LAMMPS_NS::vec_norm(fI);
    LAMMPS_NS::vec_norm(vI);
    LAMMPS_NS::vec_norm(uI);

    quatP1 = bonus[ellipsoid[iP1]].quat;
    quatP1tmp[0] = quatP1[0]; quatP1tmp[1] = quatP1[1]; quatP1tmp[2] = quatP1[2]; quatP1tmp[3] = quatP1[3];
    //fvu for site 2
    LAMMPS_NS::quat_vec_rot(fP1,f0,quatP1tmp);
    LAMMPS_NS::quat_vec_rot(vP1,v0,quatP1tmp);
    LAMMPS_NS::quat_vec_rot(uP1,u0,quatP1tmp);
    LAMMPS_NS::vec_norm(fP1);
    LAMMPS_NS::vec_norm(vP1);
    LAMMPS_NS::vec_norm(uP1);

    /* ----------------------------------------------
       Twist torque
       ---------------------------------------------- */

    twist0rad = twist0[type] * PI / 180.0;
    costwist0 = cos(twist0rad);
    sintwist0 = sin(twist0rad);

    uIdotuiP1plus1 = uI[0]*uP1[0] + uI[1]*uP1[1] + uI[2]*uP1[2] + 1.0;

    //Make sure numerical instabilities are removed
    if (uIdotuiP1plus1 < 1E-3) uIdotuiP1plus1 = 1E-3;
    if (uIdotuiP1plus1 > 2.0) uIdotuiP1plus1 = 2.0;

    inv_uIdotuiP1plus1 = 1.0/uIdotuiP1plus1;

    cosalphaiplusgammai = (fP1[0]*fI[0] + fP1[1]*fI[1] + fP1[2]*fI[2] + 
			   vP1[0]*vI[0] + vP1[1]*vI[1] + vP1[2]*vI[2]) * inv_uIdotuiP1plus1;   			   
			   
    if (cosalphaiplusgammai > 1.0) cosalphaiplusgammai = 1.0;
    if (cosalphaiplusgammai < -1.0) cosalphaiplusgammai = -1.0;
    
    alphaiplusgammai = acos(cosalphaiplusgammai);
    //printf("twits: %E \n", alphaiplusgammai);	
    cosalphaiplusgammaieq = cos(alphaiplusgammai-twist0rad);
    //printf("twits: %E \n", twist0[type]);
    //printf("twits: %E \n", cosalphaiplusgammaieq);

    fIcrossfP1[0] = fI[1]*fP1[2] - fI[2]*fP1[1];
    fIcrossfP1[1] = fI[2]*fP1[0] - fI[0]*fP1[2];
    fIcrossfP1[2] = fI[0]*fP1[1] - fI[1]*fP1[0];
    vIcrossvP1[0] = vI[1]*vP1[2] - vI[2]*vP1[1];
    vIcrossvP1[1] = vI[2]*vP1[0] - vI[0]*vP1[2];
    vIcrossvP1[2] = vI[0]*vP1[1] - vI[1]*vP1[0];
    uIcrossuP1[0] = uI[1]*uP1[2] - uI[2]*uP1[1];
    uIcrossuP1[1] = uI[2]*uP1[0] - uI[0]*uP1[2];
    uIcrossuP1[2] = uI[0]*uP1[1] - uI[1]*uP1[0];
	
    fIcrossvP1[0] = fI[1]*vP1[2] - fI[2]*vP1[1];
    fIcrossvP1[1] = fI[2]*vP1[0] - fI[0]*vP1[2];
    fIcrossvP1[2] = fI[0]*vP1[1] - fI[1]*vP1[0];
    vIcrossfP1[0] = vI[1]*fP1[2] - vI[2]*fP1[1];
    vIcrossfP1[1] = vI[2]*fP1[0] - vI[0]*fP1[2];
    vIcrossfP1[2] = vI[0]*fP1[1] - vI[1]*fP1[0];


    HI[0] = inv_uIdotuiP1plus1 * 
      (  (fIcrossfP1[0] +  vIcrossvP1[0]) * costwist0  + 
	 (vIcrossfP1[0] -  fIcrossvP1[0]) * sintwist0  -
	 cosalphaiplusgammaieq*uIcrossuP1[0]);
    HI[1] = inv_uIdotuiP1plus1 * 
      (  (fIcrossfP1[1] +  vIcrossvP1[1]) * costwist0  + 
	 (vIcrossfP1[1] -  fIcrossvP1[1]) * sintwist0  -
         cosalphaiplusgammaieq*uIcrossuP1[1]);
    HI[2] = inv_uIdotuiP1plus1 * 
      (  (fIcrossfP1[2] +  vIcrossvP1[2]) * costwist0  + 
	 (vIcrossfP1[2] -  fIcrossvP1[2]) * sintwist0  -
         cosalphaiplusgammaieq*uIcrossuP1[2]);

    kHI[0] = ktwist[type] * HI[0];
    kHI[1] = ktwist[type] * HI[1];
    kHI[2] = ktwist[type] * HI[2];

    /* ----------------------------------------------
       Allignment torque
       ---------------------------------------------- */

    tI[0] = x[iP1][0] - x[i][0]; // get vector between i and iP1
    tI[1] = x[iP1][1] - x[i][1];
    tI[2] = x[iP1][2] - x[i][2];
    dx=tI[0]; // store bead separations for virial calculation later
    dy=tI[1];
    dz=tI[2];
    bI2 = tI[0]*tI[0] + tI[1]*tI[1] + tI[2]*tI[2];
    inv_bI = 1.0/sqrt(bI2);
    tI[0]*=inv_bI;  // make tI a unit vector
    tI[1]*=inv_bI;
    tI[2]*=inv_bI;
    
    // nonzero alignment?
    //double cosphi,phi,cosphieq, align0rad;
    //twist0rad = align0[type] * PI / 180.0;

    //cosphi = uI[0]*tI[0] + uI[1]*tI[1] + uI[2]*tI[2];
    //phi = acos(cosphi);
    //cosphieq = cos(phi - twist0rad);

    

    uIcrosstI[0] = uI[1]*tI[2] - uI[2]*tI[1];
    uIcrosstI[1] = uI[2]*tI[0] - uI[0]*tI[2];
    uIcrosstI[2] = uI[0]*tI[1] - uI[1]*tI[0];

    /* ----------------------------------------------
       Total torque
       ---------------------------------------------- */
    
    torque[i][0] += kalign[type]*uIcrosstI[0] + kHI[0];  
    torque[i][1] += kalign[type]*uIcrosstI[1] + kHI[1]; 
    torque[i][2] += kalign[type]*uIcrosstI[2] + kHI[2];
 
    torque[iP1][0] -= kHI[0];  
    torque[iP1][1] -= kHI[1]; 
    torque[iP1][2] -= kHI[2]; 
    

    /* ----------------------------------------------
       Force
       ---------------------------------------------- */

    uIdottI = uI[0]*tI[0] + uI[1]*tI[1] + uI[2]*tI[2];
    
    if (uIdottI > 1.0) uIdottI = 1.0;
    if (uIdottI < -1.0) uIdottI = -1.0;

    GI[0] = inv_bI * ( uIdottI*tI[0] - uI[0] );
    GI[1] = inv_bI * ( uIdottI*tI[1] - uI[1] );
    GI[2] = inv_bI * ( uIdottI*tI[2] - uI[2] );
   
    f1[0] = kalign[type] * GI[0] ; 
    f1[1] = kalign[type] * GI[1] ;
    f1[2] = kalign[type] * GI[2] ;
      
    f[i][0] += f1[0];
    f[i][1] += f1[1];
    f[i][2] += f1[2];
    
    f[iP1][0] -= f1[0];
    f[iP1][1] -= f1[1];
    f[iP1][2] -= f1[2];
    
    f3[0] = f3[1] = f3[2]  = 0.0;  // for virial calculation

    if (eflag) eangle = kalign[type] * (1 - uIdottI) + ktwist[type] * (1 - cosalphaiplusgammaieq);

    if (evflag) // tally energy (virial=0 because force=0)
      ev_tally(i,iP1,iDummy,nlocal,newton_bond,eangle,f1,f3,
               dx,dy,dz,0.0,0.0,0.0);

  }
}

/* ---------------------------------------------------------------------- */

void AngleWLCTwist::allocate()
{
  allocated = 1;
  int n = atom->nangletypes;

  memory->create(kalign,n+1,"angle:kalign");
  memory->create(ktwist,n+1,"angle:ktwist");
  memory->create(twist0,n+1,"angle:eqtwistangle");

  memory->create(setflag,n+1,"angle:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

void AngleWLCTwist::coeff(int narg, char **arg)
{
  if (narg != 4) error->all(FLERR,"Incorrect args for angle coefficients");
  if (!allocated) allocate();

  int ilo,ihi;
  force->bounds(FLERR,arg[0],atom->nangletypes,ilo,ihi);

  double kalign_one = force->numeric(FLERR,arg[1]);
  double ktwist_one = force->numeric(FLERR,arg[2]);
  double twist0_one = force->numeric(FLERR,arg[3]);

  // convert gamma0 from degrees to radians

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    kalign[i] = kalign_one;
    ktwist[i] = ktwist_one;
    twist0[i] = twist0_one;
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR,"Incorrect args for angle coefficients");
}

/* ----------------------------------------------------------------------
   used by SHAKE
------------------------------------------------------------------------- */

double AngleWLCTwist::equilibrium_angle(int i)
{
  return 0.0;
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void AngleWLCTwist::write_restart(FILE *fp)
{
  fwrite(&kalign[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&ktwist[1],sizeof(double),atom->nangletypes,fp);
  fwrite(&twist0[1],sizeof(double),atom->nangletypes,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void AngleWLCTwist::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    fread(&kalign[1],sizeof(double),atom->nangletypes,fp);
    fread(&ktwist[1],sizeof(double),atom->nangletypes,fp);
    fread(&twist0[1],sizeof(double),atom->nangletypes,fp);
  }
  MPI_Bcast(&kalign[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&ktwist[1],atom->nangletypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&twist0[1],atom->nangletypes,MPI_DOUBLE,0,world);

  for (int i = 1; i <= atom->nangletypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   used by ComputeAngleLocal
------------------------------------------------------------------------- */

double AngleWLCTwist::single(int type, int i, int iP1, int iDummy)
{
  double **x = atom->x; // position vector
  int *ellipsoid = atom->ellipsoid;
  AtomVecEllipsoid *avec = (AtomVecEllipsoid *) atom->style_match("ellipsoid");
  AtomVecEllipsoid::Bonus *bonus = avec->bonus;

  double delx = x[iP1][0] - x[i][0];
  double dely = x[iP1][1] - x[i][1];
  double delz = x[iP1][2] - x[i][2];

  domain->minimum_image(delx,dely,delz);

  double PI, costwist0, sintwist0, twist0rad, alphaiplusgammai, cosalphaiplusgammaieq;
  double f0[3], v0[3], u0[3], *quatP1, *quatI;  
  double quatItmp[4], quatP1tmp[4];
  PI = 3.14159265359;
  twist0rad = twist0[type] * PI / 180.0;
  costwist0 = cos(twist0rad);
  sintwist0 = sin(twist0rad);

  double uI[3],uP1[3],fI[3],fP1[3],vI[3],vP1[3];

  f0[0] = 1; f0[1] = 0; f0[2] = 0;   
  v0[0] = 0; v0[1] = 1; v0[2] = 0;  
  u0[0] = 0; u0[1] = 0; u0[2] = 1; 

  quatI = bonus[ellipsoid[i]].quat;
  quatItmp[0] = quatI[0]; quatItmp[1] = quatI[1]; quatItmp[2] = quatI[2]; quatItmp[3] = quatI[3];
  //fvu for site 1
  LAMMPS_NS::quat_vec_rot(fI,f0,quatItmp);
  LAMMPS_NS::quat_vec_rot(vI,v0,quatItmp);
  LAMMPS_NS::quat_vec_rot(uI,u0,quatItmp);
  LAMMPS_NS::vec_norm(fI);
  LAMMPS_NS::vec_norm(vI);
  LAMMPS_NS::vec_norm(uI);

  quatP1 = bonus[ellipsoid[iP1]].quat;
  quatP1tmp[0] = quatP1[0]; quatP1tmp[1] = quatP1[1]; quatP1tmp[2] = quatP1[2]; quatP1tmp[3] = quatP1[3];
  //fvu for site 2
  LAMMPS_NS::quat_vec_rot(fP1,f0,quatP1tmp);
  LAMMPS_NS::quat_vec_rot(vP1,v0,quatP1tmp);
  LAMMPS_NS::quat_vec_rot(uP1,u0,quatP1tmp);
  LAMMPS_NS::vec_norm(fP1);
  LAMMPS_NS::vec_norm(vP1);
  LAMMPS_NS::vec_norm(uP1);

  double inv_uIdotuiP1plus1 = 1.0/(uI[0]*uP1[0] + uI[1]*uP1[1] + uI[2]*uP1[2] + 1);
  double cosalphaiplusgammai = (fP1[0]*fI[0] + fP1[1]*fI[1] + fP1[2]*fI[2] + 
				vP1[0]*vI[0] + vP1[1]*vI[1] + vP1[2]*vI[2]) * inv_uIdotuiP1plus1;
			   
  if (cosalphaiplusgammai > 1.0) cosalphaiplusgammai = 1.0;
  if (cosalphaiplusgammai < -1.0) cosalphaiplusgammai = -1.0;
  
  alphaiplusgammai = acos(cosalphaiplusgammai);
  cosalphaiplusgammaieq = cos(alphaiplusgammai-twist0rad);
  
  double bI2 = delx*delx + dely*dely + delz*delz;
  double cosphii = uI[0]*delx + uI[1]*dely + uI[2]*delz;
  cosphii/=sqrt(bI2);
      
  if (cosphii > 1.0) cosphii = 1.0;
  if (cosphii < -1.0) cosphii = -1.0;

  return kalign[type]*(1-cosphii) +  ktwist[type] * (1 - cosalphaiplusgammaieq); // energy
}
