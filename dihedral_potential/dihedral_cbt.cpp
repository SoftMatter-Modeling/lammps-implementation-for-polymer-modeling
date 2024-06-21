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
   Contributing author: Zhiqiang Shen (CNMS, ORNL)
   CBT potential: Ref: [1] J. Chem. Phys. 123, 114901 (2005); [2]  J. Chem. Theory Comput. 2013, 9, 3282−3292 
------------------------------------------------------------------------- */

#include "dihedral_cbt.h"
#include <mpi.h>
#include <cmath>
#include "atom.h"
#include "neighbor.h"
#include "comm.h"
#include "domain.h"
#include "force.h"
#include "update.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "utils.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define TOLERANCE 0.05
#define SMALL     0.00001

/* ---------------------------------------------------------------------- */

DihedralCbt::DihedralCbt(LAMMPS *lmp) : Dihedral(lmp)
{
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

DihedralCbt::~DihedralCbt()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(k);
    memory->destroy(a0);
    memory->destroy(a1);
    memory->destroy(a2);
    memory->destroy(a3);
  }
}

/* ---------------------------------------------------------------------- */

void DihedralCbt::compute(int eflag, int vflag)
{
  int i1,i2,i3,i4,n,type;
  double vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z;
  double ax,ay,az,bx,by,bz,cx,cy,cz;
  double ra_sq, rb_sq, rc_sq, ra, rb, rc,a_dot_b, b_dot_c, a_dot_c;
  double c1, c2, c3, cphi_part, c, s1, s1_sq, s2, s2_sq; 
  double an_cn, p; 
  double fac_theta1, fac_theta2, F1_theta1_x, F1_theta1_y, F1_theta1_z, F3_theta1_x, F3_theta1_y,F3_theta1_z; 
  double F2_theta2_x, F2_theta2_y, F2_theta2_z, F4_theta2_x, F4_theta2_y, F4_theta2_z; 
  double fac_phi, a11, a12, a13, a22, a23, a33; 
  double F1_phi_x, F1_phi_y, F1_phi_z, F4_phi_x, F4_phi_y, F4_phi_z, F23_phi_x, F23_phi_y, F23_phi_z; 
  
  double edihedral,f1[3],f2[3],f3[3],f4[3];

  edihedral = 0.0;
  ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;
  int **dihedrallist = neighbor->dihedrallist;
  int ndihedrallist = neighbor->ndihedrallist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  for (n = 0; n < ndihedrallist; n++) {
    i1 = dihedrallist[n][0];
    i2 = dihedrallist[n][1];
    i3 = dihedrallist[n][2];
    i4 = dihedrallist[n][3];
    type = dihedrallist[n][4];

    // 1st bond
    ax=x[i2][0] - x[i1][0];
    ay=x[i2][1] - x[i1][1];
    az=x[i2][2] - x[i1][2];
    
    // 2nd bond    
    bx=x[i3][0] - x[i2][0];
    by=x[i3][1] - x[i2][1];
    bz=x[i3][2] - x[i2][2];

    // 3rd bond
    cx=x[i4][0] - x[i3][0];
    cy=x[i4][1] - x[i3][1];
    cz=x[i4][2] - x[i3][2];
    
  //Consider periodic boundary conditions:
    domain->minimum_image(ax,ay,az);
    domain->minimum_image(bx,by,bz);
    domain->minimum_image(cx,cy,cz);
  
    // calculate the bond length 
    ra_sq = ax*ax+ay*ay+az*az;
    ra    = sqrt(ra_sq);
    
    rb_sq = bx*bx+by*by+bz*bz;
    rb    = sqrt(rb_sq);
    
    rc_sq = cx*cx+cy*cy+cz*cz;
    rc    = sqrt(rc_sq);
    
    //calculate the cos sin value 
    a_dot_b = ax*bx+ay*by+az*bz;
    a_dot_c = ax*cx+ay*cy+az*cz;    
    b_dot_c = bx*cx+by*cy+bz*cz;
    
    c1 = -a_dot_b/(ra*rb);
    c2 = -b_dot_c/(rb*rc);
    c3 = -a_dot_c/(ra*rc); 
    
    
   // if (c1 > 1.0)  c1 = 1.0; // numerical correction 
   // if (c1 < -1.0) c1 = -1.0;
    
   // if (c2 > 1.0)  c2 = 1.0;
   // if (c2 < -1.0) c2 = -1.0;
    
   //  if (c3 > 1.0)  c3 = 1.0;
   // if (c3 < -1.0) c3 = -1.0;
    
    
    cphi_part =c1*c2+c3; 

    s1_sq = MAX(1.0 - c1*c1,0.0);
    s1    = sqrt(s1_sq);
    if(s1<SMALL) s1=SMALL; // numerical correction 
    s1_sq = s1*s1;
     
    
    s2_sq = MAX(1.0 - c2*c2,0.0);
    s2    = sqrt(s2_sq);
    if(s2<SMALL) s2=SMALL; // numerical correction 
    s2_sq = s2*s2; 
    
    // cos phi 
    c=cphi_part/(s1*s2);
    
    // error check
     if (c > 1.0 + TOLERANCE || c < (-1.0 - TOLERANCE)) {
      int me;
      MPI_Comm_rank(world,&me);
      if (screen) {
        char str[128];
        sprintf(str,"Dihedral problem: %d " BIGINT_FORMAT " "
                TAGINT_FORMAT " " TAGINT_FORMAT " "
                TAGINT_FORMAT " " TAGINT_FORMAT,
                me,update->ntimestep,
                atom->tag[i1],atom->tag[i2],atom->tag[i3],atom->tag[i4]);
        error->warning(FLERR,str,0);
        fprintf(screen,"  1st atom: %d %g %g %g\n",
                me,x[i1][0],x[i1][1],x[i1][2]);
        fprintf(screen,"  2nd atom: %d %g %g %g\n",
                me,x[i2][0],x[i2][1],x[i2][2]);
        fprintf(screen,"  3rd atom: %d %g %g %g\n",
                me,x[i3][0],x[i3][1],x[i3][2]);
        fprintf(screen,"  4th atom: %d %g %g %g\n",
                me,x[i4][0],x[i4][1],x[i4][2]);
      }
    }
    
    if (c > 1.0) c = 1.0;
    if (c < -1.0) c = -1.0;
    
    an_cn = a0[type]+a1[type]*c+a2[type]*c*c+a3[type]*c*c*c; 
    
    p     = k[type]*pow(s1,3.0)*pow(s2,3.0)*an_cn;  
    if (eflag) edihedral = p;
    
    
    //======angle force part 
    
    fac_theta1 = 3.0*k[type]*s1*c1*s2*s2*s2*an_cn; 
    fac_theta2 = 3.0*k[type]*s1*s1*s1*s2*c2*an_cn;
    
    F1_theta1_x = bx/(ra*rb)+ax*(c1/ra_sq); 
    F1_theta1_y = by/(ra*rb)+ay*(c1/ra_sq); 
    F1_theta1_z = bz/(ra*rb)+az*(c1/ra_sq);
       
    
    F3_theta1_x = ax/(ra*rb) + bx*(c1/rb_sq);
    F3_theta1_y = ay/(ra*rb) + by*(c1/rb_sq);
    F3_theta1_z = az/(ra*rb) + bz*(c1/rb_sq);
    
    F1_theta1_x *= fac_theta1;
    F1_theta1_y *= fac_theta1;
    F1_theta1_z *= fac_theta1;
    
    F3_theta1_x *= -fac_theta1;
    F3_theta1_y *= -fac_theta1;
    F3_theta1_z *= -fac_theta1;
    
    
    F2_theta2_x = cx/(rb*rc)+bx*(c2/rb_sq); 
    F2_theta2_y = cy/(rb*rc)+by*(c2/rb_sq); 
    F2_theta2_z = cz/(rb*rc)+bz*(c2/rb_sq);
    
    F4_theta2_x = bx/(rb*rc)+cx*(c2/rc_sq); 
    F4_theta2_y = by/(rb*rc)+cy*(c2/rc_sq); 
    F4_theta2_z = bz/(rb*rc)+cz*(c2/rc_sq);

    F2_theta2_x *= fac_theta2;
    F2_theta2_y *= fac_theta2;
    F2_theta2_z *= fac_theta2;
 
    F4_theta2_x *= -fac_theta2;
    F4_theta2_y *= -fac_theta2;
    F4_theta2_z *= -fac_theta2; 
    
    //======Dihedral force part 
    
    fac_phi = -k[type]*(a1[type] + 2.0*a2[type]*c + 3.0*a3[type]*c*c); 

   /* detailed part based on the derivation, it is easier to check 
    F1_phi_x= ax*(s2_sq*cphi_part/ra_sq) + bx*(s1_sq*s2_sq*c2+s2_sq*c1*cphi_part)/ra/rb + cx*(s1_sq*s2_sq/ra/rc); 
    
    F4_phi_x= -ax*(s1_sq*s2_sq/ra/rc) -bx*(s1_sq*s2_sq*c1+s1_sq*c2*cphi_part)/rb/rc -cx*(s1_sq*cphi_part/rc_sq); 
  
    F23_phi_x =ax*(s1_sq*s2_sq*c2+s2_sq*c1*cphi_part)/ra/rb + bx*(-2.0*c3*s1_sq*s2_sq+s2_sq*cphi_part+s1_sq*cphi_part)/rb_sq+cx*(s1_sq*s2_sq*c1+s1_sq*c2*cphi_part)/rb/rc; 
   */
   
    a11 = s2_sq*cphi_part/ra_sq;
    a12 = (s1_sq*s2_sq*c2+s2_sq*c1*cphi_part)/ra/rb;
    a13 = s1_sq*s2_sq/ra/rc;   

    a22 = (-2.0*c3*s1_sq*s2_sq+s2_sq*cphi_part+s1_sq*cphi_part)/rb_sq;    
    a23 = (s1_sq*s2_sq*c1+s1_sq*c2*cphi_part)/rb/rc; 
   
    a33 = s1_sq*cphi_part/rc_sq; 
   
    F1_phi_x = a11*ax + a12*bx + a13*cx;
    F1_phi_y = a11*ay + a12*by + a13*cy;
    F1_phi_z = a11*az + a12*bz + a13*cz;

    F1_phi_x *= fac_phi; 
    F1_phi_y *= fac_phi;
    F1_phi_z *= fac_phi;
   
    F4_phi_x = -a13*ax - a23*bx - a33*cx;
    F4_phi_y = -a13*ay - a23*by - a33*cy;
    F4_phi_z = -a13*az - a23*bz - a33*cz; 
   
    F4_phi_x *= fac_phi; 
    F4_phi_y *= fac_phi;
    F4_phi_z *= fac_phi;
   
   
    F23_phi_x   = a12*ax + a22*bx + a23*cx;
    F23_phi_y   = a12*ay + a22*by + a23*cy;
    F23_phi_z   = a12*az + a22*bz + a23*cz;

    F23_phi_x *= fac_phi; 
    F23_phi_y *= fac_phi;
    F23_phi_z *= fac_phi;
   
    f1[0] = F1_theta1_x + F1_phi_x;
    f1[1] = F1_theta1_y + F1_phi_y;
    f1[2] = F1_theta1_z + F1_phi_z;

    f2[0] = -F1_theta1_x-F3_theta1_x+F2_theta2_x-F1_phi_x+F23_phi_x;
    f2[1] = -F1_theta1_y-F3_theta1_y+F2_theta2_y-F1_phi_y+F23_phi_y;
    f2[2] = -F1_theta1_z-F3_theta1_z+F2_theta2_z-F1_phi_z+F23_phi_z;

    f3[0] = F3_theta1_x-F2_theta2_x-F4_theta2_x-F4_phi_x-F23_phi_x;
    f3[1] = F3_theta1_y-F2_theta2_y-F4_theta2_y-F4_phi_y-F23_phi_y;;
    f3[2] = F3_theta1_z-F2_theta2_z-F4_theta2_z-F4_phi_z-F23_phi_z;;

    f4[0] = F4_theta2_x + F4_phi_x;
    f4[1] = F4_theta2_y + F4_phi_y;
    f4[2] = F4_theta2_z + F4_phi_z;
    
    // apply force to each of 4 atoms

    if (newton_bond || i1 < nlocal) {
      f[i1][0] += f1[0];
      f[i1][1] += f1[1];
      f[i1][2] += f1[2];
    }

    if (newton_bond || i2 < nlocal) {
      f[i2][0] += f2[0];
      f[i2][1] += f2[1];
      f[i2][2] += f2[2];
    }

    if (newton_bond || i3 < nlocal) {
      f[i3][0] += f3[0];
      f[i3][1] += f3[1];
      f[i3][2] += f3[2];
    }

    if (newton_bond || i4 < nlocal) {
      f[i4][0] += f4[0];
      f[i4][1] += f4[1];
      f[i4][2] += f4[2];
    }
    
 // for the calculation of virial pressure
 
    vb1x = x[i1][0] - x[i2][0];
    vb1y = x[i1][1] - x[i2][1];
    vb1z = x[i1][2] - x[i2][2];
        
    vb2x = x[i3][0] - x[i2][0];
    vb2y = x[i3][1] - x[i2][1];
    vb2z = x[i3][2] - x[i2][2];
     
    vb3x = x[i4][0] - x[i3][0];
    vb3y = x[i4][1] - x[i3][1];
    vb3z = x[i4][2] - x[i3][2];
    
    domain->minimum_image(vb1x,vb1y,vb1z);
    domain->minimum_image(vb2x,vb2y,vb2z);
    domain->minimum_image(vb3x,vb3y,vb3z);
    
    if (evflag)
      ev_tally(i1,i2,i3,i4,nlocal,newton_bond,edihedral,f1,f3,f4,
               vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z);
  }
}

/* ---------------------------------------------------------------------- */

void DihedralCbt::allocate()
{
  allocated = 1;
  int n = atom->ndihedraltypes;

  memory->create(k,n+1,"dihedral:k");
  memory->create(a0,n+1,"dihedral:a0");
  memory->create(a1,n+1,"dihedral:a1");
  memory->create(a2,n+1,"dihedral:a2");
  memory->create(a3,n+1,"dihedral:a3");

  memory->create(setflag,n+1,"dihedral:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one type
------------------------------------------------------------------------- */

void DihedralCbt::coeff(int narg, char **arg)
{
  if (narg != 6) error->all(FLERR,"Incorrect args for cbt coefficients");
  if (!allocated) allocate();

  int ilo,ihi;
  force->bounds(FLERR,arg[0],atom->ndihedraltypes,ilo,ihi);

  double k_one  = force->numeric(FLERR,arg[1]);
  double a0_one = force->numeric(FLERR,arg[2]);
  double a1_one = force->numeric(FLERR,arg[3]);
  double a2_one = force->numeric(FLERR,arg[4]);
  double a3_one = force->numeric(FLERR,arg[5]);
  
  // require k >= 0
  if (k_one < 0.0)
    error->all(FLERR,"Incorrect coefficient arg for dihedral coefficients");

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    k[i]  = k_one;
    a0[i] = a0_one;
    a1[i] = a1_one;
    a2[i] = a2_one;
    a3[i] = a3_one;
    
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR,"Incorrect args for dihedral coefficients");
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void DihedralCbt::write_restart(FILE *fp)
{
  fwrite(&k[1],sizeof(double),atom->ndihedraltypes,fp);
  fwrite(&a0[1],sizeof(double),atom->ndihedraltypes,fp);
  fwrite(&a1[1],sizeof(double),atom->ndihedraltypes,fp);
  fwrite(&a2[1],sizeof(double),atom->ndihedraltypes,fp);
  fwrite(&a3[1],sizeof(double),atom->ndihedraltypes,fp);

}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void DihedralCbt::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    utils::sfread(FLERR,&k[1],sizeof(double),atom->ndihedraltypes,fp,NULL,error);
    utils::sfread(FLERR,&a0[1],sizeof(double),atom->ndihedraltypes,fp,NULL,error);
    utils::sfread(FLERR,&a1[1],sizeof(double),atom->ndihedraltypes,fp,NULL,error);
    utils::sfread(FLERR,&a2[1],sizeof(double),atom->ndihedraltypes,fp,NULL,error);
    utils::sfread(FLERR,&a3[1],sizeof(double),atom->ndihedraltypes,fp,NULL,error);
    
  }
  MPI_Bcast(&k[1],atom->ndihedraltypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&a0[1],atom->ndihedraltypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&a1[1],atom->ndihedraltypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&a2[1],atom->ndihedraltypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&a3[1],atom->ndihedraltypes,MPI_DOUBLE,0,world);


  for (int i = 1; i <= atom->ndihedraltypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void DihedralCbt::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ndihedraltypes; i++)
    fprintf(fp,"%d %g %g %g %g %g\n",i,k[i],a0[i],a1[i],a2[i],a3[i]);
}
