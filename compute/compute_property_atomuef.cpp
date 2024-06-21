/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
   
   ORNL: Zhiqiang Shen 
------------------------------------------------------------------------- */

#include "compute_property_atomuef.h"
#include <cmath>
#include <cstring>
#include "math_extra.h"
#include "atom.h"
#include "fix.h"
#include "atom_vec.h"
#include "fix_nh_uef.h"
#include "modify.h"
#include "update.h"
#include "domain.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputePropertyAtomuef::ComputePropertyAtomuef(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  index(NULL), pack_choice(NULL)
{
  if (narg < 4) error->all(FLERR,"Illegal compute property/atomuef command");

  peratom_flag = 1;
  nvalues = narg - 3;
  if (nvalues == 1) size_peratom_cols = 0;
  else size_peratom_cols = nvalues;

  // parse input values
  // customize a new keyword by adding to if statement

  pack_choice = new FnPtrPack[nvalues];
  index = new int[nvalues];

  int i;
  for (int iarg = 3; iarg < narg; iarg++) {
    i = iarg-3;

    if (strcmp(arg[iarg],"xu") == 0) {
        pack_choice[i] = &ComputePropertyAtomuef::pack_xu_triclinic;
    } else if (strcmp(arg[iarg],"yu") == 0) {
        pack_choice[i] = &ComputePropertyAtomuef::pack_yu_triclinic;
    } else if (strcmp(arg[iarg],"zu") == 0) {
        pack_choice[i] = &ComputePropertyAtomuef::pack_zu_triclinic;
    // check if atom style recognizes keyword

    } else {
      index[i] = atom->avec->property_atom(arg[iarg]);
      if (index[i] < 0)
        error->all(FLERR,"Invalid keyword in compute property/atomuef command");
      pack_choice[i] = &ComputePropertyAtomuef::pack_property_atom;
    }
  }

  nmax = 0;
}

/* ---------------------------------------------------------------------- */

ComputePropertyAtomuef::~ComputePropertyAtomuef()
{
  delete [] pack_choice;
  delete [] index;
  memory->destroy(vector_atom);
  memory->destroy(array_atom);
}

/* ---------------------------------------------------------------------- */

void ComputePropertyAtomuef::init()
{
  
  //check for the uef fix 
  int i=0;
  for (i=0; i<modify->nfix; i++)
    {
      if (strcmp(modify->fix[i]->style,"nvt/uef")==0)
	break;
    }
  if (i==modify->nfix)
    error->all(FLERR,"Can't use compute property/atomuef without defining a fix nvt/uef");
  
  ifix_uef=i;
  
}

/* ---------------------------------------------------------------------- */

void ComputePropertyAtomuef::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow vector or array if necessary

  if (atom->nmax > nmax) {
    nmax = atom->nmax;
    if (nvalues == 1) {
      memory->destroy(vector_atom);
      memory->create(vector_atom,nmax,"property/atom:vector");
    } else {
      memory->destroy(array_atom);
      memory->create(array_atom,nmax,nvalues,"property/atom:array");
    }
  }

  // fill vector or array with per-atom values

  if (nvalues == 1) {
    buf = vector_atom;
    (this->*pack_choice[0])(0);
  } else {
    if (nmax) buf = &array_atom[0][0];
    else buf = NULL;
    for (int n = 0; n < nvalues; n++)
      (this->*pack_choice[n])(n);
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputePropertyAtomuef::memory_usage()
{
  double bytes = nmax*nvalues * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   one method for every keyword compute property/atom can output
   the atom property is packed into buf starting at n with stride nvalues
   customize a new keyword by adding a method
------------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */

void ComputePropertyAtomuef::pack_xu_triclinic(int n)
{
  double **x = atom->x;
  imageint *image = atom->image;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double *h = domain->h;
  int xbox,ybox,zbox;

  double rot[3][3];
  ((FixNHUef*) modify->fix[ifix_uef])->get_rot(rot);
  double p[3];
  double xn[3];
  
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      xbox = (image[i] & IMGMASK) - IMGMAX;
      ybox = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
      zbox = (image[i] >> IMG2BITS) - IMGMAX;
      p[0] = x[i][0] + h[0]*xbox + h[5]*ybox + h[4]*zbox;
      p[1] = x[i][1] + h[1]*ybox + h[3]*zbox;
      p[2] = x[i][2] + h[2]*zbox;
      p[0] -= domain->boxlo[0];
      p[1] -= domain->boxlo[1];
      p[2] -= domain->boxlo[2];
      xn[0]=rot[0][0]*p[0]+rot[1][0]*p[1]+rot[2][0]*p[2];
      xn[1]=rot[0][1]*p[0]+rot[1][1]*p[1]+rot[2][1]*p[2];
      xn[2]=rot[0][2]*p[0]+rot[1][2]*p[1]+rot[2][2]*p[2];
       
      buf[n] = xn[0];
    } else buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputePropertyAtomuef::pack_yu_triclinic(int n)
{
  double **x = atom->x;
  imageint *image = atom->image;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double *h = domain->h;
  int xbox,ybox,zbox;
  
  double rot[3][3];
  ((FixNHUef*) modify->fix[ifix_uef])->get_rot(rot);
  double p[3];
  double xn[3];
  

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      xbox = (image[i] & IMGMASK) - IMGMAX;
      ybox = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
      zbox = (image[i] >> IMG2BITS) - IMGMAX;
      p[0] = x[i][0] + h[0]*xbox + h[5]*ybox + h[4]*zbox;
      p[1] = x[i][1] + h[1]*ybox + h[3]*zbox;
      p[2] = x[i][2] + h[2]*zbox;
      p[0] -= domain->boxlo[0];
      p[1] -= domain->boxlo[1];
      p[2] -= domain->boxlo[2];
      xn[0]=rot[0][0]*p[0]+rot[1][0]*p[1]+rot[2][0]*p[2];
      xn[1]=rot[0][1]*p[0]+rot[1][1]*p[1]+rot[2][1]*p[2];
      xn[2]=rot[0][2]*p[0]+rot[1][2]*p[1]+rot[2][2]*p[2];
       
      buf[n] = xn[1];
    } else buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputePropertyAtomuef::pack_zu_triclinic(int n)
{
  double **x = atom->x;
  imageint *image = atom->image;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double *h = domain->h;
  int xbox,ybox,zbox;
  
  double rot[3][3];
  ((FixNHUef*) modify->fix[ifix_uef])->get_rot(rot);
  double p[3];
  double xn[3];

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      xbox = (image[i] & IMGMASK) - IMGMAX;
      ybox = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
      zbox = (image[i] >> IMG2BITS) - IMGMAX;
      p[0] = x[i][0] + h[0]*xbox + h[5]*ybox + h[4]*zbox;
      p[1] = x[i][1] + h[1]*ybox + h[3]*zbox;
      p[2] = x[i][2] + h[2]*zbox;
      p[0] -= domain->boxlo[0];
      p[1] -= domain->boxlo[1];
      p[2] -= domain->boxlo[2];
      xn[0]=rot[0][0]*p[0]+rot[1][0]*p[1]+rot[2][0]*p[2];
      xn[1]=rot[0][1]*p[0]+rot[1][1]*p[1]+rot[2][1]*p[2];
      xn[2]=rot[0][2]*p[0]+rot[1][2]*p[1]+rot[2][2]*p[2];
      buf[n] = xn[2];
    } else buf[n] = 0.0;
    n += nvalues;
  }
}

/* ---------------------------------------------------------------------- */

void ComputePropertyAtomuef::pack_property_atom(int n)
{
  atom->avec->pack_property_atom(index[n],&buf[n],nvalues,groupbit);
}
