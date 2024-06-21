/*
  ----------------------------------------------------------------------
  LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
  http://lammps.sandia.gov, Sandia National Laboratories
  Steve Plimpton, sjplimp@sandia.gov

  Copyright (2003) Sandia Corporation.  Under the terms of Contract
  DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
  certain rights in this software.  This software is distributed under
  the GNU General Public License.

  See the README file in the top-level LAMMPS directory.

  Contributing author: Zhiqiang Shen 
  output the rotation matrix during the simulation of uniaxial extension 
  -------------------------------------------------------------------------
*/

#include "mpi.h"
#include <cstring>
#include <cstdlib>
#include "compute_rotation_uef.h"
#include "fix_nh_uef.h"
#include "update.h"
#include "domain.h"
#include "modify.h"
#include "fix.h"
#include "error.h"

using namespace LAMMPS_NS;
/*
  ----------------------------------------------------------------------
  Default values for the ext flags
  ----------------------------------------------------------------------
*/
ComputeRotationUef::ComputeRotationUef(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  vector_flag = 1;
  extvector   = 1;
  size_vector = 12;
  timeflag    = 1;
  peflag      = 1;
  vector      = new double[12];
}

ComputeRotationUef::~ComputeRotationUef()
{
  delete [] vector;
}

/*
  ----------------------------------------------------------------------
  Check for the uef fix
  ----------------------------------------------------------------------
*/
void ComputeRotationUef::init()
{
  int i=0;
  for (i=0; i<modify->nfix; i++)
    {
      if (strcmp(modify->fix[i]->style,"nvt/uef")==0)
	break;
    }
  if (i==modify->nfix)
    error->all(FLERR,"Can't use compute rotation/uef without defining a fix nvt/uef");
  
  ifix_uef=i;
}

/* ----------------------------------------------------------------------
   Get the rotation tensor
   ------------------------------------------------------------------------- */

void ComputeRotationUef::compute_vector()
{

  ((FixNHUef*) modify->fix[ifix_uef])->get_rot(rot);
  vector[0]=rot[0][0];
  vector[1]=rot[0][1];
  vector[2]=rot[0][2];
  vector[3]=rot[1][0];
  vector[4]=rot[1][1];
  vector[5]=rot[1][2];
  vector[6]=rot[2][0];
  vector[7]=rot[2][1];
  vector[8]=rot[2][2];
  vector[9] =domain->boxlo[0];
  vector[10]=domain->boxlo[1];
  vector[11]=domain->boxlo[2];

}
