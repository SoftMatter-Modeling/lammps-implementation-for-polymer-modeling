/* -*- c++ -*- ----------------------------------------------------------
  LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
  http://lammps.sandia.gov, Sandia National Laboratories
  Steve Plimpton, sjplimp@sandia.gov

  Copyright (2003) Sandia Corporation.  Under the terms of Contract
  DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
  certain rights in this software.  This software is distributed under
  the GNU General Public License.

  See the README file in the top-level LAMMPS directory.

  Contributing author: Zhiqiang Shen 
  -------------------------------------------------------------------------
*/

#ifdef COMPUTE_CLASS

ComputeStyle(rotation/uef,ComputeRotationUef)

#else

#ifndef LMP_COMPUTE_ROTATION_UEF_H
#define LMP_COMPUTE_ROTATION_UEF_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeRotationUef : public Compute {
 public:
  ComputeRotationUef(class LAMMPS *, int, char **);
  virtual ~ComputeRotationUef();
  virtual void init();
  virtual void compute_vector();

 protected:
  int ifix_uef;
  double rot[3][3];
};

}

#endif
#endif

/*
   ERROR/WARNING messages:

   This class inherits most of the warnings from ComputePressure. The
   only addition is:

   E: Can't use compute rotation/uef without defining a fix nvt/uef

   Self-explanatory.  

*/
