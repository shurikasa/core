/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef DWR_H
#define DWR_H

/** \page dwr Dual Weighted Residuals
  * This is the main API for SCOREC's DWR library
  *
  * These functions provide capabilities to solve adjoint
  * boundary value problems using a specified functional
  * quantity of interest.
  */

/** \file dwr.h */

namespace apf {
class Mesh;
class Field;
}

/** \namespace dwr
  * \brief All DWR symbols */
namespace dwr {

/** \brief create an apf::Field using hierarchic shape functions
  * \param m the mesh over which the field is defined
  * \param name a unique name for this field
  * \param valueType the type of field data
  * \param order the polynomial order of the shape functions
  */
apf::Field* createHierarchicField(apf::Mesh* m, const char* name,
    int valueType, int order);

}

#endif