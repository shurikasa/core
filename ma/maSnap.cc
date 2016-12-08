/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#include <PCU.h>
#include "maSnap.h"
#include "maAdapt.h"
#include "maOperator.h"
#include "maSnapper.h"
#include "maLayer.h"
#include "maMatch.h"
#include <apfGeometry.h>
#include <cassert>
#include <iostream>

namespace ma {

/* this is the logic to deal with discontinuous
   periodic parametric coordinates.
   We assume that if the difference between
   the vertex coordinates of a mesh edge is more than half
   the periodic range, then the edge
   crosses over the discontinuity and we need
   to interpolate differently.
   The last parameter "mode" is there since sometimes
   this function needs to be called with mode = 0 and
   other times it needs to be called with mode = 1.
   For example, "interpolateParametricCoordinates".
   Also, this function needs to be called with with
   mode = 0 first. And if the resulting point corresponding
   to the parametric coordinates is not inside the model,
   it should be called again with mode =1 a second time*/
static double interpolateParametricCoordinate(
    double t,
    double a,
    double b,
    double range[2],
    bool isPeriodic,
    int mode)
{
  if ( ! isPeriodic)
    return (1-t)*a + t*b;
  if (range[0] > range[1])
    std::swap(range[0],range[1]);
  if (a > b)
  {
    std::swap(a,b);
    t = 1-t;
  }
  double period = range[1]-range[0];
  double span = b-a;
  if (!mode) {
    if (span < (period/2))
      return (1-t)*a + t*b;
  }
  else {
    if (span <= (period/2))
      return (1-t)*a + t*b;
  }
  a += period;
  double result = (1-t)*a + t*b;
  if (result >= range[1])
    result -= period;
  assert(result >= range[0]);
  assert(result < range[1]);
  return result;
}

void interpolateParametricCoordinates(
    apf::Mesh* m,
    Model* g,
    double t,
    Vector const& a,
    Vector const& b,
    Vector& p)
{
  double range[2];
  int dim = m->getModelType(g);
  for (int d=0; d < dim; ++d) {
    bool isPeriodic = m->getPeriodicRange(g,d,range);
    p[d] = interpolateParametricCoordinate(t,a[d],b[d],range,isPeriodic, 0);
  }

  /* check if the new point is inside the model.
   * otherwise re-run the above loop with last option
   * in "interpolateParametricCoordinae" being 1.
   * Notes
   * 1) we are assuming manifold surfaces
   * 2) we only check for faces that are periodic
   */

  // this need to be done for faces, only
  if (dim != 2)
    return;

  Vector x;
  bool ok;
  ok = m->isParamPointInsideModel(g, &p[0], x);
  if (ok)
    return;

  for (int d=0; d < dim; ++d) {
    bool isPeriodic = m->getPeriodicRange(g,d,range);
    p[d] = interpolateParametricCoordinate(t,a[d],b[d],range,isPeriodic, 1);
  }
}

static void transferParametricBetween(
    Mesh* m,
    Model* g,
    Entity* v[2],
    double t,
    Vector& p)
{
  Vector ep[2];
  for (int i=0; i < 2; ++i)
    m->getParamOn(g,v[i],ep[i]);
  ma::interpolateParametricCoordinates(m,g,t,ep[0],ep[1],p);
}

void transferParametricOnEdgeSplit(
    Mesh* m,
    Entity* e,
    double t,
    Vector& p)
{
  Model* g = m->toModel(e);
  int modelDimension = m->getModelType(g);
  if (modelDimension==m->getDimension()) return;
  Entity* ev[2];
  m->getDownward(e,0,ev);
  ma::transferParametricBetween(m, g, ev, t, p);
}

void transferParametricOnQuadSplit(
    Mesh* m,
    Entity* quad,
    Entity* v01,
    Entity* v32,
    double y,
    Vector& p)
{
  Model* g = m->toModel(quad);
  int modelDimension = m->getModelType(g);
  if (modelDimension==m->getDimension()) return;
  Entity* v[2];
  v[0] = v01; v[1] = v32;
  ma::transferParametricBetween(m, g, v, y, p);
}

static void getSnapPoint(Mesh* m, Entity* v, Vector& x)
{
  m->getPoint(v,0,x);
  Vector p;
  m->getParam(v,p);
  Model* g = m->toModel(v);
  m->snapToModel(g,p,x);
}

class SnapAll : public Operator
{
  public:
    SnapAll(Adapt* a, Tag* t, bool simple):
      snapper(a, t, simple)
    {
      adapter = a;
      tag = t;
      successCount = 0;
      didAnything = false;
      vert = 0;
    }
    int getTargetDimension() {return 0;}
    bool shouldApply(Entity* e)
    {
      if ( ! getFlag(adapter, e, SNAP))
        return false;
      vert = e;
      return true;
    }
    bool requestLocality(apf::CavityOp* o)
    {
      return snapper.setVert(vert, o);
    }
    void apply()
    {
      bool snapped = snapper.run();
      didAnything = didAnything || snapped || snapper.dug;
      if (snapped)
        ++successCount;
      clearFlag(adapter, vert, SNAP);
    }
    int successCount;
    bool didAnything;
  private:
    Adapt* adapter;
    Tag* tag;
    Entity* vert;
    Snapper snapper;
};

long tagVertsToSnap(Adapt* a, Tag*& t)
{
  Mesh* m = a->mesh;
  int dim = m->getDimension();
  t = m->createDoubleTag("ma_snap", 3);
  Entity* v;
  long n = 0;
  Iterator* it = m->begin(0);
  while ((v = m->iterate(it))) {
    int md = m->getModelType(m->toModel(v));
    if (md == dim)
      continue;
    Vector s;
    getSnapPoint(m, v, s);
    Vector x = getPosition(m, v);
    if (apf::areClose(s, x, 0.0))
      continue;
    m->setDoubleTag(v, t, &s[0]);
    if (m->isOwned(v))
      ++n;
  }
  m->end(it);
  return PCU_Add_Long(n);
}

static void markVertsToSnap(Adapt* a, Tag* t)
{
  HasTag p(a->mesh, t);
  markEntities(a, 0, p, SNAP, DONT_SNAP);
}

bool snapOneRound(Adapt* a, Tag* t, bool isSimple, long& successCount)
{
  markVertsToSnap(a, t);
  SnapAll op(a, t, isSimple);
  applyOperator(a, &op);
  successCount += PCU_Add_Long(op.successCount);
  return PCU_Or(op.didAnything);
}

long snapTaggedVerts(Adapt* a, Tag* tag)
{
  long successCount = 0;
  /* first snap all the vertices we can without digging.
     This is fast because it uses just the elements around
     the vertex and doesn't think much, it should also handle
     the vast majority of vertices */
  while (snapOneRound(a, tag, true, successCount));
  /* all the remaining vertices now need some kind of modification
     in order to snap.
     Here we turn on the "try digging before snapping" flag,
     which requires two-layer cavities so hopefully fewer vertices
     are involved here */
  while (snapOneRound(a, tag, false, successCount));
  return successCount;
}

void snap(Adapt* a)
{
  if ( ! a->input->shouldSnap)
    return;
  double t0 = PCU_Time();
  Tag* tag;
  /* we are starting to support a few operations on matched
     meshes, including snapping+UR. this should prevent snapping
     from modifying any matched entities */
  preventMatchedCavityMods(a);
  long targets = tagVertsToSnap(a, tag);
  long success = snapTaggedVerts(a, tag);
  snapLayer(a, tag);
  apf::removeTagFromDimension(a->mesh, tag, 0);
  a->mesh->destroyTag(tag);
  double t1 = PCU_Time();
  print("snapped in %f seconds: %ld targets, %ld non-layer snaps",
    t1 - t0, targets, success);
  if (a->hasLayer)
    checkLayerShape(a->mesh, "after snapping");
}

void visualizeGeometricInfo(Mesh* m, const char* name)
{
  Tag* dimensionTag = m->createIntTag("ma_geom_dim",1);
  Tag* idTag = m->createIntTag("ma_geom_id",1);
  apf::Field* field = apf::createLagrangeField(m,"ma_param",apf::VECTOR,1);
  apf::Field* targetField = apf::createLagrangeField(m,"ma_target",apf::VECTOR,1);
  Iterator* it = m->begin(0);
  Entity* v;
  while ((v = m->iterate(it)))
  {
    Model* c = m->toModel(v);
    int dimension = m->getModelType(c);
    m->setIntTag(v,dimensionTag,&dimension);
    int id = m->getModelTag(c);
    m->setIntTag(v,idTag,&id);
    Vector p;
    Vector xp = getPosition(m, v);
    m->getParam(v,p);
    if (dimension == 2 || dimension == 1) {
      Vector x;
      m->isParamPointInsideModel(c, &p[0], x);
      apf::setVector(targetField, v, 0, x - xp);
    }
    else {
      Vector x(0., 0., 0.);
      apf::setVector(targetField, v, 0, x);
    }
    apf::setVector(field,v,0,p);
  }
  m->end(it);
  apf::writeVtkFiles(name,m);
  it = m->begin(0);
  while ((v = m->iterate(it)))
  {
    m->removeTag(v,dimensionTag);
    m->removeTag(v,idTag);
  }
  m->end(it);
  m->destroyTag(dimensionTag);
  m->destroyTag(idTag);
  apf::destroyField(field);
  apf::destroyField(targetField);
}

}
