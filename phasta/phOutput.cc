#include <PCU.h>
#include "phOutput.h"
#include "phLinks.h"
#include "phAdjacent.h"
#include "phBubble.h"
#include "phAxisymmetry.h"
#include "phInterfaceCutter.h"
#include <fstream>
#include <sstream>
#include <cassert>

namespace ph {

static void getCounts(Output& o)
{
  o.nOwnedNodes = apf::countOwned(o.mesh, 0);
  o.nOverlapNodes = o.mesh->count(0);
}

static void getCoordinates(Output& o)
{
  apf::Mesh* m = o.mesh;
  int n = m->count(0);
  double* x = new double[n * 3];
  apf::MeshEntity* v;
  int i = 0;
  apf::MeshIterator* it = m->begin(0);
  while ((v = m->iterate(it))) {
    apf::Vector3 p;
    m->getPoint(v, 0, p);
    for (int j = 0; j < 3; ++j)
      x[j * n + i] = p[j]; /* FORTRAN indexing */
    ++i;
  }
  m->end(it);
  assert(i == n);
  o.arrays.coordinates = x;
}

/* so apparently old phParAdapt just used EN_id,
   and the id generator from pumi would do things
   like this. I guess PHASTA is ok with a unique
   number for each copy, regardless of part boundary
   sharing...
 
update: Michel says these global numbers are ignored
        by phasta. get rid of them when you can.
 */
static void getGlobal(Output& o)
{
  apf::Mesh* m = o.mesh;
  int n = m->count(0);
  int self = PCU_Comm_Self();
  int peers = PCU_Comm_Peers();
  int id = self + 1;
  o.arrays.globalNodeNumbers = new int[n];
  for (int i = 0; i < n; ++i) {
    o.arrays.globalNodeNumbers[i] = id;
    id += peers;
  }
}

static void getVertexLinks(Output& o, apf::Numbering* n, BCs& bcs)
{
  Links links;
  getLinks(o.mesh, 0, links, bcs);
  encodeILWORK(n, links, o.nlwork, o.arrays.ilwork);
}

static void getInterior(Output& o, BCs& bcs, apf::Numbering* n)
{
  apf::Mesh* m = o.mesh;
  Blocks& bs = o.blocks.interior;
  int*** ien     = new int**[bs.getSize()];
  int**  mattype = 0;
  if (bcs.fields.count("material type")) 
    mattype = new int* [bs.getSize()];
  apf::NewArray<int> js(bs.getSize());
  for (int i = 0; i < bs.getSize(); ++i) {
    ien    [i] = new int*[bs.nElements[i]];
    if (mattype) 
      mattype[i] = new int [bs.nElements[i]];
    js[i] = 0;
  }
  gmi_model* gm = m->getModel();
  apf::MeshEntity* e;
  apf::MeshIterator* it = m->begin(m->getDimension());
  while ((e = m->iterate(it))) {
    BlockKey k;
    getInteriorBlockKey(m, e, k);
    int nv = k.nElementVertices;
    assert(bs.keyToIndex.count(k));
    int i = bs.keyToIndex[k];
    int j = js[i];
    ien[i][j] = new int[nv];
    apf::Downward v;
    getVertices(m, e, v);
    for (int k = 0; k < nv; ++k)
      ien[i][j][k] = apf::getNumber(n, v[k], 0, 0);

    /* get material type */
    if (mattype) {
      gmi_ent* ge = (gmi_ent*)m->toModel(e);
      apf::Vector3 x;
      //m->getPoint(e, 0, x);
      x = apf::getLinearCentroid(m, e);
      std::string s("material type");
      FieldBCs& fbcs = bcs.fields[s];
      double* matval = getBCValue(gm, fbcs, ge, x);
      mattype[i][j] = *matval;
    }
    ++js[i];
  }
  m->end(it);
  for (int i = 0; i < bs.getSize(); ++i)
    assert(js[i] == bs.nElements[i]);
  o.arrays.ien     = ien;
  o.arrays.mattype = mattype;
}

static void getBoundary(Output& o, BCs& bcs, apf::Numbering* n)
{
  apf::Mesh* m = o.mesh;
  gmi_model* gm = m->getModel();
  int nbc = countNaturalBCs(*o.in);
  Blocks& bs = o.blocks.boundary;
  int*** ienb = new int**[bs.getSize()];
  int**  mattypeb = 0;
  if (bcs.fields.count("material type")) 
    mattypeb = new int*[bs.getSize()];
  int*** ibcb = new int**[bs.getSize()];
  double*** bcb = new double**[bs.getSize()];
  apf::NewArray<int> js(bs.getSize());
  for (int i = 0; i < bs.getSize(); ++i) {
    ienb[i]     = new int*[bs.nElements[i]];
    if (mattypeb) 
      mattypeb[i] = new int [bs.nElements[i]];
    ibcb[i]     = new int*[bs.nElements[i]];
    bcb[i]      = new double*[bs.nElements[i]];
    js[i] = 0;
  }
  int boundaryDim = m->getDimension() - 1;
  apf::MeshEntity* f;
  apf::MeshIterator* it = m->begin(boundaryDim);
  while ((f = m->iterate(it))) {
    apf::ModelEntity* me = m->toModel(f);
    if (m->getModelType(me) != boundaryDim)
      continue;
    apf::Matches matches;
    m->getMatches(f, matches);
    if (matches.getSize() == 1) // This prevents adding interface elements...
      continue;
    if (m->countUpward(f)>1)   // don't want interior region boundaries here...
      continue;
    gmi_ent* gf = (gmi_ent*)me;
    apf::MeshEntity* e = m->getUpward(f, 0);
    BlockKey k;
    getBoundaryBlockKey(m, e, f, k);
    assert(bs.keyToIndex.count(k));
    int i = bs.keyToIndex[k];
    int j = js[i];
    int nv = k.nElementVertices;
    apf::Downward v;
    getBoundaryVertices(m, e, f, v);
    ienb[i][j] = new int[nv];
    for (int k = 0; k < nv; ++k)
      ienb[i][j][k] = apf::getNumber(n, v[k], 0, 0);
    bcb[i][j] = new double[nbc]();
    ibcb[i][j] = new int[2](); /* <- parens initialize to zero */
    apf::Vector3 x = apf::getLinearCentroid(m, f);
    applyNaturalBCs(gm, gf, bcs, x, bcb[i][j], ibcb[i][j]);

    /* get material type */
    if (mattypeb) {
      gmi_ent* ge = (gmi_ent*)m->toModel(e);
      x = apf::getLinearCentroid(m, e);
      std::string s("material type");
      FieldBCs& fbcs = bcs.fields[s];
      double* matvalb = getBCValue(gm, fbcs, ge, x);
      mattypeb[i][j] = *matvalb;
    }
    ++js[i];
  }
  m->end(it);
  for (int i = 0; i < bs.getSize(); ++i)
    assert(js[i] == bs.nElements[i]);
  o.arrays.ienb = ienb;
  o.arrays.mattypeb = mattypeb;
  o.arrays.ibcb = ibcb;
  o.arrays.bcb = bcb;
}

bool checkInterface(Output& o, BCs& bcs) {
  o.hasDGInterface = 0;
  apf::Mesh* m = o.mesh;
  gmi_model* gm = m->getModel();

  std::string name1("DG interface");
  if (!haveBC(bcs, name1))
    return false;
  FieldBCs& fbcs1 = bcs.fields[name1];

  std::string name2("material type");
  if (!haveBC(bcs, name2))
    return false;
  FieldBCs& fbcs2 = bcs.fields[name2];

  if (PCU_Comm_Self() == 0)
    printf("Run checkInterface! Turn on hasDGInterface!\n");

  o.hasDGInterface = 1;

  int a = 0; int b = 0;
  int aID = 0;
  int bID = 1;
  int matID = 0;
  int aIDSetFlag = 0;
  int bIDSetFlag = 0;
  apf::MeshIterator* it = m->begin(m->getDimension()-1);
  apf::MeshEntity* e;
  while ((e = m->iterate(it))) {
    gmi_ent* ge = (gmi_ent*) m->toModel(e);
    if (ph::isInterface(gm, ge, fbcs1)) {
      apf::MeshEntity* eUp = m->getUpward(e, 0);
      gmi_ent* geUp = (gmi_ent*) m->toModel(eUp);
      apf::Vector3 x = apf::getLinearCentroid(m, eUp);
      double* floatID = getBCValue(gm, fbcs2, geUp, x);
      matID = (int)(*floatID+0.5);
      if (aIDSetFlag == 0) {
        aID = matID;
        aIDSetFlag = 1;
      } else if (bIDSetFlag == 0 && matID != aID) {
        bID = matID;
        bIDSetFlag = 1;
      }
      if ( matID == aID ) a++;
      if ( matID == bID ) b++;
    }
  }
  m->end(it);
  assert(aID!=bID); //assert different material ID on two sides
  assert(a==b); //assert same number of faces on each side
  return true;
}

static void getInterface
(
  Output&         o,
  BCs&            bcs,
  apf::Numbering* n
)
{
  apf::Mesh*        m  = o.mesh;
  gmi_model*        gm = m->getModel();
  BlocksInterface&  bs = o.blocks.interface;
  int***            ienif0 = new int**[bs.getSize()];
  int***            ienif1 = new int**[bs.getSize()];
  int**             mattypeif0 = 0;
  int**             mattypeif1 = 0;
  if (bcs.fields.count("material type")) {
    mattypeif0 = new int*[bs.getSize()];
    mattypeif1 = new int*[bs.getSize()];
  }
  apf::NewArray<int> js(bs.getSize());
  for (int i = 0; i < bs.getSize(); ++i) {
    ienif0[i] = new int*[bs.nElements[i]];
    ienif1[i] = new int*[bs.nElements[i]];
    if (mattypeif0) mattypeif0[i] = new int [bs.nElements[i]];
    if (mattypeif1) mattypeif1[i] = new int [bs.nElements[i]];
    js[i] = 0;
  }
  int interfaceDim = m->getDimension() - 1;
  apf::MeshEntity*   face;
  apf::MeshIterator* it = m->begin(interfaceDim);

  while ((face = m->iterate(it))) {
    apf::ModelEntity* me = m->toModel(face);
    if (getBCValue(m->getModel(), bcs.fields["DG interface"], (gmi_ent*) me) == 0)
      continue;
    if (m->getModelType(me) != interfaceDim)
      continue;
    apf::Matches matches;
    m->getMatches(face, matches);
    assert(matches.getSize() == 1);
    apf::MeshEntity* e0 = m->getUpward(face, 0);
    apf::MeshEntity* e1 = m->getUpward(matches[0].entity, 0);
    /* in order to avoid repeatation of elements */
    if (e0 > e1)
      continue;

    BlockKeyInterface k;
    getInterfaceBlockKey(m, e0, e1, face, k);
    assert(bs.keyToIndex.count(k));
    int i = bs.keyToIndex[k];
    int j = js[i];
    int nv0 = k.nElementVertices;
    int nv1 = k.nElementVertices1;
    apf::Downward v0, v1;
    getBoundaryVertices(m, e0, face, v0);
    getBoundaryVertices(m, e1, matches[0].entity, v1);
    ienif0[i][j] = new int[nv0];
    ienif1[i][j] = new int[nv1];
    for (int k = 0; k < nv0; ++k) 
      ienif0[i][j][k] = apf::getNumber(n, v0[k], 0, 0);
    for (int k = 0; k < nv1; ++k)
      ienif1[i][j][k] = apf::getNumber(n, v1[k], 0, 0);

    /* get material type */
    if (mattypeif0) {
      gmi_ent* ge0 = (gmi_ent*)m->toModel(e0);
      gmi_ent* ge1 = (gmi_ent*)m->toModel(e1);
      apf::Vector3 x0;
      apf::Vector3 x1;
      x0 = apf::getLinearCentroid(m, e0);
      x1 = apf::getLinearCentroid(m, e1);
      std::string s("material type");
      FieldBCs& fbcs = bcs.fields[s];
      double* matvalif0 = getBCValue(gm, fbcs, ge0, x0);
      double* matvalif1 = getBCValue(gm, fbcs, ge1, x1);
      mattypeif0[i][j] = *matvalif0;
      mattypeif1[i][j] = *matvalif1;
    }
    ++js[i];
  }
  m->end(it);
  for (int i = 0; i < bs.getSize(); ++i)
    assert(js[i] == bs.nElements[i]);
  o.arrays.ienif0 = ienif0;
  o.arrays.ienif1 = ienif1;
  o.arrays.mattypeif0 = mattypeif0;
  o.arrays.mattypeif1 = mattypeif1;
}  

static void getBoundaryElements(Output& o)
{
  Blocks& bs = o.blocks.boundary;
  int n = 0;
  for (int i = 0; i < bs.getSize(); ++i)
    n += bs.nElements[i];
  o.nBoundaryElements = n;
}

static void getInterfaceElements(Output& o)
{
  BlocksInterface& bs = o.blocks.interface;
  int n = 0;
  for (int i = 0; i < bs.getSize(); ++i)
    n += bs.nElements[i]; /* need to add nElementsOther as well ??? */
  o.nInterfaceElements = n;
}

static void getMaxElementNodes(Output& o)
{
  int n = 0;
  Blocks& ibs = o.blocks.interior;
  for (int i = 0; i < ibs.getSize(); ++i)
    n = std::max(n, ibs.keys[i].nElementVertices);
  Blocks& bbs = o.blocks.boundary;
  for (int i = 0; i < bbs.getSize(); ++i)
    n = std::max(n, bbs.keys[i].nElementVertices);
  BlocksInterface ifbs = o.blocks.interface;
  for (int i = 0; i < ifbs.getSize(); ++i)
    n = std::max(n, ifbs.keys[i].nElementVertices);
  o.nMaxElementNodes = n;
}

/* returns the global periodic master iff it is on this
   part, otherwise returns e */
static apf::MeshEntity* getLocalPeriodicMaster(apf::MatchedSharing* sh,
    apf::MeshEntity* e)
{
  if ( ! sh)
    return e;
  apf::Copy globalMaster = sh->getOwner(e);
  if (globalMaster.peer == PCU_Comm_Self())
    return globalMaster.entity;
  else
    return e;
}

static void getLocalPeriodicMasters(Output& o, apf::Numbering* n, BCs& bcs)
{
  apf::Mesh* m = o.mesh;
  int* iper = new int[m->count(0)];
  apf::MeshIterator* it = m->begin(0);
  apf::MeshEntity* e;
  apf::MatchedSharing* sh = m->hasMatching() ? new apf::MatchedSharing(m) : 0;
  int i = 0;
  while ((e = m->iterate(it))) {
    apf::ModelEntity* me = m->toModel(e);
    bool isDG = ph::isInterface(m->getModel(),(gmi_ent*) me,bcs.fields["DG interface"]);
    apf::MeshEntity* master = getLocalPeriodicMaster(sh, e);
    if (master == e || isDG)
      iper[i] = 0;
    else
      iper[i] = apf::getNumber(n, master, 0, 0) + 1;
    ++i;
  }
  m->end(it);
  o.arrays.iper = iper;
  delete sh;
}

static bool isMatchingSlave(apf::MatchedSharing* ms, apf::MeshEntity* v)
{
  if (!ms)
    return false;
  apf::Matches matches;
  ms->mesh->getMatches(v, matches);
  if (!matches.getSize())
    return false;
  return !ms->isOwned(v);
}

static void getEssentialBCs(BCs& bcs, Output& o)
{
  Input& in = *o.in;
  apf::Mesh* m = o.mesh;
  apf::MeshTag* angles = 0;
  apf::MatchedSharing* ms = 0;
  if (m->hasMatching())
    ms = new apf::MatchedSharing(m);
  if (in.axisymmetry)
    angles = tagAngles(m, bcs, ms);
  int nv = m->count(0);
  o.arrays.nbc = new int[nv];
  for (int i = 0; i < nv; ++i)
    o.arrays.nbc[i] = 0;
  o.arrays.ibc = new int[nv]();
  o.arrays.bc = new double*[nv];
  o.nEssentialBCNodes = 0;
  int ibc;
  int nec = countEssentialBCs(in);
  double* bc = new double[nec]();
  gmi_model* gm = m->getModel();
  int i = 0;
  int& ei = o.nEssentialBCNodes;
  apf::MeshEntity* v;
  apf::MeshIterator* it = m->begin(0);
  while ((v = m->iterate(it))) {
    gmi_ent* ge = (gmi_ent*) m->toModel(v);
    apf::Vector3 x;
    m->getPoint(v, 0, x);
    ibc = 0;
    for (int j = 0; j < nec; ++j)
      bc[j] = 0;
    //bool hasBC = applyEssentialBCs(gm, ge, bcs, x, bc, &ibc) && bcs.fields.count("material type");
    bool hasBC = applyEssentialBCs(gm, ge, bcs, x, bc, &ibc);
    /* matching introduces an iper bit */
    /* which is set for all slaves */
    if (isMatchingSlave(ms, v)) {
      hasBC = true;
      ibc |= (1<<10);
      /* axisymmetric theta for some slaves */
      if (in.axisymmetry && m->hasTag(v, angles))
        m->getDoubleTag(v, angles, &bc[11]);
    }
    if (hasBC) {
      o.arrays.nbc[i] = ei + 1;
      o.arrays.ibc[ei] = ibc;
      double* bc_ei = new double[nec];
      for (int j = 0; j < nec; ++j)
        bc_ei[j] = bc[j];
      o.arrays.bc[ei] = bc_ei;
      ++ei;
    }
    ++i;
  }
  m->end(it);
  delete [] bc;
  if (in.axisymmetry)
    m->destroyTag(angles);
  delete ms;
}

static void getInitialConditions(BCs& bcs, Output& o)
{
  Input& in = *o.in;
  if (in.solutionMigration) {
    if (!PCU_Comm_Self())
      printf("All attribute-based initial conditions, "
             "if any, "
             "are ignored due to request for SolutionMigration\n");
    return;
  }
  apf::Mesh* m = o.mesh;
  apf::NewArray<double> s(in.ensa_dof);
  apf::NewArray<double> matValue(1);
  apf::Field* f = m->findField("solution");
  assert(f);
  apf::MeshIterator* it = m->begin(3);
  apf::MeshEntity* e;
  gmi_model* gm = m->getModel();
  while ((e = m->iterate(it))) {
    gmi_ent* ge = (gmi_ent*)m->toModel(e);
    apf::Downward v;
    int nv = m->getDownward(e, 0, v);
    for (int i = 0; i < nv; ++i) {
      apf::getComponents(f, v[i], 0, &s[0]);
      apf::Vector3 x;
      m->getPoint(v[i], 0, x);
      applySolutionBCs(gm, ge, bcs, x, &s[0]);
      apf::setComponents(f, v[i], 0, &s[0]);
    }
  }
  m->end(it);
}

static void getElementGraph(Output& o, apf::Numbering* rn, BCs& bcs)
{
  if (o.in->formElementGraph) {
    o.arrays.ienneigh = formIENNEIGH(rn);
    Links links;
    getLinks(o.mesh, o.mesh->getDimension() - 1, links, bcs);
    encodeILWORKF(rn, links, o.nlworkf, o.arrays.ilworkf);
  } else {
    o.arrays.ilworkf = 0;
    o.arrays.ienneigh = 0;
  }
}

static void getEdges(Output& o, apf::Numbering* vn, apf::Numbering* rn, BCs& bcs)
{
  if (o.in->formEdges) {
    Links links;
    getLinks(o.mesh, 1, links, bcs);
    apf::Numbering* en = apf::numberOverlapDimension(o.mesh, "ph::getEdges", 1);
    encodeILWORK(en, links, o.nlworkl, o.arrays.ilworkl);
    apf::destroyNumbering(en);
  } else {
    o.arrays.ilworkl = 0;
  }
  if (o.in->formEdges) {
    apf::Mesh* m = o.mesh;
    assert(m->getDimension() == 3);
    int nelems = m->count(3);
    o.arrays.iel = new int[nelems * 6];
    apf::MeshIterator* it = m->begin(3);
    apf::MeshEntity* e;
    int i = 0;
    while ((e = m->iterate(it))) {
      apf::MeshEntity* ev[6];
      m->getDownward(e, 0, ev);
      for (int j = 0; j < 6; ++j)
        o.arrays.iel[j * nelems + i] = apf::getNumber(vn, ev[j], 0, 0) + 1;
      ++i;
    }
    m->end(it);
    assert(i == nelems);
  } else {
    o.arrays.iel = 0;
  }
  if (o.in->formEdges) {
    apf::Mesh* m = o.mesh;
    int nelems = m->count(3);
    int nedges = m->count(1);
    o.arrays.ileo = new int[nedges + 1];
    o.arrays.ile = new int[nelems * 6];
    apf::MeshIterator* it = m->begin(1);
    apf::MeshEntity* e;
    int i = 0;
    o.arrays.ileo[0] = 0;
    while ((e = m->iterate(it))) {
      apf::Adjacent adj;
      m->getAdjacent(e, 3, adj);
      int k = o.arrays.ileo[i];
      for (size_t j = 0; j < adj.getSize(); ++j)
        o.arrays.ile[k++] = apf::getNumber(rn, adj[j], 0, 0) + 1;
      o.arrays.ileo[i + 1] = k;
      ++i;
    }
    m->end(it);
    assert(i == nedges);
  } else {
    o.arrays.ileo = 0;
    o.arrays.ile = 0;
  }
}

Output::~Output()
{
  delete [] arrays.coordinates;
  delete [] arrays.ilwork;
  delete [] arrays.ilworkf;
  delete [] arrays.iper;
  delete [] arrays.globalNodeNumbers;
  Blocks& ibs = blocks.interior;
  for (int i = 0; i < ibs.getSize(); ++i) {
    for (int j = 0; j < ibs.nElements[i]; ++j)
      delete [] arrays.ien    [i][j];
    delete [] arrays.ien    [i];
    if (arrays.mattype) delete [] arrays.mattype[i];
  }
  delete [] arrays.ien;
  if (arrays.mattype) delete [] arrays.mattype;
  Blocks& bbs = blocks.boundary;
  for (int i = 0; i < bbs.getSize(); ++i) {
    for (int j = 0; j < bbs.nElements[i]; ++j) {
      delete [] arrays.ienb[i][j];
      delete [] arrays.ibcb[i][j];
      delete [] arrays.bcb[i][j];
    }
    delete [] arrays.ienb[i];
    delete [] arrays.ibcb[i];
    delete [] arrays.bcb[i];
    if (arrays.mattypeb) delete [] arrays.mattypeb[i];
  }
  delete [] arrays.ienb;
  delete [] arrays.ibcb;
  delete [] arrays.bcb;
  delete [] arrays.nbc;
  delete [] arrays.ibc;
  if (arrays.mattypeb) delete [] arrays.mattypeb;
  for (int i = 0; i < nEssentialBCNodes; ++i)
    delete [] arrays.bc[i];
  delete [] arrays.bc;
  delete [] arrays.ienneigh;
  BlocksInterface& ifbs = blocks.interface;
  for (int i = 0; i < ifbs.getSize(); ++i) {
    for (int j = 0; j < ifbs.nElements[i]; ++j) {
      delete [] arrays.ienif0[i][j];
      delete [] arrays.ienif1[i][j];
    }
    delete [] arrays.ienif0[i];
    delete [] arrays.ienif1[i];
    if (arrays.mattypeif0) delete [] arrays.mattypeif0[i];
    if (arrays.mattypeif1) delete [] arrays.mattypeif1[i];
  }
  delete [] arrays.ienif0;
  delete [] arrays.ienif1;
  if (arrays.mattypeif0) delete [] arrays.mattypeif0;
  if (arrays.mattypeif1) delete [] arrays.mattypeif1;
  delete [] arrays.ilworkl;
  delete [] arrays.iel;
  delete [] arrays.ileo;
  delete [] arrays.ile;
}

void generateOutput(Input& in, BCs& bcs, apf::Mesh* mesh, Output& o)
{
  double t0 = PCU_Time();
  o.in = &in;
  o.mesh = mesh;
  getCounts(o);
  getCoordinates(o);
  getGlobal(o);
  getAllBlocks(o.mesh, o.blocks);
  apf::Numbering* n = apf::numberOverlapNodes(mesh, "ph_local");
  apf::Numbering* rn = apf::numberElements(o.mesh, "ph_elem");
  getVertexLinks(o, n, bcs);
  getInterior(o, bcs, n);
  getBoundary(o, bcs, n);
  getInterface(o, bcs, n);
  checkInterface(o,bcs);
  getLocalPeriodicMasters(o, n, bcs);
  getEdges(o, n, rn, bcs);
  apf::destroyNumbering(n);
  getBoundaryElements(o);
  getInterfaceElements(o);
  getMaxElementNodes(o);
  getEssentialBCs(bcs, o);
  getInitialConditions(bcs, o);
  getElementGraph(o, rn, bcs);
  apf::destroyNumbering(rn);
  if (in.initBubbles)
    initBubbles(o.mesh, in);
  double t1 = PCU_Time();
  if (!PCU_Comm_Self())
    printf("generated output structs in %f seconds\n",t1 - t0);
}

}
