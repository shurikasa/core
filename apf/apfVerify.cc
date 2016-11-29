#include <PCU.h>
#include "apfMesh.h"
#include "apf.h"
#include <gmi.h>
#include <sstream>
#include <apfGeometry.h>
#include <cassert>

namespace apf {

static void intersect(
    std::set<int>& a,
    std::set<int> const& b)
{
  for (std::set<int>::iterator it = a.begin();
      it != a.end();)
  {
    if ( ! b.count(*it))
      a.erase(*(it++));
    else
      ++it;
  }
}

static Parts getCandidateParts(Mesh* m, MeshEntity* e)
{
  Downward d;
  int nd = m->getDownward(e, getDimension(m, e) - 1, d);
  Parts c;
  m->getResidence(d[0], c);
  for (int i = 1; i < nd; ++i)
  {
    Parts dr;
    m->getResidence(d[i], dr);
    intersect(c, dr);
  }
  return c;
}

static bool isSubset(Parts const& a, Parts const& b)
{
  APF_ITERATE(Parts, a, it)
    if (!b.count(*it))
      return false;
  return true;
}

static void verifyDown(Mesh* m, MeshEntity* e, int gd, int ed)
{
  Downward d;
  int nd = m->getDownward(e, ed - 1, d);
  for (int i = 0; i < nd; ++i)
  {
    int dgd = m->getModelType(m->toModel(d[i]));
    assert(dgd <= gd);
  }
  Parts r;
  m->getResidence(e, r);
  assert(isSubset(r, getCandidateParts(m, e)));
}

typedef std::map<ModelEntity*, int> UpwardCounts;
typedef std::map<ModelEntity*, bool> SideManifoldness;

static void getUpwardCounts(gmi_model* gm, int meshDimension, UpwardCounts& uc)
{
  for (int d = 0; d < meshDimension; ++d) {
    gmi_iter* it = gmi_begin(gm, d);
    gmi_ent* ge;
    while ((ge = gmi_next(gm, it))) {
      gmi_set* up = gmi_adjacent(gm, ge, d + 1);
      uc[(ModelEntity*)ge] = up->n;
      gmi_free_set(up);
    }
    gmi_end(gm, it);
  }
}

static void verifyUp(Mesh* m, UpwardCounts& guc,
    MeshEntity* e)
{
  apf::Up up;
  m->getUp(e, up);
  int upwardCount = up.n;
  int meshDimension = m->getDimension();
  int entityDimension = getDimension(m, e);
  int difference = meshDimension - entityDimension;
  assert(difference);
  ModelEntity* ge = m->toModel(e);
  int modelDimension = m->getModelType(ge);
  int modelUpwardCount = guc[ge];
  bool isOnNonManifoldFace = (modelDimension == meshDimension - 1) &&
                             (modelUpwardCount > 1);
  bool isOnManifoldBoundary = ( ! isOnNonManifoldFace) &&
                              (modelDimension < meshDimension);
  bool isShared = m->isShared(e);
  bool isMatched = false;
  if (m->hasMatching()) {
    apf::Matches matches;
    m->getMatches(e, matches);
    isMatched = (matches.getSize() > 0);
  }
  bool isExposedByModel = isOnManifoldBoundary;
  bool isExposedByMesh = isShared || isMatched;
  bool isExposed = isExposedByModel || isExposedByMesh;
  bool isOnEqualOrder = (modelDimension == entityDimension);
  int expected;
  if (isExposed)
    expected = difference;
  else
    expected = difference + 1;
  if (!isExposedByMesh && isOnEqualOrder)
    expected = std::max(expected, modelUpwardCount);
  bool okay;
  if (difference == 1)
    okay = (upwardCount == expected);
  else
    okay = (upwardCount >= expected);
  if ( ! okay) {
    std::stringstream ss;
    char const* n = apf::Mesh::typeName[m->getType(e)];
    ss << "apf::Verify: " << n << " with " << upwardCount << " adjacent "
       << dimName[entityDimension + 1] << "s\n";
    ss << "centroid: " << apf::getLinearCentroid(m, e) << '\n';
    ss << "based on the following:\n";
    ss << " - "  << n << " is classified on a model "
       << dimName[modelDimension] << '\n';
    if (modelDimension < meshDimension) {
      ss << "   "  << "which is adjacent to " << modelUpwardCount
         << " model " << dimName[modelDimension + 1] << "s\n";
    }
    if (isOnNonManifoldFace)
      ss << "   " << "and is a non-manifold " << dimName[modelDimension] << '\n';
    if (isOnManifoldBoundary)
      ss << "   " << "and is a manifold boundary " << dimName[modelDimension] << '\n';
    if (isExposedByModel)
      ss << "   making the " << n << " \"exposed\" at the geometric boundary\n";
    if (isShared)
      ss << " - " << n << " is shared in parallel\n";
    if (isMatched)
      ss << " - " << n << " has periodic matches\n";
    if (isExposedByMesh)
      ss << "   making it \"exposed\" at the part boundary\n";
    ss << "we would expect the adjacent " << dimName[entityDimension + 1] << " count to be ";
    if (difference == 1)
      ss << "exactly " << expected << '\n';
    else
      ss << "at least " << expected << '\n';

    printf("The problematic model edge number is %d\n",m->getModelTag(m->toModel(e)));
    apf::Downward de;
    int nde = m->getDownward(e, 0, de);
    printf("Vertices are:\n");
    for (int i=0; i < nde; ++i){
        Vector3 vec;
        m->getPoint(de[i], 0, vec);
        printf("vertex %d: ", i);
        for (int j = 0; j < 3; ++j){
            printf("%lf ",vec[j]);
        }
        printf("\n\n");
    }

    std::string s = ss.str();
    fail(s.c_str());
  }
  /* this is here for some spiderwebby simmetrix meshes */
  if (upwardCount >= 200) {
    std::stringstream ss;
    ss << "warning: entity of type " << m->getType(e)
      << " at " << getLinearCentroid(m, e) << " has "
      << upwardCount << " upward adjacencies\n";
    // use a stringstream to prevent output from different procs mixing
    std::string s = ss.str();
    fprintf(stderr,"%s",s.c_str());
  }
}

static void verifyResidence(Mesh* m, MeshEntity* e)
{
  Parts p;
  m->getResidence(e, p);
  assert(p.count(m->getOwner(e)));
  Copies r;
  m->getRemotes(e, r);
  assert(r.size() + 1 == p.size());
  APF_ITERATE(Copies, r, it)
    assert(p.count(it->first));
}

static void verifyEntity(Mesh* m, UpwardCounts& guc, MeshEntity* e)
{
  int ed = getDimension(m, e);
  int md = m->getDimension();
  assert(md >= ed);
  int gd = m->getModelType(m->toModel(e));
  assert(gd >= ed);
  if (ed)
    verifyDown(m, e, gd, ed);
  if (ed < md)
    verifyUp(m, guc, e);
  verifyResidence(m, e);
}

static Copies getAllCopies(Mesh* m, MeshEntity* e)
{
  Copies c;
  m->getRemotes(e, c);
  c[PCU_Comm_Self()] = e;
  return c;
}

static void verifyAllCopies(Copies& a)
{
  APF_ITERATE(Copies, a, it)
  {
    assert(it->first >= 0);
    assert(it->first < PCU_Comm_Peers());
  }
}

static void packCopies(
    int to,
    Copies& copies)
{
  int n = copies.size();
  PCU_COMM_PACK(to,n);
  APF_ITERATE(Copies,copies,it)
  {
    PCU_COMM_PACK(to,it->first);
    PCU_COMM_PACK(to,it->second);
  }
}

static void unpackCopies(
    Copies& copies)
{
  int n;
  PCU_COMM_UNPACK(n);
  for (int i=0; i < n; ++i)
  {
    int part;
    PCU_COMM_UNPACK(part);
    MeshEntity* remote;
    PCU_COMM_UNPACK(remote);
    assert(remote);
    copies[part]=remote;
  }
}

static void sendAllCopies(Mesh* m, MeshEntity* e)
{
  Copies a;
  Copies r;
  a = getAllCopies(m, e);
  verifyAllCopies(a);
  m->getRemotes(e, r);
  assert(!r.count(PCU_Comm_Self()));
  APF_ITERATE(Copies, r, it)
  {
    PCU_COMM_PACK(it->first, it->second);
    packCopies(it->first, a);
  }
}

static void receiveAllCopies(Mesh* m)
{
  MeshEntity* e;
  PCU_COMM_UNPACK(e);
  Copies a;
  unpackCopies(a);
  Copies b = getAllCopies(m, e);
  assert(a == b);
}

static void verifyRemoteCopies(Mesh* m)
{
  PCU_Comm_Begin();
  for (int d = 0; d <= m->getDimension(); ++d)
  {
    MeshIterator* it = m->begin(d);
    MeshEntity* e;
    while ((e = m->iterate(it)))
      if (m->isShared(e))
        sendAllCopies(m, e);
    m->end(it);
  }
  PCU_Comm_Send();
  while (PCU_Comm_Receive())
    receiveAllCopies(m);
}

static bool hasMatch(
    Matches& matches,
    int peer,
    MeshEntity* entity)
{
  unsigned i;
  for (i = 0; i < matches.getSize(); ++i)
    if (matches[i].peer == peer &&
        matches[i].entity == entity)
      return true;
  return false;
}

static void sendSelfToMatches(MeshEntity* e, Matches& matches)
{
  APF_ITERATE(Matches, matches, it)
  {
    assert(!((it->peer == PCU_Comm_Self())&&(it->entity == e)));
    PCU_COMM_PACK(it->peer, e);
    PCU_COMM_PACK(it->peer, it->entity);
  }
}

static void receiveMatches(Mesh* m)
{
  MeshEntity* source;
  PCU_COMM_UNPACK(source);
  MeshEntity* e;
  PCU_COMM_UNPACK(e);
  Matches matches;
  m->getMatches(e, matches);
  assert(hasMatch(matches, PCU_Comm_Sender(), source));
}

static void verifyMatches(Mesh* m)
{
  PCU_Comm_Begin();
  for (int d = 0; d <= m->getDimension(); ++d)
  {
    MeshIterator* it = m->begin(d);
    MeshEntity* e;
    while ((e = m->iterate(it))) {
      Matches matches;
      m->getMatches(e, matches);
      sendSelfToMatches(e, matches);
    }
    m->end(it);
  }
  PCU_Comm_Send();
  while (PCU_Comm_Receive())
    receiveMatches(m);
}

static void sendCoords(Mesh* m, MeshEntity* e)
{
  Vector3 x;
  m->getPoint(e, 0, x);
  Vector3 p(0,0,0);
  m->getParam(e, p);
  Copies r;
  m->getRemotes(e, r);
  APF_ITERATE(Copies, r, it)
  {
    PCU_COMM_PACK(it->first, it->second);
    PCU_COMM_PACK(it->first, x);
    PCU_COMM_PACK(it->first, p);
  }
}

static bool receiveCoords(Mesh* m)
{
  MeshEntity* e;
  Vector3 ox;
  Vector3 op;
  PCU_COMM_UNPACK(e);
  PCU_COMM_UNPACK(ox);
  PCU_COMM_UNPACK(op);
  Vector3 x;
  Vector3 p(0,0,0);
  m->getPoint(e, 0, x);
  m->getParam(e, p);
  return areClose(x, ox, 0.0) &&
         areClose(p, op, 0.0);
}

static long verifyCoords(Mesh* m)
{
  PCU_Comm_Begin();
  MeshIterator* it = m->begin(0);
  MeshEntity* e;
  while ((e = m->iterate(it)))
    if (m->isShared(e))
      sendCoords(m, e);
  m->end(it);
  PCU_Comm_Send();
  long n = 0;
  while (PCU_Comm_Receive())
    if (!receiveCoords(m))
      ++n;
  return PCU_Add_Long(n);
}

long verifyVolumes(Mesh* m, bool printVolumes)
{
  MeshIterator* it = m->begin(m->getDimension());
  MeshEntity* e;
  long n = 0;
  while ((e = m->iterate(it)))
  {
    if (!isSimplex(m->getType(e)))
      continue;
    double v = measure(m,e);
    if (v < 0) {
      if (printVolumes) {
        std::stringstream ss;
        ss << "warning: element volume " << v
          << " at " << getLinearCentroid(m, e) << '\n';
        std::string s = ss.str();
        fprintf(stdout, "%s", s.c_str());
        fflush(stdout);
      }
      ++n;
    }
  }
  m->end(it);
  return PCU_Add_Long(n);
}

static void packAlignment(Mesh* m, MeshEntity* e, MeshEntity* r, int to)
{
  PCU_COMM_PACK(to,r);
  int d = getDimension(m, e);
  Downward down;
  int nd = m->getDownward(e, d - 1, down);
  for (int i = 0; i < nd; ++i) {
    Copies remotes;
    m->getRemotes(down[i], remotes);
    MeshEntity* dr = remotes[to];
    PCU_COMM_PACK(to,dr);
  }
}

static void sendAlignment(Mesh* m, MeshEntity* e)
{
  Copies remotes;
  m->getRemotes(e, remotes);
  APF_ITERATE(Copies, remotes, it)
    packAlignment(m, e, it->second, it->first);
}

static void receiveAlignment(Mesh* m)
{
  MeshEntity* e;
  PCU_COMM_UNPACK(e);
  int d = getDimension(m, e);
  Downward down;
  int nd = m->getDownward(e, d - 1, down);
  for (int i = 0; i < nd; ++i) {
    PCU_COMM_UNPACK(e);
    assert(down[i] == e);
  }
}

static void verifyAlignment(Mesh* m)
{
  PCU_Comm_Begin();
  for (int d = 1; d <= m->getDimension(); ++d)
  {
    MeshIterator* it = m->begin(d);
    MeshEntity* e;
    while ((e = m->iterate(it)))
      if (m->isShared(e))
        sendAlignment(m, e);
    m->end(it);
  }
  PCU_Comm_Send();
  while (PCU_Comm_Receive())
    receiveAlignment(m);
}

static void verifyTags(Mesh* m)
{
  PCU_Comm_Begin();
  DynamicArray<MeshTag*> tags;
  m->getTags(tags);
  int self = PCU_Comm_Self();
  if (self) {
    int n = tags.getSize();
    PCU_COMM_PACK(self - 1, n);
    for (int i = 0; i < n; ++i)
      packTagInfo(m, tags[i], self - 1);
  }
  PCU_Comm_Send();
  while (PCU_Comm_Receive()) {
    int n;
    PCU_COMM_UNPACK(n);
    assert(tags.getSize() == (size_t)n);
    for (int i = 0; i < n; ++i) {
      std::string name;
      int type;
      int size;
      unpackTagInfo(name, type, size);
      assert(name == m->getTagName(tags[i]));
      assert(type == m->getTagType(tags[i]));
      assert(size == m->getTagSize(tags[i]));
    }
  }
}

void verify(Mesh* m)
{
  double t0 = PCU_Time();
  verifyTags(m);
  UpwardCounts guc;
  getUpwardCounts(m->getModel(), m->getDimension(), guc);
  /* got to 3 on purpose, so we can verify if
     m->getDimension is lying */
  for (int d = 0; d <= 3; ++d)
  {
    MeshIterator* it = m->begin(d);
    MeshEntity* e;
    size_t n = 0;
    while ((e = m->iterate(it)))
    {
      verifyEntity(m, guc, e);
      ++n;
    }
    m->end(it);
    assert(n == m->count(d));
    if (d > m->getDimension())
      assert(!n);
  }
  guc.clear();
  verifyRemoteCopies(m);
  verifyAlignment(m);
  verifyMatches(m);
  long n = verifyCoords(m);
  if (n && (!PCU_Comm_Self()))
    fprintf(stderr,"apf::verify fail: %ld coordinate mismatches\n", n);
  n = verifyVolumes(m);
  if (n && (!PCU_Comm_Self()))
    fprintf(stderr,"apf::verify warning: %ld negative simplex elements\n", n);
  double t1 = PCU_Time();
  if (!PCU_Comm_Self())
    printf("mesh verified in %f seconds\n", t1 - t0);
}

}
