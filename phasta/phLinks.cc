#include <PCU.h>
#include "phLinks.h"
#include "phAdjacent.h"
#include <apf.h>
#include <cassert>

namespace ph {

bool LinkKey::operator<(LinkKey const& other) const
{
  if (send != other.send)
    return send;
  return peer < other.peer;
}

typedef std::map<apf::Parts,int> PtnMdl;
typedef std::map<int,int> Owners;

/* returns the current partition model based on NormalSharing */
void getPtnMdl(apf::Mesh* m, PtnMdl& pm) {
  apf::MeshIterator* it = m->begin(0);
  apf::MeshEntity* e;
  apf::Parts res;
  while ((e = m->iterate(it))) {
    if( m->isShared(e) ) {
      m->getResidence(e,res);
      pm[res] = -1;
    }
  }
  m->end(it);
  APF_ITERATE(PtnMdl,pm,it)
    pm[it->first] = *(it->first.begin());
}

void printPtnMdl(const char* key, PtnMdl& pm) {
  APF_ITERATE(PtnMdl,pm,it) {
    PCU_Debug_Print("%s pme ", key);
    for (apf::Parts::iterator pid = it->first.begin(); pid != it->first.end(); pid++) {
      PCU_Debug_Print("%d ", *pid);
    }
    PCU_Debug_Print(" owner %d\n", it->second);
  }
}

/* exchange the number of owned partition model entities with neighbors */
void exchangeOwners(PtnMdl& pm, Owners& own) {
  const int self = PCU_Comm_Self();
  APF_ITERATE(PtnMdl,pm,it) {
    if( it->second == self )
      own[self]++;
    for (apf::Parts::iterator pid = it->first.begin(); pid != it->first.end(); pid++) {
      if( *pid != self )
        own[*pid] = 0; //need to know that the peer exists
    }
  }
  PCU_Comm_Begin();
  //send the number of owned partition model ents
  APF_ITERATE(Owners,own,it)
    PCU_COMM_PACK(it->second, own[self]);
  PCU_Comm_Send();
  //recv the neighbors owned partition model ents counts
  while (PCU_Comm_Listen()) {
    int peer = PCU_Comm_Sender();
    int ownCnt = 0;
    PCU_COMM_UNPACK(ownCnt);
    own[peer] = ownCnt;
  }
}

/* send the new owner of the partition model entity to all parts
   bounded by the entity (the resident set 'res') */
void packRes(const apf::Parts& res, int newowner) {
  int* resArray = new int[res.size()];
  int i = 0;
  for (apf::Parts::iterator pid = res.begin(); pid != res.end(); pid++)
    resArray[i++] = *pid;
  for (apf::Parts::iterator pid = res.begin(); pid != res.end(); pid++) {
    int peer = *pid;
    int sz = res.size();
    PCU_COMM_PACK(peer, sz);
    PCU_Comm_Pack(peer, resArray, sizeof(int)*res.size());
    PCU_COMM_PACK(peer, newowner);
  }
  delete [] resArray;
}

void unpackRes(apf::Parts& res, int& newowner) {
  int size = 0;
  PCU_COMM_UNPACK(size);
  int* resArray = new int[size];
  PCU_Comm_Unpack(resArray, sizeof(int)*size);
  PCU_COMM_UNPACK(newowner);
  for(int i=0; i<size; i++)
    res.insert(resArray[i]);
  delete [] resArray;
}

int poorestNeighbor(const apf::Parts& res, Owners& own) {
  int poorest = PCU_Comm_Peers();
  int poorestid = -1;
  for (apf::Parts::iterator pid = res.begin(); pid != res.end(); pid++)
     if( own[*pid] < poorest ) {
       poorest = own[*pid];
       poorestid = *pid;
     }
  assert(poorestid != -1);
  return poorestid;
}

/* returns the partition model with balanced ownership */
void balanceOwners(PtnMdl& pm) {
  double t0 = PCU_Time();
  const int self = PCU_Comm_Self();
  PCU_Debug_Open();
  Owners own;
  exchangeOwners(pm,own);
  //DEBUG {
  printPtnMdl("init",pm);
  APF_ITERATE(Owners,own,it)
    PCU_Debug_Print("init peer %d owns %d\n", it->first, it->second);
  PCU_Debug_Print("init owned %d\n", own[self]);
  //DEBUG }
  int maxOwn = PCU_Max_Int(own[self]);
  int maxPeers = PCU_Max_Int(own.size()-1);
  int iter = 0;
  const int maxIter = 100;
  const int ownershipLimit = maxPeers/2;
  if( !PCU_Comm_Self() )
    fprintf(stderr, "maxOwned %d maxPeers %d ownershipLimit %d maxIter %d\n", 
        maxOwn, maxPeers, ownershipLimit, maxIter);
  while( maxOwn > ownershipLimit && iter++ < maxIter ) {
    PCU_Comm_Begin();
    //compute and send owner changes
    APF_ITERATE(PtnMdl,pm,it) {
      if( own[self] > ownershipLimit && it->second == self ) {
        int newowner = poorestNeighbor(it->first,own);
        PCU_Debug_Print("setting newowner %d\n", newowner);
        it->second = newowner;
        packRes(it->first,newowner);
        own[self]--;
        assert(own.count(newowner));
        own[newowner]++;
      }
    }
    PCU_Comm_Send();
    //recv the owner changes
    while (PCU_Comm_Listen()) {
      int peer = PCU_Comm_Sender();
      apf::Parts res;
      int newowner = 0;
      unpackRes(res,newowner);
      PCU_Debug_Print("getting newowner %d from peer %d\n", newowner, peer);
      pm[res] = newowner;
      own[peer]--;
      assert(own.count(newowner));
      own[newowner]++;
    }
    maxOwn = PCU_Max_Int(own[self]);
    PCU_Debug_Print("owned %d\n", own[self]);
    if( !PCU_Comm_Self() )
      fprintf(stderr, "maxOwned %d\n", maxOwn);
  }
  if( !PCU_Comm_Self() )
    fprintf(stderr, "maxOwned %d iter %d time %.3f\n", 
        maxOwn, iter, PCU_Time()-t0);
  //DEBUG {
  printPtnMdl("final",pm);
  APF_ITERATE(Owners,own,it)
    PCU_Debug_Print("final peer %d owns %d\n", it->first, it->second);
  PCU_Debug_Print("final owned %d\n", own[self]);
  //DEBUG }
}

struct BalancedSharing : public apf::Sharing
{
  BalancedSharing(apf::Mesh* m) {
    mesh = m;
    helper = apf::getSharing(m);
    getPtnMdl(mesh,pm);
    balanceOwners(pm);
  }
  ~BalancedSharing() {
    delete helper;
  }
  bool isOwned(apf::MeshEntity* e) {
    apf::Parts res;
    mesh->getResidence(e,res);
    return pm[res];
  }
  void getCopies(apf::MeshEntity* e,
      apf::CopyArray& copies)
  {
    helper->getCopies(e, copies);
  }
  private:
  apf::Mesh* mesh;
  apf::Sharing* helper;
  PtnMdl pm;
};
/* the PhastaSharing class is responsible for ensuring that
   ILWORK links matched entities correctly. */

struct PhastaSharing : public apf::Sharing {
  PhastaSharing(apf::Mesh* m)
  {
    mesh = m;
    helper = apf::getSharing(m);
    BalancedSharing bs(m);
  }
  ~PhastaSharing()
  {
    delete helper;
  }
  bool isOwned(apf::MeshEntity* e)
  {
    return helper->isOwned(e);
  }
  /* this will only be called for global masters */
  void getCopies(apf::MeshEntity* e,
      apf::CopyArray& copies)
  {
    helper->getCopies(e, copies);
    if ( ! mesh->hasMatching())
      return;
    /* filter out matches which are on the same part as the global master */
    int self = PCU_Comm_Self();
    size_t i = 0;
    for (size_t j = 0; j < copies.getSize(); ++j)
      if (copies[j].peer != self)
        copies[i++] = copies[j];
    copies.setSize(i);
  }
  apf::Mesh* mesh;
  apf::Sharing* helper;
};

/* this algorithm is essential to parallel
   scalability: generate local inter-part 
   communication arrays describing shared entities.

   for every partition model face, each part
   will store an array of all entities it has
   classified on that partition model face.
   the two arrays are aligned such that
   the same index in both arrays denotes the same vertex.

   phParAdapt had an implementation which 
   used arrays of size equal to the number of processors.

   This version is copied from code in MDS
   (and/or PUMI) that does the exact same thing
   to generate SMB files, but uses space and
   time proportional only to the neighborhood size,
   thanks in part to PCU algorithms

   In a perfect world, we should only have one copy
   of this code.
*/

void getLinks(apf::Mesh* m, int dim, Links& links)
{
  PhastaSharing shr(m);
  PCU_Comm_Begin();
  apf::MeshIterator* it = m->begin(dim);
  apf::MeshEntity* v;
  while ((v = m->iterate(it))) {
/* the alignment is such that the owner part's
   array follows the order of its vertex iterator
   traversal. The owner dictates the order to the
   other part by sending remote copies */
    if ( ! shr.isOwned(v))
      continue;
    apf::CopyArray remotes;
    shr.getCopies(v, remotes);
    for (size_t i = 0; i < remotes.getSize(); ++i) {
      /* in matching we may accumulate multiple occurrences
         of the same master in the outgoing links array
         to a part that contains multiple copies of it. */
      links[LinkKey(1, remotes[i].peer)].push_back(v);
      PCU_COMM_PACK(remotes[i].peer, remotes[i].entity);
    }
  }
  m->end(it);
  PCU_Comm_Send();
  while (PCU_Comm_Listen()) {
    int peer = PCU_Comm_Sender();
    while (!PCU_Comm_Unpacked()) {
      apf::MeshEntity* v;
      PCU_COMM_UNPACK(v);
      links[LinkKey(0, peer)].push_back(v);
    }
  }
}

/* encode the local links into a big array of integers
   per phParAdapt's ILWORK array format.

   ilwork[0] = total number of links
   following that, each link is serialized as:

   tag (used to be something like (peer*peers + self), now zero)
   type (0 for send, 1 for receive)
   peer
   #entities

   then for each entity there are two integers

   local id
   #dofs on entity
*/

void encodeILWORK(apf::Numbering* n, Links& links, int& size, int*& a)
{
  size = 1; // total links
  APF_ITERATE(Links, links, it) {
    size += 4; //link header
    size += it->second.size() * 2; //entity entries
  }
  a = new int[size];
  a[0] = links.size();
  int i = 1;
  APF_ITERATE(Links, links, it) {
    LinkKey k = it->first;
    Link& l = it->second;
    a[i++] = 0;
    a[i++] = k.send;
    a[i++] = k.peer + 1; /* peers numbered from 1 */
    a[i++] = l.size();
    APF_ITERATE(Link, l, lit) {
      /* entities also numbered from 1 */
      a[i++] = apf::getNumber(n, *lit, 0, 0) + 1;
      a[i++] = 1;
    }
  }
  assert(i == size);
}

static apf::MeshEntity* getSideElement(apf::Mesh* m, apf::MeshEntity* s)
{
  return m->getUpward(s, 0);
}

void encodeILWORKF(apf::Numbering* n, Links& links, int& size, int*& a)
{
  apf::Mesh* m = apf::getMesh(n);
  size = 1; // total links
  APF_ITERATE(Links, links, it) {
    size += 2; //link header
    size += it->second.size(); //entity entries
  }
  a = new int[size];
  a[0] = links.size();
  int i = 1;
  APF_ITERATE(Links, links, it) {
    LinkKey k = it->first;
    Link& l = it->second;
    a[i++] = k.peer + 1; /* peers numbered from 1 */
    a[i++] = l.size();
    APF_ITERATE(Link, l, lit) {
      /* entities also numbered from 1 */
      apf::MeshEntity* e = getSideElement(m, *lit);
      a[i++] = apf::getNumber(n, e, 0, 0) + 1;
    }
  }
  assert(i == size);
}

static apf::MeshEntity* getOtherElem(apf::Mesh* m, apf::MeshEntity* elem,
    apf::MeshEntity* face)
{
  apf::Up up;
  m->getUp(face, up);
  if (up.n == 2)
    return up.e[1 - apf::findIn(up.e, 2, elem)];
  if (!m->hasMatching())
    return 0;
  apf::Matches matches;
  m->getMatches(face, matches);
  int self = PCU_Comm_Self();
  for (size_t i = 0; i < matches.getSize(); ++i)
    if (matches[i].peer == self)
      return m->getUpward(matches[i].entity, 0);
  return 0;
}

int* formIENNEIGH(apf::Numbering* ln)
{
  apf::Mesh* m = getMesh(ln);
  int dim = m->getDimension();
  int sideDim = dim - 1;
  int type = getFirstType(m, dim);
  int nsides = apf::Mesh::adjacentCount[type][sideDim];
  size_t nelem = m->count(dim);
  int* ienneigh = new int[nelem * nsides];
  apf::MeshIterator* it = m->begin(dim);
  apf::MeshEntity* e;
  int i = 0;
  while ((e = m->iterate(it))) {
    assert(m->getType(e) == type);
    assert(face_apf2ph[type]);
    apf::Downward sides;
    m->getDownward(e, sideDim, sides);
    for (int j = 0; j < nsides; ++j) {
      apf::MeshEntity* oe = getOtherElem(m, e, sides[j]);
      int oj = face_apf2ph[type][j];
      int oi = oe ? (getNumber(ln, oe, 0, 0) + 1) : 0;
      ienneigh[oj * nelem + i] = oi;
    }
    ++i;
  }
  m->end(it);
  return ienneigh;
}

}
