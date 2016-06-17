#include <PCU.h>
#include "phLinks.h"
#include "phAdjacent.h"
#include <apf.h>
#include <cassert>
#include <cstdlib>

namespace ph {

bool LinkKey::operator<(LinkKey const& other) const
{
  if (send != other.send)
    return send;
  return peer < other.peer;
}

typedef std::map<apf::Parts,int> PtnMdl;
typedef std::map<int,int> Tasks;
typedef std::map<int,std::set<int> > PeerTasks;

void printPtnMdlEnt(const char* key, const apf::Parts& pme) {
  PCU_Debug_Print("%s pme ", key);
  for (apf::Parts::iterator pid = pme.begin(); pid != pme.end(); pid++) {
    PCU_Debug_Print("%d ", *pid);
  }
}

void printPeerTask(const char* key, std::set<int>& tasks) {
  PCU_Debug_Print("%s sends-to-peer ", key);
  for (apf::Parts::iterator pid = tasks.begin(); pid != tasks.end(); pid++)
    PCU_Debug_Print("%d ", *pid);
}

void printPtnMdl(const char* key, PtnMdl& pm) {
  APF_ITERATE(PtnMdl,pm,it) {
    printPtnMdlEnt(key,it->first);
    PCU_Debug_Print(" owner %d\n", it->second);
  }
}

/* returns the current partition model based on NormalSharing */
void getPtnMdl(apf::Mesh* m, PtnMdl& pm) {
  apf::MeshIterator* it = m->begin(0);
  apf::MeshEntity* e;
  while ((e = m->iterate(it))) {
    if( m->isShared(e) ) {
      apf::Parts res;
      m->getResidence(e,res);
      pm[res] = m->getOwner(e);
    }
  }
}

void initTasks(PtnMdl& pm, Tasks& in, Tasks& out, Tasks& send, Tasks& recv) {
  send.clear();
  recv.clear();
  const int self = PCU_Comm_Self();
  APF_ITERATE(PtnMdl,pm,it)
    for (apf::Parts::iterator pid = it->first.begin(); pid != it->first.end(); pid++)
      in[*pid] = out[*pid] = 0; //need to know that the peer exists
  //count the outbound comm tasks
  std::set<int> peers;
  APF_ITERATE(PtnMdl,pm,it) {
    if( it->second == self )
      for (apf::Parts::iterator pid = it->first.begin(); pid != it->first.end(); pid++)
        if( *pid != self ) {
          peers.insert(*pid);
          send[*pid]++;
        }
  }
  out[self] = peers.size();
  //count the inbound comm tasks
  peers.clear();
  APF_ITERATE(PtnMdl,pm,it) {
    if( it->second != self ) {
      peers.insert(it->second);
      recv[it->second]++;
    }
  }
  in[self] = peers.size();
  PCU_Debug_Print("in %d out %d\n",in[self],out[self]);
}

void printTasks(Tasks& tasks, const char* key="") {
  APF_ITERATE(Tasks,tasks,it)
    PCU_Debug_Print("%s peer %d over %d links\n", key, it->first, it->second);
}

/* exchange the number of owned partition model entities with neighbors */
void exchangeTasks(Tasks& tasks) {
  const int self = PCU_Comm_Self();
  PCU_Comm_Begin();
  //send the number of owned partition model ents
  APF_ITERATE(Tasks,tasks,it) {
    PCU_COMM_PACK(it->first, tasks[self]);
  }
  PCU_Comm_Send();
  //recv the neighbors owned partition model ents counts
  while (PCU_Comm_Listen()) {
    int peer = PCU_Comm_Sender();
    int ownCnt = 0;
    PCU_COMM_UNPACK(ownCnt);
    tasks[peer] = ownCnt;
  }
}

void printPeerTasks(PeerTasks& pt, const char* key="") {
  APF_ITERATE(PeerTasks,pt,it) {
    PCU_Debug_Print("%s peer %d sends-to-peer", key, it->first);
    std::set<int>& task = it->second;
    for (std::set<int>::iterator pid = task.begin(); pid != task.end(); pid++)
      PCU_Debug_Print(" %d ", *pid);
    PCU_Debug_Print("\n");
  }
}

/* exchange the list of outbound comm tasks */
void exchangeSendTasks(Tasks& out, Tasks& send, PeerTasks& peerOut) {
  peerOut.clear();
  int* peersSentTo = new int[send.size()];
  int i = 0;
  APF_ITERATE(Tasks,send,it)
    peersSentTo[i++] = it->first;
  PCU_Comm_Begin();
  //send the list of outbound comm tasks
  APF_ITERATE(Tasks,out,it) {
    int peer = it->first;
    if(peer == PCU_Comm_Self() ) continue;
    int size = send.size();
    PCU_COMM_PACK(peer, size);
    PCU_Comm_Pack(peer, peersSentTo, sizeof(int)*send.size());
  }
  PCU_Comm_Send();
  //recv the neighbors outbound comm task list
  while (PCU_Comm_Listen()) {
    int peer = PCU_Comm_Sender();
    int size = 0;
    PCU_COMM_UNPACK(size);
    int* outcomm = new int[size];
    PCU_Comm_Unpack(outcomm, sizeof(int)*size);
    for(int i=0; i<size; i++)
      peerOut[peer].insert(outcomm[i]);
    delete [] outcomm;
  }
  delete [] peersSentTo;
}


/* send the new owner of the partition model entity to all parts
   bounded by the entity (the resident set 'res') */
void packRes(const apf::Parts& res, int newowner) {
  int* resArray = new int[res.size()];
  int i = 0;
  for (apf::Parts::iterator pid = res.begin(); pid != res.end(); pid++)
    resArray[i++] = *pid;
  for (apf::Parts::iterator pid = res.begin(); pid != res.end(); pid++) {
    if( *pid != PCU_Comm_Self() ) {
      int peer = *pid;
      int sz = res.size();
      printPtnMdlEnt("packres sending",res);
      PCU_Debug_Print(" to peer %d newowner %d\n",peer,newowner);
      PCU_COMM_PACK(peer, sz);
      PCU_Comm_Pack(peer, resArray, sizeof(int)*res.size());
      PCU_COMM_PACK(peer, newowner);
    }
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

/* check to see if there exists a peer with outbound tasks s.t. assigning the
 * peer ownership of the ptn mdl ent does not create additional outbound tasks*/
int poorestNeighbor(const apf::Parts& res, PeerTasks& pt) {
  const int self = PCU_Comm_Self();
  int newowner = -1;
  //loop over possible new owners ('candidate') for the ptn mdl ent
  for (apf::Parts::iterator candidate = res.begin();
      candidate != res.end();
      candidate++) {
    size_t found = 0;
    //the new owner can't be the current owner
    if( *candidate == self ) continue;
    if( ! pt.count(*candidate) ) continue;
    PCU_Debug_Print("candidate %d\n", *candidate);
    printPtnMdlEnt("checking", res);
    printPeerTask("against", pt.at(*candidate));
    PCU_Debug_Print("\n");
    //check if the candidate has the needed outbound tasks
    for (apf::Parts::iterator pid = res.begin();
        pid != res.end();
        pid++) {
      found += pt.at(*candidate).count(*pid);
    }
    printPtnMdlEnt("checking", res);
    printPtnMdlEnt("against", pt.at(*candidate));
    PCU_Debug_Print("candidate %d found %lu of %lu\n", *candidate, found, res.size()-1);
    if( found == res.size()-1 ) {
      newowner = *candidate;
      break;
    }
  }
  if( newowner != -1 )
    PCU_Debug_Print("newowner %d\n", newowner);
  return newowner;
}

void getSortedSendPeers(Tasks& send, int* sorted) {
  typedef std::multimap<int,int> mmii;
  mmii sendSort;
  APF_ITERATE(Tasks,send,s)
    sendSort.insert(std::pair<int,int>(s->second,s->first));
  int i=0;
  PCU_Debug_Print("sorted sends");
  APF_ITERATE(mmii,sendSort,s) {
    sorted[i++] = s->second;
    PCU_Debug_Print(" %d,%d ", s->second, s->first);
  }
  PCU_Debug_Print("\n");
}

/* returns the partition model with balanced ownership */
void balanceOwners(PtnMdl& pm) {
  double t0 = PCU_Time();
  const int self = PCU_Comm_Self();
  printPtnMdl("init",pm);
  Tasks in, out, send, recv;
  initTasks(pm,in,out,send,recv);
  printTasks(send,"sending to");
  exchangeTasks(in);
  exchangeTasks(out);
  PeerTasks pt;
  exchangeSendTasks(out,send,pt);
  printPeerTasks(pt);
  int tasks = in[self]+out[self];
  int totTasks = PCU_Add_Int(tasks);
  int maxTasks = PCU_Max_Int(tasks);
  int maxIn = PCU_Max_Int(in[self]);
  int maxOut = PCU_Max_Int(out[self]);
  const int peers = in.size()-1;
  const int maxPeers = PCU_Max_Int(peers);
  const float idealTasks = maxPeers/2;
  int iter = 0;
  const int maxIter = 5; //maxPeers*2;
  if( !PCU_Comm_Self() )
    fprintf(stderr, "maxIn %d maxOut %d idealTasks %.3f maxIter %d\n",
        maxIn, maxOut, idealTasks, maxIter);
  while( (maxIn > idealTasks || maxOut > idealTasks) && iter++ < maxIter ) {
    PCU_Debug_Print("startiter %d\n", iter);
    if( !PCU_Comm_Self() )
      fprintf(stderr, "totTasks %d maxTasks %d maxIn %d maxOut %d iter %d\n",
          totTasks, maxTasks, maxIn, maxOut, iter);
    PCU_Comm_Begin();
    if( out.at(self) == maxOut ) {
      int* sortedSendPeers = new int[send.size()];
      getSortedSendPeers(send,sortedSendPeers);
      for(size_t i=0; i<send.size(); i++) {
        int minSendPeer = sortedSendPeers[i];
        int minSends = send[sortedSendPeers[i]];
        PCU_Debug_Print("out %d minSendPeer %d minSends %d\n", out.at(self), minSendPeer, minSends);
        int sent = 0;
        PtnMdl toSend;
        APF_ITERATE(PtnMdl,pm,it) {
          if( it->second == self && it->first.count(minSendPeer) ) {
            int newowner = poorestNeighbor(it->first,pt);
            if( newowner != -1 ) {
              assert(newowner != -1);
              toSend[it->first] = newowner;
              sent++;
            }
          }
        }
        PCU_Debug_Print("minSends %d sent %d\n", minSends, sent);
        if( sent == minSends ) {
          APF_ITERATE(PtnMdl,toSend,it) {
            printPtnMdlEnt("setting", it->first);
            PCU_Debug_Print(" newowner %d\n", it->second);
            pm[it->first] = it->second; //change owner
            packRes(it->first,it->second);
          }
          break;
        }
      }
      delete [] sortedSendPeers;
    }
    PCU_Comm_Send();
    //recv the owner changes
    while (PCU_Comm_Listen()) {
      int peer = PCU_Comm_Sender();
      while (!PCU_Comm_Unpacked()) {
        apf::Parts res;
        int newowner = -1;
        unpackRes(res,newowner);
        printPtnMdlEnt("getting", res);
        PCU_Debug_Print(" newowner %d from peer %d\n", newowner, peer);
        assert(newowner != -1);
        assert(res.size() > 1);
        assert(pm.count(res));
        pm[res] = newowner;
      }
    }
    printPtnMdl("update",pm);
    initTasks(pm,in,out,send,recv);
    printTasks(send,"sending to");
    exchangeTasks(in);
    exchangeTasks(out);
    exchangeSendTasks(out,send,pt);
    printPeerTasks(pt);
    maxIn = PCU_Max_Int(in[self]);
    maxOut = PCU_Max_Int(out[self]);
    tasks = in[self]+out[self];
    totTasks = PCU_Add_Int(tasks);
    maxTasks = PCU_Max_Int(tasks);
    PCU_Debug_Print("enditer %d\n", iter);
  }
  if( !PCU_Comm_Self() )
    fprintf(stderr, "totTasks %d maxTasks %d maxIn %d maxOut %d iter %d time %.3f\n",
        totTasks, maxTasks, maxIn, maxOut, iter, PCU_Time()-t0);
}

struct BalancedSharing : public apf::Sharing
{
  BalancedSharing(apf::Mesh* m) {
    mesh = m;
    helper = apf::getSharing(m);
    PCU_Debug_Open();
    getPtnMdl(mesh,pm);
    balanceOwners(pm);
  }
  ~BalancedSharing() {
    delete helper;
  }
  bool isOwned(apf::MeshEntity* e) {
    apf::Parts res;
    mesh->getResidence(e,res);
    if( res.size() == 1 ) {
      return true;
    } else {
      return (pm[res] == PCU_Comm_Self());
    }
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
  BalancedSharing shr(m);
  //PhastaSharing shr(m);
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
