#include <gmi_mesh.h>
#include <gmi_null.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apfConvert.h>
#include <apf.h>
#include <PCU.h>
#include <cassert>

int main(int argc, char** argv)
{
  assert(argc==3);
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  gmi_register_mesh();
  gmi_register_null();
  int* conn;
  double* coords;
  int nelem;
  int etype;
  int nverts;

  apf::Mesh2* m = apf::loadMdsMesh(argv[1],argv[2]);
  int dim = m->getDimension();

  int nconn = m->count(dim) * 4;

  extractCoords(m, coords, nverts);
  destruct(m, conn, nelem, etype);
  m->destroyNative();
  apf::destroyMesh(m);

  printf("%d: number of connections: %d\n", PCU_Comm_Self(), nconn);
  int psize;
  PCU_Comm_Size(&psize);
  for (int pid = 0; pid < psize; ++pid) {
      if (PCU_Comm_Self() != pid)
          continue;
      printf("%d\n", PCU_Comm_Self());
      for (int i = 0; i < nconn; ++i)
          printf("%d ", conn[i]);
      printf("\n\n");
      PCU_Barrier();
  }


  gmi_model* model = gmi_load(".null");
  m = apf::makeEmptyMdsMesh(model, dim, false);
  apf::GlobalToVert outMap;
  apf::construct(m, conn, nelem, etype, outMap);
  delete [] conn;
  apf::alignMdsRemotes(m);
  apf::deriveMdsModel(m);
  apf::setCoords(m, coords, nverts, outMap);
  delete [] coords;
  outMap.clear();
  m->verify();

  //apf::writeVtkFiles("after", m);

  m->destroyNative();
  apf::destroyMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}

