#include <apf.h>
#include <gmi_null.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <PCU.h>
#include <cstdlib>

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  if ( argc != 3 ) {
    if ( !PCU_Comm_Self() )
      printf("Usage: %s <in path> <out .smb>\n", argv[0]);
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  gmi_register_null();
  gmi_register_mesh();
  apf::Mesh2* m = apf::loadMdsFromVtk(gmi_load(".null"), argv[1]);
  m->writeNative(argv[2]);
  writeVtkFiles("out_new", m);
  m->destroyNative();
  apf::destroyMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}

