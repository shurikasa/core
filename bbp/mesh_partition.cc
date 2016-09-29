#include <apf.h>
#include <pumi.h>
#include <gmi_null.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <PCU.h>
#include <cstdlib>

int main(int argc, char** argv)
{
    int num_in_part = argv[3];

    /// Initialization
    MPI_Init(&argc,&argv);
    PCU_Comm_Init();
    if ( argc != 5 ) {
        if ( !PCU_Comm_Self() )
            printf("Usage: %s <model> <in .smb> <in num_parts_per_partition> <out .smb>\n", argv[0]);
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }

    /// Init model and load the mesh on the specified number of parts
    gmi_register_null();
    pGeom g = gModel(gmi_load(argv[1]));
    pumi_mesh_load(g, argv[2], num_in_part, "mds");
    m->verify();

    /// Write the mesh
    m->writeNative(argv[4]);

    /// Clean up and exit
    m->destroyNative();    
    apf::destroyMesh(m);
    PCU_Comm_Free();
    MPI_Finalize();
}

