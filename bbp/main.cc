#include <apf.h>
#include <gmi_null.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <PCU.h>
#include <cstdlib>

int main(int argc, char** argv)
{
    /// Initialization
    MPI_Init(&argc,&argv);
    PCU_Comm_Init();
    if (argc != 3) {
        if (!PCU_Comm_Self())
            printf("Usage: %s <model> <mesh>\n", argv[0]);
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }

    /// Init model and load the (parallel) mesh
    gmi_register_null();
//    gmi_register_mesh();

    apf::Mesh2* m = apf::loadMdsMesh(argv[1], argv[2]);
    m->verify();

    /// Iterate over vertices and get the remote copy information
    apf::MeshIterator* it = m->begin(0);
    apf::MeshEntity* e;
    while ((e = m->iterate(it))) {
        apf::Copies remotes;
        m->getRemotes(e, remotes);

        /// Skip if there are no remote copies
        if (remotes.empty())
            continue;

        /// Iterate through remote copies
        for (apf::Copies::iterator it = remotes.begin(); it != remotes.end(); ++it)
        {
            ;
        }
    }
    m->end(it);

    /// Clean up and finalize
    apf::destroyMesh(m);
    PCU_Comm_Free();
    MPI_Finalize();
}

