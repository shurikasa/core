#include <apf.h>
#include <gmi_null.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <PCU.h>
#include <apfNumbering.h>
#include <cstdlib>

typedef std::set<apf::MeshEntity*> EntitySet;
typedef std::set<apf::ModelEntity*> ModelEntitySet;
std::vector<ModelEntitySet> neurons_model_ents;


/// Get faces and vertices on the cube geometrical boundary
void getFacesVertsOnCubeBdry(apf::Mesh2* m, EntitySet& faces_on_cube_bdry, EntitySet& verts_on_cube_bdry, ModelEntitySet& model_edges_on_cube) {
    apf::MeshEntity* e;
    apf::MeshIterator* it;

    /// Iterate over faces
    it = m->begin(2);
    while ((e = m->iterate(it))) {
        apf::ModelEntity* me = m->toModel(e);
        /// If face is classified on a model boundary
        if (m->getModelType(me) == 2) {
            /// Add face to the entity set and its vertices also
            faces_on_cube_bdry.insert(e);

            apf::Downward de;
            int nde = m->getDownward(e, 0, de);
            for (int i=0; i < nde; ++i)
                verts_on_cube_bdry.insert(de[i]);

            /// Also remember all the model edges of cube model faces
            /// to avoid adding them later when dealing with neuron
            nde = m->getDownward(e, 1, de);
            for (int i=0; i < nde; ++i){
                me = m->toModel(de[i]);
                if (m->getModelType(me) == 1)
                    model_edges_on_cube.insert(me);
            }
        }
    }
    m->end(it);
}


void getEdgesVertsOnNeuronBdry(apf::Mesh2* m, EntitySet& edges_on_neuron_bdry, EntitySet& verts_on_neuron_bdry, ModelEntitySet& model_edges_on_cube) {
    apf::MeshEntity* e;
    apf::MeshIterator* it;

    /// Iterate over edges
    it = m->begin(1);
    while ((e = m->iterate(it))) {
        apf::ModelEntity* me = m->toModel(e);
        /// If edge is classified on a neuron-only  model boundary
        if (m->getModelType(me) == 1 && model_edges_on_cube.count(me) == 0) {
            /// Add edge to the entity set and its vertices also
            edges_on_neuron_bdry.insert(e);

            apf::Downward de;
            int nde = m->getDownward(e, 0, de);
            for (int i=0; i < nde; ++i)
                verts_on_neuron_bdry.insert(de[i]);
        }
    }
    m->end(it);
}


/// Get internal entities
void getEntsInternal(apf::Mesh2* m, EntitySet& ents_internal, int dim) {
    apf::MeshEntity* e;
    apf::MeshIterator* it;

    /// Iterate over vertices
    it = m->begin(dim);
    while ((e = m->iterate(it))) {
        apf::ModelEntity* me = m->toModel(e);
        /// Add interior vertices to the entity set
        if (m->getModelType(me) == 3)
            ents_internal.insert(e);
    }
    m->end(it);    
}


/// Get entities on the partition boundary based on entities from a given entity set
void getEntsOnPartBdry(apf::Mesh2* m, EntitySet& ents, EntitySet& ents_on_part_bdry) {
    if (PCU_Comm_Peers() == 1)
        return;

    apf::MeshEntity* e;
    EntitySet::iterator it;

    /// Iterate through entities on the model boundary from the corresponding entity set
    for (it = ents.begin(); it != ents.end(); ++it) {
        e = *it;
        /// If the face is on the part boundary (unlikely unless something is ), add it to the entity set
        if (m->isShared(e))
           ents_on_part_bdry.insert(e);
    }
}


/// Check if the neuron field corresponds correctly to the designated neuron
/// Ideally it should be neuron's id, where for each neuron id there's a list of geometrical entities.
void verifyNeuronField(apf::Mesh2*m, const char* fname, ModelEntitySet& neuron_model_ents) {
    apf::MeshEntity* v;
    apf::MeshIterator* it;
    apf::Field* f = m->findField(fname);

    /// Iterate over vertices
    it = m->begin(0);
    while ((v = m->iterate(it))) {
       /// Check the field corresponds to its corresponding neuron group
       apf::getScalar(f, v, 0);
       
    }
    m->end(it);

    neuron_model_ents.size();
}



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
    gmi_register_mesh();

    apf::Mesh2* m = apf::loadMdsMesh(gmi_load(argv[1]), argv[2]);
//    m->verify();

    /// Entity sets with requested types of entities
    EntitySet regions_internal, verts_internal, faces_on_cube_bdry, verts_on_cube_bdry, edges_on_neuron_bdry, verts_on_neuron_bdry;
    EntitySet verts_internal_part_bdry, verts_on_cube_bdry_part_bdry, edges_on_neuron_bdry_part_bdry, verts_on_neuron_bdry_part_bdry;
    ModelEntitySet model_edges_on_cube;

    /// Get faces and vertices on the cube geometrical boundary
    getFacesVertsOnCubeBdry(m, faces_on_cube_bdry, verts_on_cube_bdry, model_edges_on_cube);

    /// Get vertices on the cube geometrical boundary and partition boundary
    getEntsOnPartBdry(m, verts_on_cube_bdry, verts_on_cube_bdry_part_bdry);


    /// Get edges and vertices on the neuron geometrical boundary
    getEdgesVertsOnNeuronBdry(m, edges_on_neuron_bdry, verts_on_neuron_bdry, model_edges_on_cube);

    /// Get edges and vertices on the neuron geometrical boundary and partition boundary
    getEntsOnPartBdry(m, edges_on_neuron_bdry, edges_on_neuron_bdry_part_bdry);
    getEntsOnPartBdry(m, verts_on_neuron_bdry, verts_on_neuron_bdry_part_bdry);


    /// Get internal vertices
    getEntsInternal(m, verts_internal, 0);

    /// Get internal vertices on partition boundary
    getEntsOnPartBdry(m, verts_internal, verts_internal_part_bdry);


    /// Get regions
    getEntsInternal(m, regions_internal, 3);


    /// Setting the neuron field on vertices classified on neuron edges
    const char* field_name = "neuron";
    apf::Field* f = apf::createLagrangeField(m, field_name, apf::SCALAR, 1);
    apf::zeroField(f);
    apf::MeshEntity* e;
    EntitySet::iterator it;

    /// Iterate through entities on the model boundary from the corresponding entity set
    for (it = verts_on_neuron_bdry.begin(); it != verts_on_neuron_bdry.end(); ++it) {
        e = *it;
        apf::setScalar(f, e, 0, 1);
    }
    apf::synchronize(f);

    /// Check if the neuron field corresponds correctly to the designated neuron
//    verifyNeuronField(m, field_name);

    if (PCU_Comm_Self() == 0)
        printf("Mesh statistics for %d parts:\n", PCU_Comm_Peers());

    printf("%d: %d regions\n", PCU_Comm_Self(), (int)regions_internal.size());
    printf("%d: %d internal vertices\n", PCU_Comm_Self(), (int)verts_internal.size());
    if (PCU_Comm_Peers() > 1)
        printf("%d: %d internal vertices on part boundary\n", PCU_Comm_Self(), (int)verts_internal_part_bdry.size());
    printf("%d: %d faces classified on cube faces\n", PCU_Comm_Self(), (int)faces_on_cube_bdry.size());
    printf("%d: %d vertices classified on cube model entities\n", PCU_Comm_Self(), (int)verts_on_cube_bdry.size());
    if (PCU_Comm_Peers() > 1)
        printf("%d: %d vertices classified on cube model entities on part boundary\n", PCU_Comm_Self(), (int)verts_on_cube_bdry_part_bdry.size());
    printf("%d: %d edges classified on neuron edges\n", PCU_Comm_Self(), (int)edges_on_neuron_bdry.size());
    if (PCU_Comm_Peers() > 1)
        printf("%d: %d edges classified on neuron edges and on part boundary\n", PCU_Comm_Self(), (int)edges_on_neuron_bdry_part_bdry.size());
    printf("%d: %d vertices classified on neuron model entities\n", PCU_Comm_Self(), (int)verts_on_neuron_bdry.size());
    if (PCU_Comm_Peers() > 1)
        printf("%d: %d vertices classified on neuron model entities on part boundary\n", PCU_Comm_Self(), (int)verts_on_neuron_bdry_part_bdry.size());

    /// Writing the vtk format
    writeVtkFiles("out", m);

    /// Clean up and finalize
    apf::destroyField(f);
    apf::destroyMesh(m);
    PCU_Comm_Free();
    MPI_Finalize();
}

