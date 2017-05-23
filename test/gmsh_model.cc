#include <sstream>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <map>
#include <set>
#include <vector>

#include <apf.h>
#include <gmi_null.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <PCU.h>

void CreateDmgFromGmsh(const char* filename)
{
    /// Declaring files for read and write
    std::ifstream mshfile(filename);
    std::ofstream mdlfile("model.dmg");

    std::string line;
    /// Skip firsti several lines until we find lines with topological info
    while (std::getline(mshfile, line)) {
        if (line.compare("$Entities") == 0)
            break;
    }    

    /// Exit if there is no topological info
    if (mshfile.eof())
    {
        std::cerr << "Mesh file has no topological information, exiting..." << std::endl;
        exit(-1);
    }

    /// Get the number of topological entities and write it to the model file
    long nents[4];
//    std::getline(mshfile, line);
    for (int i = 0; i < 4; ++i)
        mshfile >> nents[i];
    mdlfile << nents[3] << " " << nents[2] << " " << nents[1] << " " << nents[0] << "\n0 0 0\n0 0 0\n";

std::cout << "Number of geom entities: " << nents[0] << ", " << nents[1] << ", " << nents[2] << ", " << nents[3] << std::endl;

    /// Get vertices, geometrical and mesh correspondance
    long gvtag, nphys, phtag;
    std::set<long> setverts;
    for (long i = 0; i < nents[0]; ++i)
    {
        mshfile >> gvtag;
        mshfile >> nphys;

        if (!nphys)
            continue;

        mshfile >> phtag;
        setverts.insert(gvtag);
    }

    /// Get edges
    long gedge, nverts;
    std::map<long, std::vector<long> > mapedges;
    for (long i = 0; i < nents[1]; ++i)
    {
        std::vector<long> verts(2);
        mshfile >> gedge;
        mshfile >> nverts;
        for (int j = 0; j < 2; ++j)
            mshfile >> verts[j];
        mshfile >> nphys;
        mshfile >> phtag;
        
        if (!nphys)
            continue;
        mapedges[gedge] = verts;
    }

    /// Get faces
    long gface, nedges;
    std::map<long, std::vector<long> > mapfaces;
    for (long i = 0; i < nents[2]; ++i)
    {
        mshfile >> gface;
        mshfile >> nedges;
        std::vector<long> edges(nedges);
        for (int j = 0; j < nedges; ++j)
            mshfile >> edges[j];
        mshfile >> nphys;
        mshfile >> phtag;

        if (!nphys)
            continue;
        mapfaces[gface] = edges;
    }

    /// Get regions
    long gregion, nfaces;
    std::map<long, std::vector<long> > mapregions;
    for (long i = 0; i < nents[3]; ++i)
    {
        mshfile >> gregion;
        mshfile >> nfaces;
        std::vector<long> faces(nfaces);
        for (int j = 0; j < nfaces; ++j)
            mshfile >> faces[j];
        mshfile >> nphys;

        if (!nphys)
            continue;

        mshfile >> phtag;
        mapregions[gregion] = faces;
    }

    /// Reading and writing vertices
    while (std::getline(mshfile, line)) {
        if (line.compare("$Elements") == 0)
        {
            std::getline(mshfile, line);
            break;
        }
    }


    /// Associate geometrical vertices with mesh ones
    long nprocessed = 0;
    long mptag, mvtype, nattr, mvtag;
    std::map<long, long> mapverts;
    std::set<long>::iterator svit;
    while(nprocessed < (long)setverts.size())
    {
        mshfile >> mptag;
        mshfile >> mvtype;
        mshfile >> nattr;
        mshfile >> phtag;
        mshfile >> gvtag;
        mshfile >> mvtag;

        if (setverts.count(gvtag))
        {
            mapverts[mvtag] = gvtag;
            nprocessed++;
        }
    }


    mshfile.close();
    mshfile.open(filename, std::fstream::in);

    while (std::getline(mshfile, line)) {
        if (line.compare("$Nodes") == 0)
            break;
    }

    /// Number of mesh vertices
    long nmvtags;
    mshfile >> nmvtags;

    nprocessed = 0;
    double coords[3];
    std::map<long, long>::iterator mvit;
    for (long i = 0; i < nmvtags; ++i)
    {
        mshfile >> mvtag;
        for (int j = 0; j < 3; ++j)
            mshfile >> coords[j];
        mvit = mapverts.find(mvtag);
        if (mvit != mapverts.end())
        {
            mdlfile << mvit->second << " " << coords[0] << " " << coords[1] << " " << coords[2] << "\n";
            ++nprocessed;
            if (nprocessed == nents[0])
                break;
        }
    }

    if (nprocessed != nents[0])
    {
        std::cerr << "Number of mesh vertices corresponding to geometrical ones (" << nprocessed << ") is smaller than the number of geometrical vertices: " << nents[0] << std::endl;
        exit(-1);
    }

    mshfile.close();

    /// Writing edges, faces, regions
    std::map<long, std::vector<long> >::iterator mentit;
    for (mentit = mapedges.begin(); mentit != mapedges.end(); ++mentit)
    {
        mdlfile << mentit->first << " " << mentit->second[0] << " " << mentit->second[1] << "\n";
    }
    for (mentit = mapfaces.begin(); mentit != mapfaces.end(); ++mentit)
    {
        mdlfile << mentit->first << " 1\n ";
        mdlfile << mentit->second.size() << "\n";
        for (size_t j = 0; j < mentit->second.size(); ++j)
        {
             if (mentit->second[j] < 0)
                 mdlfile << "  " << -mentit->second[j] << " 0\n";
             else
                 mdlfile << "  " << mentit->second[j] << " 1\n";
        }
    }
    for (mentit = mapregions.begin(); mentit != mapregions.end(); ++mentit)
    {
        mdlfile << mentit->first << " 1\n ";
        mdlfile << mentit->second.size() << "\n";
        for (size_t j = 0; j < mentit->second.size(); ++j)
        {
             if (mentit->second[j] < 0)
                 mdlfile << "  " << -mentit->second[j] << " 0\n";
             else
                 mdlfile << "  " << mentit->second[j] << " 1\n";
        }
    }

    mdlfile.close();
}


int main(int argc, char** argv)
{
    /// Initiate communicators
    MPI_Init(&argc,&argv);
    PCU_Comm_Init();
    if ( argc != 3 ) {
        if ( !PCU_Comm_Self() )
            printf("Usage: %s <in .msh> <out .smb>\n", argv[0]);
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }

    /// Create the model from the topological information in the mesh
    CreateDmgFromGmsh(argv[1]);

    std::cout << "Created .dmg geometry, converting the mesh" << std::endl;

    /// Setting up empty models and meshes
    gmi_register_null();
    gmi_register_mesh();

    /// Load the mesh with newly created model
    apf::Mesh2* m = apf::loadMdsFromGmsh(gmi_load("model.dmg"), argv[1]);

    /// Verify mesh
    m->verify();

    /// Write out the mesh
    m->writeNative(argv[2]);

    /// Write out vtk files
    writeVtkFiles("out", m);

    /// Clean up and exit
    m->destroyNative();
    apf::destroyMesh(m);
    PCU_Comm_Free();
    MPI_Finalize();
}

