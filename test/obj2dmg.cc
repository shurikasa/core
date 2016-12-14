#include <sstream>
#include <fstream>
#include <vector>
#include <map>
#include <string>

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
            printf("Usage: %s <in .obj> <out .dmg>\n", argv[0]);
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }

    /// Declaring files for read and write
    std::ifstream ifile(argv[1]);
    std::ofstream ofile(argv[2]);

    std::vector< std::vector<double> > verts;
    std::map< std::vector<long>, long > edges;
    std::vector< std::vector<long> > faces;

    std::string line;
    /// Skip firsti several lines until we find lines with vertices
    while (std::getline(ifile, line)) {
        if (line.substr(0, 2).compare("v ") == 0)
            break; 
    }

    /// Getting coordinates
    double* coords = new double[3];
    do {
        if (line.substr(0, 2).compare("v ") == 0) {
            std::vector<double> xyz(3);
            sscanf(line.c_str(), "v %lf %lf %lf\n", coords, coords+1, coords+2);
            for (int i = 0; i < 3; ++i)
                xyz[i] = coords[i];
            verts.push_back(xyz);
        }
        else
            break;
    }
    while (std::getline(ifile, line));

std::cout << "number of verts: " << verts.size() << std::endl;

    /// Skip firsti several lines until we find lines with faces
    while (std::getline(ifile, line)) {
        if (line.substr(0, 2).compare("f ") == 0) 
            break;
    }


    /// Getting faces by vertex indices
    long* fvind = new long[3];
    do {
        if (line.substr(0, 2).compare("f ") == 0) {
            std::vector<long> ind(3);
            sscanf(line.c_str(), "f %ld %ld %ld\n", fvind, fvind+1, fvind+2);
            for (int i = 0; i < 3; ++i)
                ind[i] = fvind[i]-1;
            faces.push_back(ind);
        }
        else
            break;
    }
    while (std::getline(ifile, line));

    ifile.close();

std::cout << "number of faces: " << faces.size() << std::endl;

    /// Construct edges from face to vertex information
    size_t tag_counter = 0;
    for (size_t i = 0; i < faces.size(); ++i) {
        for (size_t j = 0; j < 3; ++j) {
            std::vector<long> ind(2);
            if (faces[i][j] > faces[i][(j+1)%3]) {
                ind[0] = faces[i][(j+1)%3];
                ind[1] = faces[i][j];
            }
            else {
                ind[0] = faces[i][j];
                ind[1] = faces[i][(j+1)%3];
            }
            if (edges.count(ind))
                continue;

            edges[ind] = tag_counter + verts.size();
            ++tag_counter;
        }
    }

    /// Writing dmg file
    ofile << "1 " << faces.size() << " " << edges.size() << " " << verts.size() << "\n0 0 0\n0 0 0\n";
    tag_counter = 0;
    /// Writing vertices
    for (size_t i = 0; i < verts.size(); ++i, ++tag_counter)
        ofile << tag_counter << " " << verts[i][0] << " " << verts[i][1] << " " << verts[i][2] << "\n";

    /// Writing edges
    for (std::map< std::vector<long>, long >::iterator mapit = edges.begin(); mapit != edges.end(); ++mapit, ++tag_counter)
        ofile << mapit->second << " " << mapit->first[0] << " " << mapit->first[1] << "\n";

    /// Writing faces
    size_t face_tag_counter = tag_counter;
    for (size_t i = 0; i < faces.size(); ++i, ++tag_counter) {
        ofile << tag_counter << " 1\n 3\n";
        for (size_t j = 0; j < 3; ++j) {
            std::vector<long> ind(2);
            int dir;
            if (faces[i][j] > faces[i][(j+1)%3]) {
                ind[0] = faces[i][(j+1)%3];
                ind[1] = faces[i][j];
                dir = 0;
            }
            else {
                ind[0] = faces[i][j];
                ind[1] = faces[i][(j+1)%3];
                dir = 1;
            }
            ofile << "  " << edges[ind]  << " " << dir << "\n";
        }
    }

    /// Write the volume tag and all face tags comprising it
    ofile << tag_counter << " 1\n " << faces.size() << "\n";
    for (size_t i = 0; i < faces.size(); ++i, ++face_tag_counter)
        ofile << "  " << face_tag_counter << " 1\n";

    ofile.close();

    PCU_Comm_Free();
    MPI_Finalize();
}

