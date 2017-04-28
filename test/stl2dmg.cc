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
    if ( argc != 4 ) {
        if ( !PCU_Comm_Self() )
            printf("Usage: %s <in .obj> <out .dmg> <out .geo>\n", argv[0]);
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }

    /// Declaring files for read and write
    std::ifstream ifile(argv[1]);
    std::ofstream odfile(argv[2]);
    std::ofstream ogfile(argv[3]);

    std::map< std::vector<double>, long > verts;
    std::map< std::vector<long>, long > edges;
    std::vector< std::vector<long> > faces;

    std::string line;
    /// Skip first several lines until we find lines with vertices
    std::getline(ifile, line);
    double* coords = new double[3];
    long vcounter = 0;
    while (std::getline(ifile, line) ) {
        if (line.find("endsolid") != std::string::npos)
            break;

        std::getline(ifile, line);

        /// Go through vertices and assign an index or check if one already exists
        /// Assign faces vertex indices
        std::vector<long> ind(3);
        for (size_t i = 0; i < 3; ++i) {
            std::getline(ifile, line);
            /// Read vertex coordinates
            sscanf(line.c_str(), "   vertex %lf %lf %lf\n", coords, coords+1, coords+2);
            std::vector<double> xyz(3);
            for (size_t j = 0; j < 3; ++j)
                xyz[j] = coords[j];

            /// Identify if the vertex already exists
            std::map< std::vector<double>, long >::iterator mapit = verts.find(xyz);
            if (mapit != verts.end()) {
                ind[i] = mapit->second;
            }
            else {
                verts[xyz] = vcounter;
                ind[i] = vcounter;
                ++vcounter;
            }
        }
        faces.push_back(ind);

        std::getline(ifile, line);
        std::getline(ifile, line);
    }

std::cout << "number of verts: " << verts.size() << std::endl;
std::cout << "number of faces: " << faces.size() << std::endl;

    ifile.close();

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
    odfile << "1 " << faces.size() << " " << edges.size() << " " << verts.size() << "\n0 0 0\n0 0 0\n";
    tag_counter = 0;
    /// Writing vertices
    for (std::map< std::vector<double>, long >::iterator mapit = verts.begin(); mapit != verts.end(); ++mapit, ++tag_counter) {
        odfile << mapit->second << " " << mapit->first[0] << " " << mapit->first[1] << " " << mapit->first[2] << "\n";
        ogfile << "Point(" << mapit->second << ") = {" << mapit->first[0] << "," << mapit->first[1] << "," << mapit->first[2] << "};\n";
        ogfile << "Physical Point(" << mapit->second << ") = {" << mapit->second << "};\n";
    }

    /// Writing edges
    for (std::map< std::vector<long>, long >::iterator mapit = edges.begin(); mapit != edges.end(); ++mapit, ++tag_counter) {
        odfile << mapit->second << " " << mapit->first[0] << " " << mapit->first[1] << "\n";
        ogfile << "Line(" << mapit->second << ") = {" << mapit->first[0] << "," << mapit->first[1] << "};\n";
        ogfile << "Physical Line(" << mapit->second << ") = {" << mapit->second << "};\n";
    }

    /// Writing faces
    size_t face_tag_counter = tag_counter;
    for (size_t i = 0; i < faces.size(); ++i, ++tag_counter) {
        odfile << tag_counter << " 1\n 3\n";
        ogfile << "Line Loop(" << tag_counter << ") = {";
        for (size_t j = 0; j < 3; ++j) {
            std::vector<long> ind(2);
            /// Identify direction based on vertex sequence
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
            odfile << "  " << edges[ind]  << " " << dir << "\n";
            if (dir)
                ogfile << edges[ind];
            else
                ogfile << -edges[ind];
            if (j < 2)
                ogfile << ",";
        }
        ogfile << "};\nRuled Surface(" << tag_counter << ") = {" << tag_counter << "};\n";
        ogfile << "Physical Surface(" << tag_counter << ") = {" << tag_counter << "};\n";
    }

    /// Write the volume tag and all face tags comprising it
    odfile << tag_counter << " 1\n " << faces.size() << "\n";
    ogfile << "Surface Loop(" << tag_counter << ") = {";
    for (size_t i = 0; i < faces.size(); ++i, ++face_tag_counter) {
        odfile << "  " << face_tag_counter << " 1\n";
        ogfile << face_tag_counter;
        if (i < faces.size() - 1)
            ogfile << ",";
    }
    ogfile << "};\nVolume(" << tag_counter << ") = {" << tag_counter << "};\n";
    ogfile << "Physical Volume(" << tag_counter << ") = {" << tag_counter << "};\n";

    odfile.close();
    ogfile.close();

    PCU_Comm_Free();
    MPI_Finalize();
}

