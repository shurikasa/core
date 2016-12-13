#include <sstream>
#include <fstream>
#include <vector>
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
                xyz.push_back(coords[i]);
            verts.push_back(xyz);
        }
        else
            break;
    }
    while (std::getline(ifile, line));

std::cout << "number of verts: " << verts.size() << std::endl;

    /// Skip firsti several lines until we find lines with faces
    while (std::getline(ifile, line)) {
        std::getline(ifile, line);
        if (line.substr(0, 2).compare("f ") == 0) 
            break;
    }


    /// Getting coordinates
    long*  = new double[3];
    do {
        if (line.substr(0, 2).compare("v ") == 0) {
            std::vector<double> xyz(3);
            sscanf(line.c_str(), "v %lf %lf %lf\n", coords, coords+1, coords+2);
            for (int i = 0; i < 3; ++i)
                xyz.push_back(coords[i]);
            verts.push_back(xyz);
        }
        else
            break;
    }
    while (std::getline(ifile, line));


    PCU_Comm_Free();
    MPI_Finalize();
}

