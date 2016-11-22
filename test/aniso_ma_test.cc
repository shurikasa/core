#include "ma.h"
#include <apf.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apfShape.h>
#include <PCU.h>
#include <cassert>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <algorithm>

#include "mthQR.h"

class AnIso : public ma::AnisotropicFunction
{
public:
    AnIso(ma::Mesh* m)
    {
        mesh = m;
        average = ma::getAverageEdgeLength(m);
        ma::getBoundingBox(m,lower,upper);
    }
    virtual void getValue(ma::Entity* v, ma::Matrix& R, ma::Vector& H)
    {
        ma::Vector p = ma::getPosition(mesh,v);
        double x = p[0], y = p[1], z = p[2];
        double Pe = 1.5;

        double uxx = 9*Pe*Pe* exp(Pe*(3*x-y));
        double uxy =-3*Pe*Pe* exp(Pe*(3*x-y));
        double uyy =  Pe*Pe* exp(Pe*(3*x-y));
        ma::Matrix r(uxx, uxy, 0.,
                uxy, uyy, 0.,
                0., 0., 1.);

        mth::Matrix<double,3,3> A, Q, Lambda;
        for (int i=0; i<3; ++i) for (int j=0; j<3; ++j) A(i,j) = r[i][j];

        // Decomposition: A = Q'*Lambda*Q
        mth::eigenQR(A, Lambda, Q, 100);

        // Rotation Matrix
        for (int i=0; i<3; ++i) 
            for (int j=0; j<3; ++j) 
                R[i][j] = Q(i,j);
        // Deformation Matrix
        double epsilon    = 0.1;
        double dim        = 2;
        double c_d        = 1./2. * std::pow(dim/(dim+1),2);
        double h_min      = 0.001;
        double h_max      = 0.5;
        for (int i=0; i<3; ++i) 
            H[i] = std::min( std::max(10./(1.e-22 + std::abs(Lambda(i,i))), h_min), h_max ); 

    }
private:
    ma::Mesh* mesh;
    double average;
    ma::Vector lower;
    ma::Vector upper;
};



int main(int argc, char** argv)
{
    assert(argc==3);
    const char* modelFile = argv[1];
    const char* meshFile = argv[2];
    MPI_Init(&argc,&argv);
    PCU_Comm_Init();
    gmi_register_mesh();
    ma::Mesh* m = apf::loadMdsMesh(modelFile,meshFile);
    m->verify();
    apf::writeVtkFiles("aniso_before",m);
    AnIso sf(m);
    ma::Input* in = ma::configure(m, &sf);
    in->shouldRunPreZoltan = true;
    in->shouldRunMidParma = true;
    in->shouldRunPostParma = true;
    in->shouldRefineLayer = true;
    ma::adapt(in);
    m->verify();
    apf::writeVtkFiles("aniso_after",m);
    m->destroyNative();
    apf::destroyMesh(m);
    PCU_Comm_Free();
    MPI_Finalize();
}

