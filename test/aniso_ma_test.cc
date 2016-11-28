#include "ma.h"
#include <apf.h>
#include <apfNumbering.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apfShape.h>
#include <PCU.h>
#include <cassert>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <algorithm>

#include <spr.h>

#include "mthQR.h"

class AnIso : public ma::AnisotropicFunction
{
public:
    
    /** \brief C-tor for mesh adaptor
     * \param m mesh
     * \param hess_fld recovered Hessian field
     */
    AnIso(apf::Mesh2* m, apf::Field* hess_fld) :
        mesh(m),
        hess(hess_fld)
    { }
    
    /** \brief get prescribed rotation and deformation at a given mesh vertex 
     * \param v mesh vertex
     * \param R rotation matrix
     * \param S deformation vector
     */
    virtual void getValue(apf::MeshEntity* v, apf::Matrix3x3& R, apf::Vector3& S)
    {
        // Get exact Hessian at vertex v
        //apf::Vector3 p = ma::getPosition(mesh,v);
        //double x = p[0], y = p[1], z = p[2];
        //double Pe = 1.5;
        //double uxx = 9*Pe*Pe*std::exp(Pe*(3*x-y));
        //double uxy =-3*Pe*Pe*std::exp(Pe*(3*x-y));
        //double uyy =  Pe*Pe*std::exp(Pe*(3*x-y));
        //double uxz = 0.;
        //double uyz = 0.;
        //double uzz = 0.;
        //apf::Matrix3x3 H(uxx, uxy, uxz,
        //        uxy, uyy, uyz,
        //        uxz, uyz, uzz);

        // Get Hessian tensor at vertex v
        apf::Matrix3x3 H;
        getMatrix(hess, v, 0, H);

        // Decomposition of Hessian: A = Q'*Lambda*Q
        mth::Matrix<double,3,3> A, Q, Lambda;
        A = *reinterpret_cast<mth::Matrix<double,3,3>* >(&H);
        mth::eigenQR(A, Lambda, Q, 100);

        // Rotation Matrix
        R = *reinterpret_cast<ma::Matrix*>(&Q);
        
        // Deformation Matrix
        double epsilon    = 10;
        double dim        = 2;
        double c_d        = 1./2. * std::pow(dim/(dim+1),2);
        double h_min      = 0.001;
        double h_max      = 0.5;
        for (int i=0; i<3; ++i)  
            S[i] = std::min( std::max(epsilon/(1.e-22 + std::abs(Lambda(i,i))), h_min), h_max ); 

    }
private:
    apf::Mesh2* mesh;
    apf::Field* hess;
};


double fun(apf::Vector3 p) {
    double x=p[0], y=p[1], z=p[2];
    double Pe = 1.5;
    return std::exp(Pe*(3*x-y));
    //return std::pow(x,2)-x;
}


int main(int argc, char** argv)
{
    assert(argc==3);
    const char* modelFile = argv[1];
    const char* meshFile = argv[2];
    MPI_Init(&argc,&argv);
    PCU_Comm_Init();
    gmi_register_mesh();
    apf::Mesh2* m = apf::loadMdsMesh(modelFile,meshFile);
    m->verify();
    
    
    // I. Fields
    apf::Field* fld = apf::createLagrangeField(m, "v", apf::SCALAR, 1);
    apf::zeroField(fld); 
    apf::GlobalNumbering*   owned = makeGlobal(apf::numberOwnedNodes(m, "owned", apf::getShape(fld)));
    apf::GlobalNumbering*  shared = makeGlobal(apf::numberOwnedNodes(m, "shared", apf::getShape(fld)));
    apf::synchronize(shared);
    apf::DynamicArray<apf::Node> nodes;
    apf::getNodes(owned, nodes);
    for (size_t i=0; i < nodes.getSize(); ++i) {
        apf::Vector3 p;
        m->getPoint(nodes[i].entity, nodes[i].node, p);
        apf::setScalar(fld, nodes[i].entity, nodes[i].node, fun(p));
    }
    apf::synchronize(fld);
    apf::Field* g = spr::getGradIPField(fld, "g", 1);
    apf::synchronize(g);
    apf::Field* gstar = spr::recoverField(g);
    apf::synchronize(gstar);
    apf::Field* H = spr::getGradIPField(gstar, "H", 1);
    apf::synchronize(H);
    apf::Field* Hstar = spr::recoverField(H);
    apf::synchronize(Hstar);
    apf::destroyField(g);
    apf::destroyField(gstar);
    apf::destroyField(H);

    apf::writeVtkFiles("aniso_before",m);


    // II. Mesh Adaptation    
    AnIso sf(m, Hstar);
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

