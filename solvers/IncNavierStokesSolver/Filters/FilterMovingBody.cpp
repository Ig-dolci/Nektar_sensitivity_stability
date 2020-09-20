///////////////////////////////////////////////////////////////////////////////
//
// File FilterMovingBody.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: Output moving body forces during time-stepping.
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <LibUtilities/BasicUtils/ParseUtils.h>
#include <iomanip>
#include <LocalRegions/Expansion1D.h>
#include <LocalRegions/Expansion2D.h>
#include <LocalRegions/Expansion3D.h>
#include <IncNavierStokesSolver/Filters/FilterMovingBody.h>

using namespace std;

namespace Nektar
{

std::string FilterMovingBody::className = SolverUtils::GetFilterFactory().
        RegisterCreatorFunction("MovingBody",
                                FilterMovingBody::create,
                                "Moving Body Filter");

NekDouble bf_drag_t, bf_lift_t, bf_drag_p, bf_lift_p;
NekDouble adj_drag_t, adj_lift_t, adj_drag_p, adj_lift_p;
NekDouble bfgrad_drag_t, bfgrad_lift_t, bfgrad_drag_p, bfgrad_lift_p;
int counter = 0.;
/**
 *
 */
FilterMovingBody::FilterMovingBody(
        const LibUtilities::SessionReaderSharedPtr         &pSession,
        const std::weak_ptr<SolverUtils::EquationSystem> &pEquation,
        const ParamMap &pParams)
    : Filter(pSession, pEquation)
{
    // OutputFile
    auto it = pParams.find("OutputFile");
    if (it == pParams.end())
    {
        m_outputFile_fce = pSession->GetSessionName();
        m_outputFile_mot = pSession->GetSessionName();
    }
    else
    {
        ASSERTL0(it->second.length() > 0, "Missing parameter 'OutputFile'.");

        m_outputFile_fce = it->second;
        m_outputFile_mot = it->second;
    }
    if (!(m_outputFile_fce.length() >= 4 &&
          m_outputFile_fce.substr(m_outputFile_fce.length() - 4) == ".fce"))
    {
        m_outputFile_fce += ".fce";
    }
    if (!(m_outputFile_mot.length() >= 4 &&
          m_outputFile_mot.substr(m_outputFile_mot.length() - 4) == ".mot"))
    {
        m_outputFile_mot += ".mot";
    }

    // OutputFrequency
    it = pParams.find("OutputFrequency");
    if (it == pParams.end())
    {
        m_outputFrequency = 1;
    }
    else
    {
        LibUtilities::Equation equ(m_session->GetInterpreter(), it->second);
        m_outputFrequency = round(equ.Evaluate());
    }

    // pSession->MatchSolverInfo("Homogeneous", "1D", m_isHomogeneous1D, false);
    // ASSERTL0(m_isHomogeneous1D, "Moving Body implemented just for 3D "
    //                             "Homogeneous 1D discetisations.");

    // Boundary (to calculate the forces)
    it = pParams.find("Boundary");
    ASSERTL0(it != pParams.end(),     "Missing parameter 'Boundary'.");
    ASSERTL0(it->second.length() > 0, "Missing parameter 'Boundary'.");
    m_BoundaryString = it->second;
    bf_drag_t = 0.;
    bf_lift_p = 0.;
    bf_drag_p = 0.;
    bf_lift_t = 0.;
    bfgrad_drag_t = 0.;
    bfgrad_lift_p = 0.;
    bfgrad_drag_p = 0.;
    bfgrad_lift_t = 0.;
    
       
}


/**
 *
 */
FilterMovingBody::~FilterMovingBody()
{

}


/**
 *
 */
void FilterMovingBody::v_Initialise(
    const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    m_index_f = 0;
    m_index_m = 0;
    m_outputStream =  Array<OneD, std::ofstream>(2);
    // Parse the boundary regions into a list.
    std::string::size_type FirstInd = m_BoundaryString.find_first_of('[') + 1;
    std::string::size_type LastInd  = m_BoundaryString.find_last_of(']') - 1;

    ASSERTL0(FirstInd <= LastInd,
            (std::string("Error reading boundary region definition:") +
             m_BoundaryString).c_str());

    std::string IndString = m_BoundaryString.substr(FirstInd,
                                                    LastInd - FirstInd + 1);

    bool parseGood = ParseUtils::GenerateSeqVector(IndString,
                                                   m_boundaryRegionsIdList);

    ASSERTL0(parseGood && !m_boundaryRegionsIdList.empty(),
             (std::string("Unable to read boundary regions index "
              "range for FilterAeroForces: ") + IndString).c_str());

    // determine what boundary regions need to be considered
    unsigned int numBoundaryRegions
                    = pFields[0]->GetBndConditions().size();

    m_boundaryRegionIsInList.insert(m_boundaryRegionIsInList.end(),
                                    numBoundaryRegions, 0);

    SpatialDomains::BoundaryConditions bcs(m_session,pFields[0]->GetGraph());

    const SpatialDomains::BoundaryRegionCollection &bregions
                    = bcs.GetBoundaryRegions();

    int cnt = 0;
    for (auto &it : bregions)
    {
        if ( std::find(m_boundaryRegionsIdList.begin(),
                       m_boundaryRegionsIdList.end(), it.first) !=
                m_boundaryRegionsIdList.end() )
        {
            m_boundaryRegionIsInList[cnt] = 1;
        }
        ++cnt;
    }

    LibUtilities::CommSharedPtr vComm = pFields[0]->GetComm();

    std::string solver_type = m_session->GetSolverInfo("SolverType");
    if(boost::iequals(solver_type, "Mapping"))
    {
        if (vComm->GetRank() == 0)
        {
            // Open output stream for cable forces
            m_outputStream[0].open(m_outputFile_fce.c_str());
            m_outputStream[0] << "#";
            m_outputStream[0].width(7);
            m_outputStream[0] << "Time";
            m_outputStream[0].width(15);
            m_outputStream[0] << "z";
            m_outputStream[0].width(15);
            m_outputStream[0] << "Fx (press)";
            m_outputStream[0].width(15);
            m_outputStream[0] << "Fx (visc)";
            m_outputStream[0].width(15);
            m_outputStream[0] << "Fx (tot)";
            m_outputStream[0].width(15);
            m_outputStream[0] << "Fy (press)";
            m_outputStream[0].width(15);
            m_outputStream[0] << "Fy (visc)";
            m_outputStream[0].width(15);
            m_outputStream[0] << "Fy (tot)";
            m_outputStream[0] << endl;

            // Open output stream for cable motions
            m_outputStream[1].open(m_outputFile_mot.c_str());
            m_outputStream[1] << "#";
            m_outputStream[1].width(7);
            m_outputStream[1] << "Time";
            m_outputStream[1].width(15);
            m_outputStream[1] << "z";
            m_outputStream[1].width(15);
            m_outputStream[1] << "Disp_x";
            m_outputStream[1].width(15);
            m_outputStream[1] << "Vel_x";
            m_outputStream[1].width(15);
            m_outputStream[1] << "Acel_x";
            m_outputStream[1].width(15);
            m_outputStream[1] << "Disp_y";
            m_outputStream[1].width(15);
            m_outputStream[1] << "Vel_y";
            m_outputStream[1].width(15);
            m_outputStream[1] << "Acel_y";
            m_outputStream[1] << endl;
        }
    }
    else if(boost::iequals(solver_type, "VelocityCorrectionScheme"))
    {
      if (vComm->GetRank() == 0)
      {
          // Open output stream for cable forces
          m_outputStream[0].open(m_outputFile_fce.c_str());
          m_outputStream[0] << "#";
          m_outputStream[0].width(7);
          m_outputStream[0] << "Time";
          m_outputStream[0].width(15);
          m_outputStream[0] << "Fx (press)";
          m_outputStream[0].width(15);
          m_outputStream[0] << "Fx (visc)";
          m_outputStream[0].width(15);
          m_outputStream[0] << "Fx (tot)";
          m_outputStream[0].width(15);
          m_outputStream[0] << "Fy (press)";
          m_outputStream[0].width(15);
          m_outputStream[0] << "Fy (visc)";
          m_outputStream[0].width(15);
          m_outputStream[0] << "Fy (tot)";
          m_outputStream[0] << endl;

          // Open output stream for cable motions
          m_outputStream[1].open(m_outputFile_mot.c_str());
          m_outputStream[1] << "#";
          m_outputStream[1].width(7);
          m_outputStream[1] << "Time";
          m_outputStream[1].width(15);
          m_outputStream[1] << "Disp_x";
          m_outputStream[1].width(15);
          m_outputStream[1] << "Vel_x";
          m_outputStream[1].width(15);
          m_outputStream[1] << "Acel_x";
          m_outputStream[1].width(15);
          m_outputStream[1] << "Disp_y";
          m_outputStream[1].width(15);
          m_outputStream[1] << "Vel_y";
          m_outputStream[1].width(15);
          m_outputStream[1] << "Acel_y";
          m_outputStream[1] << endl;
      }
    }
    
    int nt  = pFields[0]->GetNpoints();
    int dim = pFields.size()-1;
    m_baseflow = Array<OneD, Array<OneD, NekDouble> >(dim+1);
    for(int j = 0; j < dim+1; ++j)
    {
        m_baseflow[j] = Array<OneD, const NekDouble> (nt, 0.0);
    }
    
    // Getting the base flow fields
    string file = m_session->GetFunctionFilename("Baseflow", 0);
    ImportFldBound_base(file,pFields);


    
}


/**
 *
 */
void FilterMovingBody::UpdateForce(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
              Array<OneD, NekDouble> &Aeroforces,
        const NekDouble &time)
{
    int n, cnt, elmtid, nq, offset, boundary;
    int nt  = pFields[0]->GetNpoints();
    int dim = pFields.size()-1;
    int local_planes, Num_z_pos;

    StdRegions::StdExpansionSharedPtr elmt;
    Array<OneD, int> BoundarytoElmtID;
    Array<OneD, int> BoundarytoTraceID;
    Array<OneD, unsigned int> ZIDs;
    Array<OneD, MultiRegions::ExpListSharedPtr>  BndExp;

    Array<OneD, const NekDouble> P(nt);
    Array<OneD, const NekDouble> U(nt);
    Array<OneD, const NekDouble> V(nt);
    Array<OneD, const NekDouble> W(nt);

    adj_drag_t = 0.;
    adj_lift_t = 0.;
    adj_lift_p = 0.;
    adj_drag_p = 0.;
  
    Array<OneD, Array<OneD, NekDouble> > gradU(dim);
    Array<OneD, Array<OneD, NekDouble> > gradV(dim);
    Array<OneD, Array<OneD, NekDouble> > gradW(dim);
    Array<OneD, Array<OneD, NekDouble> > gradUbase(dim);
    Array<OneD, Array<OneD, NekDouble> > gradUbase1(dim*dim);
    Array<OneD, Array<OneD, NekDouble> > gradVbase(dim);
    Array<OneD, Array<OneD, NekDouble> > gradVbase1(dim*dim);
    Array<OneD, Array<OneD, NekDouble> > gradPbase(dim);
 

    Array<OneD, Array<OneD, NekDouble> > fgradU(dim);
    Array<OneD, Array<OneD, NekDouble> > fgradV(dim);
    Array<OneD, Array<OneD, NekDouble> > fgradW(dim);
    Array<OneD, Array<OneD, NekDouble> > fgradUbase(dim);
    Array<OneD, Array<OneD, NekDouble> > fgradVbase(dim);
    Array<OneD, Array<OneD, NekDouble> > fgradUbase1(dim*dim);
    Array<OneD, Array<OneD, NekDouble> > fgradVbase1(dim*dim);

    LibUtilities::CommSharedPtr vComm     = pFields[0]->GetComm();
    LibUtilities::CommSharedPtr vRowComm  = vComm->GetRowComm();
    LibUtilities::CommSharedPtr vColComm  = vComm->GetColumnComm();
    LibUtilities::CommSharedPtr vCColComm = vColComm->GetColumnComm();

    // set up storage space for forces on all the planes (z-positions)
    // on each processors some of them will remain empty as we may
    // have just few planes per processor
    std::string solver_type = m_session->GetSolverInfo("SolverType");
    std::string evol_operator = m_session->GetSolverInfo("EvolutionOperator");
    if(boost::iequals(solver_type, "Mapping"))
    {
        Num_z_pos = pFields[0]->GetHomogeneousBasis()->GetNumModes();
    }
    else
    {
        Num_z_pos = 1;
    }
    
    Array<OneD, NekDouble> Fx(Num_z_pos,0.0);
    Array<OneD, NekDouble> Fxp(Num_z_pos,0.0);
    Array<OneD, NekDouble> Fxv(Num_z_pos,0.0);
    Array<OneD, NekDouble> Fy(Num_z_pos,0.0);
    Array<OneD, NekDouble> Fyp(Num_z_pos,0.0);
    Array<OneD, NekDouble> Fyv(Num_z_pos,0.0);


    NekDouble rho = (pSession->DefinesParameter("rho"))
            ? (pSession->GetParameter("rho"))
            : 1;
    NekDouble mu = rho*pSession->GetParameter("Kinvis");

    if(boost::iequals(solver_type, "Mapping"))
    {

        for(int i = 0; i < pFields.size(); ++i)
        {
            pFields[i]->SetWaveSpace(false);
            pFields[i]->BwdTrans(pFields[i]->GetCoeffs(),
                                 pFields[i]->UpdatePhys());
            pFields[i]->SetPhysState(true);
        }



        // Get the number of local planes on the process and their IDs
        // to properly locate the forces in the Fx, Fy etc. vectors.

        ZIDs = pFields[0]->GetZIDs();
        local_planes = ZIDs.size();

        // Homogeneous 1D case  Compute forces on all WALL boundaries
        // This only has to be done on the zero (mean) Fourier mode.
        for(int plane = 0 ; plane < local_planes; plane++)
        {
            pFields[0]->GetPlane(plane)->GetBoundaryToElmtMap(BoundarytoElmtID,
                                                              BoundarytoTraceID);
            BndExp = pFields[0]->GetPlane(plane)->GetBndCondExpansions();
            StdRegions::StdExpansionSharedPtr bc;

            // loop over the types of boundary conditions
            for(cnt = n = 0; n < BndExp.size(); ++n)
            {
                if(m_boundaryRegionIsInList[n] == 1)
                {
                    for(int i = 0; i <  BndExp[n]->GetExpSize(); ++i, cnt++)
                    {
                        // find element of this expansion.
                        elmtid = BoundarytoElmtID[cnt];
                        elmt   = pFields[0]->GetPlane(plane)->GetExp(elmtid);
                        nq     = elmt->GetTotPoints();
                        offset = pFields[0]->GetPlane(plane)
                                           ->GetPhys_Offset(elmtid);

                        // Initialise local arrays for the velocity
                        // gradients size of total number of quadrature
                        // points for each element (hence local).
                        for(int j = 0; j < dim; ++j)
                        {
                            gradU[j] = Array<OneD, NekDouble>(nq,0.0);
                            gradV[j] = Array<OneD, NekDouble>(nq,0.0);
                            gradW[j] = Array<OneD, NekDouble>(nq,0.0);
                        }

                        // identify boundary of element
                        boundary = BoundarytoTraceID[cnt];

                        // Extract  fields
                        U = pFields[0]->GetPlane(plane)->GetPhys() + offset;
                        V = pFields[1]->GetPlane(plane)->GetPhys() + offset;
                        P = pFields[3]->GetPlane(plane)->GetPhys() + offset;

                        // compute the gradients
                        elmt->PhysDeriv(U,gradU[0],gradU[1]);
                        elmt->PhysDeriv(V,gradV[0],gradV[1]);

                        // Get face 1D expansion from element expansion
                        bc = BndExp[n]->GetExp(i)->as<LocalRegions::Expansion1D>();

                        // number of points on the boundary
                        int nbc = bc->GetTotPoints();

                        // several vectors for computing the forces
                        Array<OneD, NekDouble> Pb(nbc,0.0);

                        for(int j = 0; j < dim; ++j)
                        {
                            fgradU[j] = Array<OneD, NekDouble>(nbc,0.0);
                            fgradV[j] = Array<OneD, NekDouble>(nbc,0.0);
                        }

                        Array<OneD, NekDouble>  drag_t(nbc,0.0);
                        Array<OneD, NekDouble>  lift_t(nbc,0.0);
                        Array<OneD, NekDouble>  drag_p(nbc,0.0);
                        Array<OneD, NekDouble>  lift_p(nbc,0.0);
                        Array<OneD, NekDouble>  temp(nbc,0.0);
                        Array<OneD, NekDouble>  temp2(nbc,0.0);

                        // identify boundary of element .
                        boundary = BoundarytoTraceID[cnt];

                        // extraction of the pressure and wss on the
                        // boundary of the element
                        elmt->GetEdgePhysVals(boundary,bc,P,Pb);

                        for(int j = 0; j < dim; ++j)
                        {
                            elmt->GetEdgePhysVals(boundary,bc,gradU[j],fgradU[j]);
                            elmt->GetEdgePhysVals(boundary,bc,gradV[j],fgradV[j]);
                        }

                        //normals of the element
                        const Array<OneD, Array<OneD, NekDouble> > &normals
                                                = elmt->GetEdgeNormal(boundary);

                        //
                        // Compute viscous tractive forces on wall from
                        //
                        //  t_i  = - T_ij * n_j  (minus sign for force
                        //                        exerted BY fluid ON wall),
                        //
                        // where
                        //
                        //  T_ij = viscous stress tensor (here in Cartesian
                        //         coords)
                        //                          dU_i    dU_j
                        //       = RHO * KINVIS * ( ----  + ---- ) .
                        //                          dx_j    dx_i

                        //a) DRAG TERMS
                        //-rho*kinvis*(2*du/dx*nx+(du/dy+dv/dx)*ny)

                        Vmath::Vadd(nbc, fgradU[1], 1, fgradV[0],  1, drag_t,    1);
                        Vmath::Vmul(nbc, drag_t,    1, normals[1], 1, drag_t,    1);

                        Vmath::Smul(nbc, 2.0,          fgradU[0],  1, fgradU[0], 1);
                        Vmath::Vmul(nbc, fgradU[0], 1, normals[0], 1, temp2,     1);
                        Vmath::Smul(nbc, 0.5,          fgradU[0],  1, fgradU[0], 1);

                        Vmath::Vadd(nbc, temp2,     1, drag_t,     1, drag_t,    1);
                        Vmath::Smul(nbc, -mu,          drag_t,     1, drag_t,    1);

                        //zero temporary storage vector
                        Vmath::Zero(nbc, temp,  0.0);
                        Vmath::Zero(nbc, temp2, 0.0);


                        //b) LIFT TERMS
                        //-rho*kinvis*(2*dv/dy*ny+(du/dy+dv/dx)*nx)

                        Vmath::Vadd(nbc, fgradU[1], 1, fgradV[0],  1, lift_t,    1);
                        Vmath::Vmul(nbc, lift_t,    1, normals[0], 1, lift_t,    1);

                        Vmath::Smul(nbc, 2.0,          fgradV[1],  1, fgradV[1], 1);
                        Vmath::Vmul(nbc, fgradV[1], 1, normals[1], 1, temp2,     1);
                        Vmath::Smul(nbc, -0.5,         fgradV[1],  1, fgradV[1], 1);


                        Vmath::Vadd(nbc, temp2,     1, lift_t,     1, lift_t,    1);
                        Vmath::Smul(nbc, -mu,          lift_t,     1, lift_t,    1);

                        // Compute normal tractive forces on all WALL
                        // boundaries

                        Vmath::Vvtvp(nbc, Pb,       1, normals[0], 1, drag_p,    1,
                                                                      drag_p,    1);
                        Vmath::Vvtvp(nbc, Pb,       1, normals[1], 1, lift_p,    1,
                                                                      lift_p,    1);

                        //integration over the boundary
                        Fxv[ZIDs[plane]] += bc->Integral(drag_t);
                        Fyv[ZIDs[plane]] += bc->Integral(lift_t);

                        Fxp[ZIDs[plane]] += bc->Integral(drag_p);
                        Fyp[ZIDs[plane]] += bc->Integral(lift_p);
                    }
                }
                else
                {
                    cnt += BndExp[n]->GetExpSize();
                }
            }
        }

        for(int i = 0; i < pFields.size(); ++i)
        {
            pFields[i]->SetWaveSpace(true);
            pFields[i]->BwdTrans(pFields[i]->GetCoeffs(),
                                 pFields[i]->UpdatePhys());
            pFields[i]->SetPhysState(false);
        }

        // In case we are using an hybrid parallelisation we need
        // to reduce the forces on the same plane coming from
        // different mesh partitions.
        // It is quite an expensive communication, therefore
        // we check to make sure it is actually required.

        if(vComm->GetRowComm()->GetSize() > 0)
        {
            // NOTE 1: We can eventually sum the viscous and pressure
            // component before doing the communication, thus
            // reducing by a factor 2 the communication.
            // NOTE 2: We may want to set up in the Comm class an AllReduce
            // routine wich can handle arrays more efficiently
            for(int plane = 0 ; plane < local_planes; plane++)
            {
                vRowComm->AllReduce(Fxp[ZIDs[plane]], LibUtilities::ReduceSum);
                vRowComm->AllReduce(Fxv[ZIDs[plane]], LibUtilities::ReduceSum);
                vRowComm->AllReduce(Fyp[ZIDs[plane]], LibUtilities::ReduceSum);
                vRowComm->AllReduce(Fyv[ZIDs[plane]], LibUtilities::ReduceSum);
            }
        }

        // At this point on rank (0,n) of the Mesh partion communicator we have
        // the total areo forces on the planes which are on the same column
        // communicator. Since the planes are scattered on different processors
        // some of the entries of the vector Fxp, Fxp etc. are still zero.
        // Now we need to reduce the values on a single vector on rank (0,0) of the
        // global communicator.
        if(!pSession->DefinesSolverInfo("HomoStrip"))
        {
            if(vComm->GetRowComm()->GetRank() == 0)
            {
                for(int z = 0 ; z < Num_z_pos; z++)
                {
                    vColComm->AllReduce(Fxp[z], LibUtilities::ReduceSum);
                    vColComm->AllReduce(Fxv[z], LibUtilities::ReduceSum);
                    vColComm->AllReduce(Fyp[z], LibUtilities::ReduceSum);
                    vColComm->AllReduce(Fyv[z], LibUtilities::ReduceSum);
                }
            }
        }
        else
        {
            if(vComm->GetRowComm()->GetRank() == 0)
            {
                for(int z = 0 ; z < Num_z_pos; z++)
                {
                    vCColComm->AllReduce(Fxp[z], LibUtilities::ReduceSum);
                    vCColComm->AllReduce(Fxv[z], LibUtilities::ReduceSum);
                    vCColComm->AllReduce(Fyp[z], LibUtilities::ReduceSum);
                    vCColComm->AllReduce(Fyv[z], LibUtilities::ReduceSum);
                }
            }
        }

        if(!pSession->DefinesSolverInfo("HomoStrip"))
        {
            //set the forces imparted on the cable's wall
            for(int plane = 0 ; plane < local_planes; plane++)
            {
                Aeroforces[plane]                = Fxp[ZIDs[plane]]
                                                 + Fxv[ZIDs[plane]];
                Aeroforces[plane + local_planes] = Fyp[ZIDs[plane]]
                                                 + Fyv[ZIDs[plane]];
            }

            // Only output every m_outputFrequency.
            if ((m_index_f++) % m_outputFrequency)
            {
                return;
            }

            // At thi point in rank (0,0) we have the full vectors
            // containing Fxp,Fxv,Fyp and Fyv where different positions
            // in the vectors correspond to different planes.
            // Here we write it to file. We do it just on one porcess

            Array<OneD, NekDouble> z_coords(Num_z_pos,0.0);
            Array<OneD, const NekDouble> pts
                                = pFields[0]->GetHomogeneousBasis()->GetZ();

            NekDouble LZ;
            pSession->LoadParameter("LZ", LZ);
            Vmath::Smul(Num_z_pos,LZ/2.0,pts,1,z_coords,1);
            Vmath::Sadd(Num_z_pos,LZ/2.0,z_coords,1,z_coords,1);
            if (vComm->GetRank() == 0)
            {

                Vmath::Vadd(Num_z_pos,Fxp,1,Fxv,1,Fx,1);
                Vmath::Vadd(Num_z_pos,Fyp,1,Fyv,1,Fy,1);

                for(int i = 0 ; i < Num_z_pos; i++)
                {
                    m_outputStream[0].width(8);
                    m_outputStream[0] << setprecision(6) << time;

                    m_outputStream[0].width(15);
                    m_outputStream[0] << setprecision(6) << z_coords[i];

                    m_outputStream[0].width(15);
                    m_outputStream[0] << setprecision(8) << Fxp[i];
                    m_outputStream[0].width(15);
                    m_outputStream[0] << setprecision(8) << Fxv[i];
                    m_outputStream[0].width(15);
                    m_outputStream[0] << setprecision(8) << Fx[i];

                    m_outputStream[0].width(15);
                    m_outputStream[0] << setprecision(8) << Fyp[i];
                    m_outputStream[0].width(15);
                    m_outputStream[0] << setprecision(8) << Fyv[i];
                    m_outputStream[0].width(15);
                    m_outputStream[0] << setprecision(8) << Fy[i];
                    m_outputStream[0] << endl;
                }
            }
        }
        else
        {
            // average the forces over local planes for each processor
            Array<OneD, NekDouble> fces(6,0.0);
            for(int plane = 0 ; plane < local_planes; plane++)
            {
                fces[0] += Fxp[ZIDs[plane]] + Fxv[ZIDs[plane]];
                fces[1] += Fyp[ZIDs[plane]] + Fyv[ZIDs[plane]];
                fces[2] += Fxp[ZIDs[plane]] ;
                fces[3] += Fyp[ZIDs[plane]] ;
                fces[4] += Fxv[ZIDs[plane]] ;
                fces[5] += Fyv[ZIDs[plane]] ;
            }

            fces[0] = fces[0]/local_planes;
            fces[1] = fces[1]/local_planes;
            fces[2] = fces[2]/local_planes;
            fces[3] = fces[3]/local_planes;
            fces[4] = fces[4]/local_planes;
            fces[5] = fces[5]/local_planes;

            // average the forces over communicators within each strip
            vCColComm->AllReduce(fces[0], LibUtilities::ReduceSum);
            vCColComm->AllReduce(fces[1], LibUtilities::ReduceSum);
            vCColComm->AllReduce(fces[2], LibUtilities::ReduceSum);
            vCColComm->AllReduce(fces[3], LibUtilities::ReduceSum);
            vCColComm->AllReduce(fces[4], LibUtilities::ReduceSum);
            vCColComm->AllReduce(fces[5], LibUtilities::ReduceSum);

            int npts = vComm->GetColumnComm()->GetColumnComm()->GetSize();

            fces[0] = fces[0]/npts;
            fces[1] = fces[1]/npts;
            fces[2] = fces[2]/npts;
            fces[3] = fces[3]/npts;
            fces[4] = fces[4]/npts;
            fces[5] = fces[5]/npts;

            for(int plane = 0 ; plane < local_planes; plane++)
            {
                Aeroforces[plane]              = fces[0];
                Aeroforces[plane+local_planes] = fces[1];
            }

            // Only output every m_outputFrequency.
            if ((m_index_f) % m_outputFrequency)
            {
                return;
            }

            int colrank = vColComm->GetRank();
            int nstrips;

            NekDouble DistStrip;

            pSession->LoadParameter("Strip_Z", nstrips);
            pSession->LoadParameter("DistStrip", DistStrip);

            Array<OneD, NekDouble> z_coords(nstrips);
            for(int i = 0; i < nstrips; i++)
            {
                z_coords[i] = i * DistStrip;
            }

            if(colrank == 0)
            {
                m_outputStream[0].width(8);
                m_outputStream[0] << setprecision(6) << time;

                m_outputStream[0].width(15);
                m_outputStream[0] << setprecision(6) << z_coords[0];

                m_outputStream[0].width(15);
                m_outputStream[0] << setprecision(8) << fces[2];
                m_outputStream[0].width(15);
                m_outputStream[0] << setprecision(8) << fces[4];
                m_outputStream[0].width(15);
                m_outputStream[0] << setprecision(8) << fces[0];

                m_outputStream[0].width(15);
                m_outputStream[0] << setprecision(8) << fces[3];
                m_outputStream[0].width(15);
                m_outputStream[0] << setprecision(8) << fces[5];
                m_outputStream[0].width(15);
                m_outputStream[0] << setprecision(8) << fces[1];
                m_outputStream[0] << endl;

                for(int i = 1; i < nstrips; i++)
                {
                    vColComm->Recv(i, fces);

                    m_outputStream[0].width(8);
                    m_outputStream[0] << setprecision(6) << time;

                    m_outputStream[0].width(15);
                    m_outputStream[0] << setprecision(6) << z_coords[i];

                    m_outputStream[0].width(15);
                    m_outputStream[0] << setprecision(8) << fces[2];
                    m_outputStream[0].width(15);
                    m_outputStream[0] << setprecision(8) << fces[4];
                    m_outputStream[0].width(15);
                    m_outputStream[0] << setprecision(8) << fces[0];

                    m_outputStream[0].width(15);
                    m_outputStream[0] << setprecision(8) << fces[3];
                    m_outputStream[0].width(15);
                    m_outputStream[0] << setprecision(8) << fces[5];
                    m_outputStream[0].width(15);
                    m_outputStream[0] << setprecision(8) << fces[1];
                    m_outputStream[0] << endl;
                }
            }
            else
            {
                for(int i = 1; i < nstrips; i++)
                {
                    if(colrank == i)
                    {
                        vColComm->Send(0, fces);
                    }
                }
            }
        }
        

    }
    // Only twodimensional
    else if(boost::iequals(solver_type, "VelocityCorrectionScheme"))
    {

        for(int i = 0; i < pFields.size(); ++i)
        {
            pFields[i]->BwdTrans(pFields[i]->GetCoeffs(),
                                 pFields[i]->UpdatePhys());
            pFields[i]->SetPhysState(true);
        }

        pFields[0]->GetBoundaryToElmtMap(BoundarytoElmtID,
                                                          BoundarytoTraceID);
        BndExp = pFields[0]->GetBndCondExpansions();
        StdRegions::StdExpansionSharedPtr bc;
        
        MultiRegions::ExpListSharedPtr                   locExpList;
        locExpList = BndExp[1];
        int b_points = BndExp[1]->GetNpoints();

        // loop over the types of boundary conditions
        for(cnt = n = 0; n < BndExp.size(); ++n)
        {
            if(m_boundaryRegionIsInList[n] == 1)
            {
                for(int i = 0; i <  BndExp[n]->GetExpSize(); ++i, cnt++)
                {
                    // find element of this expansion.
                    elmtid = BoundarytoElmtID[cnt];
                    elmt   = pFields[0]->GetExp(elmtid);
                    nq     = elmt->GetTotPoints();
                    offset = pFields[0]->GetPhys_Offset(elmtid);

                    // Initialise local arrays for the velocity
                    // gradients size of total number of quadrature
                    // points for each element (hence local).
                    for(int j = 0; j < dim; ++j)
                    {
                        gradU[j] = Array<OneD, NekDouble>(nq,0.0);
                        gradV[j] = Array<OneD, NekDouble>(nq,0.0);
    
                    }

                    // identify boundary of element
                    boundary = BoundarytoTraceID[cnt];

                    // Extract  fields
                    U = pFields[0]->GetPhys() + offset;
                    V = pFields[1]->GetPhys() + offset;
                    P = pFields[2]->GetPhys() + offset;

                    
                    // compute the gradients
                    elmt->PhysDeriv(U,gradU[0],gradU[1]);
                    elmt->PhysDeriv(V,gradV[0],gradV[1]);
                    

                    // Get face 1D expansion from element expansion
                    bc = BndExp[n]->GetExp(i)->as<LocalRegions::Expansion1D>();


                    // number of points on the boundary
                    int nbc = bc->GetTotPoints();

                    // several vectors for computing the forces
                    Array<OneD, NekDouble> Pb(nbc,0.0);

                    for(int j = 0; j < dim; ++j)
                    {
                        fgradU[j] = Array<OneD, NekDouble>(nbc,0.0);
                        fgradV[j] = Array<OneD, NekDouble>(nbc,0.0);
            
                    }

                    Array<OneD, NekDouble>  drag_t(nbc,0.0);
                    Array<OneD, NekDouble>  lift_t(nbc,0.0);
                    Array<OneD, NekDouble>  drag_p(nbc,0.0);
                    Array<OneD, NekDouble>  lift_p(nbc,0.0);
                    Array<OneD, NekDouble>  temp(nbc,0.0);
                    Array<OneD, NekDouble>  temp2(nbc,0.0);
                    Array<OneD, Array<OneD, NekDouble> >       fp(dim);
                    Array<OneD, Array<OneD, NekDouble> >       fv(dim);


                    // identify boundary of element .
                    boundary = BoundarytoTraceID[cnt];

                    // extraction of the pressure and wss on the
                    // boundary of the element
                    elmt->GetEdgePhysVals(boundary,bc,P,Pb);

                    for(int j = 0; j < dim; ++j)
                    {
                        elmt->GetEdgePhysVals(boundary,bc,gradU[j],fgradU[j]);
                        elmt->GetEdgePhysVals(boundary,bc,gradV[j],fgradV[j]);
            
                    }

                    //normals of the element
                    const Array<OneD, Array<OneD, NekDouble> > &normals
                                            = elmt->GetEdgeNormal(boundary);

                    // Compute viscous tractive forces on wall from
                    //
                    //  t_i  = - T_ij * n_j  (minus sign for force
                    //                        exerted BY fluid ON wall),
                    //
                    // where
                    //
                    //  T_ij = viscous stress tensor (here in Cartesian
                    //         coords)
                    //                          dU_i    dU_j
                    //       = RHO * KINVIS * ( ----  + ---- ) .
                    //                          dx_j    dx_i

                    //a) DRAG TERMS
                    //-rho*kinvis*(2*du/dx*nx+(du/dy+dv/dx)*ny)

                    Vmath::Vadd(nbc, fgradU[1], 1, fgradV[0],  1, drag_t,    1);
                    Vmath::Vmul(nbc, drag_t,    1, normals[1], 1, drag_t,    1);

                    Vmath::Smul(nbc, 2.0,          fgradU[0],  1, fgradU[0], 1);
                    Vmath::Vmul(nbc, fgradU[0], 1, normals[0], 1, temp2,     1);
                    Vmath::Smul(nbc, 0.5,          fgradU[0],  1, fgradU[0], 1);

                    Vmath::Vadd(nbc, temp2,     1, drag_t,     1, drag_t,    1);
                    Vmath::Smul(nbc, -mu,          drag_t,     1, drag_t,    1);

                    //zero temporary storage vector
                    Vmath::Zero(nbc, temp,  0.0);
                    Vmath::Zero(nbc, temp2, 0.0);
        

                    //b) LIFT TERMS
                    //-rho*kinvis*(2*dv/dy*ny+(du/dy+dv/dx)*nx)

                    Vmath::Vadd(nbc, fgradU[1], 1, fgradV[0],  1, lift_t,    1);
                    Vmath::Vmul(nbc, lift_t,    1, normals[0], 1, lift_t,    1);

                    Vmath::Smul(nbc, 2.0,          fgradV[1],  1, fgradV[1], 1);
                    Vmath::Vmul(nbc, fgradV[1], 1, normals[1], 1, temp2,     1);
                    Vmath::Smul(nbc, -0.5,         fgradV[1],  1, fgradV[1], 1);


                    Vmath::Vadd(nbc, temp2,     1, lift_t,     1, lift_t,    1);
                    Vmath::Smul(nbc, -mu,          lift_t,     1, lift_t,    1);

                    // Compute normal tractive forces on all WALL
                    // boundaries
                    
                    Vmath::Vmul(nbc, Pb,       1, normals[0], 1, drag_p,    1);
                    Vmath::Vmul(nbc, Pb,       1, normals[1], 1, lift_p,    1);

                    
                    // aeroforces
                    Fxv[0] += bc->Integral(drag_t);
                    Fyv[0] += bc->Integral(lift_t);

                    Fxp[0] += bc->Integral(drag_p);
                    Fyp[0] += bc->Integral(lift_p);

                    if (boost::iequals(evol_operator, "Adjoint"))
                    {
                        // Array<OneD, NekDouble>  drag_t(nbc,0.0);
                        // Array<OneD, NekDouble>  lift_t(nbc,0.0);
                        // Array<OneD, NekDouble>  drag_p(nbc,0.0);
                        // Array<OneD, NekDouble>  lift_p(nbc,0.0);

                        //zero temporary storage vector
                        Vmath::Zero(nbc, temp,  0.0);
                        Vmath::Zero(nbc, temp2, 0.0);
                       

                        Array<OneD, NekDouble>  UElmt(nt);
                        Array<OneD, NekDouble>  VElmt(nt);
                        
                        for(int j = 0; j < dim; ++j)
                        {
                            gradUbase[j] = Array<OneD, NekDouble>(nq,0.0);
                            gradVbase[j] = Array<OneD, NekDouble>(nq,0.0);
                            
                            fgradUbase[j] = Array<OneD, NekDouble>(nbc,0.0);
                            fgradVbase[j] = Array<OneD, NekDouble>(nbc,0.0);
                        }

                        for(int j = 0; j < dim*dim; ++j)
                        {
                            gradUbase1[j] = Array<OneD, NekDouble>(nq,0.0);
                            gradVbase1[j] = Array<OneD, NekDouble>(nq,0.0);
                            fgradUbase1[j] = Array<OneD, NekDouble>(nbc,0.0);
                            fgradVbase1[j] = Array<OneD, NekDouble>(nbc,0.0);
                        }

                        UElmt = m_baseflow[0] + offset; 
                        VElmt = m_baseflow[1] + offset; 
                        elmt->PhysDeriv(UElmt,gradUbase[0],gradUbase[1]);
                        elmt->PhysDeriv(VElmt,gradVbase[0],gradVbase[1]);

                        
                        for(int j = 0; j < dim; ++j)
                        {
                            elmt->GetEdgePhysVals(boundary,bc,gradUbase[j],fgradUbase[j]);
                            elmt->GetEdgePhysVals(boundary,bc,gradVbase[j],fgradVbase[j]);
                        }


                         //a) Adj DRAG TERMS
                        // dU/dx*(drag_t + drag_p) + dV/dx*(lift_t + lift_p)             
                        for (int k = 0; k < nbc; ++k)
                        {

                            temp[k]  = (drag_t[k]+drag_p[k])*fgradUbase[0][k];
                            temp2[k] = (lift_t[k]+lift_p[k])*fgradVbase[0][k];

                        }
                        
                         // aeroforces
                        adj_drag_t += bc->Integral(temp);
                        adj_drag_p += bc->Integral(temp2);
                        
                        //zero temporary storage vector
                        Vmath::Zero(nbc, temp,  0.0);
                        Vmath::Zero(nbc, temp2, 0.0);
                        
                        for (int k = 0; k < nbc; ++k)
                        {

                            temp[k]  = (drag_t[k]+drag_p[k])*fgradUbase[1][k];
                            temp2[k] = (lift_t[k]+lift_p[k])*fgradVbase[1][k];
                            
                        }
                        // aeroforces
                        adj_lift_t += bc->Integral(temp);
                        adj_lift_p += bc->Integral(temp2);
                        

                    }
        
                    
                    if ((boost::iequals(evol_operator, "Direct") || boost::iequals(evol_operator, "Adjoint"))  && counter==0)
                    {   
                        Array<OneD, NekDouble>  pElmt(nt);
                        Array<OneD, NekDouble>  UElmt(nt);
                        Array<OneD, NekDouble>  VElmt(nt);
                        
                        for(int j = 0; j < dim; ++j)
                        {
                      
                            gradUbase[j] = Array<OneD, NekDouble>(nq,0.0);
                            gradVbase[j] = Array<OneD, NekDouble>(nq,0.0);
                            gradPbase[j] = Array<OneD, NekDouble>(nq,0.0);
                            fgradUbase[j] = Array<OneD, NekDouble>(nbc,0.0);
                            fgradVbase[j] = Array<OneD, NekDouble>(nbc,0.0);
                        }

                        for(int j = 0; j < dim*dim; ++j)
                        {
                      
                            gradUbase1[j] = Array<OneD, NekDouble>(nq,0.0);
                            gradVbase1[j] = Array<OneD, NekDouble>(nq,0.0);
                            fgradUbase1[j] = Array<OneD, NekDouble>(nbc,0.0);
                            fgradVbase1[j] = Array<OneD, NekDouble>(nbc,0.0);
                        }

                        
                        pElmt = m_baseflow[2] + offset;
                        UElmt = m_baseflow[0] + offset; 
                        VElmt = m_baseflow[1] + offset; 
                        elmt->PhysDeriv(UElmt,gradUbase[0],gradUbase[1]);
                        elmt->PhysDeriv(VElmt,gradVbase[0],gradVbase[1]);
                        elmt->PhysDeriv(pElmt,gradPbase[0],gradPbase[1]);

                        
                        for(int j = 0; j < dim; ++j)
                        {
                   
                            elmt->GetEdgePhysVals(boundary,bc,gradUbase[j],fgradUbase[j]);
                            elmt->GetEdgePhysVals(boundary,bc,gradVbase[j],fgradVbase[j]);
                        }

                        // several vectors for computing the forces
                        Array<OneD, NekDouble> Pb_base(nbc,0.0);
                        Array<OneD, NekDouble> Pbgradx_base(nbc,0.0);
                        Array<OneD, NekDouble> Pbgrady_base(nbc,0.0);
                        Array<OneD, NekDouble>  drag_t(nbc,0.0);
                        Array<OneD, NekDouble>  lift_t(nbc,0.0);
                        Array<OneD, NekDouble>  drag_p(nbc,0.0);
                        Array<OneD, NekDouble>  lift_p(nbc,0.0);
                        Array<OneD, NekDouble>  temp(nbc,0.0);
                        Array<OneD, NekDouble>  temp2(nbc,0.0);
                       // several vectors for computing the forces
                        
                        // int offset1 = i*nbc;

                        elmt->GetEdgePhysVals(boundary,bc,pElmt,Pb_base);

                         //a) DRAG TERMS
                        //-rho*kinvis*(2*du/dx*nx+(du/dy+dv/dx)*ny)

                        Vmath::Vadd(nbc, fgradUbase[1], 1, fgradVbase[0],  1, drag_t,    1);
                        Vmath::Vmul(nbc, drag_t,    1, normals[1], 1, drag_t,    1);

                        Vmath::Smul(nbc, 2.0,          fgradUbase[0],  1, fgradUbase[0], 1);
                        Vmath::Vmul(nbc, fgradUbase[0], 1, normals[0], 1, temp2,     1);
                        Vmath::Smul(nbc, 0.5,          fgradUbase[0],  1, fgradUbase[0], 1);

                        Vmath::Vadd(nbc, temp2,     1, drag_t,     1, drag_t,    1);
                        Vmath::Smul(nbc, -mu,          drag_t,     1, drag_t,    1);

                        //zero temporary storage vector
                        Vmath::Zero(nbc, temp,  0.0);
                        Vmath::Zero(nbc, temp2, 0.0);


                        //b) LIFT TERMS
                        //-rho*kinvis*(2*dv/dy*ny+(du/dy+dv/dx)*nx)

                        Vmath::Vadd(nbc, fgradUbase[1], 1, fgradVbase[0],  1, lift_t,    1);
                        Vmath::Vmul(nbc, lift_t,    1, normals[0], 1, lift_t,    1);

                        Vmath::Smul(nbc, 2.0,          fgradVbase[1],  1, fgradVbase[1], 1);
                        Vmath::Vmul(nbc, fgradVbase[1], 1, normals[1], 1, temp2,     1);
                        Vmath::Smul(nbc, -0.5,         fgradVbase[1],  1, fgradVbase[1], 1);


                        Vmath::Vadd(nbc, temp2,     1, lift_t,     1, lift_t,    1);
                        Vmath::Smul(nbc, -mu,          lift_t,     1, lift_t,    1);


                        // boundaries
                        Vmath::Vmul(nbc, Pb_base,       1, normals[0], 1, drag_p,    1);
                        Vmath::Vmul(nbc, Pb_base,       1, normals[1], 1, lift_p,    1);             
                    
                        
                        bf_drag_t += bc->Integral(drag_t);   
                        bf_lift_t += bc->Integral(lift_t);
                        bf_drag_p += bc->Integral(drag_p);   
                        bf_lift_p += bc->Integral(lift_p);

                       

                        

                        Vmath::Zero(nbc, temp,   0.0);
                        Vmath::Zero(nbc, temp2,  0.0);
                        Vmath::Zero(nbc, drag_t, 0.0);
                        Vmath::Zero(nbc, drag_p, 0.0);
                        Vmath::Zero(nbc, lift_t, 0.0);
                        Vmath::Zero(nbc, lift_p, 0.0);  

                        elmt->PhysDeriv(gradUbase[0],gradUbase1[0],gradUbase1[1]);
                        elmt->PhysDeriv(gradUbase[1],gradUbase1[2],gradUbase1[3]);
                        elmt->PhysDeriv(gradVbase[0],gradVbase1[0],gradVbase1[1]);
                        elmt->PhysDeriv(gradVbase[1],gradVbase1[2],gradVbase1[3]);
                        
                        for(int j = 0; j < dim*dim; ++j)
                        {
                            elmt->GetEdgePhysVals(boundary,bc,gradUbase1[j],fgradUbase1[j]);
                            elmt->GetEdgePhysVals(boundary,bc,gradVbase1[j],fgradVbase1[j]);
                        }

                        elmt->GetEdgePhysVals(boundary,bc,gradPbase[0],Pbgradx_base);
                        elmt->GetEdgePhysVals(boundary,bc,gradPbase[1],Pbgrady_base);

                       

                        //a) DRAG TERMS
                        //-rho*kinvis*(2*d2u/dx2*nx+(d2u/dy2+d2v/dydx)*ny)

                        Vmath::Vadd(nbc, fgradUbase1[3], 1, fgradVbase1[1],  1, drag_t,    1);
                        Vmath::Vmul(nbc, drag_t,    1, normals[1], 1, drag_t,    1);

                        Vmath::Smul(nbc, 2.0,          fgradUbase1[0],  1, fgradUbase1[0], 1);
                        Vmath::Vmul(nbc, fgradUbase1[0], 1, normals[0], 1, temp2,     1);
                        Vmath::Smul(nbc, 0.5,          fgradUbase1[0],  1, fgradUbase1[0], 1);

                        Vmath::Vadd(nbc, temp2,     1, drag_t,     1, drag_t,    1);
                        Vmath::Smul(nbc, -mu,          drag_t,     1, drag_t,    1);

                        //zero temporary storage vector
                        Vmath::Zero(nbc, temp,  0.0);
                        Vmath::Zero(nbc, temp2, 0.0);


                        //b) LIFT TERMS
                        //-rho*kinvis*(2*d2v/dy2*ny+(d2u/dydx+d2v/dx2)*nx)

                        Vmath::Vadd(nbc, fgradUbase1[1], 1, fgradVbase1[0],  1, lift_t,    1);
                        Vmath::Vmul(nbc, lift_t,    1, normals[0], 1, lift_t,    1);

                        Vmath::Smul(nbc, 2.0,          fgradVbase1[3],  1, fgradVbase1[3], 1);
                        Vmath::Vmul(nbc, fgradVbase1[3], 1, normals[1], 1, temp2,     1);
                        Vmath::Smul(nbc, -0.5,         fgradVbase1[3],  1, fgradVbase1[3], 1);


                        Vmath::Vadd(nbc, temp2,     1, lift_t,     1, lift_t,    1);
                        Vmath::Smul(nbc, -mu,          lift_t,     1, lift_t,    1);


                        // boundaries
                        Vmath::Vmul(nbc, Pbgradx_base, 1, normals[0], 1, drag_p,    1);
                        Vmath::Vmul(nbc, Pbgrady_base, 1, normals[1], 1, lift_p,    1);     
                            
                        

                        bfgrad_drag_t += bc->Integral(drag_t);   
                        bfgrad_lift_t += bc->Integral(lift_t);
                        bfgrad_drag_p += bc->Integral(drag_p);   
                        bfgrad_lift_p += bc->Integral(lift_p);

                      

                    }
                   

                        

                   



                
                }
            }
            else
            {
                cnt += BndExp[n]->GetExpSize();
            }

            
            if(boost::iequals(evol_operator, "Adjoint") )
            {
                Aeroforces[0] = Fxv[0] + Fxp[0];
                Aeroforces[1] = Fyv[0] + Fyp[0];
                Aeroforces[2] = bf_drag_t + bf_drag_p;
                Aeroforces[3] = bf_lift_t + bf_lift_p;
                Aeroforces[4] = bfgrad_drag_t + bfgrad_drag_p;
                Aeroforces[5] = bfgrad_lift_t + bfgrad_lift_p;
                Aeroforces[6] = adj_drag_t + adj_drag_p;
                Aeroforces[7] = adj_lift_t + adj_lift_p;
                
                
                

            }
            else if(boost::iequals(evol_operator, "Direct") )
            {
                Aeroforces[0] = Fxv[0] + Fxp[0];
                Aeroforces[1] = Fyv[0] + Fyp[0];
                Aeroforces[2] = bf_drag_t + bf_drag_p;
                Aeroforces[3] = bf_lift_t + bf_lift_p;
                Aeroforces[4] = bfgrad_drag_t + bfgrad_drag_p;
                Aeroforces[5] = bfgrad_lift_t + bfgrad_lift_p;
               
        
            }
            else
            {
                Aeroforces[0] = Fxv[0] + Fxp[0];
                Aeroforces[1] = Fyv[0] + Fyp[0];
            }
            

        }

    }
    else
    {
          ASSERTL0(false,
                  "This option is not enable.");
    }

    ++counter;
  

    
}


void FilterMovingBody::ImportFldBound_base(
std::string                                  pInfile,
const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields)
{
    std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef;
    std::vector<std::vector<NekDouble> >                 FieldData;

    int nqtot = m_baseflow[0].size();

    Array<OneD, NekDouble> tmp_coeff(pFields[0]->GetNcoeffs(), 0.0);


    int numexp = pFields[0]->GetExpSize();
    Array<OneD,int> ElementGIDs(numexp);

    // Define list of global element ids
    for(int i = 0; i < numexp; ++i)
    {
        ElementGIDs[i] = pFields[0]->GetExp(i)->GetGeom()->GetGlobalID();
    }

    //Get Homogeneous
    LibUtilities::FieldIOSharedPtr fld = LibUtilities::FieldIO::CreateForFile(
        m_session, pInfile);
    fld->Import(pInfile, FieldDef, FieldData,
                LibUtilities::NullFieldMetaDataMap,
                ElementGIDs);

    
    
    // fld->Import(pInfile, FieldDef, FieldData);
 
    
    int f;
    for (int j = 0; j < 3; ++j)
    {
        
        for(int i = 0; i < FieldDef.size(); ++i)
        {
            // f = 2*j + 4;
            bool flag = FieldDef[i]->m_fields[0] ==
                m_session->GetVariable(0);

            ASSERTL0(flag, (std::string("Order of ") + pInfile
                            + std::string(" data and that defined in "
                            "the session file differs")).c_str());

            pFields[j]->ExtractDataToCoeffs(
                                FieldDef[i],
                                FieldData[i],
                                FieldDef[i]->m_fields[j],
                                tmp_coeff);
            
            
           
            
        }

       
        
        pFields[0]->BwdTrans(tmp_coeff, m_baseflow[j]);
       
           
        

    }

   


}
/**
/**
 *
 */
void FilterMovingBody::UpdateMotion(
        const LibUtilities::SessionReaderSharedPtr              &pSession,
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
              Array<OneD, NekDouble>                            &MotionVars,
        const NekDouble                                         &time)
{
    // Only output every m_outputFrequency.
    if ((m_index_m++) % m_outputFrequency)
    {
        return;
    }

    int npts;
    // Length of the cable
    NekDouble Length;
    std::string evol_operator = m_session->GetSolverInfo("EvolutionOperator");
    std::string solver_type = m_session->GetSolverInfo("SolverType");
    
    if(boost::iequals(solver_type, "Mapping"))
    {
        if(!pSession->DefinesSolverInfo("HomoStrip"))
        {
            pSession->LoadParameter("LZ", Length);
            npts = m_session->GetParameter("HomModesZ");
        }
        else
        {
            pSession->LoadParameter("LC", Length);
            npts = m_session->GetParameter("HomStructModesZ");
        }
    }
    else if(boost::iequals(solver_type, "VelocityCorrectionScheme"))
    {
        npts = 1;
    }
    else
    {
          ASSERTL0(false,
                  "This option is not valid.");
    }


    NekDouble z_coords;
    for(int n = 0; n < npts; n++)
    {
        z_coords = Length/npts*n;
        m_outputStream[1].width(8);
        m_outputStream[1] << setprecision(6) << time;
        m_outputStream[1].width(15);
        // m_outputStream[1] << setprecision(6) << z_coords;
        // m_outputStream[1].width(15);
        m_outputStream[1] << setprecision(8) << MotionVars[n];
        m_outputStream[1].width(15);
        m_outputStream[1] << setprecision(8) << MotionVars[npts+n];
        m_outputStream[1].width(15);
        m_outputStream[1] << setprecision(8) << MotionVars[2*npts+n];
        m_outputStream[1].width(15);
        m_outputStream[1] << setprecision(8) << MotionVars[3*npts+n];
        m_outputStream[1].width(15);
        m_outputStream[1] << setprecision(8) << MotionVars[4*npts+n];
        m_outputStream[1].width(15);
        m_outputStream[1] << setprecision(8) << MotionVars[5*npts+n];
        m_outputStream[1] << endl;
    }
}


/**
 *
 */
void FilterMovingBody::v_Finalise(
        const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
        const NekDouble                                         &time)
{
    if (pFields[0]->GetComm()->GetRank() == 0)
    {
        m_outputStream[0].close();
        m_outputStream[1].close();
    }
}


/**
 *
 */
bool FilterMovingBody::v_IsTimeDependent()
{
    return true;
}
}
