///////////////////////////////////////////////////////////////////////////////
//
// File: ForcingMovingBody.cpp
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
// Description: Moving Body m_forcing (movement of a body in a domain is achieved
// via a m_forcing term which is the results of a coordinate system motion)
//
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/Forcing/ForcingMovingBody.h>
#include <MultiRegions/ExpList.h>
#include <SolverUtils/DriverModifiedArnoldi.h>
#include <SolverUtils/Filters/FilterAeroForces.h>
#include <SolverUtils/Forcing/ForcingBody.h>
using namespace std;

namespace Nektar
{
std::string ForcingMovingBody::className = SolverUtils::GetForcingFactory().
            RegisterCreatorFunction("MovingBody",
                                    ForcingMovingBody::create,
                                    "Moving Body Forcing");
NekDouble ForcingMovingBody::AdamsBashforth_coeffs[3][3] = {
        { 1.0       , 0.0       , 0.0     },
        { 3.0/2.0   ,-1.0/2.0   , 0.0     },
        { 23.0/12.0 ,-4.0/3.0   , 5.0/12.0}};
NekDouble ForcingMovingBody::AdamsMoulton_coeffs[3][3] = {
        { 1.0       ,  0.0      , 0.0     },
        { 1.0/2.0   ,  1.0/2.0  , 0.0     },
        { 5.0/12.0  ,  2.0/3.0  ,-1.0/12.0}};

int comp;
int c;
NekDouble y, y1, y_dagger, y1_dagger;
ForcingMovingBody::ForcingMovingBody(
                const LibUtilities::SessionReaderSharedPtr         &pSession,
                const std::weak_ptr<SolverUtils::EquationSystem> &pEquation)
    : Forcing(pSession, pEquation)
{
}

void ForcingMovingBody::v_InitObject(
        const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields,
        const unsigned int& pNumForcingFields,
        const TiXmlElement* pForce)
{
    
    std::string evol_operator = m_session->GetSolverInfo("EvolutionOperator");
    std::string solver_type = m_session->GetSolverInfo("SolverType");
    
    
    // At this point we know in the xml file where those quantities
    // are declared (equation or file) - via a function name which is now
    // stored in funcNameD etc. We need now to fill in with this info the
    // m_zta and m_eta vectors (actuallythey are matrices) Array to control if
    // the motion is determined by an equation or is from a file.(not Nektar++)
    // check if we need to load a file or we have an equation
    
    CheckIsFromFile(pForce);

    // Initialise movingbody filter
    InitialiseFilter(m_session, pFields, pForce);
    
    // Initialise the cable model

    InitialiseCableModel(m_session, pFields);
    int phystot = pFields[0]->GetTotPoints();
    if(boost::iequals(solver_type, "VCSMapping"))
    {
        // Load mapping
        m_mapping = GlobalMapping::Mapping::Load(m_session, pFields);
        m_mapping->SetTimeDependent( true );

        if(m_vdim > 0)
        {
            m_mapping->SetFromFunction( false );
        }
    }
  
    // Determine time integration order
    std::string intMethod = m_session->GetSolverInfo("TIMEINTEGRATIONMETHOD");

    if(boost::iequals(intMethod, "IMEXOrder1"))
    {
        m_intSteps = 1;
    }
    else if (boost::iequals(intMethod, "IMEXOrder2"))
    {
        m_intSteps = 2;
    }
    else if (boost::iequals(intMethod, "IMEXOrder3"))
    {
        m_intSteps = 3;
    }
    else
    {
        ASSERTL0(false, "Time integration method not supported.");
    }


    m_zta = Array<OneD, Array< OneD, NekDouble> > (3);
    m_eta = Array<OneD, Array< OneD, NekDouble> > (3);
    m_force = Array<OneD, Array<OneD,NekDouble> > (m_intSteps);
    m_force1 = Array<OneD, Array<OneD,NekDouble> > (m_intSteps);
    m_velocity     = Array<OneD, Array<OneD,NekDouble> > (m_intSteps);
    m_displacement = Array<OneD,NekDouble> (m_vdim, 0.0);
    m_previousDisp = Array<OneD,NekDouble> (m_vdim, 0.0);
    m_Aeroforces   = Array<OneD,NekDouble> (2, 0.0);
    if(boost::iequals(evol_operator, "Adjoint") || boost::iequals(evol_operator, "TransientGrowth"))
    {
        m_Aeroforces   = Array<OneD,NekDouble> (8, 0.0);
    }
    else if(boost::iequals(evol_operator, "Direct") || (boost::iequals(evol_operator, "TransientGrowth") && c==0))
    {
        m_Aeroforces   = Array<OneD,NekDouble> (6, 0.0);
    }
    else
    {
        m_Aeroforces   = Array<OneD,NekDouble> (2, 0.0);
    }
    
    // What are this bi-dimensional vectors ------------------------------------------
    // m_zta[0] = m_zta                     |  m_eta[0] = m_eta                      |
    // m_zta[1] = d(m_zta)/dt               |  m_eta[1] = d(m_eta)/dt                |
    // m_zta[2] = dd(m_zta)/ddtt            |  m_eta[2] = dd(m_eta)/ddtt             |
    //--------------------------------------------------------------------------------
    

    for(int i = 0; i < m_zta.size(); i++)
    {
        m_zta[i] = Array<OneD, NekDouble>(phystot,0.0);
        m_eta[i] = Array<OneD, NekDouble>(phystot,0.0);
    }
    forces_vel = Array<OneD, Array<OneD, NekDouble> > (3);

    for (int i = 0; i < 3; ++i)
    {
        forces_vel[i] = Array<OneD, NekDouble>(1);
    }
    comp = 0;
   

}

void ForcingMovingBody::v_Apply(
        const Array<OneD, MultiRegions::ExpListSharedPtr>&  pFields,
        const Array<OneD, Array<OneD, NekDouble> >&         inarray,
              Array<OneD, Array<OneD, NekDouble> >&         outarray,
        const NekDouble&                                    time)
{
    
    NekDouble aux0, aux1, alpha;
    std::string evol_operator = m_session->GetSolverInfo("EvolutionOperator");
    SolverUtils::DriverArnoldi::ReturnStructVector(aux0, aux1, c);
    SolverUtils::DriverArnoldi::Getalpha(alpha);
    std::string driver = m_session->GetSolverInfo("Driver");
    int       tot_step = m_session->GetParameter("NumSteps");

    int num_step = time/m_timestep;
    if(time==0)
    {
        if (comp==0)
        {
            m_MotionVars[1][0] = 0.;
            m_MotionVars[1][1] = 0.;
            
            cout << "displ: " << m_MotionVars[1][0] << ", vel: " << m_MotionVars[1][1] << endl;
        }
        else if(comp>0 && c==0){
            
            m_MotionVars[1][0] = aux0/m_structrho;
            m_MotionVars[1][1] = aux1/m_structrho;
            
                   
        }
        else if(comp>0 && c==1){
            m_MotionVars[1][0] = y_dagger;
            m_MotionVars[1][1] = y1_dagger;
            
            cout << "displ adj: " << m_MotionVars[1][0] << ", vel adj: " << m_MotionVars[1][1] << endl;
        
        }
        else{
            cout << "get dir... " << "displ: " << m_MotionVars[1][0] << ", vel: " << m_MotionVars[1][1] << endl;

        }
        
    }
      
    // Update the forces from the calculation of fluid field, which is
    // implemented in the movingbody filter
    
    std::string ode_solver    = m_session->GetSolverInfo("ODESolver");
    // Array<OneD, NekDouble> Hydroforces (2*m_np,0.0);
    // SolverUtils::FilterAeroForces::GetTotForces(Hydroforces);
    if(boost::iequals(ode_solver, "Newmark") && (boost::iequals(evol_operator, "Adjoint") || (boost::iequals(evol_operator, "TransientGrowth") && c==1)))
    {
        int cn;
        NekDouble dif, dif1, integral;
        m_MovBodyfilter->UpdateForce(m_session, pFields, m_Aeroforces, time);
        for(int n = 0, cn = 1; n < m_vdim; n++, cn--)
        {           
            RollOver(forces_vel);

            forces_vel[0][0] =  m_Aeroforces[cn];

            dif = (3*forces_vel[0][0] - 4*forces_vel[1][0] + forces_vel[2][0])/(2.*m_timestep);
            m_Aeroforces[cn] = - dif + m_Aeroforces[4+cn];

            
        }
        
    }
    else
    {
        
        m_MovBodyfilter->UpdateForce(m_session, pFields, m_Aeroforces, time);
        
    }
   
    
    std::string solver_type = m_session->GetSolverInfo("SolverType");
    // for "free" type (m_vdim = 2), the cable vibrates both in streamwise and crossflow
    // dimections, for "constrained" type (m_vdim = 1), the cable only vibrates in crossflow
    // dimection, and for "forced" type (m_vdim = 0), the calbe vibrates specifically along
    // a given function or file
    if(m_vdim == 1 || m_vdim == 2)
    {

         // For free vibration case, displacements, velocities and acceleartions
         // are obtained through solving structure dynamic model
        
        EvaluateStructDynModel(pFields, m_Aeroforces, time);
        
        
        if(boost::iequals(solver_type, "VCSMapping"))
        {
            // Convert result to format required by mapping
            int physTot = pFields[0]->GetTotPoints();
            Array< OneD, Array< OneD, NekDouble> >  coords(3);
            Array< OneD, Array< OneD, NekDouble> >  coordsVel(3);
            for(int i =0; i<3; i++)
            {
                coords[i] = Array< OneD, NekDouble> (physTot, 0.0);
                coordsVel[i] = Array< OneD, NekDouble> (physTot, 0.0);
            }
            // Get original coordinates
            pFields[0]->GetCoords(coords[0], coords[1], coords[2]);

            // Add displacement to coordinates
            NekDouble aux = m_MotionVars[0][0];
            Vmath::Sadd(physTot, aux, coords[0], 1, coords[0], 1);
            aux = m_MotionVars[1][0];
            Vmath::Sadd(physTot, aux, coords[1], 1, coords[1], 1);
            
            // fill coordsVel
            aux = m_MotionVars[0][1];
            Vmath::Fill(physTot, aux, coordsVel[0], 1);
            aux = m_MotionVars[1][1];
            Vmath::Fill(physTot, aux, coordsVel[1], 1);

            // Update mapping
            m_mapping->UpdateMapping(time, coords, coordsVel);
        }

    }
    else if(m_vdim == 0)
    {
        // For forced vibration case, load from specified file or function
        int cnt = 0;
        for(int j = 0; j < m_funcName.size(); j++)
        {
            if(m_IsFromFile[cnt] && m_IsFromFile[cnt+1])
            {
                ASSERTL0(false, "Motion loading from file needs specific "
                                "implementation: Work in Progress!");
            }
            else
            {
                GetFunction(pFields, m_session, m_funcName[j], true)->Evaluate(m_motion[0], m_zta[j], time);
                GetFunction(pFields, m_session, m_funcName[j], true)->Evaluate(m_motion[1], m_eta[j], time);
                cnt = cnt + 2;
            }
        }

        // Update mapping
        m_mapping->UpdateMapping(time);

        // Convert result from mapping
        int physTot = pFields[0]->GetTotPoints();
        Array< OneD, Array< OneD, NekDouble> >  coords(3);
        Array< OneD, Array< OneD, NekDouble> >  coordsVel(3);
        for(int i =0; i<3; i++)
        {
            coords[i] = Array< OneD, NekDouble> (physTot, 0.0);
            coordsVel[i] = Array< OneD, NekDouble> (physTot, 0.0);
        }
        // Get original coordinates
        pFields[0]->GetCoords(coords[0], coords[1], coords[2]);

        // Get Coordinates and coord velocity from mapping
        m_mapping->GetCartesianCoordinates(m_zta[0], m_eta[0], coords[2]);
        m_mapping->GetCoordVelocity(coordsVel);

        // Calculate displacement
        Vmath::Vsub(physTot, m_zta[0], 1, coords[0], 1, m_zta[0], 1);
        Vmath::Vsub(physTot, m_eta[0], 1, coords[1], 1, m_eta[0], 1);

        Vmath::Vcopy(physTot, coordsVel[0], 1, m_zta[1], 1);
        Vmath::Vcopy(physTot, coordsVel[1], 1, m_eta[1], 1);

        for(int var = 0; var < 3; var++)
        {
            for(int plane = 0; plane < m_np; plane++)
            {
                int n = pFields[0]->GetPlane(plane)->GetTotPoints();
                int offset  = plane * n;
                int Offset = var * m_np+plane;

                m_MotionVars[0][Offset] = m_zta[var][offset];
                m_MotionVars[1][Offset] = m_eta[var][offset];
            }
        }
    }
    else
    {
        ASSERTL0(false,
                 "Unrecogized vibration type for cable's dynamic model");
    }
    
    LibUtilities::CommSharedPtr vcomm = pFields[0]->GetComm();
    int colrank = vcomm->GetColumnComm()->GetRank();
    // Pass the variables of the cable's motion to the movingbody filter
    if(colrank == 0)
    {
        int n = m_MotionVars[0].size();
        Array<OneD, NekDouble> tmpArray(2*n),tmp(n);
        Vmath::Vcopy(n,m_MotionVars[0],1,tmpArray,1);
        Vmath::Vcopy(n,m_MotionVars[1],1,tmp=tmpArray+n,1);
        m_MovBodyfilter->UpdateMotion(m_session, pFields, tmpArray, time);
    }
    
    
    if(boost::iequals(solver_type, "VelocityCorrectionScheme"))
    {
       
        MappingBndConditions(pFields, inarray, time);
        
        std::string evol_operator = m_session->GetSolverInfo("EvolutionOperator");
        if(boost::iequals(evol_operator, "Nonlinear"))
        {
            Vmath::Sadd(outarray[0].size(), -m_MotionVars[0][2], outarray[0], 1, outarray[0], 1);
            Vmath::Sadd(outarray[1].size(), -m_MotionVars[1][2], outarray[1], 1, outarray[1], 1);

        }

        
    }
    if((boost::iequals(evol_operator, "TransientGrowth") && c==0))
    {
        // cout << "da0" << endl;
        y_dagger  = m_MotionVars[1][0];
        y1_dagger = m_MotionVars[1][1];
        // SolverUtils::DriverArnoldi::GetStructVector(m_MotionVars[1][0], m_MotionVars[1][1]*m_structrho);
    }
    if((boost::iequals(evol_operator, "TransientGrowth") && c==1))
    {
        // cout << "da0" << endl;
        y  = m_MotionVars[1][0];
        y1 = m_MotionVars[1][1];
        // SolverUtils::DriverArnoldi::GetStructVector(m_MotionVars[1][0], m_MotionVars[1][1]*m_structrho);
    }
    // else
    // {
    SolverUtils::DriverArnoldi::GetStructVector(m_MotionVars[1][0], m_MotionVars[1][1]);
        // cout << "da1" << endl;
    // }
}



/**
 *
 */
 void ForcingMovingBody::EvaluateStructDynModel(
         const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
               Array<OneD, NekDouble> &Hydroforces,
               NekDouble  time)
 {
    LibUtilities::CommSharedPtr vcomm = pFields[0]->GetComm();
    int colrank = vcomm->GetColumnComm()->GetRank();
    int nproc   = vcomm->GetColumnComm()->GetSize();

    bool homostrip;
    bool m_isHomogeneous1D;
    // For 3DH1D
    m_session->MatchSolverInfo("Homogeneous", "1D",
                               m_isHomogeneous1D, false);
                               
    m_session->MatchSolverInfo("HomoStrip","True",homostrip,false);
    std::string evol_operator = m_session->GetSolverInfo("EvolutionOperator");
    std::string solver_type = m_session->GetSolverInfo("SolverType");
    

     //number of structral modes and number of strips
     int npts, nstrips;
    if(boost::iequals(solver_type, "VCSMapping") && m_isHomogeneous1D)
    {
        if(!homostrip) //full resolutions
        {
            npts = m_session->GetParameter("HomModesZ");
        }
        else
        {
            m_session->LoadParameter("HomStructModesZ", npts);
            m_session->LoadParameter("Strip_Z", nstrips);

        }
    }
    else
    {
        npts = 1;
    }


    
    //the hydrodynamic forces
    Array<OneD, Array <OneD, NekDouble> > fces(2);
    //forces in x-direction
    fces[0] = Array <OneD, NekDouble> (npts,0.0);
    //forces in y-direction
    fces[1] = Array <OneD, NekDouble> (npts,0.0);

     //fill the force vectors
    if(colrank == 0)
    {
    
        if(!homostrip) //full resolutions
        {
            Vmath::Vcopy(m_np, Hydroforces,      1, fces[0], 1);
            Vmath::Vcopy(m_np, Hydroforces+m_np, 1, fces[1], 1);
            
        }
        else //strip modelling
        {
            fces[0][0] = Hydroforces[0];
            fces[1][0] = Hydroforces[m_np];
             
        }

        if(boost::iequals(solver_type, "VCSMapping") && m_isHomogeneous1D)
        {
            if(!homostrip) //full resolutions
            {
                Array<OneD, NekDouble> tmp(2*m_np);
                for (int i = 1; i < nproc; ++i)
                {
                    vcomm->GetColumnComm()->Recv(i, tmp);
                    for(int n = 0; n < m_np; n++)
                    {
                        for(int j = 0; j < 2; ++j)
                        {
                            fces[j][i*m_np + n] = tmp[j*m_np + n];
                        }
                    }
                }
            }
            else //strip modelling
            //if the body is submerged partly, the fces are filled partly
            //by the flow induced forces
            {
                Array<OneD, NekDouble> tmp(2);
                for(int i = 1; i < nstrips; ++i)
                {
                    vcomm->GetColumnComm()->Recv(i, tmp);

                    for(int j = 0 ; j < 2; ++j)
                    {
                        fces[j][i] = tmp[j];
                    }
                }
            }
        }
    }
    else
    {
        if(!homostrip) //full resolutions
        {
            vcomm->GetColumnComm()->Send(0, Hydroforces);
        }
        else //strip modelling
        {
            for(int i = 1; i < nstrips; ++i)
            {
                if(colrank == i)
                {
                    Array<OneD, NekDouble> tmp(2);
                    tmp[0] = Hydroforces[0];
                    tmp[1] = Hydroforces[m_np];
                    vcomm->GetColumnComm()->Send(0, tmp);
                }
            }
        }
    }

    if(colrank == 0)
    {
        // Fictitious mass method used to stablize the explicit coupling at
        // relatively lower mass ratio
        bool fictmass;
        m_session->MatchSolverInfo("FictitiousMassMethod", "True",
                                    fictmass, false);
        if(fictmass)
        {
            NekDouble fictrho, fictdamp;
            m_session->LoadParameter("FictMass", fictrho);
            m_session->LoadParameter("FictDamp", fictdamp);
        
            static NekDouble Betaq_Coeffs[2][2] =
                                {{1.0,  0.0},{2.0, -1.0}};
        
            // only consider second order approximation for fictitious variables
            int  intOrder= 2;
            int  nint    = min(m_movingBodyCalls+1,intOrder);
            int  nlevels = m_fV[0].size();
        
            for(int i = 0; i < m_motion.size(); ++i)
            {
                RollOver(m_fV[i]);
                RollOver(m_fA[i]);
        
                int Voffset = npts;
                int Aoffset = 2*npts;
        
                Vmath::Vcopy(npts, m_MotionVars[i]+Voffset, 1, m_fV[i][0], 1);
                Vmath::Vcopy(npts, m_MotionVars[i]+Aoffset, 1, m_fA[i][0], 1);
        
                // Extrapolate to n+1
                Vmath::Smul(npts,
                            Betaq_Coeffs[nint-1][nint-1],
                            m_fV[i][nint-1],    1,
                            m_fV[i][nlevels-1], 1);
                Vmath::Smul(npts,
                            Betaq_Coeffs[nint-1][nint-1],
                            m_fA[i][nint-1],    1,
                            m_fA[i][nlevels-1], 1);
        
                for(int n = 0; n < nint-1; ++n)
                {
                    Vmath::Svtvp(npts,
                                Betaq_Coeffs[nint-1][n],
                                m_fV[i][n],1,m_fV[i][nlevels-1],1,
                                m_fV[i][nlevels-1],1);
                    Vmath::Svtvp(npts,
                                Betaq_Coeffs[nint-1][n],
                                m_fA[i][n],1,m_fA[i][nlevels-1],1,
                                m_fA[i][nlevels-1],1);
                }
        
                // Add the fictitious forces on the RHS of the equation
                Vmath::Svtvp(npts, fictdamp,m_fV[i][nlevels-1],1,
                            fces[i],1,fces[i],1);
                Vmath::Svtvp(npts, fictrho, m_fA[i][nlevels-1],1,
                            fces[i],1,fces[i],1);
            }
        }
        
    }

    //structural solver is implemented on the root process
    if(colrank == 0)
    {
        
        //Tensioned cable model is evaluated in wave space
        for(int n = 0, cn = 1; n < m_vdim; n++, cn--)
        {   std::string evol_operator = m_session->GetSolverInfo("EvolutionOperator");
      
            if (comp==0. && (boost::iequals(evol_operator, "Adjoint") & boost::iequals(evol_operator, "Direct")))
            {   
                m_structstiff = m_structstiff - m_Aeroforces[2 + cn];
                SetDynEqCoeffMatrix(pFields,cn);
               
            }
            ++comp;
                           
            StructureSolver(pFields, cn, fces[cn], m_MotionVars[cn]);

         
        }

    }


    if(boost::iequals(solver_type, "VCSMapping") && m_isHomogeneous1D)
    {
        Array<OneD, NekDouble> Motvars(2*2*m_np);
        // send physical coeffients to all planes of each processor
        if(!homostrip)//full resolutions
        {
            Array<OneD, NekDouble> tmp(2*2*m_np);

            if(colrank != 0)
            {
                vcomm->GetColumnComm()->Recv(0, tmp);
                Vmath::Vcopy(2*2*m_np, tmp, 1, Motvars, 1);
            }
            else
            {
                for (int i = 1; i < nproc; ++i)
                {
                    for(int j = 0; j < 2; j++) //moving dimensions
                    {
                        for(int k = 0; k < 2; k++) //disp. and vel.
                        {
                            for (int n = 0; n < m_np; n++)
                            {
                                tmp[j*2*m_np+k*m_np+n] = m_MotionVars[j][k*npts+i*m_np+n];
                            }
                        }
                    }
                    vcomm->GetColumnComm()->Send(i, tmp);
                }

                for(int j = 0; j < 2; j++)
                {
                    for(int k = 0; k < 2; k++)
                    {
                        for(int n = 0; n < m_np; n++)
                        {
                            tmp[j*2*m_np+k*m_np+n] = m_MotionVars[j][k*npts+n];
                        }
                        Vmath::Vcopy(2*2*m_np, tmp, 1, Motvars, 1);
                    }
                }
            }
        }
        
        else //strip modelling
        {
            Array<OneD, NekDouble> tmp(2*2);

            if(colrank != 0)
            {
                for (int j = 1; j < nproc/nstrips; j++)
                {
                    if(colrank == j*nstrips)
                    {
                        vcomm->GetColumnComm()->Recv(0, tmp);

                        for(int plane = 0; plane < m_np; plane++)
                        {
                            for(int var = 0; var < 2; var++)
                            {
                                for(int k = 0; k < 2; k++)
                                {
                                    Motvars[var*2*m_np+k*m_np+plane]= tmp[var*2+k];
                                }
                            }
                        }
                    }
                }

                for(int i = 1; i < nstrips; i++)
                {
                    for (int j = 0; j < nproc/nstrips; j++)
                    {
                        if(colrank == i+j*nstrips)
                        {
                            vcomm->GetColumnComm()->Recv(0, tmp);

                            for(int plane = 0; plane < m_np; plane++)
                            {
                                for(int var = 0; var < 2; var++)
                                {
                                    for(int k = 0; k < 2; k++)
                                    {
                                        Motvars[var*2*m_np+k*m_np+plane] = tmp[var*2+k];
                                    }
                                }
                            }
                        }
                    }
                }
            }
            else
            {
                for(int j = 0; j < 2; ++j)
                {
                    for(int k = 0; k < 2; ++k)
                    {
                        tmp[j*2+k] = m_MotionVars[j][k*npts];
                    }
                }

                for (int j = 1; j < nproc/nstrips; j++)
                {
                    vcomm->GetColumnComm()->Send(j*nstrips, tmp);
                }

                for(int plane = 0; plane < m_np; plane++)
                {
                    for(int j = 0; j < 2; j++)
                    {
                        for(int k = 0; k < 2; ++k)
                        {
                            Motvars[j*2*m_np+k*m_np+plane] = m_MotionVars[j][k*npts];
                        }
                    }
                }

                for(int i = 1; i < nstrips; ++i)
                {
                    for(int j = 0; j < 2; ++j)
                    {
                        for(int k = 0; k < 2; ++k)
                        {
                            tmp[j*2+k] = m_MotionVars[j][i+k*npts];
                        }
                    }

                    for (int j = 0; j < nproc/nstrips; j++)
                    {
                        vcomm->GetColumnComm()->Send(i+j*nstrips, tmp);
                    }
                }
            }
        }

        // Set the m_forcing term based on the motion of the cable
        for(int var = 0; var < 2; var++)
        {
            for(int plane = 0; plane < m_np; plane++)
            {
                int n = pFields[0]->GetPlane(plane)->GetTotPoints();

                Array<OneD, NekDouble> tmp;

                int offset  = plane * n;
                int xoffset = var * m_np+plane;
                int yoffset = 2*m_np + xoffset;

                Vmath::Fill(n, Motvars[xoffset], tmp = m_zta[var] + offset, 1);
                Vmath::Fill(n, Motvars[yoffset], tmp = m_eta[var] + offset, 1);
            }
        }
     }
    


 }


/**
 *
 */
void ForcingMovingBody::StructureSolver(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
              int cn, Array<OneD, NekDouble> &HydroForces,
              Array<OneD, NekDouble> &BodyMotions)
{
    std::string evol_operator = m_session->GetSolverInfo("EvolutionOperator");
    std::string ode_solver    = m_session->GetSolverInfo("ODESolver");
    int m_subSteps = 1;
    static int intCalls = 0;
    ++intCalls;
    int order = min(intCalls, m_intSteps);

    // Array<OneD, NekDouble> Hydroforces (2*m_np,0.0);
    // SolverUtils::FilterAeroForces::GetTotForces(Hydroforces);
    if(boost::iequals(ode_solver, "Newmark"))
    {
    
        int nrows = 3;
        
        Array<OneD, NekDouble> tmp0,tmp1,tmp2;
        tmp0 = Array<OneD, NekDouble> (3,0.0);
        tmp1 = Array<OneD, NekDouble> (3,0.0);
        tmp2 = Array<OneD, NekDouble> (3,0.0);

        for(int var = 0; var < 3; var++)
        {
            tmp0[var] = BodyMotions[var];
        }
    
        tmp2[0] = m_Aeroforces[cn];
        
        Blas::Dgemv('N', nrows, nrows, 1.0,
                    &(m_CoeffMat_B[0]->GetPtr())[0],
                    nrows, &tmp0[0], 1,
                    0.0,   &tmp1[0], 1);
        Blas::Dgemv('N', nrows, nrows, 1.0/m_structrho,
                    &(m_CoeffMat_A[0]->GetPtr())[0],
                    nrows, &tmp2[0], 1,
                    1.0,   &tmp1[0], 1);
        if(boost::iequals(evol_operator, "Adjoint") || (boost::iequals(evol_operator, "TransientGrowth") && c==1))
        {
            m_displacement[0] = m_Aeroforces[cn] - m_structrho*tmp1[1];
            if(time==0){
                cout << "adj displ. " << m_displacement[0] << endl;
            }
            
            // BodyMotions[0] = m_displacement[0];
            // BodyMotions[1] = tmp1[1];
        }
        // else
        // {
        for(int i = 0; i < BodyMotions.size(); i++)
        {
            BodyMotions[i] = tmp1[i];
        }
        // }
        // std::string solver_type = m_session->GetSolverInfo("SolverType");
    
    }
    if(boost::iequals(ode_solver, "Adams"))
    {

        // Get correct time considering m_subSteps parameter
        // NekDouble   newTime = time;


        static int totalCalls = -1;
        ++totalCalls;
        if( !(totalCalls % m_subSteps) )
        {

            int expdim = pFields[0]->GetGraph()->GetMeshDimension();

            // Rotate force storage
            Array<OneD, NekDouble> tmp = m_force[m_intSteps-1];
            for(int n = m_intSteps-1; n > 0; --n)
            {
                m_force[n] = m_force[n-1];
            }
            m_force[0] = tmp;

            // Calculate total force

            m_displacement[0] = BodyMotions[0]; //1dof

            m_velocity[0][0]  = BodyMotions[1];
            if(boost::iequals(evol_operator, "Adjoint") || (boost::iequals(evol_operator, "TransientGrowth") && c==1))
            {
                m_force[0][0] =  -m_Aeroforces[cn] -
                        m_displacement[0] + m_structdamp * m_velocity[0][0]/m_structrho;
            }
            else
            {
                m_force[0][0] = (m_Aeroforces[cn] -
                    m_structstiff * m_displacement[0] - m_structdamp * m_velocity[0][0])/m_structrho;

            }

            

            // Rotate velocity storage, keeping value of velocity[0]
            Vmath::Vcopy(m_vdim, m_velocity[0], 1, m_velocity[m_intSteps-1], 1);
            tmp = m_velocity[m_intSteps-1];
            for(int n = m_intSteps-1; n > 0; --n)
            {
                m_velocity[n] = m_velocity[n-1];
            }
            m_velocity[0] = tmp;
        
            // Update velocity
            for(int j = 0; j < order; ++j)
            {
                m_velocity[0][0] += m_subSteps * m_timestep *
                    AdamsBashforth_coeffs[order-1][j] * m_force[j][0];
            }
            // Update position
            tmp = m_force1[m_intSteps-1];
            for(int n = m_intSteps-1; n > 0; --n)
            {
                m_force1[n] = m_force1[n-1];
            }
            m_force1[0] = tmp;

            if(boost::iequals(evol_operator, "Adjoint") || (boost::iequals(evol_operator, "TransientGrowth") && c==1))
            {
                m_force1[0][0] = -m_Aeroforces[4+cn] +
                        m_velocity[0][0]*m_structstiff/m_structrho;
            }
            for(int j = 0; j < order; ++j)
            {
                if(boost::iequals(evol_operator, "Adjoint") || (boost::iequals(evol_operator, "TransientGrowth") && c==1))
                {
                    m_displacement[0] += m_subSteps * m_timestep *
                        AdamsMoulton_coeffs[order-1][j] * m_force1[j][0];
                }
                else
                {
                    m_displacement[0] += m_subSteps * m_timestep *
                    AdamsMoulton_coeffs[order-1][j] * m_velocity[j][0];
                }
            }
        
            BodyMotions[0] = m_displacement[0];
            BodyMotions[1] = m_velocity[0][0];
            BodyMotions[2] = (m_Aeroforces[cn] - m_structstiff * BodyMotions[0] - m_structdamp * BodyMotions[1])/m_structrho;

    }

    }

}

/**
 *
 */
void ForcingMovingBody::InitialiseCableModel(
        const LibUtilities::SessionReaderSharedPtr& pSession,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields)
{
    m_movingBodyCalls = 0;
    m_session->LoadParameter("Kinvis",m_kinvis);
    m_session->LoadParameter("TimeStep", m_timestep, 0.01);

    LibUtilities::CommSharedPtr vcomm = pFields[0]->GetComm();
    int colrank = vcomm->GetColumnComm()->GetRank();
    int nproc   = vcomm->GetColumnComm()->GetSize();
    Array<OneD, unsigned int> ZIDs;

    //number of structral modes
    int npts, nplanes;
    bool homostrip, m_isHomogeneous1D;
    // For 3DH1D
    m_session->MatchSolverInfo("Homogeneous", "1D",
                               m_isHomogeneous1D, false);
    m_session->MatchSolverInfo("HomoStrip","True",homostrip,false);
    std::string evol_operator = m_session->GetSolverInfo("EvolutionOperator");
    std::string solver_type = m_session->GetSolverInfo("SolverType");
    if(boost::iequals(solver_type, "VCSMapping") && m_isHomogeneous1D)
    {
       if(!homostrip) //full resolutions
       {
            npts = m_session->GetParameter("HomModesZ");
       }
       else
       {
            m_session->LoadParameter("HomStructModesZ", npts);
       }

    }
    else 
    {
        npts = 1;
    }
    

    m_MotionVars = Array<OneD, Array<OneD, NekDouble> > (2);
    m_MotionVars[0] = Array<OneD, NekDouble>(3*npts,0.0);
    m_MotionVars[1] = Array<OneD, NekDouble>(3*npts,0.0);

    if(boost::iequals(solver_type, "VCSMapping") && m_isHomogeneous1D)
    {

        ZIDs = pFields[0]->GetZIDs();
        m_np = ZIDs.size();
    }
    else 
    {
        m_np = 1;
    }


    std::string vibtype = m_session->GetSolverInfo("VibrationType");

    if(boost::iequals(vibtype, "Constrained"))
    {
        m_vdim = 1;
    }
    else if (boost::iequals(vibtype, "Free"))
    {
        m_vdim = 2;
    }
    else if (boost::iequals(vibtype, "Forced"))
    {
        m_vdim = 0;
        return;
    }

    if(boost::iequals(solver_type, "VCSMapping") && m_isHomogeneous1D)
    {
        
          if(!homostrip)
          {
              m_session->LoadParameter("LZ", m_lhom);
              nplanes = m_session->GetParameter("HomModesZ");
              m_FFT =
                  LibUtilities::GetNektarFFTFactory().CreateInstance(
                                                  "NekFFTW", nplanes);
          }
          else
          {
              int nstrips;
              NekDouble DistStrip;

              m_session->LoadParameter("DistStrip", DistStrip);
              m_session->LoadParameter("Strip_Z", nstrips);
              m_lhom = nstrips * DistStrip;
              m_FFT =
                  LibUtilities::GetNektarFFTFactory().CreateInstance(
                                                  "NekFFTW", nstrips);
          }
    }
    else 
    {
        nplanes = 1;
    }

    // load the structural dynamic parameters from xml file
    m_session->LoadParameter("StructRho",  m_structrho);
    m_session->LoadParameter("StructStiff", m_structstiff,  0.0);
    m_session->LoadParameter("StructDamp", m_structdamp, 0.0);

    // Identify whether the fictitious mass method is active for explicit
    // coupling of fluid solver and structural dynamics solver
    bool fictmass;
    m_session->MatchSolverInfo("FictitiousMassMethod", "True",
                                fictmass, false);
    if(fictmass)
    {
        NekDouble fictrho, fictdamp;
        m_session->LoadParameter("FictMass", fictrho);
        m_session->LoadParameter("FictDamp", fictdamp);
        m_structrho  += fictrho;
        m_structdamp += fictdamp;

        // Storage array of Struct Velocity and Acceleration used for
        // extrapolation of fictitious force
        m_fV = Array<OneD, Array<OneD, Array<OneD, NekDouble> > > (2);
        m_fA = Array<OneD, Array<OneD, Array<OneD, NekDouble> > > (2);
        for (int i = 0; i < m_motion.size(); ++i)
        {
            m_fV[i] = Array<OneD, Array<OneD, NekDouble> > (2);
            m_fA[i] = Array<OneD, Array<OneD, NekDouble> > (2);

            for(int n = 0; n < 2; ++n)
            {
                m_fV[i][n] = Array<OneD, NekDouble>(npts, 0.0);
                m_fA[i][n] = Array<OneD, NekDouble>(npts, 0.0);
            }
        }
    }
    //Tensioned cable model is evaluated in wave space
    for(int n = 0, cn = 1; n < m_vdim; n++, cn--)
    {
        SetDynEqCoeffMatrix(pFields,cn);
    }
    // Setting the coefficient matrices for solving structural dynamic ODEs
 
    
    // Set initial condition for cable's motion
    int cnt = 0;
    
    for(int j = 0; j < m_funcName.size(); j++)
    {
        // loading from the specified files through inputstream
        if(m_IsFromFile[cnt] && m_IsFromFile[cnt+1])
        {
            std::ifstream inputStream;
            int nzpoints = 1;


                if(homostrip)
                {
                    m_session->LoadParameter("HomStructModesZ", nzpoints);
                }
                else
                {
                    nzpoints = pFields[0]->GetHomogeneousBasis()->GetNumModes();
                }



            if (vcomm->GetRank() == 0)
            {
                std::string filename = m_session->GetFunctionFilename(
                    m_funcName[j], m_motion[0]);

                // Open intputstream for cable motions
                inputStream.open(filename.c_str());

                // Import the head string from the file
                Array<OneD, std::string> tmp(9);
                for(int n = 0; n < tmp.size(); n++)
                {
                    inputStream >> tmp[n];
                }

                NekDouble time, z_cds;
                // Import the motion variables from the file
                for (int n = 0; n < nzpoints; n++)
                {
                    inputStream >> setprecision(6) >> time;
                    inputStream >> setprecision(6) >> z_cds;
                    inputStream >> setprecision(8) >> m_MotionVars[0][n];
                    inputStream >> setprecision(8) >> m_MotionVars[0][n+nzpoints];
                    inputStream >> setprecision(8) >> m_MotionVars[0][n+2*nzpoints];
                    inputStream >> setprecision(8) >> m_MotionVars[1][n];
                    inputStream >> setprecision(8) >> m_MotionVars[1][n+nzpoints];
                    inputStream >> setprecision(8) >> m_MotionVars[1][n+2*nzpoints];
                }
                // Close inputstream for cable motions
                inputStream.close();
            }
            cnt = cnt + 2;
        }
        else //Evaluate from the functions specified in xml file
        {   
            if(boost::iequals(solver_type, "VCSMapping") && m_isHomogeneous1D)
            {
                if(!homostrip)
                {
                    int nzpoints = pFields[0]->GetHomogeneousBasis()->GetNumModes();
                    Array<OneD, NekDouble> z_coords(nzpoints,0.0);
                    Array<OneD, const NekDouble> pts
                                        = pFields[0]->GetHomogeneousBasis()->GetZ();

                    Vmath::Smul(nzpoints,m_lhom/2.0,pts,1,z_coords,1);
                    Vmath::Sadd(nzpoints,m_lhom/2.0,z_coords,1,z_coords,1);

                    Array<OneD, NekDouble> x0(m_np, 0.0);
                    Array<OneD, NekDouble> x1(m_np, 0.0);
                    Array<OneD, NekDouble> x2(m_np, 0.0);
                    for (int plane = 0; plane < m_np; plane++)
                    {
                        x2[plane] = z_coords[ZIDs[plane]];
                    }

                    Array<OneD, NekDouble> tmp (m_np,0.0);
                    Array<OneD, NekDouble> tmp0(m_np,0.0);
                    Array<OneD, NekDouble> tmp1(m_np,0.0);
                    LibUtilities::EquationSharedPtr ffunc0,ffunc1;

                    ffunc0 = m_session->GetFunction(m_funcName[j], m_motion[0]);
                    ffunc1 = m_session->GetFunction(m_funcName[j], m_motion[1]);
                    ffunc0->Evaluate(x0, x1, x2, 0.0, tmp0);
                    ffunc1->Evaluate(x0, x1, x2, 0.0, tmp1);

                    int offset = j*npts;
                    Vmath::Vcopy(m_np, tmp0, 1, tmp = m_MotionVars[0]+offset,1);
                    Vmath::Vcopy(m_np, tmp1, 1, tmp = m_MotionVars[1]+offset,1);

                    if(colrank == 0)
                    {
                        for (int i = 1; i < nproc; ++i)
                        {
                            vcomm->GetColumnComm()->Recv(i, tmp0);
                            vcomm->GetColumnComm()->Recv(i, tmp1);
                            Vmath::Vcopy(m_np, tmp0, 1, tmp = m_MotionVars[0]+offset+i*m_np,1);
                            Vmath::Vcopy(m_np, tmp1, 1, tmp = m_MotionVars[1]+offset+i*m_np,1);
                        }
                    }
                    else
                    {
                        vcomm->GetColumnComm()->Send(0, tmp0);
                        vcomm->GetColumnComm()->Send(0, tmp1);
                    }
                }
            }    
            else
            {
                
                if(colrank == 0)
                {
                    int nstrips;
                    Array<OneD, NekDouble> x0(npts, 0.0);
                    Array<OneD, NekDouble> x1(npts, 0.0);
                    Array<OneD, NekDouble> x2(npts, 0.0);
                    Array<OneD, NekDouble> tmp (npts,0.0);
                    Array<OneD, NekDouble> tmp0(npts,0.0);
                    Array<OneD, NekDouble> tmp1(npts,0.0);
                    
                    if(boost::iequals(solver_type, "VCSMapping") && m_isHomogeneous1D)
                    {
                        m_session->LoadParameter("Strip_Z", nstrips);

                        ASSERTL0(m_session->DefinesSolverInfo("USEFFT"),
                                "Fourier transformation of cable motion is currently "
                                "implemented only for FFTW module.");

                        NekDouble DistStrip;
                        m_session->LoadParameter("DistStrip", DistStrip);
                        for (int plane = 0; plane < npts; plane++)
                        {
                            x2[plane] = plane*DistStrip;
                        }
                    }
                    
                    
                
                    LibUtilities::EquationSharedPtr ffunc0,ffunc1;
                    ffunc0 = m_session->GetFunction(m_funcName[j], m_motion[0]);
                    ffunc1 = m_session->GetFunction(m_funcName[j], m_motion[1]);
                    ffunc0->Evaluate(x0, x1, x2, 0.0, tmp0);
                    ffunc1->Evaluate(x0, x1, x2, 0.0, tmp1);

                    int offset = j*npts;
                    Vmath::Vcopy(npts, tmp0, 1, tmp = m_MotionVars[0]+offset,1);
                    Vmath::Vcopy(npts, tmp1, 1, tmp = m_MotionVars[1]+offset,1);
                }
            }

            cnt = cnt + 2;
        }
    }


    
    // m_MotionVars[1][0] = 0.001;

}


/**
 *
 */
void ForcingMovingBody::SetDynEqCoeffMatrix(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields, int cn)
{
    int nplanes = 1;
    int nel = 3;

    NekDouble tmp1, tmp2, tmp3;
    NekDouble tmp4, tmp5, tmp6, tmp7;

    // load the structural dynamic parameters from xml file
    NekDouble cabletension;
    NekDouble bendingstiff;
  
    m_session->LoadParameter("CableTension", cabletension, 0.0);
    m_session->LoadParameter("BendingStiff", bendingstiff, 0.0);
    
    tmp1 =   m_timestep * m_timestep;
    tmp2 =  m_structstiff / m_structrho;
    tmp3 = m_structdamp / m_structrho;
    tmp4 = cabletension / m_structrho;
    tmp5 = bendingstiff / m_structrho;
    
    m_CoeffMat_A = Array<OneD, DNekMatSharedPtr>(nplanes);
    m_CoeffMat_B = Array<OneD, DNekMatSharedPtr>(nplanes);

    m_CoeffMat_A[0]
            = MemoryManager<DNekMat>::AllocateSharedPtr(nel,nel);
    m_CoeffMat_B[0]
            = MemoryManager<DNekMat>::AllocateSharedPtr(nel,nel);

    (*m_CoeffMat_A[0])(0,0) = tmp2;
    (*m_CoeffMat_A[0])(0,1) = tmp3;
    (*m_CoeffMat_A[0])(0,2) = 1.0;
    (*m_CoeffMat_A[0])(1,0) = 0.0;
    (*m_CoeffMat_A[0])(1,1) = 1.0;
    (*m_CoeffMat_A[0])(1,2) =-m_timestep/2.0;
    (*m_CoeffMat_A[0])(2,0) = 1.0;
    (*m_CoeffMat_A[0])(2,1) = 0.0;
    (*m_CoeffMat_A[0])(2,2) =-tmp1/4.0;
    (*m_CoeffMat_B[0])(0,0) = 0.0;
    (*m_CoeffMat_B[0])(0,1) = 0.0;
    (*m_CoeffMat_B[0])(0,2) = 0.0;
    (*m_CoeffMat_B[0])(1,0) = 0.0;
    (*m_CoeffMat_B[0])(1,1) = 1.0;
    (*m_CoeffMat_B[0])(1,2) = m_timestep/2.0;
    (*m_CoeffMat_B[0])(2,0) = 1.0;
    (*m_CoeffMat_B[0])(2,1) = m_timestep;
    (*m_CoeffMat_B[0])(2,2) = tmp1/4.0;


    m_CoeffMat_A[0]->Invert();
    (*m_CoeffMat_B[0]) =
        (*m_CoeffMat_A[0]) * (*m_CoeffMat_B[0]);

}


 
void ForcingMovingBody::MappingBndConditions(
            const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
            const Array<OneD, Array<OneD, NekDouble> >        &fields,
                  NekDouble  time)
{


    std::string evol_operator = m_session->GetSolverInfo("EvolutionOperator");
    std::string driver = m_session->GetSolverInfo("Driver");
    // Declare variables
    Array<OneD, MultiRegions::ExpListSharedPtr> BndExp;
    Array<OneD, const SpatialDomains::BoundaryConditionShPtr> BndConds;
    m_baseflow = Array<OneD, Array<OneD, NekDouble> >(2);
    int nbnds    = pFields[0]->GetBndConditions().size();
    for (int n = 0; n < nbnds; ++n)
    {
        
        for ( int dim = 0; dim < m_motion.size(); dim++)
        {
            
            BndConds   = pFields[dim]->GetBndConditions();
            BndExp     = pFields[dim]->GetBndCondExpansions();
            // cout << BndExp[0] << endl;
            MultiRegions::ExpListSharedPtr                   locExpList;
            locExpList = BndExp[n];
                   
            if (BndConds[n]->GetUserDefined() =="MovingBody")
            {   
                
                // cout << BndExp[n]->UpdatePhys()[0] << endl; 
                
                // Number of points on this boundary
                int nPts = BndExp[n]->GetTotPoints();
                Array<OneD, NekDouble> x0(nPts, 0.0);
                Array<OneD, NekDouble> x1(nPts, 0.0);
                Array<OneD, NekDouble> x2(nPts, 0.0);
                Array<OneD, NekDouble> tmp(nPts,0.0);
                NekDouble x2_in =0.;
                // Homogeneous input case for x2.
                if (x2_in == NekConstants::kNekUnsetDouble)
                {
                    BndExp[n]->GetCoords(x0,x1,x2);
                }
                else
                {
                    BndExp[n]->GetCoords(x0, x1, x2);
                    Vmath::Fill(nPts, x2_in, x2, 1);
                }


                LibUtilities::Equation condition =
                std::static_pointer_cast<
                    SpatialDomains::DirichletBoundaryCondition>
                        (BndConds[n])->
                            m_dirichletCondition;
                condition.Evaluate(x0, x1, x2, time,
                                        BndExp[n]->UpdatePhys());

                if((boost::iequals(evol_operator, "Direct") || (boost::iequals(evol_operator, "TransientGrowth") && c==0)) 
                && (boost::iequals(driver, "ModifiedArnoldi") || boost::iequals(driver, "Arpack")))
                {
                
                    m_baseflow[dim] = Array<OneD, NekDouble>(nPts, 0.0);

                    string file = m_session->GetFunctionFilename("BaseFlow_BC", 0);
                    
                    int m_slices = 1;
                    ImportFldBase(file,locExpList, m_slices, dim);
                
                    Vmath::Vcopy(nPts, m_baseflow[dim], 1, tmp, 1);

                    Vmath::Smul(nPts, m_MotionVars[1][0], tmp, 1, tmp, 1);

                    Vmath::Neg(nPts, tmp, 1);

                    Vmath::Sadd(nPts, m_MotionVars[dim][1], tmp, 1, tmp, 1);

                    Vmath::Vadd(nPts, BndExp[n]->UpdatePhys(), 1,
                                                 tmp,                          1,
                                                 BndExp[n]->UpdatePhys(), 1);
                }
                if((boost::iequals(evol_operator, "Adjoint") || (boost::iequals(evol_operator, "TransientGrowth") && c==1)) 
                && (boost::iequals(driver, "ModifiedArnoldi") || boost::iequals(driver, "Arpack")))
                {
                    Vmath::Fill(nPts, m_MotionVars[dim][1]/m_structrho, tmp, 1);

                    Vmath::Vsub(nPts, BndExp[n]->UpdatePhys(), 1,
                                                 tmp,                          1,
                                                 BndExp[n]->UpdatePhys(), 1);
                }

                else
                {

                    Vmath::Fill(nPts, m_MotionVars[dim][1], tmp, 1);

                    Vmath::Vsub(nPts, BndExp[n]->UpdatePhys(), 1,
                                                 tmp,                          1,
                                                 BndExp[n]->UpdatePhys(), 1);
     

                }

                
                // Update coefficients at the boundary
                BndExp[n]->FwdTrans_BndConstrained(BndExp[n]->GetPhys(),
                                                BndExp[n]->UpdateCoeffs());
              
            }
            
        }

    }

}

void ForcingMovingBody::ImportFldBase(
std::string                                  pInfile,
MultiRegions::ExpListSharedPtr            &locExpList,
        int                                          pSlice, int f)
{
    std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef;
    std::vector<std::vector<NekDouble> >                 FieldData;

    int nqtot = m_baseflow[0].size();

    Array<OneD, NekDouble> tmp_coeff(locExpList->GetNcoeffs(), 0.0);
    Array<OneD, NekDouble> tmp_coeff1(locExpList->GetNcoeffs(), 0.0);

    LibUtilities::FieldIOSharedPtr fld = LibUtilities::FieldIO::CreateForFile(
                                                m_session, pInfile);

    fld->Import(pInfile, FieldDef, FieldData);

    int nSessionConvVar = 1;
    int nFileVar        = FieldDef[0]->m_fields.size();
    int nFileConvVar    = nFileVar - 1; // Ignore pressure


    int v = f*2 + 4;
    for(int i = 0; i < FieldDef.size(); ++i)
    {
        bool flag = FieldDef[i]->m_fields[0] ==
            m_session->GetVariable(0);

        ASSERTL0(flag, (std::string("Order of ") + pInfile
                        + std::string(" data and that defined in "
                        "the session file differs")).c_str());

        locExpList->ExtractDataToCoeffs(
                            FieldDef[i],
                            FieldData[i],
                            FieldDef[i]->m_fields[v],
                            tmp_coeff);

    }
 
    locExpList->BwdTrans(tmp_coeff, m_baseflow[f]);

}
/**
 * Function to roll time-level storages to the next step layout.
 * The stored data associated with the oldest time-level
 * (not required anymore) are moved to the top, where they will
 * be overwritten as the solution process progresses.
 */
void ForcingMovingBody::RollOver(Array<OneD, Array<OneD, NekDouble> > &input)
{
    int nlevels = input.size();

    Array<OneD, NekDouble> tmp;
    tmp = input[nlevels-1];

    for(int n = nlevels-1; n > 0; --n)
    {
        input[n] = input[n-1];
    }

    input[0] = tmp;
}



/**
 *
 */
void ForcingMovingBody::CheckIsFromFile(const TiXmlElement* pForce)
{

    m_funcName = Array<OneD, std::string> (3);
    m_motion = Array<OneD, std::string> (2);
    m_motion[0] = "x";
    m_motion[1] = "y";

    m_IsFromFile = Array<OneD, bool> (6);
    // Loading the x-dispalcement (m_zta) and the y-displacement (m_eta)
    // Those two variables are bith functions of z and t and the may come
    // from an equation (forced vibration) or from another solver which, given
    // the aerodynamic forces at the previous step, calculates the
    // displacements.

    //Get the body displacement: m_zta and m_eta
    const TiXmlElement* funcNameElmt_D
                    = pForce->FirstChildElement("DISPLACEMENTS");
    ASSERTL0(funcNameElmt_D,
             "MOVINGBODYFORCE tag has to specify a function name which "
             "prescribes the body displacement as d(z,t).");

    m_funcName[0] = funcNameElmt_D->GetText();
    ASSERTL0(m_session->DefinesFunction(m_funcName[0]),
             "Function '" + m_funcName[0] + "' not defined.");

    //Get the body velocity of movement: d(m_zta)/dt and d(m_eta)/dt
    const TiXmlElement* funcNameElmt_V
                    = pForce->FirstChildElement("VELOCITIES");
    ASSERTL0(funcNameElmt_D,
             "MOVINGBODYFORCE tag has to specify a function name which "
             "prescribes the body velocity of movement as v(z,t).");

    m_funcName[1] = funcNameElmt_V->GetText();
    ASSERTL0(m_session->DefinesFunction(m_funcName[1]),
             "Function '" + m_funcName[1] + "' not defined.");


    //Get the body acceleration: dd(m_zta)/ddt and dd(m_eta)/ddt
    const TiXmlElement* funcNameElmt_A
                    = pForce->FirstChildElement("ACCELERATIONS");
    ASSERTL0(funcNameElmt_A,
             "MOVINGBODYFORCE tag has to specify a function name which "
             "prescribes the body acceleration as a(z,t).");

    m_funcName[2] = funcNameElmt_A->GetText();
    ASSERTL0(m_session->DefinesFunction(m_funcName[2]),
             "Function '" + m_funcName[2] + "' not defined.");

    LibUtilities::FunctionType vType;

    // Check Displacement x
    vType = m_session->GetFunctionType(m_funcName[0],m_motion[0]);
    if(vType == LibUtilities::eFunctionTypeExpression)
    {
        m_IsFromFile[0] = false;
    }
    else if (vType == LibUtilities::eFunctionTypeFile)
    {
        m_IsFromFile[0] = true;
    }
    else
    {
        ASSERTL0(false, "The displacements in x must be specified via an "
                        "equation <E> or a file <F>");
    }

    // Check Displacement y
    vType = m_session->GetFunctionType(m_funcName[0],m_motion[1]);
    if(vType == LibUtilities::eFunctionTypeExpression)
    {
        m_IsFromFile[1] = false;
    }
    else if (vType == LibUtilities::eFunctionTypeFile)
    {
        m_IsFromFile[1] = true;
    }
    else
    {
        ASSERTL0(false, "The displacements in y must be specified via an "
                        "equation <E> or a file <F>");
    }

    // Check Velocity x
    vType = m_session->GetFunctionType(m_funcName[1],m_motion[0]);
    if (vType == LibUtilities::eFunctionTypeExpression)
    {
        m_IsFromFile[2] = false;
    }
    else if (vType == LibUtilities::eFunctionTypeFile)
    {
        m_IsFromFile[2] = true;
    }
    else
    {
        ASSERTL0(false, "The velocities in x must be specified via an equation "
                        "<E> or a file <F>");
    }

    // Check Velocity y
    vType = m_session->GetFunctionType(m_funcName[1],m_motion[1]);
    if (vType == LibUtilities::eFunctionTypeExpression)
    {
        m_IsFromFile[3] = false;
    }
    else if (vType == LibUtilities::eFunctionTypeFile)
    {
        m_IsFromFile[3] = true;
    }
    else
    {
        ASSERTL0(false, "The velocities in y must be specified via an equation "
                        "<E> or a file <F>");
    }

    // Check Acceleration x
    vType = m_session->GetFunctionType(m_funcName[2],m_motion[0]);
    if (vType == LibUtilities::eFunctionTypeExpression)
    {
        m_IsFromFile[4] = false;
    }
    else if (vType == LibUtilities::eFunctionTypeFile)
    {
        m_IsFromFile[4] = true;
    }
    else
    {
        ASSERTL0(false, "The accelerations in x must be specified via an "
                        "equation <E> or a file <F>");
    }

    // Check Acceleration y
    vType = m_session->GetFunctionType(m_funcName[2],m_motion[1]);
    if (vType == LibUtilities::eFunctionTypeExpression)
    {
        m_IsFromFile[5] = false;
    }
    else if (vType == LibUtilities::eFunctionTypeFile)
    {
        m_IsFromFile[5] = true;
    }
    else
    {
        ASSERTL0(false, "The accelerations in y must be specified via an "
                        "equation <E> or a file <F>");
    }
}

/**
 *
 */
void ForcingMovingBody::InitialiseFilter(
        const LibUtilities::SessionReaderSharedPtr& pSession,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
        const TiXmlElement* pForce)
{
    // Get the outputfile name, output frequency and
    // the boundary's ID for the cable's wall
    std::string typeStr = pForce->Attribute("TYPE");
    std::map<std::string, std::string> vParams;

    const TiXmlElement *param = pForce->FirstChildElement("PARAM");
    while (param)
    {
        ASSERTL0(param->Attribute("NAME"),
                 "Missing attribute 'NAME' for parameter in filter "
                 + typeStr + "'.");
        std::string nameStr = param->Attribute("NAME");

        ASSERTL0(param->GetText(), "Empty value string for param.");
        std::string valueStr = param->GetText();

        vParams[nameStr] = valueStr;

        param = param->NextSiblingElement("PARAM");
    }

    // Creat the filter for MovingBody, where we performed the calculation of
    // fluid forces and write both the aerodynamic forces and motion variables
    // into the output files
    m_MovBodyfilter = MemoryManager<FilterMovingBody>::
                                    AllocateSharedPtr(pSession, m_equ, vParams);

    // Initialise the object of MovingBody filter
    m_MovBodyfilter->Initialise(pFields, 0.0);

}

}
