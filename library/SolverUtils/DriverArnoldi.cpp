///////////////////////////////////////////////////////////////////////////////
//
// File DriverArnoldi.cpp
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
// Description: Base Driver class for the stability solver
//
///////////////////////////////////////////////////////////////////////////////
#include <iomanip>

#include <SolverUtils/DriverArnoldi.h>

using namespace std;

namespace Nektar
{
namespace SolverUtils
{

/**
 * Constructor
 */
DriverArnoldi::DriverArnoldi(
    const LibUtilities::SessionReaderSharedPtr pSession,
    const SpatialDomains::MeshGraphSharedPtr pGraph)
    : Driver(pSession, pGraph)
{
    m_session->LoadParameter("IO_InfoSteps", m_infosteps, 1);
}

int c01;
Array<OneD, NekDouble> m_y1(2);
NekDouble m_alpha0;
/**
 * Destructor
 */
DriverArnoldi::~DriverArnoldi()
{
}


/**
 * Arnoldi driver initialisation
 */
void DriverArnoldi::v_InitObject(ostream &out)
{
    Driver::v_InitObject(out);
    m_session->MatchSolverInfo("SolverType",
                               "VelocityCorrectionScheme",
                               m_timeSteppingAlgorithm, false);

    if (m_timeSteppingAlgorithm)
    {
        m_period  = m_session->GetParameter("TimeStep")
                  * m_session->GetParameter("NumSteps");
        m_nfields = m_equ[0]->UpdateFields().size() - 1;
        nq        = m_equ[0]->UpdateFields()[0]->GetNcoeffs();
        n_tot     = (m_nfields+2)*nq;

    }
    else
    {
        m_period  = 1.0;
        m_nfields = m_equ[0]->UpdateFields().size();
    }

    if(m_session->DefinesSolverInfo("ModeType") &&
       (boost::iequals(m_session->GetSolverInfo("ModeType"),
                       "SingleMode")||
        boost::iequals(m_session->GetSolverInfo("ModeType"),
                       "HalfMode")))
    {
        for(int i = 0; i < m_nfields; ++i)
        {
            m_equ[0]->UpdateFields()[i]->SetWaveSpace(true);
        }
    }
    m_negatedOp = m_equ[0]->v_NegatedOp();

    m_session->LoadParameter("kdim",  m_kdim,  16);
    m_session->LoadParameter("nvec",  m_nvec,  2);
    m_session->LoadParameter("nits",  m_nits,  500);
    m_session->LoadParameter("evtol", m_evtol, 1e-06);

    ASSERTL0( m_kdim >= m_nvec, "nvec cannot be larger than kdim.");

    m_session->LoadParameter("realShift", m_realShift, 0.0);
    m_equ[0]->SetLambda(m_realShift);

    m_session->LoadParameter("imagShift",m_imagShift,0.0);

    m_y1[0] = 0.;
    m_y1[1] = 0.;
    m_alpha = 0.;
    // The imaginary shift is applied at system assembly
    // Only if using HOMOGENEOUS expansion and ModeType set to SingleMode
    if(m_imagShift != 0.0)
    {
        if(!m_session->DefinesSolverInfo("HOMOGENEOUS")&&!m_session->DefinesSolverInfo("ModeType"))
        {
            NEKERROR(ErrorUtil::efatal, "Imaginary shift only supported with HOMOGENEOUS expansion and ModeType set to SingleMode");
        }
        else if(!boost::iequals(m_session->GetSolverInfo("ModeType"),"SingleMode"))
        {
            NEKERROR(ErrorUtil::efatal, "Imaginary shift only supported with HOMOGENEOUS expansion and ModeType set to SingleMode");
        }
    }

}

void DriverArnoldi::ArnoldiSummary(std::ostream &out)
{
    if (m_comm->GetRank() == 0)
    {
        if(m_session->DefinesSolverInfo("ModeType") &&
           boost::iequals(m_session->GetSolverInfo("ModeType"),
                          "SingleMode"))
        {
            out << "\tSingle Fourier mode    : true " << endl;
            ASSERTL0(m_session->DefinesSolverInfo("Homogeneous"),
                     "Expected a homogeneous expansion to be defined "
                     "with single mode");
        }
        else
        {
            out << "\tSingle Fourier mode    : false " << endl;
        }
        if(m_session->DefinesSolverInfo("BetaZero"))
        {
            out << "\tBeta set to Zero       : true (overrides LHom)"
                << endl;
        }
        else
        {
            out << "\tBeta set to Zero       : false " << endl;
        }

        if(m_timeSteppingAlgorithm)
        {
            out << "\tEvolution operator     : "
                << m_session->GetSolverInfo("EvolutionOperator")
                << endl;
        }
        else
        {
            out << "\tShift (Real,Imag)      : " << m_realShift
                << "," << m_imagShift <<  endl;
        }
        out << "\tKrylov-space dimension : " << m_kdim << endl;
        out << "\tNumber of vectors      : " << m_nvec << endl;
        out << "\tMax iterations         : " << m_nits << endl;
        out << "\tEigenvalue tolerance   : " << m_evtol << endl;
        out << "======================================================"
            << endl;
    }
}

/**
 * Copy Arnoldi array to field variables which depend from
 * either the m_fields or m_forces
 */
void DriverArnoldi::CopyArnoldiArrayToField(Array<OneD, NekDouble> &array)
{

    Array<OneD, MultiRegions::ExpListSharedPtr>& fields = m_equ[0]->UpdateFields();
    int nq = fields[0]->GetNcoeffs();

    for (int k = 0; k < m_nfields; ++k)
    {
        Vmath::Vcopy(nq, &array[k*nq], 1, &fields[k]->UpdateCoeffs()[0], 1);
        fields[k]->SetPhysState(false);
    }
    m_y1[0] = array[nq*m_nfields];
    m_y1[1] = array[nq*(m_nfields+1)];
    m_alpha0 = m_alpha;
    cout << "disp0: "<< m_y1[0]<< ", vel0: " << m_y1[1] << endl;

}

/**
 * Copy field variables which depend from either the m_fields
 * or m_forces array the Arnoldi array
 */
void DriverArnoldi::CopyFieldToArnoldiArray(Array<OneD, NekDouble> &array)
{
    c01 = 0;

    Array<OneD, MultiRegions::ExpListSharedPtr> fields;

    if (m_EvolutionOperator == eAdaptiveSFD)
    {
        // This matters for the Adaptive SFD method because
        // m_equ[1] is the nonlinear problem with non
        // homogeneous BCs.
        fields = m_equ[0]->UpdateFields();
    }
    else
    {
        fields = m_equ[m_nequ-1]->UpdateFields();
    }
    
    int nq = fields[0]->GetNcoeffs();
    for (int k = 0; k < m_nfields; ++k)
    {
        Vmath::Vcopy(nq,  &fields[k]->GetCoeffs()[0], 1, &array[k*nq], 1);
        fields[k]->SetPhysState(false);

        // cout << "array size: " << array.size() << ", nq: " << m_nfields*nq << endl;

    }
    
    // array[n_tot-2] = m_y1[0];
    // array[n_tot-1] = m_y1[1];    
    Array<OneD, NekDouble> tmp1(nq, m_y1[0]);
    Array<OneD, NekDouble> tmp2(nq, m_y1[1]);
    Vmath::Vcopy(nq, &tmp1[0], 1, &array[m_nfields*nq], 1);
    Vmath::Vcopy(nq, &tmp2[0], 1, &array[(m_nfields+1)*nq], 1); 
    cout << "disp:"<< m_y1[0]<< ", vel: " << m_y1[1] << endl;

      
    // Vmath::Zero(2, m_y1, 1);
}


/**
 * Initialisation for the transient growth
 */
void DriverArnoldi::CopyFwdToAdj()
{
    c01=1;
    Array<OneD, MultiRegions::ExpListSharedPtr> fields;

    ASSERTL0(m_timeSteppingAlgorithm,
            "Transient Growth non available for Coupled Solver");

    fields = m_equ[0]->UpdateFields();
    int nq = fields[0]->GetNcoeffs();

    for (int k=0 ; k < m_nfields; ++k)
    {
        Vmath::Vcopy(nq,
                     &fields[k]->GetCoeffs()[0],                      1,
                     &m_equ[1]->UpdateFields()[k]->UpdateCoeffs()[0], 1);
    }
}

void DriverArnoldi::WriteFld(std::string file, std::vector<Array<OneD, NekDouble> > coeffs)
{

    std::vector<std::string>  variables(m_nfields);

    ASSERTL1(coeffs.size() >= m_nfields, "coeffs is not of the correct length");
    for(int i = 0; i < m_nfields; ++i)
    {
        variables[i] = m_equ[0]->GetVariable(i);
    }

    m_equ[0]->WriteFld(file,m_equ[0]->UpdateFields()[0], coeffs, variables);
}


void DriverArnoldi::WriteFld(std::string file, Array<OneD, NekDouble> coeffs)
{

    std::vector<std::string>  variables(m_nfields);
    std::vector<Array<OneD, NekDouble> > fieldcoeffs(m_nfields);

    int ncoeffs = m_equ[0]->UpdateFields()[0]->GetNcoeffs();
    ASSERTL1(coeffs.size() >= ncoeffs*m_nfields,"coeffs is not of sufficient size");

    for(int i = 0; i < m_nfields; ++i)
    {
        variables[i] = m_equ[0]->GetVariable(i);
        fieldcoeffs[i] = coeffs + i*ncoeffs;
    }

    m_equ[0]->WriteFld(file,m_equ[0]->UpdateFields()[0], fieldcoeffs, variables);
}

void DriverArnoldi::WriteEvs(
        ostream &evlout,
        const int i,
        const NekDouble re_ev,
        const NekDouble im_ev,
        NekDouble resid,
        bool DumpInverse)
{
    if (m_timeSteppingAlgorithm)
    {
        NekDouble abs_ev = hypot (re_ev, im_ev);
        NekDouble ang_ev = atan2 (im_ev, re_ev);

        evlout << "EV: " << setw(2)  << i
               << setw(12) << abs_ev
               << setw(12) << ang_ev
               << setw(12) << log (abs_ev) / m_period
               << setw(12) << ang_ev       / m_period;

        if(resid != NekConstants::kNekUnsetDouble)
        {
            evlout << setw(12) << resid;
        }
        evlout << endl;
    }
    else
    {
        NekDouble invmag = 1.0/(re_ev*re_ev + im_ev*im_ev);
        NekDouble sign;
        if(m_negatedOp)
        {
            sign = -1.0;
        }
        else
        {
            sign = 1.0;
        }

        evlout << "EV: " << setw(2)  <<  i
               << setw(14) <<  sign*re_ev
               << setw(14) <<  sign*im_ev;

        if(DumpInverse)
        {
            evlout << setw(14) <<  sign*re_ev*invmag + m_realShift
                   << setw(14) <<  sign*im_ev*invmag + m_imagShift;
        }

        if(resid != NekConstants::kNekUnsetDouble)
        {
            evlout << setw(12) << resid;
        }
        evlout << endl;
    }
}

void DriverArnoldi::ReturnStructVector(NekDouble & y, NekDouble& y1, int& c)
{
    y  = m_y1[0];
    y1 = m_y1[1];
    c  = c01;
}

void DriverArnoldi::GetStructVector(NekDouble y, NekDouble y1)
{
    m_y1[0] = y;
    m_y1[1] = y1;
}

void DriverArnoldi::Getalpha(NekDouble &alpha)
{
    alpha = m_alpha0;

}
}
}
