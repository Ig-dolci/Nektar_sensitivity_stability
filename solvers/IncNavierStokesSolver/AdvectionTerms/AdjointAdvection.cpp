///////////////////////////////////////////////////////////////////////////////
//
// File AdjointAdvection.cpp
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
// Description: Evaluation of the adjoint advective term
//
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/AdvectionTerms/AdjointAdvection.h>

using namespace std;

namespace Nektar
{

string AdjointAdvection::className = SolverUtils
        ::GetAdvectionFactory().RegisterCreatorFunction("Adjoint",
                                                    AdjointAdvection::create);

/**
 *
 */
AdjointAdvection::AdjointAdvection():
    LinearisedAdvection()
{
}

AdjointAdvection::~AdjointAdvection()
{
}

void AdjointAdvection::v_Advect(
    const int nConvectiveFields,
    const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
    const Array<OneD, Array<OneD, NekDouble> >        &advVel,
    const Array<OneD, Array<OneD, NekDouble> >        &inarray,
    Array<OneD, Array<OneD, NekDouble> >              &outarray,
    const NekDouble                                   &time,
    const Array<OneD, Array<OneD, NekDouble> >        &pFwd,
    const Array<OneD, Array<OneD, NekDouble> >        &pBwd)
{
    ASSERTL1(nConvectiveFields == inarray.size(),
             "Number of convective fields and Inarray are not compatible");

    int nPointsTot  = fields[0]->GetNpoints();
    int ndim        = advVel.size();
    int nBaseDerivs = (m_halfMode || m_singleMode) ? 2 : m_spacedim;
    int nDerivs     = (m_halfMode) ? 2 : m_spacedim;

    Array<OneD, Array<OneD, NekDouble> > velocity(ndim);
    for(int i = 0; i < ndim; ++i)
    {
        if(fields[i]->GetWaveSpace() && !m_singleMode && !m_halfMode)
        {
            velocity[i] = Array<OneD, NekDouble>(nPointsTot,0.0);
            fields[i]->HomogeneousBwdTrans(advVel[i],velocity[i]);
        }
        else
        {
            velocity[i] = advVel[i];
        }
    }

    Array<OneD, Array<OneD, NekDouble> > grad (nDerivs);
    for( int i = 0; i < nDerivs; ++i)
    {
        grad[i] = Array<OneD, NekDouble> (nPointsTot);
    }

    // Evaluation of the base flow for periodic cases
    // Evaluation of the base flow for periodic cases
    if (m_slices > 1 && m_session->GetSolverInfo("BaseFlowform") == "Periodic")
    {
        for (int i = 0; i < ndim; ++i)
        {
            UpdateBase(m_slices, m_interp[i], m_baseflow[i],
                       time, m_period);
            UpdateGradBase(i, fields[i]);
        }
    }
    if (m_slices > 1 && m_session->GetSolverInfo("BaseFlowform") == "TimeDependent")
    {
        //BaseFlow from file
        if (m_session->GetFunctionType("BaseFlow", m_session->GetVariable(0))
            == LibUtilities::eFunctionTypeFile)
        {

            string file = m_session->GetFunctionFilename("BaseFlow", 0);
            int step = contador;
            // time/m_session->GetParameter("TimeStep");
            int num_step = m_session->GetParameter("NumSteps");
            
            int aux =time;
            int period = m_session->GetParameter("period");
            if (contador%m_slices==0)
            {  
                contador=0;

            }
            int adj_step = m_slices - contador;
            
            size_t found = file.find("%d");
            ASSERTL0(found != string::npos && file.find("%d", found+1) == string::npos,
                    "Since N_slices is specified, the filename provided for function "
                    "'BaseFlow' must include exactly one instance of the format "
                    "specifier '%d', to index the time-slices.");
            char* buffer = new char[file.length() + 8];

            sprintf(buffer, file.c_str(), adj_step);
            ImportFldBase(buffer,fields,1);
            
            delete[] buffer;
        }

        for (int i = 0; i < ndim; ++i)
        {
            UpdateGradBase(i, fields[i]);
        }
    }

    ++contador;
    //Evaluate the linearised advection term
    for( int i = 0; i < ndim; ++i)
    {
        // Calculate gradient
        switch(nDerivs)
        {
            case 1:
            {
                fields[i]->PhysDeriv(inarray[i], grad[0]);
            }
            break;
            case 2:
            {
                fields[i]->PhysDeriv(inarray[i], grad[0], grad[1]);
            }
            break;
            case 3:
            {
                fields[i]->PhysDeriv(inarray[i], grad[0], grad[1], grad[2]);
                if(m_multipleModes)
                {
                    // transform gradients into physical Fourier space
                    fields[i]->HomogeneousBwdTrans(grad[0], grad[0]);
                    fields[i]->HomogeneousBwdTrans(grad[1], grad[1]);
                    fields[i]->HomogeneousBwdTrans(grad[2], grad[2]);
                }
            }
            break;
        }
        
        Array<OneD, NekDouble> aux(nPointsTot, 0.0);
        NekDouble integ;

        Vmath::Vmul(nPointsTot,m_baseflow[0], 1, velocity[0], 1, aux, 1);
        Vmath::Vvtvp(nPointsTot,m_baseflow[1], 1, velocity[1], 1,aux, 1, aux, 1);
        
        
        integ = fields[0]->Integral(aux);
        // cout<< integ << endl;
        // Calculate -U_j du'_i/dx_j
        Vmath::Vmul(nPointsTot,grad[0], 1, m_baseflow[0], 1, outarray[i], 1);
        for( int j = 1; j < nDerivs; ++j)
        {
            // Vmath::Smul(nPointsTot,1./integ, velocity[j], 1, velocity[j], 1);
            Vmath::Vvtvp(nPointsTot,grad[j], 1,
                                    m_baseflow[j], 1,
                                    outarray[i], 1,
                                    outarray[i], 1);
        }
        Vmath::Neg(nPointsTot,outarray[i],1);

        // Add u'_j U_j/ dx_i
        int lim = (m_halfMode) ? 2 : ndim;
        if ( (m_halfMode || m_singleMode) && i==2)
        {
            lim = 0;
        }
        for( int j = 0; j < lim; ++j)
        {
            Vmath::Vvtvp(nPointsTot,m_gradBase[j*nBaseDerivs + i], 1,
                                    velocity[j], 1,
                                    outarray[i], 1,
                                    outarray[i], 1);
        }

        if(m_multipleModes)
        {
            fields[i]->HomogeneousFwdTrans(outarray[i],outarray[i]);
        }
        Vmath::Neg(nPointsTot,outarray[i],1);
    }
    
}

} //end of namespace

