///////////////////////////////////////////////////////////////////////////////
//
// File: ForcingBody.cpp
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
// Description: Body forcing
//
///////////////////////////////////////////////////////////////////////////////

#include <SolverUtils/Forcing/ForcingBody.h>
#include <MultiRegions/ExpList.h>

using namespace std;

namespace Nektar
{
namespace SolverUtils
{

    std::string ForcingBody::classNameBody = GetForcingFactory().
        RegisterCreatorFunction("Body",
                                ForcingBody::create,
                                "Body Forcing");
    std::string ForcingBody::classNameField = GetForcingFactory().
        RegisterCreatorFunction("Field",
                                ForcingBody::create,
                                "Field Forcing");

    ForcingBody::ForcingBody(
            const LibUtilities::SessionReaderSharedPtr &pSession,
            const std::weak_ptr<EquationSystem>      &pEquation)
        : Forcing(pSession, pEquation),
          m_hasTimeFcnScaling(false)
    {
    }
    NekDouble y, ac;
    void ForcingBody::v_InitObject(
                                   const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields,
                                   const unsigned int& pNumForcingFields,
                                   const TiXmlElement* pForce)
    {
        m_NumVariable = pNumForcingFields;
        int nq         = pFields[0]->GetTotPoints();
        y=0.;
        ac=0.;
        const TiXmlElement* funcNameElmt = pForce->FirstChildElement("BODYFORCE");
        if(!funcNameElmt)
        {
            funcNameElmt = pForce->FirstChildElement("FIELDFORCE");

            ASSERTL0(funcNameElmt, "Requires BODYFORCE or FIELDFORCE tag "
                     "specifying function name which prescribes body force.");
        }

        m_funcName = funcNameElmt->GetText();
        ASSERTL0(m_session->DefinesFunction(m_funcName),
                 "Function '" + m_funcName + "' not defined.");

        bool singleMode, halfMode;
        m_session->MatchSolverInfo("ModeType","SingleMode",singleMode,false);
        m_session->MatchSolverInfo("ModeType","HalfMode",  halfMode,  false);
        bool homogeneous = pFields[0]->GetExpType() == MultiRegions::e3DH1D ||
                           pFields[0]->GetExpType() == MultiRegions::e3DH2D;
        m_transform = (singleMode || halfMode || homogeneous);

        // Time function is optional
        funcNameElmt = pForce->FirstChildElement("BODYFORCETIMEFCN");
        if(!funcNameElmt)
        {
            funcNameElmt = pForce->FirstChildElement("FIELDFORCETIMEFCN");
        }

        // Load time function if specified
        if(funcNameElmt)
        {
            std::string funcNameTime = funcNameElmt->GetText();

            ASSERTL0(!funcNameTime.empty(),
                     "Expression must be given in BODYFORCETIMEFCN or "
                     "FIELDFORCETIMEFCN.");

            m_session->SubstituteExpressions(funcNameTime);
            m_timeFcnEqn = MemoryManager<LibUtilities::Equation>
                ::AllocateSharedPtr(m_session->GetInterpreter(),funcNameTime);

            m_hasTimeFcnScaling = true;
        }

        m_Forcing = Array<OneD, Array<OneD, NekDouble> > (m_NumVariable);
        for (int i = 0; i < m_NumVariable; ++i)
        {
            m_Forcing[i] = Array<OneD, NekDouble> (pFields[0]->GetTotPoints(), 0.0);
        }

        baseflow_sens = Array<OneD, Array<OneD, NekDouble> >(2);

        for (int i = 0; i < 2; ++i)
        {
            baseflow_sens[i] = Array<OneD, NekDouble>(pFields[i]->GetTotPoints(), 0.0);
        }
        int nvar=2;
        std::string control = m_session->GetSolverInfo("Control");
        if(m_session->DefinesParameter("SteadyF_sens"))
        {

            Array<OneD, Array<OneD, NekDouble> > veldr, veldi, velar, velai, deltaU;
            veldr = Array<OneD, Array<OneD, NekDouble> >(nvar);
            veldi = Array<OneD, Array<OneD, NekDouble> >(nvar);
            velar = Array<OneD, Array<OneD, NekDouble> >(nvar);
            velai = Array<OneD, Array<OneD, NekDouble> >(nvar);
            deltaU = Array<OneD, Array<OneD, NekDouble> >(nvar);
            Array<OneD, NekDouble> aux(nq, 0.0);
            Array<OneD, NekDouble> aux1(nq, 0.0);

            for (int i = 0; i < nvar; ++i)
            {
                veldr[i] = Array<OneD, NekDouble>(pFields[i]->GetTotPoints(), 0.0);
                veldi[i] = Array<OneD, NekDouble>(pFields[i]->GetTotPoints(), 0.0);
                velar[i] = Array<OneD, NekDouble>(pFields[i]->GetTotPoints(), 0.0);
                velai[i] = Array<OneD, NekDouble>(pFields[i]->GetTotPoints(), 0.0);
                deltaU[i] = Array<OneD, NekDouble>(pFields[i]->GetTotPoints(), 0.0);
            }
            string file1 = m_session->GetFunctionFilename("Direct_Real", 0);
            string file2 = m_session->GetFunctionFilename("Direct_Imag", 0);
            string file3 = m_session->GetFunctionFilename("Adjoint_Real", 0);
            string file4 = m_session->GetFunctionFilename("Adjoint_Imag", 0);
            string file5 = m_session->GetFunctionFilename("dU", 0);

            std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef0, FieldDef1, FieldDef2, FieldDef3, FieldDef4;
            std::vector<std::vector<NekDouble> > FieldData0, FieldData1, FieldData2, FieldData3, FieldData4;

            int nq = pFields[0]->GetTotPoints();
            Array<OneD, NekDouble> tmp_coeff0(pFields[0]->GetNcoeffs(), 0.0);
            Array<OneD, NekDouble> tmp_coeff1(pFields[0]->GetNcoeffs(), 0.0);
            Array<OneD, NekDouble> tmp_coeff2(pFields[0]->GetNcoeffs(), 0.0);
            Array<OneD, NekDouble> tmp_coeff3(pFields[0]->GetNcoeffs(), 0.0);
            Array<OneD, NekDouble> tmp_coeff4(pFields[0]->GetNcoeffs(), 0.0);

            int numexp = pFields[0]->GetExpSize();
            Array<OneD,int> ElementGIDs(numexp);

            // Define list of global element ids
            for(int i = 0; i < numexp; ++i)
            {
                ElementGIDs[i] = pFields[0]->GetExp(i)->GetGeom()->GetGlobalID();
            }

            LibUtilities::FieldIOSharedPtr fld = LibUtilities::FieldIO::CreateForFile(
                m_session, file1);

            // fld->Import(file1, FieldDef, FieldData,
            //             LibUtilities::NullFieldMetaDataMap,
            //             ElementGIDs);
            fld->Import(file1, FieldDef0, FieldData0, LibUtilities::NullFieldMetaDataMap,
            ElementGIDs);
            fld->Import(file2, FieldDef1, FieldData1, LibUtilities::NullFieldMetaDataMap,
            ElementGIDs);
            fld->Import(file3, FieldDef2, FieldData2, LibUtilities::NullFieldMetaDataMap,
            ElementGIDs);
            fld->Import(file4, FieldDef3, FieldData3, LibUtilities::NullFieldMetaDataMap,
            ElementGIDs);
            fld->Import(file5, FieldDef4, FieldData4, LibUtilities::NullFieldMetaDataMap,
            ElementGIDs);



            for(int j = 0; j < 2; ++j)
            {

                for(int i = 0; i < FieldDef0.size(); ++i)
                {

                    // ASSERTL0(flag, (std::string("Order of ") + file
                    //                 + std::string(" data and that defined in "
                    //                 "the session file differs")).c_str());

                    pFields[j]->ExtractDataToCoeffs(
                                        FieldDef0[i],
                                        FieldData0[i],
                                        FieldDef0[i]->m_fields[j],
                                        tmp_coeff0);
                    pFields[j]->ExtractDataToCoeffs(
                                        FieldDef1[i],
                                        FieldData1[i],
                                        FieldDef1[i]->m_fields[j],
                                        tmp_coeff1);
                    pFields[j]->ExtractDataToCoeffs(
                                        FieldDef2[i],
                                        FieldData2[i],
                                        FieldDef2[i]->m_fields[j],
                                        tmp_coeff2);
                    pFields[j]->ExtractDataToCoeffs(
                                        FieldDef3[i],
                                        FieldData3[i],
                                        FieldDef3[i]->m_fields[j],
                                        tmp_coeff3);
                    pFields[j]->ExtractDataToCoeffs(
                                        FieldDef4[i],
                                        FieldData4[i],
                                        FieldDef4[i]->m_fields[j],
                                        tmp_coeff4);

                }
                bool oldwavespace = pFields[j]->GetWaveSpace();
                pFields[j]->SetWaveSpace(false);
                pFields[j]->BwdTrans(tmp_coeff0, veldr[j]);
                pFields[j]->BwdTrans(tmp_coeff1, veldi[j]);
                pFields[j]->BwdTrans(tmp_coeff2, velar[j]);
                pFields[j]->BwdTrans(tmp_coeff3, velai[j]);
                pFields[j]->BwdTrans(tmp_coeff4, deltaU[j]);
                pFields[j]->SetWaveSpace(oldwavespace);


            }

            // //   /* code */
            // for (int j = 0; j < velai[0].size(); ++j)
            // {
            //     std::complex<double> ud(veldr[0][j], veldi[0][j]);
            //     std::complex<double> vd(veldr[1][j], veldi[1][j]);
            //     std::complex<double> ua(velar[0][j], velai[0][j]);
            //     std::complex<double> va(velar[1][j], velai[1][j]);
            //     aux[j] = real(ud*conj(ud) + vd*conj(vd));
            //     aux1[j] = real(ua*conj(ua) + va*conj(va));

            // }

            // NekDouble norm_d = sqrt(pFields[0]->Integral(aux));
            // NekDouble norm_a = sqrt(pFields[0]->Integral(aux1));
            // cout << norm_d << endl;
            // cout << norm_a << endl;

            // for (int j = 0; j < velai[0].size(); ++j)
            // {
            //     veldr[0][j] = veldr[0][j]/norm_d;
            //     veldr[1][j] = veldr[1][j]/norm_d;
            //     veldi[0][j] = veldi[0][j]/norm_d;
            //     veldi[1][j] = veldi[1][j]/norm_d;
            //     velar[0][j] = velar[0][j]/(norm_a);
            //     velar[1][j] = velar[1][j]/(norm_a);
            //     velai[0][j] = velai[0][j]/(norm_a);
            //     velai[1][j] = velai[1][j]/(norm_a);
            // }
            //   /* code */
            for (int j = 0; j < velai[0].size(); ++j)
            {
                std::complex<double> ud(veldr[0][j], veldi[0][j]);
                std::complex<double> vd(veldr[1][j], veldi[1][j]);
                std::complex<double> ua(velar[0][j], velai[0][j]);
                std::complex<double> va(velar[1][j], velai[1][j]);
                aux[j] = real(ua*ud + va*vd);
                aux1[j] = imag(ua*ud + va*vd);
                // aux[j] = real(conj(ua)*ud + conj(va)*vd);
                // aux1[j] = imag(conj(ua)*ud + conj(va)*vd);

            }
        

            NekDouble int_real = pFields[0]->Integral(aux);
            NekDouble int_imag = pFields[0]->Integral(aux1);
            
            std::complex<double> mycomplex(int_real, int_imag);
            cout<<mycomplex<< endl;
            Array<OneD, NekDouble> u_xdr(nq, 0.0);
            Array<OneD, NekDouble> v_xdr(nq, 0.0);
            Array<OneD, NekDouble> u_ydr(nq, 0.0);
            Array<OneD, NekDouble> v_ydr(nq, 0.0);

            Array<OneD, NekDouble> u_xdi(nq, 0.0);
            Array<OneD, NekDouble> v_xdi(nq, 0.0);
            Array<OneD, NekDouble> u_ydi(nq, 0.0);
            Array<OneD, NekDouble> v_ydi(nq, 0.0);

            Array<OneD, NekDouble> u_xd(nq, 0.0);
            Array<OneD, NekDouble> v_xd(nq, 0.0);
            Array<OneD, NekDouble> u_yd(nq, 0.0);
            Array<OneD, NekDouble> v_yd(nq, 0.0);



            Array<OneD, NekDouble> u_xar(nq, 0.0);
            Array<OneD, NekDouble> v_xar(nq, 0.0);
            Array<OneD, NekDouble> u_yar(nq, 0.0);
            Array<OneD, NekDouble> v_yar(nq, 0.0);

            Array<OneD, NekDouble> u_xai(nq, 0.0);
            Array<OneD, NekDouble> v_xai(nq, 0.0);
            Array<OneD, NekDouble> u_yai(nq, 0.0);
            Array<OneD, NekDouble> v_yai(nq, 0.0);


            Array<OneD, NekDouble> u_xa(nq, 0.0);
            Array<OneD, NekDouble> v_xa(nq, 0.0);
            Array<OneD, NekDouble> u_ya(nq, 0.0);
            Array<OneD, NekDouble> v_ya(nq, 0.0);

            pFields[0]->PhysDeriv(veldr[0],u_xdr,u_ydr);
            pFields[0]->PhysDeriv(veldr[1],v_xdr,v_ydr);
            pFields[0]->PhysDeriv(veldi[0],u_xdi,u_ydi);
            pFields[0]->PhysDeriv(veldi[1],v_xdi,v_ydi);

            pFields[0]->PhysDeriv(velar[0],u_xar,u_yar);
            pFields[0]->PhysDeriv(velar[1],v_xar,v_yar);
            pFields[0]->PhysDeriv(velai[0],u_xai,u_yai);
            pFields[0]->PhysDeriv(velai[1],v_xai,v_yai);



            for (int j = 0; j < velai[0].size(); ++j)
            {
                
                /* code */
                std::complex<double> ud(veldr[0][j], veldi[0][j]);
                std::complex<double> vd(veldr[1][j], veldi[1][j]);
                std::complex<double> ua(velar[0][j], velai[0][j]);
                std::complex<double> va(velar[1][j], velai[1][j]);

                std::complex<double> u_xd(u_xdr[j] ,  u_xdi[j]);
                std::complex<double> v_xd(v_xdr[j] ,  v_xdi[j]);
                std::complex<double> u_yd(u_ydr[j] ,  u_ydi[j]);
                std::complex<double> v_yd(v_ydr[j] ,  v_ydi[j]);
                std::complex<double> u_xa(u_xar[j] ,  u_xai[j]);
                std::complex<double> v_xa(v_xar[j] ,  v_xai[j]);
                std::complex<double> u_ya(u_yar[j] ,  u_yai[j]);
                std::complex<double> v_ya(v_yar[j] ,  v_yai[j]);

                ua = ua/mycomplex;
                va = va/mycomplex;
                u_xa = u_xa/mycomplex;
                v_xa = v_xa/mycomplex;
                u_ya = u_ya/mycomplex;
                v_ya = v_ya/mycomplex;

                // std::complex<double> grad_sigma1 = (- conj(u_xa)*ud - conj(u_ya)*vd + u_xd*conj(ua) + v_xd*conj(va));
                // std::complex<double> grad_sigma2 = (- conj(v_xa)*ud - conj(v_ya)*vd + u_yd*conj(ua) + v_yd*conj(va));

                // // std::complex<double> grad_sigma1 = (- u_xa*conj(ud) - u_ya*conj(vd) + conj(u_xd)*ua + conj(v_xd)*va);
                // std::complex<double> grad_sigma2 = (- v_xa*conj(ud) - v_ya*conj(vd) + conj(u_yd)*ua + conj(v_yd)*va);
                std::complex<double> grad_sigma1 = (+ u_xa*ud + u_ya*vd - u_xd*ua - v_xd*va);
                std::complex<double> grad_sigma2 = (+ v_xa*ud + v_ya*vd - u_yd*ua - v_yd*va);
            
                std::string evol_operator = m_session->GetSolverInfo("ForcingSensitivity");

                if(boost::iequals(evol_operator, "Frequency"))
                {
                    baseflow_sens[0][j] = imag(grad_sigma1);
                    baseflow_sens[1][j] = imag(grad_sigma2);
    
                    // baseflow_sens[0][j] = imag(grad_sigma1)*deltaU[0][j];
                    // baseflow_sens[1][j] = imag(grad_sigma2)*deltaU[1][j];
                }
                else if(boost::iequals(evol_operator, "Growth"))
                {
                    
                    baseflow_sens[0][j] = real(grad_sigma1);
                    baseflow_sens[1][j] = real(grad_sigma2);
                    // baseflow_sens[0][j] = real(grad_sigma1)*deltaU[0][j];
                    // baseflow_sens[1][j] = real(grad_sigma2)*deltaU[1][j];

                    

                }

                aux[j]  = real(grad_sigma1)*deltaU[0][j];
                aux1[j] = real(grad_sigma2)*deltaU[1][j];


            }
            NekDouble v = pFields[0]->PhysIntegral(aux);
            NekDouble u = pFields[0]->PhysIntegral(aux1);

            // NekDouble v = pFields[0]->PhysIntegral(baseflow_sens[0]);
            // NekDouble u = pFields[0]->PhysIntegral(baseflow_sens[1]);
            cout << v + u << endl;



        }
        else if(boost::iequals(control, "Activate"))
        {
            /* code */
            Array<OneD, Array<OneD, NekDouble> > tmp(pFields.size());
            for (int i = 0; i < pFields.size(); ++i)
            {
                tmp[i] = pFields[i]->GetPhys();
            }
            Update(pFields, tmp, 0.0);
            /* code */
            Array<OneD, NekDouble> aux(nq, 0.0);
            Array<OneD, Array<OneD, NekDouble> > Vel, grad;
            Vel = Array<OneD, Array<OneD, NekDouble> >(nvar);
            grad = Array<OneD, Array<OneD, NekDouble> >(nvar);
            for (int i = 0; i < nvar; ++i)
            {
                Vel[i] = Array<OneD, NekDouble>(pFields[i]->GetTotPoints(), 0.0);
                grad[i] = Array<OneD, NekDouble>(pFields[i]->GetTotPoints(), 0.0);
            }
            string file1 = m_session->GetFunctionFilename("BaseField", 0);
            string file2 = m_session->GetFunctionFilename("grad", 0);
            // NekDouble value = index/
            std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef, FieldDef1;
            std::vector<std::vector<NekDouble> > FieldData, FieldData1;
            // std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef;
            // std::vector<std::vector<NekDouble> > FieldData;
            
            Array<OneD, NekDouble> tmp_coeff0(pFields[0]->GetNcoeffs(), 0.0);
            Array<OneD, NekDouble> tmp_coeff1(pFields[0]->GetNcoeffs(), 0.0);
            int numexp = pFields[0]->GetExpSize();
            Array<OneD,int> ElementGIDs(numexp);

            // Define list of global element ids
            for(int i = 0; i < numexp; ++i)
            {
                ElementGIDs[i] = pFields[0]->GetExp(i)->GetGeom()->GetGlobalID();

            }


            LibUtilities::FieldIOSharedPtr fld = LibUtilities::FieldIO::CreateForFile(
                m_session, file1);
          
            fld->Import(file1, FieldDef, FieldData,
                        LibUtilities::NullFieldMetaDataMap,
                        ElementGIDs);
            fld->Import(file2, FieldDef1, FieldData1,
                        LibUtilities::NullFieldMetaDataMap,
                        ElementGIDs);

            printf("control\n");
            for(int j = 0; j < 2; ++j)
            {
                
                for(int i = 0; i < FieldDef.size(); ++i)
                {

                    pFields[j]->ExtractDataToCoeffs(
                                            FieldDef[i],
                                            FieldData[i],
                                            FieldDef[i]->m_fields[j],
                                            tmp_coeff0);
                    pFields[j]->ExtractDataToCoeffs(
                                            FieldDef1[i],
                                            FieldData1[i],
                                            FieldDef1[i]->m_fields[j],
                                            tmp_coeff1);

                }

                pFields[j]->SetWaveSpace(false);
                pFields[j]->BwdTrans(tmp_coeff0, Vel[j]);
                pFields[j]->BwdTrans(tmp_coeff1, grad[j]);

            }
            Array<OneD, NekDouble> x0(Vel[0].size());
            Array<OneD, NekDouble> x1(Vel[0].size());
            pFields[0]->GetCoords(x0, x1);
            NekDouble e;
            NekDouble xp, yp;
            m_session->LoadParameter("xp",  xp);
            m_session->LoadParameter("yp", yp);
            NekDouble sigma;          
            m_session->LoadParameter("sigma", sigma);

            for (int j = 0; j < Vel[0].size(); ++j)
            {
                // ||U||
                aux[j] = sqrt(Vel[0][j]*Vel[0][j] + Vel[1][j]*Vel[1][j])*sqrt(Vel[0][j]*Vel[0][j] + Vel[1][j]*Vel[1][j]);
                
               
            }
            
            for (int i = 0; i < nvar; ++i)
            {
                Vmath::Vmul(Vel[i].size(), aux, 1, Vel[i], 1,  Vel[i], 1);
                Vmath::Vmul(Vel[i].size(), m_Forcing[i], 1, Vel[i], 1,  m_Forcing[i], 1);
            }
            
            for (int j = 0; j < Vel[0].size(); ++j)
            {
                e = exp(-((x0[j]-xp)*(x0[j]-xp)/sigma + (x1[j]-yp)*(x1[j]-yp)/sigma));
                aux[j] = e*(grad[0][j]*m_Forcing[0][j] + grad[1][j]*m_Forcing[1][j]);
            }
            NekDouble sens = pFields[0]->Integral(aux);
            cout<<sens<< endl;
            printf("control\n");

        }
        else
        {
            /* code */
            Array<OneD, Array<OneD, NekDouble> > tmp(pFields.size());
            for (int i = 0; i < pFields.size(); ++i)
            {
                tmp[i] = pFields[i]->GetPhys();
            }
            Update(pFields, tmp, 0.0);
        }
        
        
    }


    void ForcingBody::Update(
            const Array< OneD, MultiRegions::ExpListSharedPtr > &pFields,
            const Array<OneD, Array<OneD, NekDouble> > &inarray,
            const NekDouble &time)
    {
        LibUtilities::EquationSharedPtr eqn = m_session->GetFunction(
            m_funcName, m_session->GetVariable(0));

        if (!boost::iequals(eqn->GetVlist(), "x y z t"))
        {
            // Coupled forcing
            int nq = pFields[0]->GetNpoints();
            Array<OneD, NekDouble> xc(nq), yc(nq), zc(nq), t(nq, time);
            std::string varstr = "x y z";
            std::vector<Array<OneD, const NekDouble>> fielddata = {
                xc, yc, zc, t};

            for (int i = 0; i < pFields.size(); ++i)
            {
                varstr += " " + m_session->GetVariable(i);
                fielddata.push_back(inarray[i]);
            }

            // Evaluate function
            for (int i = 0; i < m_NumVariable; ++i)
            {
                m_session->GetFunction(m_funcName, m_session->GetVariable(i))->
                    Evaluate(fielddata, m_Forcing[i]);
            }

            return;
        }

        for (int i = 0; i < m_NumVariable; ++i)
        {
            std::string  s_FieldStr   = m_session->GetVariable(i);
            ASSERTL0(m_session->DefinesFunction(m_funcName, s_FieldStr),
                     "Variable '" + s_FieldStr + "' not defined.");
            GetFunction(pFields, m_session, m_funcName, true)->Evaluate(
                s_FieldStr, m_Forcing[i], time);
        }

        // If singleMode or halfMode, transform the forcing term to be in
        // physical space in the plane, but Fourier space in the homogeneous
        // direction
        if (m_transform)
        {
            for (int i = 0; i < m_NumVariable; ++i)
            {
                pFields[0]->HomogeneousFwdTrans(m_Forcing[i], m_Forcing[i]);
            }
        }
    }

    void ForcingBody::GetDislpAcel(NekDouble disp, NekDouble acel)
    {
        y= disp;
        ac=acel;
    }
    void ForcingBody::v_Apply(
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> > &inarray,
            Array<OneD, Array<OneD, NekDouble> > &outarray,
            const NekDouble &time)
    {
        std::string control = m_session->GetSolverInfo("Control");
        if(m_hasTimeFcnScaling)
        {
            Array<OneD, NekDouble>  TimeFcn(1);

            for (int i = 0; i < m_NumVariable; i++)
            {
                EvaluateTimeFunction(time, m_timeFcnEqn, TimeFcn);

                Vmath::Svtvp(outarray[i].size(), TimeFcn[0],
                             m_Forcing[i], 1,
                             outarray[i],  1,
                             outarray[i],  1);
            }
        }
        else if(m_session->DefinesParameter("SteadyF_sens"))
        {

            for (int i = 0; i < m_NumVariable; i++)
            {
                  Vmath::Vadd(outarray[i].size(), outarray[i], 1,
                      baseflow_sens[i], 1, outarray[i], 1);
            }

        }
        else if(boost::iequals(control, "Activate"))
        {
            // cout<<"aqui"<<endl;
            int physTot = fields[0]->GetTotPoints();
            Array<OneD, NekDouble> x0(physTot);
            Array<OneD, NekDouble> x1(physTot);
            fields[0]->GetCoords(x0, x1);
            NekDouble sigma;
            
            m_session->LoadParameter("sigma", sigma);
        
            NekDouble xa, ya, xp, yp;
            m_session->LoadParameter("xp",  xp);
            m_session->LoadParameter("yp", yp);
           
            for(int i = 0; i < m_Forcing.size(); ++i)
            {   
                for (int j = 0; j < m_Forcing[0].size(); ++j)
                {       
                    xa = x0[j];
                    ya = x1[j];
                    xp = xp;
                    yp = yp;
                    
                    outarray[i][j] = exp(-((xa-xp)*(xa-xp)/sigma + (ya-yp)*(ya-yp)/sigma))*m_Forcing[i][j] + outarray[i][j]-ac;

                    
                              
                }
            
            }
            
            
        
            // for (int i = 0; i < m_NumVariable; i++)
            // {
            //     Vmath::Vadd(outarray[i].size(), outarray[i], 1,
            //                 m_Forcing[i], 1, outarray[i], 1);
            // }

        }
        else
        {
            
            Update(fields, inarray, time);

            for (int i = 0; i < m_NumVariable; i++)
            {
                Vmath::Vadd(outarray[i].size(), outarray[i], 1,
                            m_Forcing[i], 1, outarray[i], 1);
            }
        }
    }

}
}
