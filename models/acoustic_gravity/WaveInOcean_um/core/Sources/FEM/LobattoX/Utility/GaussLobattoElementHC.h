#pragma once

 
namespace OndoMathX {

 
    
    namespace
    {
        constexpr Real sqrt_twentyone = 4.582575694955840;

        constexpr Real a = (7.0-sqrt_twentyone)/14.0;
        constexpr Real w[3] ={1 - 3*a + 2*a*a, 4*a - 4*a*a, 2*a*a - a};
        constexpr Real z[2] = {1-a, a};


        constexpr Real dx3[5] = {1.6180339887498944584,
                                 0.61803398874989456946,
                                 8.0901699437494745126,
                                 3.0901699437494727363,
                                 2.2360679774997902491};

        constexpr Real dx4[7] = {2.4819805060619613002,
                                 0.51801949393803370381,
                                 13.513004977448479949,
                                 2.673169155390906937,
                                 1.5275252316519465356,
                                 2.8203283558848517565,
                                 3.4914862437758804603/(16.0/3.0)};
    
    } 

    inline void InterpolationQ21Q41(Real *inQ21, Real * outQ41) //30sec vs 0 sec c'est possible100
    {
        outQ41[0] = inQ21[0];
        outQ41[2] = inQ21[1];
        outQ41[4] = inQ21[2];
        outQ41[5] = inQ21[3];
        outQ41[7] = inQ21[4];
        outQ41[9] = inQ21[5];          

        outQ41[1] = inQ21[0]*w[0] + inQ21[1]*w[1] + inQ21[2]*w[2];
        outQ41[3] = inQ21[0]*w[2] + inQ21[1]*w[1] + inQ21[2]*w[0];
        outQ41[6] = inQ21[3]*w[0] + inQ21[4]*w[1] + inQ21[5]*w[2];
        outQ41[8] = inQ21[3]*w[2] + inQ21[4]*w[1] + inQ21[5]*w[0];
    }

    inline void TransposeInterpolationQ21Q41(Real *inQ41, Real * outQ21) //18sec vs 0 sec
    {
        outQ21[0] = inQ41[0];
        outQ21[1] = inQ41[2];
        outQ21[2] = inQ41[4];
        outQ21[3] = inQ41[5];
        outQ21[4] = inQ41[7];
        outQ21[5] = inQ41[9];  

        outQ21[0] += inQ41[1]*w[0] + inQ41[3]*w[2];
        outQ21[1] += w[1]*(inQ41[1] + inQ41[3]);
        outQ21[2] += inQ41[1]*w[2] + inQ41[3]*w[0];

        outQ21[3] += inQ41[6]*w[0] + inQ41[8]*w[2];
        outQ21[4] += w[1]*(inQ41[6] + inQ41[8]);
        outQ21[5] += inQ41[6]*w[2] + inQ41[8]*w[0];
    }
    
    inline void InterpolationQ12Q14(Real *inQ12, Real * outQ14) //Tested around 30x plus vite
    {
        outQ14[0] = inQ12[0];
        outQ14[1] = inQ12[1];
        outQ14[4] = inQ12[2];
        outQ14[5] = inQ12[3];
        outQ14[8] = inQ12[4];
        outQ14[9] = inQ12[5];          

        outQ14[2] = inQ12[0]*w[0] + inQ12[2]*w[1] + inQ12[4]*w[2];
        outQ14[6] = inQ12[0]*w[2] + inQ12[2]*w[1] + inQ12[4]*w[0];
        outQ14[3] = inQ12[1]*w[0] + inQ12[3]*w[1] + inQ12[5]*w[2];
        outQ14[7] = inQ12[1]*w[2] + inQ12[3]*w[1] + inQ12[5]*w[0];
    }

    inline void TransposeInterpolationQ12Q14(Real *inQ14, Real * outQ12) //Tested much faster 0 sec vs 18 ec
    {
        outQ12[0] = inQ14[0];
        outQ12[1] = inQ14[1];
        outQ12[2] = inQ14[4];
        outQ12[3] = inQ14[5];
        outQ12[4] = inQ14[8];
        outQ12[5] = inQ14[9];        

        outQ12[0] += inQ14[2]*w[0] + inQ14[6]*w[2];
        outQ12[2] += w[1]*(inQ14[2] + inQ14[6]);
        outQ12[4] += inQ14[2]*w[2] + inQ14[6]*w[0];

        outQ12[1] += inQ14[3]*w[0] + inQ14[7]*w[2];
        outQ12[3] += w[1]*(inQ14[3] + inQ14[7]);
        outQ12[5] += inQ14[3]*w[2] + inQ14[7]*w[0];
    }
    
    inline void TransposeInterpolationQ2Q4Bndy(Real * inQ4,Real *outQ2) 
    {
        outQ2[4] = 0.0;

        outQ2[0] = inQ4[0];
        outQ2[1] = inQ4[2];
        outQ2[2] = inQ4[4];
        outQ2[3] = inQ4[10];
        outQ2[5] = inQ4[14];
        outQ2[6] = inQ4[20];                     
        outQ2[7] = inQ4[22];
        outQ2[8] = inQ4[24];

        outQ2[0] += (inQ4[1]+inQ4[5])*w[0] + (inQ4[3]+inQ4[15])*w[2];
        outQ2[1] += (inQ4[1]+inQ4[3])*w[1];
        outQ2[2] += (inQ4[3]+inQ4[9])*w[0] + (inQ4[1]+inQ4[19])*w[2];
        outQ2[3] += (inQ4[5]+inQ4[15])*w[1];
        outQ2[5] += (inQ4[9]+inQ4[19])*w[1];
        outQ2[6] += (inQ4[15]+inQ4[21])*w[0] + (inQ4[5]+inQ4[23])*w[2];
        outQ2[7] += (inQ4[21]+inQ4[23])*w[1];
        outQ2[8] += (inQ4[19]+inQ4[23])*w[0] + (inQ4[9]+inQ4[21])*w[2];
    }

    inline void InterpolationQ2Q4Bndy(Real *inQ2, Real * outQ4)   
    {
        outQ4[6]  = 0.0;
        outQ4[7]  = 0.0;
        outQ4[8]  = 0.0;
        outQ4[11] = 0.0;
        outQ4[12] = 0.0;
        outQ4[13] = 0.0;                  
        outQ4[16] = 0.0;
        outQ4[17] = 0.0;
        outQ4[18] = 0.0;

        outQ4[0]  = inQ2[0];
        outQ4[2]  = inQ2[1];
        outQ4[4]  = inQ2[2];
        outQ4[10] = inQ2[3];
        outQ4[14] = inQ2[5];
        outQ4[20] = inQ2[6];                     
        outQ4[22] = inQ2[7];
        outQ4[24] = inQ2[8];

        outQ4[1]  = inQ2[0]*w[0] + inQ2[1]*w[1] + inQ2[2]*w[2];
        outQ4[3]  = inQ2[2]*w[0] + inQ2[1]*w[1] + inQ2[0]*w[2];
        outQ4[5]  = inQ2[0]*w[0] + inQ2[3]*w[1] + inQ2[6]*w[2];
        outQ4[15] = inQ2[6]*w[0] + inQ2[3]*w[1] + inQ2[0]*w[2];
        outQ4[9]  = inQ2[2]*w[0] + inQ2[5]*w[1] + inQ2[8]*w[2];
        outQ4[19] = inQ2[8]*w[0] + inQ2[5]*w[1] + inQ2[2]*w[2];
        outQ4[21] = inQ2[6]*w[0] + inQ2[7]*w[1] + inQ2[8]*w[2];
        outQ4[23] = inQ2[8]*w[0] + inQ2[7]*w[1] + inQ2[6]*w[2];
    }


    inline void InterpolationQ14Q4(Real *inQ14, Real * outQ4) //Tested 30x plus vite
    {
        outQ4[0]  = inQ14[0];
        outQ4[4]  = inQ14[1];
        outQ4[5]  = inQ14[2];
        outQ4[9]  = inQ14[3];
        outQ4[10] = inQ14[4];
        outQ4[14] = inQ14[5];                     
        outQ4[15] = inQ14[6];
        outQ4[19] = inQ14[7];
        outQ4[20] = inQ14[8];
        outQ4[24] = inQ14[9];      

        outQ4[1]  = inQ14[0]*z[0] + inQ14[1]*z[1];
        outQ4[2]  = 0.5*(inQ14[0] + inQ14[1]);
        outQ4[3]  = inQ14[0]*z[1] + inQ14[1]*z[0];

        outQ4[6]  = inQ14[2]*z[0] + inQ14[3]*z[1];
        outQ4[7]  = 0.5*(inQ14[2] + inQ14[3]);
        outQ4[8]  = inQ14[2]*z[1] + inQ14[3]*z[0];

        outQ4[11] = inQ14[4]*z[0] + inQ14[5]*z[1];
        outQ4[12] = 0.5*(inQ14[4] + inQ14[5]);
        outQ4[13] = inQ14[4]*z[1] + inQ14[5]*z[0];


        outQ4[16] = inQ14[6]*z[0] + inQ14[7]*z[1];
        outQ4[17] = 0.5*(inQ14[6] + inQ14[7]);
        outQ4[18] = inQ14[6]*z[1] + inQ14[7]*z[0];


        outQ4[21] = inQ14[8]*z[0] + inQ14[9]*z[1];
        outQ4[22] = 0.5*(inQ14[8] + inQ14[9]);
        outQ4[23] = inQ14[8]*z[1] + inQ14[9]*z[0];


    }

    inline void TransposeInterpolationQ14Q4(Real *inQ4, Real * outQ14) //Tested 18x plus vite
    {
        outQ14[0]  = inQ4[0]; 
        outQ14[1]  = inQ4[4];
        outQ14[2]  = inQ4[5];
        outQ14[3]  = inQ4[9];
        outQ14[4]  = inQ4[10];
        outQ14[5]  = inQ4[14];                     
        outQ14[6]  = inQ4[15];
        outQ14[7]  = inQ4[19];
        outQ14[8]  = inQ4[20];
        outQ14[9]  = inQ4[24];      

        outQ14[0]  += inQ4[1]*z[0] + 0.5*inQ4[2] + inQ4[3]*z[1];
        outQ14[1]  += inQ4[1]*z[1] + 0.5*inQ4[2] + inQ4[3]*z[0];

        outQ14[2]  += inQ4[6]*z[0] + 0.5*inQ4[7] + inQ4[8]*z[1];
        outQ14[3]  += inQ4[6]*z[1] + 0.5*inQ4[7] + inQ4[8]*z[0];

        outQ14[4]  += inQ4[11]*z[0] + 0.5*inQ4[12] + inQ4[13]*z[1];
        outQ14[5]  += inQ4[11]*z[1] + 0.5*inQ4[12] + inQ4[13]*z[0];

        outQ14[6]  += inQ4[16]*z[0] + 0.5*inQ4[17] + inQ4[18]*z[1];
        outQ14[7]  += inQ4[16]*z[1] + 0.5*inQ4[17] + inQ4[18]*z[0];

        outQ14[8]  += inQ4[21]*z[0] + 0.5*inQ4[22] + inQ4[23]*z[1];
        outQ14[9]  += inQ4[21]*z[1] + 0.5*inQ4[22] + inQ4[23]*z[0];

    }



    inline void InterpolationQ41Q4(Real *inQ41, Real * outQ4)  //Tested -> 25x plus vite
    {
        outQ4[0]  = inQ41[0];
        outQ4[1]  = inQ41[1];
        outQ4[2]  = inQ41[2];
        outQ4[3]  = inQ41[3];
        outQ4[4]  = inQ41[4];
        outQ4[20] = inQ41[5];                     
        outQ4[21] = inQ41[6];
        outQ4[22] = inQ41[7];
        outQ4[23] = inQ41[8];
        outQ4[24] = inQ41[9];      

        outQ4[5]  = inQ41[0]*z[0] + inQ41[5]*z[1];
        outQ4[10] = 0.5*(inQ41[0] + inQ41[5]);
        outQ4[15] = inQ41[0]*z[1] + inQ41[5]*z[0];

        outQ4[6]  = inQ41[1]*z[0] + inQ41[6]*z[1];
        outQ4[11] = 0.5*(inQ41[1] + inQ41[6]);
        outQ4[16] = inQ41[1]*z[1] + inQ41[6]*z[0];

        outQ4[7]  = inQ41[2]*z[0] + inQ41[7]*z[1];
        outQ4[12] = 0.5*(inQ41[2] + inQ41[7]);
        outQ4[17] = inQ41[2]*z[1] + inQ41[7]*z[0];

        outQ4[8]  = inQ41[3]*z[0] + inQ41[8]*z[1];
        outQ4[13] = 0.5*(inQ41[3] + inQ41[8]);
        outQ4[18] = inQ41[3]*z[1] + inQ41[8]*z[0];

        outQ4[9]  = inQ41[4]*z[0] + inQ41[9]*z[1];
        outQ4[14] = 0.5*(inQ41[4] + inQ41[9]);
        outQ4[19] = inQ41[4]*z[1] + inQ41[9]*z[0];
    }
    
    inline void TransposeInterpolationQ41Q4(Real *inQ4, Real * outQ41)  //Tested -> 11-12x plus vite
    {
        outQ41[0] = inQ4[0];
        outQ41[1] = inQ4[1];
        outQ41[2] = inQ4[2];
        outQ41[3] = inQ4[3];
        outQ41[4] = inQ4[4];
        outQ41[5] = inQ4[20];                     
        outQ41[6] = inQ4[21];
        outQ41[7] = inQ4[22];
        outQ41[8] = inQ4[23];
        outQ41[9] = inQ4[24];    

        outQ41[0]+=z[0]*inQ4[5]+inQ4[15]*z[1]+0.5*inQ4[10];
        outQ41[5]+=z[1]*inQ4[5]+inQ4[15]*z[0]+0.5*inQ4[10];

        outQ41[1]+=z[0]*inQ4[6]+inQ4[16]*z[1]+0.5*inQ4[11];
        outQ41[6]+=z[1]*inQ4[6]+inQ4[16]*z[0]+0.5*inQ4[11];


        outQ41[2]+=z[0]*inQ4[7]+inQ4[17]*z[1]+0.5*inQ4[12];
        outQ41[7]+=z[1]*inQ4[7]+inQ4[17]*z[0]+0.5*inQ4[12];

        outQ41[3]+=z[0]*inQ4[8]+inQ4[18]*z[1]+0.5*inQ4[13];
        outQ41[8]+=z[1]*inQ4[8]+inQ4[18]*z[0]+0.5*inQ4[13];

        outQ41[4]+=z[0]*inQ4[9]+inQ4[19]*z[1]+0.5*inQ4[14];
        outQ41[9]+=z[1]*inQ4[9]+inQ4[19]*z[0]+0.5*inQ4[14];

    }

    inline void InterpolationQ22Q24( Real *inQ22,  Real *outQ24  )
    {
        outQ24[0] = inQ22[0];
        outQ24[1] = inQ22[1];
        outQ24[2] = inQ22[2];

        outQ24[6] = inQ22[3];
        outQ24[7] = inQ22[4];
        outQ24[8] = inQ22[5];

        outQ24[12] = inQ22[6];
        outQ24[13] = inQ22[7];
        outQ24[14] = inQ22[8];

        outQ24[3]  = inQ22[0]*w[0] + inQ22[3]*w[1] + inQ22[6]*w[2];
        outQ24[9]  = inQ22[0]*w[2] + inQ22[3]*w[1] + inQ22[6]*w[0];

        outQ24[4]  = inQ22[1]*w[0] + inQ22[4]*w[1] + inQ22[7]*w[2];
        outQ24[10] = inQ22[1]*w[2] + inQ22[4]*w[1] + inQ22[7]*w[0];

        outQ24[5]  = inQ22[2]*w[0] + inQ22[5]*w[1] + inQ22[8]*w[2];
        outQ24[11] = inQ22[2]*w[2] + inQ22[5]*w[1] + inQ22[8]*w[0];    

    }

    inline void TransposeInterpolationQ22Q24( Real *inQ24,  Real *outQ22  )
    {
        outQ22[0] = inQ24[0];
        outQ22[1] = inQ24[1];
        outQ22[2] = inQ24[2];

        outQ22[3] = inQ24[6];
        outQ22[4] = inQ24[7];
        outQ22[5] = inQ24[8];

        outQ22[6] = inQ24[12];
        outQ22[7] = inQ24[13];
        outQ22[8] = inQ24[14];

        outQ22[0] += inQ24[3]*w[0]   + inQ24[9]*w[2];
        outQ22[3] += w[1]* (inQ24[3] + inQ24[9]);
        outQ22[6] += inQ24[3]*w[2]   + inQ24[9]*w[0];

        outQ22[1] += inQ24[4]*w[0]   + inQ24[10]*w[2];
        outQ22[4] += w[1]* (inQ24[4] + inQ24[10]);
        outQ22[7] += inQ24[4]*w[2]   + inQ24[10]*w[0];

        outQ22[2] += inQ24[5]*w[0]   + inQ24[11]*w[2];
        outQ22[5] += w[1]* (inQ24[5] + inQ24[11]);
        outQ22[8] += inQ24[5]*w[2]   + inQ24[11]*w[0];

    }


    inline void InterpolationQ24Q4 (Real *inQ24, Real *outQ4)
    {
        outQ4[0] = inQ24[0];
        outQ4[2] = inQ24[1];
        outQ4[4] = inQ24[2];
        outQ4[5] = inQ24[3];
        outQ4[7] = inQ24[4];
        outQ4[9] = inQ24[5];
        outQ4[10] = inQ24[6];
        outQ4[12] = inQ24[7];
        outQ4[14] = inQ24[8];
        outQ4[15] = inQ24[9];
        outQ4[17] = inQ24[10];
        outQ4[19] = inQ24[11];
        outQ4[20] = inQ24[12];
        outQ4[22] = inQ24[13];
        outQ4[24] = inQ24[14];

        outQ4[1] = inQ24[0]*w[0] + inQ24[1]*w[1] + inQ24[2]*w[2];
        outQ4[3] = inQ24[0]*w[2] + inQ24[1]*w[1] + inQ24[2]*w[0];

        outQ4[6] = inQ24[3]*w[0] + inQ24[4]*w[1] + inQ24[5]*w[2];
        outQ4[8] = inQ24[3]*w[2] + inQ24[4]*w[1] + inQ24[5]*w[0];

        outQ4[11] = inQ24[6]*w[0] + inQ24[7]*w[1] + inQ24[8]*w[2];
        outQ4[13] = inQ24[6]*w[2] + inQ24[7]*w[1] + inQ24[8]*w[0];

        outQ4[16] = inQ24[9]*w[0] + inQ24[10]*w[1] + inQ24[11]*w[2];
        outQ4[18] = inQ24[9]*w[2] + inQ24[10]*w[1] + inQ24[11]*w[0];

        outQ4[21] = inQ24[12]*w[0] + inQ24[13]*w[1] + inQ24[14]*w[2];
        outQ4[23] = inQ24[12]*w[2] + inQ24[13]*w[1] + inQ24[14]*w[0];

    }


    inline void TransposeInterpolationQ24Q4( Real *inQ4,  Real *outQ24  )
    {
        outQ24[0] = inQ4[0];
        outQ24[1] = inQ4[2];
        outQ24[2] = inQ4[4];
        outQ24[3] = inQ4[5];
        outQ24[4] = inQ4[7];
        outQ24[5] = inQ4[9];
        outQ24[6] = inQ4[10];
        outQ24[7] = inQ4[12];
        outQ24[8] = inQ4[14];
        outQ24[9]  = inQ4[15];
        outQ24[10] = inQ4[17];
        outQ24[11] = inQ4[19];
        outQ24[12] = inQ4[20];
        outQ24[13] = inQ4[22];
        outQ24[14] = inQ4[24];

        outQ24[0] += inQ4[1]*w[0]   + inQ4[3]*w[2];
        outQ24[1] += w[1]* (inQ4[1] + inQ4[3]);
        outQ24[2] += inQ4[1]*w[2]   + inQ4[3]*w[0];

        outQ24[3] += inQ4[6]*w[0]   + inQ4[8]*w[2];
        outQ24[4] += w[1]* (inQ4[6] + inQ4[8]);
        outQ24[5] += inQ4[6]*w[2]   + inQ4[8]*w[0];

        outQ24[6] += inQ4[11]*w[0]   + inQ4[13]*w[2];
        outQ24[7] += w[1]* (inQ4[11] + inQ4[13]);
        outQ24[8] += inQ4[11]*w[2]   + inQ4[13]*w[0];

        outQ24[9]  += inQ4[16]*w[0]   + inQ4[18]*w[2];
        outQ24[10] += w[1]* (inQ4[16] + inQ4[18]);
        outQ24[11] += inQ4[16]*w[2]   + inQ4[18]*w[0];

        outQ24[12] += inQ4[21]*w[0]   + inQ4[23]*w[2];
        outQ24[13] += w[1]* (inQ4[21] + inQ4[23]);
        outQ24[14] += inQ4[21]*w[2]   + inQ4[23]*w[0];

    }


    /*----------------Gradient---------------*/

    inline void DxQ4Q4(Real *inQ4,  Real *outQ4, Index stride=1)
    {
        Real a = *inQ4;                 inQ4+=stride;
        Real b = *inQ4;                 inQ4+=stride;
        Real c = (16.0/3.0)*(*inQ4);    inQ4+=stride;
        Real d = *inQ4;                 inQ4+=stride;
        Real e = *inQ4; 

        *outQ4 = dx4[2]*b-10*a-c+dx4[5]*d-e;            outQ4+=stride; 
        *outQ4 = dx4[1]*e-dx4[0]*a+dx4[6]*c-dx4[4]*d;   outQ4+=stride; 
        *outQ4 = 0.75*(a-e)+dx4[3]*(d-b);               outQ4+=stride; 
        *outQ4 = dx4[4]*b-dx4[1]*a-dx4[6]*c+dx4[0]*e;   outQ4+=stride; 
        *outQ4 = a-dx4[5]*b+c-dx4[2]*d+10*e;    
    }

    inline void DxQ3Q3(Real *inQ3,  Real *outQ3, Index stride=1)
    {
        Real a = *inQ3; inQ3+=stride;
        Real b = *inQ3; inQ3+=stride;
        Real c = *inQ3; inQ3+=stride;
        Real d = *inQ3; inQ3+=stride;

        *outQ3 = dx3[2]*b-dx3[3]*c+d-6*a;      outQ3+=stride; 
        *outQ3 = dx3[4]*c-dx3[0]*a-dx3[1]*d;   outQ3+=stride;
        *outQ3 = dx3[1]*a-dx3[4]*b+dx3[0]*d;   outQ3+=stride;
        *outQ3 = dx3[3]*b-dx3[2]*c-a+6*d;       
    }

    inline void DxQ3Q333(Real *inQ3,  Real *outQ333 )
    {
        for (Index i=0;i<16;++i)
        {
            DxQ3Q3(inQ3,outQ333);
            inQ3+=4;
            outQ333+=4;
        }
    }

    inline void DxQ4Q444(Real *inQ4,  Real *outQ444 )
    {
        for (Index i=0;i<25;++i)
        {
            DxQ4Q4(inQ4,outQ444);
            inQ4+=5;
            outQ444+=5;
        }
    }

    inline void DyQ4Q444(Real *inQ4,  Real *outQ444 )
    {
        for (Index k=0;k<5;++k)
        {
            for (Index i=0;i<5;++i)
            {
                DxQ4Q4(inQ4,outQ444,5);
                inQ4+=1;
                outQ444+=1;
            }
            inQ4+=20;
            outQ444+=20;
        }
    }


    inline void DxQ2Q12( Real *inQ2,  Real *outQ12 )
    {
        outQ12[0] = -3*inQ2[0] + 4*inQ2[1]   - inQ2[2];
        outQ12[1] =    inQ2[0] - 4*inQ2[1] + 3*inQ2[2];
        outQ12[2] = -3*inQ2[3] + 4*inQ2[4]   - inQ2[5];
        outQ12[3] =    inQ2[3] - 4*inQ2[4] + 3*inQ2[5];
        outQ12[4] = -3*inQ2[6] + 4*inQ2[7]   - inQ2[8];
        outQ12[5] =    inQ2[6] - 4*inQ2[7] + 3*inQ2[8];
    }
    

    inline void DyQ2Q21( Real *inQ2,  Real *outQ21 )
    {
        outQ21[0] = -3*inQ2[0] + 4*inQ2[3]   - inQ2[6];
        outQ21[1] = -3*inQ2[1] + 4*inQ2[4]   - inQ2[7];
        outQ21[2] = -3*inQ2[2] + 4*inQ2[5]   - inQ2[8];
        outQ21[3] =    inQ2[0] - 4*inQ2[3] + 3*inQ2[6];
        outQ21[4] =    inQ2[1] - 4*inQ2[4] + 3*inQ2[7];
        outQ21[5] =    inQ2[2] - 4*inQ2[5] + 3*inQ2[8];
    }

    inline void TransposeDxQ2Q12( Real *inQ12,  Real *outQ2 )
    {
        outQ2[0] =  -3*inQ12[0] +   inQ12[1];
        outQ2[1] =   4*inQ12[0] - 4*inQ12[1];
        outQ2[2] =    -inQ12[0] + 3*inQ12[1];

        outQ2[3] =  -3*inQ12[2] +   inQ12[3];
        outQ2[4] =   4*inQ12[2] - 4*inQ12[3];
        outQ2[5] =    -inQ12[2] + 3*inQ12[3];

        outQ2[6] =  -3*inQ12[4] +   inQ12[5];
        outQ2[7] =   4*inQ12[4] - 4*inQ12[5];
        outQ2[8] =    -inQ12[4] + 3*inQ12[5];

    }


    inline void TransposeDyQ2Q21( Real *inQ21,  Real *outQ2 )
    {
        outQ2[0] = -3*inQ21[0] +   inQ21[3];
        outQ2[3] =  4*inQ21[0] - 4*inQ21[3];
        outQ2[6] =   -inQ21[0] + 3*inQ21[3];
        
        outQ2[1] = -3*inQ21[1] +   inQ21[4];
        outQ2[4] =  4*inQ21[1] - 4*inQ21[4];
        outQ2[7] =   -inQ21[1] + 3*inQ21[4];

        outQ2[2] = -3*inQ21[2] +   inQ21[5];
        outQ2[5] =  4*inQ21[2] - 4*inQ21[5];
        outQ2[8] =   -inQ21[2] + 3*inQ21[5];

    }

    /*------------- 3D Gradient-------------*/ 
    inline void DxQ2Q122( Real *inQ2,  Real *outQ122 )
    {
        DxQ2Q12( inQ2,  outQ122 );
        DxQ2Q12( inQ2+9, outQ122+6 );
        DxQ2Q12( inQ2+18, outQ122+12 );
    }

    inline void TransposeDxQ2Q122( Real *inQ122,  Real *outQ2 )
    {
        TransposeDxQ2Q12( inQ122,  outQ2 );
        TransposeDxQ2Q12( inQ122+6, outQ2+9 );
        TransposeDxQ2Q12( inQ122+12, outQ2+18 );
    }


    inline void DyQ2Q212( Real *inQ2,  Real *outQ212 )
    {
        DyQ2Q21( inQ2,  outQ212 );
        DyQ2Q21( inQ2+9, outQ212+6 );
        DyQ2Q21( inQ2+18, outQ212+12 );
    }

    inline void TransposeDyQ2Q212( Real *inQ212,  Real *outQ2 )
    {
        TransposeDyQ2Q21( inQ212,  outQ2 );
        TransposeDyQ2Q21( inQ212+6, outQ2+9 );
        TransposeDyQ2Q21( inQ212+12, outQ2+18 );
    }

    inline void DzQ2Q221 ( Real *inQ2,  Real *outQ221 )
    {
        outQ221[0]  = -3*inQ2[0] + 4*inQ2[9]  -   inQ2[18];
        outQ221[9]  =    inQ2[0] - 4*inQ2[9]  + 3*inQ2[18];
        outQ221[1]  = -3*inQ2[1] + 4*inQ2[10] -   inQ2[19];
        outQ221[10] =    inQ2[1] - 4*inQ2[10] + 3*inQ2[19];
        outQ221[2]  = -3*inQ2[2] + 4*inQ2[11] -   inQ2[20];
        outQ221[11] =    inQ2[2] - 4*inQ2[11] + 3*inQ2[20];
        outQ221[3]  = -3*inQ2[3] + 4*inQ2[12] -   inQ2[21];
        outQ221[12] =    inQ2[3] - 4*inQ2[12] + 3*inQ2[21];
        outQ221[4]  = -3*inQ2[4] + 4*inQ2[13] -   inQ2[22];
        outQ221[13] =    inQ2[4] - 4*inQ2[13] + 3*inQ2[22];
        outQ221[5]  = -3*inQ2[5] + 4*inQ2[14] -   inQ2[23];
        outQ221[14] =    inQ2[5] - 4*inQ2[14] + 3*inQ2[23];
        outQ221[6]  = -3*inQ2[6] + 4*inQ2[15] -   inQ2[24];
        outQ221[15] =    inQ2[6] - 4*inQ2[15] + 3*inQ2[24];
        outQ221[7]  = -3*inQ2[7] + 4*inQ2[16] -   inQ2[25];
        outQ221[16] =    inQ2[7] - 4*inQ2[16] + 3*inQ2[25];
        outQ221[8]  = -3*inQ2[8] + 4*inQ2[17] -   inQ2[26];
        outQ221[17] =    inQ2[8] - 4*inQ2[17] + 3*inQ2[26];
    }

    inline void TransposeDzQ2Q221 ( Real *inQ221,  Real *outQ2 )
    {
        outQ2[0]  = -3*inQ221[0] +   inQ221[9];
        outQ2[9]  = 4*(inQ221[0] -   inQ221[9]); 
        outQ2[18] =   -inQ221[0] + 3*inQ221[9];

        outQ2[1]  = -3*inQ221[1] +   inQ221[10];
        outQ2[10] = 4*(inQ221[1] -   inQ221[10]); 
        outQ2[19] =   -inQ221[1] + 3*inQ221[10];

        outQ2[2]  = -3*inQ221[2] +   inQ221[11];
        outQ2[11] = 4*(inQ221[2] -   inQ221[11]); 
        outQ2[20] =   -inQ221[2] + 3*inQ221[11];

        outQ2[3]  = -3*inQ221[3] +   inQ221[12];
        outQ2[12] = 4*(inQ221[3] -   inQ221[12]); 
        outQ2[21] =   -inQ221[3] + 3*inQ221[12];

        outQ2[4]  = -3*inQ221[4] +   inQ221[13];
        outQ2[13] = 4*(inQ221[4] -   inQ221[13]); 
        outQ2[22] =   -inQ221[4] + 3*inQ221[13];

        outQ2[5]  = -3*inQ221[5] +   inQ221[14];
        outQ2[14] = 4*(inQ221[5] -   inQ221[14]); 
        outQ2[23] =   -inQ221[5] + 3*inQ221[14];

        outQ2[6]  = -3*inQ221[6] +   inQ221[15];
        outQ2[15] = 4*(inQ221[6] -   inQ221[15]); 
        outQ2[24] =   -inQ221[6] + 3*inQ221[15];

        outQ2[7]  = -3*inQ221[7] +   inQ221[16];
        outQ2[16] = 4*(inQ221[7] -   inQ221[16]); 
        outQ2[25] =   -inQ221[7] + 3*inQ221[16];

        outQ2[8]  = -3*inQ221[8] +   inQ221[17];
        outQ2[17] = 4*(inQ221[8] -   inQ221[17]); 
        outQ2[26] =   -inQ221[8] + 3*inQ221[17];
    }



    /*--------------------------------------------------------*/
    /*-------------------- 3D Interpolation ------------------*/
    /*--------------------------------------------------------*/    
    // When you pass an array as a parameter to a function it decays 
    //to a pointer to the first element of the array. 
    //So there is normally never a need to pass a pointer to an array.

    /*--------- Boundary Inerpolation Q2->Q442->Q4 ---------*/

    inline void InterpolationQ222Q444Bndy (Real *inQ2,Real *outQ4)
    {
         Real tempQ24[15], tempQ442[75];
        // 1st layer
        InterpolationQ22Q24( inQ2, tempQ24 );
        InterpolationQ24Q4( tempQ24, tempQ442);
        // 2nd layer
        InterpolationQ2Q4Bndy( inQ2+9, tempQ442+25 );
        // 3rd layer
        InterpolationQ22Q24( inQ2+18, tempQ24 );
        InterpolationQ24Q4( tempQ24, tempQ442+50);

        for (Index i = 0; i < 25; i++)
        {
            outQ4[i] = tempQ442[i];
            outQ4[i+50] = tempQ442[i+25];
            outQ4[i+100] = tempQ442[i+50];
        }
        
        outQ4[25] = tempQ442[0]*w[0] + tempQ442[25]*w[1] + tempQ442[50]*w[2];
        outQ4[26] = tempQ442[1]*w[0] + tempQ442[26]*w[1] + tempQ442[51]*w[2];
        outQ4[27] = tempQ442[2]*w[0] + tempQ442[27]*w[1] + tempQ442[52]*w[2];
        outQ4[28] = tempQ442[3]*w[0] + tempQ442[28]*w[1] + tempQ442[53]*w[2];
        outQ4[29] = tempQ442[4]*w[0] + tempQ442[29]*w[1] + tempQ442[54]*w[2];
        outQ4[30] = tempQ442[5]*w[0] + tempQ442[30]*w[1] + tempQ442[55]*w[2];
        outQ4[31] = 0.0;
        outQ4[32] = 0.0;
        outQ4[33] = 0.0;
        outQ4[34] = tempQ442[9]*w[0] + tempQ442[34]*w[1] + tempQ442[59]*w[2];
        outQ4[35] = tempQ442[10]*w[0] + tempQ442[35]*w[1] + tempQ442[60]*w[2];
        outQ4[36] = 0.0;
        outQ4[37] = 0.0;
        outQ4[38] = 0.0;
        outQ4[39] = tempQ442[14]*w[0] + tempQ442[39]*w[1] + tempQ442[64]*w[2];
        outQ4[40] = tempQ442[15]*w[0] + tempQ442[40]*w[1] + tempQ442[65]*w[2];
        outQ4[41] = 0.0;
        outQ4[42] = 0.0;
        outQ4[43] = 0.0;
        outQ4[44] = tempQ442[19]*w[0] + tempQ442[44]*w[1] + tempQ442[69]*w[2];
        outQ4[45] = tempQ442[20]*w[0] + tempQ442[45]*w[1] + tempQ442[70]*w[2];
        outQ4[46] = tempQ442[21]*w[0] + tempQ442[46]*w[1] + tempQ442[71]*w[2];
        outQ4[47] = tempQ442[22]*w[0] + tempQ442[47]*w[1] + tempQ442[72]*w[2];
        outQ4[48] = tempQ442[23]*w[0] + tempQ442[48]*w[1] + tempQ442[73]*w[2];
        outQ4[49] = tempQ442[24]*w[0] + tempQ442[49]*w[1] + tempQ442[74]*w[2];

     

        outQ4[75] = tempQ442[0]*w[2] + tempQ442[25]*w[1] + tempQ442[50]*w[0];
        outQ4[76] = tempQ442[1]*w[2] + tempQ442[26]*w[1] + tempQ442[51]*w[0];
        outQ4[77] = tempQ442[2]*w[2] + tempQ442[27]*w[1] + tempQ442[52]*w[0];
        outQ4[78] = tempQ442[3]*w[2] + tempQ442[28]*w[1] + tempQ442[53]*w[0];
        outQ4[79] = tempQ442[4]*w[2] + tempQ442[29]*w[1] + tempQ442[54]*w[0];
        outQ4[80] = tempQ442[5]*w[2] + tempQ442[30]*w[1] + tempQ442[55]*w[0];
        outQ4[81] = 0.0;
        outQ4[82] = 0.0;
        outQ4[83] = 0.0;
        outQ4[84] = tempQ442[9]*w[2] + tempQ442[34]*w[1] + tempQ442[59]*w[0];
        outQ4[85] = tempQ442[10]*w[2] + tempQ442[35]*w[1] + tempQ442[60]*w[0];
        outQ4[86] = 0.0;
        outQ4[87] = 0.0;
        outQ4[88] = 0.0;
        outQ4[89] = tempQ442[14]*w[2] + tempQ442[39]*w[1] + tempQ442[64]*w[0];
        outQ4[90] = tempQ442[15]*w[2] + tempQ442[40]*w[1] + tempQ442[65]*w[0];
        outQ4[91] = 0.0;
        outQ4[92] = 0.0;
        outQ4[93] = 0.0;
        outQ4[94] = tempQ442[19]*w[2] + tempQ442[44]*w[1] + tempQ442[69]*w[0];
        outQ4[95] = tempQ442[20]*w[2] + tempQ442[45]*w[1] + tempQ442[70]*w[0];
        outQ4[96] = tempQ442[21]*w[2] + tempQ442[46]*w[1] + tempQ442[71]*w[0];
        outQ4[97] = tempQ442[22]*w[2] + tempQ442[47]*w[1] + tempQ442[72]*w[0];
        outQ4[98] = tempQ442[23]*w[2] + tempQ442[48]*w[1] + tempQ442[73]*w[0];
        outQ4[99] = tempQ442[24]*w[2] + tempQ442[49]*w[1] + tempQ442[74]*w[0];

    }



    inline void TransposeInterpolationQ222Q444Bndy (Real *inQ4,Real *outQ2)
    {
        Real tempQ24[15], tempQ442[75];

        for (Index i = 0; i < 25; i++)
        {
            tempQ442[i] = inQ4[i];
            tempQ442[i+25] = inQ4[i+50];
            tempQ442[i+50] = inQ4[i+100];
        }
        


        tempQ442[0]  +=  inQ4[25]*w[0]  + inQ4[75]*w[2];
        tempQ442[25] +=  w[1]*(inQ4[25] + inQ4[75]);
        tempQ442[50] +=  inQ4[25]*w[2]  + inQ4[75]*w[0];

        tempQ442[1]  +=  inQ4[26]*w[0]  + inQ4[76]*w[2];
        tempQ442[26] +=  w[1]*(inQ4[26] + inQ4[76]);
        tempQ442[51] +=  inQ4[26]*w[2]  + inQ4[76]*w[0];

        tempQ442[2]  +=  inQ4[27]*w[0]  + inQ4[77]*w[2];
        tempQ442[27] +=  w[1]*(inQ4[27] + inQ4[77]);
        tempQ442[52] +=  inQ4[27]*w[2]  + inQ4[77]*w[0];

        tempQ442[3]  +=  inQ4[28]*w[0]  + inQ4[78]*w[2];
        tempQ442[28] +=  w[1]*(inQ4[28] + inQ4[78]);
        tempQ442[53] +=  inQ4[28]*w[2]  + inQ4[78]*w[0];

        tempQ442[4]  +=  inQ4[29]*w[0]  + inQ4[79]*w[2];
        tempQ442[29] +=  w[1]*(inQ4[29] + inQ4[79]);
        tempQ442[54] +=  inQ4[29]*w[2]  + inQ4[79]*w[0];

        tempQ442[5]  +=  inQ4[30]*w[0]  + inQ4[80]*w[2];
        tempQ442[30] +=  w[1]*(inQ4[30] + inQ4[80]);
        tempQ442[55] +=  inQ4[30]*w[2]  + inQ4[80]*w[0];

        tempQ442[31]  = 0.0;
        tempQ442[32]  = 0.0;
        tempQ442[33]  = 0.0;

        tempQ442[9]  +=  inQ4[34]*w[0]  + inQ4[84]*w[2];
        tempQ442[34] +=  w[1]*(inQ4[34] + inQ4[84]);
        tempQ442[59] +=  inQ4[34]*w[2]  + inQ4[84]*w[0];
        tempQ442[10] +=  inQ4[35]*w[0]  + inQ4[85]*w[2];
        tempQ442[35] +=  w[1]*(inQ4[35] + inQ4[85]);
        tempQ442[60] +=  inQ4[35]*w[2]  + inQ4[85]*w[0];

        tempQ442[36]  = 0.0;
        tempQ442[37]  = 0.0;
        tempQ442[38]  = 0.0;

        tempQ442[14] +=  inQ4[39]*w[0]  + inQ4[89]*w[2];
        tempQ442[39] +=  w[1]*(inQ4[39] + inQ4[89]);
        tempQ442[64] +=  inQ4[39]*w[2]  + inQ4[89]*w[0];
        tempQ442[15] +=  inQ4[40]*w[0]  + inQ4[90]*w[2];
        tempQ442[40] +=  w[1]*(inQ4[40] + inQ4[90]);
        tempQ442[65] +=  inQ4[40]*w[2]  + inQ4[90]*w[0];

        tempQ442[41]  = 0.0;
        tempQ442[42]  = 0.0;
        tempQ442[43]  = 0.0;

        tempQ442[19] +=  inQ4[44]*w[0]  + inQ4[94]*w[2];
        tempQ442[44] +=  w[1]*(inQ4[44] + inQ4[94]);
        tempQ442[69] +=  inQ4[44]*w[2]  + inQ4[94]*w[0];
        tempQ442[20] +=  inQ4[45]*w[0]  + inQ4[95]*w[2];
        tempQ442[45] +=  w[1]*(inQ4[45] + inQ4[95]);
        tempQ442[70] +=  inQ4[45]*w[2]  + inQ4[95]*w[0];

        tempQ442[21] +=  inQ4[46]*w[0]  + inQ4[96]*w[2];
        tempQ442[46] +=  w[1]*(inQ4[46] + inQ4[96]);
        tempQ442[71] +=  inQ4[46]*w[2]  + inQ4[96]*w[0];

        tempQ442[22] +=  inQ4[47]*w[0]  + inQ4[97]*w[2];
        tempQ442[47] +=  w[1]*(inQ4[47] + inQ4[97]);
        tempQ442[72] +=  inQ4[47]*w[2]  + inQ4[97]*w[0];

        tempQ442[23] +=  inQ4[48]*w[0]  + inQ4[98]*w[2];
        tempQ442[48] +=  w[1]*(inQ4[48] + inQ4[98]);
        tempQ442[73] +=  inQ4[48]*w[2]  + inQ4[98]*w[0];

        tempQ442[24] +=  inQ4[49]*w[0]  + inQ4[99]*w[2];
        tempQ442[49] +=  w[1]*(inQ4[49] + inQ4[99]);
        tempQ442[74] +=  inQ4[49]*w[2]  + inQ4[99]*w[0];


        TransposeInterpolationQ24Q4(tempQ442, tempQ24);
        TransposeInterpolationQ22Q24(tempQ24, outQ2); 
        
        TransposeInterpolationQ2Q4Bndy( tempQ442+25, outQ2+9 );

        TransposeInterpolationQ24Q4( tempQ442+50, tempQ24);        
        TransposeInterpolationQ22Q24( tempQ24, outQ2+18 );

    }

    inline void InterpolationQ2Q442Bndy( Real *inQ2, Real *outQ442)
    {
        Real tempQ24[15];
        // 1st layer
        InterpolationQ22Q24( inQ2, tempQ24 );
        InterpolationQ24Q4( tempQ24, outQ442);
        // 2nd layer
        InterpolationQ2Q4Bndy( inQ2+9, outQ442+25 );
        // 3rd layer
        InterpolationQ22Q24( inQ2+18, tempQ24 );
        InterpolationQ24Q4( tempQ24, outQ442+50);
    }

    inline void TransposeInterpolationQ2Q442Bndy( Real *inQ442, Real *outQ2)
    {
        Real tempQ24[15];

        TransposeInterpolationQ24Q4(inQ442, tempQ24);
        TransposeInterpolationQ22Q24(tempQ24, outQ2); 
        
        TransposeInterpolationQ2Q4Bndy( inQ442+25, outQ2+9 );

        TransposeInterpolationQ24Q4( inQ442+50, tempQ24);        
        TransposeInterpolationQ22Q24( tempQ24, outQ2+18 );

    }




    inline void InterpolationQ442Q4Bndy( Real *inQ442, Real *outQ4)
    {

        for (Index i = 0; i < 25; i++)
        {
            outQ4[i] = inQ442[i];
            outQ4[i+50] = inQ442[i+25];
            outQ4[i+100] = inQ442[i+50];
        }
        
        outQ4[25] = inQ442[0]*w[0] + inQ442[25]*w[1] + inQ442[50]*w[2];
        outQ4[26] = inQ442[1]*w[0] + inQ442[26]*w[1] + inQ442[51]*w[2];
        outQ4[27] = inQ442[2]*w[0] + inQ442[27]*w[1] + inQ442[52]*w[2];
        outQ4[28] = inQ442[3]*w[0] + inQ442[28]*w[1] + inQ442[53]*w[2];
        outQ4[29] = inQ442[4]*w[0] + inQ442[29]*w[1] + inQ442[54]*w[2];
        outQ4[30] = inQ442[5]*w[0] + inQ442[30]*w[1] + inQ442[55]*w[2];
        outQ4[31] = 0.0;
        outQ4[32] = 0.0;
        outQ4[33] = 0.0;
        outQ4[34] = inQ442[9]*w[0] + inQ442[34]*w[1] + inQ442[59]*w[2];
        outQ4[35] = inQ442[10]*w[0] + inQ442[35]*w[1] + inQ442[60]*w[2];
        outQ4[36] = 0.0;
        outQ4[37] = 0.0;
        outQ4[38] = 0.0;
        outQ4[39] = inQ442[14]*w[0] + inQ442[39]*w[1] + inQ442[64]*w[2];
        outQ4[40] = inQ442[15]*w[0] + inQ442[40]*w[1] + inQ442[65]*w[2];
        outQ4[41] = 0.0;
        outQ4[42] = 0.0;
        outQ4[43] = 0.0;
        outQ4[44] = inQ442[19]*w[0] + inQ442[44]*w[1] + inQ442[69]*w[2];
        outQ4[45] = inQ442[20]*w[0] + inQ442[45]*w[1] + inQ442[70]*w[2];
        outQ4[46] = inQ442[21]*w[0] + inQ442[46]*w[1] + inQ442[71]*w[2];
        outQ4[47] = inQ442[22]*w[0] + inQ442[47]*w[1] + inQ442[72]*w[2];
        outQ4[48] = inQ442[23]*w[0] + inQ442[48]*w[1] + inQ442[73]*w[2];
        outQ4[49] = inQ442[24]*w[0] + inQ442[49]*w[1] + inQ442[74]*w[2];

     

        outQ4[75] = inQ442[0]*w[2] + inQ442[25]*w[1] + inQ442[50]*w[0];
        outQ4[76] = inQ442[1]*w[2] + inQ442[26]*w[1] + inQ442[51]*w[0];
        outQ4[77] = inQ442[2]*w[2] + inQ442[27]*w[1] + inQ442[52]*w[0];
        outQ4[78] = inQ442[3]*w[2] + inQ442[28]*w[1] + inQ442[53]*w[0];
        outQ4[79] = inQ442[4]*w[2] + inQ442[29]*w[1] + inQ442[54]*w[0];
        outQ4[80] = inQ442[5]*w[2] + inQ442[30]*w[1] + inQ442[55]*w[0];
        outQ4[81] = 0.0;
        outQ4[82] = 0.0;
        outQ4[83] = 0.0;
        outQ4[84] = inQ442[9]*w[2] + inQ442[34]*w[1] + inQ442[59]*w[0];
        outQ4[85] = inQ442[10]*w[2] + inQ442[35]*w[1] + inQ442[60]*w[0];
        outQ4[86] = 0.0;
        outQ4[87] = 0.0;
        outQ4[88] = 0.0;
        outQ4[89] = inQ442[14]*w[2] + inQ442[39]*w[1] + inQ442[64]*w[0];
        outQ4[90] = inQ442[15]*w[2] + inQ442[40]*w[1] + inQ442[65]*w[0];
        outQ4[91] = 0.0;
        outQ4[92] = 0.0;
        outQ4[93] = 0.0;
        outQ4[94] = inQ442[19]*w[2] + inQ442[44]*w[1] + inQ442[69]*w[0];
        outQ4[95] = inQ442[20]*w[2] + inQ442[45]*w[1] + inQ442[70]*w[0];
        outQ4[96] = inQ442[21]*w[2] + inQ442[46]*w[1] + inQ442[71]*w[0];
        outQ4[97] = inQ442[22]*w[2] + inQ442[47]*w[1] + inQ442[72]*w[0];
        outQ4[98] = inQ442[23]*w[2] + inQ442[48]*w[1] + inQ442[73]*w[0];
        outQ4[99] = inQ442[24]*w[2] + inQ442[49]*w[1] + inQ442[74]*w[0];

        
    }

    inline void TransposeInterpolationQ442Q4Bndy( Real *inQ4, Real *outQ442)
    {
        for (Index i = 0; i < 25; i++)
        {
            outQ442[i] = inQ4[i];
            outQ442[i+25] = inQ4[i+50];
            outQ442[i+50] = inQ4[i+100];
        }
        
      

        outQ442[0]  +=  inQ4[25]*w[0]  + inQ4[75]*w[2];
        outQ442[25] +=  w[1]*(inQ4[25] + inQ4[75]);
        outQ442[50] +=  inQ4[25]*w[2]  + inQ4[75]*w[0];

        outQ442[1]  +=  inQ4[26]*w[0]  + inQ4[76]*w[2];
        outQ442[26] +=  w[1]*(inQ4[26] + inQ4[76]);
        outQ442[51] +=  inQ4[26]*w[2]  + inQ4[76]*w[0];

        outQ442[2]  +=  inQ4[27]*w[0]  + inQ4[77]*w[2];
        outQ442[27] +=  w[1]*(inQ4[27] + inQ4[77]);
        outQ442[52] +=  inQ4[27]*w[2]  + inQ4[77]*w[0];

        outQ442[3]  +=  inQ4[28]*w[0]  + inQ4[78]*w[2];
        outQ442[28] +=  w[1]*(inQ4[28] + inQ4[78]);
        outQ442[53] +=  inQ4[28]*w[2]  + inQ4[78]*w[0];

        outQ442[4]  +=  inQ4[29]*w[0]  + inQ4[79]*w[2];
        outQ442[29] +=  w[1]*(inQ4[29] + inQ4[79]);
        outQ442[54] +=  inQ4[29]*w[2]  + inQ4[79]*w[0];

        outQ442[5]  +=  inQ4[30]*w[0]  + inQ4[80]*w[2];
        outQ442[30] +=  w[1]*(inQ4[30] + inQ4[80]);
        outQ442[55] +=  inQ4[30]*w[2]  + inQ4[80]*w[0];

        outQ442[31]  = 0.0;
        outQ442[32]  = 0.0;
        outQ442[33]  = 0.0;

        outQ442[9]  +=  inQ4[34]*w[0]  + inQ4[84]*w[2];
        outQ442[34] +=  w[1]*(inQ4[34] + inQ4[84]);
        outQ442[59] +=  inQ4[34]*w[2]  + inQ4[84]*w[0];
        outQ442[10] +=  inQ4[35]*w[0]  + inQ4[85]*w[2];
        outQ442[35] +=  w[1]*(inQ4[35] + inQ4[85]);
        outQ442[60] +=  inQ4[35]*w[2]  + inQ4[85]*w[0];

        outQ442[36]  = 0.0;
        outQ442[37]  = 0.0;
        outQ442[38]  = 0.0;

        outQ442[14] +=  inQ4[39]*w[0]  + inQ4[89]*w[2];
        outQ442[39] +=  w[1]*(inQ4[39] + inQ4[89]);
        outQ442[64] +=  inQ4[39]*w[2]  + inQ4[89]*w[0];
        outQ442[15] +=  inQ4[40]*w[0]  + inQ4[90]*w[2];
        outQ442[40] +=  w[1]*(inQ4[40] + inQ4[90]);
        outQ442[65] +=  inQ4[40]*w[2]  + inQ4[90]*w[0];

        outQ442[41]  = 0.0;
        outQ442[42]  = 0.0;
        outQ442[43]  = 0.0;

        outQ442[19] +=  inQ4[44]*w[0]  + inQ4[94]*w[2];
        outQ442[44] +=  w[1]*(inQ4[44] + inQ4[94]);
        outQ442[69] +=  inQ4[44]*w[2]  + inQ4[94]*w[0];
        outQ442[20] +=  inQ4[45]*w[0]  + inQ4[95]*w[2];
        outQ442[45] +=  w[1]*(inQ4[45] + inQ4[95]);
        outQ442[70] +=  inQ4[45]*w[2]  + inQ4[95]*w[0];

        outQ442[21] +=  inQ4[46]*w[0]  + inQ4[96]*w[2];
        outQ442[46] +=  w[1]*(inQ4[46] + inQ4[96]);
        outQ442[71] +=  inQ4[46]*w[2]  + inQ4[96]*w[0];

        outQ442[22] +=  inQ4[47]*w[0]  + inQ4[97]*w[2];
        outQ442[47] +=  w[1]*(inQ4[47] + inQ4[97]);
        outQ442[72] +=  inQ4[47]*w[2]  + inQ4[97]*w[0];

        outQ442[23] +=  inQ4[48]*w[0]  + inQ4[98]*w[2];
        outQ442[48] +=  w[1]*(inQ4[48] + inQ4[98]);
        outQ442[73] +=  inQ4[48]*w[2]  + inQ4[98]*w[0];

        outQ442[24] +=  inQ4[49]*w[0]  + inQ4[99]*w[2];
        outQ442[49] +=  w[1]*(inQ4[49] + inQ4[99]);
        outQ442[74] +=  inQ4[49]*w[2]  + inQ4[99]*w[0];


    }


    /*--------- Q122->Q142->Q144->Q4 ---------*/

    inline void InterpolationQ122Q142( Real *inQ122,  Real *outQ142 )
    {
        InterpolationQ12Q14(inQ122, outQ142);
        InterpolationQ12Q14(inQ122+6, outQ142+10);
        InterpolationQ12Q14(inQ122+12, outQ142+20);
    }


    inline void InterpolationQ142Q144( Real *inQ142,  Real *outQ144 )
    {
        Real b[10];

        outQ144[0] = inQ142[0];
        outQ144[1] = inQ142[1];
        outQ144[2] = inQ142[2];
        outQ144[3] = inQ142[3];
        outQ144[4] = inQ142[4];
        outQ144[5] = inQ142[5];
        outQ144[6] = inQ142[6];
        outQ144[7] = inQ142[7];
        outQ144[8] = inQ142[8];
        outQ144[9] = inQ142[9];

        outQ144[20] = inQ142[10];
        outQ144[21] = inQ142[11];
        outQ144[22] = inQ142[12];
        outQ144[23] = inQ142[13];
        outQ144[24] = inQ142[14];
        outQ144[25] = inQ142[15];
        outQ144[26] = inQ142[16];
        outQ144[27] = inQ142[17];
        outQ144[28] = inQ142[18];
        outQ144[29] = inQ142[19];
 
        outQ144[40] = inQ142[20];
        outQ144[41] = inQ142[21];
        outQ144[42] = inQ142[22];
        outQ144[43] = inQ142[23];
        outQ144[44] = inQ142[24];
        outQ144[45] = inQ142[25];
        outQ144[46] = inQ142[26];
        outQ144[47] = inQ142[27];
        outQ144[48] = inQ142[28];
        outQ144[49] = inQ142[29];

        b[0] = inQ142[10]*w[1];
        b[1] = inQ142[11]*w[1];
        b[2] = inQ142[12]*w[1];
        b[3] = inQ142[13]*w[1];
        b[4] = inQ142[14]*w[1];
        b[5] = inQ142[15]*w[1];
        b[6] = inQ142[16]*w[1];
        b[7] = inQ142[17]*w[1];
        b[8] = inQ142[18]*w[1];
        b[9] = inQ142[19]*w[1];

        outQ144[10] = inQ142[0]*w[0] + b[0] +  inQ142[20]*w[2];
        outQ144[11] = inQ142[1]*w[0] + b[1] +  inQ142[21]*w[2];
        outQ144[12] = inQ142[2]*w[0] + b[2] +  inQ142[22]*w[2];
        outQ144[13] = inQ142[3]*w[0] + b[3] +  inQ142[23]*w[2];
        outQ144[14] = inQ142[4]*w[0] + b[4] +  inQ142[24]*w[2];
        outQ144[15] = inQ142[5]*w[0] + b[5] +  inQ142[25]*w[2];
        outQ144[16] = inQ142[6]*w[0] + b[6] +  inQ142[26]*w[2];
        outQ144[17] = inQ142[7]*w[0] + b[7] +  inQ142[27]*w[2];
        outQ144[18] = inQ142[8]*w[0] + b[8] +  inQ142[28]*w[2];
        outQ144[19] = inQ142[9]*w[0] + b[9] +  inQ142[29]*w[2];
       
        outQ144[30] = inQ142[0]*w[2] + b[0] +  inQ142[20]*w[0];
        outQ144[31] = inQ142[1]*w[2] + b[1] +  inQ142[21]*w[0];
        outQ144[32] = inQ142[2]*w[2] + b[2] +  inQ142[22]*w[0];
        outQ144[33] = inQ142[3]*w[2] + b[3] +  inQ142[23]*w[0];
        outQ144[34] = inQ142[4]*w[2] + b[4] +  inQ142[24]*w[0];
        outQ144[35] = inQ142[5]*w[2] + b[5] +  inQ142[25]*w[0];
        outQ144[36] = inQ142[6]*w[2] + b[6] +  inQ142[26]*w[0];
        outQ144[37] = inQ142[7]*w[2] + b[7] +  inQ142[27]*w[0];
        outQ144[38] = inQ142[8]*w[2] + b[8] +  inQ142[28]*w[0];
        outQ144[39] = inQ142[9]*w[2] + b[9] +  inQ142[29]*w[0];        

    }

    inline void InterpolationQ144Q4( Real *inQ144,  Real *outQ4 )
    {
        InterpolationQ14Q4(inQ144, outQ4);
        InterpolationQ14Q4(inQ144+10, outQ4+25);
        InterpolationQ14Q4(inQ144+20, outQ4+50);
        InterpolationQ14Q4(inQ144+30, outQ4+75);
        InterpolationQ14Q4(inQ144+40, outQ4+100);
    }

    
    /*---------  Q212->Q412->Q414->Q4 ---------*/
    
    inline void InterpolationQ212Q412( Real *inQ212,  Real *outQ412 )
    {
        InterpolationQ21Q41(inQ212, outQ412);
        InterpolationQ21Q41(inQ212+6, outQ412+10);
        InterpolationQ21Q41(inQ212+12, outQ412+20);
    }

    inline void InterpolationQ412Q414( Real *inQ412,  Real *outQ414 )
    {
        InterpolationQ142Q144( inQ412,  outQ414 );
    }

    inline void InterpolationQ414Q4( Real *inQ414,  Real *outQ4 )
    {
        InterpolationQ41Q4(inQ414, outQ4);
        InterpolationQ41Q4(inQ414+10, outQ4+25);
        InterpolationQ41Q4(inQ414+20, outQ4+50);
        InterpolationQ41Q4(inQ414+30, outQ4+75);
        InterpolationQ41Q4(inQ414+40, outQ4+100);
    }


    /*---------  Q221->Q241->Q441->Q4 ---------*/
    inline void InterpolationQ221Q241( Real *inQ221,  Real *outQ241 )
    {
        InterpolationQ22Q24(inQ221, outQ241);
        InterpolationQ22Q24(inQ221+9, outQ241+15);
    }

    inline void InterpolationQ241Q441( Real *inQ241,  Real *outQ441 )
    {
        InterpolationQ24Q4( inQ241, outQ441 );
        InterpolationQ24Q4( inQ241+15, outQ441+25 );

    }

    inline void InterpolationQ441Q4( Real *inQ441, Real *outQ4 )
    { 
        for (Index i = 0; i < 25; i++)
        {
            outQ4[i] = inQ441[i];
            outQ4[i+100] = inQ441[i+25];
        }
        

        outQ4[25] = inQ441[0]*z[0] + inQ441[25]*z[1] ;
        outQ4[50] = 0.5 * (inQ441[0] + inQ441[25]);
        outQ4[75] = inQ441[0]*z[1] + inQ441[25]*z[0];

        outQ4[26] = inQ441[1]*z[0] + inQ441[26]*z[1] ;
        outQ4[51] = 0.5*(inQ441[1] + inQ441[26]);
        outQ4[76] = inQ441[1]*z[1] + inQ441[26]*z[0];

        outQ4[27] = inQ441[2]*z[0] + inQ441[27]*z[1] ;
        outQ4[52] = 0.5*(inQ441[2] + inQ441[27]);
        outQ4[77] = inQ441[2]*z[1] + inQ441[27]*z[0];

        outQ4[28] = inQ441[3]*z[0] + inQ441[28]*z[1] ;
        outQ4[53] = 0.5*(inQ441[3] + inQ441[28]);
        outQ4[78] = inQ441[3]*z[1] + inQ441[28]*z[0];

        outQ4[29] = inQ441[4]*z[0] + inQ441[29]*z[1] ;
        outQ4[54] = 0.5*(inQ441[4] + inQ441[29]);
        outQ4[79] = inQ441[4]*z[1] + inQ441[29]*z[0];

        outQ4[30] = inQ441[5]*z[0] + inQ441[30]*z[1] ;
        outQ4[55] = 0.5*(inQ441[5] + inQ441[30]);
        outQ4[80] = inQ441[5]*z[1] + inQ441[30]*z[0];

        outQ4[31] = inQ441[6]*z[0] + inQ441[31]*z[1] ;
        outQ4[56] = 0.5*(inQ441[6] + inQ441[31]);
        outQ4[81] = inQ441[6]*z[1] + inQ441[31]*z[0];

        outQ4[32] = inQ441[7]*z[0] + inQ441[32]*z[1] ;
        outQ4[57] = 0.5*(inQ441[7] + inQ441[32]);
        outQ4[82] = inQ441[7]*z[1] + inQ441[32]*z[0];

        outQ4[33] = inQ441[8]*z[0] + inQ441[33]*z[1] ;
        outQ4[58] = 0.5*(inQ441[8] + inQ441[33]);
        outQ4[83] = inQ441[8]*z[1] + inQ441[33]*z[0];

        outQ4[34] = inQ441[9]*z[0] + inQ441[34]*z[1] ;
        outQ4[59] = 0.5*(inQ441[9] + inQ441[34]);
        outQ4[84] = inQ441[9]*z[1] + inQ441[34]*z[0];

        outQ4[35] = inQ441[10]*z[0] + inQ441[35]*z[1] ;
        outQ4[60] = 0.5*(inQ441[10] + inQ441[35]);
        outQ4[85] = inQ441[10]*z[1] + inQ441[35]*z[0];

        outQ4[36] = inQ441[11]*z[0] + inQ441[36]*z[1] ;
        outQ4[61] = 0.5*(inQ441[11] + inQ441[36]);
        outQ4[86] = inQ441[11]*z[1] + inQ441[36]*z[0];

        outQ4[37] = inQ441[12]*z[0] + inQ441[37]*z[1] ;
        outQ4[62] = 0.5*(inQ441[12] + inQ441[37]);
        outQ4[87] = inQ441[12]*z[1] + inQ441[37]*z[0];

        outQ4[38] = inQ441[13]*z[0] + inQ441[38]*z[1] ;
        outQ4[63] = 0.5*(inQ441[13] + inQ441[38]);
        outQ4[88] = inQ441[13]*z[1] + inQ441[38]*z[0];

        outQ4[39] = inQ441[14]*z[0] + inQ441[39]*z[1] ;
        outQ4[64] = 0.5*(inQ441[14] + inQ441[39]);
        outQ4[89] = inQ441[14]*z[1] + inQ441[39]*z[0];

        outQ4[40] = inQ441[15]*z[0] + inQ441[40]*z[1] ;
        outQ4[65] = 0.5*(inQ441[15] + inQ441[40]);
        outQ4[90] = inQ441[15]*z[1] + inQ441[40]*z[0];

        outQ4[41] = inQ441[16]*z[0] + inQ441[41]*z[1] ;
        outQ4[66] = 0.5*(inQ441[16] + inQ441[41]);
        outQ4[91] = inQ441[16]*z[1] + inQ441[41]*z[0];

        outQ4[42] = inQ441[17]*z[0] + inQ441[42]*z[1] ;
        outQ4[67] = 0.5*(inQ441[17] + inQ441[42]);
        outQ4[92] = inQ441[17]*z[1] + inQ441[42]*z[0];

        outQ4[43] = inQ441[18]*z[0] + inQ441[43]*z[1] ;
        outQ4[68] = 0.5*(inQ441[18] + inQ441[43]);
        outQ4[93] = inQ441[18]*z[1] + inQ441[43]*z[0];

        outQ4[44] = inQ441[19]*z[0] + inQ441[44]*z[1] ;
        outQ4[69] = 0.5*(inQ441[19] + inQ441[44]);
        outQ4[94] = inQ441[19]*z[1] + inQ441[44]*z[0];

        outQ4[45] = inQ441[20]*z[0] + inQ441[45]*z[1] ;
        outQ4[70] = 0.5*(inQ441[20] + inQ441[45]);
        outQ4[95] = inQ441[20]*z[1] + inQ441[45]*z[0];

        outQ4[46] = inQ441[21]*z[0] + inQ441[46]*z[1] ;
        outQ4[71] = 0.5*(inQ441[21] + inQ441[46]);
        outQ4[96] = inQ441[21]*z[1] + inQ441[46]*z[0];

        outQ4[47] = inQ441[22]*z[0] + inQ441[47]*z[1] ;
        outQ4[72] = 0.5*(inQ441[22] + inQ441[47]);
        outQ4[97] = inQ441[22]*z[1] + inQ441[47]*z[0];

        outQ4[48] = inQ441[23]*z[0] + inQ441[48]*z[1] ;
        outQ4[73] = 0.5*(inQ441[23] + inQ441[48]);
        outQ4[98] = inQ441[23]*z[1] + inQ441[48]*z[0];

        outQ4[49] = inQ441[24]*z[0] + inQ441[49]*z[1] ;
        outQ4[74] = 0.5*(inQ441[24] + inQ441[49]);
        outQ4[99] = inQ441[24]*z[1] + inQ441[49]*z[0];






    }

    /*--------- TRANSPOSE Q122->Q142->Q144->Q4 ---------*/

    inline void TransposeInterpolationQ122Q142( Real *inQ142,  Real *outQ122 )
    {
        TransposeInterpolationQ12Q14(inQ142, outQ122 );
        TransposeInterpolationQ12Q14(inQ142+10, outQ122+6 );
        TransposeInterpolationQ12Q14(inQ142+20, outQ122+12 );

        
    }


    inline void TransposeInterpolationQ142Q144( Real *inQ144,  Real *outQ142 )
    {
        Real b[10];

        outQ142[0] = inQ144[0];
        outQ142[1] = inQ144[1];
        outQ142[2] = inQ144[2];
        outQ142[3] = inQ144[3];
        outQ142[4] = inQ144[4];
        outQ142[5] = inQ144[5];
        outQ142[6] = inQ144[6];
        outQ142[7] = inQ144[7];
        outQ142[8] = inQ144[8];
        outQ142[9] = inQ144[9];

        outQ142[10] = inQ144[20];
        outQ142[11] = inQ144[21];
        outQ142[12] = inQ144[22];
        outQ142[13] = inQ144[23];
        outQ142[14] = inQ144[24];
        outQ142[15] = inQ144[25];
        outQ142[16] = inQ144[26];
        outQ142[17] = inQ144[27];
        outQ142[18] = inQ144[28];
        outQ142[19] = inQ144[29];
 
        outQ142[20] = inQ144[40];
        outQ142[21] = inQ144[41];
        outQ142[22] = inQ144[42];
        outQ142[23] = inQ144[43];
        outQ142[24] = inQ144[44];
        outQ142[25] = inQ144[45];
        outQ142[26] = inQ144[46];
        outQ142[27] = inQ144[47];
        outQ142[28] = inQ144[48];
        outQ142[29] = inQ144[49];

        outQ142[0]  +=  inQ144[10]*w[0]  + inQ144[30]*w[2];
        outQ142[10] +=  w[1]*(inQ144[10] + inQ144[30]);
        outQ142[20] +=  inQ144[10]*w[2]  + inQ144[30]*w[0];

        outQ142[1]  +=  inQ144[11]*w[0]  + inQ144[31]*w[2];
        outQ142[11] +=  w[1]*(inQ144[11] + inQ144[31]);
        outQ142[21] +=  inQ144[11]*w[2]  + inQ144[31]*w[0];

        outQ142[2]  +=  inQ144[12]*w[0]  + inQ144[32]*w[2];
        outQ142[12] +=  w[1]*(inQ144[12] + inQ144[32]);
        outQ142[22] +=  inQ144[12]*w[2]  + inQ144[32]*w[0];

        outQ142[3]  +=  inQ144[13]*w[0]  + inQ144[33]*w[2];
        outQ142[13] +=  w[1]*(inQ144[13] + inQ144[33]);
        outQ142[23] +=  inQ144[13]*w[2]  + inQ144[33]*w[0];

        outQ142[4]  +=  inQ144[14]*w[0]  + inQ144[34]*w[2];
        outQ142[14] +=  w[1]*(inQ144[14] + inQ144[34]);
        outQ142[24] +=  inQ144[14]*w[2]  + inQ144[34]*w[0];
        
        outQ142[5]  +=  inQ144[15]*w[0]  + inQ144[35]*w[2];
        outQ142[15] +=  w[1]*(inQ144[15] + inQ144[35]);
        outQ142[25] +=  inQ144[15]*w[2]  + inQ144[35]*w[0];

        outQ142[6]  +=  inQ144[16]*w[0]  + inQ144[36]*w[2];
        outQ142[16] +=  w[1]*(inQ144[16] + inQ144[36]);
        outQ142[26] +=  inQ144[16]*w[2]  + inQ144[36]*w[0];

        outQ142[7]  +=  inQ144[17]*w[0]  + inQ144[37]*w[2];
        outQ142[17] +=  w[1]*(inQ144[17] + inQ144[37]);
        outQ142[27] +=  inQ144[17]*w[2]  + inQ144[37]*w[0];
        
        outQ142[8]  +=  inQ144[18]*w[0]  + inQ144[38]*w[2];
        outQ142[18] +=  w[1]*(inQ144[18] + inQ144[38]);
        outQ142[28] +=  inQ144[18]*w[2]  + inQ144[38]*w[0];
        
        outQ142[9]  +=  inQ144[19]*w[0]  + inQ144[39]*w[2];
        outQ142[19] +=  w[1]*(inQ144[19] + inQ144[39]);
        outQ142[29] +=  inQ144[19]*w[2]  + inQ144[39]*w[0];
    }

    inline void TransposeInterpolationQ144Q4( Real *inQ4,  Real *outQ144 )
    {
        TransposeInterpolationQ14Q4(inQ4, outQ144);
        TransposeInterpolationQ14Q4(inQ4+25, outQ144+10);
        TransposeInterpolationQ14Q4(inQ4+50, outQ144+20);
        TransposeInterpolationQ14Q4(inQ4+75, outQ144+30);
        TransposeInterpolationQ14Q4(inQ4+100, outQ144+40);
    }

     /*--------- TRANSPOSE Q212->Q412->Q414->Q4 ---------*/
    inline void TransposeInterpolationQ212Q412( Real *inQ412,  Real *outQ212 )
    {
        TransposeInterpolationQ21Q41(inQ412, outQ212);
        TransposeInterpolationQ21Q41(inQ412+10, outQ212+6);
        TransposeInterpolationQ21Q41(inQ412+20, outQ212+12);
    }

    inline void TransposeInterpolationQ412Q414( Real *inQ414,  Real *outQ412 )
    {
        TransposeInterpolationQ142Q144( inQ414,  outQ412 );
    }

    inline void TransposeInterpolationQ414Q4( Real *inQ4,  Real *outQ414 )
    {
        TransposeInterpolationQ41Q4(inQ4, outQ414);
        TransposeInterpolationQ41Q4(inQ4+25, outQ414+10);
        TransposeInterpolationQ41Q4(inQ4+50, outQ414+20);
        TransposeInterpolationQ41Q4(inQ4+75, outQ414+30);
        TransposeInterpolationQ41Q4(inQ4+100, outQ414+40);
    }


     /*---------  TRANSPOSE Q221->Q241->Q441->Q4 ---------*/

    inline void TransposeInterpolationQ221Q241( Real *inQ241,  Real *outQ221 )
    {
        TransposeInterpolationQ22Q24(inQ241, outQ221);
        TransposeInterpolationQ22Q24(inQ241+15, outQ221+9);
    }

    inline void TransposeInterpolationQ241Q441( Real *inQ441,  Real *outQ241 )
    {
        TransposeInterpolationQ24Q4( inQ441, outQ241 );
        TransposeInterpolationQ24Q4( inQ441+25, outQ241+15 );
    }

    inline void TransposeInterpolationQ441Q4( Real *inQ4, Real *outQ441 )
     {
     
        for (Index i = 0; i < 25; i++)
        {
            outQ441[i] = inQ4[i];
            outQ441[i+25] = inQ4[i+100];
        }
        
        outQ441[0]  += inQ4[25]*z[0] + inQ4[50]*0.5 + inQ4[75]*z[1]; 
        outQ441[25] += inQ4[25]*z[1] + inQ4[50]*0.5 + inQ4[75]*z[0]; 
        outQ441[1]  += inQ4[26]*z[0] + inQ4[51]*0.5 + inQ4[76]*z[1]; 
        outQ441[26] += inQ4[26]*z[1] + inQ4[51]*0.5 + inQ4[76]*z[0]; 
        outQ441[2]  += inQ4[27]*z[0] + inQ4[52]*0.5 + inQ4[77]*z[1]; 
        outQ441[27] += inQ4[27]*z[1] + inQ4[52]*0.5 + inQ4[77]*z[0]; 
        outQ441[3]  += inQ4[28]*z[0] + inQ4[53]*0.5 + inQ4[78]*z[1]; 
        outQ441[28] += inQ4[28]*z[1] + inQ4[53]*0.5 + inQ4[78]*z[0]; 
        outQ441[4]  += inQ4[29]*z[0] + inQ4[54]*0.5 + inQ4[79]*z[1]; 
        outQ441[29] += inQ4[29]*z[1] + inQ4[54]*0.5 + inQ4[79]*z[0]; 
        outQ441[5]  += inQ4[30]*z[0] + inQ4[55]*0.5 + inQ4[80]*z[1]; 
        outQ441[30] += inQ4[30]*z[1] + inQ4[55]*0.5 + inQ4[80]*z[0]; 
        outQ441[6]  += inQ4[31]*z[0] + inQ4[56]*0.5 + inQ4[81]*z[1]; 
        outQ441[31] += inQ4[31]*z[1] + inQ4[56]*0.5 + inQ4[81]*z[0]; 
        outQ441[7]  += inQ4[32]*z[0] + inQ4[57]*0.5 + inQ4[82]*z[1]; 
        outQ441[32] += inQ4[32]*z[1] + inQ4[57]*0.5 + inQ4[82]*z[0]; 
        outQ441[8]  += inQ4[33]*z[0] + inQ4[58]*0.5 + inQ4[83]*z[1]; 
        outQ441[33] += inQ4[33]*z[1] + inQ4[58]*0.5 + inQ4[83]*z[0]; 
        outQ441[9]  += inQ4[34]*z[0] + inQ4[59]*0.5 + inQ4[84]*z[1]; 
        outQ441[34] += inQ4[34]*z[1] + inQ4[59]*0.5 + inQ4[84]*z[0]; 
        outQ441[10] += inQ4[35]*z[0] + inQ4[60]*0.5 + inQ4[85]*z[1]; 
        outQ441[35] += inQ4[35]*z[1] + inQ4[60]*0.5 + inQ4[85]*z[0]; 
        outQ441[11] += inQ4[36]*z[0] + inQ4[61]*0.5 + inQ4[86]*z[1]; 
        outQ441[36] += inQ4[36]*z[1] + inQ4[61]*0.5 + inQ4[86]*z[0]; 
    
        outQ441[12] += inQ4[37]*z[0] + inQ4[62]*0.5 + inQ4[87]*z[1]; 
        outQ441[37] += inQ4[37]*z[1] + inQ4[62]*0.5 + inQ4[87]*z[0]; 
        outQ441[13] += inQ4[38]*z[0] + inQ4[63]*0.5 + inQ4[88]*z[1]; 
        outQ441[38] += inQ4[38]*z[1] + inQ4[63]*0.5 + inQ4[88]*z[0]; 
        outQ441[14] += inQ4[39]*z[0] + inQ4[64]*0.5 + inQ4[89]*z[1]; 
        outQ441[39] += inQ4[39]*z[1] + inQ4[64]*0.5 + inQ4[89]*z[0]; 
        outQ441[15] += inQ4[40]*z[0] + inQ4[65]*0.5 + inQ4[90]*z[1]; 
        outQ441[40] += inQ4[40]*z[1] + inQ4[65]*0.5 + inQ4[90]*z[0]; 
        outQ441[16] += inQ4[41]*z[0] + inQ4[66]*0.5 + inQ4[91]*z[1]; 
        outQ441[41] += inQ4[41]*z[1] + inQ4[66]*0.5 + inQ4[91]*z[0]; 
        outQ441[17] += inQ4[42]*z[0] + inQ4[67]*0.5 + inQ4[92]*z[1]; 
        outQ441[42] += inQ4[42]*z[1] + inQ4[67]*0.5 + inQ4[92]*z[0]; 
        outQ441[18] += inQ4[43]*z[0] + inQ4[68]*0.5 + inQ4[93]*z[1]; 
        outQ441[43] += inQ4[43]*z[1] + inQ4[68]*0.5 + inQ4[93]*z[0]; 
        outQ441[19] += inQ4[44]*z[0] + inQ4[69]*0.5 + inQ4[94]*z[1]; 
        outQ441[44] += inQ4[44]*z[1] + inQ4[69]*0.5 + inQ4[94]*z[0]; 
        outQ441[20] += inQ4[45]*z[0] + inQ4[70]*0.5 + inQ4[95]*z[1]; 
        outQ441[45] += inQ4[45]*z[1] + inQ4[70]*0.5 + inQ4[95]*z[0]; 
        outQ441[21] += inQ4[46]*z[0] + inQ4[71]*0.5 + inQ4[96]*z[1]; 
        outQ441[46] += inQ4[46]*z[1] + inQ4[71]*0.5 + inQ4[96]*z[0]; 
        outQ441[22] += inQ4[47]*z[0] + inQ4[72]*0.5 + inQ4[97]*z[1]; 
        outQ441[47] += inQ4[47]*z[1] + inQ4[72]*0.5 + inQ4[97]*z[0]; 
        outQ441[23] += inQ4[48]*z[0] + inQ4[73]*0.5 + inQ4[98]*z[1]; 
        outQ441[48] += inQ4[48]*z[1] + inQ4[73]*0.5 + inQ4[98]*z[0]; 
        outQ441[24] += inQ4[49]*z[0] + inQ4[74]*0.5 + inQ4[99]*z[1]; 
        outQ441[49] += inQ4[49]*z[1] + inQ4[74]*0.5 + inQ4[99]*z[0]; 
    }

    

   
} // OndoMathX


 
