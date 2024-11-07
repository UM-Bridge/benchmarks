#include <iostream>
#include <vector>
#include <ctime>
#include <algorithm>

#include "OndoMathX.h"

using namespace OndoMathX;


int main(int argc, char *argv[])
{
 
 

    clog<< "------------------ TEST IN 2D ------------------"<<endl;


    GaussLobattoElement GLE_Q12(2,3,0);
    GaussLobattoElement GLE_Q21(3,2,0);
    GaussLobattoElement GLE_Q22(3,3,0);
    GaussLobattoElement GLE_Q44(5,5,0);

    GaussLobattoElement GLE_Q221(3,3,2);
    GaussLobattoElement GLE_Q212(3,2,3);
    GaussLobattoElement GLE_Q122(2,3,3);
    GaussLobattoElement GLE_Q222(3,3,3);
    GaussLobattoElement GLE_Q444(5,5,5);


         

    GaussLobattoInterpolator GLI_Q12_Q44(GLE_Q12,GLE_Q44);
    GaussLobattoInterpolator GLI_Q21_Q44(GLE_Q21,GLE_Q44);
    GaussLobattoInterpolator GLI_Q22_Q44(GLE_Q22,GLE_Q44);
    GaussLobattoInterpolator GLI_Q221_Q444(GLE_Q221,GLE_Q444);
    GaussLobattoInterpolator GLI_Q212_Q444(GLE_Q212,GLE_Q444);
    GaussLobattoInterpolator GLI_Q122_Q444(GLE_Q122,GLE_Q444);
    GaussLobattoInterpolator GLI_Q222_Q444(GLE_Q222,GLE_Q444);


   //Testing  Q12->Q4 Interpolation
    {
        std::clog << "---------------------------------------------"<<std::endl;


        std::cout << "Testing  Q12->Q4 Interpolation " << std::endl;
        // Testing correctness
        vector<Real> Vec(6);
        std::srand(unsigned(std::time(nullptr)));
        std::generate(Vec.begin(), Vec.end(), std::rand);

        vector<Real> InterpolationVec_Generic(25);
        GLI_Q12_Q44.Interpolation(Vec,InterpolationVec_Generic);

        Real InterpolationVec_Q14[10];
        Real InterpolationVec_HardCoded[25];
        InterpolationQ12Q14(Vec.data(),InterpolationVec_Q14);
        InterpolationQ14Q4(InterpolationVec_Q14, InterpolationVec_HardCoded );

        Real error_square = 0.0;
        Real norm_square = 0.0;
        Real diff = 0.0;

        for (Index i=0;i<25;++i)
        {
            diff = InterpolationVec_Generic[i]-InterpolationVec_HardCoded[i];
            error_square += diff*diff;
            norm_square +=  InterpolationVec_Generic[i]*InterpolationVec_Generic[i];
        }

        cout << "Error in the interpolation " << sqrt(error_square/norm_square) << std::endl;

        time_t start, end;

        time(&start);
        for (Index i=0;i<200000000;++i)
        {
            GLI_Q12_Q44.Interpolation(Vec,InterpolationVec_Generic);
            
        }
        time(&end);
    
        std::cout << "Elapsed time Generic Interpolation :" << difftime(end, start) << " seconds " << std::endl;
        

        time(&start);
        for (Index i=0;i<200000000;++i)
        {
            InterpolationQ12Q14(Vec.data(),InterpolationVec_Q14);
            InterpolationQ14Q4(InterpolationVec_Q14, InterpolationVec_HardCoded );

        }
        time(&end);
    
        std::cout << "Elapsed time Hard Coded Interpolation :" << difftime(end, start) << " seconds " << std::endl;

    }

 



   //Testing  Q12->Q4 TransposeInterpolation
    {
        std::clog << "---------------------------------------------"<<std::endl;

        
        std::cout << "Testing  Q12->Q4 TransposeInterpolation " << std::endl;
        // Testing correctness
        vector<Real> Vec(25);
        std::srand(unsigned(std::time(nullptr)));
        std::generate(Vec.begin(), Vec.end(), std::rand);

        vector<Real> InterpolationVec_Generic(6);
        GLI_Q12_Q44.TransposeInterpolation(Vec,InterpolationVec_Generic);

        Real InterpolationVec_HardCoded[6];
        Real InterpolationVec_Q14[10];
        TransposeInterpolationQ14Q4(Vec.data(),InterpolationVec_Q14);
        TransposeInterpolationQ12Q14(InterpolationVec_Q14, InterpolationVec_HardCoded);

        Real error_square = 0.0;
        Real norm_square = 0.0;
        Real diff = 0.0;

        for (Index i=0;i<6;++i)
        {
            diff = InterpolationVec_Generic[i]-InterpolationVec_HardCoded[i];
            error_square += diff*diff;
            norm_square +=  InterpolationVec_Generic[i]*InterpolationVec_Generic[i];
        }

        cout << "Error in the interpolation " << sqrt(error_square/norm_square) << std::endl;


        time_t start, end;

        time(&start);
        for (Index i=0;i<200000000;++i)
        {
            GLI_Q12_Q44.TransposeInterpolation(Vec,InterpolationVec_Generic);
            
        }
        time(&end);
    
        std::cout << "Elapsed time Generic Interpolation :" << difftime(end, start) << " seconds " << std::endl;
        

        time(&start);
        for (Index i=0;i<200000000;++i)
        {
            TransposeInterpolationQ14Q4(Vec.data(),InterpolationVec_Q14);
            TransposeInterpolationQ12Q14(InterpolationVec_Q14, InterpolationVec_HardCoded);
        }
        time(&end);
    
        std::cout << "Elapsed time Hard Coded Interpolation :" << difftime(end, start) << " seconds " << std::endl;


    }


   


    //Testing  Q21->Q4 Interpolation
    {
        std::clog << "---------------------------------------------"<<std::endl;

  
        std::cout << "Testing  Q21->Q4 Interpolation " << std::endl;
        // Testing correctness
        vector<Real> Vec(6);
        std::srand(unsigned(std::time(nullptr)));
        std::generate(Vec.begin(), Vec.end(), std::rand);

        vector<Real> InterpolationVec_Generic(25);
        GLI_Q21_Q44.Interpolation(Vec,InterpolationVec_Generic);

        Real InterpolationVec_HardCoded[25];
        Real InterpolationVec_Q41[10];
        InterpolationQ21Q41(Vec.data(),InterpolationVec_Q41);
        InterpolationQ41Q4(InterpolationVec_Q41, InterpolationVec_HardCoded);

        Real error_square = 0.0;
        Real norm_square = 0.0;
        Real diff = 0.0;
 
        for (Index i=0;i<25;++i)
        {
            diff = InterpolationVec_Generic[i]-InterpolationVec_HardCoded[i];
            error_square += diff*diff;
            norm_square +=  InterpolationVec_Generic[i]*InterpolationVec_Generic[i];
        }

        cout << "Error in the interpolation " << sqrt(error_square/norm_square) << std::endl;


        time_t start, end;

        time(&start);
        for (Index i=0;i<200000000;++i)
        {
            GLI_Q21_Q44.Interpolation(Vec,InterpolationVec_Generic);
            
        }
        time(&end);
    
        std::cout << "Elapsed time Generic Interpolation :" << difftime(end, start) << " seconds " << std::endl;
        

        time(&start);
        for (Index i=0;i<200000000;++i)
        {
            InterpolationQ21Q41(Vec.data(),InterpolationVec_Q41);
            InterpolationQ41Q4(InterpolationVec_Q41, InterpolationVec_HardCoded);
        }
        time(&end);
    
        std::cout << "Elapsed time Hard Coded Interpolation :" << difftime(end, start) << " seconds " << std::endl;
    }

    
    


    {
        std::clog << "---------------------------------------------"<<std::endl;

        //Testing  Q21->Q4 TransposeInterpolation (Q4 ->Q41 -> Q21)
      
        std::cout << "Testing  Q21->Q4 TransposeInterpolation " << std::endl;
        // Testing correctness
        vector<Real> Vec(25);
        std::srand(unsigned(std::time(nullptr)));
        std::generate(Vec.begin(), Vec.end(), std::rand);

        vector<Real> InterpolationVec_Generic(6);
        GLI_Q21_Q44.TransposeInterpolation(Vec,InterpolationVec_Generic);

        Real InterpolationVec_VecQ41[10];
        Real InterpolationVec_HardCoded[6];
        
        TransposeInterpolationQ41Q4(Vec.data() , InterpolationVec_VecQ41 );
        TransposeInterpolationQ21Q41(InterpolationVec_VecQ41 ,InterpolationVec_HardCoded);

        Real error_square = 0.0;
        Real norm_square = 0.0;
        Real diff = 0.0;
 
        for (Index i=0;i<6;++i)
        {
            diff = InterpolationVec_Generic[i]-InterpolationVec_HardCoded[i];
            error_square += diff*diff;
            norm_square +=  InterpolationVec_Generic[i]*InterpolationVec_Generic[i];
        }

        cout << "Error in the interpolation " << sqrt(error_square/norm_square) << std::endl;


        time_t start, end;

        time(&start);
        for (Index i=0;i<200000000;++i)
        {
            GLI_Q21_Q44.TransposeInterpolation(Vec,InterpolationVec_Generic);
            
        }
        time(&end);
    
        std::cout << "Elapsed time Generic Interpolation :" << difftime(end, start) << " seconds " << std::endl;
        

        time(&start);
        for (Index i=0;i<200000000;++i)
        {
            TransposeInterpolationQ41Q4(Vec.data() , InterpolationVec_VecQ41 );
            TransposeInterpolationQ21Q41(InterpolationVec_VecQ41 ,InterpolationVec_HardCoded);
        }
        time(&end);
    
        std::cout << "Elapsed time Hard Coded Interpolation :" << difftime(end, start) << " seconds " << std::endl;
    }




    {
        std::clog << "---------------------------------------------"<<std::endl;

       
        std::cout << "Testing  Q2->Q4 Boundary Interpolation " << std::endl;
        // Testing correctness
        vector<Real> Vec(9);
        std::srand(unsigned(std::time(nullptr)));
        std::generate(Vec.begin(), Vec.end(), std::rand);

        vector<Real> InterpolationVec_Generic(25);
        GLI_Q22_Q44.Interpolation(Vec,InterpolationVec_Generic);

        Real InterpolationVec_HardCoded[25];
        InterpolationQ2Q4Bndy(Vec.data(),InterpolationVec_HardCoded);
     

        Real error_square = 0.0;
        Real norm_square = 0.0;
        Real diff = 0.0;
 
        for (Index i=0;i<25;++i)
        {
            if(InterpolationVec_HardCoded[i]>0)
            {
                diff = InterpolationVec_Generic[i]-InterpolationVec_HardCoded[i];
                error_square += diff*diff;
                norm_square +=  InterpolationVec_Generic[i]*InterpolationVec_Generic[i];
            }
            
        }
        cout << "Error in the interpolation " << sqrt(error_square/norm_square) << std::endl;


        time_t start, end;
        time(&start);
        for (Index i=0;i<200000000;++i)
        {
            GLI_Q22_Q44.Interpolation(Vec,InterpolationVec_Generic);
            
        }
        time(&end);
        std::cout << "Elapsed time Generic Interpolation :" << difftime(end, start) << " seconds " << std::endl;
        

        time(&start);
        for (Index i=0;i<200000000;++i)
        {
            InterpolationQ2Q4Bndy(Vec.data(),InterpolationVec_HardCoded);
        }
        time(&end);
        std::cout << "Elapsed time Hard Coded Interpolation :" << difftime(end, start) << " seconds " << std::endl;
    }






    {
        std::clog << "---------------------------------------------"<<std::endl;

        
        std::cout << "Testing  Q2->Q4 Boundary  Transpose Interpolation " << std::endl;
        // Testing correctness
        vector<Real> Vec(25);
        std::srand(unsigned(std::time(nullptr)));
        std::generate(Vec.begin(), Vec.end(), std::rand);

        Vec[6] = 0.0;
        Vec[7] = 0.0;
        Vec[8] = 0.0;

        Vec[11] = 0.0;
        Vec[12] = 0.0;
        Vec[13] = 0.0;

        Vec[16] = 0.0;
        Vec[17] = 0.0;
        Vec[18] = 0.0;


        vector<Real> InterpolationVec_Generic(9);
        GLI_Q22_Q44.TransposeInterpolation(Vec,InterpolationVec_Generic);

        Real InterpolationVec_HardCoded[9];
        TransposeInterpolationQ2Q4Bndy(Vec.data(),InterpolationVec_HardCoded);
     

        Real error_square = 0.0;
        Real norm_square = 0.0;
        Real diff = 0.0;
 
        for (Index i=0;i<9;++i)
        {
            if(InterpolationVec_HardCoded[i]>0)
            {
                diff = InterpolationVec_Generic[i]-InterpolationVec_HardCoded[i];
                error_square += diff*diff;
                norm_square +=  InterpolationVec_Generic[i]*InterpolationVec_Generic[i];
            }
            
        }
        cout << "Error in the interpolation " << sqrt(error_square/norm_square) << std::endl;


        time_t start, end;
        time(&start);
        for (Index i=0;i<200000000;++i)
        {
            GLI_Q22_Q44.TransposeInterpolation(Vec,InterpolationVec_Generic);
            
        }
        time(&end);
        std::cout << "Elapsed time Generic Interpolation :" << difftime(end, start) << " seconds " << std::endl;
        

        time(&start);
        for (Index i=0;i<200000000;++i)
        {
            TransposeInterpolationQ2Q4Bndy(Vec.data(),InterpolationVec_HardCoded);
        }
        time(&end);
        std::cout << "Elapsed time Hard Coded Interpolation :" << difftime(end, start) << " seconds " << std::endl;
    }

    






    
    

  



   



    clog<< "\t"<<endl;
    clog<< "------------------ TEST IN 3D ------------------"<<endl;


    {   
        std::clog << "---------------------------------------------"<<std::endl;
        cout<<" --------- Interpolation of Dz  --------- "<< endl;
        //Testing  Q221->Q4 Interpolation
       
        std::cout << "Testing Q221->Q4 Interpolation " << std::endl;

        // Testing correctness
        vector<Real> Vec(18);
        std::srand(unsigned(std::time(nullptr)));
        std::generate(Vec.begin(), Vec.end(), std::rand);

        vector<Real> InterpolationVec_Generic(125);
        GLI_Q221_Q444.Interpolation(Vec,InterpolationVec_Generic);



        Real InterpolationVec_Q241[30];
        Real InterpolationVec_Q441[50];
        Real InterpolationVec_HardCoded[125];


        InterpolationQ221Q241(Vec.data(),InterpolationVec_Q241);
        InterpolationQ241Q441(InterpolationVec_Q241,InterpolationVec_Q441);
        InterpolationQ441Q4(InterpolationVec_Q441,InterpolationVec_HardCoded);

        Real error_square = 0.0;
        Real norm_square = 0.0;
        Real diff = 0.0;
 
        for (Index i=0;i<125;++i)
        {
            diff = InterpolationVec_Generic[i]-InterpolationVec_HardCoded[i];
            error_square += diff*diff;
            norm_square +=  InterpolationVec_Generic[i]*InterpolationVec_Generic[i];
        }

        cout << "Error in the interpolation " << sqrt(error_square/norm_square) << std::endl;

        time_t start, end;
        time(&start);
        for (Index i=0;i<40000000;++i)
        {
            GLI_Q221_Q444.Interpolation(Vec,InterpolationVec_Generic);
            
        }
        time(&end);
        std::cout << "Elapsed time Generic Interpolation :" << difftime(end, start) << " seconds " << std::endl;
        time(&start);
        for (Index i=0;i<40000000;++i)
        {
            InterpolationQ221Q241(Vec.data(),InterpolationVec_Q241);
            InterpolationQ241Q441(InterpolationVec_Q241,InterpolationVec_Q441);
            InterpolationQ441Q4(InterpolationVec_Q441,InterpolationVec_HardCoded);
        }
        time(&end);
    
        std::cout << "Elapsed time Hard Coded Interpolation :" << difftime(end, start) << " seconds " << std::endl; 
    }
 


    

 {   
        std::clog << "---------------------------------------------"<<std::endl;
        cout<<" --------- Transpose Interpolation of Dz  --------- "<< endl;

        //Testing  Q221->Q4 Interpolation
        std::cout << "Testing Q221->Q4 Transpose Interpolation " << std::endl;

        // Testing correctness
        vector<Real> Vec(125);
        std::srand(unsigned(std::time(nullptr)));
        std::generate(Vec.begin(), Vec.end(), std::rand);

        vector<Real> InterpolationVec_Generic(18);
        GLI_Q221_Q444.TransposeInterpolation(Vec,InterpolationVec_Generic);



        Real InterpolationVec_Q241[30];
        Real InterpolationVec_Q441[50];
        Real InterpolationVec_HardCoded[18];

        TransposeInterpolationQ441Q4(Vec.data(), InterpolationVec_Q441);
        TransposeInterpolationQ241Q441(InterpolationVec_Q441, InterpolationVec_Q241);
        TransposeInterpolationQ221Q241(InterpolationVec_Q241, InterpolationVec_HardCoded);
        
        

        Real error_square = 0.0;
        Real norm_square = 0.0;
        Real diff = 0.0 ;
 
        for (Index i=0;i<18;++i)
        {
            
            diff = InterpolationVec_Generic[i]-InterpolationVec_HardCoded[i];
            error_square += diff*diff;
            norm_square +=  InterpolationVec_Generic[i]*InterpolationVec_Generic[i];
        }

        cout << "Error in the interpolation " << sqrt(error_square/norm_square) << std::endl;

        time_t start, end;
        time(&start);
        for (Index i=0;i<20000000;++i)
        {
            GLI_Q221_Q444.TransposeInterpolation(Vec,InterpolationVec_Generic);
            
        }
        time(&end);
        std::cout << "Elapsed time Generic Interpolation :" << difftime(end, start) << " seconds " << std::endl;
        time(&start);
        for (Index i=0;i<20000000;++i)
        {
            TransposeInterpolationQ441Q4(Vec.data(), InterpolationVec_Q441);
            TransposeInterpolationQ241Q441(InterpolationVec_Q441, InterpolationVec_Q241);
            TransposeInterpolationQ221Q241(InterpolationVec_Q241, InterpolationVec_HardCoded);
        }
        time(&end);
    
        std::cout << "Elapsed time Hard Coded Interpolation :" << difftime(end, start) << " seconds " << std::endl; 
    }





    


    {   
        std::clog << "---------------------------------------------"<<std::endl;
        cout<<" --------- Interpolation of Dy  --------- "<< endl;

        //Testing  Q212->Q412->Q414->Q4 Interpolation
 
        std::cout << "Testing Q212->Q4 Interpolation " << std::endl;

        // Testing correctness
        vector<Real> Vec(18);
        std::srand(unsigned(std::time(nullptr)));
        std::generate(Vec.begin(), Vec.end(), std::rand);

        vector<Real> InterpolationVec_Generic(125);
        GLI_Q212_Q444.Interpolation(Vec,InterpolationVec_Generic);



        Real InterpolationVec_Q412[30];
        Real InterpolationVec_Q414[50];
        Real InterpolationVec_HardCoded[125];


        InterpolationQ212Q412(Vec.data(),InterpolationVec_Q412);
        InterpolationQ412Q414(InterpolationVec_Q412,InterpolationVec_Q414);
        InterpolationQ414Q4(InterpolationVec_Q414,InterpolationVec_HardCoded);

        Real error_square = 0.0;
        Real norm_square = 0.0;
        Real diff = 0.0;
 
        for (Index i=0;i<125;++i)
        {
            diff = InterpolationVec_Generic[i]-InterpolationVec_HardCoded[i];
            error_square += diff*diff;
            norm_square +=  InterpolationVec_Generic[i]*InterpolationVec_Generic[i];
        }

        cout << "Error in the interpolation " << sqrt(error_square/norm_square) << std::endl;

        time_t start, end;
        time(&start);
        for (Index i=0;i<20000000;++i)
        {
            GLI_Q212_Q444.Interpolation(Vec,InterpolationVec_Generic);
            
        }
        time(&end);
        std::cout << "Elapsed time Generic Interpolation :" << difftime(end, start) << " seconds " << std::endl;
        time(&start);
        for (Index i=0;i<20000000;++i)
        {
            InterpolationQ212Q412(Vec.data(),InterpolationVec_Q412);
            InterpolationQ412Q414(InterpolationVec_Q412,InterpolationVec_Q414);
            InterpolationQ414Q4(InterpolationVec_Q414,InterpolationVec_HardCoded);
        }
        time(&end);
    
        std::cout << "Elapsed time Hard Coded Interpolation :" << difftime(end, start) << " seconds " << std::endl; 

    }
 


    

    {   
        std::clog << "---------------------------------------------"<<std::endl;
        cout<<" --------- Transpose Interpolation of Dy  --------- "<< endl;


         //Testing  Q212->Q412->Q414->Q4 Interpolation
       
        std::cout << "Testing Q212->Q4 TransposeInterpolation " << std::endl;

        // Testing correctness
        vector<Real> Vec(125);
        std::srand(unsigned(std::time(nullptr)));
        std::generate(Vec.begin(), Vec.end(), std::rand);

        vector<Real> InterpolationVec_Generic(18);
        GLI_Q212_Q444.TransposeInterpolation(Vec,InterpolationVec_Generic);



        Real InterpolationVec_Q412[30];
        Real InterpolationVec_Q414[50];
        Real InterpolationVec_HardCoded[18];

        TransposeInterpolationQ414Q4(Vec.data(), InterpolationVec_Q414);
        TransposeInterpolationQ412Q414(InterpolationVec_Q414, InterpolationVec_Q412);
        TransposeInterpolationQ212Q412(InterpolationVec_Q412, InterpolationVec_HardCoded);
        
        

        Real error_square = 0.0;
        Real norm_square = 0.0;
        Real diff = 0.0 ;
 
        for (Index i=0;i<18;++i)
        {
            
            diff = InterpolationVec_Generic[i]-InterpolationVec_HardCoded[i];
            error_square += diff*diff;
            norm_square +=  InterpolationVec_Generic[i]*InterpolationVec_Generic[i];
        }

        cout << "Error in the interpolation " << sqrt(error_square/norm_square) << std::endl;

         time_t start, end;
        time(&start);
        for (Index i=0;i<20000000;++i)
        {
            GLI_Q212_Q444.TransposeInterpolation(Vec,InterpolationVec_Generic);
            
        }
        time(&end);
        std::cout << "Elapsed time Generic Interpolation :" << difftime(end, start) << " seconds " << std::endl;
        time(&start);
        for (Index i=0;i<20000000;++i)
        {
            TransposeInterpolationQ414Q4(Vec.data(), InterpolationVec_Q414);
            TransposeInterpolationQ412Q414(InterpolationVec_Q414, InterpolationVec_Q412);
            TransposeInterpolationQ212Q412(InterpolationVec_Q412, InterpolationVec_HardCoded);
        
        }
        time(&end);
    
        std::cout << "Elapsed time Hard Coded Interpolation :" << difftime(end, start) << " seconds " << std::endl; 

      
    }



    

    {   
        std::clog << "---------------------------------------------"<<std::endl;
        cout<<" --------- Interpolation of Dx  --------- "<< endl;
        //Testing  Q122->Q142->Q144->Q4 Interpolation
  
        std::cout << "Testing Q122->Q4 Interpolation " << std::endl;

        // Testing correctness
        vector<Real> Vec(18);
        std::srand(unsigned(std::time(nullptr)));
        std::generate(Vec.begin(), Vec.end(), std::rand);

        vector<Real> InterpolationVec_Generic(125);
        GLI_Q122_Q444.Interpolation(Vec,InterpolationVec_Generic);



        Real InterpolationVec_Q142[30];
        Real InterpolationVec_Q144[50];
        Real InterpolationVec_HardCoded[125];


        InterpolationQ122Q142(Vec.data(),InterpolationVec_Q142);
        InterpolationQ142Q144(InterpolationVec_Q142,InterpolationVec_Q144);
        InterpolationQ144Q4(InterpolationVec_Q144,InterpolationVec_HardCoded);

        Real error_square = 0.0;
        Real norm_square = 0.0;
        Real diff = 0.0;
 
        for (Index i=0;i<125;++i)
        {
            diff = InterpolationVec_Generic[i]-InterpolationVec_HardCoded[i];
            error_square += diff*diff;
            norm_square +=  InterpolationVec_Generic[i]*InterpolationVec_Generic[i];
        }

        cout << "Error in the interpolation " << sqrt(error_square/norm_square) << std::endl;

        time_t start, end;
        time(&start);
        for (Index i=0;i<20000000;++i)
        {
            GLI_Q122_Q444.Interpolation(Vec,InterpolationVec_Generic);
            
        }
        time(&end);
        std::cout << "Elapsed time Generic Interpolation :" << difftime(end, start) << " seconds " << std::endl;
        time(&start);
        for (Index i=0;i<20000000;++i)
        {
            InterpolationQ122Q142(Vec.data(),InterpolationVec_Q142);
            InterpolationQ142Q144(InterpolationVec_Q142,InterpolationVec_Q144);
            InterpolationQ144Q4(InterpolationVec_Q144,InterpolationVec_HardCoded);
        }
        time(&end);
    
        std::cout << "Elapsed time Hard Coded Interpolation :" << difftime(end, start) << " seconds " << std::endl; 


     
    }
 


   

 {   
        std::clog << "---------------------------------------------"<<std::endl;
        cout<<" --------- Transpose Interpolation of Dx  --------- "<< endl;

        //Testing  Q122->Q142->Q144->Q4 Interpolation
       
        std::cout << "Testing Q122->Q4 TransposeInterpolation " << std::endl;

        // Testing correctness
        vector<Real> Vec(125);
        std::srand(unsigned(std::time(nullptr)));
        std::generate(Vec.begin(), Vec.end(), std::rand);

        vector<Real> InterpolationVec_Generic(18);
        GLI_Q122_Q444.TransposeInterpolation(Vec,InterpolationVec_Generic);


        Real InterpolationVec_Q142[30];
        Real InterpolationVec_Q144[50];
        Real InterpolationVec_HardCoded[18];

        TransposeInterpolationQ144Q4(Vec.data(), InterpolationVec_Q144);
        TransposeInterpolationQ142Q144(InterpolationVec_Q144, InterpolationVec_Q142);
        TransposeInterpolationQ122Q142(InterpolationVec_Q142, InterpolationVec_HardCoded);
        


        Real error_square = 0.0;
        Real norm_square = 0.0;
        Real diff = 0.0 ;
 
        for (Index i=0;i<18;++i)
        {
            
            diff = InterpolationVec_Generic[i]-InterpolationVec_HardCoded[i];
            error_square += diff*diff;
            norm_square +=  InterpolationVec_Generic[i]*InterpolationVec_Generic[i];
        }

        cout << "Error in the interpolation " << sqrt(error_square/norm_square) << std::endl;


        time_t start, end;
        time(&start);
        for (Index i=0;i<20000000;++i)
        {
            GLI_Q122_Q444.TransposeInterpolation(Vec,InterpolationVec_Generic);
            
        }
        time(&end);
        std::cout << "Elapsed time Generic Interpolation :" << difftime(end, start) << " seconds " << std::endl;
        time(&start);
        for (Index i=0;i<20000000;++i)
        {
            TransposeInterpolationQ144Q4(Vec.data(), InterpolationVec_Q144);
            TransposeInterpolationQ142Q144(InterpolationVec_Q144, InterpolationVec_Q142);
            TransposeInterpolationQ122Q142(InterpolationVec_Q142, InterpolationVec_HardCoded);
        
        }
        time(&end);
    
        std::cout << "Elapsed time Hard Coded Interpolation :" << difftime(end, start) << " seconds " << std::endl; 

    }


 
  

    {
        cout << "---------------------------------------------"<<std::endl;
        cout<< "Testing Boundary Interpolation" <<endl;

     
        std::cout << "Testing Q222->Q444 Interpolation on Boundary " << std::endl;

        // Testing correctness
        vector<Real> Vec(27);
        std::srand(unsigned(std::time(nullptr)));
        std::generate(Vec.begin(), Vec.end(), std::rand);

        vector<Real> InterpolationVec_Generic(125);
        GLI_Q222_Q444.Interpolation(Vec,InterpolationVec_Generic);
       
        // Real InterpolationVec_Q442[75];
        Real InterpolationVec_HardCoded[125];
        // InterpolationQ2Q442Bndy(Vec.data(), InterpolationVec_Q442);
        // InterpolationQ442Q4Bndy(InterpolationVec_Q442, InterpolationVec_HardCoded);
        InterpolationQ222Q444Bndy(Vec.data(),  InterpolationVec_HardCoded);
    
        Real error_square = 0.0;
        Real norm_square = 0.0;
        Real diff = 0.0 ;
        for (Index i=0;i<125;++i)
        {
            if (InterpolationVec_HardCoded[i] > 0.0)
            {
                diff = InterpolationVec_Generic[i]-InterpolationVec_HardCoded[i];
                error_square += diff*diff;
                norm_square +=  InterpolationVec_Generic[i]*InterpolationVec_Generic[i];
            }

        }
        cout << "Error in the interpolation " << sqrt(error_square/norm_square) << std::endl;

        time_t start, end;
        time(&start);
        for (Index i=0;i<20000000;++i)
        {
            GLI_Q222_Q444.Interpolation(Vec,InterpolationVec_Generic);
            
        }
        time(&end);
        std::cout << "Elapsed time Generic Interpolation :" << difftime(end, start) << " seconds " << std::endl;
        time(&start);
        for (Index i=0;i<20000000;++i)
        {
            InterpolationQ222Q444Bndy(Vec.data(),  InterpolationVec_HardCoded);
        }
        time(&end);
    
        std::cout << "Elapsed time Hard Coded Interpolation :" << difftime(end, start) << " seconds " << std::endl; 

        
 
    }
        

    
    {
        cout << "---------------------------------------------"<<std::endl;
        cout<< "Testing Boundary TransposeInterpolation" <<endl; 

         
        std::cout << "Testing Q2->Q4 TransposeInterpolation on Boundary " << std::endl;

        // Testing correctness
        vector<Real> Vec(125);
        std::srand(unsigned(std::time(nullptr)));
        std::generate(Vec.begin(), Vec.end(), std::rand);

        Vec[31]=0.0;
        Vec[32]=0.0;
        Vec[33]=0.0;

        Vec[36]=0.0;
        Vec[37]=0.0;
        Vec[38]=0.0;

        Vec[41]=0.0;
        Vec[42]=0.0;
        Vec[43]=0.0;

        Vec[56]=0.0;
        Vec[57]=0.0;
        Vec[58]=0.0;

        Vec[61]=0.0;
        Vec[62]=0.0;
        Vec[63]=0.0;

        Vec[66]=0.0;
        Vec[67]=0.0;
        Vec[68]=0.0;


        Vec[81]=0.0;
        Vec[82]=0.0;
        Vec[83]=0.0;

        Vec[86]=0.0;
        Vec[87]=0.0;
        Vec[88]=0.0;

        Vec[91]=0.0;
        Vec[92]=0.0;
        Vec[93]=0.0;

        vector<Real> InterpolationVec_Generic(27);
        GLI_Q222_Q444.TransposeInterpolation(Vec,InterpolationVec_Generic);
       
        // Real InterpolationVec_Q442[75];
        Real InterpolationVec_HardCoded[27];
        // TransposeInterpolationQ442Q4Bndy(Vec.data(), InterpolationVec_Q442);
        // TransposeInterpolationQ2Q442Bndy(InterpolationVec_Q442, InterpolationVec_HardCoded);
        TransposeInterpolationQ222Q444Bndy(Vec.data(),  InterpolationVec_HardCoded);
        Real error_square = 0.0;
        Real norm_square = 0.0;
        Real diff = 0.0;
       
        for (Index i=0;i<27;++i)
        {
            
            diff = InterpolationVec_Generic[i]-InterpolationVec_HardCoded[i];
            error_square += diff*diff;
            norm_square +=  InterpolationVec_Generic[i]*InterpolationVec_Generic[i];
        }

        cout << "Error in the interpolation " << sqrt(error_square/norm_square) << std::endl;


        time_t start, end;
        time(&start);
        for (Index i=0;i<40000000;++i)
        {
            GLI_Q222_Q444.TransposeInterpolation(Vec,InterpolationVec_Generic);
            
        }
        time(&end);
        std::cout << "Elapsed time Generic Interpolation :" << difftime(end, start) << " seconds " << std::endl;
        time(&start);
        for (Index i=0;i<40000000;++i)
        {
           TransposeInterpolationQ222Q444Bndy(Vec.data(),  InterpolationVec_HardCoded);
        }
        time(&end);
    
        std::cout << "Elapsed time Hard Coded Interpolation :" << difftime(end, start) << " seconds " << std::endl; 

    }

  
 
 
    


}