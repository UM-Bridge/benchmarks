#pragma once


// -----------------------------------------------------------------------------------------//
// Class definition
namespace OndoMathX {

namespace GGL {


template <
size_t PX,
size_t PY,
size_t PZ = 0,
>
class Laplace {
    
public:
    
    // ---------------------------------------------------------------------------------//
    enum BC
    {
        Dirichlet = -1,
        Periodic  = 0,
        Neumann   = 1
    };
    
    // ---------------------------------------------------------------------------------//
    /*! \brief Deleted default constructor.
     */
    Laplace() = delete;
    
    
    
    /*! \brief Constructor
     \param aFESpace is a pointer to a finite element space the mass operator is built upon.
     \param Density is a density used to define the mass operator.
     */
    Laplace(double alpha,
            double LX,
            Index NX,
            BC BC_X,
            double LY,
            Index NY,
            BC BC_Y,
            double LZ=1.0,
            Index NZ=0,
            BC BC_Z=Periodic) :  _alpha(alpha),
                                  _LX(LX), _NX(NX), _BC_X(BC_X),
                                  _LY(LY), _NY(NY), _BC_Y(BC_Y),
                                  _LZ(LZ), _NZ(NZ), _BC_Z(BC_Z)
    {
        //Double the size of the computationnal domain to take into account B.C.
        if (BC_X != Periodic) {_NX*=2;_LX*=2;}
        if (BC_Y != Periodic) {_NY*=2;_LY*=2;}
        if (BC_Z != Periodic) {_NZ*=2;_LZ*=2;}
        
        //Compute the symbols for each frequencies
        Line<PX, 1> unitLine_X(1);
        Line<PY, 1> unitLine_Y(1);
        Line<PZ, 1> unitLine_Z(1);
        
        Field::Identity Id;
        
        LAL::SparseMatrix stiffnessMatrixX;
        LAL::SparseMatrix stiffnessMatrixY;
        LAL::SparseMatrix stiffnessMatrixZ;
        
        LAL::DiagonalMatrix massMatrixX;
        LAL::DiagonalMatrix massMatrixY;
        LAL::DiagonalMatrix massMatrixZ;
        
        LAL::ComplexDenseMatrix symbolMassX;
        LAL::ComplexDenseMatrix symbolMassY;
        LAL::ComplexDenseMatrix symbolMassZ;
        
        LAL::DiagonalMatrix symbolMass;
        
        std::vector<LAL::ComplexDenseMatrix> symbolX(NX);
        std::vector<LAL::ComplexDenseMatrix> symbolY(NY);
        std::vector<LAL::ComplexDenseMatrix> symbolZ(NZ);
        
        std::vector<std:vector<ComplexVector>> eigenvectorsX(NX);
        std::vector<std:vector<ComplexVector>> eigenvectorsY(NY);
        std::vector<std:vector<ComplexVector>> eigenvectorsZ(NZ);
        
        std::vector<std:vector<Real>> eigenvaluesX(NX);
        std::vector<std:vector<Real>> eigenvaluesY(NY);
        std::vector<std:vector<Real>> eigenvaluesZ(NZ);
        
        
        CreateStiffness(unitLine_X, Id, stiffnessMatrixX, 1.0);
        CreateStiffness(unitLine_Y, Id, stiffnessMatrixY, 1.0);
        if (PZ != 0) CreateStiffness(unitLine_Z, Id, stiffnessMatrixZ, 1.0);
        
        CreateMass(unitLine_X, Id, massMatrixX, 1.0);
        CreateMass(unitLine_Y, Id, massMatrixY, 1.0);
        if (PZ != 0) CreateMass(unitLine_Z, Id, massMatrixZ, 1.0);
        
        LAL::Allocate(symbolMassX,PX);
        LAL::Allocate(symbolMassY,PY);
        if (PZ != 0) LAL::Allocate(symbolMassZ,PZ);
        
        // symbol of the mass matrix in the X direction
        for (Index i = 0; i < PX; i++)
        {
            Real mii = LAL::getValue(MassMatrixX,i,i);
            
            if (i==0)
                LAL::setValue(symbolMassX,i,i,2*mii);
            else
                LAL::setValue(symbolMassX,i,i,mii);
        }
        
        // symbol of the mass matrix in the Y direction
        for (Index i = 0; i < PY; i++)
        {
            Real mii = LAL::getValue(MassMatrixY,i,i);
            
            if (i==0)
                LAL::setValue(symbolMassY,i,i,2*mii);
            else
                LAL::setValue(symbolMassY,i,i,mii);
        }
        
        // symbol of the mass matrix in the Z direction
        if (PZ != 0) for (Index i = 0; i < PZ; i++)
        {
            Real mii = LAL::getValue(MassMatrixZ,i,i);
            
            if (i==0)
                LAL::setValue(symbolMassZ,i,i,2*mii);
            else
                LAL::setValue(symbolMassZ,i,i,mii);
        }
        
    
        Complex value;
    
        // Symbol of the Laplace operator in frequency domain in the X direction
        for (Index f = 0; f < NX; f++ )
        {
            LAL::Allocate(symbolX[f],PX,PX);
            
            Real theta = f*2.0*M_PI/NX;
            Real cos_theta = cos(theta)
            Real sin_theta = sin(theta)
            Complex ei_theta = cos_theta + Complex(0.0,1.0)*sin_theta;

            value = 2.0*LAL::getValue(stiffnessMatrixX,0,0) + 2.0*LAL::getValue(stiffnessMatrixX,0,PX)*cos_theta;
            
            setValue(symbolX[f],0,0,value);
            
            for (Index i = 1; i < PX; i++)
            {
                value = LAL::getValue(stiffnessMatrixX,0,i)+LAL::getValue(stiffnessMatrixX,PX,i)*ei_theta;
                
                LAL::setValue(symbolX[f],0,i,value);
                
                LAL::setValue(symbolX[f],i,0,conj(value));
            }
            
            for (Index i = 1; i < PX; i++)
                for (Index j = 1; j < PX; j++)
                    LAL::setValue(symbolX[f],i,j,LAL::getValue(stiffnessMatrixX,i,j));
        }
            
        // Symbol of the Laplace operator in frequency domain in the Y direction
        for (Index f = 0; f < NY; f++ )
        {
            LAL::Allocate(symbolY[f],PY,PY);
            
            Real theta = f*2.0*M_PI/NY;
            Real cos_theta = cos(theta)
            Real sin_theta = sin(theta)
            Complex ei_theta = cos_theta + Complex(0.0,1.0)*sin_theta;

            value = 2.0*LAL::getValue(stiffnessMatrixY,0,0) + 2.0*LAL::getValue(stiffnessMatrixY,0,PY)*cos_theta;
            
            setValue(symbolY[f],0,0,value);
            
            for (Index i = 1; i < PY; i++)
            {
                value = LAL::getValue(stiffnessMatrixY,0,i)+LAL::getValue(stiffnessMatrixY,PY,i)*ei_theta;
                
                LAL::setValue(symbolY[f],0,i,value);
                
                LAL::setValue(symbolY[f],i,0,conj(value));
            }
            
            for (Index i = 1; i < PY; i++)
                for (Index j = 1; j < PY; j++)
                    LAL::setValue(symbolY[f],i,j,LAL::getValue(stiffnessMatrixY,i,j));
        }
        
        // Symbol of the Laplace operator in frequency domain in the Z direction
        if (PZ != 0) for (Index f = 0; f < NZ; f++ )
        {
            LAL::Allocate(symbolZ[f],PZ,PZ);
            
            Real theta = f*2.0*M_PI/NZ;
            Real cos_theta = cos(theta)
            Real sin_theta = sin(theta)
            Complex ei_theta = cos_theta + Complex(0.0,1.0)*sin_theta;

            value = 2.0*LAL::getValue(stiffnessMatrixZ,0,0) + 2.0*LAL::getValue(stiffnessMatrixZ,0,PZ)*cos_theta;
            
            setValue(symbolZ[f],0,0,value);
            
            for (Index i = 1; i < PZ; i++)
            {
                value = LAL::getValue(stiffnessMatrixZ,0,i)+LAL::getValue(stiffnessMatrixZ,PZ,i)*ei_theta;
                
                LAL::setValue(symbolZ[f],0,i,value/LAL::getValue(symbolMassMatrixZ,0,0));
                
                LAL::setValue(symbolZ[f],i,0,conj(value)/LAL::getValue(symbolMassMatrixZ,i,i));
            }
            
            for (Index i = 1; i < PZ; i++)
                for (Index j = 1; j < PZ; j++)
                    LAL::setValue(symbolZ[f],i,j,LAL::getValue(stiffnessMatrixZ,i,j));
        }
        
        
        // Compute eigenvalues and eigenvectors of the symbols
        for (Index f = 0; f < NX; f++ )
            LAL::GeneralizedSelfAdjointEigenvalueDecomposition(symbolX[f],symbolMassX,eigenvaluesX[f],eigenvectorsX[f]);
    
        for (Index f = 0; f < NY; f++ )
            LAL::GeneralizedSelfAdjointEigenvalueDecomposition(symbolY[f],symbolMassY,eigenvaluesY[f],eigenvectorsY[f]);
    
        if (PZ != 0) for (Index f = 0; f < NZ; f++ )
            LAL::GeneralizedSelfAdjointEigenvalueDecomposition(symbolZ[f],symbolMassZ,eigenvaluesZ[f],eigenvectorsZ[f]);
    
        //Construct the symbol of the 2d or 3d mass matrix
        if (PZ != 0)
        {
            LAL::Allocate(symbolMass, PX*PY*PZ);
            
            for (Index i = 0; i < PX; i++)
                for (Index j = 0; j < PY; j++)
                    for (Index k = 0; j < PZ; k++)
                    {
                        Index m = k*PX*PY+j*PX+i;
                    
                        LAL::setValue(symbolMass,m,m, LAL::getValue(symbolMassX,i)
                                                     *LAL::getValue(symbolMassY,j)
                                                     *LAL::getValue(symbolMassZ,k));
                    }
            
            
        }
        else
        {
            LAL::Allocate(symbolMass, PX*PY);
            
            for (Index i = 0; i < PX; i++)
                for (Index j = 0; j < PY; j++)
                {
                    Index m = j*PX+i;
                    
                    LAL::setValue(symbolMass,m,m, LAL::getValue(symbolMassX,i)
                                                 *LAL::getValue(symbolMassY,j));
                }
        }
        
        

        
        
        
        
    }
    
    
    
    
    /*! \brief Destructor
     */
    virtual ~Laplace() = default;
    // ---------------------------------------------------------------------------------//
    
    //virtual void Assemble_SymbolM_2D(size_t PX, Symbol1D_M,Symbol2D_M) // define type
    virtual void Init_for_Solve()
    {
        

        Symbol1D<PX, LAL> symbX(_LX,_NX,1);
        Symbol1D<PX, LAL> symbY(_LY,_NY,1);
        
        //compute index ilocMin = IndexMinX + PX*IndexMinY such that associated eigenvalue is equal to zero (frequency K = (0,0) )
        int IndexMinX = 0;
        for (Index iX = 1; iX < PX; iX ++)
        {
            if (real(LAL::Base::Value(_Symbol1DX->GetEigenvalues(0), iX)) < real(LAL::Base::Value(_Symbol1DX->GetEigenvalues(0), IndexMinX)))
                IndexMinX = iX;
        }
        int IndexMinY = 0;
        for (Index iY = 1; iY < PX; iY ++)
        {
            if (real(LAL::Base::Value(_Symbol1DY->GetEigenvalues(0), iY)) < real(LAL::Base::Value(_Symbol1DY->GetEigenvalues(0), IndexMinY)))
                IndexMinY = iY;
        }
        Index ilocMin = IndexMinX + PX*IndexMinY;
        
        // compute eigenvectors U and eigenvalues in 2D or 3D,
        // compute _lambda = eigenval * U'*Symbol2DM*U for symbol inversion
        ComplexVector U_tmp_1;
        LAL::Base::Allocate(U_tmp_1,PX*PX);
        Index iFreq = 0;
         _lambda.resize(_NX*_NY);
        for (Index iFreqY = 0; iFreqY < _NY; iFreqY ++ )
        {
            for (Index iFreqX = 0; iFreqX < _NX; iFreqX ++ )
            {
                iFreq = iFreqX + _NX*iFreqY;
                LAL::Base::Allocate(_lambda[iFreq],PX*PX); //lambda = eigenval * U'*Symbol2DM*U, where eigenval = eigx+eigy+alpha
                Complex * _lambda_data = LAL::Base::GetData(_lambda[iFreq]);
                Complex * U_tmp_1_data = LAL::Base::GetData(U_tmp_1);
                Real * Mass_symbol = LAL::Base::GetData(_Symbol2DMassMatrix);
                
                std::complex<double> temp_eigenvalue_y = 0.;
                std::complex<double> temp_eigenvalue_xy = 0.;
                
                for (Index i1 = 0; i1 < PX; ++i1)
                    for (Index j1 = 0; j1 < PX; ++j1)
                    {
                        Index iloc = i1 + PX*j1;
                        U_tmp_1_data[iloc]=0.0;
                        _lambda_data[iloc] = 0.0;
                    }
                // First projection.
                for (Index i1 = 0; i1 < PX; ++i1)
                {
                    Complex * eigenvector_x = LAL::Base::GetData(_Symbol1DX->GetEigenvectors(iFreqX,i1));
                    
                    for (Index j1 = 0; j1 < PX; ++j1)
                        for (Index i2 = 0; i2 < PX; ++i2)
                            U_tmp_1_data[i1 + PX*j1] += Mass_symbol[i2 + PX*j1] * conj(eigenvector_x[i2])*eigenvector_x[i2];
                    
                }
                
                // Second projection.
                for (Index j1 = 0; j1 < PX; ++j1)
                {
                    Complex * eigenvector_y = LAL::Base::GetData(_Symbol1DY->GetEigenvectors(iFreqY,j1));
                    
                    for (Index i1 = 0; i1 < PX; ++i1)
                        for (Index j2 = 0; j2 < PX; ++j2)
                            _lambda_data[i1 + PX*j1] += U_tmp_1_data[i1 + PX*j2] * conj(eigenvector_y[j2])*eigenvector_y[j2];
                }
                
                //Scale the coefficient of the decomposition
                for (Index j1 = 0; j1 < PX; ++j1)
                {
                    temp_eigenvalue_y = _alpha * 1.0 + LAL::Base::Value(_Symbol1DY->GetEigenvalues(iFreqY), j1);// we solve (_alpha *Id - Delta ) u = f
                    for (Index i1 = 0; i1 < PX; ++i1)
                    {
                        temp_eigenvalue_xy = temp_eigenvalue_y + LAL::Base::Value(_Symbol1DX->GetEigenvalues(iFreqX), i1);
                        Index iloc = i1 + PX*j1;
//                                double temp = real(LAL::Base::Value(_Symbol1DY->GetEigenvalues(iFreqY), j1)+ LAL::Base::Value(_Symbol1DX->GetEigenvalues(iFreqX), i1));
//                                std::cout<< iFreqX << " " << iFreqY << " " << temp << std::endl;
                        if (iFreq == 0 && iloc == ilocMin && _alpha == 0)
                        {
//                                     std::cout <<"temp eig: " << temp_eigenvalue_xy <<std::endl;
                            _lambda_data[iloc] = 0; // no contribution if alpha = 0 and for frequency = (0,0) and (i1,j1) such that eigenvalue_{k,iloc} =0 -> average of the solution set to zero
                        }
                        else
                        {
                            _lambda_data[iloc] = 1.0/ (_lambda_data[iloc]*temp_eigenvalue_xy);
                        }
//                                std::cout <<"temp eig: " << temp_eigenvalue_xy <<std::endl;
//                                std::cout <<"temp eig: " << temp_eigenvalue_xy << " " <<"lambda: " << _lambda_data[iloc] << std::endl;

                    }
                }
            }
        }
    }
    
    
    
    
    virtual void Solve(Vector & refU) // define type
    {
        //cout << "Total NDof: " <<  PX*PX*_NX * _NY << endl;
        Vector U_temp(PX*PX*_NX * _NY);
        std::vector<Vector> RealUbyOrder(PX*PX);
        std::vector<Vector> ImagUbyOrder(PX*PX);
        _HODFTSolutionVector.resize(PX*PX);
        
        //DST of right-hand-side by order
        Index iloc = 0;
        //size_t NElem = 0;
        //size_t index = 0;
        for (Index ilocy = 0; ilocy < PX; ilocy++)
        {
            for (Index ilocx = 0; ilocx < PX; ilocx++)
            {
                iloc = ilocx + PX * ilocy ;
                LAL::Base::Allocate(RealUbyOrder[iloc],_NX * _NY);
                LAL::Base::Allocate(ImagUbyOrder[iloc],_NX * _NY);
                LAL::Base::Allocate(_HODFTSolutionVector[iloc],_NX * _NY);
            }
        }
        switch(_BC)
        {
            case -1 :
            { // Dirichlet
                LAL::SymF(refU,U_temp,_NX/2,_NY/2,PX,-1);
                break;
            }
            case 0 :
            { // Periodic
                LAL::PerF(refU,U_temp,_NX,_NY,PX);
                break;
            }
            case 1 :
            { // Neumann
                //Vector F_sym;
                LAL::SymF(refU,U_temp,_NX/2,_NY/2,PX,1);
                break;
            }
            default :  std::cout << "Boundary condition not correct\n";
                break;
        }
        
        
        // Compute DFT in 2D (by order)
        _ThreadWiseDFT2D(0, PX*PX,  U_temp, RealUbyOrder, ImagUbyOrder);

        // Compute inverse symbol of the solution in Fourier space
        _ThreadWise2DSymbol(0,  _NX*_NY, RealUbyOrder, ImagUbyOrder, _HODFTSolutionVector);
    
        // Compute inverse DFT in 2D (by order)
        _ThreadWiseIDFT2D(0, PX*PX,  _HODFTSolutionVector, U_temp);

        
        switch(_BC)
        {
            case -1 :
            { // Dirichlet
                LAL::GetSol(U_temp,refU,_NX/2,_NY/2,PX,-1);
                break;
            }
            case 0 :
            { // Periodic
                LAL::GetSol(U_temp,refU,_NX,_NY,PX,0);
                break;
            }
            case 1 :
            { // Neumann
                LAL::GetSol(U_temp,refU,_NX/2,_NY/2,PX,1);
                break;
            }
            default :  std::cout << "Boundary condition not correct\n";
                break;
        }
        
        
        
    }

        
    // ---------------------------------------------------------------------------------//
    void _ThreadWise2DSymbol(
                             Index iFreqBegin,
                             Index iFreqEnd,
                             std::vector<Vector> & RealUbyOrder,
                             std::vector<Vector> & ImagUbyOrder,
                             std::vector<ComplexVector> & SolutionVector)
    {
        
        
        
        //compute DST(solution) = U' * Mfr * U / _lambda
        ComplexVector U_tmp_1;
        ComplexVector U_tmp_2;
        ComplexVector RHS;
        
        LAL::Base::Allocate(U_tmp_1,PX*PX);
        LAL::Base::Allocate(U_tmp_2,PX*PX);
        LAL::Base::Allocate(RHS,PX*PX);
        
        for (Index iFreq = iFreqBegin; iFreq < iFreqEnd; iFreq++)
        {
            Index iFreqX = iFreq % _NX;
            Index iFreqY = iFreq / _NX;
            
            Complex * RHS_data = LAL::Base::GetData(RHS);
            Complex * U_tmp_1_data = LAL::Base::GetData(U_tmp_1);
            Complex * U_tmp_2_data = LAL::Base::GetData(U_tmp_2);
            Complex * lambda_data = LAL::Base::GetData(_lambda[iFreq]);
            
            Real * Symbol2DMassMatrix_data = LAL::Base::GetData(_Symbol2DMassMatrix);
            
            //Construct the complex vector from real and imaginary part and apply Symbol2DM * DST(RHS)
            for (Index i1 = 0; i1 < PX; ++i1)
                for (Index j1 = 0; j1 < PX; ++j1)
                {
                    Index iloc = i1 + PX*j1;
                    Complex tmp =  complex<double>(LAL::Base::Value(RealUbyOrder[iloc],iFreq), LAL::Base::Value(ImagUbyOrder[iloc],iFreq));
                    RHS_data[iloc] = Symbol2DMassMatrix_data[iloc] * tmp;
                    U_tmp_1_data[iloc]=0.0;
                    U_tmp_2_data[iloc]=0.0;
                }
            
            
            // First projection.
            for (Index i1 = 0; i1 < PX; ++i1)
            {
                Complex * eigenvector_x = LAL::Base::GetData(_Symbol1DX->GetEigenvectors(iFreqX,i1));
                
                for (Index j1 = 0; j1 < PX; ++j1)
                    for (Index i2 = 0; i2 < PX; ++i2)
                        U_tmp_1_data[i1 + PX*j1] += RHS_data[i2 + PX*j1] * conj(eigenvector_x[i2]);
                
            }
            
            // Second projection.
            for (Index j1 = 0; j1 < PX; ++j1)
            {
                Complex * eigenvector_y = LAL::Base::GetData(_Symbol1DY->GetEigenvectors(iFreqY,j1));
                
                for (Index i1 = 0; i1 < PX; ++i1)
                    for (Index j2 = 0; j2 < PX; ++j2)
                        U_tmp_2_data[i1 + PX*j1] += U_tmp_1_data[i1 + PX*j2] * conj(eigenvector_y[j2]);
            }
            
            
            
            //Scale the coefficient of the decomposition
            for (Index i1 = 0; i1 < PX; ++i1)
                for (Index j1 = 0; j1 < PX; ++j1)
                {
                    Index iloc = i1 + PX*j1;
                    U_tmp_2_data[iloc]*=lambda_data[iloc];
                    U_tmp_1_data[iloc]=0.0;
                    RHS_data[iloc]=0.0;
                }
            
            for (Index i1 = 0; i1 < PX; ++i1)
                for (Index j1 = 0; j1 < PX; ++j1)
                {
                    
                }
            
            // Second projection transpose.
            for (Index j2 = 0; j2 < PX; ++j2)
            {
                Complex * eigenvector_y = LAL::Base::GetData(_Symbol1DY->GetEigenvectors(iFreqY,j2));
                
                for (Index i1 = 0; i1 < PX; ++i1)
                    for (Index j1 = 0; j1 < PX; ++j1)
                        U_tmp_1_data[i1 + PX*j1] += U_tmp_2_data[i1 + PX*j2] * eigenvector_y[j1];
                
            }
            
            
            // First projection transpose.
            for (Index i2 = 0; i2 < PX; ++i2)
            {
                Complex * eigenvector_x = LAL::Base::GetData(_Symbol1DX->GetEigenvectors(iFreqX,i2));
                
                for (Index j1 = 0; j1 < PX; ++j1)
                    for (Index i1 = 0; i1 < PX; ++i1)
                        RHS_data[i1 + PX*j1] += U_tmp_1_data[i2 + PX*j1] * eigenvector_x[i1];
            }
            
            // Store solution by frequence.
            for (Index i1 = 0; i1 < PX; ++i1)
                for (Index j1 = 0; j1 < PX; ++j1)
                {
                    Index iloc = i1 + PX*j1;
                    LAL::Base::Value(SolutionVector[iloc],iFreq) = RHS_data[iloc];
                }
            
        }
    }
    // ---------------------------------------------------------------------------------//
    void _ThreadWise2DSymbol(
                             Index iFreqBegin,
                             Index iFreqEnd,
                             std::vector<ComplexVector> & UbyOrder,
                             std::vector<ComplexVector> & SolutionVector)
    {
        
        
        
        //compute DST(solution) = U' * Mfr * U / _lambda
        ComplexVector U_tmp_1;
        ComplexVector U_tmp_2;
        ComplexVector RHS;
        
        LAL::Base::Allocate(U_tmp_1,PX*PX);
        LAL::Base::Allocate(U_tmp_2,PX*PX);
        LAL::Base::Allocate(RHS,PX*PX);
        
        for (Index iFreq = iFreqBegin; iFreq < iFreqEnd; iFreq++)
        {
            Index iFreqX = iFreq % _NX;
            Index iFreqY = iFreq / _NX;
            
            Complex * RHS_data = LAL::Base::GetData(RHS);
            Complex * U_tmp_1_data = LAL::Base::GetData(U_tmp_1);
            Complex * U_tmp_2_data = LAL::Base::GetData(U_tmp_2);
            Complex * lambda_data = LAL::Base::GetData(_lambda[iFreq]);
            
            Real * Symbol2DMassMatrix_data = LAL::Base::GetData(_Symbol2DMassMatrix);
            
            //Construct the complex vector from real and imaginary part and apply Symbol2DM * DST(RHS)
            for (Index i1 = 0; i1 < PX; ++i1)
                for (Index j1 = 0; j1 < PX; ++j1)
                {
                    Index iloc = i1 + PX*j1;
                    Complex tmp = LAL::Base::Value(UbyOrder[iloc],iFreq);
                    RHS_data[iloc] = Symbol2DMassMatrix_data[iloc] * tmp;
                    U_tmp_1_data[iloc]=0.0;
                    U_tmp_2_data[iloc]=0.0;
                }
            
            
            // First projection.
            for (Index i1 = 0; i1 < PX; ++i1)
            {
                Complex * eigenvector_x = LAL::Base::GetData(_Symbol1DX->GetEigenvectors(iFreqX,i1));
                
                for (Index j1 = 0; j1 < PX; ++j1)
                    for (Index i2 = 0; i2 < PX; ++i2)
                        U_tmp_1_data[i1 + PX*j1] += RHS_data[i2 + PX*j1] * conj(eigenvector_x[i2]);
                
            }
            
            // Second projection.
            for (Index j1 = 0; j1 < PX; ++j1)
            {
                Complex * eigenvector_y = LAL::Base::GetData(_Symbol1DY->GetEigenvectors(iFreqY,j1));
                
                for (Index i1 = 0; i1 < PX; ++i1)
                    for (Index j2 = 0; j2 < PX; ++j2)
                        U_tmp_2_data[i1 + PX*j1] += U_tmp_1_data[i1 + PX*j2] * conj(eigenvector_y[j2]);
            }
            
            
            
            //Scale the coefficient of the decomposition
            for (Index i1 = 0; i1 < PX; ++i1)
                for (Index j1 = 0; j1 < PX; ++j1)
                {
                    Index iloc = i1 + PX*j1;
                    U_tmp_2_data[iloc]*=lambda_data[iloc];
                    U_tmp_1_data[iloc]=0.0;
                    RHS_data[iloc]=0.0;
                }
            
            for (Index i1 = 0; i1 < PX; ++i1)
                for (Index j1 = 0; j1 < PX; ++j1)
                {
                    
                }
            
            // Second projection transpose.
            for (Index j2 = 0; j2 < PX; ++j2)
            {
                Complex * eigenvector_y = LAL::Base::GetData(_Symbol1DY->GetEigenvectors(iFreqY,j2));
                
                for (Index i1 = 0; i1 < PX; ++i1)
                    for (Index j1 = 0; j1 < PX; ++j1)
                        U_tmp_1_data[i1 + PX*j1] += U_tmp_2_data[i1 + PX*j2] * eigenvector_y[j1];
                
            }
            
            
            // First projection transpose.
            for (Index i2 = 0; i2 < PX; ++i2)
            {
                Complex * eigenvector_x = LAL::Base::GetData(_Symbol1DX->GetEigenvectors(iFreqX,i2));
                
                for (Index j1 = 0; j1 < PX; ++j1)
                    for (Index i1 = 0; i1 < PX; ++i1)
                        RHS_data[i1 + PX*j1] += U_tmp_1_data[i2 + PX*j1] * eigenvector_x[i1];
            }
            
            // Store solution by frequence.
            for (Index i1 = 0; i1 < PX; ++i1)
                for (Index j1 = 0; j1 < PX; ++j1)
                {
                    Index iloc = i1 + PX*j1;
                    LAL::Base::Value(SolutionVector[iloc],iFreq) = RHS_data[iloc];
                }
            
        }
    }
  
    void _ThreadWiseIDFT2D(
                           Index iLocBegin,
                           Index iLocEnd,
                           std::vector<ComplexVector> & HODFTSolVector,
                           ComplexVector & SolVector)
    {
        for (Index iLoc = iLocBegin; iLoc < iLocEnd; iLoc ++)
        {
            Index iLocX = iLoc % PX;
            Index iLocY = iLoc / PX;
            
            LAL::IDFT2D_Complex(HODFTSolVector[iLoc],_NX,_NY);
            Index NElem = 0;
            for (Index ny = 0; ny < _NY; ny++)
            {
                for (Index nx = 0; nx < _NX; nx++)
                {
                    Index index = (ny*PX + iLocY) * _NX*PX + (nx * PX + iLocX); // global dof
                    LAL::Base::Value(SolVector,index) = LAL::Base::Value(HODFTSolVector[iLoc],NElem) ;
                    NElem++;
                }
            }
        }
    }
    
    void _ThreadWiseDFT2D(
                          Index iLocBegin,
                          Index iLocEnd,
                          ComplexVector & SolVector,
                          std::vector<ComplexVector> & DFTSolVector
                          )
    {
        for (Index iLoc = iLocBegin; iLoc < iLocEnd; iLoc ++)
        {
            Index iLocX = iLoc % PX;
            Index iLocY = iLoc / PX;
            
            Index NElem = 0;
            for (Index ny = 0; ny < _NY; ny++)
            {
                for (Index nx = 0; nx < _NX; nx++)
                {
                    Index index = (ny*PX + iLocY) * _NX*PX + (nx * PX + iLocX); // global dof
                    LAL::Base::Value(DFTSolVector[iLoc],NElem) = LAL::Base::Value(SolVector,index);
                    NElem++;
                }
            }
            LAL::DFT2D(DFTSolVector[iLoc],_NX,_NY);
        }
    }
    
    // ---------------------------------------------------------------------------------//
    
private:
    
    /*! \brief Number of element in the x-direction. */
    double _LX;
    double _LY;
    Index _NX;
    Index _NY;
    double _alpha;
    int _BC;
    
    std::unique_ptr<Symbol1D<PX,LAL>> _Symbol1DX;
    std::unique_ptr<Symbol1D<PX,LAL>> _Symbol1DY;
    
    Vector _Symbol2DMassMatrix ;
    
    std::vector<ComplexVector> _lambda ;
    
    std::vector<ComplexVector> _HODFTSolutionVector ;
};

 
// -----------------------------------------------------------------------------------------//



} //GGL

} //OndoMathX
// -----------------------------------------------------------------------------------------//


