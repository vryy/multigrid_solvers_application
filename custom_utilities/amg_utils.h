/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
 */


/* *********************************************************
 *
 *   Last Modified by:    $Author: hbui $
 *   Date:                $Date: 2013 9 Jan 17:51:00 $
 *   Revision:            $Revision: 1.1 $
 *
 * ***********************************************************/


#if !defined(KRATOS_MULTIGRID_SOLVERS_APP_AMG_UTILS )
#define  KRATOS_MULTIGRID_SOLVERS_APP_AMG_UTILS


/* System includes */
#include <algorithm>
//#include <random>


/* External includes */
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include "external_includes/pyamg/ruge_stuben.h"
#include "external_includes/pyamg/smoothed_aggregation.h"
#include "external_includes/pyamg/graph.h"
#include "external_includes/scipy/sparse/sparsetools/csr.h"


/* Project includes */
#include "includes/define.h"
#include "utilities/timer.h"


namespace Kratos
{

/**@name Kratos Globals */
/*@{ */


/*@} */
/**@name Type Definitions */
/*@{ */


/*@} */


/**@name  Enum's */
/*@{ */


/*@} */
/**@name  Functions */
/*@{ */



/*@} */
/**@name Kratos Classes */
/*@{ */

/**
 * This class defines utility for aggregation of node clusters to be used in deflated solvers.
 */
template<class TSparseSpaceType>
class AMGUtils
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(AMGUtils);
    
    
    /**@name constant Definitions */
    /*@{ */

    static const int COLORING_NONE = 0;
    static const int COLORING_MIS  = 1;
    static const int COLORING_JP   = 2;
    static const int COLORING_LDF  = 3;

    /**@name Type Definitions */
    /*@{ */

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::MatrixPointerType SparseMatrixPointerType;

    typedef std::size_t SizeType;

    typedef std::size_t IndexType;

    typedef double ValueType;

    typedef boost::numeric::ublas::unbounded_array<IndexType> IndexContainerType;

    typedef boost::numeric::ublas::unbounded_array<int> IntegerContainerType;

    typedef boost::numeric::ublas::unbounded_array<ValueType> ValueContainerType;

    typedef boost::numeric::ublas::coordinate_matrix<ValueType> LocalSparseMatrixType;

//    typedef boost::shared_ptr<LocalSparseMatrixType> LocalSparseMatrixPointerType;

//    typedef std::vector<IndexType> IndexVectorType;
//    typedef boost::numeric::ublas::vector<IndexType> IndexVectorType;
    typedef boost::numeric::ublas::vector<int> IndexVectorType; // this must be int since this is the only one support by kratos python

//    typedef boost::shared_ptr<IndexVectorType> IndexVectorPointerType;

    /*@} */
    /**@name Life Cycle
     */
    /*@{ */

    /** Constructor.
     */


    /** Destructor.
     */

    /*@} */
    /**@name Operators
     */


    //////////////////////////////////////////////////////////////////
    //      MATRIX GALLERY
    //////////////////////////////////////////////////////////////////
    // helper function to generate the matrix from Poisson problem
    // The matrix represents a finite Difference approximation to the 
    // Poisson problem on a regular n-dimensional grid with unit grid 
    // spacing and Dirichlet boundary conditions.
    // reference: pyamg/laplacian.py/poisson
    //this matrix is m x m matrix of n x n blocks
    static SparseMatrixPointerType Poisson(const SizeType m, const SizeType n)
    {
        LocalSparseMatrixType A(m*n, m*n, (n+2*(n-1))*m + n*(m-1)*2 );

        for(SizeType i = 0; i < m; i++)
        {
            // center blocks
            SizeType iglobal = i * n;
            SizeType jglobal = i * n;
            for(SizeType j = 0; j < n; j++)
            {
                A(iglobal + j, jglobal + j) = 4.00; // for now this operation generates a stupid compilation error with gcc
                if(j > 0)
                    A(iglobal + j, jglobal + j - 1) = -1.00;
                if(j < n-1)
                    A(iglobal + j, jglobal + j + 1) = -1.00;
            }

            // lower blocks
            if(i > 0)
            {
                iglobal = i * n;
                jglobal = (i-1) * n;
                for(SizeType j = 0; j < n; j++)
                    A(iglobal + j, jglobal + j) = -1.00;
            }

            // upper blocks
            if(i < m-1)
            {
                iglobal = i * n;
                jglobal = (i+1) * n;
                for(SizeType j = 0; j < n; j++)
                    A(iglobal + j, jglobal + j) = -1.00;
            }
        }

        SparseMatrixPointerType pA = SparseMatrixPointerType(new SparseMatrixType(A)); // error compile with gcc

        return pA;
    }

    /*@} */
    /**@name Operations */
    /*@{ */

    //////////////////////////////////////////////////////////////////
    //      STRENGTH OF CONNECTION
    //////////////////////////////////////////////////////////////////
    static void ClassicalStrengthOfConnection(SparseMatrixType& C, const SparseMatrixType& A, const double theta)
    {
        SizeType m = TSparseSpaceType::Size1(A);
        SizeType n = TSparseSpaceType::Size2(A);
//        SizeType nnz = A.value_data().size();
        SizeType nnz = A.filled2();

        IndexContainerType Cp(m + 1);
        IndexContainerType Ci(nnz, 0);
        ValueContainerType Cv(nnz, 0.00);

        classical_strength_of_connection(static_cast<IndexType>(m), theta, 
                A.index1_data(), A. index2_data(), A.value_data(), Cp, Ci, Cv);
        
//        boost::timer::cpu_timer filling_time;

//        Timer::Start("Fill");
        
        //SparseMatrixType C(CreateCoordinateMatrixFromTriplet(Cp, Ci, Cv));
        
        FillCompressedMatrixFromTriplet(C, Cp, Ci, Cv); // this was faster than filling from coordinate matrix
        
//        Timer::Stop("Fill");
        
//        boost::timer::cpu_times elapsed_times(filling_time.elapsed());
//        boost::timer::nanosecond_type elapsed(elapsed_times.system + elapsed_times.user);
//        std::cout << "Filling Time : " << elapsed << std::endl;
        
    }

    static void SymmetricStrengthOfConnection(SparseMatrixType& C, const SparseMatrixType& A, const double theta)
    {
        SizeType m = TSparseSpaceType::Size1(A);
        SizeType n = TSparseSpaceType::Size2(A);
//        SizeType nnz = A.value_data().size();
        SizeType nnz = A.filled2();

        IndexContainerType Cp(m + 1);
        IndexContainerType Ci(nnz, 0);
        ValueContainerType Cv(nnz, 0.00);

        symmetric_strength_of_connection(static_cast<IndexType>(m), theta, 
                A.index1_data(), A. index2_data(), A.value_data(), Cp, Ci, Cv);
        
//        boost::timer::cpu_timer filling_time;

//        Timer::Start("Fill");
        
        //SparseMatrixType C(CreateCoordinateMatrixFromTriplet(Cp, Ci, Cv));
        
        FillCompressedMatrixFromTriplet(C, Cp, Ci, Cv); // this was faster than filling from coordinate matrix
        
//        Timer::Stop("Fill");
        
//        boost::timer::cpu_times elapsed_times(filling_time.elapsed());
//        boost::timer::nanosecond_type elapsed(elapsed_times.system + elapsed_times.user);
//        std::cout << "Filling Time : " << elapsed << std::endl;
        
    }


    //////////////////////////////////////////////////////////////////
    //      INTERPOLATION METHOD
    //////////////////////////////////////////////////////////////////
    static void DirectInterpolation(SparseMatrixType& P, const SparseMatrixType& A, 
                const SparseMatrixType& C, const IndexVectorType& splitting)
    {
        SizeType m = TSparseSpaceType::Size1(A);

        IndexContainerType Pp(m + 1);
        
//        KRATOS_WATCH("before pass1");
//        std::cout << "splitting:" << std::endl;
//        std::copy(splitting.begin(), splitting.end(), std::ostream_iterator<IndexType>(std::cout, " "));
//        std::cout << std::endl;

        rs_direct_interpolation_pass1(static_cast<IndexType>(m), C.index1_data(), C.index2_data(), splitting,  Pp);

        SizeType nnz = Pp[m];
        IndexContainerType Pi(nnz, 0);
        ValueContainerType Pv(nnz, 0.00);

//        KRATOS_WATCH(nnz);
//        std::cout << "Pp:" << std::endl;
//        std::copy(Pp.begin(), Pp.end(), std::ostream_iterator<IndexType>(std::cout, " "));
//        std::cout << std::endl;
//        KRATOS_WATCH("before pass2");

        rs_direct_interpolation_pass2(static_cast<IndexType>(m), A.index1_data(), A.index2_data(), A.value_data(),
                C.index1_data(), C.index2_data(), C.value_data(), splitting,
                Pp, Pi, Pv);

//        KRATOS_WATCH_INT_LIST(Pp);
//        KRATOS_WATCH_INT_LIST(Pi);
//        KRATOS_WATCH_INT_LIST(Pv);
        
//        KRATOS_WATCH(Pp.size());
//        KRATOS_WATCH(Pi.size());
//        std::cout << "Pp:" << std::endl;
//        std::copy(Pp.begin(), Pp.end(), std::ostream_iterator<IndexType>(std::cout, " "));
//        std::cout << std::endl;
//        std::cout << "Pi:" << std::endl;
//        std::copy(Pi.begin(), Pi.end(), std::ostream_iterator<IndexType>(std::cout, " "));
//        std::cout << std::endl;
        
        SizeType n = *std::max_element(Pi.begin(), Pi.end());
//        KRATOS_WATCH(n+1);

        
//        SparseMatrixType P(m, n+1, nnz);
        if( (TSparseSpaceType::Size1(P) != m) || (TSparseSpaceType::Size2(P) != n+1) )
            TSparseSpaceType::Resize(P, m, n+1);

//        KRATOS_WATCH("before fill");
//        KRATOS_WATCH(TSparseSpaceType::Size1(P));
//        KRATOS_WATCH(TSparseSpaceType::Size2(P));
        
        FillCompressedMatrixFromTriplet(P, Pp, Pi, Pv); // this was faster than filling from coordinate matrix
        
//        KRATOS_WATCH("after fill");

//        KRATOS_WATCH(TSparseSpaceType::Size1(P));
//        KRATOS_WATCH(TSparseSpaceType::Size2(P));
//        
//        std::cout << "P.index1_data():" << std::endl;
//        std::copy(P.index1_data().begin(), P.index1_data().end(), std::ostream_iterator<IndexType>(std::cout, " "));
//        std::cout << std::endl;
//        std::cout << "P.index2_data():" << std::endl;
//        std::copy(P.index2_data().begin(), P.index2_data().end(), std::ostream_iterator<IndexType>(std::cout, " "));
//        std::cout << std::endl;
    }


    //////////////////////////////////////////////////////////////////
    //      SPLITTING METHOD
    //////////////////////////////////////////////////////////////////
    static void RS(IndexVectorType& splitting, const SparseMatrixType& S)
    {
        SizeType m = TSparseSpaceType::Size1(S);
        SizeType n = TSparseSpaceType::Size2(S);
//        SizeType nnz = S.value_data().size();
        SizeType nnz = S.filled2();
        
        IndexContainerType Tp(n + 1, 0);
        IndexContainerType Ti(nnz, 0);
        
        TripletTranspose(m, n, S.index1_data(), S.index2_data(), Tp, Ti);
        
//        KRATOS_WATCH("TripletTranspose completed");
        
//        KRATOS_WATCH_INT_LIST(S.index1_data());
//        KRATOS_WATCH_INT_LIST(S.index2_data());
//        KRATOS_WATCH_INT_LIST(Tp);
//        KRATOS_WATCH_INT_LIST(Ti);
        
        rs_cf_splitting(static_cast<IndexType>(m), S.index1_data(), S.index2_data(), 
            Tp, Ti, splitting);
            
//        KRATOS_WATCH_INT_LIST(splitting);

    }

    static void PMIS(IndexVectorType& splitting, const SparseMatrixType& S)
    {
        SizeType m = TSparseSpaceType::Size1(S);
        SizeType n = TSparseSpaceType::Size2(S);
//        SizeType nnz = S.value_data().size();
        SizeType nnz = S.filled2();
        SparseMatrixType G(m, n, nnz);
        ValueContainerType weights(m);
        preprocess_ruge_stuben_splitting(S, COLORING_NONE, G, weights);
        MIS(splitting, G, weights);
    }
    
    static void PMISc(IndexVectorType& splitting, const SparseMatrixType& S, const std::string& coloring_method)
    {
        SizeType m = TSparseSpaceType::Size1(S);
        SizeType n = TSparseSpaceType::Size2(S);
//        SizeType nnz = S.value_data().size();
        SizeType nnz = S.filled2();
        SparseMatrixType G(m, n, nnz);
        ValueContainerType weights(m);
                
        if(coloring_method.compare("MIS"))
        {
            preprocess_ruge_stuben_splitting(S, COLORING_MIS, G, weights);
        }
        else if(coloring_method.compare("JP"))
        {
            preprocess_ruge_stuben_splitting(S, COLORING_JP, G, weights);
        }
        else if(coloring_method.compare("LDF"))
        {
            preprocess_ruge_stuben_splitting(S, COLORING_LDF, G, weights);
        }
        else
            KRATOS_THROW_ERROR(std::logic_error, "unrecognized coloring method", "");
                
        MIS(splitting, G, weights);
        
    }
    
    static void CLJP(IndexVectorType& splitting, const SparseMatrixType& S)
    {
        SizeType m = TSparseSpaceType::Size1(S);
        SizeType n = TSparseSpaceType::Size2(S);
//        SizeType nnz = S.value_data().size();
        SizeType nnz = S.filled2();

        IndexContainerType Tp(n + 1, 0);
        IndexContainerType Ti(nnz, 0);

        TripletTranspose(m, n, S.index1_data(), S.index2_data(), Tp, Ti);

        cljp_naive_splitting(m, S.index1_data(), S.index2_data(), Tp, Ti, splitting, 0);
        
    }
    
    static void CLJPc(IndexVectorType& splitting, const SparseMatrixType& S)
    {
        SizeType m = TSparseSpaceType::Size1(S);
        SizeType n = TSparseSpaceType::Size2(S);
//        SizeType nnz = S.value_data().size();
        SizeType nnz = S.filled2();

        IndexContainerType Tp(n + 1, 0);
        IndexContainerType Ti(nnz, 0);

        TripletTranspose(m, n, S.index1_data(), S.index2_data(), Tp, Ti);

        cljp_naive_splitting(m, S.index1_data(), S.index2_data(), Tp, Ti, splitting, 1);
        
    }
    
    
    
    

    //////////////////////////////////////////////////////////////////
    //      FAST LINEAR ALGEBRA
    //////////////////////////////////////////////////////////////////
    // create the transpose of a sparse matrix rB = rA^T
    static void Transpose(SparseMatrixType& rA, SparseMatrixType& rB)
    {
//        return TSparseSpaceType::Transpose(rA); //very slow for big matrix
        int m = TSparseSpaceType::Size1(rA);
        int n = TSparseSpaceType::Size2(rA);
//        int nnz = rA.value_data().size();
        int nnz = rA.filled2();
        
        //checking
        if( TSparseSpaceType::Size1(rB) != n || TSparseSpaceType::Size2(rB) != m)
            KRATOS_THROW_ERROR(std::logic_error, "The matrix dimension is not compatible", "");
        
        IndexContainerType Tp(n+1, 0);
        IndexContainerType Ti(nnz, 0);
        ValueContainerType Tx(nnz, 0.00);
        
//        KRATOS_WATCH("at Transpose: before triplet transpose");
        TripletTranspose(m, n, rA.index1_data(), rA.index2_data(), rA.value_data(), Tp, Ti, Tx);
        
//        csr_tocsc(m, n, rA.index1_data(), rA.index2_data(), rA.value_data(), Tp, Ti, Tx);
        
//        KRATOS_WATCH("at Transpose: before fill");
        FillCompressedMatrixFromTriplet(rB, Tp, Ti, Tx);
        
    }

    // perform multiplication of 2 sparse matrices by SMMP algorithm (adapted from scipy), rC = rA*rB
    static void Mult(SparseMatrixType& rA, SparseMatrixType& rB, SparseMatrixType& rC)
    {
        int m = TSparseSpaceType::Size1(rA);
        int n = TSparseSpaceType::Size2(rA);
        int k = TSparseSpaceType::Size2(rB);
        
        if(n != TSparseSpaceType::Size1(rB))
            KRATOS_THROW_ERROR(std::logic_error, "Matrix dimension is not compatible", "");
            
        IndexContainerType Cp(m+1, 0);
        
        csr_matmat_pass1(m, k, rA.index1_data(), rA.index2_data(), 
                            rB.index1_data(), rB.index2_data(), Cp);
        
        SizeType nnz = Cp[m]; // this nnz is the symbolic entries. Some symbolic entries may be zero and will be detected in pass 2
        
        IndexContainerType Ci(nnz, 0);
        ValueContainerType Cx(nnz, 0.00);
        
        csr_matmat_pass2(m, k, rA.index1_data(), rA.index2_data(), rA.value_data(),
                            rB.index1_data(), rB.index2_data(), rB.value_data(), Cp, Ci, Cx);

        FillCompressedMatrixFromTriplet(rC, Cp, Ci, Cx); // this was faster than filling from coordinate matrix
        
    }
    
    
    //////////////////////////////////////////////////////////////////
    //      HELPERS (public domain)
    //////////////////////////////////////////////////////////////////
    
    static void AddToDiagonal(SparseMatrixType& rA, const double v)
    {
        SizeType n = std::min(TSparseSpaceType::Size1(rA), TSparseSpaceType::Size2(rA));
        for(IndexType i = 0; i < n; i++)
            rA(i, i) += v;
    }
    
    /*@} */
    /**@name Access */
    /*@{ */


    /*@} */
    /**@name Inquiry */
    /*@{ */


    /*@} */
    /**@name Friends */
    /*@{ */


    /*@} */

private:
    /**@name Static Member Variables */
    /*@{ */


    /*@} */
    /**@name Member Variables */
    /*@{ */


    /*@} */
    /**@name Private Operators*/
    /*@{ */


    /*@} */
    /**@name Private Operations*/
    /*@{ */

//    LocalSparseMatrixType CreateCoordinateMatrixFromTriplet(IndexContainerType& Cp, IndexContainerType& Ci, ValueContainerType& Cv)
//    {
//        SizeType m = Cp.size() - 1;
//        SizeType nnz = Cv.size();
//        
//        LocalSparseMatrixType C(m, m, nnz);
//        
//        for(IndexType i = 0; i < m; i++)
//        {
//            SizeType local_nnz = Cp[i+1] - Cp[i];
//            for(IndexType j = 0; j < local_nnz; j++)
//            {
//                IndexType idx = Cp[i] + j;
//                C.insert_element(i, Ci[idx], Cv[idx]);
//            }
//        }
//        
//        return C;
//    }
    
    
    //////////////////////////////////////////////////////////////////
    //      HELPERS
    //////////////////////////////////////////////////////////////////
    //input: Cp, Ci, Cv
    //output: C
    static void FillCompressedMatrixFromTriplet(SparseMatrixType& C, IndexContainerType& Cp, IndexContainerType& Ci, ValueContainerType& Cv)
    {
        SizeType m = Cp.size() - 1;
//        SizeType m = TSparseSpaceType::Size1(C);
        
        for(IndexType i = 0; i < m; i++)
        {
            for(IndexType j = Cp[i]; j < Cp[i+1]; j++)
            {
                C(i, Ci[j]) = Cv[j];
//                C.push_back(i, Ci[j], Cv[j]); // note that this only works if the column indices are sorted. It is not the case if the matrix is output from SMMP method
            }
        }
        
        C.complete_index1_data(); // to ensure zero rows are handled correctly
        
    }

    //input: Cp, Ci, v (value to fill)
    //output: C
    static void FillCompressedMatrixFromTriplet(SparseMatrixType& C, IndexContainerType& Cp, IndexContainerType& Ci, double v)
    {
        SizeType m = Cp.size() - 1;
        
        for(IndexType i = 0; i < m; i++)
        {
            SizeType local_nnz = Cp[i+1] - Cp[i];
            for(IndexType j = 0; j < local_nnz; j++)
            {
                IndexType idx = Cp[i] + j;
                C(i, Ci[idx]) = v;
//                C.push_back(i, Ci[idx], v); // note that this only works if the column indices are sorted. It is not the case if the matrix is output from SMMP method
            }
        }
        
        C.complete_index1_data(); // to ensure zero rows are handled correctly
        
    }
    
    // transpose a row pointer and column index to column poiner and row index
    //input Sp, Si
    //output Tp, Ti
    static void TripletTranspose(const IndexType m, const IndexType n, 
            const IndexContainerType& Sp, const IndexContainerType& Si, 
            IndexContainerType& Tp, IndexContainerType& Ti)
    {
        SizeType nnz = Si.size();
        
        std::vector<std::vector<IndexType> > colMap;
        
        colMap.resize(n);
        SizeType estimated_nnz = nnz/n;
        for(IndexType i = 0; i < n; i++)
            colMap[i].reserve(estimated_nnz);
        
        for(IndexType i = 0; i < m; i++)
        {
            SizeType local_nnz = Sp[i+1] - Sp[i];
            for(IndexType j = 0; j < local_nnz; j++)
            {
                IndexType idx = Sp[i] + j;
                colMap[Si[idx]].push_back(i);
            }
        }
        
//        for(IndexType i = 0; i < m; i++)
//        {
////            KRATOS_WATCH_LIST(colMap[i]);
//            std::copy(colMap[i].begin(), colMap[i].end(), std::ostream_iterator<IndexType>(std::cout, " "));
//            std::cout << std::endl;
//        }
            
//        std::cout << "colMap initialization completed" << std::endl;
        
//        for(IndexType i = 0; i < m; i++)
//            std::sort(colMap[i].begin(), colMap[i].end());
//            
//        for(IndexType i = 0; i < m; i++)
//        {
//            std::copy(colMap[i].begin(), colMap[i].end(), std::ostream_iterator<IndexType>(std::cout, " "));
//            std::cout << std::endl;
//        }

//        std::cout << "colMap sort completed" << std::endl;
            
        Tp[0] = 0;
        for(IndexType i = 0; i < n; i++)
        {
            Tp[i+1] = Tp[i] + colMap[i].size();
            for(IndexType j = 0; j < colMap[i].size(); j++)
            {
                Ti[Tp[i] + j] = colMap[i][j];
            }
        }
        
//        std::copy(Tp.begin(), Tp.end(), std::ostream_iterator<IndexType>(std::cout, " "));
//        std::cout << std::endl;
//        std::copy(Ti.begin(), Ti.end(), std::ostream_iterator<IndexType>(std::cout, " "));
//        std::cout << std::endl;
//        
//        std::cout << "insert completed" << std::endl;
    }
    
    // transpose a csr sparse matrix
    // input:  Sp, Si, Sx
    // output: Tp, Ti, Tx
    static void TripletTranspose(const IndexType m, const IndexType n, 
            const IndexContainerType& Sp, const IndexContainerType& Si, const ValueContainerType& Sx, 
            IndexContainerType& Tp, IndexContainerType& Ti, ValueContainerType& Tx)
    {
        SizeType nnz = Si.size();
        
        std::vector<std::vector<IndexType> > colMap;
        std::vector<std::vector<ValueType> > valMap;
        
        colMap.resize(n);
        valMap.resize(n);
        SizeType estimated_nnz = nnz/n;
        
//        KRATOS_WATCH("at TripletTranspose: before reserve");
        
        for(IndexType i = 0; i < n; i++)
        {
            colMap[i].reserve(estimated_nnz);
            valMap[i].reserve(estimated_nnz);
        }
        
//        KRATOS_WATCH("at TripletTranspose: before push_back");
//        KRATOS_WATCH(m);
//        KRATOS_WATCH(n);
//        KRATOS_WATCH(Sp.size());
//        KRATOS_WATCH(Si.size());
//        
//        std::cout << "Sp: " << std::endl;
//        std::copy(Sp.begin(), Sp.end(), std::ostream_iterator<IndexType>(std::cout, " "));
//        std::cout << std::endl;
//        std::cout << "Si: " << std::endl;
//        std::copy(Si.begin(), Si.end(), std::ostream_iterator<IndexType>(std::cout, " "));
//        std::cout << std::endl;


        for(IndexType i = 0; i < m; i++)
        {
            SizeType local_nnz = Sp[i+1] - Sp[i];
            for(IndexType j = 0; j < local_nnz; j++)
            {
                IndexType idx = Sp[i] + j;
                colMap[Si[idx]].push_back(i); //Si[idx] is col index
                valMap[Si[idx]].push_back(Sx[idx]);
            }
        }
        
//        KRATOS_WATCH("at TripletTranspose: before fill");
        
        Tp[0] = 0;
        for(IndexType i = 0; i < n; i++)
        {
            Tp[i+1] = Tp[i] + colMap[i].size();
            for(IndexType j = 0; j < colMap[i].size(); j++)
            {
                Ti[Tp[i] + j] = colMap[i][j]; //colMap[i][j] is row index (i is col index)
                Tx[Tp[i] + j] = valMap[i][j];
            }
        }
        
    }
    
    // fill the row pointer and column index of matrix G = S+S^T
    // input : Sp, Si
    // output: Gp, Gi, nnz (number of entries of the sum)
    static void TripletSumSandST(const IndexType m, const IndexType n, 
            const IndexContainerType& Sp, const IndexContainerType& Si,
            IndexContainerType& Gp, IndexContainerType &Gi, SizeType &nnz)
    {
        assert(m==n); //must be square matrix
        
        nnz = Si.size();
        
        std::vector<std::set<IndexType> > colMap;
        typedef std::set<IndexType>::iterator set_iterator;
        
        colMap.resize(n);
        
        for(IndexType i = 0; i < m; i++)
        {
            SizeType local_nnz = Sp[i+1] - Sp[i];
            for(IndexType j = 0; j < local_nnz; j++)
            {
                IndexType idx = Sp[i] + j;
                // i is row, Si[idx] is col
                colMap[Si[idx]].insert(i);
                colMap[i].insert(Si[idx]);
            }
        }
        
        Gp[0] = 0;
        nnz = 0;
        for(IndexType i = 0; i < n; i++)
        {
            Gp[i+1] = Gp[i] + colMap[i].size();
            IndexType cnt = 0;
            nnz += colMap[i].size();
            for(set_iterator j = colMap[i].begin(); j != colMap[i].end(); ++j)
            {
                Gi[Gp[i] + cnt] = *j;
                cnt++;
            }
        }
    }
    
    
    // generate an array of random number in [0,1]
    static ValueContainerType GenerateRandomArray(IndexType n)
    {
        typedef typename boost::random::mt19937                            ENG;    // Mersenne Twister
        typedef typename boost::random::uniform_real_distribution<double> DIST;   // Uniform Distribution
        typedef typename boost::random::variate_generator<ENG,DIST>        GEN;    // Variate generator
        
        ENG  eng;
        DIST dist(0.00, 1.00);
        GEN  gen(eng,dist);
        
        ValueContainerType v(n);
        
        for(IndexType i = 0; i < n; i++)
            v[i] = gen();
        
        return v;
    }
    
    
    //////////////////////////////////////////////////////////////////
    //      HELPER FUNCTIONS FOR SPLITTING
    //////////////////////////////////////////////////////////////////
    // G must be of size m x m and weights is of size m
    static void preprocess_ruge_stuben_splitting(const SparseMatrixType &S, const IndexType coloring_method, //input
            SparseMatrixType& G, ValueContainerType& weights)                                     //output
    {
    
        SizeType m = TSparseSpaceType::Size1(S);
//        SizeType nnz = S.value_data().size();
        SizeType nnz = S.filled2();
        
        IndexContainerType Gp(m+1);
        IndexContainerType Gi(2*nnz);
        
        const IndexContainerType& Sp = S.index1_data();
        const IndexContainerType& Si = S.index2_data();
        
        TripletSumSandST(m, m, Sp, Si, Gp, Gi, nnz);
        
//        G.resize(m, m, false);
        FillCompressedMatrixFromTriplet(G, Gp, Gi, 1.00);
        
//        weights.resize(m, false);
        
        // sum along the column of S
        for(IndexType i = 0; i < m; i++)
        {
            SizeType local_nnz = Sp[i+1] - Sp[i];
            for(IndexType j = 0; j < local_nnz; j++)
            {
                IndexType idx = Sp[i] + j;
                // i is row, Si[idx] is col
                weights[Si[idx]] += 1.00;
            }
        }
        
        
        if(coloring_method == COLORING_NONE)
        {
            // generate random number
//            std::uniform_real_distribution<double> unif(0.00, 1.00);
//            std::default_random_engine re;
//            
//            for(IndexType i = 0; i < m; i++)
//                weights[i] += unif(re);

            ValueContainerType ra = GenerateRandomArray(m);
            for(IndexType i = 0; i < m; i++)
                weights[i] += ra[i];
                
        }
        else if(coloring_method == COLORING_MIS || coloring_method == COLORING_JP 
                    || coloring_method == COLORING_LDF)
        {
            IntegerContainerType coloring(m);
            
            
            if(coloring_method == COLORING_MIS)
                vertex_coloring_mis(m, Gp, Gi, coloring);
            else if(coloring_method == COLORING_JP)
            {
                ValueContainerType tmp = GenerateRandomArray(m);
                vertex_coloring_jones_plassmann(m, Gp, Gi, coloring, tmp);
            }
            else if(coloring_method == COLORING_LDF)
            {
                ValueContainerType tmp = GenerateRandomArray(m);
                vertex_coloring_LDF(m, Gp, Gi, coloring, tmp);
            }
            
            
            IndexType num_colors = 1 + *std::max(coloring.begin(), coloring.end());

//            std::uniform_real_distribution<double> unif(0.00, 1.00);
//            std::default_random_engine re;
//            
//            for(IndexType i = 0; i < m; i++)
//                weights[i] += (unif(re) + coloring[i])/num_colors;

            ValueContainerType ra = GenerateRandomArray(m);
            for(IndexType i = 0; i < m; i++)
                weights[i] += (ra[i] + coloring[i])/num_colors;
        }
        else
            KRATOS_THROW_ERROR(std::logic_error, "undetected coloring method", "");
        
    }
    
    
    //////////////////////////////////////////////////////////////////
    //      GRAPH METHOD
    //////////////////////////////////////////////////////////////////
    
    static void MIS(IndexVectorType& mis, const SparseMatrixType& G, const ValueContainerType& weights, int maxiter = -1)
    {
        int m = TSparseSpaceType::Size1(G);
        
        std::fill(mis.begin(), mis.end(), -1);
        
        if(maxiter < 0)
        {
            maximal_independent_set_parallel(m, G.index1_data(), G.index2_data(), -1, 1, 0, mis, weights);
        }
        else
        {
            maximal_independent_set_parallel(m, G.index1_data(), G.index2_data(), -1, 1, 0, mis, weights, maxiter);
        }
    }
    
    
    /*@} */
    /**@name Private  Acces */
    /*@{ */


    /*@} */
    /**@name Private Inquiry */
    /*@{ */


    /*@} */
    /**@name Un accessible methods */
    /*@{ */




    /*@} */

}; /* Class ClassName */

/*@} */

/**@name Type Definitions */
/*@{ */


/*@} */

} /* namespace Kratos.*/

#endif /* AMG_UTILS  defined */

