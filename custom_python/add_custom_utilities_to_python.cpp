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

//
//   Project Name:        Kratos
//   Last modified by:    $Author: hbui $
//   Date:                $Date: Jan 9, 2013 $
//   Revision:            $Revision: 1.1 $
//
//


// System includes
#include <string>

// External includes
#include <boost/python.hpp>

// Project includes
#include "includes/define.h"
#include "spaces/ublas_space.h"
#include "custom_utilities/amg_utils.h"
#include "custom_utilities/amg_level.h"
#include "custom_utilities/multilevel_solver_factory.h"
#include "custom_python/add_custom_utilities_to_python.h"

namespace Kratos
{

namespace Python
{

typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;

typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

typedef AMGUtils<SparseSpaceType> AMGUtilsType;

typedef typename AMGUtilsType::SparseMatrixType SparseMatrixType;

typedef typename AMGUtilsType::SparseMatrixPointerType SparseMatrixPointerType;

typedef typename AMGUtilsType::IndexVectorType IndexVectorType;

//typedef typename AMGUtilsType::IndexVectorPointerType IndexVectorPointerType;

typedef AMGUtilsType::SizeType SizeType;

typedef AMGLevel<SparseSpaceType, LocalSpaceType> AMGLevelType;

typedef LinearSolver<SparseSpaceType, LocalSpaceType> LinearSolverType;

//SparseMatrixPointerType Poisson(AMGUtilsType& dummy, const SizeType m, const SizeType n)
//{
//    return dummy.Poisson(m, n);
//}

SparseMatrixPointerType ClassicalStrengthOfConnection(AMGUtilsType& dummy, const SparseMatrixType& A, const double theta)
{
    SizeType m = SparseSpaceType::Size1(A);
    SparseMatrixPointerType C = SparseMatrixPointerType(new SparseMatrixType(m, m));
    dummy.ClassicalStrengthOfConnection(*C, A, theta);
    return C;
}

SparseMatrixPointerType SymmetricStrengthOfConnection(AMGUtilsType& dummy, const SparseMatrixType& A, const double theta)
{
    SizeType m = SparseSpaceType::Size1(A);
    SparseMatrixPointerType C = SparseMatrixPointerType(new SparseMatrixType(m, m));
    dummy.SymmetricStrengthOfConnection(*C, A, theta);
    return C;
}

SparseMatrixPointerType DirectInterpolation(AMGUtilsType& dummy, const SparseMatrixType& A,
                const SparseMatrixType& C, const IndexVectorType& splitting)
{
    SparseMatrixPointerType P = SparseMatrixPointerType(new SparseMatrixType(0, 0));
    dummy.DirectInterpolation(*P, A, C, splitting);
    return P;
}

IndexVectorType RS(AMGUtilsType& dummy, const SparseMatrixType& S)
{
    SizeType m = SparseSpaceType::Size1(S);
    IndexVectorType splitting(m);
    dummy.RS(splitting, S);
    return splitting;
}

IndexVectorType PMIS(AMGUtilsType& dummy, const SparseMatrixType& S)
{
    SizeType m = SparseSpaceType::Size1(S);
    IndexVectorType splitting(m);
    dummy.PMIS(splitting, S);
    return splitting;
}

IndexVectorType PMISc(AMGUtilsType& dummy, const SparseMatrixType& S, const std::string& coloring_method)
{
    SizeType m = SparseSpaceType::Size1(S);
    IndexVectorType splitting(m);
    dummy.PMISc(splitting, S, coloring_method);
    return splitting;
}

IndexVectorType CLJP(AMGUtilsType& dummy, const SparseMatrixType& S)
{
    SizeType m = SparseSpaceType::Size1(S);
    IndexVectorType splitting(m);
    dummy.CLJP(splitting, S);
    return splitting;
}

IndexVectorType CLJPc(AMGUtilsType& dummy, const SparseMatrixType& S)
{
    SizeType m = SparseSpaceType::Size1(S);
    IndexVectorType splitting(m);
    dummy.CLJPc(splitting, S);
    return splitting;
}

SparseMatrixPointerType Transpose(AMGUtilsType& dummy, SparseMatrixType& rA)
{
    SizeType m = SparseSpaceType::Size1(rA);
    SizeType n = SparseSpaceType::Size2(rA);
    SparseMatrixPointerType T = SparseMatrixPointerType(new SparseMatrixType(m, n));
    dummy.Transpose(rA, *T);
    return T;
}

SparseMatrixPointerType Mult(AMGUtilsType& dummy, SparseMatrixType& rA, SparseMatrixType& rB)
{
    SizeType m = SparseSpaceType::Size1(rA);
    SizeType n = SparseSpaceType::Size2(rB);
    SparseMatrixPointerType C = SparseMatrixPointerType(new SparseMatrixType(m, n));
    dummy.Mult(rA, rB, *C);
    return C;
}

void MultigridSolversApp_AddUtilitiesToPython()
{
    using namespace boost::python;

    class_<AMGUtilsType, AMGUtilsType::Pointer, boost::noncopyable>("AMGUtils", init<>())
//    .def("deref", &AMGUtilsType::deref)
//    .def("Poisson", Poisson)
    .def("ClassicalStrengthOfConnection", ClassicalStrengthOfConnection)
    .def("SymmetricStrengthOfConnection", SymmetricStrengthOfConnection)
    .def("RS", RS)
    .def("PMIS", PMIS)
    .def("PMISc", PMISc)
    .def("CLJP", CLJP)
    .def("CLJPc", CLJPc)
    .def("DirectInterpolation", DirectInterpolation)
    .def("Transpose", Transpose)
    .def("Mult", Mult)
    ;

//    //****************************************************************************************************
//    //multilevel factory and amglevel
//    //****************************************************************************************************

    class_<AMGLevelType, AMGLevelType::Pointer, boost::noncopyable >( "AMGLevel", init<LinearSolverType::Pointer, LinearSolverType::Pointer >())
    .def(self_ns::str(self))
    .def("ApplyPreSmoother", &AMGLevelType::ApplyPreSmoother)
    .def("ApplyPostSmoother", &AMGLevelType::ApplyPostSmoother)
    .def("ApplyRestriction", &AMGLevelType::ApplyRestriction)
    .def("ApplyProlongation", &AMGLevelType::ApplyProlongation)
    .def("SetPreSmoother", &AMGLevelType::SetPreSmoother)
    .def("SetPostSmoother", &AMGLevelType::SetPostSmoother)
    .def("SetRestrictionOperator", &AMGLevelType::SetRestrictionOperator)
    .def("SetProlongationOperator", &AMGLevelType::SetProlongationOperator)
    .def("SetCoarsenMatrix", &AMGLevelType::SetCoarsenMatrix)
//    .def("ComputeCoarsenMatrix", &AMGLevelType::ComputeCoarsenMatrix)
    .def("GetCoarsenMatrix", &AMGLevelType::GetCoarsenMatrix)
    ;
}

}  // namespace Python.

} // Namespace Kratos

