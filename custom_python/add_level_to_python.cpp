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
#include "custom_utilities/mg_level.h"
#include "custom_utilities/matrix_based_mg_level.h"
#include "custom_utilities/mg_projector.h"
#include "custom_utilities/matrix_based_mg_projector.h"
#include "custom_utilities/structured_mg_prolongator.h"
#include "custom_utilities/structured_mg_restrictor.h"
#include "custom_python/add_level_to_python.h"

namespace Kratos
{

namespace Python
{

typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;

typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

typedef MGLevel<SparseSpaceType, LocalSpaceType> MGLevelType;

typedef MatrixBasedMGLevel<SparseSpaceType, LocalSpaceType> MatrixBasedMGLevelType;

typedef LinearSolver<SparseSpaceType, LocalSpaceType> LinearSolverType;

typedef MGProjector<SparseSpaceType> MGProjectorType;

typedef MatrixBasedMGProjector<SparseSpaceType> MatrixBasedMGProjectorType;

typedef StructuredMGProlongator<SparseSpaceType, 2> StructuredMGProlongator2DType;

typedef StructuredMGRestrictor<SparseSpaceType, 2> StructuredMGRestrictor2DType;

void MultigridSolversApp_AddLevelToPython()
{
    using namespace boost::python;

   //****************************************************************************************************
   // level definition
   //****************************************************************************************************

    class_<MGLevelType, MGLevelType::Pointer, boost::noncopyable>
    ( "MGLevel", init<const typename MGLevelType::IndexType&>())
    .def(self_ns::str(self))
    .def("SetPreSmoother", &MGLevelType::SetPreSmoother)
    .def("SetPostSmoother", &MGLevelType::SetPostSmoother)
    .def("SetRestrictionOperator", &MGLevelType::SetRestrictionOperator)
    .def("SetProlongationOperator", &MGLevelType::SetProlongationOperator)
    .def("ApplyPreSmoother", &MGLevelType::ApplyPreSmoother)
    .def("ApplyPostSmoother", &MGLevelType::ApplyPostSmoother)
    .def("ApplyRestriction", &MGLevelType::ApplyRestriction)
    .def("ApplyProlongation", &MGLevelType::ApplyProlongation)
    ;

    class_<MatrixBasedMGLevelType, MatrixBasedMGLevelType::Pointer, bases<MGLevelType>, boost::noncopyable>
    ( "MatrixBasedMGLevel", init<const typename MGLevelType::IndexType&>())
    .def(self_ns::str(self))
    .def("GetCoarseMatrix", &MatrixBasedMGLevelType::GetCoarseMatrix)
    .def("GetCoarseUpdateVector", &MatrixBasedMGLevelType::GetCoarseUpdateVector)
    .def("GetCoarseVector", &MatrixBasedMGLevelType::GetCoarseVector)
    ;

    //****************************************************************************************************
    // projection definition
    //****************************************************************************************************

    class_<MGProjectorType, MGProjectorType::Pointer, boost::noncopyable>
    ( "MGProjector", init<>())
    .def(self_ns::str(self))
    .def("Initialize", &MGProjectorType::Initialize)
    .def("Apply", &MGProjectorType::Apply)
    ;

    class_<MatrixBasedMGProjectorType, MatrixBasedMGProjectorType::Pointer, bases<MGProjectorType>, boost::noncopyable>
    ( "MatrixBasedMGProjector", init<>())
    .def("SetOperator", &MatrixBasedMGProjectorType::SetOperator)
    .def(self_ns::str(self))
    ;

    class_<StructuredMGProlongator2DType, StructuredMGProlongator2DType::Pointer, bases<MGProjectorType>, boost::noncopyable>
    ( "StructuredMGProlongator2D", init<ModelPart::Pointer, ModelPart::Pointer>())
    .def(self_ns::str(self))
    ;

    class_<StructuredMGRestrictor2DType, StructuredMGRestrictor2DType::Pointer, bases<MGProjectorType>, boost::noncopyable>
    ( "StructuredMGRestrictor2D", init<ModelPart::Pointer, ModelPart::Pointer>())
    .def(self_ns::str(self))
    ;

}

}  // namespace Python.

} // Namespace Kratos

