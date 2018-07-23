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
#include "custom_utilities/mesh_based_mg_projector.h"
#include "custom_utilities/structured_mesh_mg_projector.h"
#include "custom_utilities/structured_mesh_mg_interpolator.h"
#include "custom_utilities/structured_mesh_mg_transpose_interpolator.h"
#include "custom_utilities/structured_mesh_matrix_based_mg_interpolator.h"
#include "custom_utilities/structured_mesh_matrix_based_mg_transpose_interpolator.h"
#include "custom_python/add_level_to_python.h"

namespace Kratos
{

namespace Python
{

typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;

typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

typedef LinearSolver<SparseSpaceType, LocalSpaceType> LinearSolverType;

void MultigridSolversApp_AddLevelToPython()
{
    using namespace boost::python;

   //****************************************************************************************************
   // level definition
   //****************************************************************************************************

    typedef MGLevel<SparseSpaceType, LocalSpaceType> MGLevelType;
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

    typedef MatrixBasedMGLevel<SparseSpaceType, LocalSpaceType> MatrixBasedMGLevelType;
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

    typedef MGProjector<SparseSpaceType> MGProjectorType;
    class_<MGProjectorType, MGProjectorType::Pointer, boost::noncopyable>
    ("MGProjector", init<>())
    .def(self_ns::str(self))
    .def("Initialize", &MGProjectorType::Initialize)
    .def("Apply", &MGProjectorType::Apply)
    ;

    typedef MatrixBasedMGProjector<SparseSpaceType> MatrixBasedMGProjectorType;
    class_<MatrixBasedMGProjectorType, MatrixBasedMGProjectorType::Pointer, bases<MGProjectorType>, boost::noncopyable>
    ("MatrixBasedMGProjector", init<>())
    .def("SetOperator", &MatrixBasedMGProjectorType::SetOperator)
    .def(self_ns::str(self))
    ;

    typedef MeshBasedMGProjector<SparseSpaceType> MeshBasedMGProjectorType;
    class_<MeshBasedMGProjectorType, MeshBasedMGProjectorType::Pointer, bases<MGProjectorType>, boost::noncopyable>
    ("MeshBasedMGProjector", init<>())
    .def("SetBlockSize", &MeshBasedMGProjectorType::SetBlockSize)
    .def("SetCoarseModelPart", &MeshBasedMGProjectorType::SetCoarseModelPart)
    .def("SetFineModelPart", &MeshBasedMGProjectorType::SetFineModelPart)
    .def(self_ns::str(self))
    ;

    typedef StructuredMeshMGProjector<SparseSpaceType, 2> StructuredMeshMGProjector2DType;
    class_<StructuredMeshMGProjector2DType, StructuredMeshMGProjector2DType::Pointer, bases<MeshBasedMGProjectorType>, boost::noncopyable>
    ("StructuredMeshMGProjector2D", init<>())
    .def("SetFineMeshSize", &StructuredMeshMGProjector2DType::SetFineMeshSize)
    .def("SetCoarseMeshSize", &StructuredMeshMGProjector2DType::SetCoarseMeshSize)
    .def(self_ns::str(self))
    ;

    typedef StructuredMeshMGInterpolator<SparseSpaceType, 2> StructuredMeshMGInterpolator2DType;
    class_<StructuredMeshMGInterpolator2DType, StructuredMeshMGInterpolator2DType::Pointer, bases<StructuredMeshMGProjector2DType>, boost::noncopyable>
    ("StructuredMeshMGInterpolator2D", init<>())
    .def(self_ns::str(self))
    ;

    typedef StructuredMeshMGTransposeInterpolator<SparseSpaceType, 2> StructuredMeshMGTransposeInterpolator2DType;
    class_<StructuredMeshMGTransposeInterpolator2DType, StructuredMeshMGTransposeInterpolator2DType::Pointer, bases<StructuredMeshMGProjector2DType>, boost::noncopyable>
    ("StructuredMeshMGTransposeInterpolator2D", init<>())
    .def(self_ns::str(self))
    ;

    typedef StructuredMeshMatrixBasedMGInterpolator<SparseSpaceType, 2> StructuredMeshMatrixBasedMGInterpolator2DType;
    class_<StructuredMeshMatrixBasedMGInterpolator2DType, StructuredMeshMatrixBasedMGInterpolator2DType::Pointer, bases<StructuredMeshMGProjector2DType, MatrixBasedMGProjectorType>, boost::noncopyable>
    ("StructuredMeshMatrixBasedMGInterpolator2D", init<>())
    .def(self_ns::str(self))
    ;

    typedef StructuredMeshMatrixBasedMGTransposeInterpolator<SparseSpaceType, 2> StructuredMeshMatrixBasedMGTransposeInterpolator2DType;
    class_<StructuredMeshMatrixBasedMGTransposeInterpolator2DType, StructuredMeshMatrixBasedMGTransposeInterpolator2DType::Pointer, bases<StructuredMeshMGProjector2DType, MatrixBasedMGProjectorType>, boost::noncopyable>
    ("StructuredMeshMatrixBasedMGTransposeInterpolator2D", init<>())
    .def(self_ns::str(self))
    ;
}

}  // namespace Python.

} // Namespace Kratos

