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
#include <boost/foreach.hpp>
#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>
#include <boost/python/operators.hpp>

// Project includes
#include "includes/define.h"
#include "spaces/ublas_space.h"
#include "custom_utilities/mg_level.h"
#include "custom_utilities/matrix_based_mg_level.h"
#include "custom_utilities/mg_projector.h"
#include "custom_utilities/mg_null_projector.h"
#include "custom_utilities/mg_transpose_projector.h"
#include "custom_utilities/matrix_based_mg_projector.h"
#include "custom_utilities/mesh_based_mg_projector.h"
#include "custom_utilities/index_based_mg_projector.h"
#include "custom_utilities/structured_mesh_mg_projector.h"
#include "custom_utilities/structured_mesh_mg_interpolator.h" // prolongation
#include "custom_utilities/structured_mesh_mg_transpose_interpolator.h" // restriction
#include "custom_utilities/structured_mesh_matrix_based_mg_interpolator.h" // prolongation
#include "custom_utilities/structured_mesh_matrix_based_mg_transpose_interpolator.h" // restriction
#include "custom_utilities/structured_mesh_mg_injector.h" // restriction
#include "custom_python/add_level_to_python.h"

namespace Kratos
{

namespace Python
{

template<class TProjectorType>
void MGProjector_AssembleOperator(TProjectorType& rDummy,
        boost::python::list rows, boost::python::list cols, const Matrix& values)
{
    std::vector<std::size_t> _rows(boost::python::len(rows));
    std::vector<std::size_t> _cols(boost::python::len(cols));

    typedef boost::python::stl_input_iterator<int> iterator_value_type;

    std::size_t cnt = 0;
    BOOST_FOREACH(const typename iterator_value_type::value_type& v, std::make_pair(iterator_value_type(rows), iterator_value_type() ) )
    {
        _rows[cnt++] = static_cast<std::size_t>(v);
    }

    cnt = 0;
    BOOST_FOREACH(const typename iterator_value_type::value_type& v, std::make_pair(iterator_value_type(cols), iterator_value_type() ) )
    {
        _cols[cnt++] = static_cast<std::size_t>(v);
    }

    rDummy.AssembleOperator(_rows, _cols, values);
}

template<class TProjectorType>
void MGProjector_AssembleOperatorByBlock(TProjectorType& rDummy,
        boost::python::list rows, boost::python::list cols, const Matrix& values,
        const std::size_t& block_size)
{
    std::vector<std::size_t> _rows(boost::python::len(rows));
    std::vector<std::size_t> _cols(boost::python::len(cols));

    typedef boost::python::stl_input_iterator<int> iterator_value_type;

    std::size_t cnt = 0;
    BOOST_FOREACH(const typename iterator_value_type::value_type& v, std::make_pair(iterator_value_type(rows), iterator_value_type() ) )
    {
        _rows[cnt++] = static_cast<std::size_t>(v);
    }

    cnt = 0;
    BOOST_FOREACH(const typename iterator_value_type::value_type& v, std::make_pair(iterator_value_type(cols), iterator_value_type() ) )
    {
        _cols[cnt++] = static_cast<std::size_t>(v);
    }

    rDummy.AssembleOperator(_rows, _cols, values, block_size);
}

template<class SparseSpaceType, std::size_t TDim>
void MultigridSolversApp_AddStructuredMeshMGProjectorToPython()
{
    using namespace boost::python;

    typedef MeshBasedMGProjector<SparseSpaceType> MeshBasedMGProjectorType;
    typedef MatrixBasedMGProjector<SparseSpaceType> MatrixBasedMGProjectorType;

    std::stringstream ss;
    ss << "StructuredMeshMGProjector" << TDim << "D";
    typedef StructuredMeshMGProjector<SparseSpaceType, TDim> StructuredMeshMGProjectorType;
    class_<StructuredMeshMGProjectorType, typename StructuredMeshMGProjectorType::Pointer, bases<MeshBasedMGProjectorType>, boost::noncopyable>
    (ss.str().c_str(), init<>())
    .def("SetFineMeshSize", &StructuredMeshMGProjectorType::SetFineMeshSize)
    .def("SetCoarseMeshSize", &StructuredMeshMGProjectorType::SetCoarseMeshSize)
    .def(self_ns::str(self))
    ;

    ss.str("");
    ss << "StructuredMeshMGInterpolator" << TDim << "D";
    typedef StructuredMeshMGInterpolator<SparseSpaceType, TDim> StructuredMeshMGInterpolatorType;
    class_<StructuredMeshMGInterpolatorType, typename StructuredMeshMGInterpolatorType::Pointer, bases<StructuredMeshMGProjectorType>, boost::noncopyable>
    (ss.str().c_str(), init<>())
    ;

    ss.str("");
    ss << "StructuredMeshMGTransposeInterpolator" << TDim << "D";
    typedef StructuredMeshMGTransposeInterpolator<SparseSpaceType, TDim> StructuredMeshMGTransposeInterpolatorType;
    class_<StructuredMeshMGTransposeInterpolatorType, typename StructuredMeshMGTransposeInterpolatorType::Pointer, bases<StructuredMeshMGProjectorType>, boost::noncopyable>
    (ss.str().c_str(), init<>())
    ;

    ss.str("");
    ss << "StructuredMeshMatrixBasedMGInterpolator" << TDim << "D";
    typedef StructuredMeshMatrixBasedMGInterpolator<SparseSpaceType, TDim> StructuredMeshMatrixBasedMGInterpolatorType;
    class_<StructuredMeshMatrixBasedMGInterpolatorType, typename StructuredMeshMatrixBasedMGInterpolatorType::Pointer, bases<StructuredMeshMGProjectorType, MatrixBasedMGProjectorType>, boost::noncopyable>
    (ss.str().c_str(), init<>())
    ;

    ss.str("");
    ss << "StructuredMeshMatrixBasedMGTransposeInterpolator" << TDim << "D";
    typedef StructuredMeshMatrixBasedMGTransposeInterpolator<SparseSpaceType, TDim> StructuredMeshMatrixBasedMGTransposeInterpolatorType;
    class_<StructuredMeshMatrixBasedMGTransposeInterpolatorType, typename StructuredMeshMatrixBasedMGTransposeInterpolatorType::Pointer, bases<StructuredMeshMGProjectorType, MatrixBasedMGProjectorType>, boost::noncopyable>
    (ss.str().c_str(), init<>())
    ;

    ss.str("");
    ss << "StructuredMeshMGInjector" << TDim << "D";
    typedef StructuredMeshMGInjector<SparseSpaceType, TDim> StructuredMeshMGInjectorType;
    class_<StructuredMeshMGInjectorType, typename StructuredMeshMGInjectorType::Pointer, bases<StructuredMeshMGProjectorType>, boost::noncopyable>
    (ss.str().c_str(), init<>())
    ;
}

void MultigridSolversApp_AddLevelToPython()
{
    using namespace boost::python;

    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;

    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

    typedef LinearSolver<SparseSpaceType, LocalSpaceType> LinearSolverType;

   //****************************************************************************************************
   // level definition
   //****************************************************************************************************

    typedef MGLevel<SparseSpaceType, LocalSpaceType> MGLevelType;
    class_<MGLevelType, MGLevelType::Pointer, boost::noncopyable>
    ( "MGLevel", init<const typename MGLevelType::IndexType&>())
    .def(self_ns::str(self))
    .def("PreSmoother", &MGLevelType::PreSmoother)
    .def("PostSmoother", &MGLevelType::PostSmoother)
    .def("SetPreSmoother", &MGLevelType::SetPreSmoother)
    .def("SetPostSmoother", &MGLevelType::SetPostSmoother)
    .def("SetRestrictionOperator", &MGLevelType::SetRestrictionOperator)
    .def("SetTransferOperator", &MGLevelType::SetTransferOperator)
    .def("SetProlongationOperator", &MGLevelType::SetProlongationOperator)
    .def("ApplyPreSmoother", &MGLevelType::ApplyPreSmoother)
    .def("ApplyPostSmoother", &MGLevelType::ApplyPostSmoother)
    .def("ApplyRestriction", &MGLevelType::ApplyRestriction)
    .def("ApplyTransfer", &MGLevelType::ApplyTransfer)
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
    .def("ApplyTranspose", &MGProjectorType::ApplyTranspose)
    .def("GetBaseSize", &MGProjectorType::GetBaseSize)
    .def("GetProjectedSize", &MGProjectorType::GetProjectedSize)
    ;

    typedef MGNullProjector<SparseSpaceType> MGNullProjectorType;
    class_<MGNullProjectorType, MGNullProjectorType::Pointer, bases<MGProjectorType>, boost::noncopyable>
    ("MGNullProjector", init<>())
    .def(self_ns::str(self))
    ;

    typedef MGTransposeProjector<SparseSpaceType> MGTransposeProjectorType;
    class_<MGTransposeProjectorType, MGTransposeProjectorType::Pointer, bases<MGProjectorType>, boost::noncopyable>
    ("MGTransposeProjector", init<typename MGProjectorType::Pointer>())
    .def(self_ns::str(self))
    ;

    typedef IndexBasedMGProjector<SparseSpaceType> IndexBasedMGProjectorType;
    class_<IndexBasedMGProjectorType, IndexBasedMGProjectorType::Pointer, bases<MGProjectorType>, boost::noncopyable>
    ("IndexBasedMGProjector", init<const std::size_t&, const std::size_t&>())
    .def("SetStride", &IndexBasedMGProjectorType::SetStride)
    .def("AssembleOperator", &MGProjector_AssembleOperator<IndexBasedMGProjectorType>)
    .def(self_ns::str(self))
    ;

    typedef MatrixBasedMGProjector<SparseSpaceType> MatrixBasedMGProjectorType;
    class_<MatrixBasedMGProjectorType, MatrixBasedMGProjectorType::Pointer, bases<MGProjectorType>, boost::noncopyable>
    ("MatrixBasedMGProjector", init<>())
    .def(init<const std::size_t&, const std::size_t&>())
    .def("SetOperator", &MatrixBasedMGProjectorType::SetOperator)
    .def("AssembleOperator", &MGProjector_AssembleOperator<MatrixBasedMGProjectorType>)
    .def("AssembleOperator", &MGProjector_AssembleOperatorByBlock<MatrixBasedMGProjectorType>)
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

    MultigridSolversApp_AddStructuredMeshMGProjectorToPython<SparseSpaceType, 2>();
    MultigridSolversApp_AddStructuredMeshMGProjectorToPython<SparseSpaceType, 3>();
}

}  // namespace Python.

} // Namespace Kratos
