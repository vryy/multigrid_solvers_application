//
//   Project Name:        Kratos
//   Last modified by:    $Author: hbui $
//   Date:                $Date: Apr 19, 2012 $
//   Revision:            $Revision: 1.0 $
//
//


// System includes

// External includes
#include <boost/python.hpp>

// Project includes
#include "includes/define.h"
#include "spaces/ublas_space.h"
#include "linear_solvers/preconditioner.h"
#include "linear_solvers/linear_solver.h"
#include "linear_solvers/direct_solver.h"
#include "linear_solvers/iterative_solver.h"
#include "custom_linear_solvers/multilevel_solver.h"
#include "custom_linear_solvers/multilevel_preconditioner.h"
#include "custom_linear_solvers/gauss_seidel_iterative_solver.h"
#include "custom_linear_solvers/jacobi_iterative_solver.h"
#include "custom_utilities/multilevel_solver_factory.h"
#include "custom_utilities/ruge_stueben_solver_factory.h"
// #include "custom_utilities/smoothed_aggregation_solver_factory.h"
#include "custom_utilities/gmg_structured_solver_factory.h"


namespace Kratos
{

template<class TParameterListType>
void SetDoubleValue(TParameterListType& dummy, const typename TParameterListType::KeyType& name, double value)
{
    dummy.set(name, value);
}

template<class TParameterListType>
void SetIntValue(TParameterListType& dummy, const typename TParameterListType::KeyType& name, int value)
{
    dummy.set(name, value);
}

template<class TParameterListType>
void SetStringValue(TParameterListType& dummy, const typename TParameterListType::KeyType& name, const std::string value)
{
    dummy.set(name, value);
}

template<class TParameterListType>
void SetBoolValue(TParameterListType& dummy, const typename TParameterListType::KeyType& name, bool value)
{
    dummy.set(name, value);
}

template<class TParameterListType>
TParameterListType& SubList(TParameterListType& dummy, const typename TParameterListType::KeyType& name)
{
    return dummy.sublist(name);
}

template<class TMultilevelSolverType>
typename TMultilevelSolverType::LevelType::Pointer MultilevelSolver_pGetLevel(TMultilevelSolverType& rDummy, int lvl)
{
    return rDummy.pGetLevel(lvl);
}

namespace Python
{

    void MultigridSolversApp_AddLinearSolversToPython()
    {
        typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
        typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

        typedef LinearSolver<SparseSpaceType, LocalSpaceType> LinearSolverType;
        typedef DirectSolver<SparseSpaceType, LocalSpaceType> DirectSolverType;
        typedef IterativeSolver<SparseSpaceType, LocalSpaceType> IterativeSolverType;
        typedef Preconditioner<SparseSpaceType, LocalSpaceType> PreconditionerType;

        using namespace boost::python;

        //***************************************************************************
        //parameter list
        //***************************************************************************

        typedef ParameterList<std::string> ParameterListType;
        class_<ParameterListType, ParameterListType::Pointer, boost::noncopyable>("MGParameterList", init<>())
        .def("SetDoubleValue", SetDoubleValue<ParameterListType>)
        .def("SetIntValue", SetIntValue<ParameterListType>)
        .def("SetStringValue", SetStringValue<ParameterListType>)
        .def("SetBoolValue", SetBoolValue<ParameterListType>)
        .def("SubList", SubList<ParameterListType>, return_internal_reference<>())
        ;

        //***************************************************************************
        //linear solvers
        //***************************************************************************

        typedef GaussSeidelIterativeSolver<SparseSpaceType, LocalSpaceType> GaussSeidelIterativeSolverType;
        class_<GaussSeidelIterativeSolverType, GaussSeidelIterativeSolverType::Pointer, bases<LinearSolverType> >( "GaussSeidelIterativeSolver")
        .def(init< >())
        .def(init<unsigned int, std::string >())
        .def(self_ns::str(self))
        ;

        typedef JacobiIterativeSolver<SparseSpaceType, LocalSpaceType> JacobiIterativeSolverType;
        class_<JacobiIterativeSolverType, JacobiIterativeSolverType::Pointer, bases<LinearSolverType> >( "JacobiIterativeSolver")
        .def(init< >())
        .def(init<unsigned int, double >())
        .def(self_ns::str(self))
        ;

        /* multilevel solver */
        typedef MultilevelSolver<SparseSpaceType, LocalSpaceType> MultilevelSolverType;
        class_<MultilevelSolverType, MultilevelSolverType::Pointer, bases<LinearSolverType> >("MultilevelSolver", init<>())
        .def(init<LinearSolverType::Pointer >())
        .def(init<LinearSolverType::Pointer, std::string >())
        .def(init<double, unsigned int, LinearSolverType::Pointer, std::string >())
        .def(self_ns::str(self))
        .def("AdditionalPhysicalDataIsNeeded", &MultilevelSolverType::AdditionalPhysicalDataIsNeeded)
        .def("ProvideAdditionalData", &MultilevelSolverType::ProvideAdditionalData)
        .def("AddPreSmoother", &MultilevelSolverType::AddPreSmoother)
        .def("ChangePreSmoother", &MultilevelSolverType::ChangePreSmoother)
        .def("AddPostSmoother", &MultilevelSolverType::AddPostSmoother)
        .def("ChangePostSmoother", &MultilevelSolverType::ChangePostSmoother)
        .def("SetUpSmoothers", &MultilevelSolverType::SetUpSmoothers)
        .def("AddLevel", &MultilevelSolverType::AddLevel)
        .def("GetLevel", MultilevelSolver_pGetLevel<MultilevelSolverType>)
        .def("GetNumberOfLevels", &MultilevelSolverType::GetNumberOfLevels)
        .def("SetCycle", &MultilevelSolverType::SetCycle)
        .def("SetCoarseSolver", &MultilevelSolverType::SetCoarseSolver)
        .def("SetTolerance", &MultilevelSolverType::SetTolerance)
        .def("SetMaxIterationsNumber", &MultilevelSolverType::SetMaxIterationsNumber)
        .def("SetMaxLevels", &MultilevelSolverType::SetMaxLevels)
        .def("SetMaxCoarseSize", &MultilevelSolverType::SetMaxCoarseSize)
        .def("SetFactory", &MultilevelSolverType::SetFactory)
        ;

        //****************************************************************************************************
        //preconditioners
        //****************************************************************************************************

        typedef MultilevelPreconditioner<SparseSpaceType, LocalSpaceType> MultilevelPreconditionerType;
        class_<MultilevelPreconditionerType, MultilevelPreconditionerType::Pointer, bases<PreconditionerType> >("MultilevelPreconditioner", init<MultilevelSolverType::Pointer>())
        .def(self_ns::str(self))
        .def("SetMultilevelSolver",&MultilevelPreconditionerType::SetMultilevelSolver)
        ;

        //***************************************************************************
        //factories
        //***************************************************************************

        typedef MultilevelSolverFactory<SparseSpaceType, LocalSpaceType> MultilevelSolverFactoryType;
        class_<MultilevelSolverFactoryType, MultilevelSolverFactoryType::Pointer, boost::noncopyable >
        ( "MultilevelSolverFactory", init<ParameterListType&>())
        .def(self_ns::str(self))
        .def("SetMute", &MultilevelSolverFactoryType::SetMute)
        .def("InitializeMultilevelSolver", &MultilevelSolverFactoryType::InitializeMultilevelSolver)
        .def("GenerateMultilevelSolver", &MultilevelSolverFactoryType::GenerateMultilevelSolver)
        ;

        typedef RugeStuebenSolverFactory<SparseSpaceType, LocalSpaceType> RugeStuebenSolverFactoryType;
        class_<RugeStuebenSolverFactoryType, RugeStuebenSolverFactoryType::Pointer, bases<MultilevelSolverFactoryType>, boost::noncopyable >
        ( "RugeStuebenSolverFactory", init<ParameterListType&>())
        .def(self_ns::str(self))
        ;

        // typedef SmoothedAggregationSolverFactory<SparseSpaceType, LocalSpaceType> SmoothedAggregationSolverFactoryType;
//        class_<SmoothedAggregationSolverFactoryType, SmoothedAggregationSolverFactoryType::Pointer, bases<MultilevelSolverFactoryType>, boost::noncopyable >
//        ( "SmoothedAggregationSolverFactory", init<ParameterListType& >())
//        ;

        typedef GMGStructuredSolverFactory<SparseSpaceType, LocalSpaceType, 2> GMGStructuredSolverFactory2DType;
        class_<GMGStructuredSolverFactory2DType, GMGStructuredSolverFactory2DType::Pointer, bases<MultilevelSolverFactoryType>, boost::noncopyable >
        ( "GMGStructuredSolverFactory2D", init<ParameterListType&>())
        .def("SetModelPart", &GMGStructuredSolverFactory2DType::SetModelPart)
        .def("SetRestrictionOperator", &GMGStructuredSolverFactory2DType::SetRestrictionOperator)
        .def("SetTransferOperator", &GMGStructuredSolverFactory2DType::SetTransferOperator)
        .def("SetProlongationOperator", &GMGStructuredSolverFactory2DType::SetProlongationOperator)
        .def(self_ns::str(self))
        ;

    }

}  // namespace Python.

} // Namespace Kratos

