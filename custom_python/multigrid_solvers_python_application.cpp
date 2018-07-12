//
//   Project Name:        Kratos
//   Last modified by:    $Author: hbui $
//   Date:                $Date: Nov 2, 2014 $
//   Revision:            $Revision: 1.0 $
//
//


// System includes


// External includes
#if defined(KRATOS_PYTHON)
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "multigrid_solvers_application.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_python/add_linear_solvers_to_python.h"

namespace Kratos
{

namespace Python
{

    using namespace boost::python;
    BOOST_PYTHON_MODULE(KratosMultigridSolversApplication)
    {

        class_<KratosMultigridSolversApplication, KratosMultigridSolversApplication::Pointer,
               bases<KratosApplication>, boost::noncopyable>
               ("KratosMultigridSolversApplication");

        MultigridSolversApp_AddUtilitiesToPython();

        MultigridSolversApp_AddLinearSolversToPython();
    }

} // namespace Python.

} // namespace Kratos.

#endif // KRATOS_PYTHON

