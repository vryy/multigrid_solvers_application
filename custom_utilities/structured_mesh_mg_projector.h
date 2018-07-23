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
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 23/7/2018 $
//   Revision:            $Revision: 1.0 $
//
//


#if !defined(KRATOS_MULTIGRID_SOLVERS_APP_STRUCTURED_MESH_MG_PROJECTOR_H_INCLUDED )
#define  KRATOS_MULTIGRID_SOLVERS_APP_STRUCTURED_MESH_MG_PROJECTOR_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <sstream>


// External includes


// Project includes
#include "includes/define.h"
#include "custom_utilities/mesh_based_mg_projector.h"

// #define CHECK_SIZE

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * Implementation of prolongation operator for geometric multigrid. The fine model_part must be double of the coarse model_part.
 * The prolongator assumes that the Dirichlet BC is not removed from the linear system matrix.
 */
template<class TSpaceType, std::size_t TDim>
class StructuredMeshMGProjector : public MeshBasedMGProjector<TSpaceType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of StructuredMeshMGProjector
    KRATOS_CLASS_POINTER_DEFINITION(StructuredMeshMGProjector);

    typedef MeshBasedMGProjector<TSpaceType> BaseType;

    typedef typename BaseType::MatrixType MatrixType;

    typedef typename BaseType::MatrixPointerType MatrixPointerType;

    typedef typename BaseType::VectorType VectorType;

    typedef typename BaseType::VectorPointerType VectorPointerType;

    typedef typename BaseType::SizeType SizeType;

    typedef typename BaseType::IndexType IndexType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Empty constructor.
    StructuredMeshMGProjector() : BaseType()
    {}

    /// Default constructor.
    StructuredMeshMGProjector(ModelPart::Pointer p_model_part_coarse, ModelPart::Pointer p_model_part_fine)
    : BaseType(p_model_part_coarse, p_model_part_fine)
    {}

    /// Default constructor.
    StructuredMeshMGProjector(ModelPart::Pointer p_model_part_coarse, ModelPart::Pointer p_model_part_fine, const std::size_t& block_size)
    : BaseType(p_model_part_coarse, p_model_part_fine, block_size)
    {}

    /// Destructor.
    virtual ~StructuredMeshMGProjector()
    {}

    /// Copy constructor
    StructuredMeshMGProjector(const StructuredMeshMGProjector& rOther)
    : BaseType(rOther)
    , m_fine_mesh_size(rOther.m_fine_mesh_size)
    , m_coarse_mesh_size(rOther.m_coarse_mesh_size)
    {}

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator. It's also important like the Copy constructor
    StructuredMeshMGProjector& operator= (const StructuredMeshMGProjector& rOther)
    {
        BaseType::operator=(rOther);
        this->m_fine_mesh_size = rOther.m_fine_mesh_size;
        this->m_coarse_mesh_size = rOther.m_coarse_mesh_size;
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    /// Create the new instance of the mesh based projector
    virtual typename MeshBasedMGProjector<TSpaceType>::Pointer Create(
        ModelPart::Pointer p_model_part_coarse, ModelPart::Pointer p_model_part_fine) const
    {
        return typename MeshBasedMGProjector<TSpaceType>::Pointer(
            new StructuredMeshMGProjector<TSpaceType, TDim>(p_model_part_coarse, p_model_part_fine));
    }

    /// Create the new instance of the mesh based projector
    virtual typename MeshBasedMGProjector<TSpaceType>::Pointer Create(
        ModelPart::Pointer p_model_part_coarse, ModelPart::Pointer p_model_part_fine,
        const std::size_t& block_size) const
    {
        return typename MeshBasedMGProjector<TSpaceType>::Pointer(
            new StructuredMeshMGProjector<TSpaceType, TDim>(p_model_part_coarse, p_model_part_fine, block_size));
    }

    /// Set the number of division on a specific direction of the fine mesh
    void SetFineMeshSize(const std::size_t& dim, const std::size_t& num_division)
    {
        if (dim < TDim)
            m_fine_mesh_size[dim] = num_division;
    }

    /// Set the number of division on a specific direction of the coarse mesh
    void SetCoarseMeshSize(const std::size_t& dim, const std::size_t& num_division)
    {
        if (dim < TDim)
            m_coarse_mesh_size[dim] = num_division;
    }

    ///@}
    ///@name Access
    ///@{

    /// Get the dimension of the structured mesh
    static const std::size_t Dimension() const {return TDim;}

    /// Get the fine mesh size
    const boost::array<std::size_t, TDim>& FineMeshSize() const {return m_fine_mesh_size;}

    /// Get the coarse mesh size
    const boost::array<std::size_t, TDim>& CoarseMeshSize() const {return m_coarse_mesh_size;}

    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream ss;
        ss << "StructuredMeshMGProjector<" << TDim << ">";
        return ss.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        BaseType::PrintInfo(rOStream);
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
    }


    ///@}
    ///@name Friends
    ///@{


    ///@}


protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}


private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    boost::array<std::size_t, TDim> m_fine_mesh_size;
    boost::array<std::size_t, TDim> m_coarse_mesh_size;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{


    ///@}

};

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template<class TSpaceType, std::size_t TDim>
inline std::istream& operator >> (std::istream& IStream, StructuredMeshMGProjector<TSpaceType, TDim>& rThis)
{
    return IStream;
}

/// output stream function
template<class TSpaceType, std::size_t TDim>
inline std::ostream& operator << (std::ostream& rOStream, const StructuredMeshMGProjector<TSpaceType, TDim>& rThis)
{
    rOStream << rThis.Info();
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

} // namespace Kratos.

#undef CHECK_SIZE

#endif // KRATOS_MULTIGRID_SOLVERS_APP_STRUCTURED_MESH_BASED_MG_PROJECTOR_H_INCLUDED  defined

