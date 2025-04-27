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
//   Date:                $Date: 16/7/2018 $
//   Revision:            $Revision: 1.0 $
//
//


#if !defined(KRATOS_MULTIGRID_SOLVERS_APP_STRUCTURED_MESH_MG_INJECTOR_H_INCLUDED )
#define  KRATOS_MULTIGRID_SOLVERS_APP_STRUCTURED_MESH_MG_INJECTOR_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <sstream>


// External includes


// Project includes
#include "includes/define.h"
#include "custom_utilities/structured_mesh_mg_projector.h"

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
 * Implementation of restriction operator for geometric multigrid by simple injection.
 * The fine model_part must be double of the coarse model_part.
 * The restrictor assumes that the Dirichlet BC is not removed from the linear system matrix.
 */
template<class TSpaceType, std::size_t TDim>
class StructuredMeshMGInjector : public StructuredMeshMGProjector<TSpaceType, TDim>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of StructuredMeshMGInjector
    KRATOS_CLASS_POINTER_DEFINITION(StructuredMeshMGInjector);

    typedef StructuredMeshMGProjector<TSpaceType, TDim> BaseType;

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
    StructuredMeshMGInjector() : BaseType()
    {}

    /// Default constructor.
    StructuredMeshMGInjector(ModelPart::Pointer p_model_part_coarse, ModelPart::Pointer p_model_part_fine)
    : BaseType(p_model_part_coarse, p_model_part_fine)
    {}

    /// Default constructor.
    StructuredMeshMGInjector(ModelPart::Pointer p_model_part_coarse, ModelPart::Pointer p_model_part_fine, const std::size_t& block_size)
    : BaseType(p_model_part_coarse, p_model_part_fine, block_size)
    {}

    /// Destructor.
    ~StructuredMeshMGInjector() override
    {}

    /// Copy constructor
    StructuredMeshMGInjector(const StructuredMeshMGInjector& rOther)
    : BaseType(rOther)
    {}

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator. It's also important like the Copy constructor
    StructuredMeshMGInjector& operator= (const StructuredMeshMGInjector& rOther)
    {
        BaseType::operator=(rOther);
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    /// Create the new instance of the mesh based projector
    typename MeshBasedMGProjector<TSpaceType>::Pointer Create(
        ModelPart::Pointer p_model_part_coarse, ModelPart::Pointer p_model_part_fine) const override
    {
        return typename MeshBasedMGProjector<TSpaceType>::Pointer(
            new StructuredMeshMGInjector<TSpaceType, TDim>(p_model_part_coarse, p_model_part_fine));
    }

    /// Create the new instance of the mesh based projector
    typename MeshBasedMGProjector<TSpaceType>::Pointer Create(
        ModelPart::Pointer p_model_part_coarse, ModelPart::Pointer p_model_part_fine,
        const std::size_t& block_size) const override
    {
        return typename MeshBasedMGProjector<TSpaceType>::Pointer(
            new StructuredMeshMGInjector<TSpaceType, TDim>(p_model_part_coarse, p_model_part_fine, block_size));
    }

    /// Set the number of division on a specific direction of the coarse mesh
    void SetDivision(const std::size_t& dim, const std::size_t& num_division)
    {
        this->CoarseMeshSize()[dim] = num_division;
        this->FineMeshSize()[dim] = num_division*2;
    }

    /// Apply the projection by injection
    int Apply(VectorType& rX, VectorType& rY) const override
    {
        #ifdef CHECK_SIZE
        KRATOS_WATCH(rX.size())
        KRATOS_WATCH(rY.size())
        #endif

        TSpaceType::SetToZero(rY);

        std::size_t row, col;

        // loop through coarse nodes
        for(ModelPart::NodeIterator it_node = this->pCoarseModelPart()->NodesBegin();
            it_node != this->pCoarseModelPart()->NodesEnd(); ++it_node)
        {
            // get the coarse multi index
            std::size_t coarse_node_id = it_node->Id();
            MultiIndex<TDim> coarse_indices = MultiIndex<TDim>::NodeIdToMultiIndex(coarse_node_id, this->CoarseMeshSize());

            // get the fine multi index
            MultiIndex<TDim> fine_indices = coarse_indices*2;

            // transfer from fine to coarse
            std::size_t fine_node_id = MultiIndex<TDim>::MultiIndexToNodeId(fine_indices, this->FineMeshSize());
            for (unsigned int ib = 0; ib < this->BlockSize(); ++ib)
            {
                row = (coarse_node_id-1)*this->BlockSize() + ib;
                col = (fine_node_id-1)*this->BlockSize() + ib;
                rY[row] = rX[col];
            }
        }

        return 0;
    }

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{

    /// Get the size of the base space
    SizeType GetBaseSize() const override
    {
        if (this->pFineModelPart() != NULL)
            return this->pFineModelPart()->NumberOfNodes() * static_cast<SizeType>(this->BlockSize());
        else
            return 0;
    }

    /// Get the size of the projected space
    SizeType GetProjectedSize() const override
    {
        if (this->pCoarseModelPart() != NULL)
            return this->pCoarseModelPart()->NumberOfNodes() * static_cast<SizeType>(this->BlockSize());
        else
            return 0;
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream ss;
        ss << "StructuredMeshMGInjector<" << TDim << ">";
        return ss.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
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


///@}

} // namespace Kratos.

#undef CHECK_SIZE

#endif // KRATOS_MULTIGRID_SOLVERS_APP_STRUCTURED_MESH_MG_INJECTOR_H_INCLUDED  defined
