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
//   Date:                $Date: 20/7/2018 $
//   Revision:            $Revision: 1.0 $
//
//


#if !defined(KRATOS_MULTIGRID_SOLVERS_APP_STRUCTURED_MESH_MATRIX_BASED_MG_TRANSPOSE_INTERPOLATOR_H_INCLUDED )
#define  KRATOS_MULTIGRID_SOLVERS_APP_STRUCTURED_MESH_MATRIX_BASED_MG_TRANSPOSE_INTERPOLATOR_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <sstream>


// External includes


// Project includes
#include "includes/define.h"
#include "custom_utilities/matrix_based_mg_projector.h"
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
 * Implementation of restriction operator for geometric multigrid by transpose of the corresponding prolongator as interpolator.
 * The fine model_part must be double of the coarse model_part.
 * The prolongator assumes that the Dirichlet BC is not removed from the linear system matrix.
 */
template<class TSpaceType, std::size_t TDim>
class StructuredMeshMatrixBasedMGTransposeInterpolator : public MatrixBasedMGProjector<TSpaceType>, public StructuredMeshMGProjector<TSpaceType, TDim>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of StructuredMeshMatrixBasedMGTransposeInterpolator
    KRATOS_CLASS_POINTER_DEFINITION(StructuredMeshMatrixBasedMGTransposeInterpolator);

    typedef MatrixBasedMGProjector<TSpaceType> BaseType1;

    typedef StructuredMeshMGProjector<TSpaceType, TDim> BaseType2;

    typedef typename BaseType1::MatrixType MatrixType;

    typedef typename BaseType1::MatrixPointerType MatrixPointerType;

    typedef typename BaseType1::VectorType VectorType;

    typedef typename BaseType1::VectorPointerType VectorPointerType;

    typedef typename BaseType1::SizeType SizeType;

    typedef typename BaseType1::IndexType IndexType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Empty constructor
    StructuredMeshMatrixBasedMGTransposeInterpolator(): BaseType1(), BaseType2()
    {}

    /// Default constructor.
    StructuredMeshMatrixBasedMGTransposeInterpolator(ModelPart::Pointer p_model_part_coarse, ModelPart::Pointer p_model_part_fine)
    : BaseType1(), BaseType2(p_model_part_coarse, p_model_part_fine)
    {}

    /// Default constructor.
    StructuredMeshMatrixBasedMGTransposeInterpolator(ModelPart::Pointer p_model_part_coarse, ModelPart::Pointer p_model_part_fine, const std::size_t& block_size)
    : BaseType1(), BaseType2(p_model_part_coarse, p_model_part_fine, block_size)
    {}

    /// Destructor.
    virtual ~StructuredMeshMatrixBasedMGTransposeInterpolator()
    {}

    /// Copy constructor
    StructuredMeshMatrixBasedMGTransposeInterpolator(const StructuredMeshMatrixBasedMGTransposeInterpolator& rOther)
    : BaseType1(rOther), BaseType2(rOther)
    {}

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator. It's also important like the Copy constructor
    StructuredMeshMatrixBasedMGTransposeInterpolator& operator= (const StructuredMeshMatrixBasedMGTransposeInterpolator& rOther)
    {
        BaseType1::operator=(rOther);
        BaseType2::operator=(rOther);
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
            new StructuredMeshMatrixBasedMGTransposeInterpolator<TSpaceType, TDim>(p_model_part_coarse, p_model_part_fine));
    }

    /// Create the new instance of the mesh based projector
    virtual typename MeshBasedMGProjector<TSpaceType>::Pointer Create(
        ModelPart::Pointer p_model_part_coarse, ModelPart::Pointer p_model_part_fine,
        const std::size_t& block_size) const
    {
        return typename MeshBasedMGProjector<TSpaceType>::Pointer(
            new StructuredMeshMatrixBasedMGTransposeInterpolator<TSpaceType, TDim>(p_model_part_coarse, p_model_part_fine, block_size));
    }

    /// Set the number of division on a specific direction of the coarse mesh
    void SetDivision(const std::size_t& dim, const std::size_t& num_division)
    {
        this->CoarseMeshSize()[dim] = num_division;
        this->FineMeshSize()[dim] = num_division*2;
    }

    /// Construct the underlying operator
    /// This shall be called after all the Set...() are called
    virtual void Initialize()
    {
        MatrixPointerType pOperator = this->GetOperator();
        if (TSpaceType::Size1(*pOperator) != this->GetProjectedSize()
         || TSpaceType::Size2(*pOperator) != this->GetBaseSize())
            TSpaceType::Resize(*pOperator, this->GetProjectedSize(), this->GetBaseSize());

        // MatrixType P(this->GetBaseSize(), this->GetProjectedSize());
        // StructuredMatrixBasedMGProlongator<TSpaceType, TDim> Op(this->pCoarseModelPart(), this->pFineModelPart());
        // Op.SetBlockSize(this->BlockSize());
        // for (std::size_t dim = 0; dim < TDim; ++dim)
        //     Op.SetDivision(dim, this->FineMeshSize()[dim]);
        // Op.Initialize();
        // TSpaceType::TransposeCopy(*(Op.GetOperator()), *pOperator);

        if (this->pFineModelPart() == NULL || this->pCoarseModelPart() == NULL)
            TSpaceType::SetToZero(*pOperator);
        else
            this->ConstructMatrix(*pOperator);
    }

    /// Apply the projection
    virtual int Apply(VectorType& rX, VectorType& rY) const
    {
        return BaseType1::Apply(rX, rY);
    }

    /// Apply the transpose projection
    virtual int ApplyTranspose(VectorType& rX, VectorType& rY) const
    {
        return BaseType1::ApplyTranspose(rX, rY);
    }

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{

    /// Get the size of the base space
    virtual SizeType GetBaseSize() const
    {
        if (this->pFineModelPart() != NULL)
            return this->pFineModelPart()->NumberOfNodes() * static_cast<SizeType>(this->BlockSize());
        else
            return 0;
    }

    /// Get the size of the projected space
    virtual SizeType GetProjectedSize() const
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
    virtual std::string Info() const
    {
        std::stringstream ss;
        ss << "StructuredMeshMatrixBasedMGTransposeInterpolator";

        ss << ", model_part_coarse: ";
        if (this->pCoarseModelPart() != NULL)
            ss << this->pCoarseModelPart()->Name();
        else
            ss << "null";

        ss << ", model_part_fine: ";
        if (this->pFineModelPart() != NULL)
            ss << this->pFineModelPart()->Name();
        else
            ss << "null";

        ss << ", block_size: " << this->BlockSize();

        return ss.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
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


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    /// Construct the restriction matrix
    void ConstructMatrix(MatrixType& rOperator) const
    {
        // loop through fine nodes
        for(ModelPart::NodeIterator it_node = this->pFineModelPart()->NodesBegin();
            it_node != this->pFineModelPart()->NodesEnd(); ++it_node)
        {
            // get the fine multi index
            std::size_t fine_node_id = it_node->Id();
            #ifdef CHECK_SIZE
            KRATOS_WATCH(fine_node_id)
            KRATOS_WATCH(this->FineMeshSize()[0])
            KRATOS_WATCH(this->FineMeshSize()[1])
            if (fine_node_id*this->BlockSize() > this->GetProjectedSize())
                KRATOS_THROW_ERROR(std::logic_error, "fine node id exceeds the vector size", __FUNCTION__)
            #endif
            MultiIndex<TDim> fine_indices = MultiIndex<TDim>::NodeIdToMultiIndex(fine_node_id, this->FineMeshSize());
            #ifdef CHECK_SIZE
            KRATOS_WATCH(fine_indices)
            #endif

            std::size_t row, col;

            // get the coarse multi index
            if (fine_indices.IsOnCoarse())
            {
                // the fine multiindex coincident with a coarse multiindex
                MultiIndex<TDim> coarse_indices = fine_indices/2;

                // transfer from coarse to fine
                std::size_t coarse_node_id = MultiIndex<TDim>::MultiIndexToNodeId(coarse_indices, this->CoarseMeshSize());
                #ifdef CHECK_SIZE
                std::cout << "1> coarse_node_id: " << coarse_node_id << std::endl;
                if (coarse_node_id*this->BlockSize() > this->GetBaseSize())
                    KRATOS_THROW_ERROR(std::logic_error, "coarse node id exceeds the vector size", __FUNCTION__)
                #endif
                for (unsigned int ib = 0; ib < this->BlockSize(); ++ib)
                {
                    row = (fine_node_id-1)*this->BlockSize() + ib;
                    col = (coarse_node_id-1)*this->BlockSize() + ib;
                    rOperator(col, row) = 1.0;
                }
            }
            else
            {
                // find the neighbors
                std::vector<MultiIndex<TDim> > neighbors_indices = fine_indices.FindCoarseNeighbours();

                double fact = 1.0 / neighbors_indices.size();
                for (std::size_t i = 0; i < neighbors_indices.size(); ++i)
                {
                    std::size_t coarse_node_id = MultiIndex<TDim>::MultiIndexToNodeId(neighbors_indices[i], this->CoarseMeshSize());
                    #ifdef CHECK_SIZE
                    std::cout << "2> "; KRATOS_WATCH(neighbors_indices[i])
                    std::cout << "2> "; KRATOS_WATCH(coarse_node_id)
                    std::cout << "2> "; KRATOS_WATCH(this->CoarseMeshSize()[0])
                    std::cout << "2> "; KRATOS_WATCH(this->CoarseMeshSize()[1])
                    if (coarse_node_id*this->BlockSize() > this->GetBaseSize())
                        KRATOS_THROW_ERROR(std::logic_error, "coarse node id exceeds the vector size", __FUNCTION__)
                    #endif
                    for (unsigned int ib = 0; ib < this->BlockSize(); ++ib)
                    {
                        row = (fine_node_id-1)*this->BlockSize() + ib;
                        col = (coarse_node_id-1)*this->BlockSize() + ib;
                        rOperator(col, row) = fact;
                    }
                }
            }
        }
    }

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
inline std::istream& operator >> (std::istream& IStream, StructuredMeshMatrixBasedMGTransposeInterpolator<TSpaceType, TDim>& rThis)
{
    return IStream;
}

/// output stream function
template<class TSpaceType, std::size_t TDim>
inline std::ostream& operator << (std::ostream& rOStream, const StructuredMeshMatrixBasedMGTransposeInterpolator<TSpaceType, TDim>& rThis)
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

#endif // KRATOS_MULTIGRID_SOLVERS_APP_STRUCTURED_MESH_MATRIX_BASED_MG_TRANSPOSE_INTERPOLATOR_H_INCLUDED  defined

