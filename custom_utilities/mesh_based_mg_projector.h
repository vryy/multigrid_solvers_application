/*
see multigrid_solvers_application/LICENSE.txt
*/

//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 23/7/2018 $
//   Revision:            $Revision: 1.0 $
//
//


#if !defined(KRATOS_MULTIGRID_SOLVERS_APP_MESH_BASED_MG_PROJECTOR_H_INCLUDED )
#define  KRATOS_MULTIGRID_SOLVERS_APP_MESH_BASED_MG_PROJECTOR_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <sstream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "custom_utilities/multi_index.h"
#include "custom_utilities/mg_projector.h"

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
 * The restrictor assumes that the Dirichlet BC is not removed from the linear system matrix.
 */
template<class TSpaceType>
class MeshBasedMGProjector : public MGProjector<TSpaceType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MeshBasedMGProjector
    KRATOS_CLASS_POINTER_DEFINITION(MeshBasedMGProjector);

    typedef MGProjector<TSpaceType> BaseType;

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
    MeshBasedMGProjector() : BaseType()
    {}

    /// Default constructor.
    MeshBasedMGProjector(ModelPart::Pointer p_model_part_coarse, ModelPart::Pointer p_model_part_fine)
    : BaseType(), mp_model_part_coarse(p_model_part_coarse), mp_model_part_fine(p_model_part_fine), m_block_size(1)
    {}

    /// Default constructor.
    MeshBasedMGProjector(ModelPart::Pointer p_model_part_coarse, ModelPart::Pointer p_model_part_fine, const int& block_size)
    : BaseType(), mp_model_part_coarse(p_model_part_coarse), mp_model_part_fine(p_model_part_fine), m_block_size(block_size)
    {}

    /// Destructor.
    ~MeshBasedMGProjector() override
    {}

    /// Copy constructor
    MeshBasedMGProjector(const MeshBasedMGProjector& rOther)
    : BaseType(), mp_model_part_coarse(rOther.mp_model_part_coarse)
    , mp_model_part_fine(rOther.mp_model_part_fine), m_block_size(rOther.m_block_size)
    {}

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator. It's also important like the Copy constructor
    MeshBasedMGProjector& operator= (const MeshBasedMGProjector& rOther)
    {
        BaseType::operator=(rOther);
        this->mp_model_part_coarse = rOther.mp_model_part_coarse;
        this->mp_model_part_fine = rOther.mp_model_part_fine;
        this->m_block_size = rOther.m_block_size;
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
            new MeshBasedMGProjector<TSpaceType>(p_model_part_coarse, p_model_part_fine));
    }

    /// Create the new instance of the mesh based projector
    virtual typename MeshBasedMGProjector<TSpaceType>::Pointer Create(
        ModelPart::Pointer p_model_part_coarse, ModelPart::Pointer p_model_part_fine,
        const std::size_t& block_size) const
    {
        return typename MeshBasedMGProjector<TSpaceType>::Pointer(
            new MeshBasedMGProjector<TSpaceType>(p_model_part_coarse, p_model_part_fine, block_size));
    }

    /// Set the number of d.o.fs per node
    void SetBlockSize(const std::size_t& block_size) {m_block_size = block_size;}

    /// Set the coarse model_part
    void SetCoarseModelPart(ModelPart::Pointer p_model_part) {mp_model_part_coarse = p_model_part;}

    /// Set the fine model_part
    void SetFineModelPart(ModelPart::Pointer p_model_part) {mp_model_part_fine = p_model_part;}

    ///@}
    ///@name Access
    ///@{

    /// Get the block size
    std::size_t BlockSize() const
    {
        return m_block_size;
    }

    /// Get the coarse model_part
    ModelPart& CoarseModelPart() {return *mp_model_part_coarse;}
    const ModelPart& CoarseModelPart() const {return *mp_model_part_coarse;}
    ModelPart::Pointer pCoarseModelPart() {return mp_model_part_coarse;}
    ModelPart::Pointer pCoarseModelPart() const {return mp_model_part_coarse;}

    /// Get the fine model_part
    ModelPart& FineModelPart() {return *mp_model_part_fine;}
    const ModelPart& FineModelPart() const {return *mp_model_part_fine;}
    ModelPart::Pointer pFineModelPart() {return mp_model_part_fine;}
    ModelPart::Pointer pFineModelPart() const {return mp_model_part_fine;}

    ///@}
    ///@name Inquiry
    ///@{

    /// Get the size of the base space
    SizeType GetBaseSize() const override
    {
        KRATOS_ERROR << "Calling base class function";
    }

    /// Get the size of the projected space
    SizeType GetProjectedSize() const override
    {
        KRATOS_ERROR << "Calling base class function";
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream ss;
        ss << "MeshBasedMGProjector";
        return ss.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << ", model_part_coarse: ";
        if (mp_model_part_coarse != NULL)
            rOStream << mp_model_part_coarse->Name();
        else
            rOStream << "null";

        rOStream << ", model_part_fine: ";
        if (mp_model_part_fine != NULL)
            rOStream << mp_model_part_fine->Name();
        else
            rOStream << "null";

        rOStream << ", block_size: " << m_block_size;
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

    ModelPart::Pointer mp_model_part_coarse;
    ModelPart::Pointer mp_model_part_fine;
    std::size_t m_block_size;

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

#endif // KRATOS_MULTIGRID_SOLVERS_APP_MESH_BASED_MG_PROJECTOR_H_INCLUDED  defined
