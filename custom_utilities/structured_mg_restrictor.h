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


#if !defined(KRATOS_MULTIGRID_SOLVERS_APP_STRUCTURED_MG_RESTRICTOR_H_INCLUDED )
#define  KRATOS_MULTIGRID_SOLVERS_APP_STRUCTURED_MG_RESTRICTOR_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <sstream>


// External includes


// Project includes
#include "includes/define.h"
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
template<class TSpaceType, std::size_t TDim>
class StructuredMGRestrictor : public MGProjector<TSpaceType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of StructuredMGRestrictor
    KRATOS_CLASS_POINTER_DEFINITION(StructuredMGRestrictor);

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
    StructuredMGRestrictor() : BaseType()
    {}

    /// Default constructor.
    StructuredMGRestrictor(ModelPart::Pointer p_model_part_coarse, ModelPart::Pointer p_model_part_fine)
    : BaseType(), mp_model_part_coarse(p_model_part_coarse), mp_model_part_fine(p_model_part_fine), m_block_size(1)
    {}

    /// Destructor.
    virtual ~StructuredMGRestrictor()
    {}

    /// Copy constructor
    StructuredMGRestrictor(const StructuredMGRestrictor& rOther)
    : BaseType(), mp_model_part_coarse(rOther.mp_model_part_coarse)
    , mp_model_part_fine(rOther.mp_model_part_fine), m_block_size(rOther.m_block_size)
    {}

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator. It's also important like the Copy constructor
    StructuredMGRestrictor& operator= (const StructuredMGRestrictor& rOther)
    {
        BaseType::operator=(rOther);
        this.mp_model_part_coarse = rOther.mp_model_part_coarse;
        this.mp_model_part_fine = rOther.mp_model_part_fine;
        this.m_block_size = rOther.m_block_size;
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    /// Set the number of d.o.fs per node
    void SetBlockSize(const int& block_size)
    {
        m_block_size = block_size;
    }

    /// Set the number of division on a specific direction of the coarse mesh
    void SetDivision(const std::size_t& dim, const std::size_t& num_division)
    {
        m_coarse_mesh_size[dim] = num_division;
        m_fine_mesh_size[dim] = num_division*2;
    }

    /// Apply the projection
    virtual int Apply(VectorType& rX, VectorType& rY) const
    {
        #ifdef CHECK_SIZE
        KRATOS_WATCH(rX.size())
        KRATOS_WATCH(rY.size())
        #endif

        TSpaceType::SetToZero(rY);

        // loop through fine nodes
        for(ModelPart::NodeIterator it_node = mp_model_part_fine->NodesBegin();
            it_node != mp_model_part_fine->NodesEnd(); ++it_node)
        {
            // get the fine multi index
            std::size_t fine_node_id = it_node->Id();
            #ifdef CHECK_SIZE
            KRATOS_WATCH(fine_node_id)
            KRATOS_WATCH(m_fine_mesh_size[0])
            KRATOS_WATCH(m_fine_mesh_size[1])
            if (fine_node_id*m_block_size > this->GetProjectedSize())
                KRATOS_THROW_ERROR(std::logic_error, "fine node id exceeds the vector size", __FUNCTION__)
            #endif
            MultiIndex<TDim> fine_indices = MultiIndex<TDim>::NodeIdToMultiIndex(fine_node_id, m_fine_mesh_size);
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
                std::size_t coarse_node_id = MultiIndex<TDim>::MultiIndexToNodeId(coarse_indices, m_coarse_mesh_size);
                #ifdef CHECK_SIZE
                std::cout << "1> coarse_node_id: " << coarse_node_id << std::endl;
                if (coarse_node_id*m_block_size > this->GetBaseSize())
                    KRATOS_THROW_ERROR(std::logic_error, "coarse node id exceeds the vector size", __FUNCTION__)
                #endif
                for (unsigned int ib = 0; ib < m_block_size; ++ib)
                {
                    row = (fine_node_id-1)*m_block_size + ib;
                    col = (coarse_node_id-1)*m_block_size + ib;
                    rY[col] += rX[row];
                }
            }
            else
            {
                // find the neighbors
                std::vector<MultiIndex<TDim> > neighbors_indices = fine_indices.FindCoarseNeighbours();

                double fact = 1.0 / neighbors_indices.size();
                for (std::size_t i = 0; i < neighbors_indices.size(); ++i)
                {
                    std::size_t coarse_node_id = MultiIndex<TDim>::MultiIndexToNodeId(neighbors_indices[i], m_coarse_mesh_size);
                    #ifdef CHECK_SIZE
                    std::cout << "2> "; KRATOS_WATCH(neighbors_indices[i])
                    std::cout << "2> "; KRATOS_WATCH(coarse_node_id)
                    std::cout << "2> "; KRATOS_WATCH(m_coarse_mesh_size[0])
                    std::cout << "2> "; KRATOS_WATCH(m_coarse_mesh_size[1])
                    if (coarse_node_id*m_block_size > this->GetBaseSize())
                        KRATOS_THROW_ERROR(std::logic_error, "coarse node id exceeds the vector size", __FUNCTION__)
                    #endif
                    for (unsigned int ib = 0; ib < m_block_size; ++ib)
                    {
                        row = (fine_node_id-1)*m_block_size + ib;
                        col = (coarse_node_id-1)*m_block_size + ib;
                        rY[col] += fact*rX[row];
                    }
                }
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
    virtual SizeType GetBaseSize() const
    {
        if (mp_model_part_fine != NULL)
            return mp_model_part_fine->NumberOfNodes() * static_cast<SizeType>(m_block_size);
        else
            return 0;
    }

    /// Get the size of the projected space
    virtual SizeType GetProjectedSize() const
    {
        if (mp_model_part_coarse != NULL)
            return mp_model_part_coarse->NumberOfNodes() * static_cast<SizeType>(m_block_size);
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
        ss << "StructuredMGRestrictor";

        ss << ", model_part_coarse: ";
        if (mp_model_part_coarse != NULL)
            ss << mp_model_part_coarse->Name();
        else
            ss << "null";

        ss << ", model_part_fine: ";
        if (mp_model_part_fine != NULL)
            ss << mp_model_part_fine->Name();
        else
            ss << "null";

        ss << ", block_size: " << m_block_size;

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

    ModelPart::Pointer mp_model_part_coarse;
    ModelPart::Pointer mp_model_part_fine;
    int m_block_size;
    boost::array<std::size_t, TDim> m_coarse_mesh_size;
    boost::array<std::size_t, TDim> m_fine_mesh_size;

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
inline std::istream& operator >> (std::istream& IStream, StructuredMGRestrictor<TSpaceType, TDim>& rThis)
{
    return IStream;
}

/// output stream function
template<class TSpaceType, std::size_t TDim>
inline std::ostream& operator << (std::ostream& rOStream, const StructuredMGRestrictor<TSpaceType, TDim>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


} // namespace Kratos.

#undef CHECK_SIZE

#endif // KRATOS_MULTIGRID_SOLVERS_APP_STRUCTURED_MG_RESTRICTOR_H_INCLUDED  defined

