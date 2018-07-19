#if !defined(KRATOS_MULTIINDEX_H_INCLUDED )
#define  KRATOS_MULTIINDEX_H_INCLUDED

namespace Kratos
{

template<std::size_t TDim, std::size_t TDivision>
class MultiIndex;

template<std::size_t TDim, std::size_t TDivision>
struct MultiIndex_Helper
{
    static std::vector<MultiIndex<TDim, TDivision> > FindCoarseNeighbours(const MultiIndex<TDim, TDivision>& rIndex)
    {
        KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, "is not implemented");
    }

    static std::size_t MultiIndexToNodeId(const MultiIndex<TDim, TDivision>& mindex,
        const boost::array<std::size_t, TDim>& mesh_size)
    {
        KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, "is not implemented");
    }

    static MultiIndex<TDim, TDivision> NodeIdToMultiIndex(const std::size_t& node_id,
        const boost::array<std::size_t, TDim>& mesh_size)
    {
        KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, "is not implemented");
    }
};

/**
 * Multi-index implementation for geometric multigrid.
 * TDivision is the number of division to obtain fine mesh from coarse mesh in each direction.
 */
template<std::size_t TDim, std::size_t TDivision = 2>
class MultiIndex
{
public:
    MultiIndex() {}

    ~MultiIndex() {}

    std::size_t& operator[](const std::size_t& dim) {return m_indices[dim];}
    const std::size_t& operator[](const std::size_t& dim) const {return m_indices[dim];}

    std::size_t size() const {return TDim;}
    std::size_t division() const {return TDivision;}
    std::size_t division(const std::size_t& dim) const {return TDivision;}

    MultiIndex<TDim>& operator*=(const std::size_t& n)
    {
        for (std::size_t i = 0; i < TDim; ++i)
            this->m_indices[i] *= n;
        return *this;
    }

    MultiIndex<TDim>& operator/=(const std::size_t& n)
    {
        for (std::size_t i = 0; i < TDim; ++i)
            this->m_indices[i] /= n;
        return *this;
    }

    MultiIndex<TDim> operator*(const std::size_t& n) const
    {
        MultiIndex<TDim> res = *this;
        res *= n;
        return res;
    }

    MultiIndex<TDim> operator/(const std::size_t& n) const
    {
        MultiIndex<TDim> res = *this;
        res /= n;
        return res;
    }

    /// Check if the fine node is on coarse node at dimension dim.
    bool IsOnCoarse(const std::size_t& dim) const
    {
        return this->m_indices[dim]%TDivision == 0;
    }

    /// Check if the fine node is on coarse node. TDivision is the number of division to obtain fine mesh from coarse mesh.
    bool IsOnCoarse() const
    {
        for (std::size_t i = 0; i < TDim; ++i)
            if (!this->IsOnCoarse(i))
                return false;
        return true;
    }

    std::vector<std::size_t> FindCoarseNeighbours(const std::size_t& dim) const
    {
        std::size_t coarse = this->m_indices[dim]/TDivision;
        if (this->IsOnCoarse(dim))
            return std::vector<std::size_t>{coarse};
        else
            return std::vector<std::size_t>{coarse, coarse+1};
    }

    std::vector<MultiIndex<TDim> > FindCoarseNeighbours()
    {
        return MultiIndex_Helper<TDim, TDivision>::FindCoarseNeighbours(*this);
    }

    static MultiIndex<TDim> NodeIdToMultiIndex(const std::size_t& node_id,
        const boost::array<std::size_t, TDim>& mesh_size)
    {
        return MultiIndex_Helper<TDim, TDivision>::NodeIdToMultiIndex(node_id, mesh_size);
    }

    static std::size_t MultiIndexToNodeId(const MultiIndex<TDim>& mindex,
        const boost::array<std::size_t, TDim>& mesh_size)
    {
        return MultiIndex_Helper<TDim, TDivision>::MultiIndexToNodeId(mindex, mesh_size);
    }

private:
    boost::array<std::size_t, TDim> m_indices;
};

/// output stream function
template<std::size_t TDim, std::size_t TDivision>
inline std::ostream& operator << (std::ostream& rOStream, const MultiIndex<TDim, TDivision>& rThis)
{
    rOStream << "(" << rThis[0];
    for (std::size_t i = 1; i < TDim; ++i)
        rOStream << ", " << rThis[i];
    rOStream << ")";
    return rOStream;
}

template<std::size_t TDivision>
struct MultiIndex_Helper<2, TDivision>
{
    static std::vector<MultiIndex<2, TDivision> > FindCoarseNeighbours(const MultiIndex<2, TDivision>& rIndex)
    {
        std::vector<MultiIndex<2, TDivision> > neighbors_indices;

        std::vector<std::vector<std::size_t> > neighbors(2);
        for (std::size_t i = 0; i < 2; ++i)
            neighbors[i] = rIndex.FindCoarseNeighbours(i);

        for (std::size_t i = 0; i < neighbors[0].size(); ++i)
        {
            for (std::size_t j = 0; j < neighbors[1].size(); ++j)
            {
                MultiIndex<2, TDivision> neighbor_indices;
                neighbor_indices[0] = neighbors[0][i];
                neighbor_indices[1] = neighbors[1][j];
                neighbors_indices.push_back(neighbor_indices);
            }
        }

        return neighbors_indices;
    }

    static std::size_t MultiIndexToNodeId(const MultiIndex<2, TDivision>& mindex,
        const boost::array<std::size_t, 2>& mesh_size)
    {
        return mindex[1]*(1+mesh_size[0]) + mindex[0] + 1;
    }

    static MultiIndex<2, TDivision> NodeIdToMultiIndex(const std::size_t& node_id,
        const boost::array<std::size_t, 2>& mesh_size)
    {
        MultiIndex<2, TDivision> mindex;

        mindex[0] = (node_id-1)%(1+mesh_size[0]);
        mindex[1] = (node_id-1)/(1+mesh_size[0]);

        return mindex;
    }
};

template<std::size_t TDivision>
struct MultiIndex_Helper<3, TDivision>
{
    static std::vector<MultiIndex<3, TDivision> > FindCoarseNeighbours(const MultiIndex<3, TDivision>& rIndex)
    {
        std::vector<MultiIndex<3, TDivision> > neighbors_indices;

        std::vector<std::vector<std::size_t> > neighbors(3);
        for (std::size_t i = 0; i < 3; ++i)
            neighbors[i] = rIndex.FindCoarseNeighbours(i);

        for (std::size_t i = 0; i < neighbors[0].size(); ++i)
        {
            for (std::size_t j = 0; j < neighbors[1].size(); ++j)
            {
                for (std::size_t k = 0; k < neighbors[2].size(); ++k)
                {
                    MultiIndex<3, TDivision> neighbor_indices;
                    neighbor_indices[0] = neighbors[0][i];
                    neighbor_indices[1] = neighbors[1][j];
                    neighbor_indices[2] = neighbors[2][k];
                    neighbors_indices.push_back(neighbor_indices);
                }
            }
        }

        return neighbors_indices;
    }

    static std::size_t MultiIndexToNodeId(const MultiIndex<3>& mindex,
        const boost::array<std::size_t, 3>& mesh_size)
    {
        return mindex[2]*(1+mesh_size[0])*(1+mesh_size[1]) + mindex[1]*(1+mesh_size[0]) + mindex[0] + 1;
    }

    static MultiIndex<3, TDivision> NodeIdToMultiIndex(const std::size_t& node_id,
        const boost::array<std::size_t, 3>& mesh_size)
    {
        MultiIndex<3, TDivision> mindex;

        mindex[0] = (node_id-1)%(1+mesh_size[0]);
        mindex[2] = (node_id-1)/((1+mesh_size[0])*(1+mesh_size[1]));
        mindex[1] = (node_id-1)/(1+mesh_size[0]) - mindex[2]*(1+mesh_size[1]);

        return mindex;
    }
};

}

#endif // KRATOS_MULTIINDEX_H_INCLUDED
