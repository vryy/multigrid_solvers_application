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


/* *********************************************************
 *
 *   Last modified by:    $Author: hbui $
 *   Date:                $Date: Nov 2, 2014 $
 *   Revision:            $Revision: 1.2 $
 *
 * ***********************************************************/


#if !defined(KRATOS_PARAMETER_LIST_H_INCLUDED )
#define  KRATOS_PARAMETER_LIST_H_INCLUDED

// External includes
#include <boost/variant.hpp>

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "containers/vector_map.h"


namespace Kratos

{

typedef boost::variant<
        bool,
        int,
        unsigned int,
        long,
        long long,
        float,
        double,
        std::string,
        Vector,
        Matrix
        > KratosParameterListAcceptedType;

/**
A ParameterList class offers an alternative to Teuchos::ParameterList
*/
template<class TKeyType>
class ParameterList : public Kratos::VectorMap<TKeyType, KratosParameterListAcceptedType>
{

    class parameter_list_visitor : public boost::static_visitor<int>
    {
    public:
        int operator()(bool b) const { return 0;}
        int operator()(int i) const { return 1;}
        int operator()(unsigned int ui) const { return 2;}
        int operator()(long l) const { return 3;}
        int operator()(long long ll) const { return 4;}
        int operator()(float f) const { return 5;}
        int operator()(double d) const { return 6;}
        int operator()(std::string s) const { return 7;}
        int operator()(Vector& v) const { return 8;}
        int operator()(Matrix& m) const { return 9;}
    };

public:

    KRATOS_CLASS_POINTER_DEFINITION(ParameterList);

    typedef Kratos::VectorMap<TKeyType, KratosParameterListAcceptedType> BaseType;

    typedef typename BaseType::iterator iterator;

    typedef typename BaseType::const_iterator const_iterator;

    typedef typename BaseType::pair_iterator pair_iterator;

    typedef typename BaseType::pair_const_iterator pair_const_iterator;

    typedef typename BaseType::key_type KeyType; //TKeyType

    typedef typename BaseType::data_type DataType; //KratosParameterListAcceptedType

    typedef Kratos::VectorMap<KeyType, Kratos::ParameterList<TKeyType> > SubListType;

    inline ParameterList() : BaseType(), mSubList()
    {}

    inline ParameterList(const BaseType& rOther ) : BaseType(rOther), mSubList()
    {}

    inline ParameterList(const ParameterList& rOther) : BaseType(rOther), mSubList(rOther.mSubList)
    {}

    inline ~ParameterList()
    {}

    inline ParameterList& operator=(const ParameterList& rOther)
    {
        BaseType::operator=(rOther);
        mSubList = rOther.mSubList;
        return *this;
    }

    inline DataType& operator[] (const TKeyType& rKey)
    {
        return BaseType::operator[](rKey);
    }

    template<class TDataType>
    inline void set(const TKeyType& rKey, const TDataType& rData)
    {
        BaseType::operator[](rKey) = rData;
    }

    template<class TDataType>
    inline TDataType& get(const TKeyType& rKey)
    {
        return boost::get<TDataType&>(BaseType::operator[](rKey));
    }

    template<class TDataType>
    inline TDataType& get(const TKeyType& rKey, const TDataType& rInitialData)
    {
        iterator i = BaseType::find(rKey);

        if(i == BaseType::end())
            set(rKey, rInitialData);

        return boost::get<TDataType&>(BaseType::operator[](rKey));
    }

    template<class TOtherKeyType>
    inline std::string& get(const TOtherKeyType& rKey, const char* rInitialData)
    {
        iterator i = BaseType::find(rKey);

        if(i == BaseType::end())
            set(rKey, std::string(rInitialData));

        return boost::get<std::string&>(BaseType::operator[](rKey));
    }

    inline int type(const TKeyType& rKey)
    {
        return boost::apply_visitor(parameter_list_visitor(), BaseType::operator[](rKey));
    }

    inline ParameterList& sublist(const TKeyType& rKey)
    {
        iterator i = BaseType::find(rKey);
        if(i != BaseType::end())
        {
            KRATOS_THROW_ERROR(std::logic_error, "Existing key has been associated with value", "");
        }

        typename SubListType::iterator ii = mSubList.find(rKey);
        if(ii == mSubList.end())
        {
            ParameterList pl;
            mSubList[rKey] = pl;
        }

        return mSubList[rKey];
    }

    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "parameter list (size = " << BaseType::size() + mSubList.size() << ") : ";

        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        BaseType::PrintData(rOStream);
        for(typename SubListType::pair_const_iterator i = mSubList.pair_begin() ; i != mSubList.pair_end() ; ++i)
        {
            rOStream << "(" << (i->first) << " , " << (i->second) << ")" << std::endl;
        }
    }

    void Print(std::ostream& rOStream)
    {
        for(pair_iterator i = BaseType::pair_begin(); i != BaseType::pair_end(); ++i)
        {
            rOStream << "(" << (i->first) << ", " << (i->second) << ", type = " << type(i->first) << ")" << std::endl;
        }
        for(typename SubListType::pair_const_iterator i = mSubList.pair_begin() ; i != mSubList.pair_end() ; ++i)
            rOStream << "(" << (i->first) << " , " << (i->second) << ")" << std::endl;
    }

private:

    SubListType mSubList;

    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
        rSerializer.save("ParameterList", *this);
    }

    virtual void load(Serializer& rSerializer)
    {
        rSerializer.load("ParameterList", *this);
    }

};


template<class TKeyType>
inline std::istream& operator >> (std::istream& is, Kratos::ParameterList<TKeyType>& rThis)
{
    return is;
}


template<class TKeyType>
inline std::ostream& operator << (std::ostream& os, const Kratos::ParameterList<TKeyType>& rThis)
{
    rThis.PrintInfo(os);
    os << std::endl;
    rThis.PrintData(os);
    return os;
}


}

#endif /* KRATOS_PARAMETER_LIST_H_INCLUDED */

