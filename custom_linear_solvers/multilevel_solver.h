/*
see multigrid_solvers_application/LICENSE.txt
*/

//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 2013 Jan 9 11:44:00 $
//   Revision:            $Revision: 1.1 $
//
//


#if !defined(KRATOS_MULTIGRID_SOLVERS_APP_MULTILEVEL_SOLVER_H_INCLUDED )
#define  KRATOS_MULTIGRID_SOLVERS_APP_MULTILEVEL_SOLVER_H_INCLUDED


// System includes
#include <string>
#include <iostream>
#include <sstream>
#include <cstddef>

// External includes
#include "external_includes/pyamg/relaxation.h"

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "linear_solvers/linear_solver.h"
#include "utilities/openmp_utils.h"
#include "custom_utilities/mg_level.h"

//#define DEBUG_PRESMOOTHER


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

// forward declaration
template<class TSparseSpaceType, class TDenseSpaceType, class TModelPartType>
class MultilevelSolverFactory;

template<class TSparseSpaceType, class TDenseSpaceType, class TModelPartType>
class MultilevelPreconditioner;

template<class TSparseSpaceType, class TDenseSpaceType,
         class TModelPartType,
         class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
class MultilevelSolver : public LinearSolver<TSparseSpaceType, TDenseSpaceType, TModelPartType, TReordererType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MultilevelSolver
    KRATOS_CLASS_POINTER_DEFINITION(MultilevelSolver);

    typedef LinearSolver<TSparseSpaceType, TDenseSpaceType, TModelPartType, TReordererType> BaseType;

    typedef typename BaseType::SparseMatrixType SparseMatrixType;

    typedef typename BaseType::VectorType VectorType;

    typedef typename BaseType::DenseMatrixType DenseMatrixType;

    typedef typename BaseType::DenseVectorType DenseVectorType;

    typedef typename BaseType::SizeType SizeType;
    typedef typename BaseType::IndexType IndexType;

    typedef MultilevelSolverFactory<TSparseSpaceType, TDenseSpaceType, TModelPartType> FactoryType;

    typedef typename FactoryType::Pointer FactoryPointerType;

    typedef MGLevel<TSparseSpaceType, TDenseSpaceType, TModelPartType> LevelType;

    typedef typename BaseType::Pointer LinearSolverPointerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MultilevelSolver() :
        mCycle("V"),
        mResidualNorm(0.00),
        mIterationsNumber(0),
        mBNorm(0.00),
        mTolerance(1e-9),
        mMaxIterationsNumber(300),
        mMaxLevels(10),
        mNumPreSmooth(1),
        mNumPostSmooth(1),
        mMaxCoarseSize(500),
        mEchoLevel(0)
    {}

    MultilevelSolver(LinearSolverPointerType pCoarseSolver) :
        mpCoarseSolver(pCoarseSolver),
        mCycle("V"),
        mResidualNorm(0.00),
        mIterationsNumber(0),
        mBNorm(0.00),
        mTolerance(1e-9),
        mMaxIterationsNumber(300),
        mMaxLevels(10),
        mNumPreSmooth(1),
        mNumPostSmooth(1),
        mMaxCoarseSize(500),
        mEchoLevel(0)
    {}

    MultilevelSolver(LinearSolverPointerType pCoarseSolver, const std::string Cycle) :
        mpCoarseSolver(pCoarseSolver),
        mCycle(Cycle),
        mResidualNorm(0.00),
        mIterationsNumber(0),
        mBNorm(0.00),
        mTolerance(1e-9),
        mMaxIterationsNumber(300),
        mMaxLevels(10),
        mNumPreSmooth(1),
        mNumPostSmooth(1),
        mMaxCoarseSize(500),
        mEchoLevel(0)
    {}

    MultilevelSolver(double NewTolerance, IndexType NewMaxIterationsNumber,
            LinearSolverPointerType pCoarseSolver, const std::string Cycle) :
        mpCoarseSolver(pCoarseSolver),
        mCycle(Cycle),
        mResidualNorm(0.00),
        mIterationsNumber(0),
        mBNorm(0.00),
        mTolerance(NewTolerance),
        mMaxIterationsNumber(NewMaxIterationsNumber),
        mMaxLevels(10),
        mNumPreSmooth(1),
        mNumPostSmooth(1),
        mMaxCoarseSize(500),
        mEchoLevel(0)
    {}

    /// Copy constructor.
    MultilevelSolver(const MultilevelSolver& Other) : BaseType(Other)
    {
        mpCoarseSolver = Other.mpCoarseSolver;
        mCycle = Other.mCycle;
        mLevels = Other.mLevels;
        mTolerance = Other.mTolerance;
        mMaxIterationsNumber = Other.mMaxIterationsNumber;
        mMaxLevels = Other.mMaxLevels;
        mMaxCoarseSize = Other.mMaxCoarseSize;
        mNumPreSmooth = Other.mNumPreSmooth;
        mNumPostSmooth = Other.mNumPostSmooth;
        mResidualNorm = 0.00;
        mIterationsNumber = 0;
        mBNorm = 0.00;
        mEchoLevel = Other.mEchoLevel;
    }

    /// Destructor.
    ~MultilevelSolver() override {}

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    MultilevelSolver& operator=(const MultilevelSolver& Other)
    {
        BaseType::operator=(Other);
        mpCoarseSolver = Other.mpCoarseSolver;
        mCycle = Other.mCycle;
        mLevels = Other.mLevels;
        mTolerance = Other.mTolerance;
        mMaxIterationsNumber = Other.mMaxIterationsNumber;
        mMaxLevels = Other.mMaxLevels;
        mMaxCoarseSize = Other.mMaxCoarseSize;
        mNumPreSmooth = Other.mNumPreSmooth;
        mNumPostSmooth = Other.mNumPostSmooth;
        mResidualNorm = 0.00;
        mIterationsNumber = 0;
        mBNorm = 0.00;
        mEchoLevel = Other.mEchoLevel;
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    /// Set the echo level
    void SetEchoLevel(const int& Level)
    {
        mEchoLevel = Level;
    }

    /// Set number of pre-smoothing steps
    void SetNumPreSmooth(const SizeType& NumPreSmooth)
    {
        mNumPreSmooth = NumPreSmooth;
    }

    /// Set number of post-smoothing steps
    void SetNumPostSmooth(const SizeType& NumPostSmooth)
    {
        mNumPostSmooth = NumPostSmooth;
    }

    /** This function is designed to be called as few times as possible. It creates the data structures
     * that only depend on the connectivity of the matrix (and not on its coefficients)
     * so that the memory can be allocated once and expensive operations can be done only when strictly
     * needed
    @param rA. System matrix
    @param rX. Solution vector. it's also the initial guess for iterative linear solvers.
    @param rB. Right hand side vector.
    */
    void Initialize(SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {
        BaseType::Initialize(rA, rX, rB);
    }

    /** Some solvers may require a minimum degree of knowledge of the structure of the matrix. To make an example
     * when solving a mixed u-p problem, it is important to identify the row associated to v and p.
     * another example is the automatic prescription of rotation null-space for smoothed-aggregation solvers
     * which require knowledge on the spatial position of the nodes associated to a given dof.
     * This function tells if the solver requires such data
     */
    bool AdditionalPhysicalDataIsNeeded() override
    {
        return true;
    }

    /** Some solvers may require a minimum degree of knowledge of the structure of the matrix. To make an example
     * when solving a mixed u-p problem, it is important to identify the row associated to v and p.
     * another example is the automatic prescription of rotation null-space for smoothed-aggregation solvers
     * which require knowledge on the spatial position of the nodes associated to a given dof.
     * This function is the place to eventually provide such data
     */
    void ProvideAdditionalData(
        SparseMatrixType& rA,
        VectorType& rX,
        VectorType& rB,
        typename TModelPartType::DofsArrayType& rdof_set,
        TModelPartType& r_model_part
    ) override
    {
        // DO NOTHING
    }

    /** Some solvers may require a minimum degree of knowledge of the structure of the matrix. To make an example
     * when solving a mixed u-p problem, it is important to identify the row associated to v and p.
     * another example is the automatic prescription of rotation null-space for smoothed-aggregation solvers
     * which require knowledge on the spatial position of the nodes associated to a given dof.
     * This function is the place to eventually provide such data
     */
    void ProvideAdditionalDataForLevel(
        const std::size_t& lvl,
        SparseMatrixType& rA,
        VectorType& rX,
        VectorType& rB,
        typename TModelPartType::DofsArrayType& rdof_set,
        TModelPartType& r_model_part
    )
    {
        if (lvl < mLevels.size())
        {
            mPreSmoothers[lvl]->ProvideAdditionalData(rA, rX, rB, rdof_set, r_model_part);
            mPostSmoothers[lvl]->ProvideAdditionalData(rA, rX, rB, rdof_set, r_model_part);
        }

        if (lvl == mLevels.size()-1)
            mpCoarseSolver->ProvideAdditionalData(rA, rX, rB, rdof_set, r_model_part);
    }

    /** This function is designed to be called every time the coefficients change in the system
     * that is, normally at the beginning of each solve.
     * For example if we are implementing a direct solver, this is the place to do the factorization
     * so that then the backward substitution can be performed effectively more than once
    @param rA. System matrix
    @param rX. Solution vector. it's also the initial guess for iterative linear solvers.
    @param rB. Right hand side vector.
    */
    void InitializeSolutionStep(SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {}

    /** Normal solve method.
    Solves the linear system Ax=b and puts the result on SystemVector& rX.
    rVectorx is also th initial guess for iterative methods.
    @param rA. System matrix
    @param rX. Solution vector. it's also the initial
    guess for iterative linear solvers.
     @param rB. Right hand side vector.
    */
    bool Solve(SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {
        double start = OpenMPUtils::GetCurrentTime();

        if(this->IsNotConsistent(rA, rX, rB))
            return false;

        if(mpFactory != NULL)
        {
            mResidualNorm = 0.00;
            mIterationsNumber = 0;
            mBNorm = 0.00;

//            KRATOS_WATCH(GetNumberOfLevels());
            mpFactory->GenerateMultilevelSolver(*this, rA);

            this->SetUpSmoothers();

            std::cout << Info();
        }

        double construct_level_time = OpenMPUtils::GetCurrentTime() - start;
        std::cout << "Time to generate the multigrid hierarchy: " << construct_level_time << std::endl;
        start = OpenMPUtils::GetCurrentTime();

        mBNorm = TSparseSpaceType::TwoNorm(rB);
        if (mEchoLevel > 0)
            std::cout << "Start solving, ||B|| = " << mBNorm << ", Tolerance = " << mTolerance << std::endl;
//        KRATOS_WATCH(&rA);

        if (mpCoarseSolver == NULL)
            KRATOS_ERROR << "The coarse solver is not set";

        const SizeType size = TSparseSpaceType::Size(rX);
        VectorType r(size, 0.00);
        TSparseSpaceType::SetToZero(rX);

        do
        {
            if(this->GetNumberOfLevels() == 1)
            {
                GetLevel(0).Inverse(mpCoarseSolver, rX, rB);

                return true; // the reason to return here is that if the coarse solver is direct solver, it may not solve accurately to the tolerance required (sounds strange but it happens with multiphase_cube example).
            }
            else
            {
                this->RecursiveSolve(0, rX, rB, mCycle);
            }

            // compute the norm

            TSparseSpaceType::Mult(rA, rX, r);

            TSparseSpaceType::UnaliasedAdd(r, -1.00, rB);

            mResidualNorm = TSparseSpaceType::TwoNorm(r);

            ++mIterationsNumber;

            if (mEchoLevel > 0)
                std::cout << " iteration " << mIterationsNumber
                          << ", residual: rel = " << mResidualNorm/mBNorm
                          << ", abs = " << mResidualNorm
                          << std::endl;
        }
        while(IterationNeeded());

        double solve_time = OpenMPUtils::GetCurrentTime() - start;

        std::cout << "#### CONSTRUCT TIME: " << construct_level_time << " ####" << std::endl;
        std::cout << "#### SOLVE TIME: " << solve_time << " ####" << std::endl;
        std::cout << "#### SOLVER TIME: " << construct_level_time + solve_time << " ####" << std::endl;

        return IsConverged();
    }

    /** Multi solve method for solving a set of linear systems with same coefficient matrix.
    Solves the linear system Ax=b and puts the result on SystemVector& rX.
    rVectorx is also th initial guess for iterative methods.
    @param rA. System matrix
    @param rX. Solution vector. it's also the initial
    guess for iterative linear solvers.
     @param rB. Right hand side vector.
    */
    bool Solve(SparseMatrixType& rA, DenseMatrixType& rX, DenseMatrixType& rB) override
    {
        return false;
    }

    void SetFactory(FactoryPointerType Factory)
    {
        mpFactory = Factory;
    }

    void AddPreSmoother(LinearSolverPointerType pPreSmootherSolver)
    {
        mPreSmoothers.push_back(pPreSmootherSolver);
    }

    void ChangePreSmoother(IndexType lvl, LinearSolverPointerType pPreSmootherSolver)
    {
        GetLevel(lvl).SetPreSmoother(pPreSmootherSolver);
    }

    void AddPostSmoother(LinearSolverPointerType pPostSmootherSolver)
    {
        mPostSmoothers.push_back(pPostSmootherSolver);
    }

    void ChangePostSmoother(IndexType lvl, LinearSolverPointerType pPostSmootherSolver)
    {
        GetLevel(lvl).SetPostSmoother(pPostSmootherSolver);
    }

    void SetUpSmoothers()
    {
        if(mPreSmoothers.size() == 0 || mPostSmoothers.size() == 0)
        {
            KRATOS_ERROR << "PreSmoother/PostSmoother must be assigned before solving";
        }

        // set up presmoother
        SizeType num_pre_smoothers = mPreSmoothers.size();

        for(IndexType i = 0; i < this->GetNumberOfLevels(); i++)
        {
            if(i < num_pre_smoothers)
                GetLevel(i).SetPreSmoother(mPreSmoothers[i]);
            else
                GetLevel(i).SetPreSmoother(mPreSmoothers.back());
        }

        // set up postsmoother
        SizeType num_post_smoothers = mPostSmoothers.size();

        for(IndexType i = 0; i < this->GetNumberOfLevels(); i++)
        {
            if(i < num_post_smoothers)
                GetLevel(i).SetPostSmoother(mPostSmoothers[i]);
            else
                GetLevel(i).SetPostSmoother(mPostSmoothers.back());
        }
    }

    ///@}
    ///@name Access
    ///@{

    /// Clean all the level
    void ResetLevel()
    {
        mLevels.clear();
    }

    /// Add the new level to the multigrid hierarchy. Use CreateLevel() for more reliable and cleaner code.
    //this function is kept to be compatible with python; use CreateLevel to avoid copying the memory
    void AddLevel(typename LevelType::Pointer plevel)
    {
//        KRATOS_WATCH(&level);
        if(this->GetNumberOfLevels() < mMaxLevels)
        {
            if(plevel->LevelDepth() == this->GetNumberOfLevels())
            {
                mLevels.push_back(plevel);
            }
            else
            {
                KRATOS_ERROR << "The input level has an invalid level depth of " << plevel->LevelDepth();
            }
        }
        else
        {
            KRATOS_ERROR << "The maximum number of levels is " << mMaxLevels;
        }
    }

    /// Set a level at a specific depth
    void SetLevel(typename LevelType::Pointer plevel)
    {
//        KRATOS_WATCH(&level);
        if(plevel->LevelDepth() < this->GetNumberOfLevels())
        {
            mLevels[plevel->LevelDepth()] = plevel;
        }
        else
        {
            KRATOS_ERROR << "The input level has an exceeded level depth > " << this->GetNumberOfLevels();
        }
    }

    /// Create a new level in the multigrid hierarchy. The level starts with level depth 0.
    void CreateLevel()
    {
        if(this->GetNumberOfLevels() < mMaxLevels)
        {
            typename LevelType::Pointer plevel = typename LevelType::Pointer(new LevelType(mLevels.size()));
            mLevels.push_back(plevel);
        }
        else
        {
            KRATOS_ERROR << "The maximum number of levels is " << mMaxLevels;
        }
    }

    typename LevelType::Pointer pGetLevel(const IndexType& idx)
    {
        return mLevels[idx];
    }

    typename LevelType::Pointer pGetLevel(const IndexType& idx) const
    {
        return mLevels[idx];
    }

    LevelType& GetLevel(const IndexType& idx)
    {
        return *(mLevels[idx]);
    }

    const LevelType& GetLevel(const IndexType& idx) const
    {
        return *(mLevels[idx]);
    }

    LevelType& GetLastLevel()
    {
        return *(mLevels.back());
    }

    const LevelType& GetLastLevel() const
    {
        return *(mLevels.back());
    }

    void SetCoarseSolver(LinearSolverPointerType pCoarseSolver)
    {
        mpCoarseSolver = pCoarseSolver;
    }

    void SetCycle(const std::string Cycle)
    {
        mCycle = Cycle;
    }

    std::string Cycle() const
    {
        return mCycle;
    }

    void SetMaxIterationsNumber(IndexType NewMaxIterationsNumber)
    {
        mMaxIterationsNumber = NewMaxIterationsNumber;
    }

    void SetTolerance(double NewTolerance)
    {
        mTolerance = NewTolerance;
    }

    void SetMaxLevels(IndexType MaxLevels)
    {
        mMaxLevels = MaxLevels;
    }

    IndexType GetMaxLevels() const
    {
        return mMaxLevels;
    }

    void SetMaxCoarseSize(IndexType MaxCoarseSize)
    {
        mMaxCoarseSize = MaxCoarseSize;
    }

    IndexType GetMaxCoarseSize() const
    {
        return mMaxCoarseSize;
    }

    ///@}
    ///@name Inquiry
    ///@{

    SizeType GetNumberOfLevels() const
    {
        return mLevels.size();
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "+++++++++" << std::endl;
        buffer << "Multilevel solver with " << this->Cycle() << " cycle";
        buffer << ", number of level(s) = " << this->GetNumberOfLevels() << std::endl;
        for(SizeType i = 0; i < this->GetNumberOfLevels(); ++i)
        {
//            buffer << mLevels[i].Info() << std::endl;
//            mLevels[i].PrintData(buffer);
            buffer << this->GetLevel(i).Info();
            buffer << std::endl;
        }
        buffer << "+++++++++" << std::endl;
        return  buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
//        for(LevelIteratorType level = mLevels.begin(); level != mLevels.end(); ++level)
//        {
//            rOStream << (*level).Info();
//            (*level).PrintData(rOStream);
//        }

        rOStream << " Summary:" << std::endl;
        if (mBNorm == 0.00)
            if (mResidualNorm != 0.00)
                rOStream << "    Residual ratio : infinite" << std::endl;
            else
                rOStream << "    Residual ratio : 0" << std::endl;
        else
        {
            rOStream << "    Initial Residual ratio : " << mBNorm << std::endl;
            rOStream << "    Final Residual ratio : " << mResidualNorm << std::endl;
            rOStream << "    Residual ratio : " << mResidualNorm / mBNorm << std::endl;
            rOStream << "    Slope : " << (mResidualNorm - mBNorm) / mIterationsNumber << std::endl;
        }

        rOStream << "    Tolerance : " << mTolerance << std::endl;
        rOStream << "    Number of iterations : " << mIterationsNumber << std::endl;
        rOStream << "    Maximum number of iterations : " << mMaxIterationsNumber;
        if (mMaxIterationsNumber == mIterationsNumber)
            rOStream << std::endl << "!!!!!!!!!!!! MULTILEVEL SOLVER NOT CONVERGED WITHIN " << mMaxIterationsNumber << " ITERATIONS!!!!!!!!!!!!";
    }

    ///@}
    ///@name Friends
    ///@{

    friend class MultilevelPreconditioner<TSparseSpaceType, TDenseSpaceType, TModelPartType>;

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

    std::vector<typename LevelType::Pointer> mLevels;
    LinearSolverPointerType mpCoarseSolver;
    std::string mCycle;

    double mResidualNorm;
    SizeType mIterationsNumber;
    double mBNorm;
    double mTolerance;
    SizeType mMaxIterationsNumber;
    SizeType mMaxLevels;
    SizeType mMaxCoarseSize;
    SizeType mNumPreSmooth;
    SizeType mNumPostSmooth;

    FactoryPointerType mpFactory;

    std::vector<LinearSolverPointerType> mPreSmoothers;
    std::vector<LinearSolverPointerType> mPostSmoothers;

    int mEchoLevel;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    void RecursiveSolve(const IndexType& lvl, VectorType& rX, VectorType& rB, const std::string& Cycle)
    {
        int err;

        const SizeType size = TSparseSpaceType::Size(rX);
        VectorType r(size, 0.00);

        double aux[10];

        if (mEchoLevel > 1)
        {
            Indent(std::cout, lvl+2); std::cout << "reporting at level " << lvl << ":" << std::endl;
            ComputeResidual(r, lvl, rX, rB);
            aux[0] = norm_2(r);
            Indent(std::cout, lvl+3); std::cout << lvl << ": residual before ApplyPreSmoother: r = " << aux[0] << std::endl;
        }

        #ifdef DEBUG_PRESMOOTHER
        WriteMatrixMarketVector("r0.mm", r);
        #endif

        // Pre smoothing
        for (SizeType i = 0; i < mNumPreSmooth; ++i)
        {
            err = GetLevel(lvl).ApplyPreSmoother(rX, rB);// ErrorCheck(err, "Error with ApplyPreSmoother at", KRATOS_HERE);
            #ifdef DEBUG_PRESMOOTHER
            ComputeResidual(r, lvl, rX, rB);
            std::stringstream ss;
            ss << "rs" << (i+1) << ".mm";
            WriteMatrixMarketVector(ss.str().c_str(), r);
            #endif
        }

        #ifdef DEBUG_PRESMOOTHER
        KRATOS_ERROR "Debug Presmoother is completed. I think we have to stop.";
        #endif

        // compute residual
        ComputeResidual(r, lvl, rX, rB);

        if (mEchoLevel > 1)
        {
            aux[1] = norm_2(r);
            Indent(std::cout, lvl+3); std::cout << lvl << ": residual after ApplyPreSmoother: r = " << aux[1] << ", reduction factor = " << aux[0]/aux[1] << std::endl;
            // REMARKS: the residual here must be small or start to reduce significantly
        }

        const SizeType csize = GetLevel(lvl).GetCoarseSize();
//        KRATOS_WATCH(csize);

        VectorType cR(csize, 0.00);
        VectorType cX(csize, 0.00);

        // restriction
        err = GetLevel(lvl).ApplyRestriction(r, cR); if (err != 0) KRATOS_ERROR << "Error with ApplyRestriction";

        if (mEchoLevel > 1)
        {
            aux[2] = norm_2(cR);
            Indent(std::cout, lvl+3); std::cout << lvl << ": residual at the coarse level after ApplyRestriction: cR = " << aux[2] << std::endl;
        }

        // solve
        if(lvl == this->GetNumberOfLevels() - 2)
        {
            err = GetLevel(lvl+1).Inverse(mpCoarseSolver, cX, cR); if (err != 0) KRATOS_ERROR << "Error with Coarse Solver";
            // KRATOS_WATCH(norm_2(cX))
        }
        else
        {
            if(Cycle.compare("V") == 0)
            {
                RecursiveSolve(lvl + 1, cX, cR, "V");
            }
            else if(Cycle.compare("W") == 0)
            {
                RecursiveSolve(lvl + 1, cX, cR, Cycle);
                RecursiveSolve(lvl + 1, cX, cR, Cycle);
            }
            else if(Cycle.compare("F") == 0)
            {
                RecursiveSolve(lvl + 1, cX, cR, Cycle);
                RecursiveSolve(lvl + 1, cX, cR, "V");
            }
            else
            {
                KRATOS_ERROR << "valid multilevel cycles are 'V', 'W' or 'F'";
            }
        }

        // prolongation
        VectorType Dx(size, 0.00);
        err = GetLevel(lvl).ApplyProlongation(cX, Dx); if (err != 0) KRATOS_ERROR << "Error with ApplyProlongation";
        // KRATOS_WATCH(norm_2(Dx))
        TSparseSpaceType::UnaliasedAdd(rX, 1.00, Dx);

        if (mEchoLevel > 1)
        {
            ComputeResidual(r, lvl, rX, rB);
            aux[3] = norm_2(r);
            Indent(std::cout, lvl+3); std::cout << lvl << ": residual after coarse solve and prolong: r = " << aux[3] << std::endl;
        }

        // Post smoothing
        for (SizeType i = 0; i < mNumPostSmooth; ++i)
            err = GetLevel(lvl).ApplyPostSmoother(rX, rB);// if (err != 0) KRATOS_ERROR << "Error with ApplyPostSmoother";

        if (mEchoLevel > 1)
        {
            ComputeResidual(r, lvl, rX, rB);
            aux[4] = norm_2(r);
            Indent(std::cout, lvl+3); std::cout << lvl << ": residual after ApplyPostSmoother: r = " << aux[4] << ", reduction factor = " << aux[3]/aux[4] << std::endl;
        }
    }

    inline void Indent(std::ostream& rOStream, const std::size_t& num) const
    {
        for (std::size_t i = 0; i < num; ++i)
            rOStream << " ";
    }

    void ComputeResidual(VectorType& r, const IndexType& lvl, VectorType& rX, VectorType& rB) const
    {
        GetLevel(lvl).Apply(rX, r); // r = A*x
        TSparseSpaceType::UnaliasedAdd(r, -1.00, rB); // r = A*x - b
        TSparseSpaceType::InplaceMult(r, -1.00); // r = b - A*x
    }

    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{

    virtual bool IterationNeeded()
    {
        return (mIterationsNumber < mMaxIterationsNumber) && (mResidualNorm > mTolerance * mBNorm);
    }

    virtual bool IsConverged()
    {
        return (mResidualNorm <= mTolerance * mBNorm);
    }

    ///@}
    ///@name Un accessible methods
    ///@{


    ///@}

}; // Class MultilevelSolver

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}

}  // namespace Kratos.

#undef DEBUG_PRESMOOTHER

#endif // KRATOS_MULTILEVEL_SOLVER_H_INCLUDED  defined
