#importing the Kratos Library
from KratosMultiphysics import *
# Check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()
import os,time

### This strategy supports NR iteration for nonlinear problems.
class SolvingStrategyPython:
    #######################################################################
    def __init__( self, model_list, gmg_level_list, Parameters ):
        #save the input parameters
        self.model_list = model_list
        self.model_solver_list = []
        for model in model_list:
            self.model_solver_list.append(model.solver.solver)
        self.gmg_level_list = gmg_level_list
        self.Parameters = Parameters
        self.PrintSparsity = self.Parameters['print_sparsity_info_flag']
        #default values for some variables
        self.max_iter = 30
        self.echo_level = 1

        #local matrices and vectors
        for i in range(0, len(self.model_solver_list)):
            self.model_solver_list[i].pA = gmg_level_list[i].GetCoarseMatrix()
            self.model_solver_list[i].pDx = gmg_level_list[i].GetCoarseUpdateVector()
            self.model_solver_list[i].pb = gmg_level_list[i].GetCoarseVector()

            self.model_solver_list[i].A = (self.model_solver_list[i].pA).GetReference()
            self.model_solver_list[i].Dx = (self.model_solver_list[i].pDx).GetReference()
            self.model_solver_list[i].b = (self.model_solver_list[i].pb).GetReference()

        self.solveCounter = 0

        self.attached_processes = []

    #######################################################################
    def Initialize(self):
        for model_solver in self.model_solver_list:
            if(model_solver.scheme.SchemeIsInitialized() == False):
                model_solver.scheme.Initialize(model_solver.model_part)
            if (model_solver.scheme.ElementsAreInitialized() == False):
                model_solver.scheme.InitializeElements(model_solver.model_part)
            if (model_solver.scheme.ConditionsAreInitialized() == False):
                model_solver.scheme.InitializeConditions(model_solver.model_part)

    #######################################################################
    def Solve( self, time, from_deac, to_deac, from_reac, to_reac ):
        for model in self.model_list:
            model.deac.Reactivate( model.model_part, from_reac, to_reac )
            model.deac.Deactivate( model.model_part, from_deac, to_deac )
            model.model_part.CloneTimeStep(time)
        self.SolveStep()

    #######################################################################
    def SolveStep(self):
        #print self.model_part
        ## - storing original condition size before adding virtual conditions.
        ## - performing contact search
        ## - creating virtual link conditions for the assembling
        self.solveCounter = self.solveCounter + 1
        self.PerformNewtonRaphsonIteration()
        #finalize the solution step
        self.FinalizeSolutionStep()
        #clear if needed - deallocates memory
        for model_solver in self.model_solver_list:
            if(model_solver.ReformDofSetAtEachStep == True):
                model_solver.Clear()

    #######################################################################
    def PerformNewtonRaphsonIteration( self ):
        model0 = self.model_list[0]
        model_solver0 = self.model_solver_list[0]
        print("time = " + str(model0.model_part.ProcessInfo[TIME]))
        for model_solver in self.model_solver_list:
            #perform the operations to be performed ONCE and ensure they will not be repeated
            # elemental function "Initialize" is called here
            if(model_solver.InitializeWasPerformed == False):
                model_solver.Initialize()
                model_solver.InitializeWasPerformed = True
            #perform initializations for the current step
            #this operation implies:
            #identifying the set of DOFs that will be solved during this step
            #organizing the DOFs so to identify the dirichlet conditions
            #resizing the matrix preallocating the "structure"
            if (model_solver.SolutionStepIsInitialized == False):
                model_solver.InitializeSolutionStep()
                model_solver.SolutionStepIsInitialized = True

        #perform prediction
        self.Predict()

        #execute iteration - first iteration is ALWAYS executed
        calculate_norm = False
        self.iterationCounter = 0
        self.iterationCounter = self.iterationCounter + 1
        normDx = self.ExecuteIteration(self.echo_level,calculate_norm)

        #non linear loop
        converged = False
        it = 0
        while(it < self.max_iter and converged == False):
            #verify convergence
            converged = model_solver0.convergence_criteria.PreCriteria(model_solver0.model_part,model_solver0.builder_and_solver.GetDofSet(),model_solver0.A,model_solver0.Dx,model_solver0.b)

            #calculate iteration
            # - system is built and solved
            # - database is updated depending on the solution
            # - nodal coordinates are updated if required
            self.iterationCounter = self.iterationCounter + 1
            normDx = self.ExecuteIteration(self.echo_level,calculate_norm)

            #verify convergence
            converged = model_solver0.convergence_criteria.PostCriteria(model_solver0.model_part,model_solver0.builder_and_solver.GetDofSet(),model_solver0.A,model_solver0.Dx,model_solver0.b)

            #update iteration count
            it = it + 1

        if( it == self.max_iter and converged == False):
            time = model0.model_part.ProcessInfo[TIME]
            print("Iteration did not converge at time step " + str(time))
            if('stop_Newton_Raphson_if_not_converge' in self.Parameters):
                if(self.Parameters['stop_Newton_Raphson_if_not_converge'] == True):
                    sys.exit("Sorry, my boss does not allow me to continue. The time step did not converge at time step " + str(time) + ", it = " + str(it) + ", max_iter = " + str(self.max_iter))
                else:
                    print('However, the iteration will still be proceeded' + ", it = " + str(it) + ", max_iter = " + str(self.max_iter))
            else:
                sys.exit("Sorry, my boss does not allow me to continue. The time step did not converge at time step " + str(time) + ", it = " + str(it) + ", max_iter = " + str(self.max_iter))
        print("PerformNewtonRaphsonIteration converged after " + str(it) + " steps")

    #######################################################################
    def Predict(self):
        for model_solver in self.model_solver_list:
            model_solver.scheme.Predict(model_solver.model_part,model_solver.builder_and_solver.GetDofSet(),model_solver.A,model_solver.Dx,model_solver.b)

    #######################################################################
    def InitializeSolutionStep(self):
        for model_solver in self.model_solver_list:
            if(model_solver.builder_and_solver.GetDofSetIsInitializedFlag() == False or model_solver.ReformDofSetAtEachStep == True):
                #initialize the list of degrees of freedom to be used
                model_solver.builder_and_solver.SetUpDofSet(model_solver.scheme,model_solver.model_part)
                #reorder the list of degrees of freedom to identify fixity and system size
                model_solver.builder_and_solver.SetUpSystem(model_solver.model_part)
                #allocate memory for the system and preallocate the structure of the matrix
                model_solver.builder_and_solver.ResizeAndInitializeVectors(model_solver.pA,model_solver.pDx,model_solver.pb,model_solver.model_part.Elements,model_solver.model_part.Conditions,model_solver.model_part.ProcessInfo)
                #updating references
                model_solver.A = (model_solver.pA).GetReference()
                model_solver.Dx = (model_solver.pDx).GetReference()
                model_solver.b = (model_solver.pb).GetReference()
            if(model_solver.SolutionStepIsInitialized == False):
                model_solver.builder_and_solver.InitializeSolutionStep(model_solver.model_part,model_solver.A,model_solver.Dx,model_solver.b)
                model_solver.scheme.InitializeSolutionStep(model_solver.model_part,model_solver.A,model_solver.Dx,model_solver.b)

    #######################################################################
    def ExecuteIteration(self,echo_level,CalculateNormDxFlag):
        for model_solver in self.model_solver_list:
            #reset system matrices and vectors prior to rebuild
            model_solver.space_utils.SetToZeroMatrix(model_solver.A)
            model_solver.space_utils.SetToZeroVector(model_solver.Dx)
            model_solver.space_utils.SetToZeroVector(model_solver.b)

            model_solver.scheme.InitializeNonLinIteration(model_solver.model_part,model_solver.A,model_solver.Dx,model_solver.b)

            #build and solve the problem
            model_solver.builder_and_solver.Build(model_solver.scheme,model_solver.model_part,model_solver.A,model_solver.b)
            model_solver.builder_and_solver.ApplyDirichletConditions(model_solver.scheme,model_solver.model_part,model_solver.A,model_solver.Dx,model_solver.b)
            model_solver.dof_util.ListDofs(model_solver.builder_and_solver.GetDofSet(),model_solver.builder_and_solver.GetEquationSystemSize())

        #provide data for the preconditioner and linear solver
        model_solver0 = self.model_solver_list[0]
        if model_solver0.linear_solver.AdditionalPhysicalDataIsNeeded():
            model_solver0.linear_solver.ProvideAdditionalData(model_solver0.A,model_solver0.Dx,model_solver0.b,model_solver0.builder_and_solver.GetDofSet(),model_solver0.model_part)
        model_solver0.linear_solver.Solve(model_solver0.A,model_solver0.Dx,model_solver0.b)

        ## update Dx for each level
        for i in range(1, len(self.gmg_level_list)):
            self.gmg_level_list[i-1].ApplyRestriction(self.model_solver_list[i-1].Dx, self.model_solver_list[i].Dx)

        #perform update
        normDx = []
        for model_solver in self.model_solver_list:
            model_solver.scheme.Update(model_solver.model_part,model_solver.builder_and_solver.GetDofSet(),model_solver.A,model_solver.Dx,model_solver.b)

            #move the mesh as needed
            if(model_solver.MoveMeshFlag == True):
                model_solver.scheme.MoveMesh(model_solver.model_part.Nodes)

            #to account for prescribed displacement, the displacement at prescribed nodes need to be updated
            for node in model_solver.model_part.Nodes:
                if node.IsFixed(DISPLACEMENT_X):
                    curr_disp = node.GetSolutionStepValue(DISPLACEMENT_X)
                    delta_disp = node.GetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_X)
                    node.SetSolutionStepValue(DISPLACEMENT_X, curr_disp + delta_disp)
                    node.SetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_X, 0.0) # set the prescribed displacement to zero to avoid update in the second step
                if node.IsFixed(DISPLACEMENT_Y):
                    curr_disp = node.GetSolutionStepValue(DISPLACEMENT_Y)
                    delta_disp = node.GetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_Y)
                    node.SetSolutionStepValue(DISPLACEMENT_Y, curr_disp + delta_disp)
                    node.SetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_Y, 0.0) # set the prescribed displacement to zero to avoid update in the second step
                if node.IsFixed(DISPLACEMENT_Z):
                    curr_disp = node.GetSolutionStepValue(DISPLACEMENT_Z)
                    delta_disp = node.GetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_Z)
                    node.SetSolutionStepValue(DISPLACEMENT_Z, curr_disp + delta_disp)
                    node.SetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_Z, 0.0) # set the prescribed displacement to zero to avoid update in the second step

            model_solver.scheme.FinalizeNonLinIteration(model_solver.model_part,model_solver.A,model_solver.Dx,model_solver.b)

            #calculate the norm of the "correction" Dx
            if(CalculateNormDxFlag == True):
                normDx.append(model_solver.space_utils.TwoNorm(model_solver.Dx))
            else:
                normDx.append(0.0)

        return normDx

    #######################################################################
    def FinalizeSolutionStep(self):
        for model_solver in self.model_solver_list:
            if(model_solver.CalculateReactionsFlag == True):
                model_solver.builder_and_solver.CalculateReactions(model_solver.scheme,model_solver.model_part,model_solver.A,model_solver.Dx,model_solver.b)

            #Finalisation of the solution step,
            model_solver.scheme.FinalizeSolutionStep(model_solver.model_part,model_solver.A,model_solver.Dx,model_solver.b)
            model_solver.builder_and_solver.FinalizeSolutionStep(model_solver.model_part,model_solver.A,model_solver.Dx,model_solver.b)
            model_solver.scheme.Clean()
            #reset flags for the next step
            model_solver.SolutionStepIsInitialized = False

    #######################################################################
    def SetEchoLevel(self,model_solver,level):
        model_solver.echo_level = level
        model_solver.builder_and_solver.SetEchoLevel(level)

    #######################################################################

