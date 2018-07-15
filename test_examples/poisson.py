import sys
import os
##################################################################
#importing Kratos modules
from KratosMultiphysics import *
from KratosMultiphysics.MultigridSolversApplication import *
kernel = Kernel()   #defining kernel

util = AMGUtils()
A = util.Poisson(5, 5)
print(A.GetReference())

