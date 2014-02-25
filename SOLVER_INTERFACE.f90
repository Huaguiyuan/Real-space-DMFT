!###############################################################
!     PROGRAM  : SOLVER_INTERFACE
!     PURPOSE  : It is an interface to modules for the impurity solvers
!     AUTHORS  : Adriano Amaricci & Antonio Privitera
!###############################################################
module SOLVER_INTERFACE
  !GLOBAL VARIABLES IN COMMON TO ALL SOLVERS:
  USE SOLVER_INPUT_VARS
  !IPT:
  USE DMFT_IPT
  !ED:
  USE DMFT_ED
  implicit none
end module SOLVER_INTERFACE
