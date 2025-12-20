!==============================================================================
! Wrapper subroutine: tracefieldline
! Purpose: Legacy interface to trace_fieldlines from fieldline_tracer module
!          Uses global variables from constant module
!==============================================================================
subroutine tracefieldline()
  use constant, only: n_d, nt, np, RZ_i, TF_current, fieldline, discrete_coil
  use fieldline_tracer, only: trace_fieldlines
  implicit none

  ! Call the new module subroutine
  call trace_fieldlines(nt * np, RZ_i, discrete_coil(1:n_d, 1:4), TF_current, fieldline)

end subroutine tracefieldline









