!==============================================================================
! Module: fieldline_tracer
! Description: Module for tracing magnetic field lines in cylindrical coordinates
!              using LSODE (Livermore Solver for Ordinary Differential Equations)
! Author: Refactored from original traceline.f90
! Date: 2025
!==============================================================================

module fieldline_tracer
  implicit none

  ! Module parameters
  integer, parameter :: dp = selected_real_kind(15, 307)  ! Double precision
  real(dp), parameter :: pi = 3.14159265358979323846264338327950288420_dp
  real(dp), parameter :: mu0 = 4.0_dp * pi * 1.0e-7_dp
  real(dp), parameter :: permeability_constant = 1.0e-7_dp  ! Permeability factor

  ! LSODE solver parameters
  integer, parameter :: lsode_neq = 2                  ! Number of equations
  integer, parameter :: lsode_mf = 10                  ! Adams-Moulton with functional iteration
  integer, parameter :: lsode_lrw = 20 + 16 * lsode_neq  ! Real workspace size
  integer, parameter :: lsode_liw = 20                 ! Integer workspace size
  real(dp), parameter :: lsode_rtol = 0.0_dp           ! Relative tolerance
  real(dp), parameter :: lsode_atol = 1.0e-13_dp       ! Absolute tolerance
  integer, parameter :: lsode_itol = 1                 ! Scalar RTOL and ATOL
  integer, parameter :: lsode_itask = 1                ! Normal task
  integer, parameter :: lsode_istate = 1               ! First call
  integer, parameter :: lsode_iopt = 0                 ! No optional inputs

  ! Module-level variables for storing coil data
  ! These are accessible to all subroutines in the module
  real(dp), allocatable, save :: coils_data_global(:, :)
  integer, save :: n_coils_global = 0

  ! Flag to indicate if coil data has been loaded
  logical, save :: coils_initialized = .false.
  
  ! Verbose flag for controlling print output
  logical, save :: verbose_mode = .false.

  private
  public :: trace_fieldlines, initialize_coils, cleanup_coils
  public :: get_magnetic_field  ! 添加新的公共接口
  public :: set_verbose  ! 添加设置verbose的接口

contains

  !============================================================================
  ! Subroutine: set_verbose
  ! Purpose: Set verbose mode for print output control
  ! Input:
  !   verbose - Logical flag to enable/disable verbose output
  !============================================================================
  subroutine set_verbose(verbose)
    implicit none
    logical, intent(in) :: verbose
    verbose_mode = verbose
  end subroutine set_verbose

  !============================================================================
  ! Subroutine: initialize_coils
  ! Purpose: Initialize coil data for use in field calculations
  ! Input:
  !   coils_data - Coil segment data (n_coils x 4): [x, y, z, current]
  !============================================================================
  subroutine initialize_coils(coils_data)
    implicit none
    real(dp), intent(in) :: coils_data(:, :)
    
    integer :: n_coils
    
    ! Get number of coils
    n_coils = size(coils_data, dim=1)
    
    ! Check input dimensions
    if (size(coils_data, dim=2) /= 4) then
      error stop "Error: coils_data must have 4 columns [x, y, z, current]"
    end if
    
    ! Clean up any existing data
    if (allocated(coils_data_global)) then
      deallocate(coils_data_global)
    end if
    
    ! Allocate and store coil data
    allocate(coils_data_global(n_coils, 4))
    coils_data_global = coils_data
    n_coils_global = n_coils
    coils_initialized = .true.
    
    if (verbose_mode) then
      print *, "Coil data initialized with", n_coils_global, "segments"
    end if
    
  end subroutine initialize_coils
  
  !============================================================================
  ! Subroutine: cleanup_coils
  ! Purpose: Clean up allocated coil data
  !============================================================================
  subroutine cleanup_coils()
    implicit none
    
    if (allocated(coils_data_global)) then
      deallocate(coils_data_global)
    end if
    
    n_coils_global = 0
    coils_initialized = .false.
    
    if (verbose_mode) then
      print *, "Coil data cleaned up"
    end if
    
  end subroutine cleanup_coils

  !============================================================================
  ! Subroutine: trace_fieldlines
  ! Purpose: Trace magnetic field lines using LSODE integration
  ! Input:
  !   nturn         - Number of toroidal turns to trace
  !   nphi          - Number of points per toroidal turn
  !   initial_rz    - Initial [R, Z] coordinates in cylindrical system
  ! Output:
  !   fieldline_data - Array of field line points (nturn*nphi x 4): [x, y, z, |B|]
  !============================================================================
  subroutine trace_fieldlines(nturn, nphi, initial_rz, fieldline_data)
    implicit none

    ! Arguments
    integer, intent(in) :: nturn, nphi
    real(dp), intent(in) :: initial_rz(2)
    real(dp), intent(out) :: fieldline_data(:, :)

    ! Local variables
    integer :: i, n_points
    integer :: lsode_istate_local
    real(dp) :: phi, phi_stop
    real(dp) :: rz_current(lsode_neq)
    real(dp), allocatable :: rwork(:), fieldpos(:), field(:)
    integer, allocatable :: iwork(:)

    ! Calculate total number of points
    n_points = nturn * nphi

    ! Check if coils are initialized
    if (.not. coils_initialized) then
      error stop "Error: Coil data not initialized. Call initialize_coils first."
    end if

    ! Validate output array
    if (size(fieldline_data, dim=2) < 4) then
      error stop "Error: fieldline_data must have at least 4 columns"
    end if

    ! Allocate workspace arrays
    allocate(rwork(lsode_lrw))
    allocate(iwork(lsode_liw))
    allocate(fieldpos(3))
    allocate(field(3))

    ! Initialize LSODE state
    lsode_istate_local = lsode_istate
    rz_current = initial_rz
    phi = 0.0_dp

    ! Trace field lines
    do i = 1, n_points
      ! Calculate phi_stop for this integration step
      phi_stop = phi + 2.0_dp * pi / real(nphi, kind=dp)

      ! Convert from cylindrical (R, Z) to Cartesian (x, y, z)
      fieldline_data(i, 1) = rz_current(1) * cos(phi)
      fieldline_data(i, 2) = rz_current(1) * sin(phi)
      fieldline_data(i, 3) = rz_current(2)

      ! Calculate magnetic field magnitude at this point
      fieldpos = fieldline_data(i, 1:3)
      call calculate_total_field(fieldpos, field)
      fieldline_data(i, 4) = sqrt(field(1)**2 + field(2)**2 + field(3)**2)

      ! Integrate along field line using LSODE
      call dlsode(fieldline_ode, lsode_neq, rz_current, phi, phi_stop, &
                  lsode_itol, lsode_rtol, lsode_atol, lsode_itask, &
                  lsode_istate_local, lsode_iopt, rwork, lsode_lrw, &
                  iwork, lsode_liw, jacobian_stub, lsode_mf)

      phi = phi_stop

      ! Check LSODE status
      if (lsode_istate_local < 0) then
        write(*, '(A, I0)') "Error: LSODE solver failed with ISTATE = ", lsode_istate_local
        exit
      end if
    end do

    ! Deallocate workspace
    deallocate(rwork, iwork, fieldpos, field)

  end subroutine trace_fieldlines


  !============================================================================
  ! Subroutine: fieldline_ode
  ! Purpose: Define the ODE system for field line tracing in cylindrical coords
  ! The system is: dR/dphi = R * Br / Bphi, dZ/dphi = R * Bz / Bphi
  !============================================================================
  subroutine fieldline_ode(neq, t, v, vdot)
    implicit none

    ! Arguments (LSODE required interface)
    integer, intent(in) :: neq
    real(dp), intent(in) :: t
    real(dp), intent(in) :: v(neq)
    real(dp), intent(out) :: vdot(neq)
    
    ! Local variables
    real(dp) :: r, z, phi
    real(dp) :: fieldpos(3), field(3)
    real(dp) :: br, bz, bphi
    real(dp) :: zero_threshold

    ! Extract coordinates
    r = v(1)
    z = v(2)
    phi = t

    ! Convert cylindrical to Cartesian for field calculation
    fieldpos(1) = r * cos(phi)
    fieldpos(2) = r * sin(phi)
    fieldpos(3) = z

    ! Calculate magnetic field components
    call calculate_total_field(fieldpos, field)

    ! Project field to cylindrical components
    br = field(1) * cos(phi) + field(2) * sin(phi)
    bz = field(3)
    bphi = -field(1) * sin(phi) + field(2) * cos(phi)

    ! Avoid division by zero
    zero_threshold = 1.0e-15_dp
    if (abs(bphi) < zero_threshold) then
      vdot(1) = 0.0_dp
      vdot(2) = 0.0_dp
    else
      vdot(1) = r * br / bphi
      vdot(2) = r * bz / bphi
    end if

  end subroutine fieldline_ode


  !============================================================================
  ! Subroutine: jacobian_stub
  ! Purpose: Jacobian matrix calculation (stub - not used in functional iteration)
  !============================================================================
  subroutine jacobian_stub(neq, t, y, ml, mu, pd, nrowpd)
    implicit none
    integer, intent(in) :: neq, ml, mu, nrowpd
    real(dp), intent(in) :: t, y(neq)
    real(dp), intent(out) :: pd(nrowpd, neq)
    ! This is a stub - Jacobian not used with MF=10
    pd = 0.0_dp
    return
  end subroutine jacobian_stub

  !============================================================================
  ! Subroutine: calculate_total_field
  ! Purpose: Calculate magnetic field from discrete coil segments using global data
  ! Input:
  !   fieldpos   - Position where field is calculated [x, y, z]
  ! Output:
  !   field      - Magnetic field components [Bx, By, Bz]
  !============================================================================
  subroutine calculate_total_field(fieldpos, field)
    implicit none

    ! Arguments
    real(dp), intent(in) :: fieldpos(3)
    real(dp), intent(out) :: field(3)

    ! Local variables
    integer :: i_coil
    real(dp) :: x, y, z
    real(dp) :: x0, y0, z0
    real(dp) :: dlx, dly, dlz
    real(dp) :: current
    real(dp) :: distance, distance_cubed
    real(dp) :: dx, dy, dz
    real(dp) :: field_contrib(3)

    ! Extract field point coordinates
    x = fieldpos(1)
    y = fieldpos(2)
    z = fieldpos(3)

    ! Initialize field
    field = 0.0_dp

    ! Calculate field from discrete coil segments
    do i_coil = 1, n_coils_global - 1
      ! Calculate segment midpoint
      x0 = 0.5_dp * (coils_data_global(i_coil, 1) + coils_data_global(i_coil + 1, 1))
      y0 = 0.5_dp * (coils_data_global(i_coil, 2) + coils_data_global(i_coil + 1, 2))
      z0 = 0.5_dp * (coils_data_global(i_coil, 3) + coils_data_global(i_coil + 1, 3))

      ! Calculate segment length vector
      dlx = coils_data_global(i_coil, 1) - coils_data_global(i_coil + 1, 1)
      dly = coils_data_global(i_coil, 2) - coils_data_global(i_coil + 1, 2)
      dlz = coils_data_global(i_coil, 3) - coils_data_global(i_coil + 1, 3)

      ! Get segment current
      current = coils_data_global(i_coil, 4)

      ! Calculate distance from field point to segment midpoint
      dx = x - x0
      dy = y - y0
      dz = z - z0
      distance = sqrt(dx**2 + dy**2 + dz**2)
      distance_cubed = distance**3

      ! Avoid division by zero
      if (distance_cubed < 1.0e-20_dp) cycle

      ! Calculate Biot-Savart contribution: B = (μ₀/4π) * I * (dl × r) / r³
      field_contrib(1) = current * (dly * dz - dlz * dy) / distance_cubed
      field_contrib(2) = current * (dlz * dx - dlx * dz) / distance_cubed
      field_contrib(3) = current * (dlx * dy - dly * dx) / distance_cubed
      
      ! Accumulate field
      field = field + field_contrib
    end do

    ! Apply permeability factor (μ₀/4π equivalent)
    field = field * permeability_constant
    
  end subroutine calculate_total_field

  !============================================================================
  ! Subroutine: get_magnetic_field
  ! Purpose: Public wrapper for calculate_total_field
  ! Input:
  !   fieldpos - Position [x, y, z]
  ! Output:
  !   field    - Magnetic field [Bx, By, Bz]
  !============================================================================
  subroutine get_magnetic_field(fieldpos, field)
    implicit none
    real(dp), intent(in) :: fieldpos(3)
    real(dp), intent(out) :: field(3)
    
    if (.not. coils_initialized) then
      error stop "Error: Coil data not initialized. Call initialize_coils first."
    end if
    
    call calculate_total_field(fieldpos, field)
  end subroutine get_magnetic_field

end module fieldline_tracer