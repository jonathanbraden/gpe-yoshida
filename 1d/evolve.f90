#include "macros.h"

program Evolve_GPE
  use, intrinsic :: iso_c_binding
  use constants, only : dl, twopi
  use Simulation
  use Initial_Conditions
  use gaussianRandomField
  use Yoshida
  use Equations
  use Equations_imag
  implicit none
  
  type TimeParams
     real(dl) :: dt, dtout, alpha
     integer :: nstep, nout_step, out_step_size
  end type TimeParams
  
  !call run_trapped_background( 0.1, 10., 0.01, 128, 6)
  call run_imprinted_wave( 0.01, 0.2*twopi, 1.e-4, 0.5_dl, 2, 512, 6 )
  
contains
  
  !>@brief
  !> Run time evolution for trapped background evolution.
  !> The inputs are:
  !>  phi0  -> The initial phase displacement
  !>  g2    -> The dimensionless pairwise scattering cross section (get units)
  !>  nu_1d -> The dimensionless nu parameter (in the effective reduced 1D theory)
  !
  !>@todo - Implement automatic calculation of lSize using the Thomas-Fermi approximation for a given potential
  subroutine run_trapped_background(phi0, g2, nu_1d, nLat, ord)
    real(dl), intent(in) :: phi0, g2, nu_1d
    integer, intent(in) :: nLat, ord

    type(Lattice) :: mySim
    real(dl), parameter :: eps = 1.e-15
    real(dl) :: lSize
    integer :: nf
    real(dl) :: lperp, nu_dim
    real(dl) :: chemP, en  ! Remove these in future
    real(dl) :: period
    integer :: log_u, lat_u
    real(dl) :: dt
    integer :: out_size
    real(dl) :: alpha, dt_courant
    real(dl) :: dt_period
    integer :: bg_periods
    real(dl) :: dt_1d
    integer :: num_out
    integer :: i

    ! TBD: Automate the calculation of lSize using TF approximation
    lSize = 10._dl; nf = 2
    
    call create_lattice(mySim, nLat, lSize, nf)
    call set_model_parameters(g2, 0._dl, 0._dl, 0._dl, nf)
    call initialize_trap_potential(mySim, 32._dl, 3)

    ! Solve for the inhomogeneous background
    call imprint_gaussian(mySim,1._dl)
    call solve_background_w_grad_flow(mySim, eps, eps)

    ! Convert to parameters of 2D simulation
    chemP = chemical_potential(mySim)
    en = energy(mySim)
    lperp = g2 / ( chemP - en )  ! Change this to more direct integration
    nu_dim = g2*nu_1d / lperp

    !print*,"L_perp / L_SHO is ", lperp / sqrt(twopi) 
    
    ! Turn on nu, and rotate the condensate
    call set_model_parameters(g2,0.,nu_dim,0.,nf)
    ! Now that I've fixed up my chemical potential calculation, I think
    ! this is incorrect and the factor of 1/2 should be removed
    call set_chemical_potential(chemical_potential_full(mySim)/2._dl)
    call rotate_condensate(mySim, phi0, 2)

    ! With units used here, period is (omega_phi T)/2 * (lperp / g2) / sqrt(nu_1d).  Approximate omega_phi T = 2pi
    period = 0.5*twopi*lperp/sqrt(nu_1d)/g2
    dt_period = period / 32
    
    alpha = 8.; dt_courant = mySim%dx**2/alpha

    dt = dt_courant ! change to a minimum somewhere
    bg_periods = 5; out_size = int(period/dt) / 32
    print*,"dt_out = ",dt*out_size, dt*out_size/period
  
    num_out = 32*bg_periods
    print*,"t_final / period = ", dt*out_size*num_out/period
  
    dt_1d = dt*(g2/lperp)*sqrt(nu_1d)*2. ! Figure out what this is and re-express in 2d units

    log_u = 51; lat_u = 50
    open(unit=log_u,file='log.out')
    write(log_u,*) "# Time (sim units), Chemical Potential, Energy, Field Norm, Time (mass units)"
    write(log_u,*) 0., chemical_potential_full(mySim), energy(mySim), field_norm(mySim), 0.
  
    call write_lattice_data(mySim, lat_u)
    do i=1,num_out
       call step_lattice(mySim,dt,out_size,ord)
       write(log_u,*) dt*out_size*i, chemical_potential_full(mySim), energy(mySim), field_norm(mySim), dt_1d*out_size*i
       call write_lattice_data(mySim, lat_u)
    enddo
    close(log_u)
    close(lat_u)    
  end subroutine run_trapped_background

  subroutine run_imprinted_wave(nu, phi0, amp, floq_frac, wave_num, nLat, ord)
    real(dl), intent(in) :: nu, phi0, amp, floq_frac
    integer, intent(in) :: wave_num
    integer, intent(in) :: nLat, ord
    
    type(Lattice) :: mySim
    integer :: nf
    real(dl) :: lSize, lyap
    
    integer :: u_log, u_fld
    integer :: i
    
    ! Time-stepping params to factor out (can shove into a subroutine)
    real(dl) :: dt, dt_out
    real(dl) :: alpha, period
    integer :: out_size, out_steps
    integer :: num_periods, steps_per_period

    real(dl) :: w_tot, k_floq
    real(dl) :: lSize_heal
    
#ifndef PERIODIC
    stop "Error, trying to run preheating simulation without periodic BCs"
#endif
    
    nf = 2
    k_floq = 0.5_dl*phi0/sqrt(2._dl)
    k_floq = k_floq*floq_frac
    lyap = get_lyap(k_floq,phi0)
    lyap = lyap / (2.*sqrt(nu))                        ! Check this

    lSize = twopi/k_floq*wave_num
    lSize_heal = lSize / (2.*sqrt(nu))
    
    call create_lattice(mySim, nLat, lSize_heal, nf)
    call set_model_parameters(1._dl, 0._dl, nu, 0._dl, nf)
    call initialize_trap_potential(mySim, 0._dl, 1)

    call set_chemical_potential(1._dl-nu)

    call imprint_preheating_sine_wave(mySim, phi0, amp, wave_num, 1, 2)
    
    ! Work out time-stepping
    w_tot = nyquist_freq(mySim)
    alpha = 16.
    dt = (twopi/w_tot)/alpha
    
    period = 0.5_dl*twopi/sqrt(nu)   ! This is approximate.  Improve
    steps_per_period = 32
    out_size = int(period/dt)/steps_per_period  ! Fix this
    num_periods = 50
    out_steps = steps_per_period*num_periods 
    
    u_log = 51; u_fld = 50 ! Automate these with newunit
    open(unit=u_log, file='log.out')
    write(u_log,*) "# Time (sim units), Chemical Potential, Energy, Time (mass units), n_1, n_2"
    write(u_log,*) 0., chemical_potential_full(mySim), energy(mySim), 0., num_part(mySim,1), num_part(mySim,2)

    call write_lattice_data(mySim,u_fld)
    do i=1,out_steps  ! fix this
       call step_lattice(mySim, dt, out_size, ord)
       write(u_log,*) dt*out_size*i, chemical_potential_full(mySim), energy(mySim), 2._dl*sqrt(nu)*dt*out_size*i, num_part(mySim,1), num_part(mySim,2)
       call write_lattice_data(mySim, u_fld)
    enddo
    close(u_log); close(u_fld)
  end subroutine run_imprinted_wave

  !>@brief
  !> Return the angular frequency of the Nyquist mode for the total phase.
  !> Currently assumes symmetric condensate sitting at a minima
  real(dl) function nyquist_freq(this) result(om)
    type(Lattice), intent(in) :: this

    real(dl) :: k_nyq
    k_nyq = 0.5_dl*twopi/this%dx
    om = k_nyq*sqrt(1._dl + 0.25_dl*k_nyq**2)
  end function nyquist_freq

  real(dl) function get_lyap(k,phi0) result(lyap)
    real(dl), intent(in) :: k, phi0
    lyap = 0.25_dl*k*sqrt(phi0**2 - 4._dl*k**2)
  end function get_lyap
  
  subroutine run_preheating_symmetric(g2, nu, phi0, nLat, lSize)
    real(dl), intent(in) :: g2, nu, phi0
    integer, intent(in) :: nLat
    real(dl), intent(in) :: lSize
    
    integer :: nf
    type(Lattice) :: mySim

    nf = 2

    call create_lattice(mySim, nLat, lSize, nf)
    
    call set_model_parameters(g2, 0., nu, 0., nf)
    !call set_chemical_potential()  ! Fix this
    call rotate_condensate(mySim, phi0, 2)
    
  end subroutine run_preheating_symmetric
  
  subroutine summarise_ic(sim, g2)
    type(Lattice), intent(inout) :: sim
    real(dl), intent(in) :: g2
    
    real(dl) :: en, chemP
    real(dl) :: lperp
    
    chemP = chemical_potential(sim)
    en = energy(sim)
    print*,"chemical potential is ", chemP
    print*,"energy is ", en

    lperp = g2 / ( chemP - en )  ! Fix to have number of fields
    print*,"L_perp / L_perp,SHO is ", lperp / sqrt(twopi) 
    print*,"full one is ", chemical_potential_full(sim)
    print*,"field norm is ",field_norm(sim)
  end subroutine summarise_ic

end program Evolve_GPE
