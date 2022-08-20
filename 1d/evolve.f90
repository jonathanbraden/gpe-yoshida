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

  type(Lattice) :: mySim
  integer :: i, nf
  real(dl) :: dx, dt, alpha
  real(dl), dimension(:), allocatable :: psi_i
  integer :: ord

  ! This is all ugly temporary stuff
  real(dl) :: g2, nu_dim, nu_1d
  real(dl) :: lperp, en, chemP
  real(dl) :: dt_1d
  integer :: out_size, num_out
  real(dl) :: period
  integer :: bg_periods
  
  g2 = 10.; nu_1d = 0.01; nf = 2
  ord = 6
  
  call create_lattice(mySim,128,10._dl,nf)
  call set_model_parameters(g2,0.,0.,0._dl,nf)  
  call initialize_trap_potential(mySim, 32._dl,3)

  call imprint_gaussian(mySim,1._dl)
  call solve_background_w_grad_flow(mySim, 1.e-15, 1.e-15)

  chemP = chemical_potential(mySim)
  en = energy(mySim)
  lperp = g2 / ( chemP - en )  ! Fix to have number of fields
  nu_dim = g2*nu_1d / lperp

  ! With the normalization I'm using, I don't get the correct
  ! time-evolution.  But I think this is due to my field normalization
  !  Check what happens if I divide by the field normalization
  ! No, problem is I don't have nu in the def. of chemical potential
  !
  ! Did I fix this?

  call set_model_parameters(g2,0.,nu_dim,0.,nf)  
  call set_chemical_potential(chemical_potential_full(mySim)/2.)
  call rotate_condensate(mySim,0.1,2)

  ! With this def, period is (omega_phi T)/2 * (lperp / g2) / sqrt(nu_1d)
  period = 0.5*twopi*lperp/sqrt(nu_1d)/g2
  print*, "Approximate period is (assuming 2pi for period in mass units) ", period

  dx = mySim%dx; alpha = 8.; dt = dx**2/alpha
  bg_periods = 5; out_size = int(period/dt) / 32
  print*,"dt_out = ",dt*out_size, dt*out_size/period
  
  num_out = 32*bg_periods
  print*,"t_final / period = ", dt*out_size*num_out/period
  
  dt_1d = dt*(g2/lperp)*sqrt(nu_1d)*2. ! Figure out what this is and re-express in 2d units
  
  open(unit=51,file='log.out')
  write(51,*) "# Time (sim units), Chemical Potential, Energy, Field Norm, Time (mass units)"
  write(51,*) 0., chemical_potential_full(mySim), energy(mySim), field_norm(mySim), 0.
  
  call write_lattice_data(mySim,50)
  do i=1,num_out
     call step_lattice(mySim,dt,out_size,ord)
     write(51,*) dt*out_size*i, chemical_potential_full(mySim), energy(mySim), field_norm(mySim), dt_1d*out_size*i
     call write_lattice_data(mySim,50)
  enddo
  close(51)
  
contains

  !>@brief
  !> Run time evolution for trapped background evolution.
  !> The inputs are:
  !>  phi0  -> The initial phase displacement
  !>  g2    -> The dimensionless pairwise scattering cross section (get units)
  !>  nu_1d -> The dimensionless nu parameter (in the effective reduced 1D theory)
  !
  !>@todo - Implement automatic calculation of lSize using the Thomas-Fermi approximation for a given potential
  subroutine run_trapped_background(phi0, g2, nu_1d, nLat)
    real(dl), intent(in) :: phi0, g2, nu_1d
    integer, intent(in) :: nLat

    real(dl) :: lSize
    integer :: nf
    type(Lattice) :: mySim
    real(dl), parameter :: eps = 1.e-15
    real(dl) :: lperp, nu_dim
    real(dl) :: chemP, en  ! Remove these in future
    
    lSize = 10._dl; nf = 2

    ! TBD: Automate the calculation of lSize using TF approximation
    
    call create_lattice(mySim, nLat, lSize, nf)
    call set_model_parameters(g2, 0._dl, 0._dl, 0._dl, nf)
    call initialize_trap_potential(mySim, 32._dl,3)

    ! Solve for the inhomogeneous background
    call imprint_gaussian(mySim,1._dl)
    call solve_background_w_grad_flow(mySim, eps, eps)

    ! Convert to parameters of 2D simulation
    chemP = chemical_potential(mySim)
    en = energy(mySim)
    lperp = g2 / ( chemP - en ) 
    nu_dim = g2*nu_1d / lperp

    ! This comment is here to remind be of the normalization.
    !print*,"L_perp / L_SHO is ", lperp / sqrt(twopi) 
    
    ! Turn on nu, and rotate the condensate
    call set_model_parameters(g2,0.,nu_dim,0.,nf)
    call set_chemical_potential(chemical_potential_full(mySim)/2._dl)  ! Check my norm here
    call rotate_condensate(mySim, phi0, 2)

    ! Work out number of required time-steps, etc.
    
  end subroutine run_trapped_background

  subroutine run_imprinted_wave(amp, wave_num, nLat)
    real(dl), intent(in) :: amp
    integer, intent(in) :: wave_num
    integer, intent(in) :: nLat
    
    type(Lattice) :: mySim
    integer :: nf
    real(dl) :: lSize, nu
    integer :: u_log

    ! Time-stepping params to factor out
    real(dl) :: dt, dt_out
    integer :: out_size
    integer :: i
    
    nf = 2; lSize = 10._dl
    nu = 1.e-2
    
    call create_lattice(mySim, nLat, lSize, nf)
    call set_model_parameters(1._dl, 0._dl, 0._dl, nu, nf)
    call initialize_trap_potential(mySim, 0._dl, 1)

    ! Work on this
    chemP = 1._dl-nu
    call set_chemical_potential(chemP)

    call imprint_preheating_sine_wave(mySim, 0.125*twopi, 1.e-6, 2, 1, 2)

    ! Work out time-stepping, etc.

    u_log = 51
    open(unit=u_log, file='log.out')
    write(u_log,*) "# Time (sim units), Chemical Potential, Energy, Field Norm"
    write(u_log,*) "0., chemical_potential_full(mySim), energy(mySim), field_norm(mySim)"

    call write_lattice_data(mySim,50)
    do i=1,num_out
       call step_lattice(mySim, dt, out_size, ord)
       write(u_log,*) dt*out_size*i, chemical_potential_full(mySim), energy(mySim), field_norm(mySim)
       call write_lattice_data(mySim, 50)
    enddo
    close(u_log)
  end subroutine run_imprinted_wave

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
  
  subroutine summarise_ic(sim)
    type(Lattice), intent(inout) :: sim
    real(dl) :: en, chemP

    chemP = chemical_potential(mySim)
    en = energy(mySim)
    print*,"chemical potential is ", chemP
    print*,"energy is ", en

    lperp = g2 / ( chemP - en )  ! Fix to have number of fields
    print*,"L_perp / L_perp,SHO is ", lperp / sqrt(twopi) 
    print*,"full one is ", chemical_potential_full(mySim)
    print*,"field norm is ",field_norm(mySim)
  end subroutine summarise_ic
  
end program Evolve_GPE
