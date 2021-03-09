!> This module stores relevant particle data
!<
module particle

  use params
  
  implicit none

  integer, parameter :: ALREADYALLOCATED = params_particleModuleErrorId + 1
  integer, parameter :: NOTALLOCATED     = params_particleModuleErrorId + 2
  
  double precision, allocatable, dimension(:) :: prt_mass
  double precision, allocatable, dimension(:) :: prt_charge
  double precision, allocatable, dimension(:,:) :: prt_rpz
  double precision, allocatable, dimension(:) :: prt_mu
  double precision, allocatable, dimension(:) :: prt_rho_par

  logical :: prt_isAllocated
     
contains

  !> Updates the phase variables rpz and rho_par of the particle
  !<
  function particle_updatePhase(prt_mass, prt_charge, prt_rpz, prt_mu, prt_rho_par, prt_isAllocated, newPhase,veclen,iblock) result(err)

    integer, intent(in) :: veclen !> block size
    integer, intent(in) :: iblock !> block id  
    double precision, allocatable, dimension(:), intent(inout) :: prt_mass
    double precision, allocatable, dimension(:), intent(inout) :: prt_charge
    double precision, allocatable, dimension(:,:), intent(inout) :: prt_rpz
    double precision, allocatable, dimension(:), intent(inout) :: prt_mu
    double precision, allocatable, dimension(:), intent(inout) :: prt_rho_par
    logical, intent(inout) :: prt_isAllocated

    double precision, allocatable, dimension(:,:), intent(in) :: newPhase

    integer :: err
    integer :: iv,iglob

    err = 0
    
#ifdef _OPENMP
    !$omp target teams distribute 
#elif _OPENACC 
    !$acc parallel loop
#endif
    do iv = 1,veclen
       iglob = (iblock - 1) * veclen + iv
       prt_rpz(iglob,1)   = newPhase(iv,1)
       prt_rpz(iglob,2)   = newPhase(iv,2)
       prt_rpz(iglob,3)   = newPhase(iv,3)
       prt_rho_par(iglob) = newPhase(iv,4)
    end do
    
  end function particle_updatePhase

  !> Returns the phase variables rpz and rho_par of the particle
  !<
  function particle_getPhase(prt_mass, prt_charge, prt_rpz, prt_mu, prt_rho_par, prt_isAllocated, phase,veclen,iblock) result(err)


    integer, intent(in) :: veclen !> block size
    integer, intent(in) :: iblock !> block id  
    double precision, allocatable, dimension(:), intent(in) :: prt_mass
    double precision, allocatable, dimension(:), intent(in) :: prt_charge
    double precision, allocatable, dimension(:,:), intent(in) :: prt_rpz
    double precision, allocatable, dimension(:), intent(in) :: prt_mu
    double precision, allocatable, dimension(:), intent(in) :: prt_rho_par
    logical, intent(in) :: prt_isAllocated

    double precision, dimension(veclen,4), intent(out) :: phase

    integer :: err
    integer :: iv, iglob

    err = 0
    

#ifdef _OPENMP
    !$omp target teams distribute 
#elif _OPENACC 
    !$acc parallel loop
#endif
    do iv = 1,veclen
       iglob = (iblock - 1) * veclen + iv
       phase(iv,1)   = prt_rpz(iglob,1)
       phase(iv,2)   = prt_rpz(iglob,2)
       phase(iv,3)   = prt_rpz(iglob,3)
       phase(iv,4)   = prt_rho_par(iglob)
    end do
    
  end function particle_getPhase
  
  function particle_init(prt_mass, prt_charge, prt_rpz, prt_mu, prt_rho_par, prt_isAllocated,arraydim) result(err)

    integer, intent(in) :: arraydim !> size of arrays
    double precision, allocatable, dimension(:), intent(inout) :: prt_mass
    double precision, allocatable, dimension(:), intent(inout) :: prt_charge
    double precision, allocatable, dimension(:,:), intent(inout) :: prt_rpz
    double precision, allocatable, dimension(:), intent(inout) :: prt_mu
    double precision, allocatable, dimension(:), intent(inout) :: prt_rho_par
    logical, intent(inout) :: prt_isAllocated

    integer :: err

    if(prt_isallocated) then
       err = ALREADYALLOCATED
       return
    end if

    allocate(prt_mass(arraydim))
    allocate(prt_charge(arraydim))
    allocate(prt_mu(arraydim))

    allocate(prt_rpz(arraydim,3))
    allocate(prt_rho_par(arraydim))

#ifdef _OPENMP
    !$omp target enter data map(alloc: prt_mass, prt_charge, prt_mu, prt_rpz, prt_rho_par)
#elif _OPENACC
    !$acc enter data create(prt_mass, prt_charge, prt_mu, prt_rpz, prt_rho_par)
#endif

    prt_isAllocated = .true.
    err = 0
    
  end function particle_init

  function particle_deallocate(prt_mass, prt_charge, prt_rpz, prt_mu, prt_rho_par, prt_isAllocated) result(err)

    double precision, allocatable, dimension(:), intent(inout) :: prt_mass
    double precision, allocatable, dimension(:), intent(inout) :: prt_charge
    double precision, allocatable, dimension(:,:), intent(inout) :: prt_rpz
    double precision, allocatable, dimension(:), intent(inout) :: prt_mu
    double precision, allocatable, dimension(:), intent(inout) :: prt_rho_par
    logical, intent(inout) :: prt_isAllocated

    integer :: err

    if(.not. prt_isallocated) then
       err = NOTALLOCATED
       return
    end if

    deallocate(prt_mass)
    deallocate(prt_charge)
    deallocate(prt_mu)

    deallocate(prt_rpz)
    deallocate(prt_rho_par)

#ifdef _OPENMP
    !$omp target exit data map(release: prt_mass, prt_charge, prt_mu, prt_rpz, prt_rho_par)
#elif _OPENACC
    !$acc exit data delete(prt_mass, prt_charge, prt_mu, prt_rpz, prt_rho_par)
#endif

    prt_isallocated = .false.
    err = 0
    
  end function particle_deallocate
  
end module particle
