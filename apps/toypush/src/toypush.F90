!> This program is a simple implementation of a RK4 particle push routine to test theoretical peak performance
!> @author T. Koskela
!> @date Aug 17 2016
!<
program toypush

  use params
  use rk4, only: rk4_push, y, dy, ytmp, dyt, dym, y2, jacb, efield, bfield
! use particle, only: particle_data
  use initmodule, only : init
use grid_module, only : grid_efield, grid_mapping, grid_node, grid_tri

  implicit none  
  
  integer :: err
   double precision, allocatable, dimension(:)   :: prt_mass
   double precision, allocatable, dimension(:)   :: prt_charge
   double precision, allocatable, dimension(:,:) :: prt_rpz
   double precision, allocatable, dimension(:)   :: prt_mu
   double precision, allocatable, dimension(:)   :: prt_rho_par
   logical                                       :: prt_isAllocated

  integer :: it,iblock
  integer :: nblock
  
  integer :: pid

  integer :: num_procs, my_id, ith

  double precision :: t1,t2

  my_id = 0
  num_procs = 1

  params_nprt = params_nprtPerRank * num_procs

  if(my_id .eq. 0) write(*,*) 'program toypush started'
  if(my_id .eq. 0) write(*,*) 'veclen = ',veclen
  if(my_id .eq. 0) write(*,*)
  
  if(my_id .eq. 0) write(*,*) 'initializing particles with ',params_nprtperrank,'particles per rank.'
  if(my_id .eq. 0) write(*,*) 'initializing particles with ',params_nprt,'total particles.'
  if(my_id .eq. 0) write(*,*) 'initializing grid      with ',params_nnode,'nodes.'
  if(my_id .eq. 0) write(*,*) 'initializing grid      with ',params_ntri,'triangles.'
  err = init(prt_mass, prt_charge, prt_rpz, prt_mu, prt_rho_par, prt_isAllocated)
  if(err .ne. 0) stop
  if(my_id .eq. 0) write(*,*) 'done initialising'
  if(my_id .eq. 0) write(*,*)
  
  nblock = params_nprt / veclen  
  if(my_id .eq. 0) write(*,*) 'pushing particles in ',nblock,' blocks'

  call cpu_time(t1)
  
  !pid = 88
  !open(unit=15,file='orbit.dat',action='write')
#ifdef _OPENMP
  !$omp target update to(prt_rpz, prt_rho_par, y, bfield, efield, prt_charge, prt_mass, prt_mu, jacb, ytmp, dy, dyt, dym, y2, grid_mapping, grid_efield, grid_tri)
#elif _OPENACC
  !$acc update device(prt_rpz, prt_rho_par, y, bfield, efield, prt_charge, prt_mass, prt_mu, jacb, ytmp, dy, dyt, dym, y2, grid_mapping, grid_efield, grid_tri)
#endif
  do iblock = 1, nblock

        write(*,*) "iblock=", iblock, " of nblock=", nblock
        do it = 1, nt
           err = rk4_push(prt_mass, prt_charge, prt_rpz, prt_mu, prt_rho_par, prt_isAllocated, iblock)
           
        end do
  end do
  !close(15)
#ifdef _OPENMP
  !$omp target update from(prt_rpz, prt_rho_par, y2, y, bfield, efield, prt_charge, prt_mass, prt_mu, jacb, dy, dyt, dym, ytmp, grid_mapping, grid_efield, grid_tri)
#elif _OPENACC
  !$acc update host(prt_rpz, prt_rho_par, y2, y, bfield, efield, prt_charge, prt_mass, prt_mu, jacb, dy, dyt, dym, ytmp, grid_mapping, grid_efield, grid_tri)
#endif

  call cpu_time(t2)

  if(my_id .eq. 0) write(*,*) 'done pushing'
  if(my_id .eq. 0) write(*,*) 'spent ',t2-t1,'s'
  if(my_id .eq. 0) write(*,*)

  !call mpi_gatherv
  
  if(my_id .eq. 0)  write(*,*) 'finalizing'
  err = finalize(prt_mass, prt_charge, prt_rpz, prt_mu, prt_rho_par, prt_isAllocated)
  if(my_id .eq. 0) write(*,*) 'done finalizing'
  if(my_id .eq. 0) write(*,*)
  
  if(my_id .eq. 0) write(*,*) 'program toypush successfully completed!'
  
contains
  
  function finalize(prt_mass, prt_charge, prt_rpz, prt_mu, prt_rho_par, prt_isAllocated) result(err)

    use rk4, only: rk4_deallocate
    use particle, only : particle_deallocate
#ifdef USEIO
    use particleio, only : particleio_write
#endif
    
    implicit none

   double precision, allocatable, dimension(:)   :: prt_mass
   double precision, allocatable, dimension(:)   :: prt_charge
   double precision, allocatable, dimension(:,:) :: prt_rpz
   double precision, allocatable, dimension(:)   :: prt_mu
   double precision, allocatable, dimension(:)   :: prt_rho_par
   logical                                       :: prt_isAllocated

    integer :: err

#ifdef USEIO
    err = particleio_write(prt_mass, prt_charge, prt_rpz, prt_mu, prt_rho_par, prt_isAllocated, 'endstate.dat')
#endif
    
    err = particle_deallocate(prt_mass, prt_charge, prt_rpz, prt_mu, prt_rho_par, prt_isAllocated)
    err = rk4_deallocate()
    
    deallocate(grid_efield)
    deallocate(grid_mapping)
    deallocate(grid_node)
    deallocate(grid_tri)
#ifdef _OPENMP
    !$omp target exit data map(release: grid_efield, grid_mapping, grid_tri)
#elif _OPENACC
    !$acc exit data delete(grid_efield, grid_mapping, grid_tri)
#endif

  end function finalize
  
end program toypush
