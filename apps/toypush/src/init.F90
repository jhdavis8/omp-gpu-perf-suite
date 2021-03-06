module initmodule

  use params
  
  implicit none

contains
  
  !> Initialize the particles, grid, E-field, etc. This stuff could/should be read from an input
  !> file in a future release.
  !<
  function init(prt_mass, prt_charge, prt_rpz, prt_mu, prt_rho_par, prt_isAllocated) result(err)

    use rk4, only: rk4_init
    use particle, only : particle_init, particle_updatephase
#ifdef USEIO
    use particleIO, only: particleio_write
#endif
    use grid_module, only : get_coords, get_nodes, grid_init, grid_init_coords, grid_init_efield

    implicit none

    double precision, allocatable, dimension(:)   :: prt_mass
    double precision, allocatable, dimension(:)   :: prt_charge
    double precision, allocatable, dimension(:,:) :: prt_rpz
    double precision, allocatable, dimension(:)   :: prt_mu
    double precision, allocatable, dimension(:)   :: prt_rho_par
    logical                                       :: prt_isAllocated

    integer :: i,ir,ith,nr,nth
    integer :: err

    double precision, allocatable, dimension(:,:) :: y

    double precision, dimension(2,params_nnode) :: coords
    integer, dimension(3,params_ntri) :: tri
    double precision, dimension(3,params_nnode) :: efield

    double precision :: v

    allocate(y(params_nprt,4))
#ifdef _OPENMP
    !$omp target data map(alloc: y)
#elif _OPENACC
    !$acc data create(y)
#endif

    nth = 32
    nr = params_nprt / nth
    err = particle_init(prt_mass, prt_charge, prt_rpz, prt_mu, prt_rho_par, prt_isAllocated, params_nprt)

    v = 1D4
#ifdef _OPENMP
    !$omp target teams distribute collapse(2) 
#elif _OPENACC
    !$acc parallel loop collapse(2)
#endif
    do ir = 1,nr
       do ith = 1,nth
          i = (ith - 1) * nr + ir
          y(i,1) = rmin + (rmax-rmin) / 2D0 + &
               cos(dble(ith-1)/dble(nth-1) * twopi) * dble(ir)/dble(nr) * 0.5D0
          y(i,2) = 0.0D0
          y(i,3) = sin(dble(ith-1)/dble(nth-1) * twopi) * dble(ir)/dble(nr) * 0.5D0
          y(i,4) = protonmass / unitcharge * v
       end do
    end do

    err = particle_updatePhase(prt_mass, prt_charge, prt_rpz, prt_mu, prt_rho_par, prt_isAllocated, y,params_nprt,1)
#ifdef _OPENMP
    !$omp target update from(prt_rpz,prt_rho_par)
#elif _OPENACC
    !$acc update host(prt_rpz,prt_rho_par)
#endif

    prt_mass = 2.0D0 * protonmass
    prt_charge = unitcharge
    prt_mu = 0.5D0 * protonmass * v ** 2
    
    err = rk4_init()

    err = grid_init(params_nnode,params_ntri)
    
    coords(:,1) = [rmin,       zmin      ]
    coords(:,2) = [2D0 * rmax, zmin      ]
    coords(:,3) = [rmin,       2D0 * zmax]

    tri(:,1) = [1,2,3]

#ifdef MULTIPLEELEMENTS
    coords(:,1) = [rmin, zmin]
    coords(:,2) = [rmax, zmin]
    coords(:,3) = [rmin, zmax]
    coords(:,4) = [rmax, zmax]
    tri(:,2) = [4,2,3]
    efield(:,4) = [5D-1, 0D0, 0D0] ! node 3
#endif

    
    err = grid_init_coords(coords,tri)

    !              ER   Ep   Ez
    efield(:,1) = [ 0D0, 0D0, 0D0] ! node 1
    efield(:,2) = [ 1D0, 0D0, 0D0] ! node 2
    efield(:,3) = [-1D0, 0D0, 0D0] ! node 3

    efield = efield * 1D0
    
    err = grid_init_efield(efield)
    
#ifdef USEIO
    err = particleio_write(prt_mass, prt_charge, prt_rpz, prt_mu, prt_rho_par, prt_isAllocated, 'inistate.dat')
#endif
    deallocate(y)
#ifdef _OPENMP
    !$omp end target data 
#elif _OPENACC
    !$acc end data
#endif
  end function init

end module initmodule
