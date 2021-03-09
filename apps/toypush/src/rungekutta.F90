module rk4

  use params
  
  implicit none  
  
  double precision, allocatable, dimension(:,:) :: ytmp,dy,dyt,dym,y,y2
  double precision, allocatable, dimension(:,:) :: Efield !> Electric field vector (Er,Ephi,Ez)
  double precision, allocatable, dimension(:,:) :: Bfield !> Magnetic field vector (Br,Bphi,Bz)
  !> Magnetic field jacobian
  !> |  dBRdR,   dBRdphi,   dBRdz  |
  !> | dBphidR, dBphidphi, dBphidz |
  !> |  dBzdR,   dBzdphi,   dBzdz  |
  !<
  double precision, allocatable, dimension(:,:,:) :: jacB
  
contains
  
  function rk4_push(prt_mass, prt_charge, prt_rpz, prt_mu, prt_rho_par, prt_isAllocated, iblock) result(err)

    use eom, only : eom_eval
    use interpolate, only : b_interpol_analytic, e_interpol_tri
    use particle, only : particle_getphase, particle_updatephase
    use search_module, only : search_tr_vec
    use grid_module, only : grid_ntri, grid_mapping, grid_efield, grid_tri
    
    implicit none

    double precision, allocatable, dimension(:)   :: prt_mass
    double precision, allocatable, dimension(:)   :: prt_charge
    double precision, allocatable, dimension(:,:) :: prt_rpz
    double precision, allocatable, dimension(:)   :: prt_mu
    double precision, allocatable, dimension(:)   :: prt_rho_par
    logical                                       :: prt_isAllocated

    integer, intent(in) :: iblock
    
    double precision :: hdt, dt6
    integer :: err

    integer :: iv,iy
    integer, dimension(veclen) :: itri

#ifdef _OPENMP
    !$omp target data map(alloc: itri)
#elif _OPENACC
    !$acc data create(itri)
#endif

    hdt = dt * 0.5D0
    dt6 = dt / 6.0D0

    err = particle_getPhase(prt_mass, prt_charge, prt_rpz, prt_mu, prt_rho_par, prt_isAllocated,y,veclen,iblock)

#ifndef MULTIPLEELEMENTS
    itri = 1
#ifdef _OPENMP
    !$omp target update to(itri)
#elif _OPENACC
    !$acc update device(itri)
#endif
#endif

#ifdef DEBUG
    err = check_bounds(y)
    if(err .eq. 1) stop
#endif
    
    ! get derivs with existing E-field
#ifdef MULTIPLEELEMENTS
    err = search_tr_vec(y,itri)
#endif

    err = e_interpol_tri(y,itri,efield)

    err = b_interpol_analytic(y,bfield,jacb)

    err = eom_eval(y   ,bfield,jacb,efield,dt,dy ,prt_mu,prt_charge,prt_mass)

#ifdef _OPENMP
    !$omp target teams distribute collapse(2)
#elif _OPENACC
    !$acc parallel loop collapse(2)
#endif
    do iy = 1,4
       do iv = 1,veclen
          ytmp(iv,iy) = y(iv,iy) + hdt * dy(iv,iy)
       end do
    end do

#ifdef DEBUG
    err = check_bounds(ytmp)
    if(err .eq. 1) stop
#endif
#ifdef MULTIPLEELEMENTS
    err = search_tr_vec(ytmp,itri)
#endif

    err = e_interpol_tri(ytmp,itri,efield)

    err = b_interpol_analytic(ytmp,bfield,jacb)

    err = eom_eval(ytmp,bfield,jacb,efield,dt,dyt,prt_mu,prt_charge,prt_mass)

#ifdef _OPENMP
    !$omp target teams distribute collapse(2)
#elif _OPENACC
    !$acc parallel loop collapse(2)
#endif
    do iy = 1,4
       do iv = 1,veclen
          ytmp(iv,iy) = y(iv,iy) + hdt * dyt(iv,iy)
       end do
    end do

#ifdef DEBUG
    err = check_bounds(ytmp)
    if(err .eq. 1) stop
#endif

#ifdef MULTIPLEELEMENTS
    err = search_tr_vec(ytmp,itri)
#endif
    err = e_interpol_tri(ytmp,itri,efield)
    err = b_interpol_analytic(ytmp,bfield,jacb)

    err = eom_eval(ytmp,bfield,jacb,efield,dt,dym,prt_mu,prt_charge,prt_mass)
    
#ifdef _OPENMP
    !$omp target teams distribute collapse(2)
#elif _OPENACC
    !$acc parallel loop collapse(2)
#endif
    do iy = 1,4
       do iv = 1,veclen
          ytmp(iv,iy) = y(iv,iy) + dt * dym(iv,iy)
       end do
    end do

#ifdef DEBUG
    err = check_bounds(ytmp)
    if(err .eq. 1) stop
#endif

#ifdef _OPENMP
    !$omp target teams distribute collapse(2)
#elif _OPENACC
    !$acc parallel loop collapse(2)
#endif
    do iy = 1,4
       do iv = 1,veclen
          dym(iv,iy) = dyt(iv,iy) + dym(iv,iy)
       end do
    end do

#ifdef MULTIPLEELEMENTS
    err = search_tr_vec(ytmp,itri)
#endif
    err = e_interpol_tri(ytmp,itri,efield)

    err = b_interpol_analytic(ytmp,bfield,jacb)

    err = eom_eval(ytmp,bfield,jacb,efield,dt,dyt,prt_mu,prt_charge,prt_mass)
    
    ! Obtain new_phase
#ifdef _OPENMP
    !$omp target teams distribute collapse(2)
#elif _OPENACC
    !$acc parallel loop collapse(2)
#endif
    do iy = 1,4
       do iv = 1,veclen
          y2(iv,iy) = y(iv,iy) + dt6 * ( dy(iv,iy) + dyt(iv,iy) + 2.0D0*dym(iv,iy) )
       end do
    end do
#ifdef DEBUG
    err = check_bounds(y2)
    if(err .eq. 1) stop
#endif

    err = particle_updatePhase(prt_mass, prt_charge, prt_rpz, prt_mu, prt_rho_par, prt_isAllocated,y2,veclen,iblock)

#ifdef _OPENMP
    !$omp end target data
#elif _OPENACC
    !$acc end data
#endif

  end function rk4_push

  function rk4_init() result(err)

    implicit none

    integer :: err
    
    allocate(ytmp(veclen,4)) ! used by rk4_push, goes into eom_eval
    allocate(y(   veclen,4)) ! used by rk4_push, goes into eom_eval
    allocate(y2(  veclen,4)) ! used by rk4_push only
    allocate(dy(  veclen,4)) ! used by eom_eval
    allocate(dyt( veclen,4)) ! used by eom_eval
    allocate(dym( veclen,4)) ! used by eom_eval

    allocate(efield(veclen,3))   ! This is used by e_interpol_tri
    allocate(bfield(veclen,3))   ! This is used by b_interpol_analytic
    allocate(jacb(veclen,3,3))   ! This is used by b_interpol_analytic

#ifdef _OPENMP
    !$omp target enter data map(alloc: ytmp, y, y2, dy, dyt, dym, efield, bfield, jacb)
#elif _OPENACC
    !$acc enter data create(ytmp, y, y2, dy, dyt, dym, efield, bfield, jacb)
#endif

    err = 0
    
  end function rk4_init

  function rk4_deallocate() result(err)

    implicit none

    integer :: err
    
    deallocate(ytmp)
    deallocate(y)
    deallocate(y2)
    deallocate(dy)
    deallocate(dyt)
    deallocate(dym)

    deallocate(efield)
    deallocate(bfield)
    deallocate(jacb)

#ifdef _OPENMP
    !$omp target exit data map(release: ytmp, y, y2, dy, dyt, dym, efield, bfield, jacb)
#elif _OPENACC
    !$acc exit data delete(ytmp, y, y2, dy, dyt, dym, efield, bfield, jacb)
#endif

    err = 0

  end function rk4_deallocate

  function check_bounds(y) result(err)

    use params, only : rmin,rmax,zmin,zmax,veclen

    integer :: err
    double precision, intent(in), dimension(veclen,4) :: y
    
    integer :: iv

    do iv = 1,veclen
       if(y(iv,1) .gt. rmax .or. y(iv,1) .lt. rmin &
            .or. y(iv,3) .gt. zmax .or. y(iv,3) .lt. zmin) then
          err = 1
          write(*,*) 'particle ',iv,' is out of bounds!'
          return
       end if       
    end do
    
  end function check_bounds
  
end module rk4
