!> This program is a simple implementation of a RK4 particle push routine to test theoretical peak performance
!> @author T. Koskela
!> @date Aug 17 2016
!<

! Fortran bindings for a small subset of the NVIDIA Tools Extensions library
! Thanks to Jeff Larkin: https://gist.github.com/jefflarkin/b64f63d79bdc978a2503
#ifdef NVTX
module nvtx
  use iso_c_binding
  public :: nvtxrangepusha, nvtxrangepop
  ! public :: nvtxrangepushaargb
  interface
     ! Annotate the timeline with a message
     ! Parameters:
     ! * string : the message in a string format
     subroutine nvtxrangepusha(string) bind(C, name="nvtxRangePushA")
       use iso_c_binding , only : c_char
       character(kind=c_char) :: string(*)
     end subroutine nvtxrangepusha

     ! Annotate the timeline with both a message and an ARGB color
     ! Parameters:
     ! * string : the message in a string format
     ! * argb   : the color in argb format (example: Z'FF880000'
     ! subroutine nvtxrangepushaargb(string,argb) bind(C, name="_nvtxRangePushAARGB")
     !   use iso_c_binding , only : c_char, c_int
     !   character(kind=c_char) :: string(*)
     !   integer(kind=c_int), value  :: argb
     ! end subroutine nvtxrangepushaargb

     ! Pop the last range off the stack
     subroutine nvtxrangepop() bind(C, name="nvtxRangePop")
     end subroutine nvtxrangepop

     ! Place a mark on the timeline with a message
     ! Parameters:
     ! * string : the message in a string format
     ! NOT YET EXPOSED
     subroutine nvtxMarkA(string) bind(C, name="nvtxMarkA")
       use iso_c_binding , only : c_char
       character(kind=c_char) :: string(*)
     end subroutine nvtxMarkA

     ! Name an OS thread
     ! NOT YET EXPOSED
     subroutine nvtxNameOsThread(tid, string) bind(C, name="nvtxNameOsThread")
       use iso_c_binding , only : c_int, c_char
       integer(kind=c_int) :: tid
       character(kind=c_char) :: string(*)
     end subroutine nvtxNameOsThread
  end interface
end module nvtx
#endif


program toypush

  use params
  use rk4, only: rk4_push
  use particle, only: particle_data
  use initmodule, only : init

#ifdef OPENMP
  use omp_lib
#endif
#ifdef MPI
  use mpi
#endif
#ifdef NVTX
  use nvtx
#endif
  
  implicit none  
  
  integer :: err
  type(particle_data) :: prt
  integer :: it,iblock
  integer :: nblock
  
  integer :: pid

  integer :: num_procs, my_id, ith

  double precision :: t1,t2

  my_id = 0
  num_procs = 1

#ifdef OPENMP
  !$omp parallel
  !$omp master
  if(my_id .eq. 0) then
     write(*,*) 'number of OpenMP threads = ',omp_get_num_threads()
  end if
  !$omp end master
  !$omp end parallel
#endif

#ifdef MPI
  call mpi_init(err)
  call MPI_COMM_RANK (MPI_COMM_WORLD, my_id, err)
  call MPI_COMM_SIZE (MPI_COMM_WORLD, num_procs, err)
#endif

  params_nprt = params_nprtPerRank * num_procs

  if(my_id .eq. 0) write(*,*) 'program toypush started'
  if(my_id .eq. 0) write(*,*) 'veclen = ',veclen
  if(my_id .eq. 0) write(*,*)
  
  if(my_id .eq. 0) write(*,*) 'initializing particles with ',params_nprtperrank,'particles per rank.'
  if(my_id .eq. 0) write(*,*) 'initializing particles with ',params_nprt,'total particles.'
  if(my_id .eq. 0) write(*,*) 'initializing grid      with ',params_nnode,'nodes.'
  if(my_id .eq. 0) write(*,*) 'initializing grid      with ',params_ntri,'triangles.'
  err = init(prt)
  if(err .ne. 0) stop
  if(my_id .eq. 0) write(*,*) 'done initialising'
  if(my_id .eq. 0) write(*,*)
  
  nblock = params_nprt / veclen  
  if(my_id .eq. 0) write(*,*) 'pushing particles in ',nblock,' blocks'
  
  !pid = 88
  !open(unit=15,file='orbit.dat',action='write')

#ifdef OPENMP
  t1 = omp_get_wtime()
#else
  t1 = mpi_wtime()
#endif

  !$omp parallel do private(iblock, it)
  do iblock = 1,nblock

#ifdef MPI
     if (mod(iblock,num_procs) .eq. my_id) then
#endif
     
        do it = 1,nt
           err = rk4_push(prt, iblock)
           
#ifdef VERBOSE
#ifdef OPENMP
           !$omp critical
           if (mod(it,nt) .eq. 0) then
              ith = omp_get_thread_num()
              if(my_id .eq. 0) write(*,*) 'Thread ',ith,' completed block ',iblock
           end if
           !$omp end critical
#endif
#endif
        end do
#ifdef MPI
     end if
#endif
  end do
  !close(15)

#ifdef OPENMP
  t2 = omp_get_wtime()
#else
  t2 = mpi_wtime()
#endif

  if(my_id .eq. 0) write(*,*) 'done pushing'
  if(my_id .eq. 0) write(*,*) 'spent ',t2-t1,'s'
  if(my_id .eq. 0) write(*,*)

  !call mpi_gatherv
  
  if(my_id .eq. 0)  write(*,*) 'finalizing'
  err = finalize(prt)
#ifdef MPI
  call mpi_finalize(err)
#endif
  if(my_id .eq. 0) write(*,*) 'done finalizing'
  if(my_id .eq. 0) write(*,*)
  
  if(my_id .eq. 0) write(*,*) 'program toypush successfully completed!'
  
contains
  
  function finalize(prt) result(err)

    use rk4, only: rk4_deallocate
    use particle, only : particle_data, particle_deallocate
#ifdef USEIO
    use particleio, only : particleio_write
#endif
    
    implicit none

    type(particle_data) :: prt
    integer :: err

#ifdef USEIO
    err = particleio_write(prt, 'endstate.dat')
#endif
    
    err = particle_deallocate(prt)
    err = rk4_deallocate()
    
  end function finalize
  
end program toypush
