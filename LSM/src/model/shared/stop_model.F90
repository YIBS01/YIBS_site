#include "rundeck_opts.h"
#if ( defined USE_ESMF )  || ( defined USE_MPP )
#define USE_MPI
#endif

subroutine stop_model( message, retcode )
!@sum Aborts the execution of the program. Passes an error message and
!@+ a return code to the calling script. Should be used instead of STOP
  USE PARAM
  implicit none
!@var message an error message (reason to stop)
  character*(*), intent (in) :: message
!@var retcode return code to be passed to the calling script
  integer, intent(in) :: retcode
  integer, parameter :: iu_err = 9
  integer :: rank
#ifdef USE_MPI
  integer :: mpi_err
#  ifdef MPI_DEFS_HACK
#  include "mpi_defs.h"
#  endif
#include "mpif.h"
#endif

#ifdef USE_MPI
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, mpi_err)
#else
  rank =0
#endif
  ! skip writing status file for retcode<0
  if ( retcode >= 0 ) call write_run_status( message, retcode )
  if (rank == 0) then
    write (6,'(//2(" ",132("*")/))')
    write (6,*) ' Program terminated due to the following reason:'
    write (6,*) ' >>  ', message, '  <<'
    write (6,'(/2(" ",132("*")/))')
  endif

  call sys_flush(6)

  if ( retcode > 13 ) then
    write (0,*) 'Model crashed due to ',message
#ifdef USE_MPI
    !??? bad: the next line will prevent a job from terminating unless
    !???          all processors reach this point
    !??? bug: without it, jobs don't terminate even if
    !???          all processors reach this point
    call mpi_finalize(mpi_err)
    !??? hopefully, we can get rid of the above line soon
    call mpi_abort(MPI_COMM_WORLD, retcode, iu_err)
#else
    call sys_abort
#endif
  else
#ifdef USE_MPI
    call mpi_finalize(mpi_err)
#endif
    call exit_rc (0)
  endif

end subroutine stop_model

subroutine throwException(message, retcode)
!@sum Either invokes pFUnit exception for testing or
!@+ stop_model() for run-time testing.
!@auth T. Clune
#ifdef USE_PFUNIT
  use pFUnit
#endif
  character(len=*), intent(in) :: message
  integer, intent(in) :: retcode

#ifdef USE_PFUNIT
  call throw(Exception(message))
#else
  call stop_model(message, retcode)
#endif
end subroutine throwException
