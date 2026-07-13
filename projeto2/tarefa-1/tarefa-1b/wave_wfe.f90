---------------------------------------------------!
!  This module contains the essential functions to   !
!simulate a wave with FIXED ENDS, starting at rest.  !
! (wave_wfe == "wave with fixed ends.")              !
!----------------------------------------------------!

! subroutines: initialize (have function to specify shape), iterate (impose boundary conditions, but dont have function for it),
! simulate(specify final time), write (specify file. Later will have to save only one part of the wave for spectral analysis).

module wave_wfe  
use, intrinsic :: iso_fortran_env, only: sp=>real32 , dp=>real64
implicit none

type:: wave
  real(dp)              :: dx, c, r
  integer               :: L, npoints
  real(dp), allocatable :: sys(:,:)

  contains 

  procedure :: initialize, iterate, simulate
end type 

contains 

  subroutine initialize(self, dx, r, c, L)
    class(wave), intent(out) :: self
    real(dp), intent(in)     :: dx, r, c
    integer, intent(in)      :: L
    real(dp)                 :: x_0 , sigma
    integer                  :: i
  
    self%dx      = dx 
    self%c       = c 
    self%r       = r
    self%L       = L
    self%npoints = int(L/dx)
    allocate(self%sys(0:self%npoints, 3)) ! same wave at 3 different times. 
    
    ! INITIALIZING. MODIFY THIS...
    x_0   = real(self%npoints,dp) / 3.0_dp
    sigma = real(self%npoints,dp) / 30.0_dp
    do i = 1, self%npoints-1 
      self%sys(i,:) = dexp(-((i - x_0) / sigma)**2)
    end do
    ! BOUNDARY CONDITIONS: FIXED ENDS
    self%sys(0,:)            = 0.0_dp 
    self%sys(self%npoints,:) = 0.0_dp
  end subroutine initialize

  ! finds next position of wave.
  subroutine iterate(self)
    class(wave), intent(in out) :: self
    integer                     :: i 

    do i = 1, self%npoints-1
      self%sys(i,3) = 2.0_dp*(1.0_dp - self%r**2)*self%sys(i,2) + &
        (self%r**2)*(self%sys(i+1,2) + self%sys(i-1,2)) - self%sys(i,1)
    end do
    ! BOUNDARY CONDITIONS: FIXED ENDS
    self%sys(0,3)            = 0.0_dp 
    self%sys(self%npoints,3) = 0.0_dp

    ! get ready for next iteration 
    do i = 0, self%npoints 
      self%sys(i,1) = self%sys(i,2)
      self%sys(i,2) = self%sys(i,3)
    end do
  end subroutine iterate 
       
  ! simulates waves until time t, saving results to given file.
  subroutine simulate(self, t, lag, filename)
    class(wave), intent(in out) :: self 
    integer, intent(in)         :: t, lag   ! number of iterations and time between saves.
    character(len=*), intent(in):: filename ! file to overwrite data.
    real(dp), allocatable       :: data(:,:)! For column first, I have to save before write.
    integer                     :: i, j, io, idx
    
    allocate(data(0:self%npoints,int(t/lag)+2)) 

    ! save initial configuration
    do i = 0, self%npoints
      data(i,1) = i*self%dx
      data(i,2) = self%sys(i,2)
    end do 

    idx = 3
    do i = 1, (t - mod(t,lag)) ! avoid unnecessary iterations
      call self%iterate()
      if (mod(i,lag) .eq. 0) then 
        data(:,idx) = self%sys(:,2)
        idx = idx + 1
      end if
    end do

    open(newunit=io, file=filename, action="write")
      do i = 0, self%npoints 
        write(io,*) (data(i,j), j = 1, (int(t/lag)+2)) ! "dynamic formatting" :o
      end do
    close(io)
        
  end subroutine simulate
        
end module wave_wfe
