!----------------------------------------------------!
!  This module contains the essential functions to   !
!simulate a wave with one FIXED END and one free end,!
! starting at rest.                                  !
!----------------------------------------------------!


module wave_wfre_guitar_FT
use, intrinsic :: iso_fortran_env, only: sp=>real32 , dp=>real64
use Fourier_Transform
implicit none

type:: wave
  real(dp)              :: dx, c, r
  integer               :: L, npoints
  real(dp), allocatable :: sys(:,:)

  contains 

  procedure :: initialize, iterate, simulate_and_getFT
end type 

contains 

  subroutine initialize(self, dx, r, c, L)
    class(wave), intent(out) :: self
    real(dp), intent(in)     :: dx, r, c
    integer, intent(in)      :: L
    real(dp)                 :: x_kink
    integer                  :: i
  
    self%dx      = dx 
    self%c       = c 
    self%r       = r
    self%L       = L
    self%npoints = int(L/dx)
    allocate(self%sys(0:self%npoints, 3)) ! same wave at 3 different times. 
    
    ! INITIALIZING. MODIFY THIS...
    x_kink   = real(self%npoints,dp) / 4.0_dp
    do i = 1, int(x_kink)  ! y = x
      self%sys(i,:) = self%dx*i
    end do

    do i = int(x_kink), self%npoints-1  ! y = -x/3 + L/3
      self%sys(i,:) = -1.0_dp/3.0_dp*self%dx*i + 1.0_dp/3.0_dp*self%npoints*self%dx
    end do
    ! BOUNDARY CONDITIONS: FIXED END + FREE END
    self%sys(0,:)            = 0.0_dp 
    self%sys(self%npoints,:) = self%sys(self%npoints-1,:)
  end subroutine initialize

  ! finds next position of wave.
  subroutine iterate(self)
    class(wave), intent(in out) :: self
    integer                     :: i 

    do i = 1, self%npoints-1
      self%sys(i,3) = 2.0_dp*(1.0_dp - self%r**2)*self%sys(i,2) + &
        (self%r**2)*(self%sys(i+1,2) + self%sys(i-1,2)) - self%sys(i,1)
    end do
    ! BOUNDARY CONDITIONS: FIXED END + FREE END
    self%sys(0,3)            = 0.0_dp 
    self%sys(self%npoints,3) = self%sys(self%npoints-1,3)

    ! get ready for next iteration 
    do i = 0, self%npoints 
      self%sys(i,1) = self%sys(i,2)
      self%sys(i,2) = self%sys(i,3)
    end do
  end subroutine iterate 
       
  ! simulates waves until time t, saving Y(t) at x = x_0 (e.g., L/4)
  subroutine simulate_and_getFT(self, t, x_0, filename)
    class(wave), intent(in out) :: self 
    integer, intent(in)         :: t                   ! number of iterations.
    real(dp), intent(in)        :: x_0                 ! position to get signal.
    character(len=*), intent(in):: filename            ! file to overwrite data.
    complex(dpc), allocatable   :: data(:), FTdata(:)  ! signal and its FT.
    integer                     :: i, io, x0_idx
    real(dp)                    :: dt
  
    allocate(data(t+1)) 
    allocate(FTdata(t+1))
    x0_idx = int(self%npoints * x_0) ! index of x_0 to use in self%sys(:,:)
    dt     = self%dx*self%r/self%c

    ! 
    do i = 1, (t + 1)
      data(i) = cmplx(self%sys(x0_idx,2), kind=dp)
      call self%iterate()
    end do

    call FT(data(:), FTdata(:), t+1)

    open(newunit=io, file=filename, action="write")
      do i = 1, t+1 
        write(io,*) (i-1)/(dt*real(t+1,dp)), abs(FTdata(i)) 
      end do
    close(io)
        
  end subroutine simulate_and_getFT
        
end module wave_wfre_guitar_FT
