!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program contains subroutines to calculate
! the Fourier Transform (FT) and the Inverse 
! Fourier Transform (IFT) of a signal with equal 
! time spacing. (NOT OPTIMIZED AT ALL.)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! NOTE: could optimize by calculating exp(twopi/N) once and every iteration just multiplying by itself.

module Fourier_Transform
        ! cexp -> simple precision , zexp -> double precision (not fortran standard)
        use, intrinsic :: iso_fortran_env, only: sp=>real32 , dp=>real64
        implicit none

        integer, parameter      :: dpc = kind((1.0_dp, 1.0_dp)) ! for double precision complex numbers
        real(dp), parameter     :: twopi = 8.0_dp*atan(1.0_dp)
        complex(dpc), parameter :: im  = (0.0_dp , 1.0_dp)

        contains 

        subroutine FT(Y, res, N) 
        ! Calculates FT of Y and puts result in res(:)
        complex(dpc) , intent(in)     :: Y(:)  ! Signal separated by EQUAL time intervals.
        complex(dpc) , intent(in out) :: res(:)! FT of Y. 
        integer , intent(in)          :: N     ! length of Y and res.
        integer                       :: j, k

        do k = 0, N-1 ! for each FT component...
            res(k+1) = (0.0_dp , 0.0_dp) 
            do j = 0,N-1 ! get FT from Y...
                res(k+1) = res(k+1) + Y(j+1)*zexp(twopi * im * j * k /real(N,dp))
            end do
        end do

        end subroutine FT


        subroutine IFT(Y, res, N) 
        ! Calculates inverse FT of Y and puts result in res(:)
        complex(dpc) , intent(in)     :: Y(:)  ! Signal separated by EQUAL time intervals.
        complex(dpc) , intent(in out) :: res(:)! Inverse FT of Y. 
        integer , intent(in)          :: N     ! length of Y and res.
        integer                       :: j, k

        do j = 0, N-1 ! for each FT component...
            res(j+1) = (0.0_dp , 0.0_dp) 
            do k = 0,N-1 ! get IFT from Y...
                res(j+1) = res(j+1) + Y(k+1)*zexp(-twopi * im * j * k / real(N,dp))
            end do
            res(j+1) = res(j+1) / N
        end do

        end subroutine IFT
end module Fourier_Transform
