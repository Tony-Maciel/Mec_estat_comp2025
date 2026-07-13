!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! This program calculates the Fourier transform of a specified
! signal (given in "data.in") and writes the result to "data.out".
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! NOTE: bounds of j and k were wrong in FT. also make sure that i/dtN 
! can have good precision to distinguish w1 and w2. Else it looks like one peak.
! e.g. make N large, so i can vary a lot and make dt not too small.

! NOTE: don't forget to convert the signal to be complex when writing
! in data.in, else when you read it again, error will be:
! "Fortran runtime error: Bad repeat count in item 1 of list input", 
! because fortran will try to read a real as a complex.

use Fourier_Transform
implicit none

real(dp), parameter :: pi = 4.0_dp*atan(1.0_dp)
real(dp), parameter :: dt = 1e-2_dp   ! time interval 
integer, parameter  :: N = 1000       ! number of time intervals
real(dp), parameter :: a1 = 2.0_dp, a2 = 4.0_dp, w1 = 4.0_dp * pi, w2 = 2.5_dp * pi ! w = 2pif
complex(dpc)        :: Y(N) , FTY(N) 
real(dp)            :: read_dt
integer             :: read_N
integer             :: i, io1, io2, io3

! Initialize data to file
open(newunit=io1,file="data.in",action="write")
    write(io1,*) N, dt
    do i = 1, N 
        !Y(i) = a1*dcos(w1*i*dt) + a2*dsin(w2*i*dt)
        write(io1,*) cmplx(a1*dcos(w1*i*dt) + a2*dsin(w2*i*dt), kind=dp)
    end do
close(io1)

! Read data from file. first line has N and dt, rest is the signal from times dt to N*dt
open(newunit=io2,file="data.in",action="read")
    read(io2,*) read_N , read_dt
    do i = 1, N
        read(io2,*) Y(i)
    end do
close(io2)

! Get Fourier transform
call FT(Y, FTY, read_N)

! Save data to file (magnitude plot)
open(newunit=io3,file="data.out",action="write")
    do i = 1, read_N 
        write(io3,*) (i-1)/(read_dt*real(read_N,dp)) , abs(FTY(i))
    end do
close(io3)

end
