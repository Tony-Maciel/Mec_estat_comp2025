!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! This program verifies that the time complexity of the DFT is O(N^2)
! by returning data for a log-log plot of time taken:
! t ~ N**2 => log(t) ~ 2log(N)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! NOTE: better fit if compile with -O0 and NOT -Ofast, because optimizations
! make the times calculated vary more. (but it doesn't change overall time complexity)

use Fourier_Transform
implicit none

real(dp), parameter :: pi = 4.0_dp*atan(1.0_dp)
real(dp), parameter :: dt = 4e-2_dp   ! time interval 
integer, parameter  :: N_max = 400    ! max number of time intervals
real(dp), parameter :: a1 = 2.0_dp, a2 = 4.0_dp, w1 = 4.0_dp * pi, w2 = 2.5_dp * pi ! w = 2pif
complex(dpc)        :: Y(N_max) , FTY(N_max) 
real(dp)            :: start, fin, average
integer             :: i, io, N, total

! Save data to file
open(newunit=io,file="saida-5-12624850.out",action="write")
total = 50 ! number of averages for times
do N = 50, N_max
    do i = 1, N 
        Y(i) =cmplx(a1*dcos(w1*i*dt) + a2*dsin(w2*i*dt), kind=dp)
    end do
    
    ! time FT
    average = 0.0_dp
    do i = 1, total
        call cpu_time(start)   
        call FT(Y(1:N), FTY(1:N), N)
        call cpu_time(fin)
        average = average + (fin - start)
    end do

    write(io,*) dlog(real(N,dp)), dlog(average/total) 
end do    
close(io)
end
