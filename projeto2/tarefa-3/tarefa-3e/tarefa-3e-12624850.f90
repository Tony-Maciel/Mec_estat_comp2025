!----------------------------------------------------!
!       This program simulates waves discretizing the!
! wave equation. It assumes the wave starts at rest  !
! and with an initial shape. It records the position !
! of the wave after each time step to an output file !
!----------------------------------------------------!

program tarefa3e
!use wave_wfre 
!use wave_wfre_guitar     ! uncomment to use.
!use wave_wfre_guitar_FT
use wave_wfre_FT
use Fourier_Transform
implicit none 

type(wave) :: wave1 

! if using gaussian wave, use wave_wfre and run:
!call wave1%initialize(dx=1e-3_dp, r=1.0_dp, c = 3e2_dp, L=1)
!call wave1%simulate(t=2000, lag=200, filename="saida-3e-1-12624850.dat")  

! if using guitar wave, use wave_wfre_guitar and run:
!call wave1%initialize(dx=1e-3_dp, r=1.0_dp, c = 3e2_dp, L=1)
!call wave1%simulate(t=2000, lag=200, filename="saida-3e-2-12624850.dat")  

! if using guitar wave and FT, use wave_wfre_guitar_FT and Fourier_Transform. Run:
!call wave1%initialize(dx=1e-3_dp, r=1.0_dp, c = 3e2_dp, L=1)
!call wave1%simulate_and_getFT(t=20000, x_0=0.25_dp , filename="saida-3e-3-12624850.dat") 

! if using gaussian wave and FT, use wave_wfre_FT and Fourier_Transform. Run:
call wave1%initialize(dx=1e-3_dp, r=1.0_dp, c = 3e2_dp, L=1, x_0=1.0_dp/20.0_dp)
call wave1%simulate_and_getFT(t=20000, x_0=0.25_dp , filename="saida-3e-6-12624850.dat") 

end program tarefa3e
