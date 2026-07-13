!----------------------------------------------------!
!       This program simulates waves discretizing the!
! wave equation. It assumes the wave starts at rest  !
! and with an initial shape. It records the position !
! of the wave after each time step to an output file !
!----------------------------------------------------!

program tarefa3a
use wave_wfe_guitar_FT
use Fourier_Transform
implicit none 

type(wave) :: wave1 

call wave1%initialize(dx=1e-3_dp, r=1.0_dp, c = 3e2_dp, L=1)
call wave1%simulate_and_getFT(t=10000, x_0=0.25_dp , filename="saida-3a-12624850.dat") 
! missing frequencies k = 4, 8, 12, ...

end program tarefa3a
