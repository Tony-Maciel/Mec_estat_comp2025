!----------------------------------------------------!
!       This program simulates waves discretizing the!
! wave equation. It assumes the wave starts at rest  !
! and with an initial shape. It records the position !
! of the wave after each time step to an output file !
!----------------------------------------------------!

program tarefa2
use wave_wfe_guitar 
implicit none 

type(wave) :: wave1 

call wave1%initialize(dx=1e-3_dp, r=1.0_dp, c = 3e2_dp, L=1)
call wave1%simulate(t=2000, lag=200, filename="saida-2-12624850.dat") 
! t = 500, lag = 100 -> see propagation

end program tarefa2
