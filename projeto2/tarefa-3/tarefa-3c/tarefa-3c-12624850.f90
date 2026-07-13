!----------------------------------------------------!
!       This program simulates waves discretizing the!
! wave equation. It assumes the wave starts at rest  !
! and with an initial shape. It records the position !
! of the wave after each time step to an output file !
!----------------------------------------------------!

program tarefa3c
use wave_wfe_FT
use Fourier_Transform
implicit none 

type(wave) :: wave1 

call wave1%initialize(dx=1e-3_dp, r=1.0_dp, c = 3e2_dp, L=1, x_0=1.0_dp/3.0_dp)
call wave1%simulate_and_getFT(t=10000, x_0=0.25_dp , filename="saida-3c-12624850.dat") 

end program tarefa3c
