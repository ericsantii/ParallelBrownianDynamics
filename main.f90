!!****m* Program/BrownianDynamics
!!
!! NAME
!!  Program/BrownianDynamics
!!
!! SYNOPSIS
!! This PROGRAM runs Brownian dynamics simulations and
!! in particular microrheology experiments.  The PROGRAM is
!! highly modular, so almost nothing happens directly within
!! the PROGRAM statement.  The flow control of the PROGRAM is
!! fairly straightforward.  First the user preferences are read
!! from the file conf.in, then memory is allocated for the various
!! arrays used throughout and all the variables and initialized.  
!! A loop is started and the position of the particles is advanced
!! according to whatever dynamics are relevant.  Then the particles
!! are checked for overlap and their positions corrected.  Output
!! is generated periodiCALLy and various statistical measures are
!! taken.  When the loop exits, the final statistical measurements
!! are made and the final data is output.
!!
!! NOTES
!!  Located in file main.f90
!!
!! CREATION DATE
!!    02/28/05
!!
!! AUTHOR 
!!    James Swan
!!
!! COPYRIGHT
!!    28-Feb-2005 - James William Swan
!!******
PROGRAM BrownianDynamics
  USE Globals
  USE VECTORCLASS
  USE omp_lib
  IMPLICIT NONE
    
  ! Loop counters.
  INTEGER :: i, j, id1, id2, dIFfid
  ! Vectors used for measuring the displacement between particles.
  type ( Vector ) :: dIFf, disp1, disp2, motDisp, MinImage, ranVec, comparison
  ! Scale factors used for moving particles the correct amount.
  REAL( 8 ) :: radID1, radID2, radSum2, dIFf2
  ! DIFference between particle 1 and bath particles, reaction probability
  REAL :: dIFfx, ran0
  ! Counter for number of simulation
  INTEGER :: k, sim, l 
  real :: wtime
  real, dimension(2) :: tarray 
  integer count_0, count_1
  integer count_rate, count_max
  double precision time_init, time_final, elapsed_time
  double precision w_time



  gND = 3 !! Do not change! Dimensionality is defined with variable gNDimensionality

  !  ! Starting time
  !  call system_clock(count_0, count_rate, count_max)
  !  time_init = count_0*1.0/count_rate
  w_time = omp_get_wtime()

  CALL Initialize                    ! init.f90 

  
       
  
    
  
  DO k = 1, gNrun
     sim = k
     write( gsim, 100 ) sim
     
100  format( i4.4 )
     
     CALL OpenFiles                    !init.f90
    
     ! read data of Restart_pos_orient.out file --Luis & Ronal
     if(gRestart.eq.1) then                                           !RONAL-LUIS RESTARTING
       do i=1,gNPart                                                  !RONAL-LUIS RESTARTING
         read(2000, 2001) gT, gR(i,:), gN(i,:), gQuat(i,:)            !RONAL-LUIS RESTARTING
         2001 FORMAT(1X, E15.7E3, 3(1X, E24.15E3), 3(1X, E24.15E3), 4(1X, E24.15E3)) !RONAL-LUIS RESTARTING
       enddo                                                          !RONAL-LUIS RESTARTING
       gT=0.d0                                                        !RONAL-LUIS RESTARTING
     endif                                                            !RONAL-LUIS RESTARTING
  
     DO while ( gT .le. gTf )
        
        CALL CalculateFT              ! dyn.f90         
        CALL UpdateParticles          ! dyn.f90
        CALL OutputInter              ! out.f90
       !CALL ClusterStats
        
        gT = gT + gDT 
        
     END DO
     
     CALL FinalEquilValues
     CALL CloseFiles                   ! init.f90
     
  END DO
  ! ! Ending time
  !  call system_clock(count_1, count_rate, count_max)
  !  time_final = count_1*1.0/count_rate
  !  ! Elapsed time
  !  elapsed_time = time_final - time_init
  w_time = omp_get_wtime() - w_time

  !   ! Write elapsed time
  !  write(*,1003) int(elapsed_time),elapsed_time-int(elapsed_time)
  write(*, '(a,g14.6)') ' Elapsed seconds = ' , w_time
   

END PROGRAM BrownianDynamics
