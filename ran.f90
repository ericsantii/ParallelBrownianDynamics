!!****h* Brownian Motion
!!
!! NAME
!!  Knuth_random
!!
!! SYNOPSIS
!!  This module is handles all the random number generation.  It was
!!  originally written by Donald Knuth in Fortran 77 and then ported
!!  to Fortran 90 by Alan Miller.  The statistics were checked and
!!  confirmed statisfactory, but the details of the random number 
!!  generation was not important.  The routine rand takes an array
!!  and the size of the array and fills the array with random numbers
!!  between zero and one.      
!!
!! NOTES
!!  Located in file ran.f90
!!******
MODULE Knuth_random

! Code converted using TO_F90 by Alan Miller
! Date: 2000-09-10  Time: 16:37:48
! Latest revision - 16 January 2003

! FORTRAN 77 version of "ran_array"
! from Seminumerical Algorithms by D E Knuth, 3rd edition (1997)
!       including the MODIFICATIONS made in the 9th printing (2002)
! ********* see the book for explanations and caveats! *********
! Author: Steve Kifowit
! http://ourworld.compuserve.com/homepages/steve_kifowit
! with modifications by Alan Miller to rnarry and rnstrt based upon
! Knuth's code.

   !!*********************
   !date: 04/03/09-7:22pm
   !Hasta el momento no sabemos que es "rnarry" and "rnstrt".
   !claramente sabemos que son subrutines
   !!*********************

! For Donald Knuth's Fortran 77 versions, go to:
! http://www-cs-faculty.stanford.edu/~knuth/programs
! Look for frng.f and frngdb.f

IMPLICIT NONE
Real*8, SAVE       :: seed ! --> "seed" is a value that does not change
                           ! throgouth the program. This is an initial
                           ! value or is calculated in other subroutine ?.

CONTAINS


SUBROUTINE rand(u, n)     ! Which is the theory for this subroutine.
                          

REAL, INTENT(OUT)    :: u(:)   ! Output-this is why size is not specified.  - Luis -
INTEGER, INTENT(IN)  :: n      ! Imput
integer :: i, j
real*8 :: sum, rr, rrsq, term1, term2, term3, term4
real*8 :: d2p31m, d2p31m_inv

    parameter (d2p31m=2147483647.0)
    parameter (d2p31m_inv=1.d0/d2p31m)

    do i = 1, n

         sum = 0.d0

         do j = 1, 12                        ! why the last value of j is 12 and is not other value? 

            seed = dmod( 16807.0 * seed, d2p31m )  ! dmod stand for double-precision module - By Luis

            sum = sum + seed * d2p31m_inv

         end do

        rr = 0.25d0 * ( sum - 6.d0 )
        rrsq = rr * rr

        term1 = 0.029899776d0 * rrsq + 0.008355968d0  
        term2 = term1 * rrsq + 0.076542912d0
        term3 = term2 * rrsq + 0.252408784d0
        term4 = term3 * rrsq + 3.949846138d0 !term4=function(sum)         

        u(i) = term4 * rr                    !1.414213562d0 * term4 * rr, u(i)=function(sum)
                                             
    end do

!!*************************************
!! Write in screem the random numbers
!!************************************
!    do i=1,n
!      write(*,*) i, u(i)
!    end do
!    stop

END SUBROUTINE rand



SUBROUTINE rnstrt(seed0)          ! This subroutine assigned el name
				  ! "seed0" to "seed".
INTEGER, INTENT(IN)  :: seed0     ! yet does not know the value of "seed"
  				  !
seed = seed0			  !
  				  !
END SUBROUTINE rnstrt             !             

END MODULE Knuth_random


!!****m* BrownianDynamics/RandomVectors
!!
!! NAME
!!  RandomVectors
!!
!! SYNOPSIS
!!  This subroutine calculates the random Brownian displacements
!!  that are scaled in the UpdateParticles routine on the square root
!!  of the time step and the size of the particle.  These are uniformly
!!  distributed random numbers between -sqrt(3) and sqrt(3) for each
!!  orthonormal direction.
!!
!!  -sqrt(3)=-1.73
!!   sqrt(3)= 1.73
!!
!! NOTES
!!  Located in file ran.f90
!!
!! CREATION DATE
!!    02/28/05
!!
!! AUTHOR 
!!    James Swan
!!
!! COPYRIGHT
!!    28-Feb-2005 - James William Swan
!!***************************************************************************************
subroutine RandomVectors
use Globals
use Knuth_random
implicit none
 
    character*3  ipfi
    Integer :: i, j, hh, mm, ss, kk                                 
    Integer :: seed0, time( 3 )                             
    Real, allocatable, dimension( : ) :: u
    Real ( 8 ) norm, normo, dty, ang
    
!!***************************************************************************************
    allocate( u( 3 * gNPart ) )                           

    ! Reseed the random number generator.
    if ( mod( int( gT / gDT ), gReseed ) .eq. 0 ) then   
        call iTime( time ) 				 
	! Where is defined this subroutine		 

        seed0 = 10000 * time( 3 ) + 100 * time( 2 ) + time( 1 )
        seed0 = seed0 + time( 1 ) * time( 2 ) * ( time( 3 ) + time( 1 ) ) + time( 3 )

        call rnstrt( seed0 )
    end if

    ! Get the random number array.

!!***************************************************************************************
!!  Brownian displacement
!!***************************************************************************************

    call rand( u, 3 * gNPart )

    ! Rescale the random numbers and store them as Brownian displacements.

    do i = 1, gNMove					  
    	gBVec( i ) % x = u( 3 * ( i - 1 ) + 1 )          
        gBVec( i ) % y = u( 3 * ( i - 1 ) + 2 )           
        
        if ( g2Dproblem .eq. 1 ) then                   
            gBVec( i ) % z = 0.0                          
        else
            gBVec( i ) % z = u( 3 * ( i - 1 ) + 3 )
        end if
        
    end do
    
    ! Brownian displacement of the particle experiencing
    ! brownian forces
    do i = gNMove + 1, gNPart                             
                                                          
        gBVec( i ) % x = u( 3 * ( i - 1 ) + 1 )          
        gBVec( i ) % y = u( 3 * ( i - 1 ) + 2 )         
                                                          
        if ( g2Dproblem .eq. 1 ) then                     
        	gBVec( i ) % z = 0.0                    
        else                                              
        	gBVec( i ) % z = u( 3 * ( i - 1 ) + 3 )
        end if
        

    end do


!!***************************************************************************************
!! Brownian Rotation
!!***************************************************************************************
        if ( gRotateFlag .eq. 1 ) then			  
            call rand( u, 3 * gNPart )                   
            
            if ( u( 1 ) .lt. gRotateFactor ) then       
                do i = 1, gNMove                          

                    dty = ( 2.D0 * u( i + 1 ) - 1.D0 ) * gPi             
                    orient( i ) = Vector( cos( dty ), sin( dty ), 0 ) 

                end do

            end if


        end if
    deallocate( u ) 

end subroutine RandomVectors


      FUNCTION ran0(idum)
      INTEGER idum,IA,IM,IQ,IR,MASK
      REAL ran0,AM
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836, &
     MASK=123459876)
      INTEGER k
      idum=ieor(idum,MASK)
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      ran0=AM*idum
      idum=ieor(idum,MASK)
      return
      END
