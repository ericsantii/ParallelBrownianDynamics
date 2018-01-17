!!****m* BrownianDynamics/DiffCalcu
!!
!! NAME
!! DiffCalcu
!!
!! SYNOPSIS
!! This subroutine calculate the mean square displacement
!! of the positions and rotation of the osmotic motors
!! and calculate the diffusivity coefficient for translation
!! and rotation
!!
!! NOTES
!! Located in file diff.f90
!!
!! CREATION DATE
!!  06/27/2009
!!
!! AUTHOR
!!   Glenn Vidal
!!
!! COPYRIGHT
!!   06-April-2010 - Glenn Vidal
!!*****

program DiffCalcu
implicit none
   Integer :: i

   Real( 8 ), dimension( 3 ) :: diffur1, Sddifr1
   Real( 8 ), dimension( 3 ) :: diffur2, Sddifr2
   Real( 8 ), dimension( 3 ) :: diffur3, Sddifr3
   Real( 8 ), dimension( 3 ) :: diffur4, Sddifr4
   Real( 8 ), dimension( 3 ) :: diffur5, Sddifr5
   Real( 8 ), dimension( 3 ) :: diffur6, Sddifr6
   Real( 8 ), dimension( 3 ) :: diffur7, Sddifr7

   Real( 8 ), dimension( 7 ) :: alpha

!! input data
   open ( unit = 100, file = 'OUT/DIFFU/with_affine/Ps_0.1/diffrotat_0.0.plo', status = 'unknown' )
   open ( unit = 101, file = 'OUT/DIFFU/with_affine/Ps_0.1/diffrotat_0.01.plo', status = 'unknown' )
   open ( unit = 102, file = 'OUT/DIFFU/with_affine/Ps_0.1/diffrotat_0.1.plo', status = 'unknown' )
   open ( unit = 103, file = 'OUT/DIFFU/with_affine/Ps_0.1/diffrotat_0.5.plo', status = 'unknown' )
   open ( unit = 104, file = 'OUT/DIFFU/with_affine/Ps_0.1/diffrotat_1.5.plo', status = 'unknown' )
   open ( unit = 105, file = 'OUT/DIFFU/with_affine/Ps_0.1/diffrotat_10.0.plo', status = 'unknown' )
   open ( unit = 106, file = 'OUT/DIFFU/with_affine/Ps_0.1/diffrotat_50.0.plo', status = 'unknown' ) 

!! output data
   open ( unit = 107, file = 'OUT/DIFFU/with_affine/Ps_0.1/diffrx.plo')
   open ( unit = 108, file = 'OUT/DIFFU/with_affine/Ps_0.1/diffry.plo')
   open ( unit = 109, file = 'OUT/DIFFU/with_affine/Ps_0.1/diffrz.plo')

   do i = 1, 3
         read( 100, * ) diffur1(i), Sddifr1(i)
         read( 101, * ) diffur2(i), Sddifr2(i)
         read( 102, * ) diffur3(i), Sddifr3(i)
         read( 103, * ) diffur4(i), Sddifr4(i)
         read( 104, * ) diffur5(i), Sddifr5(i)
         read( 105, * ) diffur6(i), Sddifr6(i) 
         read( 106, * ) diffur7(i), Sddifr7(i)
   end do

   alpha(1) = 0.0
   alpha(2) = 0.01
   alpha(3) = 0.1
   alpha(4) = 0.5
   alpha(5) = 1.5
   alpha(6) = 10.0
   alpha(7) = 50.0

!! Rotational diffusivity in the x-direction as a function of alpha

   write( 107, * ) alpha(1), diffur1(1), Sddifr1(1)
   write( 107, * ) alpha(2), diffur2(1), Sddifr2(1)
   write( 107, * ) alpha(3), diffur3(1), Sddifr3(1)
   write( 107, * ) alpha(4), diffur4(1), Sddifr4(1)
   write( 107, * ) alpha(5), diffur5(1), Sddifr5(1)
   write( 107, * ) alpha(6), diffur6(1), Sddifr6(1)
   write( 107, * ) alpha(7), diffur7(1), Sddifr7(1)
   
!*****

!! Rotational diffusivity in the y-direction as a function of alpha

   write( 108, * ) alpha(1), diffur1(2), Sddifr1(2)
   write( 108, * ) alpha(2), diffur2(2), Sddifr2(2)
   write( 108, * ) alpha(3), diffur3(2), Sddifr3(2)
   write( 108, * ) alpha(4), diffur4(2), Sddifr4(2)
   write( 108, * ) alpha(5), diffur5(2), Sddifr5(2)
   write( 108, * ) alpha(6), diffur6(2), Sddifr6(2)
   write( 108, * ) alpha(7), diffur7(2), Sddifr7(2)
    
!*****

!! Rotational diffusivity in the z-direction as a function of alpha

   write( 109, * ) alpha(1), diffur1(3), Sddifr1(3)
   write( 109, * ) alpha(2), diffur2(3), Sddifr2(3)
   write( 109, * ) alpha(3), diffur3(3), Sddifr3(3)
   write( 109, * ) alpha(4), diffur4(3), Sddifr4(3)
   write( 109, * ) alpha(5), diffur5(3), Sddifr5(3)
   write( 109, * ) alpha(6), diffur6(3), Sddifr6(3)
   write( 109, * ) alpha(7), diffur7(3), Sddifr7(3)
    
!*****

end program 
