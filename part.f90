!!****h* BrownianDynamics/ParticleClass
!!
!! NAME
!!  ParticleClass
!!
!! SYNOPSIS
!!  This module is used for storing data about the particles in the simulation.
!!  Encapsulating the data this way makes the program much more clear.
!!
!! NOTES
!!  Located in file part.f90
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
module ParticleClass
use VectorClass
implicit none

    !!****t* ParticleClass/Particle
    !!
    !! NAME
    !!  ParticleClass
    !!
    !! SYNOPSIS
    !!  This derived type contains all the data concerning a particle in the
    !!  simulation cell.  This includes the particle position, the peclet
    !!  number, the radii of the particle and the initial position in the
    !!  simulation cell.
    !!
    !! NOTES
    !!  Located in file part.f90
    !!******
    type :: Particle
        type ( Vector ) pos             ! Particle position. 
        type ( Vector ) pecDTInvRad     ! The degree of external forcing.
        type ( Vector ) init            ! The initial position of the particle.
        Real( 8 ) :: rat                ! The radius of the particle.
        Real( 8 ) :: pec                ! The peclet number forcing the particle.
        Real( 8 ) :: sqrt2DTInvRad      ! Used for Brownian displacement scaling.
        Real( 8 ) :: pecDT              ! Pec * dT    - Luis

        type ( Vector ) abspos          ! Absolute position
        type ( Vector ) poscounter      ! Absolute position counter
    end type Particle

end module ParticleClass
