PROGRAM countlines
IMPLICIT NONE
  integer :: lines
  lines = 0

  OPEN ( unit = 1, file='Data/periodic_pos.out' , status="unknown")

  DO
     READ (1,*, END=10)
     lines = lines + 1
  END DO

  10 write(*,*) lines

  write(*,*) mod(11.0,4.0)
END PROGRAM countlines
