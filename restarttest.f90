PROGRAM countlines
IMPLICIT NONE
  integer :: lines
  integer :: np
  integer :: nlcut ! number of lines to cut
  CHARACTER(LEN=60) :: cmd
  CHARACTER(LEN=10) :: nlcut_str

  np = 100
  lines = 0
  nlcut = 0 
  OPEN ( unit = 1, file='testfile.dat' , status="unknown")

  DO
     READ (1,*, END=10)
     lines = lines + 1
  END DO

  10 write(*,*) lines
  
  nlcut = mod(lines,np)
  write(*,*) "Excess number of lines =", nlcut
  !1000 format('(i3)')

  nlcut = lines - nlcut
  write(nlcut_str,'(i10)') nlcut
  
  cmd = 'head -'//trim(ADJUSTL(nlcut_str))//' testfile.dat > testfile.dat_temp'
  
  CALL SYSTEM(cmd)
  
  cmd = 'rm testfile.dat ; mv testfile.dat_temp testfile.dat'
  
  CALL SYSTEM(cmd)

  !OPEN ( unit = 1, file='testfile.dat' , status="unknown", access='append')
  !write(1,*) '-----------------------------------TEST-----------------------------------'

END PROGRAM countlines
