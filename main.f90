PROGRAM main

  ! Test out distance calculations from metric tensor and gdot subroutine

  IMPLICIT NONE

  INTERFACE
    FUNCTION factorial(n)
      INTEGER :: factorial 
      INTEGER, INTENT(IN) :: n
    END FUNCTION factorial

    REAL FUNCTION gdot(a,b,gij)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION a(3), b(3), gij(3,3)
    END FUNCTION gdot
  END INTERFACE

  INTEGER :: i,j,k,l,m,n
  INTEGER :: ntimesteps, nions, totatoms, nCr
  INTEGER, PARAMETER :: r=2
  INTEGER, ALLOCATABLE, DIMENSION(:) :: natom
  REAL*8 :: disp, avg_vol
  REAL*8, DIMENSION(:) :: v1(3), v2(3)
  REAL*8, ALLOCATABLE, DIMENSION(:) :: volume
  REAL*8, DIMENSION(:,:) :: gij(3,3)
  REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: acell, xyz
  CHARACTER(len=3), ALLOCATABLE, DIMENSION(:) :: atype
  LOGICAL :: sim_type, file_exists 

  ! Find parameters: number of timesteps, number of ions in system, and total atoms
  CALL find_params(ntimesteps, nions, totatoms)

  ! Allocate arrays according to parameters
  ALLOCATE(natom(nions))
  ALLOCATE(atype(nions))
  ALLOCATE(volume(ntimesteps))
  ALLOCATE(acell(ntimesteps,3,3))
  ! ALLOCATE(gij(ntimesteps,3,3))
  ALLOCATE(xyz(ntimesteps,totatoms,3))

  ! Read in XDATCAR
  CALL read_xdat(ntimesteps,nions,totatoms,acell,volume,natom,xyz,atype, sim_type)

  ! True == NPT or NPH simulation
  ! False == NVT simulation
  PRINT*, sim_type
  ! IF (sim_type) THEN
  !   ! NPT/NPH simulation
  !   WRITE(*,*) "NPT/NPH SIMULATION. INPUT EQUILIBRIUM VOLUME: "
  !   READ(*,*) avg_vol
  ! ELSE 
  !   ! NVT simulation
  !   avg_vol = volume(1)
  !   WRITE(*,*) "NVT SIMULATION. VOLUME IS: ", avg_vol
  ! END IF

  ! Find total number of atom-pair combinations
  m = nions 
  nCr = (factorial(m) / (factorial(r) * factorial(m-r))) + nions 
  WRITE(*,*) 'Total number of atom type combinations: ', nCr

  ! Calculate radial distribution function
  CALL calc_rdf(ntimesteps, natom, atype, nions, nCr, acell, volume, xyz, sim_type)
  WRITE(*,*)' TEST2'

  ! Deallocate arrays
  DEALLOCATE(natom, atype, volume, acell, xyz)

END PROGRAM main