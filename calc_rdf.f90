SUBROUTINE calc_rdf(ntimesteps, natom, atype, nions, nCr, acell, avg_vol, xyz, sim_type)

  IMPLICIT NONE

  INTERFACE
    REAL FUNCTION gdot(a,b,gij)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION a(3), b(3), gij(3,3)
    END FUNCTION gdot
  END INTERFACE

  INTEGER :: i,j,k,l,m,n,ab,cd,ibins,nbins
  INTEGER, INTENT(IN) :: ntimesteps, nions, nCr
  INTEGER, DIMENSION(:), INTENT(IN) :: natom(3)
  INTEGER, DIMENSION(:) :: ntype(nions+1)
  REAL*8 :: disp, delta, density, vshell, rmax, r
  REAL*8, INTENT(IN) :: avg_vol
  REAL*8, DIMENSION(:) :: v1(3)
  REAL*8, DIMENSION(:,:) :: gij(3,3), gr(nCr,20000), hist(nCr,20000)
  REAL*8, DIMENSION(:,:,:), INTENT(IN) :: acell(ntimesteps,3,3), xyz(ntimesteps,SUM(natom),3)
  CHARACTER(len=3), DIMENSION(:) :: aname(nions+1)
  CHARACTER(len=3), DIMENSION(:), INTENT(IN) :: atype(nions)
  LOGICAL, INTENT(IN) :: sim_type
  REAL, PARAMETER :: pi=3.1415926535897932385

  ntype(:) = 0
  rmax = avg_vol**(1./3.)
  WRITE(*,*) 'RMAX', rmax

  ! Assign number of atoms and atom names to arrays: ntype and aname
  DO i = 1,nions
    ntype(i+1) = ntype(i+1) + ntype(i) + natom(i)
    aname(i+1) = atype(i)
    PRINT*, ntype(i+1), aname(i+1)
  END DO

  OPEN(UNIT=5,FILE='atoms',STATUS='UNKNOWN')
    ! print out atoms pairs for python plotting legend
  DO m=2,nions+1
    DO l = m,nions+1
      WRITE(5,*) aname(m), aname(l)
    END DO
  END DO
  CLOSE(5)

  ! Initialize histograms
  hist(:,:) = 0
  ibins = 0

  ! Determine whether XDATCAR is an NPT/NPH or NVT simulation
  IF (sim_type) THEN 
    PRINT*, sim_type
    PRINT*, '************************************** Calculating NPT or NPH displacements **************************************'
    
    ! Initialize RDF parameters
    delta = 0.05
    WRITE(*,*) 'Delta: ', delta
    nbins = rmax/delta
    
    ! Calculate distance vectors between two points. Then take the dot product of vectors with metric tensor Gij
    n = 0
    DO m = 2,nions+1
      DO l = m,nions+1
        n = n + 1
        WRITE(*,*)
        WRITE(*,*) 'Computing atom ', aname(m), aname(l), ' pair'
        WRITE(*,*)
        ab = ntype(m-1)+1
        cd = ntype(l-1)+1
        ! WRITE(*,*) ab, ntype(m)-1, cd, ntype(l)
        DO i = 1,ntimesteps
          gij(:,:) = 0.0
          CALL metric(acell(i,:,:),gij(:,:),1.0)
          DO j = ab,ntype(m)-1
            DO k = cd,ntype(l)
              v1(1) = (xyz(i,j,1) - xyz(i,k,1)) - ANINT((xyz(i,j,1) - xyz(i,k,1)))
              v1(2) = (xyz(i,j,2) - xyz(i,k,2)) - ANINT((xyz(i,j,2) - xyz(i,k,2)))
              v1(3) = (xyz(i,j,3) - xyz(i,k,3)) - ANINT((xyz(i,j,3) - xyz(i,k,3)))
              disp = SQRT(gdot(v1(:), v1(:), gij(:,:)))

              ibins = (disp/delta) + 1
              hist(n,ibins) = hist(n,ibins) + 1
              IF (ibins.eq.1) hist(n,ibins) = 0

            END DO
          END DO
          IF (mod(i,100) .eq. 0) write(6,'(a25,i7,a10)',ADVANCE='NO') 'Computing displacement ',i,'         '//CHAR(13)
        END DO
      END DO
    END DO

  ELSE
    PRINT*, sim_type
    PRINT*, '************************************** Calculating NVT displacements **************************************'

    ! Calculate metric tensor Gij
    CALL metric(acell(1,:,:),gij(:,:),1.0)
    ! print*, (gij(:,:))

    ! Initialize number of bins
    nbins = 250
    delta = rmax/nbins

    ! Calculate distance vectors between two points. Then take the dot product of vectors with metric tensor Gij
    n = 0
    DO m = 2,nions+1
      DO l = m,nions+1
        n = n + 1
        WRITE(*,*)
        WRITE(*,*) 'Computing atom ', aname(m), aname(l), ' pair'
        WRITE(*,*)
        ab = ntype(m-1)+1
        cd = ntype(l-1)+1
        ! WRITE(9,*) ab, cd
        DO i = 1,ntimesteps
          DO j = ab,ntype(m)-1
            DO k = cd,ntype(l)
              v1(1) = (xyz(i,j,1) - xyz(i,k,1)) - ANINT((xyz(i,j,1) - xyz(i,k,1)))
              v1(2) = (xyz(i,j,2) - xyz(i,k,2)) - ANINT((xyz(i,j,2) - xyz(i,k,2)))
              v1(3) = (xyz(i,j,3) - xyz(i,k,3)) - ANINT((xyz(i,j,3) - xyz(i,k,3)))
              disp = SQRT(gdot(v1(:), v1(:), gij(:,:)))

              ibins = (disp/delta) + 1
              hist(n,ibins) = hist(n,ibins) + 1
              IF (ibins.eq.1) hist(n,ibins) = 0

            END DO
          END DO
          IF (mod(i,100) .eq. 0) write(6,'(a25,i7,a10)',ADVANCE='NO') 'Computing displacement ',i,'         '//CHAR(13)
        END DO
      END DO
    END DO

  END IF

  WRITE(*,*) 'CALCULATING RADIAL DISTRIBUTION FUNCTION'

  OPEN(UNIT=5,FILE='atoms',STATUS='UNKNOWN')

  vshell = (4 * pi * delta)
  write(*,*) 'nbins: ', nbins, 'delta: ', delta

  130  FORMAT(f18.10, 8X, f18.10)
  131  FORMAT(A3, 8X, A3)

  k = 21
  n = 0
  DO m = 2,nions+1
    DO l = m,nions+1
      WRITE(k,130)
      k = k + 1
      n = n + 1
      density = natom(m-1)*natom(l-1)
      WRITE(*,*) k, aname(m), aname(l)
      WRITE(*,*) density
      DO i = 1, nbins
        r = (i * delta) - (delta/2.)
        IF (r.ge.rmax/2) EXIT
        gr(n,i) = hist(n,i)*avg_vol/vshell/r**2./density/ntimesteps
        WRITE(k,130) r, gr(n,i)
      END DO
    END DO
  END DO

  CALL EXECUTE_COMMAND_LINE("rm fort.21")
  CALL EXECUTE_COMMAND_LINE("paste fort.2* > rad1")
  CALL EXECUTE_COMMAND_LINE("rm fort.2*")

END SUBROUTINE calc_rdf