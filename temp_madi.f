!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE MPI_INST_PRO(OPACK,IPACK)
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      include 'header'
      INCLUDE "mpif.h"
      INCLUDE 'header_mpi'
      integer :: send_rank!,en,st
      REAL(8),DIMENSION(0:NXM,1:NY*NPROCS) :: IPACK
      REAL(8),DIMENSION(0:NXM,1:NY) :: OPACK
!      REAL(8),DIMENSION(0:NXM,1:NY,1:NPROCS) ::temp

      send_rank=0
! READ IN LOCAL DT BASED ON CFL
!      OPACK = DT

! SEND ALL LOCAL DT'S TO EVERY PROCESS
      CALL MPI_GATHER(OPACK,NY*(NXM+1),MPI_DOUBLE_PRECISION,IPACK,NY
     + *(NXM+1), MPI_DOUBLE_PRECISION,send_rank,MPI_COMM_WORLD,IERROR)

!     T = MINVAL(IPACK)
!       IF ( RANK .EQ. 0) THEN
!        temp(:,:) = 0.0d0
!        st=1
!        DO I = 1,NPROCS
!          en=st+SIZE_1-1
!          DO J =1,SIZE_2
!           temp(st:en,J) =  STAT_TMP(:,J,I)
!          ENDDO
!          st=en+1
!        ENDDO
!      ENDIF

      RETURN
      END

