C******************************************************************************|
C diablo.f -> DNS In A Box, Laptop Optimized                       VERSION 0.9
C
C This Fortran 77 code computes incompressible flow in a box.
C
C Primative variables (u,v,w,p) are used, and continuity is enforced with a
C fractional step algorithm.
C
C SPATIAL DERIVATIVES:
C   0, 1, 2, or 3 directions are taken to be periodic and handled spectrally
C   (these cases are referred to as the "periodic", "channel", "duct", and
C    "cavity" cases respectively).
C   The remaining directions are taken to be bounded by walls and handled with
C   momentum- and energy-conserving second-order central finite differences.
C
C TIME ADVANCEMENT
C   Two main approaches are implemented:
C     1. RKW3 on nonlinear terms and CN on viscous terms over each RK substep.
C     2. RKW3 on y-derivative terms and CN on other terms over each RK substep.
C
C The emphasis in this introductory code is on code simplicity:
C   -> All variables are in core.
C   -> The code is not explicitly designed for use with either MPI or SMP.
C   -> Overindexing is not used.
C A few simple high-performance programming constructs are used:
C   -> The inner 2 loops are broken out in such a way as to enable out-of-order
C      execution of these loops as much as possible, thereby leveraging
C      vector and superscalar CPU architectures.
C   -> The outer loops are fairly long (including as many operations as
C      possible inside on a single J plane of data) in order to make effective
C      use of cache.
C Multiple time advancement algorithms are implemented for the periodic,
C channel, duct, and cavity cases in order to compare their efficiency for
C various flows on different computational architectures.  In a few of the
C algorithms, the code is broken out fully into a PASS1/PASS2 architecture
C to maximize the efficient use of cache.
C
C This code was developed as a joint project in MAE 223 (CFD), taught by
C Thomas Bewley, at UC San Diego (spring 2001, spring 2005).
C Primary contributions follow:
C Thomas Bewley was the chief software architect
C John R. Taylor wrote the channel flow solvers
C******************************************************************************|
C
C This code is free software; you can redistribute it and/or modify it
C under the terms of the GNU General Public License as published by the
C Free Software Foundation; either version 2 of the License, or (at your
C option) any later version. This code is distributed in the hope that it
C will be useful, but WITHOUT ANY WARRANTY; without even the implied
C warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
C GNU General Public License for more details. You should have received a
C copy of the GNU General Public License along with this code; if not,
C write to the Free Software Foundation, Inc., 59 Temple Place - Suite
C 330, Boston, MA 02111-1307, USA.
C
C******************************************************************************|

C----|--.---------.---------.---------.---------.---------.---------.-|-------|
      PROGRAM DIABLO
      INCLUDE 'header'
      INCLUDE "mpif.h"
      INCLUDE 'header_mpi'
 
      INTEGER N
      REAL*8  int_time 
      

      WRITE(6,*) 
      WRITE(6,*) '             ****** WELCOME TO DIABLO ******'
      WRITE(6,*)
      CALL INITIALIZE
! Initialize START_TIME for run timing
      CALL DATE_AND_TIME (VALUES=TIME_ARRAY)
      START_TIME=TIME_ARRAY(5)*3600+TIME_ARRAY(6)*60
     &   +TIME_ARRAY(7)+TIME_ARRAY(8)*0.001

C A flag to determine if we are considering the first time-step
      FIRST_TIME=.TRUE.
  
c      FIRST_TIME=.FALSE. 
c      CALL SAVE_STATS(.FALSE.)

c      stop

c        int_time = TIME
c       int_time = 2.90270002252950


      DO TIME_STEP = TIME_STEP+1, TIME_STEP+N_TIME_STEPS
        IF (RANK .eq. 0) THEN
         WRITE(6,*) 'Now beginning TIME_STEP = ',TIME_STEP
        ENDIF 
       DO N=1,N_TH 
c        If ( (TIME-int_time) .LT. 2.d0*PI) then
c           RI_TAU(N) = RI_FINAL*(TIME-int_time)/(2.d0*PI)
c        ELSE
           RI_TAU(N) = RI_FINAL
c        ENDIF    
       ENDDO  

        DO RK_STEP=1,3
          IF (NUM_PER_DIR.EQ.3) THEN
            IF (TIME_AD_METH.EQ.1) CALL RK_PER_1
            IF (TIME_AD_METH.EQ.2) CALL RK_PER_2            
          ELSEIF (NUM_PER_DIR.EQ.2) THEN
            IF (TIME_AD_METH.EQ.1) CALL RK_CHAN_1
            IF (TIME_AD_METH.EQ.2) CALL RK_CHAN_2            
            IF (TIME_AD_METH.EQ.3) CALL RK_CHAN_3            
          ELSEIF (NUM_PER_DIR.EQ.1) THEN
            IF (TIME_AD_METH.EQ.1) CALL RK_DUCT_1
            IF (TIME_AD_METH.EQ.2) CALL RK_DUCT_2            
          ELSEIF (NUM_PER_DIR.EQ.0) THEN
            IF (TIME_AD_METH.EQ.1) CALL RK_CAV_1
            IF (TIME_AD_METH.EQ.2) CALL RK_CAV_2            
          END IF
        END DO
        TIME=TIME+DELTA_T
        FIRST_TIME=.FALSE.
! Save statistics to an output file
        IF (MOD(TIME_STEP,SAVE_STATS_INT).EQ.0) THEN
            CALL SAVE_STATS(.FALSE.)
        END IF
! Save the flow to a restart file
        IF (MOD(TIME_STEP,SAVE_FLOW_INT).EQ.0) THEN
          CALL SAVE_FLOW(.FALSE.)
        END IF
!adding the last_saved module to save for backup

!        IF (MOD(TIME_STEP,1000).EQ.) THEN
!            CALL SAVE_BU(TIME_STEP)
!      save flow after long time interval

        IF (MOD(TIME_STEP,5000).EQ.0) THEN
          call SAVE_FLOW_INTER
        ENDIF
       
!        write(6,*) 'I am here 4'
! Filter the scalar field
!        DO N=1,N_TH
!          IF (FILTER_TH(N)
!     &       .AND.(MOD(TIME_STEP,FILTER_INT(N)).EQ.0)) THEN
!          write(*,*) 'Filtering...'
!          CALL FILTER(N)
!          END IF 
!        END DO

        
! If we are considering a Near Wall Model, it may be necessary to
! filter the velocity field.  If so, use the following:
!          IF (FILTER_VEL
!     &       .AND.(MOD(TIME_STEP,FILTER_VEL_INT).EQ.0)) THEN
!             write(*,*) 'Filtering Velocity'
!           CALL APPLY_FILTER_VEL(CU1,1,NY)
!           CALL APPLY_FILTER_VEL(CU2,1,NY)
!           CALL APPLY_FILTER_VEL(CU3,1,NY)
!        END IF
 
      END DO

! Calculate and display the runtime for the simulation
      CALL DATE_AND_TIME (VALUES=TIME_ARRAY)
      END_TIME=TIME_ARRAY(5)*3600+TIME_ARRAY(6)*60
     &   +TIME_ARRAY(7)+TIME_ARRAY(8)*0.001
      WRITE(*,*) 'Elapsed Time (sec): ',end_time-start_time
      WRITE(*,*) 'Seconds per Iteration: '
     &     ,(end_time-start_time)/N_TIME_STEPS

      TIME_STEP=TIME_STEP-1
      CALL SAVE_FLOW(.TRUE.)
      CALL SAVE_STATS(.TRUE.)
      WRITE(6,*)
      WRITE(6,*) '        ****** Hello world!  Have a nice day! ******'
      WRITE(6,*)
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE INITIALIZE
      INCLUDE 'header'
      INCLUDE "mpif.h"
      INCLUDE 'header_mpi'


      REAL    VERSION, CURRENT_VERSION, DAMP_FACT
      logical RESET_TIME, int_kick
      INTEGER I, J, K, N
      CHARACTER*35 FGRID
      REAL*8 RNUM1,RNUM2,RNUM3
      INTEGER, DIMENSION(:), ALLOCATABLE :: seed
 
      OPEN (11,file='input.dat',form='formatted',status='old')      
      
      WRITE(6,*) 'Note that this code is distributed under the ',
     *           'GNU General Public License.'
      WRITE(6,*) 'No warranty is expressed or implied.'
      WRITE(6,*)

C Read input file.
C   (Note - if you change the following section of code, update the
C    CURRENT_VERSION number to make obsolete previous input files!)

      CURRENT_VERSION=0.9
      READ(11,*)
      READ(11,*)
      READ(11,*)
      READ(11,*)
      READ(11,*) VERSION
      IF (VERSION .NE. CURRENT_VERSION) STOP 'Wrong input data format. '
      READ(11,*)
      READ(11,*) USE_MPI
      READ(11,*)
      READ(11,*) NU, LX, LY, LZ
      READ(11,*)
      READ(11,*) N_TIME_STEPS, DELTA_T, RESET_TIME, VARIABLE_DT, CFL
     &            , UPDATE_DT
      READ(11,*)
      READ(11,*) NUM_PER_DIR, TIME_AD_METH, LES, LES_MODEL_TYPE
      READ(11,*)
      READ(11,*) VERBOSITY, SAVE_FLOW_INT, SAVE_STATS_INT, MOVIE
      READ(11,*)
      READ(11,*) CREATE_NEW_FLOW, IC_TYPE, KICK
      READ(11,*)
      READ(11,*) F_TYPE, UBULK0, PX0, OMEGA0, AMP_OMEGA0, ANG_BETA
      IF (NUM_PER_DIR.eq.3) THEN
      READ(11,*)
      READ(11,*) U_BC_XMIN, U_BC_XMIN_C1, U_BC_XMIN_C2, U_BC_XMIN_C3
      READ(11,*) 
      READ(11,*) V_BC_XMIN, V_BC_XMIN_C1, V_BC_XMIN_C2, V_BC_XMIN_C3
      READ(11,*)
      READ(11,*) W_BC_XMIN, W_BC_XMIN_C1, W_BC_XMIN_C2, W_BC_XMIN_C3
      READ(11,*)
      READ(11,*) U_BC_XMAX, U_BC_XMAX_C1, U_BC_XMAX_C2, U_BC_XMAX_C3
      READ(11,*)
      READ(11,*) V_BC_XMAX, V_BC_XMAX_C1, V_BC_XMAX_C2, V_BC_XMAX_C3
      READ(11,*)
      READ(11,*) W_BC_XMAX, W_BC_XMAX_C1, W_BC_XMAX_C2, W_BC_XMAX_C3
      END IF
      IF (NUM_PER_DIR.gt.0) THEN
      READ(11,*)
      READ(11,*) U_BC_YMIN, U_BC_YMIN_C1, U_BC_YMIN_C2, U_BC_YMIN_C3
      READ(11,*) 
      READ(11,*) V_BC_YMIN, V_BC_YMIN_C1, V_BC_YMIN_C2, V_BC_YMIN_C3
      READ(11,*)
      READ(11,*) W_BC_YMIN, W_BC_YMIN_C1, W_BC_YMIN_C2, W_BC_YMIN_C3
      READ(11,*)
      READ(11,*) U_BC_YMAX, U_BC_YMAX_C1, U_BC_YMAX_C2, U_BC_YMAX_C3
      READ(11,*)
      READ(11,*) V_BC_YMAX, V_BC_YMAX_C1, V_BC_YMAX_C2, V_BC_YMAX_C3
      READ(11,*)
      READ(11,*) W_BC_YMAX, W_BC_YMAX_C1, W_BC_YMAX_C2, W_BC_YMAX_C3
      READ(11,*)
      READ(11,*) STOCHASTIC_FORCING
! Filter for the velocity field
      READ(11,*)
      READ(11,*) FILTER_VEL, FILTER_VEL_INT
      END IF
      IF (NUM_PER_DIR.lt.2) THEN
      READ(11,*)
      READ(11,*) U_BC_ZMIN, U_BC_ZMIN_C1, U_BC_ZMIN_C2, U_BC_ZMIN_C3
      READ(11,*) 
      READ(11,*) V_BC_ZMIN, V_BC_ZMIN_C1, V_BC_ZMIN_C2, V_BC_ZMIN_C3
      READ(11,*)
      READ(11,*) W_BC_ZMIN, W_BC_ZMIN_C1, W_BC_ZMIN_C2, W_BC_ZMIN_C3
      READ(11,*)
      READ(11,*) U_BC_ZMAX, U_BC_ZMAX_C1, U_BC_ZMAX_C2, U_BC_ZMAX_C3
      READ(11,*)
      READ(11,*) V_BC_ZMAX, V_BC_ZMAX_C1, V_BC_ZMAX_C2, V_BC_ZMAX_C3
      READ(11,*)
      READ(11,*) W_BC_ZMAX, W_BC_ZMAX_C1, W_BC_ZMAX_C2, W_BC_ZMAX_C3
! Read Stochastic forcing parameters
      READ(11,*)
      READ(11,*) STOCHASTIC_FORCING
! Filter for the velocity field
      READ(5,*)
      READ(5,*) FILTER_VEL, FILTER_VEL_INT
      END IF
      READ(11,*)
      write(*,*) FILTER_VEL_INT
! Read in the parameters for the N_TH scalars
      DO N=1,N_TH
        READ(11,*)
        READ(11,*) CREATE_NEW_TH(N)
        READ(11,*)
        READ(11,*) FILTER_TH(N), FILTER_INT(N)
        READ(11,*)
        READ(11,*)RI_FINAL,PR(N),BACKGROUND_GRAD(N),DEV_BACK_TH,
     &            SPONGE_FLAG
        READ(11,*)
        READ(11,*) TH_BC_YMIN(N),TH_BC_YMIN_C1(N),TH_BC_YMIN_C2(N)
     &             ,TH_BC_YMIN_C3(N)
        READ(11,*)
        READ(11,*) TH_BC_YMAX(N),TH_BC_YMAX_C1(N),TH_BC_YMAX_C2(N)
     &             ,TH_BC_YMAX_C3(N)
      END DO

C If we are using MPI, then Initialize the MPI Variables
      IF (USE_MPI) THEN
        CALL INIT_MPI
      END IF

      delta_T_lt = delta_T
      b_fact = 60.d0 ;
      ANG_BETA = ANG_BETA*3.14159265/180.0
      Gravity_g= 10.d0 ;
      rho_0    = 1000.0d0 ;
      
      Temp_0     = 1.0d0 ;
      Sal_0    = 35.0d0 ;
      alpha_T  = 8.0d0*10.0d0**(-5.0)
      gamma_S  = 8.0d0*10.0d0**(-4.0)

      m_FP     = -0.060d0
      L_heat   = 3.35*10.0d0**5.0
      C_sp_heat = 4.184*10.0d0**3.0
      MELTING_MODEL = .true.

      Ratio_gr = Gravity_g*alpha_T ;

      IF(N_TH ==2)THEN
       Ratio_gr_Sc(1) =  Gravity_g*alpha_T ;
       Ratio_gr_Sc(2) = -Gravity_g*gamma_S ;
      ENDIF

      
      Grad_rho(1) = 0.0d0 ; !RI_FINAL*rho_0/Gravity_g
      IF(N_TH ==2)THEN
       Grad_rho(2) = 0.0d0 !(10.0)/LX ;
      ENDIF  

!      DO N=1,N_TH
!       if(N=1) then
       RI_TAU(1) = RI_TAU(1)*Gravity_g*alpha_T ;
!       else
       RI_TAU(2) = RI_TAU(2)*Gravity_g*gamma_S ;
!       endif
!      ENDDO
C Initialize grid

      IF (RANK .eq. 0) THEN
      WRITE(6,*) '#############################################'
      WRITE(6,*) 'Box size: LX =',LX,', LY =',LY,', LZ =',LZ,'.'
      WRITE(6,*) '#############################################'
      WRITE(6,*) 'Grid size: NX =',NX,', NY =',NY,', NZ =',NZ,'.'
      WRITE(6,*) 'Inclined angle beta', ANG_BETA*180.0/3.14159265
      WRITE(6,*) '#############################################'
      WRITE(*,*) 'Ceartion of a new flowfield', CREATE_NEW_FLOW
     

      DO N=1,N_TH
        WRITE(6,*) '  Scalar number: ',N
        WRITE(6,*) '  Final Richardson number: ',RI_FINAL
        WRITE(6,*) '  Prandlt number: ',PR(N)
        WRITE(6,*) '  Gravity: ', Gravity_g
        WRITE(6,*) '  Density Gradient: ', Grad_rho(N)
        WRITE(6,*) '  Ratio gravity/rho : ', Ratio_gr_Sc(N)
        WRITE(6,*) '#############################################'
      END DO
        WRITE(6,*) '  Background Temp: ', Temp_0
        WRITE(6,*) '  Background Salinity: ', Sal_0
      ENDIF

      NXM=NX-1
      NYM=NY-1
      NZM=NZ-1

      I = RANK+1

           FGRID ='grids/ygrid_'
     &        //CHAR(MOD(I,100)/10+48)
     &        //CHAR(MOD(I,10)+48) // '.txt'


      IF (NUM_PER_DIR .GT. 0) THEN
         if (rank .eq.0) then
         WRITE (6,*) 'Fourier in X'
         endif
         DO I=0,NX
           GX(I)=(I*LX)/NX
           DX(I)=LX/NX
           IF (VERBOSITY .GT. 3) WRITE(6,*) 'GX(',I,') = ',GX(I)
         END DO
      ELSE
         WRITE (6,*) 'Finite-difference in X'
         OPEN (30,file='xgrid.txt',form='formatted',status='old')
         READ (30,*) NX_T
C Check to make sure that grid file is the correct dimensions
         IF (NX_T.ne.NX) THEN
           WRITE(6,*) 'NX, NX_T',NX,NX_T
           STOP 'Error: xgrid.txt wrong dimensions'
         END IF
         DO I=1,NX+1
           READ(30,*) GX(I)
           IF (VERBOSITY .GT. 3) WRITE(6,*) 'GX(',I,') = ',GX(I)
         END DO 
         DO I=1,NX
           READ(30,*) GXF(I)
           IF (VERBOSITY .GT. 3) WRITE(6,*) 'GXF(',I,') = ',GXF(I)
         END DO
C Define ghost cells, if needed for this grid...
         GXF(0)=2.d0*GXF(1)-GXF(2)
         GXF(NX+1)=2.d0*GXF(NX)-GXF(NXM)
         GX(0)=2.d0*GX(1)-GX(2)
C Define the grid spacings 
         DO I=1,NX+1
           DX(I)=(GXF(I)-GXF(I-1))
         END DO
         DO I=1,NX
           DXF(I)=(GX(I+1)-GX(I))
         END DO
         CLOSE(30)
      END IF

      IF (NUM_PER_DIR .GT. 1) THEN
         if (rank .eq. 0) then
          WRITE (6,*) 'Fourier in Z'
         endif
         DO K=0,NZ
           GZ(K)=(K*LZ)/NZ
           DZ(K)=LZ/NZ
           IF (VERBOSITY .GT. 3) WRITE(6,*) 'GZ(',K,') = ',GZ(K)
         END DO
      ELSE
        
         WRITE (6,*) 'Finite-difference in Z'
         OPEN (30,file='zgrid.txt',form='formatted',status='old')
         READ (30,*) NZ_T
C Check to make sure that grid file is the correct dimensions
         IF (NZ_T.ne.NZ) THEN
           WRITE(6,*) 'NZ, NZ_T',NZ,NZ_T
           STOP 'Error: zgrid.txt wrong dimensions'
         END IF
         DO K=1,NZ+1
           READ(30,*) GZ(k)
           IF (VERBOSITY .GT. 3) WRITE(6,*) 'GZ(',K,') = ',GZ(K)
         END DO 
         DO K=1,NZ
           READ(30,*) GZF(k)
           IF (VERBOSITY .GT. 3) WRITE(6,*) 'GZF(',K,') = ',GZF(K)
         END DO
C Define ghost cells, if needed for this grid...
         GZF(0)=2.d0*GZF(1)-GZF(2)
         GZF(NZ+1)=2.d0*GZF(NZ)-GZF(NZM)
         GZ(0)=2.d0*GZ(1)-GZ(2)
C Define grid spacing 
         DO K=1,NZ+1
           DZ(K)=(GZF(K)-GZF(K-1))
         END DO
         DO K=1,NZ
           DZF(K)=(GZ(K+1)-GZ(K))
         END DO
         CLOSE(30)
      END IF
      
      IF (NUM_PER_DIR .GT. 2) THEN
         if (rank .eq. 0) then
         WRITE (6,*) 'Fourier in Y'
         endif 
         DO J=0,NY
           GY(J)=(J*LY)/NY
           DY(J)=LY/NY
           IF (VERBOSITY .GT. 3) WRITE(6,*) 'GY(',J,') = ',GY(J)
         END DO
      ELSE
         if (rank .eq.0)then
         WRITE (6,*) 'Finite-difference in Y'
         endif
         OPEN (30,file=FGRID,form='formatted',status='old')
         READ (30,*) NY_T
C Check to make sure that grid file is the correct dimensions
         IF (NY_T.ne.NY) THEN
           WRITE(6,*) 'NY, NY_T',NY,NY_T
           STOP 'Error: ygrid.txt wrong dimensions'
         END IF
         DO J=1,NY+1
           READ(30,*) GY(j)
           IF (VERBOSITY .GT. 3) WRITE(6,*) 'GY(',J,') = ',GY(J)
         END DO
         DO J=1,NY
           READ(30,*) GYF(j)
           IF (VERBOSITY .GT. 3) WRITE(6,*) 'GYF(',J,') = ',GYF(J)
         END DO
         CLOSE(30)
 
          
         IF (USE_MPI) THEN
           CALL GHOST_GRID_MPI
         ELSE
C Define ghost cells
           GYF(0)=2.d0*GYF(1)-GYF(2)
           GYF(NY+1)=2.d0*GYF(NY)-GYF(NYM)
           GY(0)=2.d0*GY(1)-GY(2)
         END IF

C Define grid spacing
         DO J=0,NY+2
           DY(J)=(GYF(J)-GYF(J-1))
         END DO
         DO J=0,NY+1
           DYF(J)=(GY(J+1)-GY(J))
         END DO
      END IF


       call MPI_INST_PROF(GYF(1:NY),GYF_TOT)


       IF (RANK .eq. 0)THEN
        call NETCDF_OPEN_STATS_CHAN
       ENDIF


       IF(SPONGE_FLAG)THEN
C Seting up sponge in the domain
c      y_sp_max =  57437.d0/1000000.d0
c      y_sp_min = 0.039915.d0
c      y_sp_min =  39548.d0/1000000.d0

!      IF (rank .eq.NPROCS-1) then
       y_sp_max =  GYF(NY+1)
!      ENDIF

!      write(6,*) 'Sponge before', y_sp_max,rank
      
      CALL  MPI_BCAST_REAL_L(y_sp_max,1,1)  
      if (rank .eq.0)then
      write(6,*) 'Sponge', y_sp_max
      endif

!      y_sp_min = 1.9907902464      
       y_sp_min = 4.8132616724
      

      DO j=0,NY+1
       SPONGE_SIGMA(j) = 0.d0
      ENDDO 

      I=0
c      J=0 
!      WRITE(199,*) rank,100000
    
      DO J=0,NY+1
       IF (I .eq. 0) THEN  
        IF (GYF(j).GT.y_sp_min-0.0000001d0) THEN
          NYC = j
          WRITE(410,*)NYC,GYF(NYC),GYF(NY+1),y_sp_min
          call sponge_modi
           I=1
         ENDIF
        ENDIF
      ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    CALL NEW SPONGE FOR FORCING  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
!      CALL sponge_forcing  
      write(6,*) 'I am here ??????'    
      ENDIF
c       DO J=0,NY+1
c        IF (I.EQ.0)THEN
c         IF (GYF(j).GT.y_sp_min-0.0000001d0) THEN
c          NYC = j 
c          WRITE(410,*)NYC,GYF(NYC),GYF(NY+1),y_sp_min
c          call sponge_modi
c          I = 1
c         ELSE
c          I = 0
c         ENDIF  
c        ENDIF
c      ENDDO
c      DO j =0,NY+1
c       WRITE(199,*)j,SPONGE_SIGMA(j),gyf(j)
c      ENDDO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C Initialize storage arrays.
      DO K=0,NZ+1
        DO I=0,NX+1 
          DO J=0,NY+1
            U1(I,K,J)=0.
            U3(I,K,J)=0.
            U2(I,K,J)=0.
            P (I,K,J)=0.
            R1(I,K,J)=0.
            R2(I,K,J)=0.
            R3(I,K,J)=0.
            F1(I,K,J)=0.
            F2(I,K,J)=0.
            F3(I,K,J)=0.
! Array for LES subgrid model
            NU_T(I,K,J)=0.
            DO N=1,N_TH
              KAPPA_T(I,K,J,N)=0.
              TH(I,K,J,N)     =0.
            END DO
          END DO
            U1(I,K,-1)=0.
            U3(I,K,-1)=0.
            U2(I,K,-1)=0.
            U1(I,K,NY+2)=0.
            U3(I,K,NY+2)=0.
            U2(I,K,NY+2)=0.
            DO N=1,N_TH
              TH(I,K,-1,N)  =0.
              TH(I,K,NY+2,N)=0.
            END DO
        END DO
      END DO

C Initialize FFT package (includes defining the wavenumber vectors).
      CALL INIT_FFT

C Initialize RKW3 parameters.
      H_BAR(1)=DELTA_T*(8.0/15.0)
      H_BAR(2)=DELTA_T*(2.0/15.0)
      H_BAR(3)=DELTA_T*(5.0/15.0)
      BETA_BAR(1)=1.0
      BETA_BAR(2)=25.0/8.0
      BETA_BAR(3)=9.0/4.0
      ZETA_BAR(1)=0.0
      ZETA_BAR(2)=-17.0/8.0
      ZETA_BAR(3)=-5.0/4.0
      
C Initialize case-specific packages.
      IF (NUM_PER_DIR.EQ.3) THEN
        CALL INIT_PER
      ELSEIF (NUM_PER_DIR.EQ.2) THEN
        IF (USE_MPI) THEN
          CALL INIT_CHAN_MPI
        ELSE 
          CALL INIT_CHAN
        END IF
      ELSEIF (NUM_PER_DIR.EQ.1) THEN
        CALL INIT_DUCT
      ELSEIF (NUM_PER_DIR.EQ.0) THEN
        CALL INIT_CAV
      END IF

C Initialize values for reading of scalars
      NUM_READ_TH=0
      DO N=1,N_TH
        IF (CREATE_NEW_TH(N)) THEN
          NUM_READ_TH=NUM_READ_TH 
        ELSE
          NUM_READ_TH=NUM_READ_TH + 1
          READ_TH_INDEX(NUM_READ_TH)=N
        END IF
      END DO
      IF (NUM_PER_DIR.EQ.2) THEN
        CALL CREATE_TH_CHAN 
      ELSE IF (NUM_PER_DIR.EQ.3) THEN
        CALL CREATE_TH_PER
      END IF 
     
      if (rank .eq. 0)then
        write(*,*) 'A new flowfield', CREATE_NEW_FLOW
      endif  

      IF (CREATE_NEW_FLOW) THEN
        IF (NUM_PER_DIR.EQ.3) THEN
          CALL CREATE_FLOW_PER
        ELSEIF (NUM_PER_DIR.EQ.2) THEN
          CALL CREATE_FLOW_CHAN
        ELSEIF (NUM_PER_DIR.EQ.1) THEN
          CALL CREATE_FLOW_DUCT
        ELSEIF (NUM_PER_DIR.EQ.0) THEN
          CALL CREATE_FLOW_CAV
        END IF
        if (rank .eq. 0)then
        write(*,*) 'A new flowfield has been created'
        endif

!        CALL SAVE_STATS(.FALSE.)
        !write(*,*) 5000
      ELSE
        write(*,*) 'Reading flow...'
        CALL READ_FLOW
        write(*,*) 'Done reading flow'

      int_kick = .FALSE.
      if (int_kick) then
C   case of strong stratification for low Reynolds no flow we need to perturb 
C   System 

C Initialize the random number generator
      CALL RANDOM_SEED(SIZE = K)
      Allocate (seed(1:K))
      seed(1:K)=10
      CALL RANDOM_SEED(PUT = seed)


      DO J=1,NY+1
        DO I=0,NKX
          DO K=0,TNKZ
           IF ((I .EQ.0).AND.(K .EQ.0)) THEN
           ELSE
            CU1(I,K,J)=0.d0
            CU2(I,K,J)=0.d0
            CU3(I,K,J)=0.d0
           ENDIF
          ENDDO
         ENDDO
       ENDDO

       DO I=1,NKX
        DO J=1,NY
          DO K=1,TNKZ
C Now, give the velocity field a random perturbation
            CALL RANDOM_NUMBER(RNUM1)
            CALL RANDOM_NUMBER(RNUM2)
            CALL RANDOM_NUMBER(RNUM3)
           CU1(I,K,J)=CU1(I,K,J)+(RNUM1-0.5)*KICK
           CU2(I,K,J)=CU2(I,K,J)+(RNUM2-0.5)*KICK
           CU3(I,K,J)=CU3(I,K,J)+(RNUM3-0.5)*KICK
          ENDDO
        ENDDO
       ENDDO

       CALL REM_DIV_CHAN

!       CALL POISSON_P_CHAN

      endif 

! TEMPORARY!! SCALE THE VELOCITY FLUCTIONS
c      DO J=0,NY+1
c        DO K=0,TNKZ
c          DO I=0,NKX
c            IF ((K .EQ. 0).AND.(I.EQ.0)) THEN
c                DAMP_FACT = 1.0
c            ELSE
c              IF (J .EQ. 0) THEN
c                DAMP_FACT = 1.0
c              ELSE
c                DAMP_FACT = exp(-10.d0*dble((J-1))/dble((NY-1)))
c              ENDIF   
c               CU1(I,K,J)=CU1(I,K,J)*DAMP_FACT
c               CU3(I,K,J)=CU3(I,K,J)*DAMP_FACT             
c            ENDIF
c          END DO
c        END DO
c      END DO
c     
c      DO K =0,TNKZ
c        DO I=0,NKX
c          CU2(I,K,NY+1) = 0.d0
c        ENDDO
c      ENDDO    
c
c      DO J=NY,1,-1
c        DO K=0,TNKZ
c          DO I=0,NKX
c               CU2(I,K,J)=CU2(I,K,J+1)+DYF(J)*(CIKX(I)*CU1(I,K,J)+ 
c     :                     CIKZ(K)*CU3(I,K,J))
c          END DO
c        END DO
c       END DO  

! TEMPORARY!!! REMOVE MEAN VELOCITY
c      DO J=0,NY+1
c        CU1(0,0,J)=0.d0 
c        CU2(0,0,J)=0.d0 
c        CU3(0,0,J)=0.d0 
c      END DO
c        DO J=0,NY+1
c         DO k=0,TNKZ
c          DO i=0,NKX 
c           CU1(i,k,J)=0.d0 
c           CU2(i,k,J)=0.d0 
c           CU3(i,k,J)=0.d0
c          ENDDO
c         ENDDO  
c       END DO 

C Initialize flow.
      IF (RESET_TIME .OR. CREATE_NEW_FLOW) THEN
        PREVIOUS_TIME_STEP=0
        TIME_STEP=0
        TIME=0
      END IF

!        CALL SAVE_STATS(.FALSE.)
        IF (NUM_PER_DIR.EQ.3) THEN
          CALL POISSON_P_PER
        ELSEIF (NUM_PER_DIR.EQ.2) THEN
!          CALL POISSON_P_CHAN
        ELSEIF (NUM_PER_DIR.EQ.1) THEN
          CALL POISSON_P_DUCT
        ELSEIF (NUM_PER_DIR.EQ.0) THEN
          CALL POISSON_P_CAV
        END IF

      END IF

      RETURN
      END





C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE SAVE_STATS(FINAL)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      INCLUDE "mpif.h"
      include 'header_mpi'
 
      LOGICAL FINAL

      IF (NUM_PER_DIR.EQ.3) THEN
        CALL SAVE_STATS_PER(FINAL)          
      ELSEIF (NUM_PER_DIR.EQ.2) THEN
        CALL SAVE_STATS_CHAN(FINAL)          
      ELSEIF (NUM_PER_DIR.EQ.1) THEN
        CALL SAVE_STATS_DUCT(FINAL)          
      ELSEIF (NUM_PER_DIR.EQ.0) THEN
        CALL SAVE_STATS_CAV(FINAL)          
      END IF
      if (rank .eq. 0 ) then
      write(*,*) 'done save_stats diablo'
      endif

      IF (FINAL) THEN
        IF (NUM_PER_DIR.EQ.3) THEN
          CALL VIS_FLOW_PER         
        ELSEIF (NUM_PER_DIR.EQ.2) THEN
!          CALL VIS_FLOW_CHAN         
        ELSEIF (NUM_PER_DIR.EQ.1) THEN
          CALL VIS_FLOW_DUCT          
        ELSEIF (NUM_PER_DIR.EQ.0) THEN
          CALL VIS_FLOW_CAV         
        END IF
      END IF

      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE READ_FLOW
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      INCLUDE "mpif.h"
      INCLUDE 'header_mpi' 
      CHARACTER*35 FNAME
      CHARACTER*35 FNAME_TH
      INTEGER I, J, K, N, NUM_PER_DIR_T, TIME_STEP_T
      REAL*8   TIME_T
 


c      FNAME='diablo.start'
      I = RANK+1

           FNAME ='last_saved/diablo_'
     &        //CHAR(MOD(I,100)/10+48)
     &        //CHAR(MOD(I,10)+48) // '.start'
           FNAME_TH ='last_saved/diablo_th01'
     &        //CHAR(MOD(I,100)/10+48)
     &        //CHAR(MOD(I,10)+48) // '.start'




c      DO I=1,RANK+1
c        FNAME_TH(m)='diablo_th'
c     &        //CHAR(MOD(N,100)/10+48)
c     &        //CHAR(MOD(N,10)+48) // '.start'
c      END DO

      WRITE(6,*)   'Reading flow from ',FNAME
      WRITE(6,*)   'Reading flow from ',FNAME_TH, CREATE_NEW_TH(1)

      OPEN(UNIT=10,FILE=FNAME,STATUS="OLD",FORM="UNFORMATTED")
      READ (10) NX_T, NY_T, NZ_T, NUM_PER_DIR_T, TIME, TIME_STEP


      write(*,*) 'NX_T, NY_T, NZ_T: ',NX_T,NY_T,NZ_T,NUM_PER_DIR_T,
     :             NUM_PER_DIR

      IF ((NX .NE. NX_T) .OR. (NY .NE. NY_T) .OR. (NZ .NE. NZ_T))
     *     STOP 'Error: old flowfield wrong dimensions. '
      IF (NUM_PER_DIR .NE. NUM_PER_DIR_T)
     *     STOP 'Error: old flowfield wrong NUM_PER_DIR. '

      write(*,*) 'READING FLOW'
      IF (NUM_PER_DIR.EQ.3) THEN
        READ (10) (((CU1(I,K,J),I=0,NKX),K=0,TNKZ),J=0,TNKY),
     *            (((CU2(I,K,J),I=0,NKX),K=0,TNKZ),J=0,TNKY),
     *            (((CU3(I,K,J),I=0,NKX),K=0,TNKZ),J=0,TNKY)
        DO N=1,NUM_READ_TH
! Specify in input.dat which scalars are to be read
          OPEN(UNIT=11,FILE=FNAME_TH,STATUS="OLD"
     &           ,FORM="UNFORMATTED")
          READ (11) NX_T, NY_T, NZ_T, NUM_PER_DIR_T, TIME_T, TIME_STEP_T
          READ (11) (((CTH(I,K,J,READ_TH_INDEX(N))
     &           ,I=0,NKX),K=0,TNKZ),J=0,TNKY)
         CLOSE(11)
        END DO
      ELSEIF (NUM_PER_DIR.EQ.2) THEN
        READ (10) (((CU1(I,K,J),I=0,NKX),K=0,TNKZ),J=0,NY+1),
     *            (((CU2(I,K,J),I=0,NKX),K=0,TNKZ),J=0,NY+1),
     *            (((CU3(I,K,J),I=0,NKX),K=0,TNKZ),J=0,NY+1),
     *            (((CP(I,K,J),I=0,NKX),K=0,TNKZ),J=0,NY+1)  
!        DO N=1,NUM_READ_TH
! Specify in input.dat which scalars are to be read
        IF (CREATE_NEW_TH(1)) THEN
        ELSE
          OPEN(UNIT=11,FILE=FNAME_TH,STATUS="OLD"
     &           ,FORM="UNFORMATTED")
          READ (11) NX_T, NY_T, NZ_T, NUM_PER_DIR_T, TIME_T,
     &              TIME_STEP_T
          READ (11) (((CTH(I,K,J,1),I=0,NKX),K=0,TNKZ),J=0,NY+1)
          READ (11) (((CTH(I,K,J,2),I=0,NKX),K=0,TNKZ),J=0,NY+1)
         CLOSE(11)
       ENDIF

      ELSEIF (NUM_PER_DIR.EQ.1) THEN
        READ (10) (((CU1(I,K,J),I=0,NKX),K=1,NZ),J=1,NY),
     *            (((CU2(I,K,J),I=0,NKX),K=1,NZ),J=2,NY),
     *            (((CU3(I,K,J),I=0,NKX),K=2,NZ),J=1,NY)
      ELSEIF (NUM_PER_DIR.EQ.0) THEN
        READ (10) (((U1(I,K,J),I=2,NX),K=1,NZ),J=1,NY),
     *            (((U2(I,K,J),I=1,NX),K=1,NZ),J=2,NY),
     *            (((U3(I,K,J),I=1,NX),K=2,NZ),J=1,NY)
      END IF
      CLOSE(10)
      CLOSE(11)

C Apply initial boundary conditions, set ghost cells
      IF (USE_MPI) THEN
        call APPLY_BC_VEL_MPI
      ELSE
        call APPLY_BC_VEL_LOWER
        call APPLY_BC_VEL_UPPER
      END IF

c For Scalar field end values need to be assigned 
c as boundary condition.
      
     
      DO J=0,NY+1
        TH_BACK(J,1)=Temp_0 ;
        TH_BACK(J,2)=Sal_0 ;
      END DO

       DO N=1,N_TH
        DO  K=0,TNKZ
         DO I=0,NKX 
          IF ((TH_BC_YMIN(N).EQ.0).AND.(TH_BC_YMAX(N).EQ.0)) THEN
           CTH(I,K,0,N) = CTH(I,K,1,N)
           CTH(I,K,NY+1,N) = CTH(I,K,NY,N)
          ENDIF
          
          IF ((TH_BC_YMIN(N).EQ.1).AND.(TH_BC_YMAX(N).EQ.1)) THEN
           CTH(I,K,0,N)    = CTH(I,K,1,N)-TH_BC_YMIN_C1(N)*DY(1)
           CTH(I,K,NY+1,N) = CTH(I,K,NY,N)+TH_BC_YMAX_C1(N)*DY(NY+1)
          ENDIF 

         ENDDO
        ENDDO
       ENDDO 
    
     


      RETURN
      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE SAVE_FLOW(FINAL)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      INCLUDE "mpif.h"
      INCLUDE 'header_mpi'

      CHARACTER*35 FNAME
      CHARACTER*35 FNAME_TH
      INTEGER      I, J, K, N
      LOGICAL      FINAL

      I = RANK +1 
      IF (FINAL) THEN
         FNAME='last_saved/diablo_'
     &        //CHAR(MOD(I,100)/10+48)
     &        //CHAR(MOD(I,10)+48) // '.res'


          FNAME_TH='last_saved/diablo_th01'
     &        //CHAR(MOD(I,100)/10+48)
     &        //CHAR(MOD(I,10)+48) // '.res'
      ELSE
         FNAME='last_saved/diablo_'
     &        //CHAR(MOD(I,100)/10+48)
     &        //CHAR(MOD(I,10)+48) // '.saved'
         FNAME_TH='last_saved/diablo_th01'
     &        //CHAR(MOD(I,100)/10+48)
     &        //CHAR(MOD(I,10)+48) // '.saved'

       IF (RANK .eq. 0) THEN 
       write(*,*) 'Ri value in the simulation', RI_TAU(1)
       ENDIF 
 
      END IF
      IF (RANK .eq. 0) THEN
      WRITE(6,*) 'Writing flow to ',FNAME
      ENDIF       

      OPEN(UNIT=10,FILE=FNAME,STATUS="UNKNOWN",FORM="UNFORMATTED")
      WRITE(10) NX, NY, NZ, NUM_PER_DIR, TIME, TIME_STEP


      IF (NUM_PER_DIR.EQ.3) THEN
        WRITE(10) (((CU1(I,K,J),I=0,NKX),K=0,TNKZ),J=0,TNKY),
     *            (((CU2(I,K,J),I=0,NKX),K=0,TNKZ),J=0,TNKY),
     *            (((CU3(I,K,J),I=0,NKX),K=0,TNKZ),J=0,TNKY)
        DO N=1,N_TH
          OPEN(UNIT=11,FILE=FNAME_TH,STATUS="UNKNOWN"
     &       ,FORM="UNFORMATTED")
          WRITE(11) NX, NY, NZ, NUM_PER_DIR, TIME, TIME_STEP
          WRITE(11) (((CTH(I,K,J,N),I=0,NKX),K=0,TNKZ),J=0,TNKY)
          CLOSE(11)
        END DO
      ELSEIF (NUM_PER_DIR.EQ.2) THEN
          WRITE(10) (((CU1(I,K,J),I=0,NKX),K=0,TNKZ),J=0,NY+1),
     *            (((CU2(I,K,J),I=0,NKX),K=0,TNKZ),J=0,NY+1),
     *            (((CU3(I,K,J),I=0,NKX),K=0,TNKZ),J=0,NY+1),
     *            (((CP(I,K,J),I=0,NKX),K=0,TNKZ),J=0,NY+1)
        !DO N=1,N_TH
          OPEN(UNIT=11,FILE=FNAME_TH,STATUS="UNKNOWN"
     &       ,FORM="UNFORMATTED")
          WRITE(11) NX, NY, NZ, NUM_PER_DIR, TIME, TIME_STEP
          WRITE(11) (((CTH(I,K,J,1),I=0,NKX),K=0,TNKZ),J=0,NY+1)
          WRITE(11) (((CTH(I,K,J,2),I=0,NKX),K=0,TNKZ),J=0,NY+1)
          CLOSE(11)
        !END DO
      ELSEIF (NUM_PER_DIR.EQ.1) THEN
        WRITE(10) (((CU1(I,K,J),I=0,NKX),K=1,NZ),J=1,NY),
     *            (((CU2(I,K,J),I=0,NKX),K=1,NZ),J=2,NY),
     *            (((CU3(I,K,J),I=0,NKX),K=2,NZ),J=1,NY)
      ELSEIF (NUM_PER_DIR.EQ.0) THEN
        WRITE(10) (((U1(I,K,J),I=2,NX),K=1,NZ),J=1,NY),
     *            (((U2(I,K,J),I=1,NX),K=1,NZ),J=2,NY),
     *            (((U3(I,K,J),I=1,NX),K=2,NZ),J=1,NY)
      END IF
      CLOSE(10)
      CLOSE(11)

      RETURN
      END

C----* --.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE SAVE_FLOW_INTER
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      INCLUDE "mpif.h"
      INCLUDE 'header_mpi'

      CHARACTER*35 FNAME
      CHARACTER*35 FNAME_TH
      INTEGER      I, J, K, N
      LOGICAL      FINAL

      I = RANK +1 
      j = TIME_STEP


    
         FNAME='save_int/diablo_'
     &        //CHAR(MOD(I,100)/10+48)
     &        //CHAR(MOD(I,10)+48) //'-'
     &        //CHAR(MOD(j,10000000)/1000000+48)
     &        //CHAR(MOD(j,1000000)/100000+48)
     &        //CHAR(MOD(j,100000)/10000+48)
     &        //CHAR(MOD(j,10000)/1000+48)
     &        //CHAR(MOD(j,1000)/100+48)
     &        //CHAR(MOD(j,100)/10+48)
     &        //CHAR(MOD(j,10)+48) //
     &         '.saved'
         FNAME_TH='save_int/diablo_th01'
     &        //CHAR(MOD(I,100)/10+48)
     &        //CHAR(MOD(I,10)+48) //'-'
     &        //CHAR(MOD(j,10000000)/1000000+48)
     &        //CHAR(MOD(j,1000000)/100000+48)
     &        //CHAR(MOD(j,100000)/10000+48)
     &        //CHAR(MOD(j,10000)/1000+48)
     &        //CHAR(MOD(j,1000)/100+48)
     &        //CHAR(MOD(j,100)/10+48)
     &        //CHAR(MOD(j,10)+48) //
     &           '.saved'

      IF (RANK .eq. 0)THEN
      write(*,*) 'Ri value in the simulation', RI_TAU(1)

      
      WRITE(6,*) 'Writing flow to ',FNAME
      ENDIF

      OPEN(UNIT=10,FILE=FNAME,STATUS="UNKNOWN",FORM="UNFORMATTED")
      WRITE(10) NX, NY, NZ, NUM_PER_DIR, TIME, TIME_STEP


      IF (NUM_PER_DIR.EQ.3) THEN
        WRITE(10) (((CU1(I,K,J),I=0,NKX),K=0,TNKZ),J=0,TNKY),
     *            (((CU2(I,K,J),I=0,NKX),K=0,TNKZ),J=0,TNKY),
     *            (((CU3(I,K,J),I=0,NKX),K=0,TNKZ),J=0,TNKY)
        DO N=1,N_TH
          OPEN(UNIT=11,FILE=FNAME_TH,STATUS="UNKNOWN"
     &       ,FORM="UNFORMATTED")
          WRITE(11) NX, NY, NZ, NUM_PER_DIR, TIME, TIME_STEP
          WRITE(11) (((CTH(I,K,J,N),I=0,NKX),K=0,TNKZ),J=0,TNKY)
          CLOSE(11)
        END DO
      ELSEIF (NUM_PER_DIR.EQ.2) THEN
          WRITE(10) (((CU1(I,K,J),I=0,NKX),K=0,TNKZ),J=1,NY),
     *            (((CU2(I,K,J),I=0,NKX),K=0,TNKZ),J=2,NY),
     *            (((CU3(I,K,J),I=0,NKX),K=0,TNKZ),J=1,NY)
        DO N=1,N_TH
          OPEN(UNIT=11,FILE=FNAME_TH,STATUS="UNKNOWN"
     &       ,FORM="UNFORMATTED")
          WRITE(11) NX, NY, NZ, NUM_PER_DIR, TIME, TIME_STEP
          WRITE(11) (((CTH(I,K,J,N),I=0,NKX),K=0,TNKZ),J=1,NY)
          CLOSE(11)
        END DO
      ELSEIF (NUM_PER_DIR.EQ.1) THEN
        WRITE(10) (((CU1(I,K,J),I=0,NKX),K=1,NZ),J=1,NY),
     *            (((CU2(I,K,J),I=0,NKX),K=1,NZ),J=2,NY),
     *            (((CU3(I,K,J),I=0,NKX),K=2,NZ),J=1,NY)
      ELSEIF (NUM_PER_DIR.EQ.0) THEN
        WRITE(10) (((U1(I,K,J),I=2,NX),K=1,NZ),J=1,NY),
     *            (((U2(I,K,J),I=1,NX),K=1,NZ),J=2,NY),
     *            (((U3(I,K,J),I=1,NX),K=2,NZ),J=1,NY)
      END IF
      CLOSE(10)
      CLOSE(11)

      RETURN
      END

C----*
!C--.---------.---------.---------.---------.---------.---------.-|-------|
!      SUBROUTINE SAVE_BU(TIME_STEP)
!C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
!      INCLUDE 'header'
!      INCLUDE "mpif.h"
!      INCLUDE 'header_mpi'
!
!      CHARACTER*35 FNAME
!      CHARACTER*35 FNAME_TH
!!      CHARACTER*35 DIR
!      INTEGER      I, J, K, N
!      LOGICAL      FINAL
!
!      I = RANK +1
!      j = TIME_STEP
 
!        DIR=write (i, "(A11, I4)") ,"last_saved_",j    

!         FNAME='save_int/diablo_'
!     &        //CHAR(MOD(I,100)/10+48)
!     &        //CHAR(MOD(I,10)+48) //'-'

!     &        //CHAR(MOD(j,10000000)/1000000+48)

!     &        //CHAR(MOD(j,1000000)/100000+48)
!     &        //CHAR(MOD(j,100000)/10000+48)
!     &        //CHAR(MOD(j,10000)/1000+48)
!     &        //CHAR(MOD(j,1000)/100+48)
!     &        //CHAR(MOD(j,100)/10+48)
!     &        //CHAR(MOD(j,10)+48) //
!     &         '.saved'
!         FNAME_TH='save_int/diablo_th01'
!     &        //CHAR(MOD(I,100)/10+48)
!     &        //CHAR(MOD(I,10)+48) //'-'
!     &        //CHAR(MOD(j,10000000)/1000000+48)
!     &        //CHAR(MOD(j,1000000)/100000+48)
!     &        //CHAR(MOD(j,100000)/10000+48)
!    &        //CHAR(MOD(j,10000)/1000+48)
!     &        //CHAR(MOD(j,1000)/100+48)
!     &        //CHAR(MOD(j,100)/10+48)
!     &        //CHAR(MOD(j,10)+48) //
!     &           '.saved'

!      OPEN(UNIT=10,FILE=FNAME,STATUS="UNKNOWN",FORM="UNFORMATTED")
!      WRITE(10) NX, NY, NZ, NUM_PER_DIR, TIME, TIME_STEP


!      IF (NUM_PER_DIR.EQ.3) THEN
!        WRITE(10) (((CU1(I,K,J),I=0,NKX),K=0,TNKZ),J=0,TNKY),
!     *            (((CU2(I,K,J),I=0,NKX),K=0,TNKZ),J=0,TNKY),
!     *            (((CU3(I,K,J),I=0,NKX),K=0,TNKZ),J=0,TNKY)
!        DO N=1,N_TH
!          OPEN(UNIT=11,FILE=FNAME_TH,STATUS="UNKNOWN"
!     &       ,FORM="UNFORMATTED")
!          WRITE(11) NX, NY, NZ, NUM_PER_DIR, TIME, TIME_STEP
!          WRITE(11) (((CTH(I,K,J,N),I=0,NKX),K=0,TNKZ),J=0,TNKY)
!          CLOSE(11)
!        END DO
!      ELSEIF (NUM_PER_DIR.EQ.2) THEN
!          WRITE(10) (((CU1(I,K,J),I=0,NKX),K=0,TNKZ),J=1,NY),
!     *            (((CU2(I,K,J),I=0,NKX),K=0,TNKZ),J=2,NY),
!     *            (((CU3(I,K,J),I=0,NKX),K=0,TNKZ),J=1,NY)
!        DO N=1,N_TH
!          OPEN(UNIT=11,FILE=FNAME_TH,STATUS="UNKNOWN"
!     &       ,FORM="UNFORMATTED")
!          WRITE(11) NX, NY, NZ, NUM_PER_DIR, TIME, TIME_STEP
!          WRITE(11) (((CTH(I,K,J,N),I=0,NKX),K=0,TNKZ),J=1,NY)
!          CLOSE(11)
!        END DO
!      ELSEIF (NUM_PER_DIR.EQ.1) THEN
!        WRITE(10) (((CU1(I,K,J),I=0,NKX),K=1,NZ),J=1,NY),
!     *            (((CU2(I,K,J),I=0,NKX),K=1,NZ),J=2,NY),
!     *            (((CU3(I,K,J),I=0,NKX),K=2,NZ),J=1,NY)
!      ELSEIF (NUM_PER_DIR.EQ.0) THEN
!        WRITE(10) (((U1(I,K,J),I=2,NX),K=1,NZ),J=1,NY),
!     *            (((U2(I,K,J),I=1,NX),K=1,NZ),J=2,NY),
!     *            (((U3(I,K,J),I=1,NX),K=2,NZ),J=1,NY)
!      END IF
!      CLOSE(10)
!      RETURN
!      END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE FILTER(n)
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'

      integer n

      IF (NUM_PER_DIR.EQ.3) THEN
        CALL FILTER_PER
      ELSEIF (NUM_PER_DIR.EQ.2) THEN
        CALL FILTER_CHAN(n)
      END IF

      RETURN
      END

