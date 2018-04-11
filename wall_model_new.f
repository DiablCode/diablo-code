      SUBROUTINE WALL_MODEL_LOWER
! This subroutine creates an artificial slip condition for the velocity
! field using a modification of the wall model by Schumann (1975) and
! Grotzbach (1987) to account for a constant angle between the outer
! flow and the surface stress owing to system rotation (for an Ekman layer)
! The velocity field should be in **PHYSICAL** space upon calling
! Returns the wall stress in U_BC_LOWER_NWM(0:NXM,0:NZM), etc.

      include 'header'

      integer i,j,k

! alpha0 is the angle of surface stress in degrees
! For a nonrotating flow, set = 0
      real*8 alpha_0
!      parameter (alpha_0=23.45d0)
      parameter (alpha_0=0.d0)

! Constants for the log-law equation
      real*8 ikappa,B,z0
! Smooth wall:
      parameter (ikappa=2.44d0,B=5.2d0,z0=0.d0)     
! Rough wall
!      parameter (ikappa=2.44d0,B=0.d0,z0=2.3422d-6)     

! Location at which to apply log law (in wall units)
      real*8 zplus

      real*8 U1_mean,U3_mean,U_mean
      real*8 tauwall_mean_x,tauwall_mean_z
      real*8 tauwall_x(0:NX+1,0:NZ+1)
      real*8 tauwall_z(0:NX+1,0:NZ+1)
      real*8 tauwall_z_test
      real*8 coef1,coef2
      real*8 coef_x1,coef_x2,coef_z1,coef_z2
      integer index_x1,index_x2,index_z1,index_z2
      real*8 shift_amp,shift_x,shift_z

! Variables for Newton-Raphson solver
      real*8 f,fp
      real*8 res,error
      integer iterations

! Get interpolation coefficients at location of log law
! For smooth walls      
!      zplus=20.d0
! For rough wall
!      zplus=GYF(1)/z0
!      do j=1,NY 
!        if ((gyf(j)-gy(2))/NU.lt.zplus) j1=j
!      end do
!      j2=j1+1
! Use linear interpolation to get coefficients
! For open channel flow
!      coef1=(gyf(j2)-zplus*NU)/(gyf(j2)-gyf(j1))
! For closed chanel flow
!      coef1=((gyf(j2)-gy(2))-zplus*NU)/(gyf(j2)-gyf(j1))
!      coef2=1.d0-coef1

! Or, use the first point, instead of interpolated value
      coef1=1.d0
      coef2=0.d0
      j1=2
      j2=1
! For smooth wall
      zplus=(gyf(2)-gy(2))/NU
! For rough wall
!      zplus=(gyf(2)-gy(2))/z0

! Get the plane average of the horizontal velocity at zplus     
      U1_mean=(coef1*SUM(U1(0:NXM,0:NZM,j1))
     &       +coef2*SUM(U1(0:NXM,0:NZM,j2)))
     &             /dble(NX*NZ)
      U3_mean=(coef1*SUM(U3(0:NXM,0:NZM,j1))
     &       +coef2*SUM(U3(0:NXM,0:NZM,j2)))
     &             /dble(NX*NZ)
      U_mean=sqrt(U1_mean**2.d0+U3_mean**2.d0)      

! Smooth wall
! Make a first guess at the mean utau 
      utau_mean_lower=U_mean         
      iterations=0
      error=1.d-9
      res=1. 
      do while (res.ge.error.and.iterations.lt.50)
        iterations=iterations+1
        utau_mean_lower=abs(utau_mean_lower)
        f=ikappa*utau_mean_lower*log(zplus*utau_mean_lower)
     &          +B*utau_mean_lower-U_mean
        fp=ikappa*log(zplus*utau_mean_lower)+ikappa+B
        res=abs(f/sqrt(utau_mean_lower))
        utau_mean_lower=utau_mean_lower-f/fp
      end do
      if (res.gt.error) then
        write(*,*) 'WARNING: Failed convergence in wall_model'
      end if
! Schumann model
!      utau_mean_lower=1.0d0
! Rough wall
!      utau_mean_lower=U_mean/(ikappa*log(zplus))

      tauwall_mean=utau_mean_lower**2.d0


      tauwall_mean_x=tauwall_mean*cos(alpha_0*2*pi/360.d0)
      tauwall_mean_z=tauwall_mean*sin(alpha_0*2*pi/360.d0)

! Redefine the mean velocity to be at the first gridpoint gyf(2)
      U1_mean=SUM(U1(0:NXM,0:NZM,2))
     &          /dble(NX*NZ)
      U3_mean=SUM(U3(0:NXM,0:NZM,2))
     &         /dble(NX*NZ)

      do k=0,NZM
        do i=0,NXM
! Modified Schumann-Grotzbach model:
          tauwall_x(i,k)=U1(i,k,2)*tauwall_mean_x/U1_mean
! For Ekman layer
!          tauwall_z(i,k)=U3(i,k,2)*tauwall_mean_z/U3_mean
!          tauwall_z_test=U3(i,k,2)*tauwall_mean_x/U1_mean
! Take max of rotating and nonrotating versions
!          if (abs(tauwall_z_test).gt.abs(tauwall_z(i,k))) then
!             tauwall_z(i,k)=tauwall_z_test
!          end if
! For nonrotating
          tauwall_z(i,k)=U3(i,k,2)*tauwall_mean_x/U1_mean


! MKP Model, no shift
!           tauwall_x(i,k)=tauwall_mean*U1_mean
!     &             /sqrt(U1_mean**2.d0+U3_mean**2.d0)
!           tauwall_z(i,k)=tauwall_mean*U3_mean
!     &             /sqrt(U1_mean**2.d0+U3_mean**2.d0)

! MKP Model, with stress-wise shift

!      shift_amp=GYF(2)*1.d0/tan(13*2*PI/360.d0)
!      shift_x=shift_amp*U1_mean/sqrt(U1_mean**2.d0+U3_mean**2.d0) 
!      shift_z=shift_amp*U3_mean/sqrt(U1_mean**2.d0+U3_mean**2.d0) 

!      index_x1=int(shift_x/DX(1))
!      index_x2=index_x1+1
!      index_z1=int(shift_z/DZ(1))
!      index_z2=index_z1+1

!      coef_x1=index_x2-shift_x/DX(1) 
!      coef_x2=1.d0-coef_x1 
!      coef_z1=index_z2-shift_z/DZ(1) 
!      coef_z2=1.d0-coef_z1 

!      tauwall_x(i,k)=tauwall_mean*U1_mean
!     &                /sqrt(U1_mean**2.d0+U3_mean**2.d0)
!     &       -0.1d0*utau_mean_lower
!     &        *(coef_x1*coef_z1*U1(index_x1,index_z1,2)
!     &         +coef_x1*coef_z2*U1(index_x1,index_z2,2)
!     &         +coef_x2*coef_z1*U1(index_x2,index_z1,2)
!     &         +coef_x2*coef_z2*U1(index_x2,index_z2,2)-U1_mean)
!      tauwall_z(i,k)=tauwall_mean*U3_mean
!     &             /sqrt(U1_mean**2.d0+U3_mean**2.d0)
!     &       -0.1d0*utau_mean_lower
!     &        *(coef_x1*coef_z1*U3(index_x1,index_z1,2)
!     &         +coef_x1*coef_z2*U3(index_x1,index_z2,2)
!     &         +coef_x2*coef_z1*U3(index_x2,index_z1,2)
!     &         +coef_x2*coef_z2*U3(index_x2,index_z2,2)-U3_mean)


! Ok, now all models should do the following:

! Set the velocity gradient at the wall based on the local wall stress
          U_BC_LOWER_NWM(i,k)=tauwall_x(i,k)/NU
          W_BC_LOWER_NWM(i,k)=tauwall_z(i,k)/NU


        end do 
      end do 

      write(200,*) U_mean,utau_mean_lower
     &    ,SUM(U_BC_LOWER_NWM(0:NXM,0:NZM))/dble(NX*NZ)

      return
      end 






      SUBROUTINE WALL_MODEL_UPPER
! This subroutine creates an artificial slip condition for the velocity
! field using a modification of the wall model by Schumann (1975) and
! Grotzbach (1987) to account for a constant angle between the outer
! flow and the surface stress owing to system rotation (for an Ekman layer)
! The velocity field should be in **PHYSICAL** space upon calling
! Returns the wall stress in U_BC_UPPER_NWM(0:NXM,0:NZM), etc.

      include 'header'

      integer i,k,N

! alpha0 is the angle of surface stress in degrees
! For a nonrotating flow, set = 0
      real*8 alpha_0
      parameter (alpha_0=0.d0)

! Constants for the log-law equation
      real*8 ikappa,B
      parameter (ikappa=2.44d0,B=5.2d0)     

      real*8 U1_mean,U3_mean,U_mean
      real*8 tauwall_mean_x,tauwall_mean_z
      real*8 tauwall_x(0:NX+1,0:NZ+1)
      real*8 tauwall_z(0:NX+1,0:NZ+1)
      real*8 coef1,coef2
      real*8 coef_x1,coef_x2,coef_z1,coef_z2
      integer index_x1,index_x2,index_z1,index_z2
      real*8 shift_amp,shift_x,shift_z

! Variables for Newton-Raphson solver
      real*8 f,fp
      real*8 res,error
      integer iterations

! Get the plane average of the horizontal velocity at GYF(NY-1)     
      U1_mean=SUM(U1(0:NXM,0:NZM,NY-1))/dble(NX*NZ)
      U3_mean=SUM(U3(0:NXM,0:NZM,NY-1))/dble(NX*NZ)
      U_mean=sqrt(U1_mean**2.d0+U3_mean**2.d0)      
 
! Make a first guess at the mean utau 
      utau_mean_upper=U_mean         
          
      iterations=0
      error=1.d-9
      res=1. 
      do while (res.ge.error.and.iterations.lt.50)
        iterations=iterations+1
        utau_mean_upper=abs(utau_mean_upper)
        f=ikappa*utau_mean_upper
     &         *log(abs((GY(NY)-GYF(NY-1))/NU)*utau_mean_upper)
     &          +B*utau_mean_upper-U_mean
        fp=ikappa*log(abs((GY(NY)-GYF(NY-1))/NU)*utau_mean_upper)
     &           +ikappa+B
        res=abs(f/sqrt(utau_mean_upper))
        utau_mean_upper=utau_mean_upper-f/fp
      end do
      if (res.gt.error) then
        write(*,*) 'WARNING: Failed convergence in wall_model'
      end if
      tauwall_mean=utau_mean_upper**2.d0

      tauwall_mean_x=tauwall_mean*cos(alpha_0*2*pi/360.d0)
      tauwall_mean_z=tauwall_mean*sin(alpha_0*2*pi/360.d0)

      do k=0,NZM
        do i=0,NXM
! SG Model
          tauwall_x(i,k)=U1(i,k,NY-1)*tauwall_mean_x/U1_mean
! Rotating version
!          tauwall_z(i,k)=U3(i,k,NY-1)*tauwall_mean_z/U3_mean
! Non-rotating version
          tauwall_z(i,k)=U3(i,k,NY-1)*tauwall_mean_x/U1_mean
       
 
! MKP Model, with stress-wise shift

!      shift_amp=GYF(2)*1.d0/tan(13*2*PI/360.d0)
!      shift_x=shift_amp*U1_mean/sqrt(U1_mean**2.d0+U3_mean**2.d0) 
!      shift_z=shift_amp*U3_mean/sqrt(U1_mean**2.d0+U3_mean**2.d0) 

!      index_x1=int(shift_x/DX(1))
!      index_x2=index_x1+1
!      index_z1=int(shift_z/DZ(1))
!      index_z2=index_z1+1

!      coef_x1=index_x2-shift_x/DX(1) 
!      coef_x2=1.d0-coef_x1 
!      coef_z1=index_z2-shift_z/DZ(1) 
!      coef_z2=1.d0-coef_z1 

!      tauwall_x(i,k)=tauwall_mean*U1_mean
!     &                /sqrt(U1_mean**2.d0+U3_mean**2.d0)
!     &       -0.1d0*utau_mean_upper
!     &        *(coef_x1*coef_z1*U1(index_x1,index_z1,NY-1)
!     &         +coef_x1*coef_z2*U1(index_x1,index_z2,NY-1)
!     &         +coef_x2*coef_z1*U1(index_x2,index_z1,NY-1)
!     &         +coef_x2*coef_z2*U1(index_x2,index_z2,NY-1)-U1_mean)
!      tauwall_z(i,k)=tauwall_mean*U3_mean
!     &             /sqrt(U1_mean**2.d0+U3_mean**2.d0)
!     &       -0.1d0*utau_mean_upper
!     &        *(coef_x1*coef_z1*U3(index_x1,index_z1,NY-1)
!     &         +coef_x1*coef_z2*U3(index_x1,index_z2,NY-1)
!     &         +coef_x2*coef_z1*U3(index_x2,index_z1,NY-1)
!     &         +coef_x2*coef_z2*U3(index_x2,index_z2,NY-1)-U3_mean)

! Set the velocity gradient at the wall based on the local wall stress
          U_BC_UPPER_NWM(i,k)=-1.d0*tauwall_x(i,k)/NU
          W_BC_UPPER_NWM(i,k)=-1.d0*tauwall_z(i,k)/NU

! Use a dirichlet version of the wall model
        end do 
      end do 

      write(201,*) U_mean,utau_mean_upper,U1(10,10,NY-1)
     &      ,U_BC_UPPER_NWM(10,10),(GY(NY)-GYF(NY-1))/NU,GY(NY)
     &      ,GYF(NY-1),NU

      return
      end 


      subroutine apply_filter_vel(CVEL,js,je)
C This subroutine applies a filter to the highest wavenumbers
C It should be applied with the velcoity in Fourier space
C The filter used is a sharpened raised cosine filter in the horizontal
C and a fourth order implicit compact filter in the vertical, with the
C parameter alpha determining the width of the vertical filtering window

      include 'header'

      COMPLEX*16 CVEL(0:NX/2,0:NZ+1,0:NY+1) 

      integer I,J,K,js,je,n

! Variables for horizontal filtering
      real*8 sigma(0:NKX,0:TNKZ),sigma0

! Variables for vertical filtering
      real*8 alpha
      parameter (alpha=0.48d0)
! Parameters for a larger stencil filter
      real*8 f_a,f_b,f_c

C Set the filtering constants for the horizontal direction
!      DO i=0,NKX
!       DO k=0,TNKZ
!        sigma0=0.5d0*(1.d0+
!     &       cos(sqrt((KX(i)*LX*1.d0/float(NX))**2.d0
!     &            +(KZ(k)*LZ*1.d0/float(NZ))**2.d0)))
! Apply a sharpened raised cosine filter
!        sigma(i,k)=sigma0**4.d0*(35.d0-84.d0*sigma0
!     &        +70.d0*sigma0**2.d0-20.d0*sigma0**3.d0)
!       END DO
!      END DO

C Do the spectral filtering in the horizontal
!        DO K=0,TNKZ
!          DO I=0,NKX
!            if ((I.ne.0).or.(K.ne.0)) then
!            DO J=js+1,je-1
!              CVEL(I,K,J)=CVEL(I,K,J)*sigma(i,k)
!            END DO
!            end if
!          END DO 
!        END DO

C Set the filtering constants
      f_a=(1.d0/8.d0)*(5.d0+6.d0*alpha)
      f_b=0.5d0*(1.d0+2.d0*alpha)
      f_c=(-1.d0/8.d0)*(1.d0-2.d0*alpha)      


C First, zero the tridiagonal matrix components
      DO I=0,NKX
        DO J=1,NY
          MATD_C(I,J)=1.d0
          MATL_C(I,J)=0.d0
          MATU_C(I,J)=0.d0
          VEC_C(I,J)=0.d0
        END DO
      END DO  

      IF (MOD(TIME_STEP,FILTER_VEL_INT).eq.0) THEN

C Filter the velocity in the vertical direction
      DO K=0,TNKZ
        DO I=0,NKX
C Construct the centered difference terms
          DO J=3,NY-2
            MATL_C(I,J)=alpha
            MATD_C(I,J)=1.d0
            MATU_C(I,J)=alpha
            VEC_C(I,J)=f_a*CVEL(I,K,J)
     &                +(f_b/2.d0)*(CVEL(I,K,J+1)+CVEL(I,K,J-1))
     &                +(f_c/2.d0)*(CVEL(I,K,J+2)+CVEL(I,K,J-2))
          END DO
C Now, construct the equations for the boundary nodes
C Filter explicitly with a fourth-order accurate scheme
          J=2
            MATL_C(I,J)=0.d0
            MATD_C(I,J)=1.d0
            MATU_C(I,J)=0.d0
            VEC_C(I,J)=(3.d0/4.d0)*CVEL(I,K,2)
     &          +(1.d0/16.d0)*(1.d0*CVEL(I,K,1)+6.d0*CVEL(I,K,3)
     &                        -4.d0*CVEL(I,K,4)+1.d0*CVEL(I,K,5))
          J=1
            MATL_C(I,J)=0.d0
            MATD_C(I,J)=1.d0
            MATU_C(I,J)=0.d0
            VEC_C(I,J)=(15.d0/16.d0)*CVEL(I,K,1)
     &          +(1.d0/16.d0)*(4.d0*CVEL(I,K,2)-6.d0*CVEL(I,K,3)
     &                     +4.d0*CVEL(I,K,4)-1.d0*CVEL(I,K,5))
          J=NY
            MATL_C(I,J)=0.d0
            MATD_C(I,J)=1.d0
            MATU_C(I,J)=0.d0
            VEC_C(I,J)=(15.d0/16.d0)*CVEL(I,K,NY)
     &          +(1.d0/16.d0)*(4.d0*CVEL(I,K,NY-1)-6.d0*CVEL(I,K,NY-2)
     &                     +4.d0*CVEL(I,K,NY-3)-1.d0*CVEL(I,K,NY-4))
          J=NY-1
            MATL_C(I,J)=0.d0
            MATD_C(I,J)=1.d0
            MATU_C(I,J)=0.d0
            VEC_C(I,J)=(3.d0/4.d0)*CVEL(I,K,NY-1)
     &          +(1.d0/16.d0)*(1.d0*CVEL(I,K,NY)+6.d0*CVEL(I,K,NY-2)
     &                        -4.d0*CVEL(I,K,NY-3)+1.d0*CVEL(I,K,NY-4))
         END DO
C Now, solve the tridiagonal system
         CALL THOMAS_COMPLEX(MATL_C,MATD_C,MATU_C,VEC_C,NY,NKX)
         DO I=0,NKX
           DO J=js+1,je-1
             CVEL(I,K,J)=VEC_C(I,J)
           END DO
         END DO
C END DO K  
       END DO

       end if

       return
       end 

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE APPLY_BC_NWM
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C This subroutine linearly interpolates the velocity to the ghost cells
C For use when a wall model is used.  This should be applied before calling
C les_chan
C It is important not to use the wall gradient in calculating the SGS stress
      INCLUDE 'header'
      INTEGER I,K,N

C If we are using the near wall model, first convert the local boundary 
C condition to Fourier Space
      IF (U_BC_LOWER.eq.3) THEN
C Store the boundary condition in one plane of CS1 which is sized
C correctly for the FFT calls
        DO K=0,NZM
          DO I=0,NXM
            S1(I,K,0)=U_BC_LOWER_NWM(I,K)
          END DO
        END DO
        call FFT_XZ_TO_FOURIER(S1,CS1,0,0)
        DO K=0,TNKZ
          DO I=0,NKX
            CU_BC_LOWER_NWM(I,K)=CS1(I,K,0)
          END DO
        END DO
      END IF
      IF (U_BC_UPPER.eq.3) THEN
C Store the boundary condition in one plane of CS1 which is sized
C correctly for the FFT calls
        DO K=0,NZM
          DO I=0,NXM
            S1(I,K,0)=U_BC_UPPER_NWM(I,K)
          END DO
        END DO
        call FFT_XZ_TO_FOURIER(S1,CS1,0,0)
        DO K=0,TNKZ
          DO I=0,NKX
            CU_BC_UPPER_NWM(I,K)=CS1(I,K,0)
          END DO
        END DO
      END IF
      IF (W_BC_LOWER.eq.3) THEN
C Store the boundary condition in one plane of CS1 which is sized
C correctly for the FFT calls
        DO K=0,NZM
          DO I=0,NXM
            S1(I,K,0)=W_BC_LOWER_NWM(I,K)
          END DO
        END DO
        call FFT_XZ_TO_FOURIER(S1,CS1,0,0)
        DO K=0,TNKZ
          DO I=0,NKX
            CW_BC_LOWER_NWM(I,K)=CS1(I,K,0)
          END DO
        END DO
      END IF
      IF (W_BC_UPPER.eq.3) THEN
C Store the boundary condition in one plane of CS1 which is sized
C correctly for the FFT calls
        DO K=0,NZM
          DO I=0,NXM
            S1(I,K,0)=W_BC_UPPER_NWM(I,K)
          END DO
        END DO
        call FFT_XZ_TO_FOURIER(S1,CS1,0,0)
        DO K=0,TNKZ
          DO I=0,NKX
            CW_BC_UPPER_NWM(I,K)=CS1(I,K,0)
          END DO
        END DO
      END IF

       IF (U_BC_LOWER.EQ.3) THEN
        DO K=0,TNKZ
        DO I=0,NKX
! Linear extrapolation of resolved velocity
          CU1(I,K,1)=0.5d0*(CU1(I,K,3)+CU1(I,K,2))
     &            -((CU1(I,K,3)-CU1(I,K,2))/DY(3))*(GY(3)-GYF(1))
          CU1(I,K,0)=0.5d0*(CU1(I,K,3)+CU1(I,K,2))
     &            -((CU1(I,K,3)-CU1(I,K,2))/DY(3))*(GY(3)-GYF(0))
! One half of the imposed gradient
!         CU1(I,K,1)=CU1(I,K,2)-0.5d0*DY(2)*CU_BC_LOWER_NWM(I,K)
!         CU1(I,K,0)=CU1(I,K,2)-0.5d0*(DY(2)+DY(1))*CU_BC_LOWER_NWM(I,K)
        END DO
        END DO
       END IF
       IF (U_BC_UPPER.EQ.3) THEN
        DO K=0,TNKZ
        DO I=0,NKX
! Linear extrapolation of resolved velocity
          CU1(I,K,NY)=0.5d0*(CU1(I,K,NY-2)+CU1(I,K,NY-1))
     &     +((CU1(I,K,NY-1)-CU1(I,K,NY-2))/DY(NY-2))*(GYF(NY)-GY(NY-1))
          CU1(I,K,NY+1)=0.5d0*(CU1(I,K,NY-2)+CU1(I,K,NY-1))
     &    +((CU1(I,K,NY-1)-CU1(I,K,NY-2))/DY(NY-2))*(GYF(NY+1)-GY(NY-1))
        END DO
        END DO
       END IF

       IF (W_BC_LOWER.EQ.3) THEN
        DO K=0,TNKZ
        DO I=0,NKX
! Linear extrapolation of resolved velocity
          CU3(I,K,1)=0.5d0*(CU3(I,K,3)+CU3(I,K,2))
     &            -((CU3(I,K,3)-CU3(I,K,2))/DY(3))*(GY(3)-GYF(1))
          CU3(I,K,0)=0.5d0*(CU3(I,K,3)+CU3(I,K,2))
     &            -((CU3(I,K,3)-CU3(I,K,2))/DY(3))*(GY(3)-GYF(0))
! One half of the imposed gradient
!         CU3(I,K,1)=CU3(I,K,2)-0.5d0*DY(2)*CW_BC_LOWER_NWM(I,K)
!         CU3(I,K,0)=CU3(I,K,2)-0.5d0*(DY(2)+DY(1))*CW_BC_LOWER_NWM(I,K)

        END DO
        END DO
       END IF
       IF (W_BC_UPPER.EQ.3) THEN
        DO K=0,TNKZ
        DO I=0,NKX
! Linear extrapolation of resolved velocity
          CU3(I,K,NY)=0.5d0*(CU3(I,K,NY-2)+CU3(I,K,NY-1))
     &     +((CU3(I,K,NY-1)-CU3(I,K,NY-2))/DY(NY-2))*(GYF(NY)-GY(NY-1))
          CU3(I,K,NY+1)=0.5d0*(CU3(I,K,NY-2)+CU3(I,K,NY-1))
     &    +((CU3(I,K,NY-1)-CU3(I,K,NY-2))/DY(NY-2))*(GYF(NY+1)-GY(NY-1))
        END DO
        END DO
       END IF

       IF (V_BC_LOWER.EQ.3) THEN
        DO K=0,TNKZ
        DO I=0,NKX
          CU2(I,K,2)=0.d0
          CU2(I,K,1)=0.5d0*(CU2(I,K,3))
     &            -((CU2(I,K,3))/DYF(2))*(GYF(2)-GY(1))
          CU2(I,K,0)=0.5d0*(CU2(I,K,3))
     &            -((CU2(I,K,3))/DYF(2))*(GYF(2)-GY(0))
        END DO
        END DO
       END IF
       IF (V_BC_UPPER.EQ.3) THEN
        DO K=0,TNKZ
        DO I=0,NKX
          CU2(I,K,NY)=0.d0
          CU2(I,K,NY+1)=0.5d0*(CU2(I,K,NY-1))
     &            +((0.d0-CU2(I,K,NY-1))/DYF(NY-1))*(GY(NY)-GYF(NY-1))
        END DO
        END DO
       END IF
     
       IF ((U_BC_LOWER.EQ.3).OR.(W_BC_LOWER.EQ.3)
     &     .OR.(V_BC_LOWER.EQ.3)) THEN
       DO N=1,N_TH
        DO K=0,TNKZ
        DO I=0,NKX
          CTH(I,K,1,N)=0.5d0*(CTH(I,K,3,N)+CTH(I,K,2,N))
     &            -((CTH(I,K,3,N)-CTH(I,K,2,N))/DY(3))*(GY(3)-GYF(1))
          CTH(I,K,0,N)=0.5d0*(CTH(I,K,3,N)+CTH(I,K,2,N))
     &            -((CTH(I,K,3,N)-CTH(I,K,2,N))/DY(3))*(GY(3)-GYF(0))
        END DO
        END DO
       END DO
       END IF
       IF ((U_BC_UPPER.EQ.3).OR.(W_BC_UPPER.EQ.3)
     &     .OR.(V_BC_UPPER.EQ.3)) THEN
       DO N=1,N_TH
        DO K=0,TNKZ
        DO I=0,NKX
          CTH(I,K,NY,N)=0.5d0*(CTH(I,K,NY-2,N)+CTH(I,K,NY-1,N))
     &        +((CTH(I,K,3,N)-CTH(I,K,NY-1,N))/DY(3))*(GY(NY-2)-GYF(NY))
          CTH(I,K,NY+1,N)=0.5d0*(CTH(I,K,NY-2,N)+CTH(I,K,NY-1,N))
     &          +((CTH(I,K,NY-2,N)-CTH(I,K,NY-1,N))
     &           /DY(NY-2))*(GY(NY-2)-GYF(NY+1))
        END DO
        END DO
       END DO
       END IF
         
 
       RETURN
       END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE APPLY_BC_ANWM
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C This subroutine linearly interpolates the velocity to the ghost cells
C For use when a wall model is used.  This should be applied before calling
C les_chan
C It is important not to use the wall gradient in calculating the SGS stress
      INCLUDE 'header'
      INTEGER I,K

       IF (U_BC_LOWER.EQ.3) THEN
        DO K=0,NZM
        DO I=0,NXM
          U1(I,K,1)=U1(I,K,2)-DY(2)*U_BC_LOWER_NWM(I,K)
          U1(I,K,0)=U1(I,K,1)
        END DO
        END DO
       END IF
       IF (U_BC_UPPER.EQ.3) THEN
        DO K=0,NZM
        DO I=0,NXM
          U1(I,K,NY)=U1(I,K,NY-1)+DY(NY)*U_BC_UPPER_NWM(I,K)
          U1(I,K,NY+1)=U1(I,K,NY)
        END DO
        END DO
       END IF

       IF (W_BC_LOWER.EQ.3) THEN
        DO K=0,NZM
        DO I=0,NXM
          U3(I,K,1)=U3(I,K,2)-DY(2)*W_BC_LOWER_NWM(I,K)
          U3(I,K,0)=U3(I,K,1)
        END DO
        END DO
       END IF
       IF (W_BC_UPPER.EQ.3) THEN
        DO K=0,NZM
        DO I=0,NXM
          U3(I,K,NY)=U3(I,K,NY-1)+DY(NY)*W_BC_UPPER_NWM(I,K)
          U3(I,K,NY+1)=U3(I,K,NY)
        END DO
        END DO
       END IF

       IF (V_BC_LOWER.EQ.3) THEN
        DO K=0,NZM
        DO I=0,NXM
          U2(I,K,1)=0.d0
          U2(I,K,2)=0.0d0
        END DO
        END DO
       END IF
       IF (V_BC_UPPER.EQ.3) THEN
        DO K=0,NZM
        DO I=0,NXM
          U2(I,K,NY)=0.d0
          U2(I,K,NY+1)=0.0d0
        END DO
        END DO
       END IF

       
       RETURN
       END


       SUBROUTINE STOCHASTIC_FORCING_CHAN(JF1,JF2)
       
       include 'header'
       real x2,x3,w,y1,y2,y3
       integer i,j,k
       real*8 x1
       real*8 F1_mean(0:NY+1),F2_mean(0:NY+1),F3_mean(0:NY+1)
       real*8 F1_rms(0:NY+1),F2_rms(0:NY+1),F3_rms(0:NY+1)
       real*8 F1_u1(0:NY+1),F2_u2(0:NY+1),F3_u3(0:NY+1)
       real*8 error(0:NY+1),error_int ,U1max(0:NY+1)
       
       integer jf1,jf2
         

       IF ((MOD(TIME_STEP,10).EQ.0).AND.(RK_STEP.EQ.1)) THEN
! First, we need to get the plane averaged velocity
       do j=0,NY+1
       U1_bar(j)=0.d0
       U3_bar(j)=0.d0
       do k=0,NZM
       do i=0,NXM
         U1_bar(j)=U1_bar(j)+U1(i,k,j)
         U3_bar(j)=U3_bar(j)+U3(i,k,j)
       end do
       end do
       U1_bar(j)=U1_bar(j)/dble(NX*NZ)
       U3_bar(j)=U3_bar(j)/dble(NX*NZ)
       end do
! Compute the nondimensional shear and the error from the log law
         do j=1,NY
           dudy(j)=(U1_bar(j)-U1_bar(j-1))
     &          /(GYF(j)-GYF(j-1))
           dwdy(j)=(U3_bar(j)-U3_bar(j-1))
     &          /(GYF(j)-GYF(j-1))
           shear(j)=sqrt(dudy(j)**2.d0+dwdy(j)**2.d0)
           if (gy(j).lt.0.d0) then
             error(j)=shear(j)*(gy(j)-gy(2))/(2.44d0)-1.d0
           else
             error(j)=shear(j)*(gy(NY)-gy(j))/2.44d0-1.d0
           end if
         end do
! Integrate the deviation from the log law

         do j=1,NY
! Adjust pointwise forcing function
            if ((((gy(j)-gy(2)).lt.DELTA_YF_V(5))
     &      .or.((gy(NY)-gy(j)).lt.DELTA_YF_V(5)))
!     &      .and.((j.gt.2).and.(j.lt.NY))) then
     &      .and.((j.gt.2).and.(j.lt.NY))) then
            STOCHASTIC_Y(J)=STOCHASTIC_Y(J)
     &              +error(j)*50.d0
            end if
         end do
       END IF
! Set the vertical forcing function in order to increase the resolved
! Reynolds stress

       do j=0,NY+1
       U1max(j)=0.d0
       do k=0,NZM
       do i=0,NXM
         U1max(j)=max(U1max(j),
     &     ((U1(I,K,J)-U1_bar(j))*DYF(j-1)
     &     +(U1(I,K,J-1)-U1_bar(j-1))*DYF(J))/(2.d0*DY(J)) )
!     &     +((U3(I,K,J)-U3_bar(j))*DYF(j-1)
!     &     +(U3(I,K,J-1)-U3_bar(j-1))*DYF(j))/(2.d0*DY(j)) )
       end do
       end do
       end do

       do j=0,NY+1
       do k=0,NZM
       do i=0,NXM
       call RANDOM_NUMBER(x1)
! Random number between zero and one 
       S1(I,K,J)=STOCHASTIC_Y(J)*x1
!       S1(I,K,J)=STOCHASTIC_Y(J)
! Sign should be dependent on horizontal velocity
       if (gy(j).lt.0.d0) then
! We want a negative correlation with U2 between U1 and U3
       if ( ((U1(I,K,J)-U1_bar(j))*DYF(j-1)
     &           +(U1(I,K,J-1)-U1_bar(j-1))*DYF(J))/(2.d0*DY(J))
     &     +((U3(I,K,J)-U3_bar(j))*DYF(j-1)
     &           +(U3(I,K,J-1)-U3_bar(j-1))*DYF(j))/(2.d0*DY(j))
     &       .gt.0.d0) then
         S1(I,K,J)=-1.d0*S1(I,K,J)
       else
         S1(I,K,J)=1.d0*S1(I,K,J)
       end if
       else
! We want a positive correlation with U2 between U1 and U3
       if ( ((U1(I,K,J)-U1_bar(j))*DYF(j-1)
     &           +(U1(I,K,J-1)-U1_bar(j-1))*DYF(J))/(2.d0*DY(J))
     &     +((U3(I,K,J)-U3_bar(j))*DYF(j-1)
     &           +(U3(I,K,J-1)-U3_bar(j-1))*DYF(j))/(2.d0*DY(j))
     &       .gt.0.d0) then
         S1(I,K,J)=1.d0*S1(I,K,J)
       else
         S1(I,K,J)=-1.d0*S1(I,K,J)
       end if
       end if

       end do
       end do
       end do

! Now, remove the mean of the forcing
       do j=0,NY+1
       F2_mean(j)=0.d0
       do k=0,NZM
       do i=0,NXM
         F2_mean(j)=F2_mean(j)+S1(I,K,J)
       end do
       end do
       F2_mean(j)=F2_mean(j)/dble(NX*NZ)
       end do
       do j=0,NY+1
       do k=0,NZM
       do i=0,NXM
         S1(I,K,J)=S1(I,K,J)-F2_mean(j)
       end do
       end do
       end do

       IF ((MOD(TIME_STEP,SAVE_STATS_INT).EQ.0).AND.(RK_STEP.EQ.1)) THEN
! Save information about the forcing 

       do j=1,NY
         F2_rms(j)=0.d0
         do k=0,NZ-1
         do i=0,NX-1
           F2_rms(j)=F2_rms(j)+S1(I,K,J)**2.d0
         end do
         end do
         F2_rms(j)=sqrt(F2_rms(j)/dble(NX*NZ))
       end do
       do j=1,NY
         F2_u2(j)=0.d0
         do k=0,NZ-1
         do i=0,NX-1
           F2_u2(j)=F2_u2(j)+S1(I,K,J)*U2(I,K,J)
         end do
         end do
         F2_u2(j)=F2_u2(j)/dble(NX*NZ)
       end do

       open(46,file='mean_stochastic.txt',status='unknown'
     &       ,form='formatted')
       write(46,*) TIME_STEP,TIME,DELTA_T
       do j=1,NY
         write(46,460) j,GYF(j),F2_rms(j),F2_u2(j),STOCHASTIC_Y(J)
     &             ,error(j)
       end do
460    format(I3,' ',5(F25.9,' '))


       end if
       CALL FFT_XZ_TO_FOURIER(S1,CS1,0,NY+1)
! Add forcing term as explicit Euler
       do j=3,NY-1
         do k=0,TNKZ
           do i=0,NKX
             if ((i.ne.0).and.(k.ne.0)) then
!             call RANDOM_NUMBER(x1)
!             CR2(I,K,J)=CR2(I,K,J)+H_BAR(RK_STEP)*x1*CS1(I,K,J)
             CR2(I,K,J)=CR2(I,K,J)+CS1(I,K,J)*H_BAR(RK_STEP)
             end if
           end do
         end do
       end do


       return
       end

