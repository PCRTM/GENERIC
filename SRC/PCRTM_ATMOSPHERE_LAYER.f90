MODULE PCRTM_ATMOSPHERE_LAYER
  USE PCRTM_TYPE_KIND
  USE PCRTM_ATMOSPHERE_DEFINE, ONLY    : PCRTM_ATMOSPHERE_TYPE         
  USE PCRTM_ATM_ABSORPTION_DEFINE,ONLY : PCRTM_ATM_ABS_STRUCT_TYPE
                      

  USE PCRTM_CONSTANTS,ONLY : WTDRYAIR,  &
                             WTMOL,     &
                             R_AIR,     &
                             R_GAS,     &
                             R_AVG,     &
                             G_SFC
  USE PCRTM_MATH_UTILITY,ONLY    : LOGXLOGP

  TYPE PCRTM_GEOMETRY_TYPE
     REAL(SINGLE) :: PSFC
     REAL(SINGLE) :: POBS
     INTEGER      :: LSFC
     INTEGER      :: LOBS
     INTEGER      :: NTOP
     INTEGER      :: NBOT
     REAL(SINGLE) :: SCALSFC
     REAL(SINGLE) :: SCALOBS
     REAL(SINGLE) :: SECZANG  
     REAL(SINGLE) :: SATANG   
     REAL(SINGLE) :: SOLAR_ZANG
     REAL(SINGLE) :: SAT_AZIMUTH
     REAL(SINGLE) :: SOLAR_AZIMUTH
  END TYPE PCRTM_GEOMETRY_TYPE

CONTAINS
  SUBROUTINE PCRTM_GEOMETRY_INFO(GEOMETRY,P,NLEV)
    TYPE(PCRTM_GEOMETRY_TYPE),INTENT(INOUT):: GEOMETRY
    INTEGER,      INTENT(IN)               :: NLEV
    REAL(SINGLE), INTENT(IN)               :: P(NLEV)

    DO L=1,NLEV
       IF( ABS(P(L)-GEOMETRY%POBS).LT.0.0001) THEN
          GEOMETRY%LOBS=L
          EXIT
       ELSEIF (P(L)-GEOMETRY%POBS.GT.0.0) THEN
          GEOMETRY%LOBS=L-1
          EXIT
       ENDIF
    ENDDO
    IF (GEOMETRY%LOBS.EQ.0) GEOMETRY%LOBS=1   

    GEOMETRY%SCALOBS = (P(GEOMETRY%LOBS+1)-GEOMETRY%POBS)/(P(GEOMETRY%LOBS+1)-P(GEOMETRY%LOBS))
    GEOMETRY%NTOP    = GEOMETRY%LOBS
    
    DO L=NLEV,1,-1
       IF( ABS(P(L)- GEOMETRY%PSFC).LT.0.0001) THEN
          GEOMETRY%LSFC=L
          EXIT
       ELSEIF (P(L)- GEOMETRY%PSFC.LT.0.0) THEN
          GEOMETRY%LSFC=L+1 
          EXIT
       ENDIF
    ENDDO
    IF (GEOMETRY%LSFC.GT.NLEV) GEOMETRY%LSFC=NLEV
    GEOMETRY%SCALSFC = (GEOMETRY%PSFC-P(GEOMETRY%LSFC-1))/(P(GEOMETRY%LSFC)-P(GEOMETRY%LSFC-1))
    GEOMETRY%NBOT    = GEOMETRY%LSFC-1

    GEOMETRY%SECZANG = 1/COS(GEOMETRY%SATANG*.01745329)

  END SUBROUTINE PCRTM_GEOMETRY_INFO


  SUBROUTINE LAYERAVG( GEOMETRY, ATM, PCRTM_STND )
    TYPE(PCRTM_ATMOSPHERE_TYPE),INTENT(INOUT) :: ATM
    TYPE(PCRTM_ATM_ABS_STRUCT_TYPE),INTENT(IN):: PCRTM_STND
    TYPE(PCRTM_GEOMETRY_TYPE), INTENT(INOUT)  :: GEOMETRY
 
    REAL(SINGLE) :: P(PCRTM_STND%NLEV),T(PCRTM_STND%NLEV)
    REAL(SINGLE) :: WH2O(PCRTM_STND%NLEV)
    REAL(SINGLE) :: TOBS
    REAL(SINGLE) :: SCALOBS
    REAL(SINGLE) :: TSFC
    REAL(SINGLE) :: SCALSFC
    REAL(SINGLE) :: VMR(PCRTM_STND%NLEV,PCRTM_STND%NMOL)
    REAL(SINGLE) :: WOBS
    REAL(SINGLE) :: PSFC
    REAL(SINGLE) :: POBS
    INTEGER      :: L,K
    INTEGER      :: LOBS
    INTEGER      :: LSFC

      DO L = 1, PCRTM_STND%NLEV
         P(L)    = PCRTM_STND%PBND(L)
         T(L)    = ATM%TLEV(L)
         VMR(L,:)= ATM%VMR(L,:) 
      END DO

      LSFC    = GEOMETRY%LSFC
      LOBS    = GEOMETRY%LOBS
      SCALSFC = GEOMETRY%SCALSFC
      SCALOBS = GEOMETRY%SCALOBS
      PSFC    = GEOMETRY%PSFC
      POBS    = GEOMETRY%POBS

!------------------------------------------------------------------------------
!     CALCULATE LAYER AMOUNTS AND LAYER QUANTITIES:
      
      DO L=1,PCRTM_STND%NLEV
         WH2O(L)=VMR(L,1)   !SAVE H2O VMR FOR GHITE ROUTINE (G/KG)
!     ---------------------------------
!     CONVERT H2O FROM G/KG TO VOLUMING MIXING RATIO
         VMR(L,1)=VMR(L,1)/1000.*WTDRYAIR/WTMOL(1) &
             /(1+WTDRYAIR/WTMOL(1)*WH2O(L)/1000.)
!     CONVERT O3 AND CO FROM PPMV TO MIXING RATIO
!!$         VMR(L,3)=VMR(L,3)*1.E-6  !*28.964/47.998
!!$         VMR(L,5)=VMR(L,5)*1.E-6  !*28.964/47.998
!!$         VMR(L,2)=VMR(L,2)*1.E-6
!!$         VMR(L,4)=VMR(L,4)*1.E-6
!!$         VMR(L,6)=VMR(L,6)*1.E-6
         VMR(L,2:PCRTM_STND%NMOL)=VMR(L,2:PCRTM_STND%NMOL)*1.E-6
      ENDDO
      
!---  GET BOUNDARY INDEX AND INTERPOLATE T AND VMR TO BOUNDARIES:

!      CALL GETLOBS(P,POBS,NLEV,LOBS,T,TOBS,SCALOBS)

!--INTERPOLATE T AND MIXING RATIOS:
      CALL LOGXLOGP(T(LOBS),T(LOBS+1),P(LOBS),P(LOBS+1),POBS,TOBS)
      T(LOBS)        = TOBS
      !ATM%TLEV(LOBS) = TOBS
      DO K=1,PCRTM_STND%NMOL
         IF(PCRTM_STND%MOLINDX(K).EQ.2) THEN
            CALL LOGXLOGP(VMR(LOBS,K),VMR(LOBS+1,K),P(LOBS), &
                P(LOBS+1),POBS,WOBS)
            VMR(LOBS,K)     =WOBS
         ENDIF
      ENDDO
      P(LOBS) = POBS

!--INTERPOLATE T AND MIXING RATIOS:
      A = alog(PSFC/P(LSFC))/alog(P(LSFC-1)/P(LSFC))

      TSFC = T(LSFC)*(T(LSFC-1)/T(LSFC))**A
!!$      CALL LOGXLOGP(T(LSFC),T(LSFC-1),P(LSFC),P(LSFC-1),PSFC,TSFC)
!!$      ATM%DLAYDLEV(LSFC) = ALOG(PSFC/P(LSFC-1))/ALOG(P(LSFC)/P(LSFC-1))*TSFC/T(LSFC)
!!$      if ( ATM%DLAYDLEV(LSFC) .lt. 0) then
!!$         print*, 'LSFC, ATM%DLAYDLEV(LSFC),P((LSFC-1):LSFC),PSFC'
!!$         write(*,*)LSFC, ATM%DLAYDLEV(LSFC),P((LSFC-1):LSFC),PSFC
!!$      end if
      ATM%TSFC        = TSFC
      ATM%dTSFCdTLSFC =(1-A)*(T(LSFC-1)/T(LSFC))**A
      ATM%dTSFCdTLSFC_=A*(T(LSFC-1)/T(LSFC))**(A-1) 
      T(LSFC)         = TSFC

!!$      do l = 1,LSFC
!!$         print*, l,ATM%TLEV(l),T(l)
!!$      end do
!!$      print*,ATM%TSFC , ATM%dTSFCdTLSFC, ATM%dTSFCdTLSFC_
      

      DO K=1,PCRTM_STND%NMOL
         IF(PCRTM_STND%MOLINDX(K).EQ.2) THEN
!!$            CALL LOGXLOGP(VMR(LSFC,K)+1e-20,VMR(LSFC-1,K),P(LSFC), &
!!$                P(LSFC-1),PSFC,WSFC)
            WSFC = VMR(LSFC,K)*( VMR(LSFC-1,K)/(VMR(LSFC,K)+1e-20) )**A
            VMR(LSFC,K) = WSFC
            ATM%VSFC(K) = ATM%VMR(LSFC,K)*( ATM%VMR(LSFC-1,K)/(ATM%VMR(LSFC,K)+1e-20) )**A
            ATM%dVSFCdVLSFC(K) =(1-A)*(VMR(LSFC-1,K)/(VMR(LSFC,K)+1E-20))**A
            ATM%dVSFCdVLSFC_(K)=A*(VMR(LSFC-1,K)/(VMR(LSFC,K)+1E-20))**(A-1)
         ENDIF
      ENDDO
      P(LSFC) = PSFC


!---  PERFORM LAYER INTEGRATION TO GET LAYER AMOUNTS AND
!     LAYER AVERAGE QUANTITIES:

      CALL LAYAMT(P,T,WH2O,LSFC,LOBS,ATM,PCRTM_STND%NLEV,VMR,PCRTM_STND%NMOL)
      
!     ----------------------------------------
!     CALCULATE LAYER AMOUNTS FOR ATMOSPHERIC GASES:
!     ----------------------------------------

      L = PCRTM_STND%NT                     ! SELECT STANDARD PROFILE
      DO K=1,PCRTM_STND%NMOL           
         IF (PCRTM_STND%MOLINDX(K).EQ.0) THEN
            ATM%GASPROF(:,K)     = PCRTM_STND%GASSTD(:,L,K)
            ATM%GASPROF(LSFC-1,K)= ATM%GASPROF(LSFC,K)*SCALSFC
            ATM%GASPROF(LOBS,K)  = ATM%GASPROF(LOBS,K)*SCALOBS
         ELSE IF (PCRTM_STND%MOLINDX(K).EQ.1) THEN
            ATM%GASPROF(:,K)     = PCRTM_STND%GASSTD(:,L,K)*PCRTM_STND%SCALFAC(K)
            ATM%GASPROF(LSFC-1,K)= ATM%GASPROF(LSFC,K)*SCALSFC
            ATM%GASPROF(LOBS,K)  = ATM%GASPROF(LOBS,K)*SCALOBS
         ELSE IF (PCRTM_STND%MOLINDX(K).NE.2) THEN
            PRINT*,'WRONG MOLECULE INDEX, AVAILABLE OPTIONS: 0 1 2',PCRTM_STND%MOLINDX(K)
            STOP
         ENDIF
      END DO

      RETURN
    END SUBROUTINE LAYERAVG


! --------------------------------------------------------------------------
    SUBROUTINE LAYAMT(P,T,WH2O,LSFC,LOBS,ATM,NLEV,VMR,NMOL)
        !    ------------
      IMPLICIT NONE
      TYPE(PCRTM_ATMOSPHERE_TYPE),INTENT(INOUT) :: ATM 
      INTEGER      :: LSFC
      INTEGER      :: LOBS
      REAL(SINGLE) :: P(NLEV),T(NLEV),Z(NLEV),WH2O(NLEV)
      REAL(SINGLE) :: WTAIR(NLEV)
      REAL(SINGLE) :: G(NLEV),C(NLEV),RHO_AIR(NLEV)
      INTEGER      :: NLEV, NMOL
      INTEGER      :: L
      INTEGER      :: K
      REAL(SINGLE) :: RHO1
      REAL(SINGLE) :: RHO2
      REAL(SINGLE) :: P1
      REAL(SINGLE) :: P2
      REAL(SINGLE) :: W1
      REAL(SINGLE) :: W2
      REAL(SINGLE) :: DP
      REAL(SINGLE) :: G_AVG
      REAL(SINGLE) :: W_AVG
      REAL(SINGLE) :: R_HGT
      REAL(SINGLE) :: VMR(NLEV,NMOL)
      

!-----------------------------------------------------------------------
!     REFERENCE: AIRS LAYER PACKAGE SCIENCE NOTE, S. HANNON AND L. STROW
!     UMBC (WITH MODIFICATIONS)

!     -- CALCULATE INITIAL VALUES OF PRESSURE HEIGHTS --
      CALL GPHITE(P,T,WH2O,-690.0,NLEV,1,Z)
      
      DO L=LSFC,LOBS,-1

!     CALCULATE MOLECULAR WEIGHT OF MOIST AIR (G/MOL)

         WTAIR(L)=((1.0-WH2O(L)/1000.)*WTDRYAIR)+(WH2O(L)/1000.*WTMOL(1))

!     AIR MASS DENSITY
         C(L)=0.001*WTAIR(L)/R_AIR 

         RHO_AIR(L)=C(L)*P(L)/T(L)
         
!     GRAVITY
         R_HGT=R_AVG+Z(L)   ! (M)
         G(L)=G_SFC*(R_AVG*R_AVG)/(R_HGT*R_HGT) ! (M/S^2)

      END DO
!-------------------------------------------
!    CALCULATE LAYER AIR DENSITY FIRST: 
!----------------------------------------------
      DO L=LSFC-1,LOBS,-1
         RHO1=RHO_AIR(L)
         RHO2=RHO_AIR(L+1)
         ATM%DLAYDLEV(L)=RHO1/(RHO1+RHO2)
         ATM%TLAY(L)=((RHO1*T(L))+(RHO2*T(L+1)))/(RHO1+RHO2)
         P1=P(L)
         P2=P(L+1)
         DP=P2-P1
         G_AVG = 0.5*(G(L)+G(L+1))
         ATM%AIRAMT(L)=DP*10./G_AVG*2./(WTAIR(L)+WTAIR(L+1))
         DO K=1,NMOL         
            W1=AMAX1( VMR(L,K),1.0E-20)
            W2=AMAX1( VMR(L+1,K),1.0E-20)
            W_AVG=((RHO1*W1)+(RHO2*W2))/(RHO1+RHO2)
            ATM%GASPROF(L,K)=W_AVG*ATM%AIRAMT(L)
         END DO
      ENDDO                     !END OF LAYER LOOP

      RETURN
    END SUBROUTINE LAYAMT

    SUBROUTINE GPHITE( P, T, W, Z_SFC, N_LEVELS, I_DIR, Z)
        
      IMPLICIT NONE
!
! PURPOSE:
!  ROUTINE TO COMPUTE GEOPOTENTIAL HEIGHT GIVEN THE ATMOSPHERIC STATE.
!    INCLUDES VIRTUAL TEMPERATURE ADJUSTMENT.
! CREATED:
!  19-SEP-1996 RECEIVED FROM HAL WOOLF
!  ARGUMENTS:
!     INPUT
!    --------
!       P     - REAL*4 PRESSURE ARRAY (MB)
!       T     - REAL*4 TEMPERATURE PROFILE ARRAY (K)
!       W     - REAL*4 WATER VAPOUR PROFILE ARRAY (G/KG)
!     Z_SFC   - REAL*4 SURFACE HEIGHT (M).  0.0 IF NOT KNOWN.
!    N_LEVELS - INT*4 NUMBER OF ELEMENTS USED IN PASSED ARRAYS
!     I_DIR   - INT*4 DIRECTION OF INCREASING LAYER NUMBER
!                 I_DIR = +1, LEVEL(1) == P(TOP)         } SATELLITE/AC
!                             LEVEL(N_LEVELS) == P(SFC)  }    CASE
!                 I_DIR = -1, LEVEL(1) == P(SFC)         } GROUND-BASED
!                             LEVEL(N_LEVELS) == P(TOP)  }    CASE
!     OUTPUT
!    --------
!       Z     - REAL*4 PRESSURE LEVEL HEIGHT ARRAY (M)
! COMMENTS:
!   DIMENSION OF HEIGHT ARRAY MAY NOT NOT BE THE SAME AS THAT OF THE
!     INPUT PROFILE DATA.
!
!=======================================================================
!-----------------------------------------------------------------------
!              -- PREVENT IMPLICIT TYPING OF VARIABLES --
!-----------------------------------------------------------------------
!      IMPLICIT NONE
!-----------------------------------------------------------------------
!                           -- ARGUMENTS --
!-----------------------------------------------------------------------
! -- ARRAYS
      REAL(SINGLE)    P(*), T(*), W(*), Z(*)
! -- SCALARS
      INTEGER N_LEVELS, I_DIR
      REAL(SINGLE)    Z_SFC
!-----------------------------------------------------------------------
!                         -- LOCAL VARIABLES --
!-----------------------------------------------------------------------
! -- PARAMETERS
! -- ROG=R/G/WTDRYAIR=8.3143/9.80665/28.97=0.029265538(KM/K)=29.2655(M/K)
      REAL(SINGLE), PARAMETER :: ROG = 29.2898,  FAC = 0.5 * ROG 
! -- SCALARS
      INTEGER I_START, I_END, L
      REAL(SINGLE)    V_LOWER, V_UPPER, ALGP_LOWER, ALGP_UPPER, HGT
!***********************************************************************
!                         ** EXECUTABLE CODE **
!***********************************************************************
!-----------------------------------------------------------------------
!  -- CALCULATE VIRTUAL TEMPERATURE ADJUSTMENT AND EXPONENTIAL       --
!  -- PRESSURE HEIGHT FOR LEVEL ABOVE SURFACE.  ALSO SET INTEGRATION --
!  -- LOOP BOUNDS                                                    --
!-----------------------------------------------------------------------
      IF( I_DIR .GT. 0 ) THEN   ! DATA STORED TOP DOWN
        V_LOWER = T(N_LEVELS) * ( 1.0 + ( 0.00061 * W(N_LEVELS) ) )
        ALGP_LOWER = ALOG( P(N_LEVELS) )
        I_START = N_LEVELS-1
        I_END   = 1
      ELSE                      ! DATA STORED BOTTOM UP
        V_LOWER = T(1) * ( 1.0 + ( 0.00061 * W(1) ) )
        ALGP_LOWER = ALOG( P(1) )
        I_START = 2
        I_END   = N_LEVELS
      END IF
      HGT = Z_SFC               ! ASSIGN SURFACE HEIGHT
!-----------------------------------------------------------------------
!             -- LOOP OVER LAYERS ALWAYS FROM SFC -> TOP --
!-----------------------------------------------------------------------
      DO L = I_START, I_END, -1*I_DIR
!       ----------------------------------------------------
!       APPLY VIRTUAL TEMPERATURE ADJUSTMENT FOR UPPER LEVEL
!       ----------------------------------------------------
        V_UPPER = T(L)
        IF( P(L) .GE. 300.0 ) &
         V_UPPER = V_UPPER * ( 1.0 + ( 0.00061 * W(L) ) )
!       ----------------------------------------------------- 
!       CALCULATE EXPONENTIAL PRESSURE HEIGHT FOR UPPER LAYER
!       ----------------------------------------------------- 
        ALGP_UPPER = ALOG( P(L) )
!       ----------------
!       CALCULATE HEIGHT USE EQ. D[LN(P)]=-[(WTDRYAIR*G_SFC/R_GAS)/T]*D(HGT)
!       ----------------
        HGT = HGT + ( FAC*(V_UPPER+V_LOWER)*(ALGP_LOWER-ALGP_UPPER) )
!       -------------------------------
!       OVERWRITE VALUES FOR NEXT LAYER
!       -------------------------------
        V_LOWER = V_UPPER
        ALGP_LOWER = ALGP_UPPER
!       ---------------------------------------------
!       STORE HEIGHTS IN SAME DIRECTION AS OTHER DATA
!       ---------------------------------------------
        Z(L) = HGT
      END DO

      RETURN
    END SUBROUTINE GPHITE

  END MODULE PCRTM_ATMOSPHERE_LAYER
