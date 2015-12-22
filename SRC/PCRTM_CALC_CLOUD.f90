MODULE PCRTM_CALC_CLOUD
  USE PCRTM_TYPE_KIND
  USE PCRTM_CLOUD_DEFINE, ONLY   : PCRTM_CLOUD_TYPE,CLOUD_TAB_INDEX_TYPE,INIT_PCRTM_CLOUD_PARAM
  USE PCRTM_ATM_ABSORPTION_DEFINE, ONLY : PCRTM_ATM_ABS_STRUCT_TYPE
  USE PCRTM_ATMOSPHERE_DEFINE, ONLY : PCRTM_ATMOSPHERE_TYPE
  USE PCRTM_MATH_UTILITY, ONLY : PLANCK, DPLANCKDT, LOGXLOGP1
  USE PCRTM_CLOUD_LUT_IO, ONLY : PCRTM_CLD_TABLE_DEF
  USE PCRTM_CLOUD_LUT_INTERP, ONLY : PCRTM_INTERP_IDX,PCRTM_INTERP_CLOUD_TAB
  USE PCRTM_JACOBIAN, ONLY : PCRTM_CLD_JACOBIAN_TYPE

CONTAINS
  SUBROUTINE CALC_CLD_PROPERTY( CLD,                  &
                                NCLD,                 &
                                ATM,                  &
                                PCRTM_STND,           &
                                ICE_GRID,             &
                                WAT_GRID,             &
                                SATANG,               &
                                JACOB,                &
                                K_CLD)
    
    INTEGER,                        INTENT(IN)     :: NCLD
    TYPE(PCRTM_ATMOSPHERE_TYPE),    INTENT(IN)     :: ATM
    TYPE(PCRTM_ATM_ABS_STRUCT_TYPE),INTENT(IN)     :: PCRTM_STND
    TYPE(PCRTM_CLD_TABLE_DEF),      INTENT(IN)     :: ICE_GRID
    TYPE(PCRTM_CLD_TABLE_DEF),      INTENT(IN)     :: WAT_GRID
    REAL(SINGLE),                   INTENT(IN)     :: SATANG
    LOGICAL,                        INTENT(IN)     :: JACOB
    TYPE(PCRTM_CLOUD_TYPE),         INTENT(OUT)    :: CLD(NCLD)
    TYPE(PCRTM_CLD_JACOBIAN_TYPE),  INTENT(OUT)    :: K_CLD(NCLD)    

    TYPE(CLOUD_TAB_INDEX_TYPE)                     :: CLD_TAB_INDEX(NCLD)  
    INTEGER                                        :: K

    CALL INIT_PCRTM_CLOUD_PARAM( CLD, NCLD, PCRTM_STND%NM, JACOB, K_CLD )

    CALL PCRTM_ASSIGN_CLOUD_LAYER( CLD, NCLD, ATM, PCRTM_STND )

    DO K = 1,NCLD
       IF (CLD(K)%PHASE .EQ. 1) THEN
!!$          PRINT*,'ICE CLOUD LOCATED AT',CLD(K)%PCLD, 'MB' 
!!$          PRINT*,'CLOUD VISUAL OPTICAL DEPTH', CLD(K)%VISTAU
!!$          PRINT*,'CLOUD PARTICLE SIZE', CLD(K)%DE, 'UM'
!!$          PRINT*,'CLOUD TEMPERATURE',CLD(K)%TCLD,'K'
          CALL PCRTM_INTERP_IDX(SATANG,                   &
                                ICE_GRID,                 &
                                CLD(K),                   &
                                CLD_TAB_INDEX(K))
          DO IM = 1, PCRTM_STND%NM
             CALL PCRTM_INTERP_CLOUD_TAB(CLD(K),          &
                                         ICE_GRID,        & 
                                         CLD_TAB_INDEX(K),&
                                         IM,              &
                                         JACOB)
          END DO
       ELSE IF (CLD(K)%PHASE .EQ. 2)  THEN
!!$          PRINT*,'WATER CLOUD LOCATED AT',CLD(K)%PCLD, 'MB' 
!!$          PRINT*,'CLOUD VISUAL OPTICAL DEPTH', CLD(K)%VISTAU
!!$          PRINT*,'CLOUD PARTICLE SIZE', CLD(K)%DE, 'UM'
!!$          PRINT*,'CLOUD TEMPERATURE',CLD(K)%TCLD,'K'
          CALL PCRTM_INTERP_IDX(SATANG,                   &
                                WAT_GRID,                 &
                                CLD(K),                   &
                                CLD_TAB_INDEX(K))
          DO IM = 1, PCRTM_STND%NM
             CALL PCRTM_INTERP_CLOUD_TAB(CLD(K),          &
                                         WAT_GRID,        & 
                                         CLD_TAB_INDEX(K),&
                                         IM,              &
                                         JACOB)
          END DO
       ELSE
          PRINT*,'WRONG CLOUD PHASE'
          STOP
       END IF
    END DO
    
  END SUBROUTINE CALC_CLD_PROPERTY


 SUBROUTINE PCRTM_ASSIGN_CLOUD_LAYER( CLD, NCLD, ATM, PCRTM_STND )
   INTEGER,                         INTENT(IN) :: NCLD
   TYPE(PCRTM_CLOUD_TYPE),       INTENT(INOUT) :: CLD(NCLD)
   TYPE(PCRTM_ATMOSPHERE_TYPE),     INTENT(IN) :: ATM
   TYPE(PCRTM_ATM_ABS_STRUCT_TYPE), INTENT(IN) :: PCRTM_STND

   INTEGER                                     :: K,L
   REAL(SINGLE)                                :: A

    K=0
    DO L= 1, PCRTM_STND%NLAY
       IF ( ATM%CLD_FLAG(L) .GT.0 ) THEN
          K = K+1
          CLD(K)%ILAY   = L
          CLD(K)%PHASE  = ATM%CLD_FLAG(L)
          CLD(K)%PCLD   = ATM%PCLD(L)
          CLD(K)%VISTAU = ATM%TAUCLD(L)
          CLD(K)%DE     = ATM%DECLD(L)
          CALL LOGXLOGP1(ATM%TLEV(L),ATM%TLEV(L+1),PCRTM_STND%PBND(L),&
               PCRTM_STND%PBND(L+1),CLD(K)%PCLD,CLD(K)%TCLD,CLD(K)%DTCLDDPCLD,A)
!!$          PRINT*, 'DTCLDDPCLD', CLD(K)%DTCLDDPCLD
          CLD(K)%BCLD       = PLANCK(PCRTM_STND%FRQ,CLD(K)%TCLD,PCRTM_STND%NM)
          CLD(K)%DBCLDDTCLD = DPLANCKDT(PCRTM_STND%FRQ,CLD(K)%TCLD,PCRTM_STND%NM)
          IF (ABS(PCRTM_STND%PBND(L+1)-CLD(K)%PCLD).GT.ABS(PCRTM_STND%PBND(L)-CLD(K)%PCLD)) THEN
             CLD(K)%ILEV = L
          ELSE
             CLD(K)%ILEV = L + 1 
          ENDIF
          CLD(K)%DTCLDDTLEVU = (1-A)*CLD(K)%TCLD/ATM%TLEV(L)
          CLD(K)%DTCLDDTLEVD = A*CLD(K)%TCLD/ATM%TLEV(L+1)
       ENDIF
       IF (K.EQ.NCLD) EXIT
    ENDDO

    RETURN
  END SUBROUTINE PCRTM_ASSIGN_CLOUD_LAYER



END MODULE PCRTM_CALC_CLOUD
