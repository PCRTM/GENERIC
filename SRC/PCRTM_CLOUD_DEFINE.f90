!#############################################################################
! Module cloud_parameter
!
! Defines data type CLOUD_PARAMETER_TYPE
!#############################################################################
MODULE PCRTM_CLOUD_DEFINE

  USE PCRTM_TYPE_KIND
  USE PCRTM_JACOBIAN, ONLY  : PCRTM_CLD_JACOBIAN_TYPE, &
                              INIT_PCRTM_CLOUD_JACOB,  &
                              CLEAR_PCRTM_CLOUD_JACOB

  TYPE PCRTM_CLOUD_TYPE
     REAL(SINGLE)  :: PCLD          ! CLOUD PRES
     REAL(SINGLE)  :: TCLD          ! CLOUD TEMPERATURE
     REAL(SINGLE)  :: VISTAU        ! CLOUD VISTAU
     REAL(SINGLE)  :: DE            ! CLOUD PARTICLE SIZE
     INTEGER       :: ILEV          ! CLOUD LEVEL
     INTEGER       :: ILAY          ! CLOUD LAYER
     INTEGER       :: PHASE         ! CLOUD PHASE
     REAL(SINGLE)  :: DTCLDDPCLD    ! CLOUD LAYER TEMP GRADIENT
     REAL(SINGLE)  :: DTCLDDTLEVU, DTCLDDTLEVD ! CLOUD LEVEL TEMP DERIV
     REAL(SINGLE), &
      ALLOCATABLE  :: BCLD(:)       ! CLOUD PLANCK RADIANCE
     REAL(SINGLE), &
      ALLOCATABLE  :: DBCLDDTCLD(:) ! CLOUD PLANCK DERIVATIVE
     REAL(SINGLE), &
      ALLOCATABLE  :: TRANS(:)      ! CLOUD TRANS(NM)
     REAL(SINGLE), &
      ALLOCATABLE  :: REFL(:)       ! CLOUD REFL(NM)
     REAL(SINGLE), &
      ALLOCATABLE  :: DRFDDE(:)     ! JACOBIAN CLOUD REFLECTANCE/PARTICLE SIZE
     REAL(SINGLE), &
      ALLOCATABLE  :: DRFDTAU(:)    ! JACOBIAN CLOUD REFLECTANCE/OPTICAL DEPTH
     REAL(SINGLE), &
      ALLOCATABLE  :: DTRDDE(:)     ! JACOBIAN CLOUD TRANSMITTANCE/PARTICLE SIZE 
     REAL(SINGLE), &
      ALLOCATABLE  :: DTRDTAU(:)    ! JACOBIAN CLOUD TRANSMITTANCE/OPTICAL DEPTH
  END TYPE PCRTM_CLOUD_TYPE
 
  TYPE CLOUD_TAB_INDEX_TYPE
     INTEGER       :: TAU_IDX
     INTEGER       :: TAU_IDX_
     INTEGER       :: D_IDX
     INTEGER       :: D_IDX_
     INTEGER       :: THET_IDX
     INTEGER       :: THET_IDX_
     REAL(SINGLE)  :: R_TAU
     REAL(SINGLE)  :: R_D
     REAL(SINGLE)  :: R_THET
  END TYPE CLOUD_TAB_INDEX_TYPE

CONTAINS

!##############################################################################
!---  DETERMINE THE CLOUD TEMPERATURE, LEVEL, TEMPERATURE GRADIENT AT THE 
!     CLOUD HEIGHT
!##############################################################################
  
  SUBROUTINE INIT_PCRTM_CLOUD_PARAM( CLD, NCLD, NM, JACOB, K_CLD )
    
    INTEGER,                      INTENT(IN) :: NCLD
    INTEGER,                      INTENT(IN) :: NM
    TYPE(PCRTM_CLOUD_TYPE),       INTENT(OUT):: CLD(NCLD)
    TYPE(PCRTM_CLD_JACOBIAN_TYPE),INTENT(OUT):: K_CLD(NCLD)
    LOGICAL,                      INTENT(IN) :: JACOB

    INTEGER                                  :: K
    INTEGER                                  :: ALLOC_STAT

    DO K = 1, NCLD
       ALLOCATE( CLD(K)%BCLD(NM),            &
                 CLD(K)%DBCLDDTCLD(NM),      &
                 CLD(K)%TRANS(NM),           &
                 CLD(K)%REFL(NM),            &
                 CLD(K)%DRFDDE(NM),          &       
                 CLD(K)%DRFDTAU(NM),         &       
                 CLD(K)%DTRDDE(NM),          &       
                 CLD(K)%DTRDTAU(NM),         &
                 STAT = ALLOC_STAT )
       IF ( ALLOC_STAT /= 0 ) THEN
          PRINT*,'ERROR TRYING TO ALLOCATE CLOUD PARAM'
          STOP
       ENDIF
    END DO

    IF(JACOB) THEN
       CALL INIT_PCRTM_CLOUD_JACOB( K_CLD, NCLD, NM)
    END IF

  END SUBROUTINE INIT_PCRTM_CLOUD_PARAM
    
  SUBROUTINE CLEAR_PCRTM_CLOUD_PARAM( CLD, NCLD, JACOB, K_CLD)
    TYPE(PCRTM_CLOUD_TYPE),       INTENT(INOUT) :: CLD(NCLD)  
    TYPE(PCRTM_CLD_JACOBIAN_TYPE),INTENT(INOUT) :: K_CLD(NCLD)
    LOGICAL,                      INTENT(IN)    :: JACOB

    INTEGER                                     :: K
    INTEGER                                     :: DEALLOC_STAT

    DO  K = 1, NCLD
       DEALLOCATE( CLD(K)%BCLD,                    &
                   CLD(K)%TRANS,                   &
                   CLD(K)%REFL,                    &
                   CLD(K)%DRFDDE,                  &       
                   CLD(K)%DRFDTAU,                 &       
                   CLD(K)%DTRDDE,                  &       
                   CLD(K)%DTRDTAU,                 &
                   STAT = DEALLOC_STAT )
       IF ( DEALLOC_STAT /= 0 ) THEN
          PRINT*,'ERROR TRYING TO DEALLOCATE CLOUD PARAM'
          STOP
       END IF
    END DO

 
    IF(JACOB) THEN
       CALL CLEAR_PCRTM_CLOUD_JACOB( K_CLD, NCLD )
    END IF

  END SUBROUTINE CLEAR_PCRTM_CLOUD_PARAM


END MODULE PCRTM_CLOUD_DEFINE
