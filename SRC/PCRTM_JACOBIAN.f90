MODULE PCRTM_JACOBIAN
  USE PCRTM_TYPE_KIND
  USE PCRTM_ATM_ABSORPTION_DEFINE, ONLY : PCRTM_ATM_ABS_STRUCT_TYPE

  TYPE PCRTM_NM_JACOBIAN_TYPE
     LOGICAL                   :: JACOB
     REAL(SINGLE), ALLOCATABLE :: DTAUDGAS(:,:,:)
     REAL(SINGLE), ALLOCATABLE :: DTAUDH2O(:,:)
     REAL(SINGLE), ALLOCATABLE :: DTAUDT(:,:)
     REAL(SINGLE), ALLOCATABLE :: DBDT(:,:)
     REAL(SINGLE), ALLOCATABLE :: DRUPDT(:,:) 
     REAL(SINGLE), ALLOCATABLE :: DRUPDGAS(:,:,:)
     REAL(SINGLE), ALLOCATABLE :: DRUPDH2O(:,:) 
     REAL(SINGLE), ALLOCATABLE :: DRDNCLDDT(:,:)
     REAL(SINGLE), ALLOCATABLE :: DRDNCLDDGAS(:,:,:)
     REAL(SINGLE), ALLOCATABLE :: DRDNCLDDH2O(:,:)
     REAL(SINGLE), ALLOCATABLE :: DRDNDT(:,:) 
     REAL(SINGLE), ALLOCATABLE :: DRDNDGAS(:,:,:)
     REAL(SINGLE), ALLOCATABLE :: DRDNDH2O(:,:)
     REAL(SINGLE), ALLOCATABLE :: DRUPDEM(:)
     REAL(SINGLE), ALLOCATABLE :: DRUPDTS(:)
  END TYPE PCRTM_NM_JACOBIAN_TYPE

  TYPE PCRTM_CLD_JACOBIAN_TYPE
     REAL(SINGLE), ALLOCATABLE :: DRUPDCLDRF(:)
     REAL(SINGLE), ALLOCATABLE :: DRUPDCLDTR(:)
     REAL(SINGLE), ALLOCATABLE :: DRUPDTCLD(:)
     REAL(SINGLE), ALLOCATABLE :: DRDNDCLDRF(:)
     REAL(SINGLE), ALLOCATABLE :: DRDNDCLDTR(:)
     REAL(SINGLE), ALLOCATABLE :: DRDNDTCLD(:)
     REAL(SINGLE), ALLOCATABLE :: DRDNCLDDCLDRF(:)
     REAL(SINGLE), ALLOCATABLE :: DRDNCLDDCLDTR(:)
     REAL(SINGLE), ALLOCATABLE :: DRDNCLDDTCLD(:)
     REAL(SINGLE), ALLOCATABLE :: DRDCLDRF(:)
     REAL(SINGLE), ALLOCATABLE :: DRDCLDTR(:)
     REAL(SINGLE), ALLOCATABLE :: DRDTCLD(:)
     REAL(SINGLE), ALLOCATABLE :: DRDPCLD(:)
     REAL(SINGLE), ALLOCATABLE :: DRDDECLD(:)
     REAL(SINGLE), ALLOCATABLE :: DRDTAUCLD(:)
  END TYPE PCRTM_CLD_JACOBIAN_TYPE

!!$  TYPE CH_JACOBIAN_TYPE
!!$     LOGICAL                   :: JACOB_CH
!!$     REAL(SINGLE), ALLOCATABLE :: R_TLEV(:,:), R_H2OLEV(:,:)
!!$     REAL(SINGLE), ALLOCATABLE :: R_N2OLEV(:,:), R_CH4LEV(:,:), R_CO2LEV(:,:)
!!$     REAL(SINGLE), ALLOCATABLE :: R_EM(:), R_TS(:)
!!$     REAL(SINGLE), ALLOCATABLE :: R_COLEV(:,:), R_O3LEV(:,:)
!!$  END TYPE CH_JACOBIAN_TYPE

  TYPE PCRTM_CH_JACOBIAN_TYPE
     LOGICAL                   :: JACOB_CH
     REAL(SINGLE), ALLOCATABLE :: R_TLEV(:,:)
     REAL(SINGLE), ALLOCATABLE :: R_GASLEV(:,:,:)
     REAL(SINGLE), ALLOCATABLE :: R_EM(:), R_TS(:)
     LOGICAL                   :: JACOB_BTCH
     REAL(SINGLE), ALLOCATABLE :: BT_TLEV(:,:)
     REAL(SINGLE), ALLOCATABLE :: BT_GASLEV(:,:,:)
     REAL(SINGLE), ALLOCATABLE :: BT_EM(:), BT_TS(:)
  END TYPE PCRTM_CH_JACOBIAN_TYPE

  TYPE PCRTM_PC_JACOBIAN_TYPE
     LOGICAL                   :: JACOB_PC
     REAL(SINGLE), ALLOCATABLE :: K_T(:,:),K_TLAY(:,:),K_EM(:), K_TS(:)
!!$     REAL(SINGLE), ALLOCATABLE :: K_H2O(:,:)
!!$     REAL(SINGLE), ALLOCATABLE :: K_CO2(:,:), K_CH4(:,:), K_N2O(:,:)
!!$     REAL(SINGLE), ALLOCATABLE :: K_CO(:,:), K_O3(:,:) 
!!$     REAL(SINGLE), ALLOCATABLE :: K_H2OLAY(:,:)
!!$     REAL(SINGLE), ALLOCATABLE :: K_CO2LAY(:,:), K_N2OLAY(:,:), K_CH4LAY(:,:)
!!$     REAL(SINGLE), ALLOCATABLE :: K_COLAY(:,:), K_O3LAY(:,:) 
     REAL(SINGLE), ALLOCATABLE :: K_gaslev(:,:,:),K_gaslay(:,:,:)
  END TYPE PCRTM_PC_JACOBIAN_TYPE

CONTAINS

  SUBROUTINE INIT_PCRTM_JACOB_NM( K_NM, PCRTM_STND )
    TYPE(PCRTM_NM_JACOBIAN_TYPE),   INTENT(OUT) :: K_NM
    TYPE(PCRTM_ATM_ABS_STRUCT_TYPE),INTENT(IN)  :: PCRTM_STND

    INTEGER :: NLAY,NM,NMOL
    INTEGER :: ALLOC_STAT

    K_NM%JACOB = .TRUE.

    NM    = PCRTM_STND%NM
    NMOL  = PCRTM_STND%NMOL
    NLAY  = PCRTM_STND%NLAY
    
    ALLOCATE(K_NM%DTAUDGAS(NMOL+1, NM, NLAY),  & 
             K_NM%DTAUDH2O(NM, NLAY),          &
             K_NM%DTAUDT(NM, NLAY),            &
             K_NM%DBDT(NM, NLAY),              &
             K_NM%DRUPDT(NM, NLAY),            &
             K_NM%DRUPDGAS(NMOL+1,NM,NLAY),    &
             K_NM%DRUPDH2O(NM, NLAY),          &
             K_NM%DRDNCLDDT(NM, NLAY),         &
             K_NM%DRDNCLDDGAS(NMOL+1,NM,NLAY), &
             K_NM%DRDNCLDDH2O(NM, NLAY),       &
             K_NM%DRDNDT(NM, NLAY),            &
             K_NM%DRDNDGAS(NMOL+1,NM,NLAY),    &
             K_NM%DRDNDH2O(NM, NLAY),          &
             K_NM%DRUPDTS(NM),                 &
             K_NM%DRUPDEM(NM),                 &
             STAT = ALLOC_STAT)
    
    IF ( ALLOC_STAT /= 0 ) THEN
       PRINT*,'ERROR TRYING TO ALLOCATE K_NM'
       STOP
    ENDIF
       
    K_NM%DTAUDGAS       = 0
    K_NM%DTAUDH2O       = 0
    K_NM%DTAUDT         = 0
    K_NM%DBDT           = 0
    K_NM%DRUPDT         = 0
    K_NM%DRUPDGAS       = 0
    K_NM%DRUPDH2O       = 0
    K_NM%DRDNCLDDT      = 0
    K_NM%DRDNCLDDGAS    = 0
    K_NM%DRDNCLDDH2O    = 0
    K_NM%DRDNDGAS       = 0.0
    K_NM%DRDNDH2O       = 0     
    K_NM%DRDNDT         = 0    
    K_NM%DRUPDTS        = 0
    K_NM%DRUPDEM        = 0
    
  END SUBROUTINE INIT_PCRTM_JACOB_NM

  SUBROUTINE CLEAR_PCRTM_JACOB_NM( K_NM )
    TYPE(PCRTM_NM_JACOBIAN_TYPE), INTENT(INOUT) :: K_NM
    INTEGER                                     :: DEALLOC_STAT
    
    IF(K_NM%JACOB) THEN
    
       DEALLOCATE(K_NM%DTAUDGAS,                & 
                  K_NM%DTAUDH2O,                &
                  K_NM%DTAUDT,                  &
                  K_NM%DBDT,                    &
                  K_NM%DRUPDT,                  &
                  K_NM%DRUPDGAS,                &
                  K_NM%DRUPDH2O,                &
                  K_NM%DRDNDT,                  &
                  K_NM%DRDNDGAS,                &
                  K_NM%DRDNDH2O,                &
                  K_NM%DRDNCLDDT,               &
                  K_NM%DRDNCLDDGAS,             &
                  K_NM%DRDNCLDDH2O,             &
                  K_NM%DRUPDTS,                 &
                  K_NM%DRUPDEM,                 &
                  STAT = DEALLOC_STAT)
       
       IF ( DEALLOC_STAT /= 0 ) THEN
          PRINT*,'ERROR TRYING TO DEALLOCATE K_NM'
          STOP
       ENDIF

    END IF

  END SUBROUTINE CLEAR_PCRTM_JACOB_NM

  
  SUBROUTINE INIT_PCRTM_CLOUD_JACOB( K_CLD, NCLD, NM)
    TYPE(PCRTM_CLD_JACOBIAN_TYPE),INTENT(OUT)    :: K_CLD(NCLD)
    INTEGER,                      INTENT(IN)     :: NCLD
    INTEGER,                      INTENT(IN)     :: NM

    INTEGER                                      :: K
    INTEGER                                      :: ALLOC_STAT

    DO K = 1, NCLD
       ALLOCATE( K_CLD(K)%DRUPDCLDRF(NM),     &
                 K_CLD(K)%DRUPDCLDTR(NM),     &
                 K_CLD(K)%DRUPDTCLD(NM),      &
                 K_CLD(K)%DRDNDCLDRF(NM),     &
                 K_CLD(K)%DRDNDCLDTR(NM),     &
                 K_CLD(K)%DRDNDTCLD(NM),      &
                 K_CLD(K)%DRDNCLDDCLDRF(NM),  &
                 K_CLD(K)%DRDNCLDDCLDTR(NM),  &
                 K_CLD(K)%DRDNCLDDTCLD(NM),   &
                 K_CLD(K)%DRDCLDRF(NM),       &
                 K_CLD(K)%DRDCLDTR(NM),       &
                 K_CLD(K)%DRDTCLD(NM),        &
                 K_CLD(K)%DRDPCLD(NM),        &
                 K_CLD(K)%DRDDECLD(NM),       &
                 K_CLD(K)%DRDTAUCLD(NM),      &
                 STAT = ALLOC_STAT )
       IF ( ALLOC_STAT /= 0 ) THEN
          PRINT*,'ERROR TRYING TO ALLOCATE CLOUD JACOBIAN PARAM'
          STOP
       ENDIF

       K_CLD(K)%DRUPDCLDRF   = 0
       K_CLD(K)%DRUPDCLDTR   = 0
       K_CLD(K)%DRUPDTCLD    = 0
       K_CLD(K)%DRDNDCLDRF   = 0
       K_CLD(K)%DRDNDCLDTR   = 0
       K_CLD(K)%DRDNDTCLD    = 0
       K_CLD(K)%DRDNDCLDRF   = 0
       K_CLD(K)%DRDNDCLDTR   = 0
       K_CLD(K)%DRDNCLDDTCLD = 0
       K_CLD(K)%DRDCLDRF     = 0 
       K_CLD(K)%DRDCLDTR     = 0
       K_CLD(K)%DRDTCLD      = 0
       K_CLD(K)%DRDPCLD      = 0   
       K_CLD(K)%DRDDECLD     = 0
       K_CLD(K)%DRDTAUCLD    = 0
    END DO
  END SUBROUTINE INIT_PCRTM_CLOUD_JACOB

  SUBROUTINE CLEAR_PCRTM_CLOUD_JACOB( K_CLD, NCLD)
  
    TYPE(PCRTM_CLD_JACOBIAN_TYPE),INTENT(INOUT)  :: K_CLD(NCLD)
    INTEGER,                      INTENT(IN)     :: NCLD

    INTEGER                                     :: K
    INTEGER                                     :: DEALLOC_STAT

    DO K = 1, NCLD
       DEALLOCATE( K_CLD(K)%DRUPDCLDRF,                  &
                   K_CLD(K)%DRUPDCLDTR,                  &
                   K_CLD(K)%DRUPDTCLD,                   &
                   K_CLD(K)%DRDNDCLDRF,                  &
                   K_CLD(K)%DRDNDCLDTR,                  &
                   K_CLD(K)%DRDNDTCLD,                   &
                   K_CLD(K)%DRDNCLDDCLDRF,               &
                   K_CLD(K)%DRDNCLDDCLDTR,               &
                   K_CLD(K)%DRDNCLDDTCLD,                &
                   K_CLD(K)%DRDCLDRF,                    &
                   K_CLD(K)%DRDCLDTR,                    &
                   K_CLD(K)%DRDTCLD,                     &
                   K_CLD(K)%DRDPCLD,                     &
                   K_CLD(K)%DRDDECLD,                    &
                   K_CLD(K)%DRDTAUCLD,                   &
                   STAT = DEALLOC_STAT )
       IF ( DEALLOC_STAT /= 0 ) THEN
          PRINT*,'ERROR TRYING TO ALLOCATE CLOUD JACOBIAN PARAM'
          STOP
       ENDIF
    END DO
       
  END SUBROUTINE CLEAR_PCRTM_CLOUD_JACOB

  SUBROUTINE INIT_PCRTM_JACOB_PC( K_PC, NPC, NLAY)
    
    TYPE(PCRTM_PC_JACOBIAN_TYPE), INTENT(INOUT) :: K_PC
    INTEGER,                      INTENT(IN)    :: NPC
    INTEGER,                      INTENT(IN)    :: NLAY
    
    INTEGER                                     :: NLEV
    INTEGER                                     :: ALLOC_STAT
    
    IF(K_PC%JACOB_PC) THEN

       NLEV = NLAY + 1
       
!!$       ALLOCATE( K_PC%K_T(NPC,NLEV),                         &
!!$                 K_PC%K_CO2(NPC,NLEV),                       &
!!$                 K_PC%K_H2O(NPC,NLEV),                       &
!!$                 K_PC%K_CO(NPC,NLEV),                        &
!!$                 K_PC%K_CH4(NPC,NLEV),                       &
!!$                 K_PC%K_O3(NPC,NLEV),                        &
!!$                 K_PC%K_N2O(NPC,NLEV),                       &
!!$                 K_PC%K_TLAY(NPC,NLAY),                      &
!!$                 K_PC%K_CO2LAY(NPC,NLAY),                    &
!!$                 K_PC%K_H2OLAY(NPC,NLAY),                    &
!!$                 K_PC%K_COLAY(NPC,NLAY),                     &
!!$                 K_PC%K_CH4LAY(NPC,NLAY),                    &
!!$                 K_PC%K_O3LAY(NPC,NLAY),                     &
!!$                 K_PC%K_N2OLAY(NPC,NLAY),                    &
!!$                 K_PC%K_EM(NPC),                             & 
!!$                 K_PC%K_TS(NPC),                             &
!!$                 STAT = ALLOC_STAT )

       ALLOCATE( K_PC%K_T(NPC,NLEV),                         &
                 K_PC%K_TLAY(NPC,NLAY),                      &
                 K_PC%K_GASLEV(NPC,NLEV,15),                 &
                 K_PC%K_GASLAY(NPC,NLAY,15),                 &
                 K_PC%K_EM(NPC),                             & 
                 K_PC%K_TS(NPC),                             &
                 STAT = ALLOC_STAT )

       
       IF ( ALLOC_STAT /= 0 ) THEN
          PRINT*,'ERROR TRYING TO ALLOCATE K_PC'
          STOP
       ENDIF

       K_PC%K_T      =     0
       K_PC%K_GASLEV =     0
       K_PC%K_TLAY   =     0  
       K_PC%K_GASLAY =     0
       K_PC%K_EM     =     0
       K_PC%K_TS     =     0

    END IF
    
  END SUBROUTINE INIT_PCRTM_JACOB_PC

  SUBROUTINE CLEAR_PCRTM_JACOB_PC( K_PC )
    
    TYPE(PCRTM_PC_JACOBIAN_TYPE), INTENT(INOUT) :: K_PC

    INTEGER                                     :: DEALLOC_STAT

    IF(K_PC%JACOB_PC) THEN
    
       DEALLOCATE( K_PC%K_T,                         &
                   K_PC%K_GASLEV,                    &
                   K_PC%K_TLAY,                      &
                   K_PC%K_GASLAY,                    &
                   K_PC%K_EM,                        & 
                   K_PC%K_TS,                        &
                   STAT = DEALLOC_STAT )
       
       IF ( DEALLOC_STAT /= 0 ) THEN
          PRINT*,'ERROR TRYING TO DEALLOCATE K_PC'
          STOP
       ENDIF

    END IF

  END SUBROUTINE CLEAR_PCRTM_JACOB_PC
  

!!$  SUBROUTINE INIT_PCRTM_JACOB_CH( K_CH, NCH, NLEV)
!!$    
!!$    TYPE(CH_JACOBIAN_TYPE), INTENT(OUT)   :: K_CH
!!$    INTEGER, INTENT(IN)                   :: NCH
!!$    INTEGER, INTENT(IN)                   :: NLEV
!!$    
!!$    INTEGER                               :: NPC
!!$    INTEGER                               :: ALLOC_STAT
!!$
!!$    ALLOCATE( K_CH%R_TLEV(NCH,NLEV),                         &
!!$              K_CH%R_CO2LEV(NCH,NLEV),                       &
!!$              K_CH%R_H2OLEV(NCH,NLEV),                       &
!!$              K_CH%R_COLEV(NCH,NLEV),                        &
!!$              K_CH%R_CH4LEV(NCH,NLEV),                       &
!!$              K_CH%R_O3LEV(NCH,NLEV),                        &
!!$              K_CH%R_N2OLEV(NCH,NLEV),                       &
!!$              K_CH%R_EM(NCH),                                &
!!$              K_CH%R_TS(NCH),                                &
!!$              STAT = ALLOC_STAT )
!!$
!!$    IF ( ALLOC_STAT /= 0 ) THEN
!!$       PRINT*,'ERROR TRYING TO ALLOCATE K_CH'
!!$       STOP
!!$    ENDIF
!!$
!!$  END SUBROUTINE INIT_PCRTM_JACOB_CH
!!$
!!$  SUBROUTINE CLEAR_PCRTM_JACOB_CH( K_CH )
!!$    
!!$    TYPE(CH_JACOBIAN_TYPE), INTENT(INOUT) :: K_CH
!!$
!!$    INTEGER                               :: DEALLOC_STAT
!!$
!!$    DEALLOCATE( K_CH%R_TLEV,                       &
!!$              K_CH%R_CO2LEV,                       &
!!$              K_CH%R_H2OLEV,                       &
!!$              K_CH%R_COLEV,                        &
!!$              K_CH%R_CH4LEV,                       &
!!$              K_CH%R_O3LEV,                        &
!!$              K_CH%R_N2OLEV,                       &
!!$              K_CH%R_EM,                           &
!!$              K_CH%R_TS,                           &
!!$              STAT = DEALLOC_STAT )
!!$
!!$    IF ( DEALLOC_STAT /= 0 ) THEN
!!$       PRINT*,'ERROR TRYING TO DEALLOCATE K_CH'
!!$       STOP
!!$    ENDIF
!!$
!!$  END SUBROUTINE CLEAR_PCRTM_JACOB_CH

  SUBROUTINE INIT_PCRTM_JACOB_CH( K_CH, NCH, NLEV)
    
    TYPE(PCRTM_CH_JACOBIAN_TYPE), INTENT(INOUT) :: K_CH
    INTEGER,                      INTENT(IN)    :: NCH
    INTEGER,                      INTENT(IN)    :: NLEV

    INTEGER                                     :: ALLOC_STAT

    IF(K_CH%JACOB_CH) THEN
       ALLOCATE( K_CH%R_TLEV(NCH,NLEV),                       &
                 K_CH%R_GASLEV(NCH,NLEV,15),                  &
                 K_CH%R_EM(NCH),                              &
                 K_CH%R_TS(NCH),                              &
                 STAT = ALLOC_STAT )
       IF ( ALLOC_STAT /= 0 ) THEN
          PRINT*,'ERROR TRYING TO ALLOCATE K_CH -- Radiance'
          STOP
       ENDIF
       K_CH%R_TLEV     =   0
       K_CH%R_GASLEV   =   0
       K_CH%R_EM       =   0
       K_CH%R_TS       =   0
    END IF

    IF(K_CH%JACOB_BTCH) THEN
       ALLOCATE( K_CH%BT_TLEV(NCH,NLEV),                      &
                 K_CH%BT_GASLEV(NCH,NLEV,15),                 &
                 K_CH%BT_EM(NCH),                             &
                 K_CH%BT_TS(NCH),                             &
                 STAT = ALLOC_STAT )
       IF ( ALLOC_STAT /= 0 ) THEN
          PRINT*,'ERROR TRYING TO ALLOCATE K_CH -- Brightness Temp'
          STOP
       ENDIF
       K_CH%BT_TLEV     =   0
       K_CH%BT_GASLEV   =   0
       K_CH%BT_EM       =   0
       K_CH%BT_TS       =   0       
    END IF


  END SUBROUTINE INIT_PCRTM_JACOB_CH

  SUBROUTINE CLEAR_PCRTM_JACOB_CH( K_CH )
    
    TYPE(PCRTM_CH_JACOBIAN_TYPE), INTENT(INOUT) :: K_CH

    INTEGER                                     :: DEALLOC_STAT

    IF (K_CH%JACOB_CH) THEN
       DEALLOCATE( K_CH%R_TLEV,                                 &
                   K_CH%R_GASLEV,                               &
                   K_CH%R_EM,                                   &
                   K_CH%R_TS,                                   &
                   STAT = DEALLOC_STAT )
       
       IF ( DEALLOC_STAT /= 0 ) THEN
          PRINT*,'ERROR TRYING TO DEALLOCATE K_CH -- Radiance'
          STOP
       ENDIF
    END IF

    IF(K_CH%JACOB_BTCH) THEN
       DEALLOCATE( K_CH%BT_TLEV,                           &
                   K_CH%BT_GASLEV,                         &
                   K_CH%BT_EM,                             &
                   K_CH%BT_TS,                             &
                   STAT = DEALLOC_STAT )
       IF ( DEALLOC_STAT /= 0 ) THEN
          PRINT*,'ERROR TRYING TO DEALLOCATE K_CH -- Brightness Temp'
          STOP
       ENDIF
    END IF

  END SUBROUTINE CLEAR_PCRTM_JACOB_CH

END MODULE PCRTM_JACOBIAN
