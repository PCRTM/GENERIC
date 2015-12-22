! PCRTM_ATMOSPHERE_DEFINE
!
! DEFINING Atmosphere level and layer properties

! CREATION HISTORY:
!       MODIFIED BY:     WAN WU, 27-AGU-2012



MODULE PCRTM_ATMOSPHERE_DEFINE

  USE PCRTM_TYPE_KIND

  TYPE PCRTM_ATMOSPHERE_TYPE
    REAL(SINGLE), ALLOCATABLE  :: TLEV(:)          ! TLEVELS(NLEV)
    REAL(SINGLE), ALLOCATABLE  :: VMR(:,:)         ! VMR(NLEV,NMOL)
    REAL(SINGLE), ALLOCATABLE  :: TLAY(:)          ! TLAY(NLAY)
    REAL(SINGLE), ALLOCATABLE  :: GASPROF(:,:)     ! GASPROF(NLAY,NMOL)
    REAL(SINGLE), ALLOCATABLE  :: AIRAMT(:)        ! AIRAMT(NLAY)
    REAL(SINGLE), ALLOCATABLE  :: DLAYDLEV(:)      ! DLAYDLEV(NLAY)
    INTEGER,      ALLOCATABLE  :: CLD_FLAG(:)      ! CLOUD FLAG FOR EACH LAYER 0 -CLEAR 1 - ICE 2 - WATER
    REAL(SINGLE), ALLOCATABLE  :: PCLD(:)          ! CLOUD PARAM FOR EACH LAYER
    REAL(SINGLE), ALLOCATABLE  :: TAUCLD(:)
    REAL(SINGLE), ALLOCATABLE  :: DECLD(:)
    REAL(SINGLE)               :: TSKIN
    REAL(SINGLE)               :: TSFC
    REAL(SINGLE), ALLOCATABLE  :: VSFC(:)
    REAL(SINGLE)               :: DTSFCDTLSFC,DTSFCDTLSFC_
    REAL(SINGLE), ALLOCATABLE  :: DVSFCDVLSFC(:),DVSFCDVLSFC_(:)
  END TYPE PCRTM_ATMOSPHERE_TYPE

CONTAINS


  SUBROUTINE INIT_PCRTM_ATMOSPHERE( ATMOSPHERE,NLAY,NMOL )
    TYPE(PCRTM_ATMOSPHERE_TYPE), INTENT(OUT)   :: ATMOSPHERE
    INTEGER, INTENT(IN) :: NLAY, NMOL

    INTEGER :: ALLOC_STAT

    ! PERFORM THE ALLOCATION
    ALLOCATE( ATMOSPHERE%TLEV(NLAY+1),                      &
              ATMOSPHERE%VMR(NLAY+1,NMOL),                  &
              ATMOSPHERE%TLAY(NLAY),                        &
              ATMOSPHERE%GASPROF(NLAY,NMOL+1),              &
              ATMOSPHERE%AIRAMT(NLAY),                      &
              ATMOSPHERE%DLAYDLEV(NLAY),                    &
              ATMOSPHERE%CLD_FLAG(NLAY),                    &
              ATMOSPHERE%PCLD(NLAY),                        &
              ATMOSPHERE%TAUCLD(NLAY),                      &
              ATMOSPHERE%DECLD(NLAY),                       &
              ATMOSPHERE%VSFC(NMOL),                        &
              ATMOSPHERE%DVSFCDVLSFC(NMOL),                 &
              ATMOSPHERE%DVSFCDVLSFC_(NMOL),                &
              STAT = ALLOC_STAT )
    
    IF ( ALLOC_STAT /= 0 ) THEN
       PRINT*,'ERROR TRYING TO ALLOCATE ATMOSPHERE'
       STOP
    ENDIF
    
    ATMOSPHERE%TLEV       = 0
    ATMOSPHERE%VMR        = 0
    ATMOSPHERE%TLAY       = 0
    ATMOSPHERE%GASPROF    = 0
    ATMOSPHERE%AIRAMT     = 0
    ATMOSPHERE%DLAYDLEV   = 0
    ATMOSPHERE%CLD_FLAG   = 0
    ATMOSPHERE%PCLD       = 0
    ATMOSPHERE%TAUCLD     = 0
    ATMOSPHERE%DECLD      = 0

  END SUBROUTINE INIT_PCRTM_ATMOSPHERE

  SUBROUTINE CLEAR_PCRTM_ATMOSPHERE( ATMOSPHERE)
    TYPE(PCRTM_ATMOSPHERE_TYPE), INTENT(INOUT) :: ATMOSPHERE
    INTEGER :: DEALLOC_STAT

    ! PERFORM THE ALLOCATION
    DEALLOCATE( ATMOSPHERE%TLEV,                        &
                ATMOSPHERE%VMR,                         &
                ATMOSPHERE%TLAY,                        &
                ATMOSPHERE%GASPROF,                     &
                ATMOSPHERE%AIRAMT,                      &
                ATMOSPHERE%DLAYDLEV,                    &
                ATMOSPHERE%CLD_FLAG,                    &
                ATMOSPHERE%PCLD,                        &
                ATMOSPHERE%TAUCLD,                      &
                ATMOSPHERE%DECLD,                       &
                ATMOSPHERE%VSFC,                        &
                ATMOSPHERE%DVSFCDVLSFC,                 &
                ATMOSPHERE%DVSFCDVLSFC_,                &
                STAT = DEALLOC_STAT )
    
    IF ( DEALLOC_STAT /= 0 ) THEN
       PRINT*,'ERROR TRYING TO DEALLOCATE ATMOSPHERE'
       STOP
    ENDIF
    
  END SUBROUTINE CLEAR_PCRTM_ATMOSPHERE



END MODULE PCRTM_ATMOSPHERE_DEFINE





