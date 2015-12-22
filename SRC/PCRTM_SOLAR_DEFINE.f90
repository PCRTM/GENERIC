Module PCRTM_Solar_define

  USE PCRTM_FILE_UTILITIES,only:GETLUN
  USE PCRTM_MATH_UTILITY,only:linear_order_interp
  USe pcrtm_type_kind
  USE PCRTM_ATM_ABSORPTION_DEFINE,only: PCRTM_ATM_ABS_STRUCT_TYPE

  TYPE PCRTM_SOLAR_DEF
     integer      :: StartWaveIndex
     REAL(SINGLE),ALLOCATABLE :: SolarSpectrum(:)
     REAL(DOUBLE),ALLOCATABLE :: BRDF(:),TBD(:),Tdb(:),Rdd(:) 
     REAL(SINGLE),ALLOCATABLE :: SolarRadup(:) 
  END TYPE PCRTM_SOLAR_DEF

Contains

  subroutine INIT_PCRTM_SOLAR_SOLUTION( PCRTM_STND,                        &
                                        PCRTM_SOLAR_SOLUTION)
    TYPE(PCRTM_ATM_ABS_STRUCT_TYPE),   INTENT(IN)  :: PCRTM_STND
    Type(PCRTM_SOLAR_DEF),             INTENT(out) :: PCRTM_SOLAR_SOLUTION
    
    integer  :: IO, n,im,StartWaveIndex
    real*8   :: solarspectrum(480001),frq(480001)

    integer :: alloc_stat


    do im = 1, pcrtm_stnd%nM
       IF (PCRTM_STND%FRQ(im).GT.1800) THEN
          StartWaveIndex = im 
          exit
       ENDIF
    end do
    PCRTM_SOLAR_SOLUTION%StartWaveIndex = StartWaveIndex

    nfrq = pcrtm_stnd%nM - StartWaveIndex + 1 

    allocate( PCRTM_SOLAR_SOLUTION%SolarSpectrum(NFRQ),         &
              PCRTM_SOLAR_SOLUTION%SolarRadUp(NFRQ),             &
              PCRTM_SOLAR_SOLUTION%BRDF(NFRQ),                   &
              PCRTM_SOLAR_SOLUTION%TDB(NFRQ),                    &
              PCRTM_SOLAR_SOLUTION%TBD(NFRQ),                    &
              PCRTM_SOLAR_SOLUTION%RDD(NFRQ),                    &
              STAT = ALLOC_STAT)
    IF ( ALLOC_STAT /= 0 ) THEN
       PRINT*,'ERROR TRYING TO ALLOCATE PCRTM_SOLAR_SOLUTION'
       STOP
    ENDIF
    PCRTM_SOLAR_SOLUTION%BRDF     =0
    PCRTM_SOLAR_SOLUTION%TDB      =0
    PCRTM_SOLAR_SOLUTION%TBD      =0
    PCRTM_SOLAR_SOLUTION%RDD      =0


    CALL GETLUN(IO)
    OPEN(IO,FILE='Solar_Spectrum_1800_3000.txt')
    Do n = 1,480001
       read(IO,*)frq(n),solarspectrum(n)
    end Do
    CLOSE(IO)

    CALL linear_order_interp(frq, solarspectrum, pcrtm_stnd%frq(StartWaveIndex:pcrtm_stnd%nM), &
                             PCRTM_SOLAR_SOLUTION%Solarspectrum, 480001, nfrq)
    

  end subroutine INIT_PCRTM_SOLAR_SOLUTION

  subroutine CLEAR_PCRTM_SOLAR_SOLUTION(PCRTM_SOLAR_SOLUTION)

    Type(PCRTM_SOLAR_DEF),             INTENT(inout) :: PCRTM_SOLAR_SOLUTION
    INTEGER                                          :: DEALLOC_STAT


    deallocate( PCRTM_SOLAR_SOLUTION%SolarSpectrum,          &
                PCRTM_SOLAR_SOLUTION%SolarRadUp,             &
                PCRTM_SOLAR_SOLUTION%BRDF,                   &
                PCRTM_SOLAR_SOLUTION%TDB,                    &
                PCRTM_SOLAR_SOLUTION%TBD,                    &
                PCRTM_SOLAR_SOLUTION%RDD,                    &
                STAT = DEALLOC_STAT)
    IF ( DEALLOC_STAT /= 0 ) THEN
       PRINT*,'ERROR TRYING TO DEALLOCATE PCRTM_SOLAR_SOLUTION'
       STOP
    ENDIF

  end subroutine CLEAR_PCRTM_SOLAR_SOLUTION

end Module PCRTM_SOLAR_DEFINE
