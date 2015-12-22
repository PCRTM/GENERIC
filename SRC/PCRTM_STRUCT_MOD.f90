MODULE PCRTM_STRUCT_MOD

  USE PCRTM_TYPE_KIND
  USE PCRTM_ATM_ABSORPTION_DEFINE, ONLY : PCRTM_ATM_ABS_STRUCT_TYPE, &
                                          PCRTM_ATM_ABSORPTION_TYPE, &
                                          INIT_PCRTM_ATM_ABS_STRUCT, &
                                          INIT_PCRTM_ATM_ABSORPTION
  USE PCRTM_CLOUD_LUT_IO, ONLY : PCRTM_CLD_TABLE_DEF,INIT_PCRTM_ICECLD_GRID,INIT_PCRTM_WATCLD_GRID
  USE PCRTM_PC_SOLUTION,  ONLY : PCRTM_EOF_SOLUTION_TYPE, INIT_PCRTM_EOF_SOLUTION
  USE PCRTM_MATH_UTILITY, ONLY : FIND_BINARY_INDEX

  TYPE FRQ_INDEX_TYPE
     INTEGER, ALLOCATABLE :: I_CH(:)
  END TYPE FRQ_INDEX_TYPE

CONTAINS
  
  SUBROUTINE PCRTM_STRUCT_WN_CONVERT(FRQCH,            &
                                     NCH,              &  
                                     PCRTM_STND,       &
                                     ATM_ABS_COEF,     &
                                     ICE_GRID,         &
                                     WAT_GRID,         &
                                     EOF_SOLUTION,     &
                                     PCRTM_STND_,      &
                                     ATM_ABS_COEF_,    &
                                     ICE_GRID_,        &
                                     WAT_GRID_,        &
                                     EOF_SOLUTION_)

    TYPE(PCRTM_ATM_ABS_STRUCT_TYPE), INTENT(IN)   :: PCRTM_STND
    TYPE(PCRTM_ATM_ABSORPTION_TYPE), INTENT(IN)   :: ATM_ABS_COEF
    TYPE(PCRTM_CLD_TABLE_DEF),       INTENT(IN)   :: ICE_GRID
    TYPE(PCRTM_CLD_TABLE_DEF),       INTENT(IN)   :: WAT_GRID
    TYPE(PCRTM_EOF_SOLUTION_TYPE),   &
                        ALLOCATABLE, INTENT(IN)   :: EOF_SOLUTION(:)

    TYPE(PCRTM_ATM_ABS_STRUCT_TYPE), INTENT(OUT)  :: PCRTM_STND_
    TYPE(PCRTM_ATM_ABSORPTION_TYPE), INTENT(OUT)  :: ATM_ABS_COEF_
    TYPE(PCRTM_CLD_TABLE_DEF),       INTENT(OUT)  :: ICE_GRID_
    TYPE(PCRTM_CLD_TABLE_DEF),       INTENT(OUT)  :: WAT_GRID_
    TYPE(PCRTM_EOF_SOLUTION_TYPE),   &
                        ALLOCATABLE, INTENT(OUT)  :: EOF_SOLUTION_(:)
    INTEGER,                   INTENT(IN)   :: NCH
    REAL(DOUBLE),              INTENT(IN)   :: FRQCH(NCH)

    INTEGER,                   ALLOCATABLE  :: NCHBND(:)
    INTEGER                                 :: NBND, NBND_
    INTEGER,                   ALLOCATABLE  :: BINDX(:)
    TYPE(FRQ_INDEX_TYPE),      ALLOCATABLE  :: CHINDX(:)
    INTEGER                                 :: NM
    INTEGER                                 :: NS,NT,NS_,NT_
    INTEGER                                 :: IB
    INTEGER                                 :: ICH
    INTEGER                                 :: K, K1, K2
    INTEGER                                 :: ALLOC_STAT, DEALLOC_STAT

    NBND  = SIZE(EOF_SOLUTION)

    ALLOCATE(BINDX(NBND),NCHBND(NBND),STAT = ALLOC_STAT)
    IF (ALLOC_STAT /= 0) PRINT*, 'ERROR TRYING TO ALLOCATE BND OR BND INDEX'
    ALLOCATE(CHINDX(NBND),STAT = ALLOC_STAT)
    IF (ALLOC_STAT /= 0) PRINT*, 'ERROR TRYING TO ALLOCATE CHANNNEL INDEX'
    DO IB = 1, NBND
       ALLOCATE(CHINDX(IB)%I_CH(EOF_SOLUTION(IB)%NCHBND),STAT = ALLOC_STAT)
       IF (ALLOC_STAT /= 0) THEN
          PRINT*, 'ERROR TRYING TO ALLOCATE CHANNNEL INDEX OF BND', IB
       END IF
    END DO

    BINDX  = 0
    NCHBND = 0
    NBND_  = 0
    NM     = 0

    DO IB = 1, NBND
       CHINDX(IB)%I_CH = 0
       DO ICH = 1, NCH
!!$          IF ( FRQCH(ICH) .GE. EOF_SOLUTION(IB)%FRQCH(1) .AND.           &
!!$               FRQCH(ICH) .LE.                                           & 
!!$               EOF_SOLUTION(IB)%FRQCH(EOF_SOLUTION(IB)%NCHBND) )  THEN 
!!$             NCHBND(IB) = NCHBND(IB) + 1
!!$             CALL FIND_BINARY_INDEX(EOF_SOLUTION(IB)%FRQCH,FRQCH(ICH),K)
!!$             CHINDX(IB)%I_CH(NCHBND(IB)) = K
!!$          END IF
          DO K = 1, EOF_SOLUTION(IB)%NCHBND
             IF(FRQCH(ICH) .EQ. EOF_SOLUTION(IB)%FRQCH(K)) THEN
                NCHBND(IB) = NCHBND(IB) + 1
                CHINDX(IB)%I_CH(NCHBND(IB)) = K
             END IF
          END DO
       END DO
       IF (NCHBND(IB) .GT. 0) THEN
          NBND_ = NBND_ + 1
          BINDX(NBND_)  = IB 
          NM    = NM + EOF_SOLUTION(IB)%NREG
       END IF
    END DO

    IF ( NBND_ .EQ. 0 ) THEN
       PRINT*,'ERROR - INPUT CHANNEL FREQUENCY OUTSIDE OF THE FREQUENCY RANGE DEFINED BY THE SENSOR'
       STOP
    END IF

    ALLOCATE(EOF_SOLUTION_(NBND_),STAT = ALLOC_STAT)
    IF (ALLOC_STAT /= 0) PRINT*, 'ERROR TRYING TO ALLOCATE EOF_SOLUTION_'

    DO IB = 1, NBND_
       CALL INIT_PCRTM_EOF_SOLUTION(EOF_SOLUTION_(IB),              &
                                    EOF_SOLUTION(BINDX(IB))%NREG,   &
                                    EOF_SOLUTION(BINDX(IB))%NPCBND, &
                                    NCHBND(BINDX(IB)))
       EOF_SOLUTION_(IB)%FRQM     = EOF_SOLUTION(BINDX(IB))%FRQM
       EOF_SOLUTION_(IB)%RADSTD   = EOF_SOLUTION(BINDX(IB))%RADSTD
       EOF_SOLUTION_(IB)%COEFMEAN = EOF_SOLUTION(BINDX(IB))%COEFMEAN
       EOF_SOLUTION_(IB)%REGCOEF  = EOF_SOLUTION(BINDX(IB))%REGCOEF
       DO ICH = 1, NCHBND(BINDX(IB))
          EOF_SOLUTION_(IB)%FRQCH(ICH)     = EOF_SOLUTION(BINDX(IB))%FRQCH(CHINDX(BINDX(IB))%I_CH(ICH))
          EOF_SOLUTION_(IB)%RADSTDCH(ICH)  = EOF_SOLUTION(BINDX(IB))%RADSTDCH(CHINDX(BINDX(IB))%I_CH(ICH))
          EOF_SOLUTION_(IB)%RADMEANCH(ICH) = EOF_SOLUTION(BINDX(IB))%RADMEANCH(CHINDX(BINDX(IB))%I_CH(ICH))
          EOF_SOLUTION_(IB)%PC(:,ICH)      = EOF_SOLUTION(BINDX(IB))%PC(:,CHINDX(BINDX(IB))%I_CH(ICH))
       END DO
    END DO
   
    EOF_SOLUTION_(1)%INDX    = EOF_SOLUTION(BINDX(1))%INDX 
    IF (BINDX(1) .GT. 1) THEN
       K1 = EOF_SOLUTION(BINDX(1)-1)%NREG
       EOF_SOLUTION_(1)%INDX =  EOF_SOLUTION_(1)%INDX - EOF_SOLUTION(BINDX(1)-1)%INDX(K1)
    END IF
    DO IB = 2, NBND_
       K1 = EOF_SOLUTION(BINDX(IB)-1)%NREG
       K2 = EOF_SOLUTION_(IB-1)%NREG
       EOF_SOLUTION_(IB)%INDX  = EOF_SOLUTION(BINDX(IB))%INDX - EOF_SOLUTION(BINDX(IB)-1)%INDX(K1) &
                               + EOF_SOLUTION_(IB-1)%INDX(K2)
    END DO

    CALL INIT_PCRTM_ATM_ABS_STRUCT( PCRTM_STND_,            &
                                    PCRTM_STND%NLAY,        &
                                    PCRTM_STND%NT,          &
                                    PCRTM_STND%NMOL,        &
                                    NM )
       
    CALL INIT_PCRTM_ICECLD_GRID( ICE_GRID_, NM )
    CALL INIT_PCRTM_WATCLD_GRID( WAT_GRID_, NM )

    CALL INIT_PCRTM_ATM_ABSORPTION(ATM_ABS_COEF_, PCRTM_STND_)

    PCRTM_STND_%MOLINDX  =  PCRTM_STND%MOLINDX       
    PCRTM_STND_%SCALFAC  =  PCRTM_STND%SCALFAC
    PCRTM_STND_%GASSTD   =  PCRTM_STND%GASSTD
    PCRTM_STND_%TSTD     =  PCRTM_STND%TSTD
    PCRTM_STND_%FIXAMT   =  PCRTM_STND%FIXAMT
    PCRTM_STND_%PBND     =  PCRTM_STND%PBND

    DO IB  = 1, NBND_

       NS  = EOF_SOLUTION(BINDX(IB))%INDX(1)
       NT  = EOF_SOLUTION(BINDX(IB))%INDX(EOF_SOLUTION(BINDX(IB))%NREG)
       NS_ = EOF_SOLUTION_(IB)%INDX(1)
       NT_ = EOF_SOLUTION_(IB)%INDX(EOF_SOLUTION_(IB)%NREG)

       PCRTM_STND_%FRQ(NS_:NT_)          =  PCRTM_STND%FRQ(NS:NT) 
       PCRTM_STND_%NGAS(NS_:NT_)         =  PCRTM_STND%NGAS(NS:NT)
       PCRTM_STND_%IDGAS(:,NS_:NT_)      =  PCRTM_STND%IDGAS(:,NS:NT)
       
       ICE_GRID_%R_TAB(:,:,:,NS_:NT_)    =  ICE_GRID%R_TAB(:,:,:,NS:NT)
       ICE_GRID_%T_TAB(:,:,:,NS_:NT_)    =  ICE_GRID%T_TAB(:,:,:,NS:NT)
       WAT_GRID_%R_TAB(:,:,:,NS_:NT_)    =  WAT_GRID%R_TAB(:,:,:,NS:NT)
       WAT_GRID_%T_TAB(:,:,:,NS_:NT_)    =  WAT_GRID%T_TAB(:,:,:,NS:NT)

       DO K = 1,PCRTM_STND%NMOL
          ATM_ABS_COEF_%KGAS(K,NS_:NT_,:,:) =  ATM_ABS_COEF%KGAS(K,NS:NT,:,:)
       END DO
       ATM_ABS_COEF_%KH2O1(NS_:NT_,:,:)  =  ATM_ABS_COEF%KH2O1(NS:NT,:,:)
       ATM_ABS_COEF_%KH2O2(NS_:NT_,:,:)  =  ATM_ABS_COEF%KH2O2(NS:NT,:,:)
       ATM_ABS_COEF_%KH2O3(NS_:NT_,:,:)  =  ATM_ABS_COEF%KH2O3(NS:NT,:,:)

    END DO

    DEALLOCATE(BINDX,NCHBND,STAT = DEALLOC_STAT)
    IF (DEALLOC_STAT /= 0) PRINT*, 'ERROR TRYING TO DEALLOCATE BND OR BND INDEX'
    DO IB = 1, NBND
       DEALLOCATE(CHINDX(IB)%I_CH, STAT = DEALLOC_STAT)
       IF (DEALLOC_STAT /= 0) PRINT*, 'ERROR TRYING TO DEALLOCATE CHANNNEL INDEX'
    END DO
    DEALLOCATE(CHINDX,STAT = DEALLOC_STAT)
    IF (DEALLOC_STAT /= 0) PRINT*, 'ERROR TRYING TO ALLOCATE CHANNNEL INDEX'


  END SUBROUTINE PCRTM_STRUCT_WN_CONVERT
  


END MODULE PCRTM_STRUCT_MOD
