MODULE PCRTM_PC_SOLUTION
  USE PCRTM_TYPE_KIND
  USE PCRTM_FILE_UTILITIES
  USE PCRTM_RT_SOLUTION_DEFINE,ONLY : PCRTM_RT_SOLUTION_TYPE
  USE PCRTM_JACOBIAN
  USE PCRTM_CONSTANTS,ONLY    : WTDRYAIR, WTMOL
  USE PCRTM_MATH_UTILITY

  TYPE PCRTM_EOF_SOLUTION_TYPE
     LOGICAL                    :: CHDOMAIN     !DETERMINE IF OUTPUT CHANNEL RADIANCE
     LOGICAL                    :: BT_flag      !DETERMINE IF OUTPUT BRIGHTNESS TEMPERATURE
     INTEGER                    :: NREG         !NO. OF MONO FRQ IN EACH BAND
     INTEGER                    :: NCHBND       !NO. OF CHAN FRQ IN EACH BAND
     INTEGER                    :: NPCBND       !NO. OF EOFS IN EACH BAND
     INTEGER,      ALLOCATABLE  :: INDX(:)
     REAL(DOUBLE), ALLOCATABLE  :: FRQM(:)      !MONO FRQ
     REAL(DOUBLE), ALLOCATABLE  :: FRQCH(:)     !CHANNEL FRQ 
     REAL(SINGLE), ALLOCATABLE  :: PC(:,:)      !PCS FOR EACH BANDS
     REAL(SINGLE), ALLOCATABLE  :: RADSTDCH(:)  !NORMALIZATION FACTOR
     REAL(SINGLE), ALLOCATABLE  :: RADMEANCH(:) !MEAN
     ! PC REGR COEF FOR EACH BAND AND THE MONO SCALE FACTOR FOR EACH BAND
     REAL(SINGLE), ALLOCATABLE  :: REGCOEF(:,:)
     REAL(SINGLE), ALLOCATABLE  :: RADSTD(:)
     REAL(SINGLE), ALLOCATABLE  :: COEFMEAN(:)
     REAL(SINGLE), ALLOCATABLE  :: RADCH(:)
     REAL(SINGLE), ALLOCATABLE  :: RADPC(:)
     REAL(SINGLE), ALLOCATABLE  :: BTCH(:)
  END TYPE PCRTM_EOF_SOLUTION_TYPE

CONTAINS
  
  SUBROUTINE INIT_PCRTM_EOF_SOLUTION(EOF_SOLUTION,NREG,NPCBND,NCHBND)
    TYPE(PCRTM_EOF_SOLUTION_TYPE),&
                        INTENT(OUT)  :: EOF_SOLUTION 
    INTEGER,            INTENT(IN)   :: NREG
    INTEGER,            INTENT(IN)   :: NPCBND
    INTEGER,            INTENT(IN)   :: NCHBND

    INTEGER                          :: ALLOC_STAT

    EOF_SOLUTION%CHDOMAIN = .FALSE.
    EOF_SOLUTION%NREG     = NREG
    EOF_SOLUTION%NPCBND   = NPCBND
    EOF_SOLUTION%NCHBND   = NCHBND
    ALLOCATE(EOF_SOLUTION%INDX(NREG),                      &
             EOF_SOLUTION%FRQM(NREG),                      &
             EOF_SOLUTION%FRQCH(NCHBND),                   &
             EOF_SOLUTION%RADSTD(NREG),                    &
             EOF_SOLUTION%RADSTDCH(NCHBND),                &
             EOF_SOLUTION%REGCOEF(NREG,NPCBND),            &
             EOF_SOLUTION%PC(NPCBND,NCHBND),               &
             EOF_SOLUTION%RADMEANCH(NCHBND),               &
             EOF_SOLUTION%COEFMEAN(NPCBND),                &
             EOF_SOLUTION%RADCH(NCHBND),                   &
             EOF_SOLUTION%RADPC(NPCBND),                   &
             EOF_SOLUTION%BTCH(NCHBND),                    &
             STAT = ALLOC_STAT )                                  
    
    IF ( ALLOC_STAT /= 0 ) THEN
       PRINT*,'ERROR TRYING TO ALLOCATE EOF_SOLUTION'
       STOP
    ENDIF
    EOF_SOLUTION%RADSTD        =  0
    EOF_SOLUTION%RADSTDCH      =  0
    EOF_SOLUTION%REGCOEF       =  0
    EOF_SOLUTION%PC            =  0
    EOF_SOLUTION%RADMEANCH     =  0
  END SUBROUTINE INIT_PCRTM_EOF_SOLUTION

  SUBROUTINE CLEAR_PCRTM_EOF_SOLUTION(EOF_SOLUTION)
    TYPE(PCRTM_EOF_SOLUTION_TYPE), &
                        INTENT(INOUT) :: EOF_SOLUTION
    INTEGER                           :: DEALLOC_STAT

    DEALLOCATE(   EOF_SOLUTION%INDX,                       &
                  EOF_SOLUTION%FRQM,                       &
                  EOF_SOLUTION%FRQCH,                      &
                  EOF_SOLUTION%RADSTD,                     &
                  EOF_SOLUTION%RADSTDCH,                   &
                  EOF_SOLUTION%REGCOEF,                    &
                  EOF_SOLUTION%PC,                         &
                  EOF_SOLUTION%RADMEANCH,                  &
                  EOF_SOLUTION%COEFMEAN,                   &
                  EOF_SOLUTION%RADCH,                      &
                  EOF_SOLUTION%RADPC,                      &
                  EOF_SOLUTION%BTCH,                       &
                  STAT = DEALLOC_STAT )                                  
    IF ( DEALLOC_STAT /= 0 ) THEN
       PRINT*,'ERROR TRYING TO DEALLOCATE EOF_SOLUTION'
       STOP
    ENDIF

  END SUBROUTINE CLEAR_PCRTM_EOF_SOLUTION

  SUBROUTINE TRUNCATE_PCRTM_EOF_JACOB(EOF_SOLUTION,         &
                                      JACOB,                &
                                      K_PC,                 &
                                      NPCBND)
    TYPE(PCRTM_EOF_SOLUTION_TYPE),  &
                        allocatable,INTENT(INOUT) :: EOF_SOLUTION(:)     
    TYPE(PCRTM_PC_JACOBIAN_TYPE),   &
                        ALLOCATABLE,INTENT(INOUT) :: K_PC(:)

    LOGICAL,                        INTENT(IN)    :: JACOB
    INTEGER,            ALLOCATABLE,INTENT(IN)    :: NPCBND(:)
 
    INTEGER                                       :: NLAY, NB
    
    NB = SIZE(EOF_SOLUTION)

    DO IB = 1,NB
       CALL TRUNCATE_PCRTM_EOF_SOLUTION(EOF_SOLUTION(IB), NPCBND(IB))
    END DO

    IF(JACOB) THEN
       NLAY = size(K_PC(1)%K_TLAY,2)
       DO IB = 1,NB
          CALL CLEAR_PCRTM_JACOB_PC( K_PC(IB) )
          CALL INIT_PCRTM_JACOB_PC( K_PC(IB), EOF_SOLUTION(IB)%NPCBND, NLAY )
       END DO
    END IF

  END SUBROUTINE TRUNCATE_PCRTM_EOF_JACOB

  SUBROUTINE TRUNCATE_PCRTM_EOF_SOLUTION(EOF_SOLUTION, NPCBND)
    TYPE(PCRTM_EOF_SOLUTION_TYPE),&
                        INTENT(INOUT)  :: EOF_SOLUTION     
    INTEGER,            INTENT(IN)     :: NPCBND

    TYPE(PCRTM_EOF_SOLUTION_TYPE)      :: EOF_SOLUTION_
    INTEGER :: NREG
    INTEGER :: NCHBND
    INTEGER :: DEALLOC_STAT

    IF( NPCBND .GT. EOF_SOLUTION%NPCBND) THEN
       PRINT*,'ERROR --- NUMBER OF EOFS AFTER TRUNCATION CAN NOT EXCEED THE ORIGINAL NUMBER'
       PRINT*,'MAXIMUM EOF NUMBER ALLOWED', EOF_SOLUTION%NPCBND
       STOP
    END IF

    NREG   = EOF_SOLUTION%NREG
    NCHBND = EOF_SOLUTION%NCHBND
    
    CALL INIT_PCRTM_EOF_SOLUTION(EOF_SOLUTION_,NREG,NPCBND,NCHBND)
    
    EOF_SOLUTION_%CHDOMAIN = EOF_SOLUTION%CHDOMAIN
    EOF_SOLUTION_%BT_FLAG  = EOF_SOLUTION%BT_FLAG 
    EOF_SOLUTION_%INDX     = EOF_SOLUTION%INDX
    EOF_SOLUTION_%FRQM     = EOF_SOLUTION%FRQM
    EOF_SOLUTION_%FRQCH    = EOF_SOLUTION%FRQCH 
    EOF_SOLUTION_%RADSTD   = EOF_SOLUTION%RADSTD
    EOF_SOLUTION_%RADSTDCH = EOF_SOLUTION%RADSTDCH
    EOF_SOLUTION_%REGCOEF  = EOF_SOLUTION%REGCOEF(:,1:NPCBND)
    EOF_SOLUTION_%PC       = EOF_SOLUTION%PC(1:NPCBND,:)
    EOF_SOLUTION_%RADMEANCH= EOF_SOLUTION%RADMEANCH
    EOF_SOLUTION_%COEFMEAN = EOF_SOLUTION%COEFMEAN(1:NPCBND)
    EOF_SOLUTION_%RADCH    = EOF_SOLUTION%RADCH
    EOF_SOLUTION_%RADPC    = EOF_SOLUTION%RADPC(1:NPCBND)
    EOF_SOLUTION_%BTCH     = EOF_SOLUTION%BTCH

    DEALLOCATE(EOF_SOLUTION%INDX,                       &
               EOF_SOLUTION%FRQM,                       &
               EOF_SOLUTION%FRQCH,                      &
               EOF_SOLUTION%RADSTD,                     &
               EOF_SOLUTION%RADSTDCH,                   &
               EOF_SOLUTION%REGCOEF,                    &
               EOF_SOLUTION%PC,                         &
               EOF_SOLUTION%RADMEANCH,                  &
               EOF_SOLUTION%COEFMEAN,                   &
               EOF_SOLUTION%RADCH,                      &
               EOF_SOLUTION%RADPC,                      &
               EOF_SOLUTION%BTCH,                       &
               STAT = DEALLOC_STAT )                                  
    IF ( DEALLOC_STAT /= 0 ) THEN
       PRINT*,'ERROR TRYING TO DEALLOCATE EOF_SOLUTION IN TRUNCATE_EOF_SOLUTION()'
       STOP
    ENDIF

    CALL INIT_PCRTM_EOF_SOLUTION(EOF_SOLUTION,NREG,NPCBND,NCHBND)

    EOF_SOLUTION = EOF_SOLUTION_

    DEALLOCATE(EOF_SOLUTION_%INDX,                       &
               EOF_SOLUTION_%FRQM,                       &
               EOF_SOLUTION_%FRQCH,                      &
               EOF_SOLUTION_%RADSTD,                     &
               EOF_SOLUTION_%RADSTDCH,                   &
               EOF_SOLUTION_%REGCOEF,                    &
               EOF_SOLUTION_%PC,                         &
               EOF_SOLUTION_%RADMEANCH,                  &
               EOF_SOLUTION_%COEFMEAN,                   &
               EOF_SOLUTION_%RADCH,                      &
               EOF_SOLUTION_%RADPC,                      &
               EOF_SOLUTION_%BTCH,                       &
               STAT = DEALLOC_STAT )                                  
    IF ( DEALLOC_STAT /= 0 ) THEN
       PRINT*,'ERROR TRYING TO DEALLOCATE EOF_SOLUTION_ IN TRUNCATE_EOF_SOLUTION'
       STOP
    ENDIF

  END SUBROUTINE TRUNCATE_PCRTM_EOF_SOLUTION


  SUBROUTINE RD_SENSOR_BND_INFO( EOF_SOLUTION, PARFILE )
    TYPE(PCRTM_EOF_SOLUTION_TYPE), &
                        ALLOCATABLE,INTENT(OUT) :: EOF_SOLUTION(:) 
    CHARACTER(160),                 INTENT(IN)  :: PARFILE

    INTEGER :: IUPAR
    INTEGER :: NBND
    INTEGER :: NREG, NPCBND, NCHBND
    INTEGER :: IB, ALLOC_STAT

    CALL GETLUN(IUPAR)

    OPEN(UNIT=IUPAR,FILE=PARFILE,FORM='UNFORMATTED')

    PRINT*, PARFILE
    READ(IUPAR) NBND
    PRINT*, 'NUMBER OF BANDS', NBND

    ALLOCATE(EOF_SOLUTION(NBND),STAT = ALLOC_STAT)
    IF ( ALLOC_STAT /= 0 ) THEN
       PRINT*,'ERROR TRYING TO ALLOCATE EOF_SOLUTION'
       STOP
    ENDIF
    
    DO IB=1,NBND
       READ(IUPAR) NREG, NPCBND, NCHBND
       PRINT*, IB, NREG, NPCBND, NCHBND
       CALL INIT_PCRTM_EOF_SOLUTION(EOF_SOLUTION(IB),NREG,NPCBND,NCHBND)
       READ(IUPAR)(EOF_SOLUTION(IB)%INDX(I),I=1,NREG)
       READ(IUPAR)(EOF_SOLUTION(IB)%FRQM(I),I=1,NREG)
       READ(IUPAR)(EOF_SOLUTION(IB)%FRQCH(I),I=1,NCHBND)
       READ(IUPAR)(EOF_SOLUTION(IB)%RADSTD(I),I=1,NREG)
       READ(IUPAR)(EOF_SOLUTION(IB)%RADSTDCH(I),I=1,NCHBND)
       READ(IUPAR)((EOF_SOLUTION(IB)%REGCOEF(I,J),I=1,NREG),J=1,NPCBND)
       READ(IUPAR)((EOF_SOLUTION(IB)%PC(J,I),J=1,NPCBND),I=1,NCHBND)
       READ(IUPAR)(EOF_SOLUTION(IB)%RADMEANCH(I),I=1,NCHBND)
       DO I=1,NPCBND
          EOF_SOLUTION(IB)%COEFMEAN(I)=                                     & 
          DOT_PRODUCT(EOF_SOLUTION(IB)%PC(I,:),EOF_SOLUTION(IB)%RADMEANCH)
       ENDDO
    END DO
    
    CLOSE(IUPAR)

  END SUBROUTINE RD_SENSOR_BND_INFO



  SUBROUTINE PRED2PCSCORE(RT_SOLUTION, EOF_SOLUTION)

    TYPE(PCRTM_RT_SOLUTION_TYPE), INTENT(IN)    :: RT_SOLUTION    
    TYPE(PCRTM_EOF_SOLUTION_TYPE),INTENT(INOUT) :: EOF_SOLUTION

    INTEGER :: I,K,J
    REAL(SINGLE), ALLOCATABLE :: TMP(:)
    INTEGER :: ALLOC_STAT, DEALLOC_STAT

    ALLOCATE(TMP(EOF_SOLUTION%NREG), STAT = ALLOC_STAT)
    IF(ALLOC_STAT /= 0 ) THEN
       PRINT*,'ERROR TRYING TO ALLOCATE TMP VARIABLE IN PRED2PCSCORE SUBROUTINE'
       STOP
    ENDIF

    DO K = 1,EOF_SOLUTION%NREG
       J = EOF_SOLUTION%INDX(K)
       TMP(K) = RT_SOLUTION%RADUP(J)* EOF_SOLUTION%RADSTD(K)
    END DO

    EOF_SOLUTION%RADPC = 0
    DO I=1,EOF_SOLUTION%NPCBND
       DO K=1,EOF_SOLUTION%NREG          
          EOF_SOLUTION%RADPC(I) = EOF_SOLUTION%RADPC(I)                     &
                                + TMP(K)                                    &
                                * EOF_SOLUTION%REGCOEF(K,I)          
       ENDDO
    ENDDO
    EOF_SOLUTION%RADPC = EOF_SOLUTION%RADPC - EOF_SOLUTION%COEFMEAN 

    DEALLOCATE(TMP, STAT = DEALLOC_STAT)
    IF(DEALLOC_STAT /= 0 ) THEN
       PRINT*,'ERROR TRYING TO DEALLOCATE TMP VARIABLE IN PRED2PCSCORE SUBROUTINE'
       STOP
    ENDIF    
   
  END SUBROUTINE PRED2PCSCORE


  SUBROUTINE PRED2PC_JACOB(TMP_M, EOF_SOLUTION, TMP_PC)

    REAL(SINGLE),                 INTENT(IN):: TMP_M(:)    
    TYPE(PCRTM_EOF_SOLUTION_TYPE),INTENT(IN):: EOF_SOLUTION
    REAL(SINGLE),ALLOCATABLE,  INTENT(INOUT):: TMP_PC(:)
 
    REAL(SINGLE), ALLOCATABLE               :: TMP(:)
    INTEGER :: I,K,J
    INTEGER :: NREG, NPCBND
    INTEGER :: ALLOC_STAT, DEALLOC_STAT

    NREG  = EOF_SOLUTION%NREG
    ALLOCATE(TMP(NREG),                                              &
             STAT = ALLOC_STAT)
    IF(ALLOC_STAT /= 0 ) THEN
       PRINT*,'ERROR TRYING TO ALLOCATE TMP VARIABLE IN PRED2PC_JACOB SUBROUTINE'
       STOP
    ENDIF

    DO K=1,NREG
       J = EOF_SOLUTION%INDX(K)
       TMP(K) = TMP_M(J)*EOF_SOLUTION%RADSTD(K)
    END DO
    
    NPCBND = EOF_SOLUTION%NPCBND
 
    TMP_PC = 0
    DO I=1,NPCBND
       DO K=1,NREG
          TMP_PC(I) = TMP_PC(I) + TMP(K)                              &
                                * EOF_SOLUTION%REGCOEF(K,I)           
       ENDDO
    ENDDO
    
    DEALLOCATE(TMP, STAT = DEALLOC_STAT)
    IF(DEALLOC_STAT /= 0 ) THEN
       PRINT*,'ERROR TRYING TO DEALLOCATE TMP VARIABLE IN PRED2PC_JACOB SUBROUTINE'
       STOP
    ENDIF  

  END SUBROUTINE PRED2PC_JACOB


  SUBROUTINE PCSCORE2CH(EOF_SOLUTION)
    
    TYPE(PCRTM_EOF_SOLUTION_TYPE),INTENT(INOUT) :: EOF_SOLUTION
   
    INTEGER :: I, NPCBND, NCHBND, J

    NPCBND = EOF_SOLUTION%NPCBND
    NCHBND = EOF_SOLUTION%NCHBND

    EOF_SOLUTION%RADCH = 0
    DO I=1,NCHBND
       DO J = 1, NPCBND
          EOF_SOLUTION%RADCH(I) = EOF_SOLUTION%RADCH(I) +             &
                       EOF_SOLUTION%PC(J,I)*EOF_SOLUTION%RADPC(J)
       END DO
    ENDDO

!!$    DO I=1,NCHBND
!!$       EOF_SOLUTION%RADCH(I)                                          &
!!$                      = DOT_PRODUCT( EOF_SOLUTION%PC(1:NPCBND,I),     &
!!$                                     EOF_SOLUTION%RADPC(1:NPCBND) )
!!$    ENDDO

    EOF_SOLUTION%RADCH = ( EOF_SOLUTION%RADCH                         & 
                           + EOF_SOLUTION%RADMEANCH )                 &
                           * EOF_SOLUTION%RADSTDCH

    
    IF(EOF_SOLUTION%BT_FLAG) THEN
       EOF_SOLUTION%BTCH = Bright_Temp(EOF_SOLUTION%frqch, &
                                       EOF_SOLUTION%Radch, &
                                       EOF_SOLUTION%nchbnd)
    END IF

  END SUBROUTINE PCSCORE2CH

  
  SUBROUTINE PCRTM_JACOBIAN_PC_SOLUTION(K_M,                &
                                        DLAYDLEV,           &
                                        AIRAMT,             &
                                        NTOP,               &
                                        NBOT,               &
                                        EOF_SOLUTION,       &
                                        K_PC)

    TYPE(PCRTM_EOF_SOLUTION_TYPE),  INTENT(IN)   :: EOF_SOLUTION
    TYPE(PCRTM_PC_JACOBIAN_TYPE),   INTENT(INOUT):: K_PC
    TYPE(PCRTM_NM_JACOBIAN_TYPE),   INTENT(IN)   :: K_M

    INTEGER,                  INTENT(IN)   :: NTOP,NBOT
    REAL(SINGLE),             INTENT(IN)   :: DLAYDLEV(:)
    REAL(SINGLE),             INTENT(IN)   :: AIRAMT(:)

    INTEGER                                :: NM, NREG,NPCBND

    REAL(SINGLE), ALLOCATABLE              :: TMP_M(:)
    REAL(SINGLE), ALLOCATABLE              :: TMP_PC(:)


    INTEGER                                :: ALLOC_STAT, DEALLOC_STAT
    INTEGER                                :: L, N

    NM = SIZE(K_M%DRUPDTS)
    ALLOCATE(TMP_M(NM),STAT=ALLOC_STAT)
    IF(ALLOC_STAT /= 0 ) THEN
       PRINT*,'ERROR TRYING TO ALLOCATE TMP_M VARIABLE IN JACOBIAN_PC_SOLUTION SUBROUTINE'
       STOP
    ENDIF

    NPCBND = EOF_SOLUTION%NPCBND
    ALLOCATE(TMP_PC(NPCBND),STAT=ALLOC_STAT)
    IF(ALLOC_STAT /= 0 ) THEN
       PRINT*,'ERROR TRYING TO ALLOCATE TMP_PC VARIABLE IN JACOBIAN_PC_SOLUTION SUBROUTINE'
       STOP
    ENDIF 
    
    NREG  = EOF_SOLUTION%NREG
    
    DO L=NTOP,NBOT

       TMP_M = K_M%DRUPDT(:,L)
       CALL PRED2PC_JACOB(TMP_M, EOF_SOLUTION, TMP_PC)
       K_PC%K_TLAY(:,L) = TMP_PC

       TMP_M = K_M%DRUPDH2O(:,L)
       CALL PRED2PC_JACOB(TMP_M, EOF_SOLUTION, TMP_PC)
       K_PC%K_GASLAY(:,L,1) = TMP_PC*AIRAMT(L)/1.E3*WTDRYAIR/WTMOL(1)
       
       DO N = 1,15
          TMP_M = K_M%DRUPDGAS(n,:,L)
          CALL PRED2PC_JACOB(TMP_M, EOF_SOLUTION, TMP_PC)
          K_PC%K_GASLAY(:,L,n) = TMP_PC*AIRAMT(L)/1.E6
       END DO
          
!!$       TMP_M = K_M%DRUPDGAS(3,:,L)
!!$       CALL PRED2PC_JACOB(TMP_M, EOF_SOLUTION, TMP_PC)
!!$       K_PC%K_GASLAY(:,L,3) = TMP_PC*AIRAMT(L)/1.E6
!!$
!!$       TMP_M = K_M%DRUPDGAS(2,:,L)
!!$       CALL PRED2PC_JACOB(TMP_M, EOF_SOLUTION, TMP_PC)
!!$       K_PC%K_GASLAY(:,L,2) = TMP_PC*AIRAMT(L)/1.E6
!!$
!!$       TMP_M = K_M%DRUPDGAS(4,:,L)
!!$       CALL PRED2PC_JACOB(TMP_M, EOF_SOLUTION, TMP_PC)
!!$       K_PC%K_GASLAY(:,L,4) = TMP_PC*AIRAMT(L)/1.E6
!!$       
!!$       TMP_M = K_M%DRUPDGAS(5,:,L)
!!$       CALL PRED2PC_JACOB(TMP_M, EOF_SOLUTION, TMP_PC)
!!$       K_PC%K_GASLAY(:,L,5) = TMP_PC*AIRAMT(L)/1.E6
!!$       
!!$       TMP_M = K_M%DRUPDGAS(6,:,L)
!!$       CALL PRED2PC_JACOB(TMP_M, EOF_SOLUTION, TMP_PC)
!!$       K_PC%K_GASLAY(:,L,6) = TMP_PC*AIRAMT(L)/1.E6

    END DO

    DEALLOCATE(TMP_PC, STAT = DEALLOC_STAT)
    IF(DEALLOC_STAT /= 0 ) THEN
       PRINT*,'ERROR TRYING TO DEALLOCATE TMP_PC VARIABLE IN JACOBIAN_PC_SOLUTION SUBROUTINE'
       STOP
    ENDIF

    TMP_M = K_M%DRUPDTS
    CALL PRED2PC_JACOB(TMP_M, EOF_SOLUTION, K_PC%K_TS)
    TMP_M = K_M%DRUPDEM
    CALL PRED2PC_JACOB(TMP_M, EOF_SOLUTION, K_PC%K_EM)

    K_PC%K_T(:,NTOP)    = K_PC%K_TLAY(:,NTOP)*DLAYDLEV(NTOP)
    K_PC%K_GASLEV(:,NTOP,:)  = K_PC%K_GASLAY(:,NTOP,:)*DLAYDLEV(NTOP)


    K_PC%K_T(:,NBOT+1)  = K_PC%K_TLAY(:,NBOT)*(1.0-DLAYDLEV(NBOT))
    K_PC%K_GASLEV(:,NBOT+1,:)= K_PC%K_GASLAY(:,NBOT,:)*(1.0-DLAYDLEV(NBOT))


    K_PC%K_T(:,NBOT+1)  = K_PC%K_T(:,NBOT+1)*DLAYDLEV(NBOT+1)
    K_PC%K_GASLEV(:,NBOT+1,:)= K_PC%K_GASLEV(:,NBOT+1,:)*DLAYDLEV(NBOT+1)


    DO L= NTOP+1,NBOT

       K_PC%K_T(:,L)= K_PC%K_TLAY(:,L-1)*(1.0-DLAYDLEV(L-1))+             &
               K_PC%K_TLAY(:,L)*DLAYDLEV(L) 

       DO N = 1,15
          K_PC%K_GASLEV(:,L,N)= K_PC%K_GASLAY(:,L-1,N)*(1.0-DLAYDLEV(L-1))+         &
               K_PC%K_GASLAY(:,L,N)*DLAYDLEV(L) 
       END DO
          
!!$       K_PC%K_H2O(:,L)= K_PC%K_H2OLAY(:,L-1)*(1.0-DLAYDLEV(L-1))+         &
!!$               K_PC%K_H2OLAY(:,L)*DLAYDLEV(L) 
!!$          
!!$       K_PC%K_CO2(:,L)= K_PC%K_CO2LAY(:,L-1)*(1.0-DLAYDLEV(L-1))+         &
!!$               K_PC%K_CO2LAY(:,L)*DLAYDLEV(L) 
!!$          
!!$       K_PC%K_O3(:,L)=  K_PC%K_O3LAY(:,L-1)*(1.0-DLAYDLEV(L-1))+          &
!!$               K_PC%K_O3LAY(:,L)*DLAYDLEV(L) 
!!$          
!!$       K_PC%K_N2O(:,L)= K_PC%K_N2OLAY(:,L-1)*(1.0-DLAYDLEV(L-1))+         &
!!$               K_PC%K_N2OLAY(:,L)*DLAYDLEV(L) 
!!$          
!!$       K_PC%K_CO(:,L)= K_PC%K_COLAY(:,L-1)*(1.0-DLAYDLEV(L-1))+           &
!!$                  K_PC%K_COLAY(:,L)*DLAYDLEV(L) 
!!$          
!!$       K_PC%K_CH4(:,L)= K_PC%K_CH4LAY(:,L-1)*(1.0-DLAYDLEV(L-1))+         &
!!$                  K_PC%K_CH4LAY(:,L)*DLAYDLEV(L) 
          
    ENDDO

    DEALLOCATE(TMP_M, STAT = DEALLOC_STAT)
    IF(DEALLOC_STAT /= 0 ) THEN
       PRINT*,'ERROR TRYING TO DEALLOCATE TMP VARIABLE IN JACOBIAN_PC_SOLUTION SUBROUTINE'
       STOP
    ENDIF 

  END SUBROUTINE PCRTM_JACOBIAN_PC_SOLUTION

  
  SUBROUTINE PCRTM_JACOBIAN_CH_SOLUTION(NTOP,               &
                                        NBOT,               &
                                        EOF_SOLUTION,       &
                                        K_PC,               &
                                        K_CH)

    TYPE(PCRTM_EOF_SOLUTION_TYPE),  INTENT(IN)   :: EOF_SOLUTION
    TYPE(PCRTM_CH_JACOBIAN_TYPE),   INTENT(INOUT):: K_CH
    TYPE(PCRTM_PC_JACOBIAN_TYPE),   INTENT(IN)   :: K_PC
    INTEGER,                        INTENT(IN)   :: NTOP,NBOT

    REAL(SINGLE), ALLOCATABLE              :: K_T(:,:)
    REAL(SINGLE), ALLOCATABLE              :: K_GASLEV(:,:,:)
!!$    REAL(SINGLE), ALLOCATABLE              :: K_CO2(:,:), K_CH4(:,:), K_N2O(:,:)
!!$    REAL(SINGLE), ALLOCATABLE              :: K_CO(:,:), K_O3(:,:),K_H2O(:,:)

    REAL(SINGLE), ALLOCATABLE              :: PC(:,:), tmp(:)

    INTEGER                                :: NPCBND, NCHBND

    INTEGER                                :: ALLOC_STAT, DEALLOC_STAT
    INTEGER                                :: I,L,N


    NPCBND = EOF_SOLUTION%NPCBND
    NCHBND = EOF_SOLUTION%NCHBND

    ALLOCATE(K_T  (NPCBND, NBOT+1),                                         &
             K_GASLEV(NPCBND, NBOT+1,15),                                   &
!!$             K_H2O(NPCBND, NBOT+1),                                         &
!!$             K_CO2(NPCBND, NBOT+1),                                         &
!!$             K_O3 (NPCBND, NBOT+1),                                         & 
!!$             K_CO (NPCBND, NBOT+1),                                         &
!!$             K_CH4(NPCBND, NBOT+1),                                         &
!!$             K_N2O(NPCBND, NBOT+1),                                         &
             STAT = ALLOC_STAT)
    IF(ALLOC_STAT /= 0 ) THEN
       PRINT*,'ERROR TRYING TO ALLOCATE K-JACOBIAN VARIABLES IN JACOBIAN_PC_SOLUTION SUBROUTINE'
       STOP
    ENDIF

    K_T(:,NTOP:NBOT+1)   = K_PC%K_T(:,NTOP:NBOT+1)
    K_GASLEV(:,NTOP:NBOT+1,:) = K_PC%K_GASLEV(:,NTOP:NBOT+1,:)

!!$    K_H2O(:,NTOP:NBOT+1) = K_PC%K_H2O(:,NTOP:NBOT+1)
!!$    K_CO2(:,NTOP:NBOT+1) = K_PC%K_CO2(:,NTOP:NBOT+1)
!!$    K_O3(:,NTOP:NBOT+1)  = K_PC%K_O3(:,NTOP:NBOT+1)
!!$    K_CO(:,NTOP:NBOT+1)  = K_PC%K_CO(:,NTOP:NBOT+1)
!!$    K_CH4(:,NTOP:NBOT+1) = K_PC%K_CH4(:,NTOP:NBOT+1)
!!$    K_N2O(:,NTOP:NBOT+1) = K_PC%K_N2O(:,NTOP:NBOT+1)


    ALLOCATE(PC(NPCBND, NCHBND), STAT = ALLOC_STAT)
    IF(ALLOC_STAT /= 0 ) THEN
       PRINT*,'ERROR TRYING TO ALLOCATE PC VARIABLE IN JACOBIAN_PC_SOLUTION SUBROUTINE'
       STOP
    ENDIF
    
    DO I = 1, NPCBND
       PC(I,:) = EOF_SOLUTION%PC(I,:)* EOF_SOLUTION%RADSTDCH
    END DO

    DO L =  NTOP, NBOT+1
       
       DO I=1,NCHBND


          K_CH%R_TLEV(I,L)     = DOT_PRODUCT( PC(:,I), K_T(:,L) )

          DO N = 1,15
             K_CH%R_GASLEV(I,L,N) = DOT_PRODUCT( PC(:,I), K_GASLEV(:,L,N) )
          END DO
             
!!$          K_CH%R_GASLEV(I,L,1) = DOT_PRODUCT( PC(:,I), K_H2O(:,L) )
!!$          
!!$          K_CH%R_GASLEV(I,L,2) = DOT_PRODUCT( PC(:,I), K_CO2(:,L) )
!!$          
!!$          K_CH%R_GASLEV(I,L,3) = DOT_PRODUCT( PC(:,I), K_O3(:,L) )
!!$          
!!$          K_CH%R_GASLEV(I,L,4) = DOT_PRODUCT( PC(:,I), K_N2O(:,L) )
!!$          
!!$          K_CH%R_GASLEV(I,L,5) = DOT_PRODUCT( PC(:,I), K_CO(:,L) )
!!$          
!!$          K_CH%R_GASLEV(I,L,6) = DOT_PRODUCT( PC(:,I), K_CH4(:,L) )
          
       END DO
       
    END DO
    
    
    DO I=1,NCHBND
       K_CH%R_EM(I)  = DOT_PRODUCT( PC(:,I), K_PC%K_EM )
       K_CH%R_TS(I)  = DOT_PRODUCT( PC(:,I), K_PC%K_TS )
    END DO
    
    
    
    DEALLOCATE(PC, STAT = DEALLOC_STAT)
    IF(DEALLOC_STAT /= 0 ) THEN
       PRINT*,'ERROR TRYING TO DEALLOCATE PC VARIABLE IN JACOBIAN_PC_SOLUTION SUBROUTINE'
       STOP
    ENDIF

 
    DEALLOCATE(K_T,                                           &
               K_GASLEV,                                      &
!!$               K_H2O,                                         &
!!$               K_CO2,                                         &
!!$               K_O3,                                          & 
!!$               K_CO,                                          &
!!$               K_CH4,                                         &
!!$               K_N2O,                                         &
               STAT = DEALLOC_STAT)
    IF(DEALLOC_STAT /= 0 ) THEN
       PRINT*,'ERROR TRYING TO DEALLOCATE K-JACOBIAN VARIABLES IN JACOBIAN_PC_SOLUTION SUBROUTINE'
       STOP
    ENDIF
    
    IF(K_CH%JACOB_BTCH) THEN

       ALLOCATE(tmp(NCHBND), STAT = ALLOC_STAT)
       IF(ALLOC_STAT /= 0 ) THEN
          PRINT*,'ERROR TRYING TO ALLOCATE tmp VARIABLE IN JACOBIAN_PC_SOLUTION SUBROUTINE'
          STOP
       ENDIF
       
       tmp = dBTdRad(EOF_SOLUTION%frqch,EOF_SOLUTION%Radch,Nchbnd)
      
       DO L =  NTOP, NBOT+1

          DO I=1,NCHBND

             K_CH%BT_TLEV(I,L)     = tmp(I)*K_CH%R_TLEV(I,L)
             
             DO N= 1,15
                K_CH%BT_GASLEV(I,L,N) = tmp(I)*K_CH%R_GASLEV(I,L,N)
             END DO
             
!!$             K_CH%BT_GASLEV(I,L,1) = tmp(I)*K_CH%R_GASLEV(I,L,1)
!!$             
!!$             K_CH%BT_GASLEV(I,L,2) = tmp(I)*K_CH%R_GASLEV(I,L,2)
!!$             
!!$             K_CH%BT_GASLEV(I,L,3) = tmp(I)*K_CH%R_GASLEV(I,L,3)
!!$             
!!$             K_CH%BT_GASLEV(I,L,4) = tmp(I)*K_CH%R_GASLEV(I,L,4)
!!$             
!!$             K_CH%BT_GASLEV(I,L,5) = tmp(I)*K_CH%R_GASLEV(I,L,5)
!!$             
!!$             K_CH%BT_GASLEV(I,L,6) = tmp(I)*K_CH%R_GASLEV(I,L,6)
             
          END DO
          
       END DO

           
       DO I=1,NCHBND
          K_CH%BT_EM(I)  = tmp(I)*K_CH%R_EM(I)
          K_CH%BT_TS(I)  = tmp(I)*K_CH%R_TS(I)
       END DO
       
       DEALLOCATE(tmp, STAT = DEALLOC_STAT)
       IF(DEALLOC_STAT /= 0 ) THEN
          PRINT*,'ERROR TRYING TO DEALLOCATE tmp VARIABLE IN JACOBIAN_PC_SOLUTION SUBROUTINE'
          STOP
       ENDIF

    END IF



  END SUBROUTINE PCRTM_JACOBIAN_CH_SOLUTION


END MODULE PCRTM_PC_SOLUTION
