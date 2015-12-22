MODULE PCRTM_FORWARD_MODEL

    USE PCRTM_ATMOSPHERE_LAYER
    USE PCRTM_ATM_ABSORPTION
    USE PCRTM_JACOBIAN
    USE PCRTM_RT_SOLUTION_DEFINE,ONLY : PCRTM_RT_SOLUTION_TYPE
    USE PCRTM_CALC_RAD
    USE PCRTM_CALC_CLOUD
    USE PCRTM_SOLAR_RT
    USE PCRTM_PC_SOLUTION
    USE PCRTM_TR_Solution

CONTAINS
  SUBROUTINE  PCRTM_FORWARD_RT_M( GEOMETRY,      &
                                  ATM,           &
                                  NCLD,          &
                                  CLD,           &
                                  ATM_ABS,       &
                                  PCRTM_STND,    &
                                  ICE_GRID,      &
                                  WAT_GRID,      &
                                  RT_SOLUTION,   &
                                  SOLAR_TAB,     &
                                  SOLAR_SOLUTION,&
                                  K_M,           &
                                  K_CLD)

    TYPE(PCRTM_GEOMETRY_TYPE),   INTENT(INOUT):: GEOMETRY
    TYPE(PCRTM_ATM_ABS_STRUCT_TYPE),INTENT(IN):: PCRTM_STND
    TYPE(PCRTM_ATMOSPHERE_TYPE), INTENT(INOUT):: ATM
    TYPE(PCRTM_ATM_ABSORPTION_TYPE),INTENT(IN):: ATM_ABS
    TYPE(PCRTM_NM_JACOBIAN_TYPE),INTENT(INOUT):: K_M
   

    INTEGER,                     INTENT(IN)   :: NCLD
    TYPE(PCRTM_RT_SOLUTION_TYPE),INTENT(INOUT):: RT_SOLUTION
    
    TYPE(PCRTM_SOLAR_DEF),       INTENT(INOUT):: SOLAR_SOLUTION
    type(PCRTM_SOLAR_LUT_DEF),   intent(in)    :: solar_tab(2)

    TYPE(PCRTM_CLD_TABLE_DEF),   INTENT(IN)   :: ICE_GRID
    TYPE(PCRTM_CLD_TABLE_DEF),   INTENT(IN)   :: WAT_GRID
    
    INTEGER                                   :: IM, K, I1
    TYPE(PCRTM_CLOUD_TYPE),       INTENT(OUT) :: CLD(NCLD)
    TYPE(PCRTM_CLD_JACOBIAN_TYPE),INTENT(OUT) :: K_CLD(NCLD)
    

    CALL PCRTM_GEOMETRY_INFO( GEOMETRY,             &
                        PCRTM_STND%PBND,            &
                        PCRTM_STND%NLEV )


    CALL LAYERAVG( GEOMETRY, ATM, PCRTM_STND )

    
    CALL CALCODETC( ATM,                            &
                    ATM_ABS,                        &
                    PCRTM_STND,                     &
                    GEOMETRY,                       &
                    RT_SOLUTION,                    &
                    K_M )

    IF (NCLD .GT. 0) THEN 
!!$       PRINT*, 'TOTAL NUMBER OF CLOUD LAYERS OF THE FOV:', NCLD
       CALL CALC_CLD_PROPERTY( CLD,                 &
                               NCLD,                &
                               ATM,                 &
                               PCRTM_STND,          &
                               ICE_GRID,            &
                               WAT_GRID,            &
                               GEOMETRY%SATANG,     &
                               K_M%JACOB,           &
                               K_CLD)
    END IF


    CALL PCRTM_CALCRADUP(GEOMETRY,                        &
                         PCRTM_STND,                      &
                         NCLD,                            &
                         CLD,                             &
                         K_M,                             &
                         K_CLD,                           &
                         RT_SOLUTION )
    
    
    CALL PCRTM_CALCRADDN(GEOMETRY,                        &
                         PCRTM_STND,                      &
                         NCLD,                            &
                         CLD,                             &
                         K_M,                             &
                         K_CLD,                           &
                         RT_SOLUTION)
    
    RT_SOLUTION%RADUP = RT_SOLUTION%RADUP+          &
                       (1.0-RT_SOLUTION%EMIS)*      &
                        RT_SOLUTION%RADDN+          &
                        RT_SOLUTION%RADDNCLD
    

    IF (GEOMETRY%SOLAR_ZANG .lT. 90) THEN
       CALL PCRTM_FORWARD_RT_M_Solar ( RT_SOLUTION,     &
                                       PCRTM_STND,      &
                                       NCLD,            &
                                       CLD,             &
                                       GEOMETRY,        &
                                       SOLAR_TAB,       &
                                       SOLAR_SOLUTION)

       rt_solution%radup(SOLAR_SOLUTION%StartWaveIndex:PCRTM_STND%nM) = SOLAR_SOLUTION%SolarRadUp + &
                                        rt_solution%radup(SOLAR_SOLUTION%StartWaveIndex:PCRTM_STND%nM)  

    END IF



    IF (K_M%JACOB) THEN

       K_M%DRUPDEM  = K_M%DRUPDEM - RT_SOLUTION%RADDN

       DO I1 = GEOMETRY%NTOP,GEOMETRY%NBOT
          K_M%DRUPDT(:,I1) = K_M%DRUPDT(:,I1) + (1.0-RT_SOLUTION%EMIS)*K_M%DRDNDT(:,I1) + K_M%DRDNCLDDT(:,I1)
          K_M%DRUPDH2O(:,I1) = K_M%DRUPDH2O(:,I1) + (1.0-RT_SOLUTION%EMIS)*K_M%DRDNDH2O(:,I1) + K_M%DRDNCLDDH2O(:,I1)
          DO IM = 1, PCRTM_STND%NM
             DO IGG=1,PCRTM_STND%NGAS(IM)
                IG=PCRTM_STND%IDGAS(IGG,IM)
                K_M%DRUPDGAS(IG,IM,I1) = K_M%DRUPDGAS(IG,IM,I1) + (1.0-RT_SOLUTION%EMIS(IM))*K_M%DRDNDGAS(IG,IM,I1) &
                                         + K_M%DRDNCLDDGAS(IG,IM,I1) 
             END DO
          END DO
       ENDDO

       IF (NCLD .GT. 0) THEN
          DO K = 1, NCLD
             K_CLD(K)%DRDCLDRF  = K_CLD(K)%DRUPDCLDRF + (1-RT_SOLUTION%EMIS)*K_CLD(K)%DRDNDCLDRF + K_CLD(K)%DRDNCLDDCLDRF 
             K_CLD(K)%DRDCLDTR  = K_CLD(K)%DRUPDCLDTR + (1-RT_SOLUTION%EMIS)*K_CLD(K)%DRDNDCLDTR + K_CLD(K)%DRDNCLDDCLDTR
             K_CLD(K)%DRDDECLD  = K_CLD(K)%DRDCLDRF*CLD(K)%DRFDDE  + K_CLD(K)%DRDCLDTR*CLD(K)%DTRDDE
             K_CLD(K)%DRDTAUCLD = K_CLD(K)%DRDCLDRF*CLD(K)%DRFDTAU + K_CLD(K)%DRDCLDTR*CLD(K)%DTRDTAU
             K_CLD(K)%DRDTCLD   = K_CLD(K)%DRUPDTCLD  + (1-RT_SOLUTION%EMIS)*K_CLD(K)%DRDNDTCLD + K_CLD(K)%DRDNCLDDTCLD
             K_CLD(K)%DRDPCLD   = K_CLD(K)%DRDTCLD*CLD(K)%DTCLDDPCLD
          END DO
       END IF
    END IF
    
  END SUBROUTINE PCRTM_FORWARD_RT_M


  SUBROUTINE PRED_EOF_CH_SOLUTION(ATM,              &
                                  GEOMETRY,         &
                                  RT_SOLUTION,      &
                                  EOF_SOLUTION,     &
                                  K_M,              &
                                  K_PC,             &
                                  K_CH )

    TYPE(PCRTM_RT_SOLUTION_TYPE),     INTENT(IN)    :: RT_SOLUTION
    TYPE(PCRTM_NM_JACOBIAN_TYPE),     INTENT(IN)    :: K_M
    TYPE(PCRTM_CH_JACOBIAN_TYPE), ALLOCATABLE, &
                                      INTENT(INOUT) :: K_CH(:)
    TYPE(PCRTM_PC_JACOBIAN_TYPE), ALLOCATABLE, &
                                      INTENT(INOUT) :: K_PC(:)    
    TYPE(PCRTM_EOF_SOLUTION_TYPE),&
                        ALLOCATABLE,  INTENT(INOUT) :: EOF_SOLUTION(:)    
    TYPE(PCRTM_ATMOSPHERE_TYPE),      INTENT(IN)    :: ATM
    TYPE(PCRTM_GEOMETRY_TYPE),        INTENT(IN)    :: GEOMETRY
    
    INTEGER                                         :: IB, NB

    NB = SIZE(EOF_SOLUTION)

    DO IB = 1,NB 
       CALL PRED2PCSCORE(RT_SOLUTION, EOF_SOLUTION(IB))
       IF(EOF_SOLUTION(IB)%CHDOMAIN) THEN
          CALL PCSCORE2CH(EOF_SOLUTION(IB))
       END IF
    END DO

    IF(ALLOCATED(K_PC)) THEN
       DO IB = 1, NB
          IF(K_PC(IB)%JACOB_PC) THEN
             CALL PCRTM_JACOBIAN_PC_SOLUTION(K_M,                &
                                             ATM%DLAYDLEV,       &
                                             ATM%AIRAMT,         & 
                                             GEOMETRY%NTOP,      &
                                             GEOMETRY%NBOT,      &      
                                             EOF_SOLUTION(IB),   &
                                             K_PC(IB))
          END IF
       END DO
    END IF

    IF(ALLOCATED(K_CH)) THEN
       DO IB = 1, NB
          IF(K_CH(IB)%JACOB_CH) THEN
             CALL PCRTM_JACOBIAN_CH_SOLUTION(GEOMETRY%NTOP,   &
                                             GEOMETRY%NBOT,   &      
                                             EOF_SOLUTION(IB),&
                                             K_PC(IB),        &  
                                             K_CH(IB))
          END IF
       END DO
    END IF

    
  END SUBROUTINE PRED_EOF_CH_SOLUTION

  
  SUBROUTINE  PCRTM_FORWARD_RT( GEOMETRY,       &
                                ATM,            &
                                ATM_ABS,        &
                                PCRTM_STND,     &
                                ICE_GRID,       &
                                WAT_GRID,       &
                                RT_SOLUTION,    &
                                SOLAR_TAB,      &
                                SOLAR_SOLUTION, &
                                EOF_SOLUTION,   &
                                K_M,            &
                                K_PC,           &
                                K_CH,           & 
                                TR_SOLUTION )                          
                          

    TYPE(PCRTM_GEOMETRY_TYPE),        INTENT(INOUT) :: GEOMETRY
    TYPE(PCRTM_ATM_ABS_STRUCT_TYPE),  INTENT(IN)    :: PCRTM_STND
    TYPE(PCRTM_ATMOSPHERE_TYPE),      INTENT(INOUT) :: ATM
    TYPE(PCRTM_ATM_ABSORPTION_TYPE),  INTENT(IN)    :: ATM_ABS
    TYPE(PCRTM_NM_JACOBIAN_TYPE),     INTENT(INOUT) :: K_M
   
    TYPE(PCRTM_RT_SOLUTION_TYPE),     INTENT(INOUT) :: RT_SOLUTION
    type(PCRTM_SOLAR_LUT_DEF),        intent(in)    :: solar_tab(2)
    TYPE(PCRTM_SOLAR_DEF),            INTENT(inout) :: SOLAR_SOLUTION

    TYPE(PCRTM_CLD_TABLE_DEF),        INTENT(IN)    :: ICE_GRID
    TYPE(PCRTM_CLD_TABLE_DEF),        INTENT(IN)    :: WAT_GRID
    TYPE(PCRTM_CH_JACOBIAN_TYPE), ALLOCATABLE, &
                                      INTENT(INOUT) :: K_CH(:)
    TYPE(PCRTM_PC_JACOBIAN_TYPE), ALLOCATABLE, &
                                      INTENT(INOUT) :: K_PC(:)    
    TYPE(PCRTM_EOF_SOLUTION_TYPE),&
                        ALLOCATABLE,  INTENT(INOUT) :: EOF_SOLUTION(:)  
    type(PCRTM_TR_solution_type),&
                        ALLOCATABLE,  INTENT(INOUT) :: TR_SOLUTION(:)

      
    TYPE(PCRTM_CLD_JACOBIAN_TYPE),    ALLOCATABLE   :: K_CLD(:)
    TYPE(PCRTM_CLOUD_TYPE),           ALLOCATABLE   :: CLD(:)
    INTEGER                                         :: IB,IM,NCLD
    INTEGER                                         :: NB,I,ILEV
    INTEGER                                         :: K1, K2
    INTEGER                                         :: ALLOC_STAT

    REAL(SINGLE),                     ALLOCATABLE   :: TMP_M1(:), TMP_M2(:)
    REAL(SINGLE),                     ALLOCATABLE   :: TMP_PC(:), TMP_CH(:) 

    NCLD = 0
    DO IM = 1, PCRTM_STND%NLAY
       IF (ATM%CLD_FLAG(IM) .GT. 0) THEN
          NCLD = NCLD + 1
       END IF
    END DO

    NB = SIZE(EOF_SOLUTION)

    IF (NCLD .GT. 0) THEN
       ALLOCATE(CLD(NCLD),K_CLD(NCLD),STAT = ALLOC_STAT)
       IF ( ALLOC_STAT /= 0 ) THEN
          PRINT*,'ERROR TRYING TO ALLOCATE CLOUD IN PCRTM_FORWARD_RT'
          STOP
       ENDIF
    END IF

    CALL PCRTM_FORWARD_RT_M( GEOMETRY,      &
                             ATM,           &
                             NCLD,          &
                             CLD,           &
                             ATM_ABS,       &
                             PCRTM_STND,    &
                             ICE_GRID,      &
                             WAT_GRID,      &
                             RT_SOLUTION,   &
                             SOLAR_TAB,     &
                             SOLAR_SOLUTION,&
                             K_M,           &
                             K_CLD )


    CALL PRED_EOF_CH_SOLUTION(ATM,              &
                              GEOMETRY,         &
                              RT_SOLUTION,      &
                              EOF_SOLUTION,     &
                              K_M,              &
                              K_PC,             &
                              K_CH )

    IF (NCLD .GT. 0 .AND. K_M%JACOB ) THEN
       ALLOCATE(TMP_M1(PCRTM_STND%NM),TMP_M2(PCRTM_STND%NM),STAT = ALLOC_STAT)
       IF ( ALLOC_STAT /= 0 ) THEN
          PRINT*,'ERROR TRYING TO ALLOCATE TMP_M IN PCRTM_FORWARD_RT'
          STOP
       ENDIF
       DO I = 1, NCLD
          ILEV = CLD(I)%ILAY
          TMP_M1 = K_CLD(I)%DRDTCLD*CLD(I)%DTCLDDTLEVU 
          TMP_M2 = K_CLD(I)%DRDTCLD*CLD(I)%DTCLDDTLEVD 
          DO IB = 1, NB
             ALLOCATE(TMP_PC(EOF_SOLUTION(IB)%NPCBND),STAT = ALLOC_STAT)
             IF ( ALLOC_STAT /= 0 ) THEN
                PRINT*,'ERROR TRYING TO ALLOCATE TMP_PC IN PCRTM_FORWARD_RT'
                STOP
             ENDIF
             CALL PRED2PC_JACOB(TMP_M1, EOF_SOLUTION(IB), TMP_PC)
             K_PC(IB)%K_T(:,ILEV) = K_PC(IB)%K_T(:,ILEV) + TMP_PC
!!$              DO IM = 1, EOF_SOLUTION(IB)%NPCBND
!!$                write(104,*)K_PC(IB)%K_TLAY(IM,ILEV-1:ILEV),K_PC(IB)%K_T(IM,ILEV),  TMP_PC(IM)
!!$             END DO
             CALL PRED2PC_JACOB(TMP_M2, EOF_SOLUTION(IB), TMP_PC)
             K_PC(IB)%K_T(:,ILEV+1) = K_PC(IB)%K_T(:,ILEV+1) + TMP_PC
!!$             DO IM = 1, EOF_SOLUTION(IB)%NPCBND
!!$                write(105,*)K_PC(IB)%K_TLAY(IM,ILEV:ILEV+1),K_PC(IB)%K_T(IM,ILEV+1),TMP_PC(IM)
!!$             END DO
             DEALLOCATE(TMP_PC,STAT = ALLOC_STAT)
             IF ( ALLOC_STAT /= 0 ) THEN
                PRINT*,'ERROR TRYING TO DEALLOCATE TMP_PC IN PCRTM_FORWARD_RT'
                STOP
             ENDIF
          END DO
       END DO
       DEALLOCATE(TMP_M1,TMP_M2,STAT = ALLOC_STAT)
       IF ( ALLOC_STAT /= 0 ) THEN
          PRINT*,'ERROR TRYING TO DEALLOCATE TMP_M IN PCRTM_FORWARD_RT'
          STOP
       ENDIF
    END IF

    IF(ALLOCATED(K_CH)) THEN
       DO IM = 1, NCLD
          ILEV = CLD(IM)%ILAY
          DO IB = 1, NB
             IF(K_CH(IB)%JACOB_CH) THEN
                ALLOCATE(TMP_CH(EOF_SOLUTION(ib)%NCHBND))
                DO I=1,EOF_SOLUTION(IB)%NCHBND
                   K_CH(IB)%R_TLEV(I,ILEV)   = DOT_PRODUCT( EOF_SOLUTION(IB)%PC(:,I), K_PC(IB)%K_T(:,ILEV) )
                   K_CH(IB)%R_TLEV(I,ILEV+1) = DOT_PRODUCT( EOF_SOLUTION(IB)%PC(:,I), K_PC(IB)%K_T(:,ILEV+1) )
                END DO
                K_CH(IB)%R_TLEV(:,ILEV)      =  K_CH(IB)%R_TLEV(:,ILEV)*EOF_SOLUTION(IB)%RADSTDCH
                K_CH(IB)%R_TLEV(:,ILEV+1)    =  K_CH(IB)%R_TLEV(:,ILEV+1)*EOF_SOLUTION(IB)%RADSTDCH
                IF(K_CH(IB)%JACOB_BTCH) THEN
                   TMP_CH = dBTdRad(EOF_SOLUTION(ib)%frqch,EOF_SOLUTION(ib)%Radch,EOF_SOLUTION(IB)%Nchbnd)
                   K_CH(IB)%BT_TLEV(:,ILEV)     = TMP_CH*K_CH(IB)%R_TLEV(:,ILEV)
                   K_CH(IB)%BT_TLEV(:,ILEV+1)   = TMP_CH*K_CH(IB)%R_TLEV(:,ILEV+1)
                END IF
                DEALLOCATE(TMP_CH)
             END IF
          END DO
       END DO
    END IF


    IF(ALLOCATED(TR_SOLUTION)) THEN
       DO IB = 1, SIZE(EOF_SOLUTION)
          IF(TR_SOLUTION(IB)%FLAG) THEN
             CALL PCRTM_TR_SOLUTION_SOLVER( TR_SOLUTION(IB),    &
                                            RT_SOLUTION,        &
                                            EOF_SOLUTION(IB),   &
                                            PCRTM_STND,         &
                                            GEOMETRY%NTOP,      &
                                            GEOMETRY%NBOT )
          END IF
       END DO
    END IF

    IF (NCLD .GT.0 ) THEN
       CALL CLEAR_PCRTM_CLOUD_PARAM( CLD, NCLD, K_M%JACOB, K_CLD)
       DEALLOCATE(CLD,K_CLD,STAT = ALLOC_STAT)
       IF ( ALLOC_STAT /= 0 ) THEN
          PRINT*,'ERROR TRYING TO DEALLOCATE CLOUD IN PCRTM_FORWARD_RT'
          STOP
       ENDIF
    END IF


  END SUBROUTINE PCRTM_FORWARD_RT

  
END MODULE PCRTM_FORWARD_MODEL
