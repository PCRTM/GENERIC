MODULE CLEAR_PCRTM
  USE PCRTM_ATM_ABSORPTION_DEFINE
  USE PCRTM_CLOUD_LUT_IO
  USE PCRTM_ATMOSPHERE_DEFINE
  USE PCRTM_RT_SOLUTION_DEFINE
  USE PCRTM_SOLAR_DEFINE  
  USE PCRTM_SOLAR_PARAMETER,ONLY: PCRTM_SOLAR_LUT_DEF,CLEAR_SOLAR_BRDF_TAB
  USE PCRTM_JACOBIAN
  USE PCRTM_PC_SOLUTION, ONLY : PCRTM_EOF_SOLUTION_TYPE, &
                                CLEAR_PCRTM_EOF_SOLUTION
  USE PCRTM_TR_SOLUTION, ONLY : PCRTM_TR_SOLUTION_TYPE, &
                                    CLEAR_PCRTM_TR_SOLUTION


CONTAINS

  SUBROUTINE CLEAR_PCRTM_LUT( PCRTM_STND,     &
                              ICE_GRID,     &
                              WAT_GRID,     &
                              ATM_ABS_COEF, &
                              EOF_SOLUTION )


    TYPE(PCRTM_ATM_ABS_STRUCT_TYPE), &
                    INTENT (INOUT):: PCRTM_STND
    TYPE(PCRTM_CLD_TABLE_DEF),INTENT (INOUT):: ICE_GRID
    TYPE(PCRTM_CLD_TABLE_DEF),INTENT (INOUT):: WAT_GRID
    TYPE(PCRTM_ATM_ABSORPTION_TYPE), &
                    INTENT (INOUT):: ATM_ABS_COEF
    TYPE(PCRTM_EOF_SOLUTION_TYPE),   &
        ALLOCATABLE,INTENT (INOUT):: EOF_SOLUTION(:)
    INTEGER                       :: DEALLOC_STAT
    INTEGER                       :: N


    CALL CLEAR_PCRTM_ATM_ABSORPTION( ATM_ABS_COEF )
    
    CALL CLEAR_PCRTM_ATM_ABS_STRUCT( PCRTM_STND ) 
    
    CALL CLEAR_PCRTM_CLD_GRID(ICE_GRID)
      
    CALL CLEAR_PCRTM_CLD_GRID(WAT_GRID)

    DO N = 1,SIZE(EOF_SOLUTION)
       CALL CLEAR_PCRTM_EOF_SOLUTION(EOF_SOLUTION(N))
    END DO

    DEALLOCATE( EOF_SOLUTION, STAT = DEALLOC_STAT )
    IF ( DEALLOC_STAT /= 0 ) THEN
       PRINT*,'ERQROR TRYING TO DEALLOCATE EOF_SOLUTION'
       STOP
    ENDIF

  END SUBROUTINE CLEAR_PCRTM_LUT

  SUBROUTINE CLEAR_PCRTM_PARAM( ATM,RT_SOLUTION)

    TYPE(PCRTM_ATMOSPHERE_TYPE),     INTENT(INOUT)  :: ATM
    TYPE(PCRTM_RT_SOLUTION_TYPE),    INTENT(INOUT)  :: RT_SOLUTION


    CALL CLEAR_PCRTM_ATMOSPHERE( ATM )
     
    CALL CLEAR_PCRTM_RT_SOLUTION( RT_SOLUTION ) 
       
  END SUBROUTINE CLEAR_PCRTM_PARAM

  SUBROUTINE CLEAR_PCRTM_LP( PCRTM_STND,     &
                             ICE_GRID,       &
                             WAT_GRID,       &
                             ATM_ABS_COEF,   &
                             ATM,            & 
                             RT_SOLUTION,    &
                             SOLAR_TAB,      &
                             SOLAR_SOLUTION, &
                             EOF_SOLUTION,   &
                             TR_SOLUTION)

    TYPE(PCRTM_ATM_ABS_STRUCT_TYPE), &
                              INTENT (INOUT):: PCRTM_STND
    TYPE(PCRTM_CLD_TABLE_DEF),INTENT (INOUT):: ICE_GRID
    TYPE(PCRTM_CLD_TABLE_DEF),INTENT (INOUT):: WAT_GRID
    TYPE(PCRTM_ATM_ABSORPTION_TYPE), &
                              INTENT (INOUT):: ATM_ABS_COEF 

    TYPE(PCRTM_EOF_SOLUTION_TYPE),   &
                  ALLOCATABLE,INTENT (INOUT):: EOF_SOLUTION(:)
    TYPE(PCRTM_ATMOSPHERE_TYPE), &
                              INTENT(INOUT) :: ATM
    TYPE(PCRTM_RT_SOLUTION_TYPE),      &
                              INTENT(INOUT) :: RT_SOLUTION
    Type(PCRTM_SOLAR_DEF),    INTENT(inout) :: SOLAR_SOLUTION
    type(PCRTM_SOLAR_LUT_DEF), &
                              intent(inout) :: solar_tab(2)

    TYPE(PCRTM_TR_SOLUTION_TYPE),&
                  ALLOCATABLE,INTENT(INOUT) :: TR_SOLUTION(:)

    INTEGER                                 :: IB
    INTEGER                                 :: DEALLOC_STAT
    
    CALL CLEAR_PCRTM_PARAM( ATM,RT_SOLUTION)
    
    CALL CLEAR_PCRTM_SOLAR_SOLUTION(SOLAR_SOLUTION)

    CALL clear_solar_brdf_tab(solar_tab)

    CALL CLEAR_PCRTM_LUT( PCRTM_STND,   &
                          ICE_GRID,     &
                          WAT_GRID,     &
                          ATM_ABS_COEF, &
                          EOF_SOLUTION )

    IF(ALLOCATED(TR_SOLUTION)) THEN
       DO IB = 1, SIZE(TR_SOLUTION)
          IF(TR_SOLUTION(IB)%FLAG) THEN
             CALL  CLEAR_PCRTM_TR_SOLUTION(TR_SOLUTION(IB))
          END IF
       END DO
       DEALLOCATE(TR_SOLUTION,STAT = DEALLOC_STAT)
       IF(DEALLOC_STAT/=0) THEN
          PRINT*,'ERROR TRYING TO DEALLOCATE TR_SOLUTION'
          STOP
       ENDIF
    END IF
                              
  END SUBROUTINE CLEAR_PCRTM_LP


  SUBROUTINE CLEAR_PCRTM_JACOB( K_NM, K_PC, K_CH)
    TYPE(PCRTM_NM_JACOBIAN_TYPE),      &
                    INTENT(INOUT) :: K_NM 
    TYPE(PCRTM_PC_JACOBIAN_TYPE), ALLOCATABLE, &
                    INTENT(INOUT) :: K_PC(:)
    TYPE(PCRTM_CH_JACOBIAN_TYPE), ALLOCATABLE, &
                    INTENT(INOUT) :: K_CH(:)
                    
    INTEGER                       :: IB, DEALLOC_STAT

    CALL CLEAR_PCRTM_JACOB_NM( K_NM )

    IF(ALLOCATED(K_CH)) THEN 
       DO IB = 1, SIZE(K_CH)                
          CALL CLEAR_PCRTM_JACOB_CH( K_CH(IB) )
       END DO
       DEALLOCATE(K_CH,STAT = DEALLOC_STAT)
       IF ( DEALLOC_STAT /= 0 ) THEN
          PRINT*,'ERROR TRYING TO CLEAR K_CH'
          STOP
       ENDIF
    END IF

    IF(ALLOCATED(K_PC)) THEN 
       DO IB = 1, SIZE(K_PC)                
          CALL CLEAR_PCRTM_JACOB_PC( K_PC(IB) )
       END DO
       DEALLOCATE(K_PC,STAT = DEALLOC_STAT)
       IF ( DEALLOC_STAT /= 0 ) THEN
          PRINT*,'ERROR TRYING TO CLEAR K_PC'
          STOP
       ENDIF
    END IF


  END SUBROUTINE CLEAR_PCRTM_JACOB

END MODULE CLEAR_PCRTM
