MODULE PCRTM_CLOUD_LUT_INTERP

  USE PCRTM_TYPE_KIND
  USE PCRTM_CLOUD_DEFINE, ONLY   : PCRTM_CLOUD_TYPE,CLOUD_TAB_INDEX_TYPE
  USE PCRTM_CLOUD_LUT_IO, ONLY   : PCRTM_CLD_TABLE_DEF
  USE PCRTM_MATH_UTILITY, ONLY   : FIND_RANDOM_INDEX, BLEND_103, BLEND_103_JACOB 
  USE PCRTM_CLOUD_LUT_IO, ONLY   : TAUMAX, DMAX, THETMAX

CONTAINS
!##############################################################################

 SUBROUTINE PCRTM_INTERP_IDX(THETA, TAB, CLOUD, CLD_TAB_INDEX)

    IMPLICIT NONE 
    REAL(SINGLE),               INTENT(IN):: THETA
    TYPE(PCRTM_CLOUD_TYPE),  INTENT(INOUT):: CLOUD
    TYPE(PCRTM_CLD_TABLE_DEF),  INTENT(IN):: TAB
    TYPE(CLOUD_TAB_INDEX_TYPE),INTENT(OUT):: CLD_TAB_INDEX

    INTEGER :: TAU_IDX, TAU_IDX_, D_IDX, D_IDX_, THET_IDX, THET_IDX_
    REAL(SINGLE) :: R_TAU, R_D, R_THET

    IF (CLOUD%VISTAU.GT. TAB%TAU_V(TAUMAX+1))THEN  
       PRINT*, 'INPUT TAU LARGER THAN LIMIT' 
       PRINT*, 'VISTAU, TAB%TAU_V(TAUMAX+1)',CLOUD%VISTAU, TAB%TAU_V(TAUMAX+1)
       TAU_IDX  = TAUMAX+1
       TAU_IDX_ = TAUMAX
       R_TAU    = ( CLOUD%VISTAU - TAB%TAU_V(TAU_IDX_) )/( TAB%TAU_V(TAU_IDX)- TAB%TAU_V(TAU_IDX_) ) 
    ELSE
       CALL FIND_RANDOM_INDEX(TAB%TAU_V, CLOUD%VISTAU, TAU_IDX_, TAU_IDX, R_TAU)
    END IF

       IF (CLOUD%DE.LT.TAB%D(1)) then
       PRINT*, 'INPUT DE SMALLER THAN LIMIT'
       PRINT*, 'DE, TAB%D(1)',CLOUD%DE, TAB%D(1)
       CLOUD%DE = TAB%D(1)
       endif
       
    IF (CLOUD%DE.GT.TAB%D(DMAX))THEN
       PRINT*, 'INPUT DE LARGER THAN LIMIT'
       PRINT*, 'DE, TAB%D(DMAX)',CLOUD%DE, TAB%D(DMAX)
       CLOUD%DE = TAB%D(DMAX)
       PRINT*, 'CHANGE DE TO TAB%D(DMAX)','DE', CLOUD%DE
       D_IDX = DMAX
       D_IDX_ = DMAX-1
       R_D  = ( CLOUD%DE - TAB%D(D_IDX_) )/( TAB%D(D_IDX)- TAB%D(D_IDX_) ) 
    ELSE
       CALL FIND_RANDOM_INDEX(TAB%D, CLOUD%DE, D_IDX_, D_IDX, R_D)
    END IF
    
    IF (THETA .GT. TAB%THET(THETMAX)) THEN
       PRINT*, 'INPUT THET LARGER THAN 90 dgree' 
       PRINT*, 'THETA,TAB%THET(THETMAX)',THETA,TAB%THET(THETMAX)
       THET_IDX = THETMAX
       THET_IDX_= THETMAX-1
       stop
    ELSE
       CALL FIND_RANDOM_INDEX(TAB%THET, THETA, THET_IDX_, THET_IDX, R_THET)
    END IF

    CLD_TAB_INDEX%TAU_IDX   = TAU_IDX
    CLD_TAB_INDEX%TAU_IDX_  = TAU_IDX_
    CLD_TAB_INDEX%D_IDX     = D_IDX
    CLD_TAB_INDEX%D_IDX_    = D_IDX_
    CLD_TAB_INDEX%THET_IDX  = THET_IDX
    CLD_TAB_INDEX%THET_IDX_ = THET_IDX_
    CLD_TAB_INDEX%R_THET    = R_THET
    CLD_TAB_INDEX%R_TAU     = R_TAU
    CLD_TAB_INDEX%R_D       = R_D

    
 !--------------------------------------------------------
 
    RETURN
  END SUBROUTINE PCRTM_INTERP_IDX




!##############################################################################

  SUBROUTINE PCRTM_INTERP_CLOUD_TAB(CLOUD,         &
                                    TAB,           &
                                    CLD_TAB_INDEX, &
                                    IM,            &
                                    JACOB)
    TYPE(PCRTM_CLD_TABLE_DEF), INTENT(IN) :: TAB
    TYPE(PCRTM_CLOUD_TYPE), INTENT(INOUT) :: CLOUD
    INTEGER,                   INTENT(IN) :: IM 
    LOGICAL,                   INTENT(IN) :: JACOB
    TYPE(CLOUD_TAB_INDEX_TYPE),INTENT(IN) :: CLD_TAB_INDEX

    INTEGER      :: TAU_IDX, TAU_IDX_
    INTEGER      :: D_IDX, D_IDX_ 
    INTEGER      :: THET_IDX, THET_IDX_ 
    REAL(SINGLE) :: R_THET, R_TAU, R_D 
    REAL(SINGLE) :: TMP,TMP1,TMP2

    TAU_IDX   = CLD_TAB_INDEX%TAU_IDX
    TAU_IDX_  = CLD_TAB_INDEX%TAU_IDX_
    D_IDX     = CLD_TAB_INDEX%D_IDX  
    D_IDX_    = CLD_TAB_INDEX%D_IDX_
    THET_IDX  = CLD_TAB_INDEX%THET_IDX
    THET_IDX_ = CLD_TAB_INDEX%THET_IDX_
    R_THET    = CLD_TAB_INDEX%R_THET
    R_TAU     = CLD_TAB_INDEX%R_TAU
    R_D       = CLD_TAB_INDEX%R_D


    CALL BLEND_103 ( R_D, R_TAU, R_THET,                        &  
                     TAB%R_TAB(D_IDX_,TAU_IDX_,THET_IDX_,IM),   &
                     TAB%R_TAB(D_IDX_,TAU_IDX_,THET_IDX,IM),    &
                     TAB%R_TAB(D_IDX_,TAU_IDX,THET_IDX_,IM),    &
                     TAB%R_TAB(D_IDX_,TAU_IDX,THET_IDX,IM),     &
                     TAB%R_TAB(D_IDX,TAU_IDX_,THET_IDX_,IM),    &
                     TAB%R_TAB(D_IDX,TAU_IDX_,THET_IDX,IM),     &
                     TAB%R_TAB(D_IDX,TAU_IDX,THET_IDX_,IM),     &  
                     TAB%R_TAB(D_IDX,TAU_IDX,THET_IDX,IM),      & 
                     TMP )
    CLOUD%REFL(IM)  = TMP

    CALL BLEND_103 ( R_D, R_TAU, R_THET,                        &  
                     TAB%T_TAB(D_IDX_,TAU_IDX_,THET_IDX_,IM),   &
                     TAB%T_TAB(D_IDX_,TAU_IDX_,THET_IDX,IM),    &
                     TAB%T_TAB(D_IDX_,TAU_IDX,THET_IDX_,IM),    &
                     TAB%T_TAB(D_IDX_,TAU_IDX,THET_IDX,IM),     &
                     TAB%T_TAB(D_IDX,TAU_IDX_,THET_IDX_,IM),    &
                     TAB%T_TAB(D_IDX,TAU_IDX_,THET_IDX,IM),     &
                     TAB%T_TAB(D_IDX,TAU_IDX,THET_IDX_,IM),     &
                     TAB%T_TAB(D_IDX,TAU_IDX,THET_IDX,IM),      &
                     TMP )
    CLOUD%TRANS(IM) = TMP


    IF (JACOB) THEN
       CALL BLEND_103_JACOB ( R_D, R_TAU, R_THET,               &  
                     TAB%R_TAB(D_IDX_,TAU_IDX_,THET_IDX_,IM),   &
                     TAB%R_TAB(D_IDX_,TAU_IDX_,THET_IDX,IM),    &
                     TAB%R_TAB(D_IDX_,TAU_IDX,THET_IDX_,IM),    &
                     TAB%R_TAB(D_IDX_,TAU_IDX,THET_IDX,IM),     &
                     TAB%R_TAB(D_IDX,TAU_IDX_,THET_IDX_,IM),    &
                     TAB%R_TAB(D_IDX,TAU_IDX_,THET_IDX,IM),     &
                     TAB%R_TAB(D_IDX,TAU_IDX,THET_IDX_,IM),     &
                     TAB%R_TAB(D_IDX,TAU_IDX,THET_IDX,IM),      &
                     TMP1, TMP2, TMP )
       CLOUD%DRFDTAU(IM) = TMP2/( TAB%TAU_V(TAU_IDX) - TAB%TAU_V(TAU_IDX_) )
       CLOUD%DRFDDE(IM)  = TMP1/( TAB%D(D_IDX) - TAB%D(D_IDX_) )

       CALL BLEND_103_JACOB ( R_D, R_TAU, R_THET,               &  
                     TAB%T_TAB(D_IDX_,TAU_IDX_,THET_IDX_,IM),   &
                     TAB%T_TAB(D_IDX_,TAU_IDX_,THET_IDX,IM),    &
                     TAB%T_TAB(D_IDX_,TAU_IDX,THET_IDX_,IM),    &
                     TAB%T_TAB(D_IDX_,TAU_IDX,THET_IDX,IM),     &
                     TAB%T_TAB(D_IDX,TAU_IDX_,THET_IDX_,IM),    &
                     TAB%T_TAB(D_IDX,TAU_IDX_,THET_IDX,IM),     &
                     TAB%T_TAB(D_IDX,TAU_IDX,THET_IDX_,IM),     & 
                     TAB%T_TAB(D_IDX,TAU_IDX,THET_IDX,IM),      &
                     TMP1, TMP2, TMP )
       CLOUD%DTRDTAU(IM) = TMP2/( TAB%TAU_V(TAU_IDX) - TAB%TAU_V(TAU_IDX_) )
       CLOUD%DTRDDE(IM)  = TMP1/( TAB%D(D_IDX) - TAB%D(D_IDX_) )

    END IF

  END SUBROUTINE PCRTM_INTERP_CLOUD_TAB
!##############################################################################


END MODULE PCRTM_CLOUD_LUT_INTERP
