MODULE PCRTM_CALC_RAD

  USE PCRTM_TYPE_KIND
  USE PCRTM_MATH_UTILITY
  USE PCRTM_ATM_ABSORPTION_DEFINE
  USE PCRTM_JACOBIAN
  USE PCRTM_CLOUD_DEFINE
  USE PCRTM_RT_SOLUTION_DEFINE, ONLY : PCRTM_RT_SOLUTION_TYPE
  USE PCRTM_ATMOSPHERE_LAYER, ONLY   : PCRTM_GEOMETRY_TYPE

CONTAINS
!##############################################################################
  SUBROUTINE PCRTM_CALCRADUP(GEOMETRY,       &
                             PCRTM_STND,     &
                             NCLD,           &
                             CLD,            &
                             K_NM,           &
                             K_CLD,          &
                             RT_SOLUTION )

    TYPE(PCRTM_ATM_ABS_STRUCT_TYPE),INTENT(IN)  :: PCRTM_STND
    INTEGER,                        INTENT(IN)  :: NCLD
    TYPE(PCRTM_CLOUD_TYPE),       INTENT(INOUT) :: CLD(NCLD)
    TYPE(PCRTM_NM_JACOBIAN_TYPE), INTENT(INOUT) :: K_NM
    TYPE(PCRTM_CLD_JACOBIAN_TYPE),INTENT(INOUT) :: K_CLD(NCLD)
    TYPE(PCRTM_RT_SOLUTION_TYPE), INTENT(INOUT) :: RT_SOLUTION
    TYPE(PCRTM_GEOMETRY_TYPE),    INTENT(IN)    :: GEOMETRY

    REAL(SINGLE) :: DRDTAU(PCRTM_STND%NM)
    REAL(SINGLE) :: TRAN(PCRTM_STND%NM)
    REAL(SINGLE) :: TR(PCRTM_STND%NM)
    REAL(SINGLE) :: TRT2L(PCRTM_STND%NM)
    INTEGER      :: K,I1,IG,IGG,IM
    INTEGER      :: NTOP, NBOT
    REAL(SINGLE) :: SECZANG

    NTOP    = GEOMETRY%NTOP
    NBOT    = GEOMETRY%NBOT
    SECZANG = GEOMETRY%SECZANG
      
    RT_SOLUTION%RADUP = RT_SOLUTION%EMIS*RT_SOLUTION%BS

    TRAN =1
    IF (NCLD .GT.0) THEN
       DO K = 1,NCLD
          TRAN = TRAN*CLD(K)%TRANS
       END DO
    END IF
    IF(K_NM%JACOB)THEN
       K_NM%DRUPDTS = RT_SOLUTION%EMIS * RT_SOLUTION%TRTOP2LV(:,NBOT+1)*TRAN             &
            *RT_SOLUTION%DBSDTS
       K_NM%DRUPDEM = RT_SOLUTION%TRTOP2LV(:,NBOT+1)*TRAN*RT_SOLUTION%BS
    END IF

    DO I1=NBOT,NTOP,-1

       TR=RT_SOLUTION%TRLAY(:,I1)
       DO K = 1,NCLD
          IF ( I1 .EQ. CLD(K)%ILEV-1 ) THEN
             TRAN =TRAN/(CLD(K)%TRANS+1E-20)
          END IF
       END DO
       TRT2L= RT_SOLUTION%TRTOP2LV(:,I1)*TRAN
       
       IF(K_NM%JACOB)THEN 
          DRDTAU = (RT_SOLUTION%BLAY(:,I1)- RT_SOLUTION%RADUP)*SECZANG                      &
               *TR*TRT2L 
          K_NM%DRUPDT(:,I1)  = DRDTAU*K_NM%DTAUDT(:,I1)+RT_SOLUTION%DBDT(:,I1)*             & 
               (1.0-TR)*TRT2L
          K_NM%DRUPDH2O(:,I1)= DRDTAU*K_NM%DTAUDH2O(:,I1)
          DO IM = 1, PCRTM_STND%NM
             DO IGG=1,PCRTM_STND%NGAS(IM)
                IG =PCRTM_STND%IDGAS(IGG,IM)
                K_NM%DRUPDGAS(IG,IM,I1) = DRDTAU(IM)*K_NM%DTAUDGAS(IG,IM,I1)
             ENDDO
          END DO
       ENDIF
          
       RT_SOLUTION%RADUP = RT_SOLUTION%RADUP*TR+(1.0-TR)*RT_SOLUTION%BLAY(:,I1) 
       
       IF (NCLD .GT.0) THEN
          DO K = 1,NCLD
             IF( I1 .EQ. CLD(K)%ILEV ) THEN
                IF (K_NM%JACOB) THEN
                   K_CLD(K)%DRUPDCLDTR = (RT_SOLUTION%RADUP-CLD(K)%BCLD)*TRT2L/(CLD(K)%TRANS+1E-20)
                   K_CLD(K)%DRUPDCLDRF = -CLD(K)%BCLD*TRT2L/(CLD(K)%TRANS+1E-20)
                   K_CLD(K)%DRUPDTCLD  = (1.0-CLD(K)%TRANS -CLD(K)%REFL)*TRT2L*             &
                                          CLD(K)%DBCLDDTCLD/(CLD(K)%TRANS+1E-20)
                END IF
                RT_SOLUTION%RADUP      = RT_SOLUTION%RADUP*CLD(K)%TRANS+CLD(K)%BCLD*        & 
                                         (1.0-CLD(K)%TRANS-CLD(K)%REFL)
             ENDIF
          END DO
       END IF

    ENDDO !ENDDO I1
   
    RETURN
  END SUBROUTINE PCRTM_CALCRADUP

!#######################################################################################################################

  SUBROUTINE PCRTM_CALCRADDN(GEOMETRY,        &
                             PCRTM_STND,      &
                             NCLD,            &
                             CLD,             &
                             K_NM,            &
                             K_CLD,           &
                             RT_SOLUTION)

    TYPE(PCRTM_ATM_ABS_STRUCT_TYPE),INTENT(IN)   :: PCRTM_STND
    INTEGER,                        INTENT(IN)   :: NCLD
    TYPE(PCRTM_CLOUD_TYPE),         INTENT(IN)   :: CLD(NCLD)
    TYPE(PCRTM_NM_JACOBIAN_TYPE),   INTENT(INOUT):: K_NM
    TYPE(PCRTM_CLD_JACOBIAN_TYPE),  INTENT(INOUT):: K_CLD(NCLD)
    TYPE(PCRTM_RT_SOLUTION_TYPE),   INTENT(INOUT):: RT_SOLUTION
    TYPE(PCRTM_GEOMETRY_TYPE),      INTENT(IN)   :: GEOMETRY


    REAL(SINGLE) :: DRDTAU(PCRTM_STND%NM)
    REAL(SINGLE) :: TRAN(PCRTM_STND%NM)
    REAL(SINGLE) :: TRC(PCRTM_STND%NM)
    REAL(SINGLE) :: TRSFC(PCRTM_STND%NM)
    REAL(SINGLE) :: TR(PCRTM_STND%NM)
    REAL(SINGLE) :: TRT2C2L(PCRTM_STND%NM)
    REAL(SINGLE) :: TRT2B2L(PCRTM_STND%NM)
    REAL(SINGLE) :: TRT2B2C(NCLD,PCRTM_STND%NM)
    INTEGER      :: K, I1, N, M, IM  
    REAL(SINGLE) :: DRADCLDRFDCLDRF(NCLD,PCRTM_STND%NM)
    REAL(SINGLE) :: DRADCLDRFDCLDTR(NCLD,NCLD,PCRTM_STND%NM)
    REAL(SINGLE) :: DRADCLDRFDTCLD(NCLD,NCLD,PCRTM_STND%NM)
    REAL(SINGLE) :: RADCLDRF(NCLD,PCRTM_STND%NM)
    INTEGER      :: NTOP, NBOT
    REAL(SINGLE) :: SECZANG
    NTOP    = 1
    NBOT    = GEOMETRY%NBOT
    SECZANG = GEOMETRY%SECZANG
    TRSFC   = RT_SOLUTION%TRTOP2LV(:,NBOT+1)

    RT_SOLUTION%RADDN = 0
    
    IF (NCLD .GT.0) THEN
       DO K = 1,NCLD
          TRSFC = TRSFC*CLD(K)%TRANS
       END DO
    END IF

    DO I1=NTOP,NBOT
       TRAN = 1
       DO K = 1,NCLD
          IF (I1 .GE. CLD(K)%ILEV  )THEN
             TRAN=TRAN*CLD(K)%TRANS
          END IF
       END DO

       TR      = RT_SOLUTION%TRLAY(:,I1)  
       TRT2B2L = TRSFC*TRSFC/(TRAN+1.0E-15)/(RT_SOLUTION%TRTOP2LV(:,I1)+1.0E-15)/(TR+1.0E-15)
       
       IF(K_NM%JACOB)THEN
          DRDTAU              = (RT_SOLUTION%BLAY(:,I1)-RT_SOLUTION%RADDN)*TR               &
                                 *SECZANG*TRT2B2L
          K_NM%DRDNDT(:,I1)   = (DRDTAU* K_NM%DTAUDT(:,I1)+                                 &
                                (1.0-TR)*RT_SOLUTION%DBDT(:,I1))*TRT2B2L
!!$          print*, I1, sum(K_NM%DRdnDT(:,I1)), isnan(sum(K_NM%DRdnDT(:,I1)))
!!$          print*,'****************************************'
!!$          do im = 1, 466
!!$             print*, im, K_NM%drdndt(im,I1),drdtau(im), TR(im),  TRT2B2L(im), TRSFC(im), RT_SOLUTION%TRTOP2LV(im,I1)
!!$             if (I1 .eq. 10 )stop
!!$          end do
!!$          print*,'****************************************'

          K_NM%DRDNDH2O(:,I1) = DRDTAU*K_NM%DTAUDH2O(:,I1)
          DO IM = 1, PCRTM_STND%NM
             DO IGG=1,PCRTM_STND%NGAS(IM)
                IG =PCRTM_STND%IDGAS(IGG,IM)
                K_NM%DRDNDGAS(IG,IM,I1)= DRDTAU(IM)*K_NM%DTAUDGAS(IG,IM,I1)
             ENDDO
          END DO
       ENDIF
       RT_SOLUTION%RADDN = RT_SOLUTION%RADDN*TR+(1.0-TR)*RT_SOLUTION%BLAY(:,I1)

       IF (NCLD .GT. 0) THEN
          DO K =1, NCLD
             IF(I1 .EQ. CLD(K)%ILEV-1) THEN
                TRT2B2C(K,:) = TRT2B2L/(TR+1.0E-20)
                IF (K_NM%JACOB) THEN
                   K_CLD(K)%DRDNDCLDTR= (RT_SOLUTION%RADDN-CLD(K)%BCLD)*TRT2B2C(K,:)
                   K_CLD(K)%DRDNDCLDRF= -CLD(K)%BCLD*TRT2B2C(K,:)     
                   K_CLD(K)%DRDNDTCLD = (1.0-CLD(K)%TRANS-CLD(K)%REFL)*TRT2B2C(K,:)*                 &
                        CLD(K)%DBCLDDTCLD
                   DRADCLDRFDCLDRF(K,:) = TRAN*RT_SOLUTION%TRTOP2LV(:,I1+1)*TR*RT_SOLUTION%RADDN
                END IF
                RADCLDRF(K,:) = RT_SOLUTION%RADDN*RT_SOLUTION%TRTOP2LV(:,I1+1)*TRAN*CLD(K)%REFL
                IF (K_NM%JACOB) THEN
                   IF (K .GT. 1) THEN
                      DO L = 1,K-1
                         TRC = 1
                         DO N = 1,K-L-1
                            TRC = TRC*CLD(L+N)%TRANS
                         END DO
                         TRC = TRAN/(TRT2B2C(L,:)+1E-20)*RT_SOLUTION%TRTOP2LV(:,CLD(K)%ILEV)*        &
                               CLD(K)%REFL/(TR+1.0E-20)*RT_SOLUTION%TRTOP2LV(:,CLD(K)%ILEV)          &
                               /(RT_SOLUTION%TRTOP2LV(:,CLD(L)%ILEV)+1E-20)*TRC +                    &
                               RT_SOLUTION%RADDN*RT_SOLUTION%TRTOP2LV(:,CLD(K)%ILEV)*                &
                               CLD(K)%REFL*TRAN/(CLD(L)%TRANS+1.0E-20)
                         DRADCLDRFDCLDTR(K,L,:) =   K_CLD(L)%DRDNDCLDTR*TRC
                         DRADCLDRFDTCLD(K,L,:)  =   K_CLD(L)%DRDNDTCLD*TRC  
                      END DO
                   END IF
                END IF
                RT_SOLUTION%RADDN  =  RT_SOLUTION%RADDN*CLD(K)%TRANS+                                &
                                       (1.0-CLD(K)%TRANS-CLD(K)%REFL)*CLD(K)%BCLD
             ENDIF
          END DO  ! DO K

          IF(K_NM%JACOB) THEN
             K = 0
             DO N = 1,NCLD
                IF ( I1 .LT. CLD(N)%ILEV ) THEN
                   K = N
                   EXIT
                END IF
             END DO
             IF (K .GT. 0) THEN
                TRT2C2L   = RT_SOLUTION%TRTOP2LV(:,CLD(K)%ILEV)**2*TRAN*CLD(K)%REFL         &
                     /(RT_SOLUTION%TRTOP2LV(:,I1+1)+1E-20)   
                IF ( K .LT. NCLD ) THEN
                   DO N = 1, NCLD-K
                      TRC = 1
                      DO M = K,K+N-1
                         TRC = TRC*(CLD(M)%TRANS)**2
                      END DO
                      TRT2C2L= TRT2C2L + RT_SOLUTION%TRTOP2LV(:,CLD(K+N)%ILEV)**2           &
                           /(RT_SOLUTION%TRTOP2LV(:,I1+1)+1E-20)*TRAN*CLD(K+N)%REFL*TRC   
                   END DO
                END IF
             ELSE
                TRT2C2L = 0
             END IF
             DRDTAU                 = (RT_SOLUTION%BLAY(:,I1)-RT_SOLUTION%RADDN(:))         &
                                       *TR*SECZANG*TRT2C2L 
             K_NM%DRDNCLDDT(:,I1)   = DRDTAU*K_NM%DTAUDT(:,I1)+(1.0-TR)*                    &
                                      RT_SOLUTION%DBDT(:,I1)*TRT2C2L
             K_NM%DRDNCLDDH2O(:,I1) = DRDTAU*K_NM%DTAUDH2O(:,I1)
             DO IM = 1, PCRTM_STND%NM
                DO IGG=1,PCRTM_STND%NGAS(IM)
                   IG=PCRTM_STND%IDGAS(IGG,IM)
                   K_NM%DRDNCLDDGAS(IG,IM,I1)=DRDTAU(IM)*K_NM%DTAUDGAS(IG,IM,I1)
                END DO
             END DO
          END IF
       END IF

    END DO ! DO I1
    
    RT_SOLUTION%RADDN    = RT_SOLUTION%RADDN*TRSFC
    RT_SOLUTION%RADDNCLD = 0

    IF (NCLD .GT. 0) THEN
       DO K = 1,NCLD
          RT_SOLUTION%RADDNCLD         = RT_SOLUTION%RADDNCLD+RADCLDRF(K,:)
          IF (K_NM%JACOB) THEN
             K_CLD(K)%DRDNCLDDCLDTR    = 0 
             K_CLD(K)%DRDNCLDDTCLD     = 0
             K_CLD(K)%DRDNCLDDCLDRF    = DRADCLDRFDCLDRF(K,:)
             IF(K .LE. NCLD-1)THEN 
                DO L = K+1,NCLD
                   K_CLD(K)%DRDNCLDDCLDTR=K_CLD(K)%DRDNCLDDCLDTR+DRADCLDRFDCLDTR(L,K,:)
                   K_CLD(K)%DRDNCLDDTCLD =K_CLD(K)%DRDNCLDDTCLD +DRADCLDRFDTCLD(L,K,:)
                END DO
             END IF
          END IF
       ENDDO
    END IF

    RETURN
  END SUBROUTINE PCRTM_CALCRADDN

!#######################################################################################################################



END MODULE PCRTM_CALC_RAD





