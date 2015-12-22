MODULE INIT_PCRTM

  USE PCRTM_ATM_ABSORPTION_IO    
  USE PCRTM_CLOUD_LUT_IO
  USE PCRTM_ATMOSPHERE_DEFINE
  USE PCRTM_RT_SOLUTION_DEFINE
  USE PCRTM_SOLAR_DEFINE
  USE PCRTM_SOLAR_PARAMETER,ONLY: PCRTM_SOLAR_LUT_DEF, INIT_SOLAR_BRDF_TAB
  USE PCRTM_JACOBIAN
  USE PCRTM_TR_SOLUTION, ONLY : Init_PCRTM_TR_solution,PCRTM_TR_solution_type
  USE PCRTM_PC_SOLUTION, ONLY : PCRTM_EOF_SOLUTION_TYPE, RD_SENSOR_BND_INFO
  Use PCRTM_struct_Mod, only : PCRTM_struct_WN_Convert
  Use CLEAR_PCRTM, only : Clear_PCRTM_LUT

CONTAINS

    SUBROUTINE PCRTM_INIT( SENSOR_ID,        &
                           JACOB,            &
                           CHDOMAIN,         &
                           BT_FLAG,          &
                           JACOB_CH,         &
                           JACOB_BT,         &
                           TR_FLAG,          &
                           PCRTM_STND,       &
                           ATM_ABS_COEF,     &
                           ICE_GRID,         &
                           WAT_GRID,         &
                           EOF_SOLUTION,     &
                           ATM,              &
                           RT_SOLUTION,      &
                           SOLAR_TAB,        &
                           SOLAR_SOLUTION,   &
                           K_M,              &
                           K_PC,             &
                           K_CH,             & 
                           TR_SOLUTION,      &
                           FRQCH )                  


     INTEGER,              INTENT(IN)              :: SENSOR_ID
     LOGICAL,              INTENT(IN)              :: JACOB
     LOGICAL,              INTENT(IN)              :: CHDOMAIN
     LOGICAL,              INTENT(IN)              :: JACOB_CH
     LOGICAL,              INTENT(IN)              :: BT_FLAG
     LOGICAL,              INTENT(IN)              :: JACOB_BT
     LOGICAL,              INTENT(IN)              :: TR_FLAG
     REAL(DOUBLE),OPTIONAL,allocatable,INTENT(IN)  :: FRQCH(:)

     TYPE(PCRTM_ATM_ABS_STRUCT_TYPE),   INTENT(OUT):: PCRTM_STND
     TYPE(PCRTM_ATM_ABSORPTION_TYPE),   INTENT(OUT):: ATM_ABS_COEF
     TYPE(PCRTM_CLD_TABLE_DEF),         INTENT(OUT):: ICE_GRID
     TYPE(PCRTM_CLD_TABLE_DEF),         INTENT(OUT):: WAT_GRID
     TYPE(PCRTM_EOF_SOLUTION_TYPE),   &
                          ALLOCATABLE,  INTENT(OUT):: EOF_SOLUTION(:)
     TYPE(PCRTM_RT_SOLUTION_TYPE),      INTENT(OUT):: RT_SOLUTION
     type(PCRTM_SOLAR_LUT_DEF),         INTENT(OUT):: SOLAR_TAB(2)
     TYPE(PCRTM_SOLAR_DEF),             INTENT(OUT):: SOLAR_SOLUTION
     TYPE(PCRTM_NM_JACOBIAN_TYPE),      INTENT(OUT):: K_M
     TYPE(PCRTM_CH_JACOBIAN_TYPE),    &
                          ALLOCATABLE,  INTENT(OUT):: K_CH(:)
     TYPE(PCRTM_PC_JACOBIAN_TYPE),    &
                          ALLOCATABLE,  INTENT(OUT):: K_PC(:)

     TYPE(PCRTM_ATMOSPHERE_TYPE),       INTENT(OUT):: ATM
     type(PCRTM_TR_solution_type),&
                          ALLOCATABLE,  INTENT(OUT):: TR_SOLUTION(:)


     TYPE(PCRTM_ATM_ABS_STRUCT_TYPE)               :: PCRTM_STND_
     TYPE(PCRTM_ATM_ABSORPTION_TYPE)               :: ATM_ABS_COEF_
     TYPE(PCRTM_CLD_TABLE_DEF)                     :: ICE_GRID_
     TYPE(PCRTM_CLD_TABLE_DEF)                     :: WAT_GRID_
     TYPE(PCRTM_EOF_SOLUTION_TYPE),    ALLOCATABLE :: EOF_SOLUTION_(:)     
     INTEGER                                       :: NCH
     INTEGER                                       :: IB, NBND
     INTEGER                                       :: ALLOC_STAT

 
!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------     
     if(present(FRQCH))then
!!$   INITIATE ALLOCATABLE TYPE LUT VARIABLES
     
        CALL INIT_PCRTM_LUT( SENSOR_ID,                &
                             PCRTM_STND_,              &
                             ICE_GRID_,                &
                             WAT_GRID_,                &
                             ATM_ABS_COEF_,            & 
                             EOF_SOLUTION_)
     
       
!!$ USE THE USER DEFINED CHANNEL INFO TO REDEFINE THE PCRTM MODEL STRUCTURE
        NCH = size(FRQCH)
        CALL PCRTM_STRUCT_WN_CONVERT(FRQCH,            &
                                     NCH,              &  
                                     PCRTM_STND_,      &
                                     ATM_ABS_COEF_,    &
                                     ICE_GRID_,        &
                                     WAT_GRID_,        &
                                     EOF_SOLUTION_,    &
                                     PCRTM_STND,       &
                                     ATM_ABS_COEF,     &
                                     ICE_GRID,         &
                                     WAT_GRID,         &
                                     EOF_SOLUTION)
        
        CALL CLEAR_PCRTM_LUT( PCRTM_STND_,   &
                              ICE_GRID_,     &
                              WAT_GRID_,     &
                              ATM_ABS_COEF_, &
                              EOF_SOLUTION_ )

     else
        
        CALL INIT_PCRTM_LUT( SENSOR_ID,                &
                             PCRTM_STND,               &
                             ICE_GRID,                 &
                             WAT_GRID,                 &
                             ATM_ABS_COEF,             & 
                             EOF_SOLUTION)


     end if

     NBND = SIZE(EOF_SOLUTION)
     PRINT*, 'STRUCTURE OF EOF_SOLUTION'
     PRINT*,'----------#------MONO FREQ','------PC','--------CHANNEL---------'
     DO IB = 1, NBND
        PRINT*, IB, EOF_SOLUTION(IB)%NREG, EOF_SOLUTION(IB)%NPCBND, EOF_SOLUTION(IB)%NCHBND
        EOF_SOLUTION(IB)%chdomain = CHDOMAIN
        EOF_SOLUTION(IB)%BT_FLAG  = BT_FLAG
     END DO
     PRINT*,'----------------------------------------------------------'
!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------
!!$   INITIATE ALLOCATABLE TYPE PCRTM INPUT OUTPUT SOLUTION VARIABLES

     CALL INIT_PCRTM_PARAM( PCRTM_STND,             &
                            ATM,                    &
                            RT_SOLUTION )

     CALL INIT_PCRTM_SOLAR_SOLUTION( PCRTM_STND,    &
                                     SOLAR_SOLUTION)

     CALL INIT_SOLAR_BRDF_TAB(SOLAR_TAB)

     K_M%Jacob = .false.
     IF(JACOB) THEN
        CALL INIT_PCRTM_JACOB( PCRTM_STND,          &
                               EOF_SOLUTION,        &
                               K_M,                 &
                               K_PC,                &
                               K_CH,                &
                               JACOB_CH,            &
                               JACOB_BT )
     END IF

     if(TR_FLAG)then
        
        ALLOCATE(TR_SOLUTION(NBND),STAT = ALLOC_STAT)
        IF ( ALLOC_STAT /= 0 ) THEN
           PRINT*,'ERROR TRYING TO ALLOCATE TR_SOLUTION'
           STOP
        ENDIF

        do ib = 1,NBND
           call INIT_PCRTM_TR_solution( TR_SOLUTION(ib),           &
                                        EOF_SOLUTION(ib)%nchbnd,   &
                                        PCRTM_STND%nlay )
        end do
        
     end if

     
   END SUBROUTINE PCRTM_INIT 



  SUBROUTINE INIT_PCRTM_LUT( SENSOR_ID,    &
                             PCRTM_STND,   &
                             ICE_GRID,     &
                             WAT_GRID,     &
                             ATM_ABS_COEF, &
                             EOF_SOLUTION )


!!$  INPUT ARGUMENTS:
!!$     SENSOR NAME          SENSOR_ID  
!!$    ***********************************************        
!!$      'CLARREO_0.1'          1
!!$      'CLARREO_0.5'          2
!!$      'CLARREO_1.0'          3
!!$      'IASI-6341'            4
!!$      'IASI-8461'            5
!!$      'IASI-subset'          6    
!!$      'AIRS'                 7
!!$      'AIRS-subset'          8
!!$      'CRISBLK'              9
!!$      'CRISBOX'              10
!!$      'CRISHAM'              11
!!$      'NAST-I3'              12
!!$      'NAST-I44'             13

!!$    * absorption and regression coefficient files associated with senor id 1 -13 are obtained using lblrtm v-11.7;

!!$      'IASI'                    14    
!!$      'AIRS'                    15
!!$      'High-resolution CRISBLK' 16
!!$      'High-resolution CRISBOX' 17
!!$      'High-resolution CRISHAM' 18
!!$      'NAST-I BOX'              19
!!$      'NAST-I KSR'              20
!!$      'SHIS BOX'                21
!!$      'SHIS KSR'                22
!!$      'CRIS BOX'                23
!!$      'CRIS HAM'                24
!!$      'CRIS BLK'                25

!!$    * absorption and regression coefficient files associated with senor id 13-25 are obtained using lblrtm v-12.2;

     INTEGER,                  INTENT (IN)  :: SENSOR_ID
     TYPE(PCRTM_ATM_ABS_STRUCT_TYPE), &
                               INTENT (OUT) :: PCRTM_STND
     TYPE(PCRTM_CLD_TABLE_DEF),INTENT (OUT) :: ICE_GRID
     TYPE(PCRTM_CLD_TABLE_DEF),INTENT (OUT) :: WAT_GRID
     TYPE(PCRTM_ATM_ABSORPTION_TYPE), &
                               INTENT (OUT) :: ATM_ABS_COEF
     TYPE(PCRTM_EOF_SOLUTION_TYPE),   &
                  ALLOCATABLE, INTENT (OUT) :: EOF_SOLUTION(:)

     
     CHARACTER(160) :: ABSCOEFFILE(100)
     CHARACTER(160) :: ICE_TAB(100)
     CHARACTER(160) :: WAT_TAB(100)
     CHARACTER(160) :: PARFILE(100)

     ABSCOEFFILE(1) = 'INPUTDIR/ABS_COEF/clarreo_abscoef_0.1'
     ICE_TAB(1)     = 'INPUTDIR/CLD_LUT/ice_clarreo_0.1.dat'
     WAT_TAB(1)     = 'INPUTDIR/CLD_LUT/wat_clarreo_0.1.dat'
     PARFILE(1)     = 'INPUTDIR/EOF_COEF/clarreo_0.1_pccoef'
     
     ABSCOEFFILE(2) = 'INPUTDIR/ABS_COEF/clarreo_abscoef_0.5'
     ICE_TAB(2)     = 'INPUTDIR/CLD_LUT/ice_clarreo_0.5.dat'
     WAT_TAB(2)     = 'INPUTDIR/CLD_LUT/wat_clarreo_0.5.dat'
     PARFILE(2)     = 'INPUTDIR/EOF_COEF/clarreo_0.5_pccoef'

     ABSCOEFFILE(3) = 'INPUTDIR/ABS_COEF/clarreo_abscoef_1.0'
     ICE_TAB(3)     = 'INPUTDIR/CLD_LUT/ice_clarreo_1.0.dat'
     WAT_TAB(3)     = 'INPUTDIR/CLD_LUT/wat_clarreo_1.0.dat'
     PARFILE(3)     = 'INPUTDIR/EOF_COEF/clarreo_1.0_pccoef'

     ABSCOEFFILE(4) = 'INPUTDIR/ABS_COEF/iasi_abscoef_680'
     ICE_TAB(4)     = 'INPUTDIR/CLD_LUT/ice_iasi_rt680.dat'
     WAT_TAB(4)     = 'INPUTDIR/CLD_LUT/wat_iasi_rt680.dat'
     PARFILE(4)     = 'INPUTDIR/EOF_COEF/iasi_pccoef_400'

     ABSCOEFFILE(5) = 'INPUTDIR/ABS_COEF/iasi_abscoef_873'
     ICE_TAB(5)     = 'INPUTDIR/CLD_LUT/ice_iasi_rt873.dat'
     WAT_TAB(5)     = 'INPUTDIR/CLD_LUT/wat_iasi_rt873.dat'
     PARFILE(5)     = 'INPUTDIR/EOF_COEF/iasi_pccoef_1240'  

     ABSCOEFFILE(6) = 'INPUTDIR/ABS_COEF/iasi_subset_616_abscoef_407'
     ICE_TAB(6)     = 'INPUTDIR/CLD_LUT/ice_iasi_subset_rt_407.dat'
     WAT_TAB(6)     = 'INPUTDIR/CLD_LUT/wat_iasi_subset_rt_407.dat'
     PARFILE(6)     = 'INPUTDIR/EOF_COEF/iasi_subset_616_pccoef_500'       

     ABSCOEFFILE(7) = 'INPUTDIR/ABS_COEF/airs_abscoef_425'
     ICE_TAB(7)     = 'INPUTDIR/CLD_LUT/ice_airs_rt425.dat'
     WAT_TAB(7)     = 'INPUTDIR/CLD_LUT/wat_airs_rt425.dat'
     PARFILE(7)     = 'INPUTDIR/EOF_COEF/airs_pccoef'

     ABSCOEFFILE(8) = 'INPUTDIR/ABS_COEF/airs_subset_281_abscoef_238'
     ICE_TAB(8)     = 'INPUTDIR/CLD_LUT/ice_airs_subset_rt_238.dat'
     WAT_TAB(8)     = 'INPUTDIR/CLD_LUT/wat_airs_subset_rt_238.dat'
     PARFILE(8)     = 'INPUTDIR/EOF_COEF/airs_subset_281_pccoef_160'        

     ABSCOEFFILE(9) = 'INPUTDIR/ABS_COEF/cris_abscoef_360'
     ICE_TAB(9)     = 'INPUTDIR/CLD_LUT/ice_cris_rt360.dat'
     WAT_TAB(9)     = 'INPUTDIR/CLD_LUT/wat_cris_rt360.dat'
     PARFILE(9)     = 'INPUTDIR/EOF_COEF/crisblk_pccoef'     

     ABSCOEFFILE(10)= 'INPUTDIR/ABS_COEF/cris_abscoef_360'
     ICE_TAB(10)    = 'INPUTDIR/CLD_LUT/ice_cris_rt360.dat'
     WAT_TAB(10)    = 'INPUTDIR/CLD_LUT/wat_cris_rt360.dat'
     PARFILE(10)    = 'INPUTDIR/EOF_COEF/crisbox_pccoef' 

     ABSCOEFFILE(11)= 'INPUTDIR/ABS_COEF/cris_abscoef_360'
     ICE_TAB(11)    = 'INPUTDIR/CLD_LUT/ice_cris_rt360.dat'
     WAT_TAB(11)    = 'INPUTDIR/CLD_LUT/wat_cris_rt360.dat'
     PARFILE(11)    = 'INPUTDIR/EOF_COEF/crisham_pccoef' 

     ABSCOEFFILE(12)= 'INPUTDIR/ABS_COEF/nasti_abscoef_710'
     ICE_TAB(12)    = 'INPUTDIR/CLD_LUT/ice_nasti_rt710.dat'
     WAT_TAB(12)    = 'INPUTDIR/CLD_LUT/wat_nasti_rt710.dat'
     PARFILE(12)    = 'INPUTDIR/EOF_COEF/nasti_pccoef'
     
     ABSCOEFFILE(13)= 'INPUTDIR/ABS_COEF/nasti_abscoef_44_bnds'
     ICE_TAB(13)    = 'INPUTDIR/CLD_LUT/ice_nasti_rt2065.dat'
     WAT_TAB(13)    = 'INPUTDIR/CLD_LUT/wat_nasti_rt2065.dat'
     PARFILE(13)    = 'INPUTDIR/EOF_COEF/nasti_pccoef_44_bnds'

     ABSCOEFFILE(14)= 'INPUTDIR/ABS_COEF/iasi_abscoef_solar'
     ICE_TAB(14)    = 'INPUTDIR/CLD_LUT/ice_rt_753_iasi_solar.dat'
     WAT_TAB(14)    = 'INPUTDIR/CLD_LUT/wat_rt_753_iasi_solar.dat'
     PARFILE(14)    = 'INPUTDIR/EOF_COEF/iasi_pccoef_solar'     

     ABSCOEFFILE(15)= 'INPUTDIR/ABS_COEF/airs_abscoef_solar'
     ICE_TAB(15)    = 'INPUTDIR/CLD_LUT/ice_rt_500_airs_solar.dat'
     WAT_TAB(15)    = 'INPUTDIR/CLD_LUT/wat_rt_500_airs_solar.dat'
     PARFILE(15)    = 'INPUTDIR/EOF_COEF/airs_pccoef_solar'     
     
     ABSCOEFFILE(16)= 'INPUTDIR/ABS_COEF/crishrblk_abscoef_solar'
     ICE_TAB(16)    = 'INPUTDIR/CLD_LUT/ice_rt_374_crishr_blk_solar.dat'
     WAT_TAB(16)    = 'INPUTDIR/CLD_LUT/wat_rt_374_crishr_blk_solar.dat'
     PARFILE(16)    = 'INPUTDIR/EOF_COEF/crishrblk_pccoef_solar'   

     ABSCOEFFILE(17)= 'INPUTDIR/ABS_COEF/crishrbox_abscoef_solar'
     ICE_TAB(17)    = 'INPUTDIR/CLD_LUT/ice_rt_540_crishr_box_solar.dat'
     WAT_TAB(17)    = 'INPUTDIR/CLD_LUT/wat_rt_540_crishr_box_solar.dat'
     PARFILE(17)    = 'INPUTDIR/EOF_COEF/crishrbox_pccoef_solar'      

     ABSCOEFFILE(18)= 'INPUTDIR/ABS_COEF/crishrham_abscoef_solar'
     ICE_TAB(18)    = 'INPUTDIR/CLD_LUT/ice_rt_398_crishr_ham_solar.dat'
     WAT_TAB(18)    = 'INPUTDIR/CLD_LUT/wat_rt_398_crishr_ham_solar.dat'
     PARFILE(18)    = 'INPUTDIR/EOF_COEF/crishrham_pccoef_solar'

     ABSCOEFFILE(19)= 'INPUTDIR/ABS_COEF/nastibox_abscoef_solar'
     ICE_TAB(19)    = 'INPUTDIR/CLD_LUT/ice_rt_748_nasti_box_solar.dat'
     WAT_TAB(19)    = 'INPUTDIR/CLD_LUT/wat_rt_748_nasti_box_solar.dat'
     PARFILE(19)    = 'INPUTDIR/EOF_COEF/nastibox_pccoef_solar'

     ABSCOEFFILE(20)= 'INPUTDIR/ABS_COEF/nastiksr_abscoef_solar'
     ICE_TAB(20)    = 'INPUTDIR/CLD_LUT/ice_rt_559_nasti_ksr_solar.dat'
     WAT_TAB(20)    = 'INPUTDIR/CLD_LUT/wat_rt_559_nasti_ksr_solar.dat'
     PARFILE(20)    = 'INPUTDIR/EOF_COEF/nastiksr_pccoef_solar'

     ABSCOEFFILE(21)= 'INPUTDIR/ABS_COEF/shisbox_abscoef_solar'
     ICE_TAB(21)    = 'INPUTDIR/CLD_LUT/ice_rt_647_shis_box_solar.dat'
     WAT_TAB(21)    = 'INPUTDIR/CLD_LUT/wat_rt_647_shis_box_solar.dat'
     PARFILE(21)    = 'INPUTDIR/EOF_COEF/shisbox_pccoef_solar'

     ABSCOEFFILE(22)= 'INPUTDIR/ABS_COEF/shisksr_abscoef_solar'
     ICE_TAB(22)    = 'INPUTDIR/CLD_LUT/ice_rt_647_shis_ksr_solar.dat'
     WAT_TAB(22)    = 'INPUTDIR/CLD_LUT/wat_rt_647_shis_ksr_solar.dat'
     PARFILE(22)    = 'INPUTDIR/EOF_COEF/shisksr_pccoef_solar'

     ABSCOEFFILE(23)= 'INPUTDIR/ABS_COEF/crisbox_abscoef_solar'
     ICE_TAB(23)    = 'INPUTDIR/CLD_LUT/ice_rt_485_cris_box_solar.dat'
     WAT_TAB(23)    = 'INPUTDIR/CLD_LUT/wat_rt_485_cris_box_solar.dat'
     PARFILE(23)    = 'INPUTDIR/EOF_COEF/crisbox_pccoef_solar'      

     ABSCOEFFILE(24)= 'INPUTDIR/ABS_COEF/crisham_abscoef_solar'
     ICE_TAB(24)    = 'INPUTDIR/CLD_LUT/ice_rt_384_cris_ham_solar.dat'
     WAT_TAB(24)    = 'INPUTDIR/CLD_LUT/wat_rt_384_cris_ham_solar.dat'
     PARFILE(24)    = 'INPUTDIR/EOF_COEF/crisham_pccoef_solar'    

     ABSCOEFFILE(25)= 'INPUTDIR/ABS_COEF/crisblk_abscoef_solar'
     ICE_TAB(25)    = 'INPUTDIR/CLD_LUT/ice_rt_369_cris_blk_solar.dat'
     WAT_TAB(25)    = 'INPUTDIR/CLD_LUT/wat_rt_369_cris_blk_solar.dat'
     PARFILE(25)    = 'INPUTDIR/EOF_COEF/crisblk_pccoef_solar'       

     ABSCOEFFILE(26)= 'INPUTDIR/ABS_COEF/clarreo_abscoef_0.5'
     ICE_TAB(26)    = 'INPUTDIR/CLD_LUT/ice_clarreo_0.5.dat'
     WAT_TAB(26)    = 'INPUTDIR/CLD_LUT/wat_clarreo_0.5.dat'
     PARFILE(26)    = 'INPUTDIR/EOF_COEF/clarreo_0.5_pccoef_5121'

     ABSCOEFFILE(99)= 'INPUTDIR/ABS_COEF/abscoef'
     ICE_TAB(99)    = 'INPUTDIR/CLD_LUT/ice_rt.dat'
     WAT_TAB(99)    = 'INPUTDIR/CLD_LUT/wat_rt.dat'
     PARFILE(99)    = 'INPUTDIR/EOF_COEF/pccoef'     
     

        
     CALL RD_PCRTM_ATM_ABSCOEF( ATM_ABS_COEF,                      &
                                ABSCOEFFILE(SENSOR_ID),            &
                                PCRTM_STND )
     CALL INIT_PCRTM_ICECLD_GRID(ICE_GRID,PCRTM_STND%NM)
     CALL READ_PCRTM_CLD_TAB(ICE_TAB(SENSOR_ID),ICE_GRID,PCRTM_STND%NM)
     CALL INIT_PCRTM_WATCLD_GRID(WAT_GRID,PCRTM_STND%NM)
     CALL READ_PCRTM_CLD_TAB(WAT_TAB(SENSOR_ID),WAT_GRID,PCRTM_STND%NM)

     CALL RD_SENSOR_BND_INFO(EOF_SOLUTION,PARFILE(SENSOR_ID)) 

   END SUBROUTINE INIT_PCRTM_LUT


   SUBROUTINE INIT_PCRTM_PARAM( PCRTM_STND,     &
                                ATM,            &
                                RT_SOLUTION)         
 
     TYPE(PCRTM_ATM_ABS_STRUCT_TYPE), INTENT(IN)    :: PCRTM_STND
     TYPE(PCRTM_ATMOSPHERE_TYPE),     INTENT(OUT)   :: ATM
     TYPE(PCRTM_RT_SOLUTION_TYPE),    INTENT(OUT)   :: RT_SOLUTION


     CALL INIT_PCRTM_ATMOSPHERE( ATM, PCRTM_STND%NLAY, PCRTM_STND%NMOL )
     
     CALL INIT_PCRTM_RT_SOLUTION( RT_SOLUTION, PCRTM_STND )

   END SUBROUTINE INIT_PCRTM_PARAM


   SUBROUTINE INIT_PCRTM_LP( SENSOR_ID,    &
                             PCRTM_STND,   &
                             ICE_GRID,     &
                             WAT_GRID,     &
                             ATM_ABS_COEF, & 
                             EOF_SOLUTION, &
                             ATM,          &
                             RT_SOLUTION )         
 
     INTEGER,                        INTENT (IN) :: SENSOR_ID
     TYPE(PCRTM_ATM_ABS_STRUCT_TYPE),INTENT (OUT):: PCRTM_STND
     TYPE(PCRTM_CLD_TABLE_DEF),      INTENT (OUT):: ICE_GRID
     TYPE(PCRTM_CLD_TABLE_DEF),      INTENT (OUT):: WAT_GRID
     TYPE(PCRTM_ATM_ABSORPTION_TYPE),INTENT (OUT):: ATM_ABS_COEF 
     TYPE(PCRTM_EOF_SOLUTION_TYPE),   &
                        ALLOCATABLE, INTENT (OUT):: EOF_SOLUTION(:)
     TYPE(PCRTM_ATMOSPHERE_TYPE),    INTENT (OUT):: ATM
     TYPE(PCRTM_RT_SOLUTION_TYPE),   INTENT (OUT):: RT_SOLUTION

     CALL INIT_PCRTM_LUT( SENSOR_ID,       &
                          PCRTM_STND,      &
                          ICE_GRID,        &
                          WAT_GRID,        &
                          ATM_ABS_COEF,    &
                          EOF_SOLUTION )

     CALL INIT_PCRTM_PARAM( PCRTM_STND,    &
                            ATM,           &
                            RT_SOLUTION )  
          
   END SUBROUTINE INIT_PCRTM_LP

   SUBROUTINE INIT_PCRTM_JACOB( PCRTM_STND,     &
                                EOF_SOLUTION,   &
                                K_NM,           &
                                K_PC,           &
                                K_CH,           &
                                JACOB_CH,       &
                                JACOB_BT )
                                

     TYPE(PCRTM_NM_JACOBIAN_TYPE), INTENT(OUT)  :: K_NM
     TYPE(PCRTM_PC_JACOBIAN_TYPE), ALLOCATABLE, &
                                   INTENT(OUT)  :: K_PC(:)
     TYPE(PCRTM_CH_JACOBIAN_TYPE), ALLOCATABLE, &
                                   INTENT(OUT)  :: K_CH(:)

     TYPE(PCRTM_ATM_ABS_STRUCT_TYPE), INTENT(IN):: PCRTM_STND
     TYPE(PCRTM_EOF_SOLUTION_TYPE),   &
                         ALLOCATABLE, INTENT(IN):: EOF_SOLUTION(:)
     LOGICAL,                         INTENT(IN):: JACOB_CH
     LOGICAL,                         INTENT(IN):: JACOB_BT
     
     
     INTEGER                                    :: IB, NB
     INTEGER                                    :: ALLOC_STAT

     CALL INIT_PCRTM_JACOB_NM( K_NM, PCRTM_STND)    

     NB = SIZE(EOF_SOLUTION)

     ALLOCATE(K_PC(NB),STAT = ALLOC_STAT)
     IF ( ALLOC_STAT /= 0 ) THEN
        PRINT*,'ERROR TRYING TO ALLOCATE K_PC'
        STOP
     ENDIF

     DO IB = 1, NB
        K_PC(IB)%Jacob_PC = .true.
        CALL INIT_PCRTM_JACOB_PC( K_PC(IB), EOF_SOLUTION(IB)%NPCBND,PCRTM_STND%NLAY)
     END DO

     IF(JACOB_CH) THEN

        ALLOCATE(K_CH(NB),STAT = ALLOC_STAT)
        IF ( ALLOC_STAT /= 0 ) THEN
           PRINT*,'ERROR TRYING TO ALLOCATE K_CH'
           STOP
        ENDIF
        
        DO IB = 1, NB
           K_CH(IB)%Jacob_CH = .true.
           K_CH(IB)%Jacob_BTCH = Jacob_BT
           CALL INIT_PCRTM_JACOB_CH( K_CH(IB),                  &
                                     EOF_SOLUTION(IB)%NCHBND,   &
                                     PCRTM_STND%NLEV)
        END DO
        
     END IF

   END SUBROUTINE INIT_PCRTM_JACOB
   
   
 END MODULE INIT_PCRTM
