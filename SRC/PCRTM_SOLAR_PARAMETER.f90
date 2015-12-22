MODULE PCRTM_SOLAR_PARAMETER

  USE PCRTM_FILE_UTILITIES

  INTEGER,PARAMETER:: maxskytau = 20, maxtau = 30,  maxdef = 18,  maxwvl = 121
  INTEGER,PARAMETER:: size_SZA = 18, size_VZA = 18, size_VAA = 37, size_Wave = 1
  INTEGER,PARAMETER:: STdb = 15
  INTEGER,PARAMETER:: SBRDF = 100
  INTEGER,PARAMETER:: STbd = 15
  REAL*8,PARAMETER::pi = 3.141592653589793 !acos(-1.0d0)

  type :: PCRTM_SOLAR_LUT_DEF
     !Parameters for BRDF, TDB, and TBD LUTs!
     !TDB LUTs!  
     REAL*4, allocatable :: MTDBA(:,:)
     REAL*4, allocatable :: MTDBB(:,:,:,:,:)
     !TBD LUTs!  
     REAL*4, allocatable :: MTBDA(:,:)
     REAL*4, allocatable :: MTBDB(:,:,:,:,:)
     !BRDF LUTs!
     REAL*4, allocatable :: MBRDFA(:,:,:,:)
     REAL*4, allocatable :: MBRDFB(:,:,:,:)
     REAL*4, allocatable :: Rdd_lut(:,:,:,:)
     REAL*8,allocatable::de_lut(:),tau_lut(:),bottom_sky_tau_lut(:) 
     
  end type PCRTM_SOLAR_LUT_DEF

!!!! Original LUT
!!$  DATA tau_lut /.01, .05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, &
!!$                 1.2, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, &
!!$                 9.0, 10., 11., 12., 15., 20./
!!$  
!!$  DATA def_ic_lut /10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.,&
!!$                 110.,120.,130.,140.,150.,160.,170.,180./
!!$  
!!$  DATA def_wc_lut /4.0,8.0,12.0,16.0,20.0,24.0,28.0,32.0,36.0,40.,&
!!$                   50.,60.,70.,80.,100.,120.,160.,200./
!!$
!!$  DATA bottom_sky_tau_lut /0.0, .01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,&
!!$                   1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0/

contains
  subroutine INIT_solar_brdf_tab(solar_tab)
    
    type(PCRTM_SOLAR_LUT_DEF),intent(out) :: solar_tab(2)
    character(256)                        :: lutname
    integer                               :: lun, alloc_stat, ib

    do ib = 1, 2
       allocate(solar_tab(ib)%MTDBA(maxwvl,STdb),                              &
                solar_tab(ib)%MTDBB(STdb,maxskytau,maxtau,maxdef, size_VZA),   &
                solar_tab(ib)%MTBDA(maxwvl,STbd),                              &
                solar_tab(ib)%MTBDB(STbd,maxskytau,maxtau,maxdef,size_SZA),    &
                solar_tab(ib)%MBRDFA(maxtau,size_VZA,size_VAA,SBRDF),          &
                solar_tab(ib)%MBRDFB(SBRDF,maxdef,size_SZA,maxwvl),            &
                solar_tab(ib)%RDD_lut(maxskytau,maxtau,maxdef,maxwvl),         &
                solar_tab(ib)%de_lut(maxdef),                                  &
                solar_tab(ib)%tau_lut(maxtau),                                 &
                solar_tab(ib)%bottom_sky_tau_lut(maxskytau),                   &
                stat = alloc_stat )                
       IF (ALLOC_STAT /= 0) THEN
          PRINT*,'ERROR TRYING TO ALLOCATE SOLAR_TAB'
          STOP
       ENDIF
    end do

    solar_tab(1)%tau_lut = (/.01, .05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, &
                            1.2, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, &
                            9.0, 10., 11., 12., 15., 20./)
    solar_tab(1)%bottom_sky_tau_lut = (/0.0, .01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,&
                                       1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0/)
    solar_tab(1)%de_lut  = (/10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.,&
                            110.,120.,130.,140.,150.,160.,170.,180./)
    call getlun(lun)
    ! Read BRDF !
    WRITE(lutname, '(A)')'INPUTDIR/CLD_LUT/Ice_Corrected_BRDF_M_A_1800_3000.dat'
    OPEN(lun,FILE=lutname,FORM='unformatted',ACCESS='direct',RECL=maxtau*size_VZA*size_VAA*SBRDF*4)
    READ(lun,REC=1)solar_tab(1)%MBRDFA
    CLOSE(lun)
    WRITE(lutname, '(A)')'INPUTDIR/CLD_LUT/Ice_Corrected_BRDF_M_B_1800_3000.dat'
    OPEN(lun,FILE=lutname,FORM='unformatted',ACCESS='direct',RECL=SBRDF*maxdef*size_SZA*maxwvl*4)
    READ(lun,REC=1)solar_tab(1)%MBRDFB
    CLOSE(lun)
    ! read RDD
    WRITE(lutname,'(A)')'INPUTDIR/CLD_LUT/Rdd_ic_20x30x18x121.dat'
    OPEN(lun,FILE=lutname,FORM='unformatted',ACCESS='direct',RECL=maxskytau*maxtau*maxdef*maxwvl*4)
    READ(lun,REC=1)solar_tab(1)%Rdd_lut
    CLOSE(lun)
    ! Read TBD and TDB LUTs !
    WRITE(lutname,'(A)')'INPUTDIR/CLD_LUT/Ice_Tbd_M_A_1800_3000.dat'
    OPEN(lun,FILE=lutname,FORM='unformatted',ACCESS='direct',RECL=maxwvl*STbd*4)
    READ(lun,REC=1)solar_tab(1)%MTBDA
    CLOSE(lun)   
    WRITE(lutname,'(A)')'INPUTDIR/CLD_LUT/Ice_Tbd_M_B_1800_3000.dat'
    OPEN(lun,FILE=lutname,FORM='unformatted',ACCESS='direct',RECL=STbd*maxskytau*maxtau*maxdef*size_SZA*4)
    READ(lun,REC=1)solar_tab(1)%MTBDB
    CLOSE(lun)
    WRITE(lutname,'(A)')'INPUTDIR/CLD_LUT/Ice_Tdb_M_A_1800_3000.dat'
    OPEN(lun,FILE=lutname,FORM='unformatted',ACCESS='direct',RECL=maxwvl*STdb*4)
    READ(lun,REC=1)solar_tab(1)%MTDBA
    CLOSE(lun)   
    WRITE(lutname,'(A)')'INPUTDIR/CLD_LUT/Ice_Tdb_M_B_1800_3000.dat'
    OPEN(lun,FILE=lutname,FORM='unformatted',ACCESS='direct',RECL=STdb*maxskytau*maxtau*maxdef*size_VZA*4)
    READ(lun,REC=1)solar_tab(1)%MTDBB
    CLOSE(lun)
    
    
    solar_tab(2)%tau_lut =  solar_tab(1)%tau_lut
    solar_tab(2)%bottom_sky_tau_lut = solar_tab(1)%bottom_sky_tau_lut
    solar_tab(2)%de_lut  = (/4.0,8.0,12.0,16.0,20.0,24.0,28.0,32.0,36.0,40.,&
                            50.,60.,70.,80.,100.,120.,160.,200./)
    
    
  end subroutine INIT_solar_brdf_tab

  subroutine clear_solar_brdf_tab(solar_tab)

    type(PCRTM_SOLAR_LUT_DEF), &
                         intent(inout) :: solar_tab(2)
    integer                            :: ib, dealloc_stat


    do ib = 1, 2
       deallocate(solar_tab(ib)%MTDBA,                              &
                  solar_tab(ib)%MTDBB,                              &
                  solar_tab(ib)%MTBDA,                              &
                  solar_tab(ib)%MTBDB,                              &
                  solar_tab(ib)%MBRDFA,                             &
                  solar_tab(ib)%MBRDFB,                             &
                  solar_tab(ib)%RDD_lut,                            &
                  solar_tab(ib)%de_lut,                             &
                  solar_tab(ib)%tau_lut,                            &
                  solar_tab(ib)%bottom_sky_tau_lut,                 &
                  stat = dealloc_stat )                
       IF (DEALLOC_STAT /= 0) THEN
          PRINT*,'ERROR TRYING TO DEALLOCATE SOLAR_TAB'
          STOP
       ENDIF
    end do
  end subroutine clear_solar_brdf_tab


END MODULE PCRTM_SOLAR_PARAMETER
