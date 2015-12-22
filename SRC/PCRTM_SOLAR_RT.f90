Module PCRTM_SOLAR_RT
  USE  PCRTM_SOLAR_PARAMETER
  USE  PCRTM_SOLAR_Define
  USE  PCRTM_RT_SOLUTION_define,only :  PCRTM_RT_SOLUTION_TYPE
  USE  PCRTM_CLOUD_DEFINE, ONLY   : PCRTM_CLOUD_TYPE
  USE  PCRTM_ATMOSPHERE_LAYER, only : PCRTM_GEOMETRY_TYPE
  
contains

  subroutine PCRTM_FORWARD_RT_M_Solar ( RT_SOLUTION,     &
                                        PCRTM_Stnd,      &
                                        NCLD,            &
                                        CLD,             &
                                        Geometry,        &
                                        Solar_tab,       &
                                        SOLAR_SOLUTION)
    Type(PCRTM_SOLAR_DEF),             INTENT(inout) :: SOLAR_SOLUTION
    TYPE(PCRTM_RT_SOLUTION_TYPE),      INTENT(in)    :: RT_SOLUTION
    TYPE(PCRTM_ATM_ABS_STRUCT_TYPE),   INTENT(IN)    :: PCRTM_STND
    integer,                           INTENT(in)    :: NCLD
    TYPE(PCRTM_CLOUD_TYPE),            INTENT(in)    :: CLD(NCLD)
    TYPE(PCRTM_GEOMETRY_TYPE),         INTENT(in)    :: Geometry
    type(PCRTM_SOLAR_LUT_DEF),         intent(in)    :: solar_tab(2)
    integer  ::  StartWaveIndex
    real*4   ::  wavenumber, tau_above, tau_below
    real*8   ::  view_mu, mu0, albedo

    StartWaveIndex = SOLAR_SOLUTION%StartWaveIndex
    view_mu = COS(Geometry%satang*pi/180.0)
    mu0     = COS(Geometry%solar_zang*pi/180.0)

!!$    CALL R_T_INTPOL(RT_SOLUTION%taugas_below(StartWaveIndex:),  &
!!$                    CLD(1)%VISTAU,                              &
!!$                    CLD(1)%De,                                  &
!!$                    Geometry%solar_zang,                        &
!!$                    Geometry%SATANG,                            &
!!$                    view_phi,                                   &
!!$                    cloud_type,                                 &
!!$                    solar_freq,                                 &
!!$                    start_freqgrid_index,                       &
!!$                    end_freqgrid_index,                         &
!!$                    numberfreqbins,                             &
!!$                    BRDF, Rdd, Tbd, Tdb)


!!!! Assign LUTs according to cloud type.
    IF (cld(1)%phase == 1) THEN
       CALL R_T_INTPOL(RT_SOLUTION,                        &
                       NCLD,                               &
                       CLD,                                &
                       Geometry,                           &
                       Solar_tab(1),                       &
                       SOLAR_SOLUTION,                     &
                       PCRTM_Stnd)
    ELSE IF (cld(1)%phase == 2) THEN
       CALL R_T_INTPOL(RT_SOLUTION,                        &
                       NCLD,                               &
                       CLD,                                &
                       Geometry,                           &
                       Solar_tab(2),                       &
                       SOLAR_SOLUTION,                     &
                       PCRTM_Stnd)
    END IF

    
    SOLAR_SOLUTION%SolarRadUp = 0.0
    DO iwave = 1, pcrtm_stnd%nM - StartWaveIndex + 1
       wavenumber = pcrtm_stnd%frq(iwave+StartWaveIndex-1)
       IF ((wavenumber .LT. 1800.00) .OR. (wavenumber .GT. 3000.00)) THEN
          PRINT *, 'Wavenumber must be larger than 1800 and smaller than 3000!'
          STOP
       ENDIF
       tau_above = RT_SOLUTION%Taugas_above(iwave+StartWaveIndex-1)
       tau_below = RT_SOLUTION%Taugas_below(iwave+StartWaveIndex-1)
       albedo    = (1.0 - RT_SOLUTION%EMIS(iwave+StartWaveIndex-1))*pi
       
       IF (NCLD .EQ. 0) THEN
          R_TOA = albedo*exp(-(tau_above+tau_below)/inci_mu)*&
               exp(-(tau_above+tau_below)/view_mu)
       ELSE
          R_TOA = (SOLAR_SOLUTION%BRDF(iwave) + albedo*SOLAR_SOLUTION%Tbd(iwave)&
               *SOLAR_SOLUTION%Tdb(iwave)/(1-albedo*SOLAR_SOLUTION%Rdd(iwave)/pi))*&
               exp(-tau_above/view_mu)*exp(-tau_above/mu0)
          SOLAR_SOLUTION%SolarRadUp(iwave) = 1000.0*mu0*SOLAR_SOLUTION%SolarSpectrum(iwave)*R_TOA/pi
       ENDIF
!!$       write(*,*) wavenumber, SOLAR_SOLUTION%BRDF(iwave),albedo,SOLAR_SOLUTION%Tbd(iwave),&
!!$            SOLAR_SOLUTION%Tdb(iwave),SOLAR_SOLUTION%Rdd(iwave)
    end DO


  end subroutine PCRTM_FORWARD_RT_M_SOLAR


!!$  SUBROUTINE R_T_INTPOL(bottomskytau, ctau, De, SZA, VZA, VAA,       &
!!$                      cloud_type, Wavenumber, start_freqgrid_index,&
!!$                      end_freqgrid_index, numberfreqbins,          &
!!$                      BRDF, Rdd, Tbd, Tdb)

  subroutine R_T_INTPOL(RT_SOLUTION,                    &
                        NCLD,                           &
                        CLD,                            &
                        Geometry,                       &
                        Solar_tab,                      &      
                        SOLAR_SOLUTION,                 &
                        PCRTM_Stnd)

    Type(PCRTM_SOLAR_DEF),          INTENT(inout) :: SOLAR_SOLUTION
    TYPE(PCRTM_RT_SOLUTION_TYPE),      INTENT(in) :: RT_SOLUTION
    TYPE(PCRTM_ATM_ABS_STRUCT_TYPE),   INTENT(IN) :: PCRTM_STND
    integer,                           INTENT(in) :: NCLD
    TYPE(PCRTM_CLOUD_TYPE),            INTENT(in) :: CLD(NCLD)
    TYPE(PCRTM_GEOMETRY_TYPE),         INTENT(in) :: Geometry
    type(PCRTM_SOLAR_LUT_DEF),         intent(in) :: solar_tab


!!$    REAL*8, DIMENSION(numberfreqbins), INTENT(OUT) :: BRDF, Rdd, Tbd, Tdb

    !INTERNAL

    REAL*8, allocatable :: Wavenumber(:)    !DIMENSION(numberfreqbins)
    REAL*4, allocatable :: bottomskytau(:)  !DIMENSION(numberfreqbins)
    REAL*8, allocatable :: Rdd_grid(:,:), Tbd_grid(:,:), Tdb_grid(:,:) !DIMENSION(maxskytau,end_freqgrid_index)
    REAL*8, allocatable :: BRDF_grid(:) !DIMENSION(end_freqgrid_index)

    INTEGER :: start_freqgrid_index, end_freqgrid_index, numberfreqbins,StartWaveIndex
    REAL*8  :: ctau, De
    REAL*8  :: SZA, VZA


    REAL*8 VAA !(NEED TO BE CHANGED IF > 180)

    REAL*8, DIMENSION(maxtau):: tau_ln
    REAL*8, DIMENSION(maxdef):: De_ln
    REAL*8, DIMENSION(size_SZA):: cSZA
    REAL*8, DIMENSION(size_VZA):: cVZA
    REAL*8, DIMENSION(size_VAA):: cVAA
    REAL*8, DIMENSION(maxwvl)::cWavenumber
    INTEGER index1_bottomskytau, index2_bottomskytau
    INTEGER index1_tau, index2_tau, index1_De, index2_De
    INTEGER index1_SZA, index2_SZA
    INTEGER index1_VZA, index2_VZA, index1_VAA, index2_VAA
    INTEGER index1_wave, index2_wave
    INTEGER i_unknown, i, j, ALLOC_STAT, DEALLOC_STAT
    REAL*8 fraction_bottomskytau, fraction_tau, fraction_De
    REAL*8 fraction_SZA, fraction_VZA, fraction_VAA, fraction_wave 
    REAL*8 uSZA, view_mu, Rdd_Out, Tbd_Out, Tdb_Out, BRDF_Out 


    REAL::timeSum, timex, dtime, tarray(2)
    REAL::dis_timeSum, dis_timex, dis_dtime
!!$bottomskytau, ctau, De, SZA, VZA, VAA,       &
!!$                      cloud_type, Wavenumber, start_freqgrid_index,&
!!$                      end_freqgrid_index, numberfreqbins,          &
!!$                      BRDF, Rdd, Tbd, Tdb)


    ctau = cld(1)%vistau
    De   = cld(1)%De
    SZA  = Geometry%solar_zang
    VZA  = Geometry%satang
    VAA  = Geometry%sat_azimuth-Geometry%solar_azimuth
    StartWaveIndex = SOLAR_SOLUTION%StartWaveIndex
    numberfreqbins = pcrtm_stnd%nM-StartWaveIndex+1

    allocate( bottomskytau(numberfreqbins) )
    bottomskytau = RT_SOLUTION%Taugas_below(StartWaveIndex:pcrtm_stnd%nM)

    allocate( wavenumber(numberfreqbins))
    Wavenumber = pcrtm_stnd%frq(StartWaveIndex:pcrtm_stnd%nM)
   
    start_freqgrid_index = (aint(pcrtm_stnd%frq(StartWaveIndex)/10)*10 - 1790)/10
    end_freqgrid_index = (aint(pcrtm_stnd%frq(pcrtm_stnd%nM)/10)*10 + 10 - 1790)/10

    allocate(Rdd_grid(maxskytau,end_freqgrid_index), &
             Tbd_grid(maxskytau,end_freqgrid_index), &
             Tdb_grid(maxskytau,end_freqgrid_index), &
             BRDF_grid(end_freqgrid_index))
    
!!!! grid points of the angles in LUTs
    DO j = 1, size_SZA
       cSZA(j) = REAL(COS((0.5+(j-1)*5.0)*pi/180.0)) 
    ENDDO
    DO i = 1, size_VZA
       cVZA(i) = REAL(COS((72.0 - i*4.0)*pi/180.0))
    ENDDO
    DO j = 1, size_VAA
       cVAA(j) = Real(5.0*(j-1))
    ENDDO
    
!!!! Cloud Type, Optical Depth and Particle Size Positions in Their Corresponding LUTs.
    tau_ln = log(solar_tab%tau_lut)
    IF (log(ctau) .LE. tau_ln(1)) THEN
       index1_tau = 1
       index2_tau = 2
    ELSE IF (log(ctau) .GE. tau_ln(maxtau)) THEN
       index1_tau = maxtau - 1
       index2_tau = maxtau
    ELSE
       DO i_unknown = 2, maxtau
          IF (log(ctau) .LE. tau_ln(i_unknown)) THEN
             index1_tau = i_unknown - 1
             index2_tau = i_unknown
             EXIT
          ENDIF
       ENDDO
    ENDIF
    fraction_tau = (log(ctau)-tau_ln(index1_tau)) / &
         (tau_ln(index2_tau)-tau_ln(index1_tau))
    
    De_ln(:)  = log(solar_tab%De_lut(:))
    IF (log(De) .LE. De_ln(1)) THEN
       index1_De = 1
       index2_De = 2
    ELSE IF (log(De) .GE. De_ln(maxdef)) THEN
       index1_De = maxdef - 1
       index2_De = maxdef
    ELSE
       DO i_unknown = 2, maxdef
          IF (log(De) .LE. De_ln(i_unknown)) THEN
             index1_De = i_unknown - 1
             index2_De = i_unknown
             EXIT
          ENDIF
       ENDDO
    ENDIF
    fraction_De = (log(De)-De_ln(index1_De)) / &
         (De_ln(index2_De)-De_ln(index1_De))
    
!!!! Solar Zenith Angle in TBD
    uSZA = cos(SZA*pi/180.0)
    IF (-uSZA .LE. -cSZA(1)) THEN
       index1_SZA = 1
       index2_SZA = 2
    ELSE IF (-uSZA .GE. -cSZA(size_SZA)) THEN
       index1_SZA = size_SZA - 1
       index2_SZA = size_SZA
    ELSE
       DO i_unknown = 2, size_SZA
          IF (-uSZA .LE. -cSZA(i_unknown)) THEN
             index1_SZA = i_unknown - 1
             index2_SZA = i_unknown
             EXIT
          ENDIF
       ENDDO
    ENDIF
    fraction_SZA = (-uSZA+cSZA(index1_SZA))/(-cSZA(index2_SZA)+cSZA(index1_SZA))
    
!!!! Wavenumber
    DO j = 1, maxwvl
       cWavenumber(j) = Real((j-1)*10.0+1800.0)
    EndDO
    
!!!! Satellite VZA Position in BRDF Tbd LUTs !!!!
    view_mu = cos(VZA*pi/180.0)
    IF (view_mu .LE. cVZA(1)) THEN
       index1_VZA = 1
       index2_VZA = 2
    ELSE IF (view_mu .GE. cVZA(size_VZA)) THEN
       index1_VZA = size_VZA - 1
       index2_VZA = size_VZA
    ELSE
       DO i_unknown = 2, size_VZA
          IF (view_mu .LE. cVZA(i_unknown)) THEN
             index1_VZA = i_unknown - 1
             index2_VZA = i_unknown
             EXIT
          ENDIF
       ENDDO
    ENDIF
    fraction_VZA = (view_mu-cVZA(index1_VZA))/(cVZA(index2_VZA)-cVZA(index1_VZA))
    
!!!! Satellite VAA Position
    IF (VAA .GT. 180.0) THEN
       VAA = 360.0 - VAA
    ENDIF
    IF (VAA .LE. cVAA(1)) THEN
       index1_VAA = 1
       index2_VAA = 2
    ELSE IF (VAA .GE. cVAA(size_VAA)) THEN
       index1_VAA = size_VAA - 1
       index2_VAA = size_VAA
    ELSE
       DO i_unknown = 2, size_VAA
          IF (VAA .LE. cVAA(i_unknown)) THEN
             index1_VAA = i_unknown - 1
             index2_VAA = i_unknown
             EXIT
          ENDIF
       ENDDO
    ENDIF
    fraction_VAA = (VAA-cVAA(index1_VAA))/(cVAA(index2_VAA)-cVAA(index1_VAA))
    
!!!! Start interpolation .... 
    DO i = start_freqgrid_index, end_freqgrid_index
       IF (ctau .LT. 0.01) THEN
          BRDF_grid(i) = 0.00
       ELSE
          Call  LIN_INT5(index1_tau,index2_tau,fraction_tau,&
               index1_De,index2_De,fraction_De,             &
               index1_SZA, index2_SZA, fraction_SZA,        &
               index1_VZA,index2_VZA, fraction_VZA,         &
               index1_VAA,index2_VAA, fraction_VAA,         &
               i, BRDF_Out,solar_tab)
          BRDF_grid(i) = BRDF_Out
!!$          write(*,*)i,BRDF_grid(i)
       ENDIF
    ENDDO
!!$print*, BRDF_grid
!!$stop
    
    DO j = 1, maxskytau
       DO i = start_freqgrid_index, end_freqgrid_index
          CALL LIN_INT2_0(index1_tau,index2_tau,fraction_tau, &
               index1_De, index2_De, fraction_De,             &
               j, i, Rdd_Out,solar_tab)
          Rdd_Grid(j,i) = Rdd_Out
          
          Call  LIN_INT3(index1_tau,index2_tau,fraction_tau, &
               index1_De, index2_De, fraction_De,            &
               index1_SZA, index2_SZA, fraction_SZA,         &
               j, i, 1, Tbd_Out,solar_tab)
          Tbd_Grid(j,i) = Tbd_Out
          
          CALL  LIN_INT3(index1_tau,index2_tau,fraction_tau, &
               index1_De, index2_De, fraction_De,            &
               index1_VZA,index2_VZA, fraction_VZA,          &
               j, i, 2, Tdb_Out,solar_tab)
          Tdb_Grid(j,i) = Tdb_Out
!!$          write(*,*)j,i,Rdd_Grid(j,i),Tdb_Grid(j,i),Tbd_Grid(j,i)
       ENDDO
    ENDDO

    DO i = 1, numberfreqbins
       IF (Wavenumber(i) .LE. cWavenumber(1)) THEN
          index1_wave = 1
          index2_wave = 2
       ELSE IF (Wavenumber(i) .GE. cWavenumber(maxwvl)) THEN
          index1_wave = maxwvl - 1
          index2_wave = maxwvl
       ELSE
          DO i_unknown = 2, maxwvl
             IF (Wavenumber(i) .LE. cWavenumber(i_unknown)) THEN
                index1_wave = i_unknown - 1
                index2_wave = i_unknown
                EXIT
             ENDIF
          ENDDO
       ENDIF
       fraction_wave = (Wavenumber(i)-cWavenumber(index1_wave))/(cWavenumber(index2_wave)-cWavenumber(index1_wave))
       
       SOLAR_SOLUTION%BRDF(i) = BRDF_grid(index1_wave) + (BRDF_grid(index2_wave) - BRDF_grid(index1_wave))*fraction_wave
       
       IF (bottomskytau(i) .LE. solar_tab%bottom_sky_tau_lut(1)) THEN
          index1_bottomskytau = 1
          index2_bottomskytau = 2
       ELSE IF (bottomskytau(i) .GE. solar_tab%bottom_sky_tau_lut(maxskytau)) THEN
          index1_bottomskytau = maxskytau - 1
          index2_bottomskytau = maxskytau
       ELSE
          DO i_unknown = 2, maxskytau
             IF (bottomskytau(i) .LE. solar_tab%bottom_sky_tau_lut(i_unknown)) THEN
                index1_bottomskytau = i_unknown - 1
                index2_bottomskytau = i_unknown
                EXIT
             ENDIF
          ENDDO
       ENDIF
       fraction_bottomskytau = (bottomskytau(i)-solar_tab%bottom_sky_tau_lut(index1_bottomskytau))/&
            (solar_tab%bottom_sky_tau_lut(index2_bottomskytau)-solar_tab%bottom_sky_tau_lut(index1_bottomskytau))
!!$       print*, 'fraction_bottomskytau',bottomskytau(i), index1_bottomskytau, index2_bottomskytau,fraction_bottomskytau
       
       CALL LIN_INT2(index1_bottomskytau,index2_bottomskytau,fraction_bottomskytau,&
            index1_wave, index2_wave, fraction_wave, Rdd_Out, end_freqgrid_index, Rdd_grid)
       SOLAR_SOLUTION%Rdd(i) = Rdd_Out
       CALL LIN_INT2(index1_bottomskytau,index2_bottomskytau,fraction_bottomskytau,&
            index1_wave, index2_wave, fraction_wave, Tbd_Out, end_freqgrid_index, Tbd_grid)
       SOLAR_SOLUTION%Tbd(i) = Tbd_Out
       CALL LIN_INT2(index1_bottomskytau,index2_bottomskytau,fraction_bottomskytau,&
            index1_wave, index2_wave, fraction_wave, Tdb_Out, end_freqgrid_index, Tdb_grid)
       SOLAR_SOLUTION%Tdb(i) = Tdb_Out 
!!$       print*,i, Wavenumber(i), index1_wave, index2_wave, fraction_wave, &
!!$            SOLAR_SOLUTION%BRDF(i),SOLAR_SOLUTION%Rdd(i),SOLAR_SOLUTION%Tbd(i),SOLAR_SOLUTION%Tdb(i)
    EndDO
    
    deallocate( bottomskytau )
    deallocate( wavenumber )
    deallocate( Rdd_grid, Tbd_grid, Tdb_grid, BRDF_grid)

    
  END SUBROUTINE R_T_INTPOL
  
  
  SUBROUTINE LIN_INT2_0(X1_D1,X1_D2,F1,X2_D1,X2_D2,F2,X3,X4,R, solar_tab)

    INTEGER X1_D1, X1_D2, X2_D1, X2_D2, X1, X2, X3, X4, i, j 
    REAL*8 F1, F2
    REAL*8, DIMENSION(2,2)::R2
    REAL*8, DIMENSION(2)::R1
    REAL*8  R
    type(PCRTM_SOLAR_LUT_DEF),         intent(in) :: solar_tab

    DO i = 1, 2
       if (i .eq. 1) then
          X1 = X1_D1
       else
          X1 = X1_D2
       endif
       DO j = 1, 2
          if (j .eq. 1) then
             X2 = X2_D1
          else
             X2 = X2_D2
          endif
          R2(i,j) = solar_tab%Rdd_lut(X3,X1,X2,X4)
       ENDDO
    END DO
    DO j = 1, 2
       R1(j) = R2(1,j) + (R2(2,j) - R2(1,j)) * F1
    ENDDO
    R = R1(1) + (R1(2) - R1(1)) * F2
    
  END SUBROUTINE LIN_INT2_0
  
  
  SUBROUTINE LIN_INT2(X1_D1,X1_D2,F1,X2_D1,X2_D2,F2,RT,end_freqgrid_index, RT_grid)

    INTEGER X1_D1, X1_D2, X2_D1, X2_D2, X1, X2, i, j 
    INTEGER end_freqgrid_index
    REAL*8 F1, F2
    REAL*8, allocatable :: RT_grid(:,:)
    REAL*8, DIMENSION(2,2)::R2
    REAL*8, DIMENSION(2)::R1
    REAL*8  RT
    
    DO i = 1, 2
       if (i .eq. 1) then
          X1 = X1_D1
       else
          X1 = X1_D2
       endif
       DO j = 1, 2
          if (j .eq. 1) then
             X2 = X2_D1
          else
             X2 = X2_D2
          endif
          R2(i,j) = RT_grid(X1,X2)
       ENDDO
    END DO
    
    DO j = 1, 2
       R1(j) = R2(1,j) + (R2(2,j) - R2(1,j)) * F1
    ENDDO
    RT = R1(1) + (R1(2) - R1(1)) * F2
    
    
  END SUBROUTINE LIN_INT2
  
  
  SUBROUTINE LIN_INT3(X2_D1,X2_D2,F2,X3_D1,X3_D2,F3,X4_D1,X4_D2,F4,X1,X5,id,T,solar_tab)

    IMPLICIT NONE
    INTEGER X2_D1, X2_D2, X3_D1, X3_D2, X4_D1, X4_D2
    INTEGER X1, X2, X3, X4, X5,i, j, k, l,q,id
    REAL*8 F1, F2, F3, F4
    REAL*8, DIMENSION(2,2,2,2)::T4
    REAL*8, DIMENSION(2,2,2)::T3
    REAL*8, DIMENSION(2,2)::T2
    REAL*8, DIMENSION(2)::T1
    REAL*8  T
    type(PCRTM_SOLAR_LUT_DEF),         intent(in) :: solar_tab


    DO j = 1, 2
       if (j .eq. 1) then
          X2 = X2_D1
       else
          X2 = X2_D2
       endif
       DO k = 1, 2
          if (k .eq. 1) then
             X3 = X3_D1
          else
             X3 = X3_D2
          endif
          DO l = 1, 2
             if (l.eq.1) then
                X4 = X4_D1
             else
                X4 = X4_D2
             endif
             IF (id.eq.2) THEN
                T3(j,k,l) = DOT_PRODUCT(solar_tab%MTDBA(X5,:),solar_tab%MTDBB(:,X1,X2,X3,X4))
             ELSE
                T3(j,k,l) = DOT_PRODUCT(solar_tab%MTBDA(X5,:),solar_tab%MTBDB(:,X1,X2,X3,X4))
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    
    DO k = 1, 2
       DO l = 1, 2
          T2(k,l) = T3(1,k,l) + (T3(2,k,l) - T3(1,k,l))*F2
       ENDDO
    ENDDO
    DO l = 1, 2
       T1(l) = T2(1,l) + (T2(2,l) - T2(1,l))*F3
    ENDDO
    T = T1(1) + (T1(2) - T1(1)) * F4
    
  END SUBROUTINE LIN_INT3

!!$SUBROUTINE LIN_INT4(X1_D1,X1_D2,F1,X2_D1,X2_D2,F2,X3_D1,X3_D2,F3,X4_D1,X4_D2,F4,X5,id,T)
!!$IMPLICIT NONE
!!$INTEGER X1_D1, X1_D2, X2_D1, X2_D2, X3_D1, X3_D2, X4_D1, X4_D2
!!$INTEGER X1, X2, X3, X4, X5,i, j, k, l,q,id
!!$REAL*8 F1, F2, F3, F4
!!$REAL*8, DIMENSION(2,2,2,2)::T4
!!$REAL*8, DIMENSION(2,2,2)::T3
!!$REAL*8, DIMENSION(2,2)::T2
!!$REAL*8, DIMENSION(2)::T1
!!$REAL*8  T
!!$
!!$DO i = 1, 2
!!$   if (i .eq. 1) then
!!$      X1 = X1_D1
!!$   else
!!$      X1 = X1_D2
!!$   endif
!!$   DO j = 1, 2
!!$      if (j .eq. 1) then
!!$         X2 = X2_D1
!!$      else
!!$         X2 = X2_D2
!!$      endif
!!$      DO k = 1, 2
!!$         if (k .eq. 1) then
!!$            X3 = X3_D1
!!$         else
!!$            X3 = X3_D2
!!$         endif
!!$         DO l = 1, 2
!!$            if (l.eq.1) then
!!$               X4 = X4_D1
!!$            else
!!$               X4 = X4_D2
!!$            endif
!!$            IF (id.eq.2) THEN
!!$               T4(i,j,k,l) = DOT_PRODUCT(MTDBA(X5,:),MTDBB(:,X1,X2,X3,X4))
!!$            ELSE
!!$               T4(i,j,k,l) = DOT_PRODUCT(MTBDA(X5,:),MTBDB(:,X1,X2,X3,X4))
!!$            ENDIF
!!$         ENDDO
!!$      ENDDO
!!$   ENDDO
!!$ENDDO
!!$
!!$DO j = 1, 2
!!$   DO k = 1, 2
!!$      DO l = 1, 2
!!$         T3(j,k,l) = T4(1,j,k,l) + (T4(2,j,k,l) - T4(1,j,k,l))*F1
!!$      ENDDO
!!$   ENDDO
!!$ENDDO
!!$DO k = 1, 2
!!$   DO l = 1, 2
!!$      T2(k,l) = T3(1,k,l) + (T3(2,k,l) - T3(1,k,l))*F2
!!$   ENDDO
!!$ENDDO
!!$DO l = 1, 2
!!$   T1(l) = T2(1,l) + (T2(2,l) - T2(1,l))*F3
!!$ENDDO
!!$T = T1(1) + (T1(2) - T1(1)) * F4
!!$
!!$END



  SUBROUTINE LIN_INT5(X1_D1,X1_D2,F1,X2_D1,X2_D2,F2,X3_D1,X3_D2,F3,X4_D1,X4_D2,F4,&
       X5_D1,X5_D2,F5,X6,BRDF,solar_tab)
    

    INTEGER X1_D1, X1_D2, X2_D1, X2_D2, X3_D1, X3_D2, X4_D1, X4_D2, X5_D1,X5_D2
    INTEGER X1, X2, X3, X4, X5, X6, i, j, k, l, m 
    REAL*8 F1, F2, F3, F4, F5
    REAL*8, DIMENSION(2,2,2,2,2)::R5
    REAL*8, DIMENSION(2,2,2,2)::R4
    REAL*8, DIMENSION(2,2,2)::R3
    REAL*8, DIMENSION(2,2)::R2
    REAL*8, DIMENSION(2)::R1
    REAL*8  BRDF
    type(PCRTM_SOLAR_LUT_DEF),         intent(in) :: solar_tab

    DO i = 1, 2
       if (i .eq. 1) then
          X1 = X1_D1
       else
          X1 = X1_D2
       endif
       DO j = 1, 2
          if (j .eq. 1) then
             X2 = X2_D1
          else
             X2 = X2_D2
          endif
          DO k = 1, 2
             if (k .eq. 1) then
                X3 = X3_D1
             else
                X3 = X3_D2
             endif
             DO l = 1, 2
                if (l .eq. 1) then
                   X4 = X4_D1
                else
                   X4 = X4_D2
                endif
                DO m = 1, 2
                   if (m .eq. 1) then
                      X5 = X5_D1
                   else
                      X5 = X5_D2
                   endif
                   R5(i,j,k,l,m) = DOT_PRODUCT(solar_tab%MBRDFA(X1,X4,X5,1:100),solar_tab%MBRDFB(1:100,X2,X3,X6))
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    
!!!! Interpolation
    DO j = 1, 2
       DO k = 1, 2
          DO l = 1, 2
             DO m = 1, 2
                R4(j,k,l,m) = R5(1,j,k,l,m) + (R5(2,j,k,l,m) - R5(1,j,k,l,m))*F1
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    DO k = 1, 2
       DO l = 1, 2
          DO m = 1, 2
             R3(k,l,m) = R4(1,k,l,m) + (R4(2,k,l,m) - R4(1,k,l,m))*F2
          ENDDO
       ENDDO
    ENDDO
!!$print*, F2
!!$stop
    DO l = 1, 2
       DO m = 1, 2
          R2(l,m) = R3(1,l,m) + (R3(2,l,m) - R3(1,l,m))*F3
       ENDDO
    ENDDO
    DO m = 1, 2
       R1(m) = R2(1,m) + (R2(2,m) - R2(1,m))*F4
    ENDDO
    BRDF = R1(1) + (R1(2) - R1(1)) * F5
!!$print*, F1, F2, F3, F4, F5,BRDF

    
  END SUBROUTINE LIN_INT5
  
end Module PCRTM_SOLAR_RT
