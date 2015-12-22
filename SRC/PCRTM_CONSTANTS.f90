Module PCRTM_constants

!Module ATM_constants
  use PCRTM_Type_Kind

!!$    absorption Molecule Name =(/'H2O','CO2','O3 ','N2O','CO ','CH4', &
!!$                                'O2', 'NO', 'SO2', 'NO2', 'NH3', 'HNO3' &
!!$                                'CCL4', 'CCl3F', 'CCl2F2'/)

  integer, parameter            :: mxmol    = 15 
  real(single), dimension(mxMol):: WtMol    = (/18.015,  44.010,  47.998,  44.01, &
                                               28.011,  16.043,  31.999,  30.01, &
                                               64.06,   46.01,   17.03,   63.01, &
                                               53.82,  137.37,  120.91 /)
  real(single), parameter       :: WtDryAir = 28.97
  real(single), parameter       :: R_gas=8.3143
  real(single), parameter       :: R_air=0.9975*R_gas

   real(single), parameter :: Cplanck = 6.626176e-27 !Planck Constant (mW)
   real(single), parameter :: Vlight = 2.997925e+10 ! Velocity of (cm/s)
   real(single), parameter :: Baltz = 1.380662e-16 !Baltzman Constant
   real(single), parameter :: c1 = 2.0*Cplanck*Vlight*Vlight
   real(single), parameter :: c2 = Cplanck*Vlight/Baltz

   REAL(single), PARAMETER :: R_EQUATOR = 6.378388e+06
   REAL(single), PARAMETER :: R_POLAR = 6.356911e+06
   REAL(single), PARAMETER :: R_AVG = 0.5*(R_EQUATOR+R_POLAR)
   REAL(single), PARAMETER :: g_sfc = 9.80665
   REAl(single), dimension(101) :: pbnd_101 =  &
                            (/   0.0050,    0.0161,    0.0384,    0.0769,    0.1370,   &
                                 0.2244,    0.3454,    0.5064,    0.7140,    0.9753,   &
                                 1.2972,    1.6872,    2.1526,    2.7009,    3.3398,   &
                                 4.0770,    4.9204,    5.8776,    6.9567,    8.1655,   &
                                 9.5119,   11.0038,   12.6492,   14.4559,   16.4318,   &
                                 18.5847,   20.9224,   23.4526,   26.1829,   29.1210,  &
                                 32.2744,   35.6505,   39.2566,   43.1001,   47.1882,  &
                                 51.5278,   56.1260,   60.9895,   66.1253,   71.5398,  &
                                 77.2396,   83.2310,   89.5204,   96.1138,  103.0172,  &
                                 110.2366,  117.7775,  125.6456,  133.8462,  142.3848, &
                                 151.2664,  160.4959,  170.0784,  180.0183,  190.3203, &
                                 200.9887,  212.0277,  223.4415,  235.2338,  247.4085, &
                                 259.9691,  272.9191,  286.2617,  300.0000,  314.1369, &
                                 328.6753,  343.6176,  358.9665,  374.7241,  390.8926, &
                                 407.4738,  424.4698,  441.8819,  459.7118,  477.9607, &
                                 496.6298,  515.7200,  535.2322,  555.1669,  575.5248, &
                                 596.3062,  617.5112,  639.1398,  661.1920,  683.6673, &
                                 706.5654,  729.8857,  753.6275,  777.7897,  802.3714, &
                                 827.3713,  852.7880,  878.6201,  904.8659,  931.5236, &
                                 958.5911,  986.0666, 1013.9476, 1042.2319, 1070.9170, &
                                 1100.0000 /)        
  

 end Module PCRTM_constants
