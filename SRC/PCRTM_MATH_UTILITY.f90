Module PCRTM_math_utility
  use PCRTM_Type_Kind
  use PCRTM_constants, ONLY : c1, c2


contains
!#############################################################################
! SUBROUTINE LogXLogP
!#############################################################################      
  subroutine LogXLogP(x1,x2,y1,y2,y,x)

    IMPLICIT NONE

!---  Subroutine arguments
    REAL(SINGLE), INTENT(IN)  :: x1
    REAL(SINGLE), INTENT(IN)  :: x2
    REAL(SINGLE), INTENT(IN)  :: y1
    REAL(SINGLE), INTENT(IN)  :: y2
    REAL(SINGLE), INTENT(IN)  :: y
    REAL(SINGLE), INTENT(OUT) :: x
    
!---  Local variable
    REAL    :: A

!ed. alog has been superceeded by LOG, intrinsic natural log function
    A = alog(y/y1)/alog(y2/y1)
    x = x1*(x2/x1)**A
    
    RETURN
  END subroutine LogXLogP


!#############################################################################
! SUBROUTINE LogXLogP1
!#############################################################################      
  subroutine  LogXLogP1(x1,x2,y1,y2,y,x,dxdy,A)

    IMPLICIT NONE
      
!---  Subroutine arguments
    REAL(SINGLE), INTENT(IN)  :: x1
    REAL(SINGLE), INTENT(IN)  :: x2
    REAL(SINGLE), INTENT(IN)  :: y1
    REAL(SINGLE), INTENT(IN)  :: y2
    REAL(SINGLE), INTENT(IN)  :: y
    REAL(SINGLE), INTENT(OUT) :: x
    REAL(SINGLE), INTENT(OUT) :: dxdy
    REAL(SINGLE), INTENT(OUT) :: A

!ed. alog has been superceeded by LOG, intrinsic natural log function
    A = alog(y/y1)/alog(y2/y1)
    x = x1*(x2/x1)**A
    dxdy = x*alog(x2/x1)/alog(y2/y1)/y
    
    RETURN
  END subroutine LogXLogP1

!#############################################################################  
! SUBROUTINE LinXLinP
!#############################################################################       
  subroutine linXlinP(x1,x2,y1,y2,y,x)
      IMPLICIT NONE
      
!---  Subroutine arguments
      REAL(SINGLE), INTENT(IN)  :: x1
      REAL(SINGLE), INTENT(IN)  :: x2
      REAL(SINGLE), INTENT(IN)  :: y1
      REAL(SINGLE), INTENT(IN)  :: y2
      REAL(SINGLE), INTENT(IN)  :: y
      REAL(SINGLE), INTENT(OUT) :: x

!---  Local variable
      REAL(SINGLE)    :: A

      A = (y-y1)/(y2-y1)
      x = x1+(x2-x1)*A
      
      RETURN
    END subroutine linXlinP
      

!#############################################################################  
! SUBROUTINE LinXLinP1
!#############################################################################       
   subroutine linXlinP1(x1,x2,y1,y2,y,x,dxdy)
      IMPLICIT NONE
      
!---  Subroutine arguments
      REAL(SINGLE), INTENT(IN)  :: x1
      REAL(SINGLE), INTENT(IN)  :: x2
      REAL(SINGLE), INTENT(IN)  :: y1
      REAL(SINGLE), INTENT(IN)  :: y2
      REAL(SINGLE), INTENT(IN)  :: y
      REAL(SINGLE), INTENT(OUT) :: x
      REAL(SINGLE), INTENT(OUT) :: dxdy

!---  Local variable
      REAL(SINGLE)    :: A

      A = (y-y1)/(y2-y1)
      x = x1+(x2-x1)*A
      dxdy = 1/(y2-y1)
      
      RETURN
    END subroutine linXlinP1
      


!#############################################################################
! FUNCTION Planck
!
! Convert temperature (K) to planck radiance (mW*cm/m2/sr) at frq (cm-1)
!#############################################################################
    function Planck(frq,T,N)
      
      IMPLICIT NONE
      
      !---  Function arguments
      integer,        intent(IN) :: N
      REAL(double),   INTENT(IN) :: frq(N)
      REAL(SINGLE),   INTENT(IN) :: T
      
      !---  Function return value
      REAL(SINGLE)               :: Planck(N)
      
      !---  Local variables 
      REAL(SINGLE)               :: tmp
      REAL(SINGLE)               :: tmp1(N)
      REAL(SINGLE)               :: tmp2(N)
      
      tmp    = T
      tmp1   = c1*frq**3
      tmp2   = EXP(c2*frq/tmp)
      Planck = tmp1/(tmp2-1.d0)
      
      RETURN
    end function Planck

    function Bright_Temp(frq,Radch,N)
      
      IMPLICIT NONE
      
      !---  Function arguments
      integer,        intent(IN) :: N
      REAL(double),   INTENT(IN) :: frq(N)
      REAL(SINGLE),   INTENT(IN) :: Radch(N)
      
      !---  Function return value
      REAL(SINGLE)               :: Bright_Temp(N)
      
      !---  Local variables 
      REAL(SINGLE)               :: tmp1(N)

      tmp1   = c1*frq**3
      Bright_Temp = c2*frq/log(tmp1/Radch + 1)

    END function Bright_Temp

    Function dBTdRad(frq,Rad,N)
      integer,        intent(IN) :: N
      REAL(double),   INTENT(IN) :: frq(N)
      REAL(SINGLE),   INTENT(IN) :: Rad(N)      

      !---  Function return value
      REAL(SINGLE)               :: dBTdRad(N)

      !---  Local variables 
      REAL(SINGLE)               :: tmp1(N)      
      
      tmp1    = c1*frq**3
      dBTdRad = c2*frq*tmp1/(log(tmp1/Rad + 1))**2/(tmp1 + Rad)/Rad

    END Function dBTdRad


!#############################################################################
      function dPlanckdT(frq,T,N)
      implicit none
      integer     :: N
      real(double):: frq(N)
      real(single):: T
      real(single):: dPlanckdT(N)
      real(single):: Tmp,Tmp1(N)
 
      Tmp=1.d0/T
      Tmp1=exp(c2*frq*Tmp)
      dPlanckdT=c1*c2*Tmp1*(frq*frq*Tmp/(Tmp1-1.d0))**2
      return
      end function dPlanckdT

!#############################################################################

      subroutine blend_103 ( r, s, t, x000, x001, x010, x011, x100, x101, x110, &
           x111, x )
    
!*******************************************************************************
!
!! BLEND_103 extends scalar point data into a cube.
!
!
!  Diagram:
!
!    011--------------111 
!      |               |
!      |               | 
!      |               |
!      |               |
!      |               |
!    001--------------101
!
!
!      *---------------*
!      |               |
!      |               |
!      |      rst      |
!      |               |
!      |               |
!      *---------------*
!
!
!    010--------------110
!      |               |
!      |               |
!      |               |
!      |               | 
!      |               |
!    000--------------100 
!
!
!  Formula:
!
!    Written as a polynomial in R, S and T, the interpolation map has the 
!    form:
!
!      X(R,S,T) =
!        1         * ( + x000 )
!      + r         * ( - x000 + x100 )
!      +     s     * ( - x000        + x010 )
!      +         t * ( - x000               + x001 )
!      + r * s     * ( + x000 - x100 - x010                       + x110 )
!      + r     * t * ( + x000 - x100        - x001        + x101 )
!      +     s * t * ( + x000        - x010 - x001 + x011 )
!      + r * s * t * ( - x000 + x100 + x010 + x001 - x011 - x101 - x110 + x111 )
!
!  Reference:
!
!    William Gordon,
!    Blending-Function Methods of Bivariate and Multivariate Interpolation
!      and Approximation,
!    SIAM Journal on Numerical Analysis,
!    Volume 8, Number 1, March 1971, pages 158-177.
!
!    William Gordon and Charles Hall,
!    Transfinite Element Methods: Blending-Function Interpolation over
!      Arbitrary Curved Element Domains,
!    Numerische Mathematik,
!    Volume 21, Number 1, 1973, pages 109-129.
!
!    William Gordon and Charles Hall,
!    Construction of Curvilinear Coordinate Systems and Application to
!      Mesh Generation,
!    International Journal of Numerical Methods in Engineering,
!    Volume 7, 1973, pages 461-477.
!
!    Joe Thompson, Bharat Soni, Nigel Weatherill,
!    Handbook of Grid Generation,
!    CRC Press, 1999.
!
!  Modified:
!
!    11 October 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real R, S, T, the coordinates where an interpolated value
!    is desired.
!
!    Input, real(single) X000, X001, X010, X011, X100, X101, X110, X111, the
!    data values at the corners.
!
!    Output, real X, the interpolated data value at (R,S,T).
!
        implicit none
        !
        real(single) :: r
        real(single) :: s
        real(single) :: t
        real(single) :: x
        real(single) :: x000
        real(single) :: x001
        real(single) :: x010
        real(single) :: x011
        real(single) :: x100
        real(single) :: x101
        real(single) :: x110
        real(single) :: x111
!
!  Interpolate the interior point.
!
        x = &
             1.0E+00     * ( + x000 ) &
             + r         * ( - x000 + x100 ) &
             +     s     * ( - x000        + x010 ) &
             +         t * ( - x000               + x001 ) &
             + r * s     * ( + x000 - x100 - x010                      + x110 ) &
             + r     * t * ( + x000 - x100        - x001        + x101 ) &
             +     s * t * ( + x000        - x010 - x001 + x011 ) &
             + r * s * t * ( - x000 + x100 + x010 + x001 - x011 - x101 - x110 + x111 )
        
        return
      end subroutine blend_103

!#############################################################################

      subroutine blend_103_jacob ( r, s, t, x000, x001, x010, x011, x100, x101, x110, &
           x111, dxdr, dxds, dxdt )
        
        real(single) :: r
        real(single) :: s
        real(single) :: t
        real(single) :: dxdr
        real(single) :: dxds
        real(single) :: dxdt
        real(single) :: x000
        real(single) :: x001
        real(single) :: x010
        real(single) :: x011
        real(single) :: x100
        real(single) :: x101
        real(single) :: x110
        real(single) :: x111
        
        dxdr =  - x000 + x100  &
             + s     * ( + x000 - x100 - x010                      + x110 ) &
             + t * ( + x000 - x100        - x001        + x101 ) &
             + s * t * ( - x000 + x100 + x010 + x001 - x011 - x101 - x110 + x111 )
        
        
        dxds =   - x000        + x010  &
             + r *  ( + x000 - x100 - x010                      + x110 ) &
             + t * ( + x000        - x010 - x001 + x011 ) &
             + r * t * ( - x000 + x100 + x010 + x001 - x011 - x101 - x110 + x111 )
        
        dxdt =    - x000               + x001  &
             + r * ( + x000 - x100        - x001        + x101 ) &
             +    s  * ( + x000        - x010 - x001 + x011 ) &
             + r * s * ( - x000 + x100 + x010 + x001 - x011 - x101 - x110 + x111 )
        
        return
      end subroutine blend_103_jacob
      
      
!#############################################################################

      subroutine find_random_index(x, x_int, i1, i2, r)
        REAL(single), INTENT(IN)  :: x(:)
        REAL(single), INTENT(IN)  :: x_int
        REAL(single), INTENT(out) :: r
        INTEGER , INTENT(OUT) :: i1, i2
        INTEGER :: k, n
        integer,parameter :: order = 1
        n = SIZE(x)
        DO k=1,n
           IF (x_int <= x(k) ) EXIT
        END DO
        i1 = MIN(MAX(1,k-1-(ORDER/2)),n-ORDER)
        i2 = i1 + ORDER
        r  = (x_int - x(i1))/(x(i2)-x(i1)) 
      END SUBROUTINE find_random_index
      
!#############################################################################
      subroutine find_binary_index(x, x_int, index)
        REAL(double), INTENT(IN)  :: x(:)
        REAL(double), INTENT(IN)  :: x_int
        INTEGER , INTENT(OUT)     :: index
      
        integer                   :: i1, i2
        integer                   :: mid
        real                      :: tmp
        real(double)              :: diff

        i1 = 1
        i2 = size(x)
        index = 0
        do while (i1 .lt. i2)
           tmp = (i1 + i2)/2
           mid = NINT(tmp)    
           diff = x(mid)-x_int
           if (diff .eq. 0) then
              index = mid 
              exit
           elseif (diff .lt.0 )then
              i1 = i1+1
           else        
              i2 = mid-1
           end if
        end do
        
        if(index .eq.0) then
           if (i1 .eq. size(x)) then
              index = i1
           else if(abs(x_int-x(i1)) .gt. abs(x_int-x(i1+1))) then
              index = i1+1
           else
              index = i1  
           end if
        end if
              
        
      end subroutine find_binary_index
 
!#############################################################################
      subroutine linear_order_interp(x1, y1, x2, y2,n1,n2)
        REAL(double), INTENT(IN)  :: x1(n1),x2(n2),y1(n1)
        REAL(single), INTENT(out) :: y2(n2)
        INTEGER , INTENT(in) :: n1, n2
        INTEGER :: l, n
!!$ x1, x2 ,y1 must be arragned in ascending order

        i2=1
        DO l=1,n2
           if(x2(l) .lt. x1(1)) THEN
              i1 =1 
              i2 =2
           else if(x2(l) .gt. x1(n1)) then
              i1 = n1-1
              i2 = n1
           else
              Do n = i2,n1
                 IF (x1(n) >= x2(l) ) then
                    i1 = n-1
                    i2 = n
                    EXIT
                 end IF
              end Do
           end if
           y2(l)  = y1(i1)+(x2(l) - x1(i1))/(x1(i2)-x1(i1))*(y1(i2)-y1(i1))
        END DO

      END SUBROUTINE linear_order_interp
!#############################################################################
      
    end Module PCRTM_math_utility
