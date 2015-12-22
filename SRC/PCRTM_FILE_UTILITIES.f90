MODULE PCRTM_FILE_UTILITIES

CONTAINS
!#############################################################################
!  SUBROUTINE GETLUN(LUN)
!
! GET A FREE LOGICAL UNIT FOR FILE ACCESS
! !#############################################################################
  SUBROUTINE GETLUN(LUN)
    INTEGER, INTENT(OUT) :: LUN
    INTEGER              :: I
    LOGICAL              :: OPENUNIT
    
    DO I = 10,99
       INQUIRE(I,OPENED = OPENUNIT)
       IF( .NOT. OPENUNIT) EXIT
    ENDDO
    
    LUN = I
    
  END SUBROUTINE GETLUN
  
  
  
END MODULE PCRTM_FILE_UTILITIES
