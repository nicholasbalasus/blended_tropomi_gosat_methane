!! cakeCut_m
!! Credit:
!! Algorithm and implementation by Kai Yang (kaiyang@umd.edu)
!! April 2, 1998
!!
!! Usage: 
!! A polygon is represented by a list of vertices - (x,y) pairs, stored
!! in a two dimensions array. The parameter VMAX determines the 
!! maximum number of vertices that can be stored in
!! a polygon, and may be increased or decreased.
!!
!! The vertices stored in the Polygon should be in the clockwise
!! order when looking down at the Polygon. 
!! 
!! HcakCut performs horizontal cut, namely cutting along the x direction
!! at the integer grid lines (constant integer y).
!! Given an input polygon Pin, this function cuts Pin
!! into a number of output polygons Pout. 
!! Example:
!! npouts = HcakeCut( Pin, Pout)
!! 
!! VcakCut performs vertical cut, namely cutting along the y direction
!! at the integer grid lines (constant integer x).
!! Given an input polygon Pin, this function cuts Pin 
!! into a number of output polygons Pout.  
!! Example:
!! npouts = VcakeCut( Pin, Pout)
!!
!! area2
!! This function returns the 2 times the area of a polygon.
!!

MODULE cakeCut_m
    IMPLICIT NONE
    INTEGER, PARAMETER :: VMAX = 10 
    TYPE Polygon
      INTEGER :: nv 
      REAL(KIND=8), DIMENSION(VMAX,2) :: vList
    END TYPE Polygon
    
    INTERFACE  append
     MODULE PROCEDURE append1
     MODULE PROCEDURE append2
    END INTERFACE

    PUBLIC :: HcakeCut
    PUBLIC :: VcakeCut
    PUBLIC :: area2

    CONTAINS  

      SUBROUTINE append1( p, x, y )
        TYPE( Polygon ), INTENT(INOUT) :: p
        REAL(KIND=8), INTENT(IN) :: x, y 
        INTEGER :: iv
        IF( p%nv < 0 .OR. p%nv >= VMAX ) THEN
           WRITE(*,*) 'number of vertices exceeds maximum.'
           STOP
        ENDIF

        iv = p%nv + 1;
        p%vList(iv,1) = x
        p%vList(iv,2) = y
        p%nv = iv
        RETURN 
      END SUBROUTINE append1

      SUBROUTINE append2( p, v )
        TYPE( Polygon ), INTENT(INOUT) :: p
        REAL(KIND=8), DIMENSION(2), INTENT(IN) :: v
        INTEGER :: iv
        IF( p%nv < 0 .OR. p%nv >= VMAX ) THEN
           WRITE(*,*) 'number of vertices exceeds maximum.'
           STOP
        ENDIF

        iv = p%nv + 1;
        p%vList(iv,1:2) = v(1:2)
        p%nv = iv
        RETURN 
      END SUBROUTINE append2

      FUNCTION HcakeCut( Pin, Pouts ) RESULT (npouts)
        TYPE( Polygon ), INTENT(IN) :: Pin
        TYPE( Polygon ), DIMENSION(:), INTENT(OUT) :: Pouts
        INTEGER :: npouts 
        !! local variable
        REAL(KIND=8), DIMENSION(Pin%nv,2) :: vList 
        INTEGER :: nv, jmin, jmax, i, irowMax, irowMin
        REAL(KIND=8) ::  ymin, ymax, x, y
        INTEGER :: sign, ip, j, j1 
       
        nv = Pin%nv
        vList(1:nv,1:2) = Pin%vList(1:nv,1:2)
        IF( nv < 3 ) THEN
           npouts = 0
           RETURN
        ENDIF

        !! find ymin, ymax and their corresponding indices 

        jmin = 1;  jmax = 1
        ymin = vList(jmin, 2); ymax = vList(jmax, 2)

        DO i = 2, nv
          IF( vList(i, 2) < ymin ) THEN
             jmin = i; ymin = vList(jmin, 2)
          ENDIF
          IF( vList(i, 2) > ymax ) THEN
             jmax = i; ymax = vList(jmax, 2)
          ENDIF
        ENDDO
   
        irowMax = CEILING( ymax ); irowMin = FLoOR( ymin ) 
        npouts = irowMax - irowMin;
 
        IF( npouts <= 0 ) RETURN
      
        IF( npouts == 1 ) THEN
           Pouts(1)%nv = nv 
           Pouts(1)%vList(1:nv,1:2) = vList(1:nv,1:2)
           RETURN
        ENDIF
 
        DO i = 1, npouts
          Pouts(i)%nv = 0
          Pouts(i)%vList = 0.d0
        ENDDO

        i  = 0; ip = 1; sign = 1
        j  = jmin 
        j1 = MOD( j, nv) + 1
        y  = FLOOR( ymin ) + 1

        DO WHILE( i < nv ) 
          IF( sign*vList(j1, 2) > sign* y  ) THEN !true, for intersecting y 
             x = vList( j,1 ) + (y-vList(j,2)) &
               *(vList(j1,1 ) -  vList(j,1))/(vList(j1,2) - vList(j,2))  
             CALL append( Pouts(ip), x, y )
             ip = ip + sign;
             CALL append( Pouts(ip), x, y )
             y  =  y + DBLE(sign)  
          ELSE
             CALL append( Pouts(ip), vList(j1,1:2) )
             IF( sign*vList(j1, 2) >= sign* y  ) THEN !! > not possible, for == 
                IF( 1 <= ip + sign .AND. ip + sign <= npouts ) THEN 
                   ip = ip + sign 
                   y  =  y + DBLE(sign) 
                   CALL append( Pouts(ip), vList(j1,1:2) )
                ENDIF
             ENDIF
             j = j1; j1 = MOD(j, nv) + 1
             IF( j == jmax ) THEN
                sign = -1; y = y + DBLE(sign)
             ENDIF
             i = i + 1
          ENDIF
        END DO
        RETURN
      END FUNCTION HcakeCut

      FUNCTION VcakeCut( Pin, Pouts, pindex, id_sub ) RESULT (npouts)
        TYPE( Polygon ), INTENT(IN) :: Pin
        TYPE( Polygon ), DIMENSION(:), INTENT(OUT) :: Pouts
        INTEGER :: pindex
        INTEGER :: id_sub
        INTEGER :: npouts 

        !! Local
        REAL(KIND=8), DIMENSION(Pin%nv,2) :: vList 
        INTEGER :: nv, jmin, jmax, i, icolMax, icolMin
        REAL(KIND=8) ::  xmin, xmax, x, y
        INTEGER :: sign, ip, j, j1 
       
        nv = Pin%nv
        vList(1:nv,1:2) = Pin%vList(1:nv,1:2)
        IF( nv < 3 ) THEN
           npouts = 0
           RETURN
        ENDIF

       ! find ymin, ymax and their corresponding indices 
       jmin = 1; jmax = 1
       xmin = vList(jmin, 1)
       xmax = vList(jmax, 1)

       DO i = 2, nv
         IF( vList(i, 1) < xmin ) THEN
            jmin = i; xmin = vList(jmin, 1)
         ENDIF
         IF( vList(i, 1) > xmax ) THEN
            jmax = i; xmax = vList(jmax, 1)
         ENDIF 
       ENDDO

       icolMax = CEILING( xmax ); icolMin = FLOOR( xmin )
       npouts = icolMax - icolMin

        IF( npouts <= 0 ) RETURN
      
        IF( npouts == 1 ) THEN
           Pouts(1)%nv = nv 
           Pouts(1)%vList(1:nv,1:2) = vList(1:nv,1:2)
           RETURN
        ENDIF
 
        DO i = 1, npouts
          Pouts(i)%nv = 0
          !IF( (pindex.GE.1792284) .AND. (id_sub==5) .AND. (i.GE.520) ) THEN
          !   PRINT*, i, npouts
          !   PRINT*, "---------------"
          !   PRINT*, SHAPE(Pouts)
          !   PRINT*, "---------------"
          !   PRINT*, Pouts(i)%nv
          !   PRINT*, "---------------"
          !   PRINT*, Pouts(i)%vList
          !   PRINT*, vList
          !ENDIF
          Pouts(i)%vList = 0.d0
        ENDDO
 
        i  = 0; ip = 1; sign = 1
        j  = jmin
        j1 = MOD( j, nv) + 1
        x  = FLOOR( xmin ) + 1;

        DO WHILE( i < nv )
          IF( sign*vList(j1, 1) > sign* x ) THEN  !! true, for intersecting x 
             y = vList( j,2 ) + (x-vList(j,1)) &
               *(vList(j1,2 ) - vList(j,2))/(vList(j1,1) - vList(j,1))
             CALL append( Pouts(ip), x, y )
             ip = ip + sign;
             CALL append( Pouts(ip), x, y ) 
             x  =  x + DBLE(sign)
          ELSE
             CALL append( Pouts(ip), vList(j1,1:2) )
             IF( sign*vList(j1, 1) >= sign* x  ) THEN   !! > not possible, for == 
                IF( 1 <= ip + sign .AND. ip + sign <= npouts ) THEN
                   ip =  ip + sign 
                   x  =  x  + DBLE(sign)
                   CALL append( Pouts(ip), vList(j1,1:2) )
                ENDIF
             ENDIF
             j = j1
             j1 = MOD(j, nv) + 1
             IF( j == jmax ) THEN
                sign = -1; x =  x + DBLE(sign)
             ENDIF
             i = i + 1 
          ENDIF       
        ENDDO
        RETURN
      END FUNCTION VcakeCut

      FUNCTION area2( p ) 
        TYPE( Polygon ), INTENT(IN ) :: p
        REAL(KIND=8) ::  area2  

        INTEGER :: nv,  i,i1 

        area2 = 0.d0 
        nv = p%nv
        IF( nv < 3 ) RETURN 
        
        DO i = 2, nv-1
          i1 = i+1 
          area2 = area2 + (p%vList(i1,1)-p%vList(1,1))*(p%vList(i,2)-p%vList(1,2)) &
                        - (p%vList(i1,2)-p%vList(1,2))*(p%vList(i,1)-p%vList(1,1)) 
        ENDDO
        RETURN
      END FUNCTION area2
END MODULE cakeCut_m
