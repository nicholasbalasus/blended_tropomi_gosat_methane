PROGRAM RegridPixels

        !==============================================================================
        ! Cakecut module is from Kai Yang (kaiyang@umd.edu)
        ! Applied by Lei Zhu (leizhu@fas.harvard.edu) to regrid pixles
        ! 01/15/15

        ! Now save Pixels_count in the output
        ! lei, 10/27/16

        ! Adapted for SCIA oversampling - JDM 10/28/2016

        !==============================================================================
        ! How to use:
        ! Compile: gfortran -o RegridPixels.x cakecut_m.f90 tools_m.f90 RegridPixels.f90
        ! Change user inputs in run.sh
        ! Run: ./run.sh

        !==============================================================================
        ! Modules
        USE cakecut_m
        USE tools_m

        IMPLICIT NONE

        !==============================================================================
        ! Variables

        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! Input and output dirs
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! Input and output file dirs
        CHARACTER*200 Output_Dir
        ! Input and output file names
        CHARACTER*200 Input_Filename, Output_Filename
        ! Input file unit
        INTEGER*4 :: In_Unit

        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! Satellite pixels information read from the input file
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! Total number of pixels in the input file
        INTEGER*4 :: NPIXELS
        ! Lat and Lon read from the input file
        ! 4 corners and the center
        REAL,    ALLOCATABLE :: LAT(:,:), LON(:,:)
        ! Columns read from the input
        REAL,    ALLOCATABLE :: VCD(:), VCD_UNC(:)
        !REAL,    ALLOCATABLE :: tropday(:), tropqa(:)
        ! Other information read from the input file, discarded
        REAL    :: temp_ID 
        !, temp_ID2, temp_ID3
        !REAL    :: tempval1
        !, tempval2, tempval3, tempval4
        !REAL    :: tempval5, tempval6, tempval7, tempval8
        !REAL    :: tempval9, tempval10, tempval11
	! We need temp_ID to be a real since Matlab (JDM)

        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! Single pixel
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! The area of the pixel
        REAL    :: A
        ! Lat and Lon of 4 corners read from the input 
        REAL    :: Lat_o(4), Lon_o(4)
        ! Sorted Lat and Lon of the 4 corners, so that it's in clockwise order
        INTEGER :: id(4)
        REAL    :: Lat_r(4), Lon_r(4)

        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! Related to cakecut scripts
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! Single pixel
        TYPE( POLYGON )                 :: pixel 
        ! Results from Horizontal and Veritical cut 
        TYPE( POLYGON ), DIMENSION(500) :: sub_pixels, final_pixels 
        INTEGER :: id_all, id_sub, id_final, n_sub_pixels, n_final_pixels
        ! Overlapped area
        REAL    :: temp_area
        ! Total area of the pixel, will be compared with A
        ! to test area conservativeness
        REAL    :: pixel_area

        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! Domain parameters
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! 4 edges of the domain
        REAL    :: Lat_low, Lat_up, Lon_left, Lon_right
        ! Resolution, use 0.02 
        REAL    :: Res
        ! # of rows (lat) and cols (lon) of the domain
        INTEGER :: NROWS, NCOLS
        ! Lat and Lon of the out grids
        REAL,    ALLOCATABLE :: Lat_grid(:), Lon_grid(:) 

        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! Results
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! Overlapped area in each grid
        REAL,    ALLOCATABLE :: Sub_Area(:,:)
        ! Sum at each step (pixel)
        REAL,    ALLOCATABLE :: Sum_Above(:,:), Sum_Below(:,:)
        ! Count, how many pixels have contribution to this cell
        INTEGER, ALLOCATABLE :: Pixels_count(:,:)
        ! Overlapped area and error weighted final results
        REAL,    ALLOCATABLE :: Average(:,:)

        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! Loop index and other
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++
        INTEGER*4 :: i, row, col, p

        !==============================================================================
        !==============================================================================
        !==============================================================================
        ! User input, read from run.sh 
        !==============================================================================
        !==============================================================================
        !==============================================================================
        
        READ ( 5, '(A)' ) Output_Dir
        READ ( 5, '(A)' ) Input_Filename
        READ ( 5, '(A)' ) Output_Filename
        READ ( 5, *     ) Res

        !==============================================================================
        !==============================================================================
        !==============================================================================
        ! Read satellite pixels
        !==============================================================================
        !==============================================================================
        !==============================================================================

        ! First, get the number of total pixels

        In_Unit = 100
        p = 0

        ! Open pixels file
        CALL Check_File(Input_Filename, In_Unit)
11      CONTINUE
        p = p + 1
        READ(In_Unit, *, END=12, ERR=12)
        GOTO 11
12      CONTINUE

        NPIXELS = p - 1

        PRINT*, " ----------------------------------------------------------"
        PRINT*, " ----------------------------------------------------------"
        PRINT*, "  Total pixels : ", NPIXELS

        ! Allocate arraies, change to NPIXELS instead of NPIXELS+1 (NB)
        ALLOCATE ( LAT     (NPIXELS, 5) )
        ALLOCATE ( LON     (NPIXELS, 5) )
        ALLOCATE ( VCD     (NPIXELS) )
        ALLOCATE ( VCD_UNC (NPIXELS) )
        !ALLOCATE ( tropday     (NPIXELS+1) )
        !ALLOCATE ( tropqa     (NPIXELS+1) )

        ! Re-read the file, and store data fields
        p = 0
        REWIND(In_Unit)
13      CONTINUE
        p = p + 1
        ! Input file format
        ! temp_ID: pixel ID, not used in the code
        ! LAT(p,1:5): lat of the 4 corners (1:4) and center (5) of pixel p
        ! LON(p,1:5): lon of the 4 corners (1:4) and center (5) of pixel p
        ! temp_SZA, temp_CF, temp_AMF: not used in the code
        ! VCD(p): Vertical column density of pixel p
        ! VCD_Unc(p): Error associated with VCD(p)
	! Free format reading (JDM)
        !names = [,lat0,lat1,lat2,lat3,lat,lon0,lon1,lon2,lon3,lon,xch4,xch4_unc,qa,weekday,
        !year,day,month,albedo_swir,aerosol_op,aerosol_num,ch4col,surfpres],
        !READ(In_Unit,"(I8,13F15.6,2E15.6)",end=15,err=14) temp_ID,                          &
        !READ(In_Unit,*,end=15,err=14) temp_ID,                         &
        !READ(In_Unit,"(I8,12F16.6)",end=15,err=14) temp_ID, &
        READ(In_Unit,*,end=15,err=14) temp_ID, &
             LAT(p,1:5), LON(p,1:5), VCD(p), VCD_UNC(p)
                !tempval1, &
                !tempval2, tempval3, tempval4, tempval5, tempval6,      &
                !tempval7, tempval8, tempval9, tempval10, &
                !tempval11

        GOTO 13
14      CONTINUE
        PRINT*, "  Error in reading: ", TRIM(Input_Filename)
        PRINT*, "  At line: ", p
        READ(In_Unit,"(I8,10F15.6,2E15.6)") temp_ID, LAT(p,1:5), &
             LON(p,1:5), VCD(p), VCD_UNC(p)
        PRINT*, "  Check the format"
        STOP
15      CONTINUE
        CLOSE(In_Unit)

        !==============================================================================
        !==============================================================================
        !==============================================================================
        ! Get domain parameters based on the input information
        !==============================================================================
        !==============================================================================
        !==============================================================================

        ! Get 4 edges of the domain
        Lat_low   = FLOOR(MINVAL(LAT))-7
        Lat_up    = CEILING(MAXVAL(LAT))+7
        Lon_left  = FLOOR(MINVAL(LON))-7
        Lon_right = CEILING(MAXVAL(LON))+7

        ! Get # rows and cols of the domain
        NROWS = ( Lat_up    - Lat_low  ) / Res
        NCOLS = ( Lon_right - Lon_left ) / Res

        ! Get the output grid coordinate, center of grid
        ALLOCATE ( Lat_grid (NROWS) )
        ALLOCATE ( Lon_grid (NCOLS) )
        DO row = 1,NROWS
          Lat_grid(row) = Lat_low  + row*Res - 0.5*Res
        ENDDO
        DO col = 1,NCOLS
          Lon_grid(col) = Lon_left + col*Res - 0.5*Res
        ENDDO

        ! Print domain parameters
        PRINT*, " ----------------------------------------------------------"
        PRINT*, "  Domain parameters"        
        PRINT*, "  Lat range    : ", Lat_low,  Lat_up
        PRINT*, "  Lon range    : ", Lon_left, Lon_right
        PRINT*, "  Resolution   : ", Res
        PRINT*, "  Rows         : ", NROWS
        PRINT*, "  Cols         : ", NCOLS

        ! Allocate arraies and initializ them
        ALLOCATE ( Sub_Area(NROWS,     NCOLS) )
        ALLOCATE ( Sum_Above(NROWS,    NCOLS) )
        ALLOCATE ( Sum_Below(NROWS,    NCOLS) )
        ALLOCATE ( Pixels_count(NROWS, NCOLS) )
        ALLOCATE ( Average(NROWS,      NCOLS) )

        Sub_Area     = 0.0D0
        Sum_Above    = 0.0D0
        Sum_Below    = 0.0D0
        Pixels_count = 0
        Average      = 0.0D0

        !==============================================================================
        !==============================================================================
        !==============================================================================
        ! Loop pixels, get overlapped area using Cakecut functions
        !==============================================================================
        !==============================================================================

        PRINT*, " ----------------------------------------------------------"

        DO p = 1, NPIXELS

        ! IF( p.GE.2838660 ) THEN
        !    PRINT*, "First: ", p
        ! ENDIF        

        ! Print process
        IF(MOD(p, 100000)==0) THEN
          PRINT*,"  Finish (%): ", NINT(REAL(p)/NPIXELS*100)
        ENDIF

        Lat_r = 0.0D0
        Lon_r = 0.0D0

        ! Get 4 corners from the satellite
        Lon_o = LON(p,1:4)
        Lat_o = LAT(p,1:4)

        IF( MAXVAL(ABS(Lon_o)) .GT. 178 ) CYCLE
        IF( MAXVAL(ABS(Lat_o)) .GT. 88 ) CYCLE

        ! IF( p.GE.2838659 ) THEN
        !IF(p.GE.1792284) THEN
        !   PRINT*, "o: ", Lon_o, Lat_o
        !   PRINT*, "r: ", Lon_r, Lat_r
        !   PRINT*, id
        !ENDIF

        ! Weekday selector
        !IF( tropday(p) .GT. 4) CYCLE

        ! Require quality flag > 0
        !IF( tropqa(p) .LT. 0.3) CYCLE
       
        ! Arrange the corners so that it's in colockwise order
        id = -9999
        CALL Clockwise_Sort( Lon_o, Lat_o, id )

        ! Check for bad pixels
        IF( SUM(id) /= 10 ) THEN
          PRINT*, "Bad pixel: ", p
          PRINT*, Lon_o, Lat_o
          PRINT*, Lon_r, Lat_r
          PRINT*, id
          STOP
        ENDIF

        ! Get the colockwise-ordered 4 corners
        DO i = 1, 4
          Lon_r(i) = Lon_o(id(i))
          Lat_r(i) = Lat_o(id(i))
        ENDDO

        ! Arrange the pixel corners in the grids
        pixel%nv = 4

        DO i = 1, 4
          pixel%vList(i,1:2) = (/ (Lon_r(i)-Lon_left)/Res, (Lat_r(i)-Lat_low)/Res /) 
        ENDDO

        ! Get the area of the pixel
        CALL Cal_Quard_Area(Lon_r,Lat_r,A)

        ! IF( p.GE.2838659 ) THEN
        !    PRINT*, "After Cal_Quard_Area: ", p
        ! ENDIF

        ! Used for the next step
        id_all = 0
        pixel_area = 0.0d0

        ! Perform Horizontal cut first at the integer grid lines
        n_sub_pixels = HcakeCut( pixel, sub_pixels)
        ! IF( p.GE.2838659 ) THEN
        !   PRINT*, "After HcakeCut: ", n_sub_pixels
        ! ENDIF        

        ! Then perform Vertical cut for each sub pixel obtainted 
        ! from the Horizontal cut at the integer grid lines
        DO id_sub = 1, n_sub_pixels
          !IF( (p.GE.1792284) .AND. (id_sub==5) ) THEN
          !   PRINT*, sub_pixels(id_sub)
          !   PRINT*, "-----"
          !   PRINT*, final_pixels
          !   PRINT*, "-----"
          !ENDIF
          n_final_pixels = VcakeCut( sub_pixels(id_sub), final_pixels, p, id_sub )
          !IF( p.GE.2838659 ) THEN
          !     PRINT*, "In n_sub_pixels DO: ", id_sub, n_sub_pixels
          !ENDIF
          DO id_final = 1, n_final_pixels
            id_all = id_all + 1
            temp_area = area2(final_pixels(id_final))*0.5d0*Res*Res
            pixel_area = pixel_area + temp_area
            !DO i = 1, final_pixels(id_final)%nv
            !  WRITE(*,'(A,F10.3,A,F10.3,A)') '{', final_pixels(id_final)%vList(i,1),',', &
            !                                      final_pixels(id_final)%vList(i,2),'},' 
            !ENDDO
            row = FLOOR(MINVAL(final_pixels(id_final)%vList(1:final_pixels(id_final)%nv,2))) + 1
            col = FLOOR(MINVAL(final_pixels(id_final)%vList(1:final_pixels(id_final)%nv,1))) + 1
            !WRITE(*,*) 'id =', id_all, ', aera=',temp_area, ', row=', row, ',col=',col

            ! Get the overlaped area between the pixel and each cell
            Sub_Area(row,col) = temp_area

            ! Sum weighted value and weights
            ! Here use temp_area/A/VCD_Unc(p) as averaging weight, meaning that
            ! averaging weight is assumed to be proportional to the ratio of the overlap area (temp_area) to the
            ! pixel size (A) and inversely proportional to the error standard deviation (VCD_Unc(p)).
            ! If you just want fraction of overlap area as averaging weight, use: temp_area/A
            ! If you just want area weighted average, use: temp_area
            ! Use VCD_Unc (NB)
            Sum_Above(row,col) = Sum_Above(row,col) + temp_area/A/VCD_Unc(p)*VCD(p)
            Sum_Below(row,col) = Sum_Below(row,col) + temp_area/A/VCD_Unc(p)
            ! Sum_Above(row,col) = Sum_Above(row,col) + temp_area/A*VCD(p)
            ! Sum_Below(row,col) = Sum_Below(row,col) + temp_area/A
            Pixels_count(row,col) = Pixels_count(row,col) + 1
          ENDDO
        ENDDO
        
        ! Check area consvertiveness
        !IF( ABS(A-pixel_area)/A >=0.05  ) THEN
        !  PRINT*, "------------------------------------------------------"
        !  PRINT*, "  -Area not conservative at pixel:", p, A, pixel_area    
        !  PRINT*, Lon_o, Lat_o
        !  PRINT*, Lon_r, Lat_r
        !ENDIF
        
        ! End loop pixels
        ENDDO 

        !IF( p.GE.2838660 ) THEN
        !   PRINT*, "After pixel DO: ", p
        !ENDIF

        !==============================================================================
        !==============================================================================
        !==============================================================================
        ! Get weigthed average and save output
        !==============================================================================
        !==============================================================================
        !==============================================================================

        ! Open output file
        Output_Filename = TRIM(Output_Dir)//TRIM(Output_Filename)
        OPEN(99,file=TRIM(Output_Filename))

        ! Save
        DO row = 1, NROWS
        DO col = 1, NCOLS
          ! Change to just disqualifying Sum_Below==0 since I am regridding surface altitude with a lot of 0s (NB)
          ! IF( Sum_Above(row,col)/=0 .AND. Sum_Below(row,col)/=0  ) THEN
          IF( Sum_Below(row,col)/=0  ) THEN
            Average(row,col)= Sum_Above(row,col)/Sum_Below(row,col)     
            WRITE(99,'(2I6,2f12.4,E15.6,I6)') row,col,Lat_grid(row),Lon_grid(col),      &
                                              Average(row,col),Pixels_count(row,col)
          ENDIF
        ENDDO
        ENDDO

        ! Close output
        CLOSE(99)

        PRINT*, " ----------------------------------------------------------"
        PRINT*, "  Output saved at: ", TRIM(Output_Filename)

        !==============================================================================
        !==============================================================================
        !==============================================================================
        ! Free memories
        !==============================================================================
        !==============================================================================
        !==============================================================================

        IF ( ALLOCATED( LAT          ) ) DEALLOCATE( LAT          )
        IF ( ALLOCATED( LON          ) ) DEALLOCATE( LON          )
        IF ( ALLOCATED( VCD          ) ) DEALLOCATE( VCD          )
        IF ( ALLOCATED( VCD_UNC      ) ) DEALLOCATE( VCD_UNC      )
        IF ( ALLOCATED( Lat_grid     ) ) DEALLOCATE( Lat_grid     )
        IF ( ALLOCATED( Lon_grid     ) ) DEALLOCATE( Lon_grid     )
        IF ( ALLOCATED( Sub_Area     ) ) DEALLOCATE( Sub_Area     )
        IF ( ALLOCATED( Sum_Above    ) ) DEALLOCATE( Sum_Above    )
        IF ( ALLOCATED( Sum_Below    ) ) DEALLOCATE( Sum_Below    )
        IF ( ALLOCATED( Pixels_count ) ) DEALLOCATE( Pixels_count )
        IF ( ALLOCATED( Average      ) ) DEALLOCATE( Average      )

END PROGRAM RegridPixels
