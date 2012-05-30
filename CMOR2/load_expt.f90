!!$program driver
!!$  use exp_info
!!$  implicit none
!!$  integer::i
!!$  character(len=256)::exp_file
!!$  !
!!$  exp(:)%case(1:)          = ' '
!!$  exp(:)%model_id(1:)      = ' '
!!$  exp(:)%loc(1:)           = ' '
!!$  exp(:)%expt_id(1:)       = ' '
!!$  exp(:)%rip_code(1:)      = ' '
!!$  exp(:)%cmip(1:)          = ' '
!!$  exp(:)%run_refcase(1:)   = ' '
!!$  exp(:)%run_refdate(1:)   = ' '
!!$  exp(:)%begin_end(1:)     = ' '
!!$  exp(:)%grid(1:)          = ' '
!!$  exp(:)%compset(1:)       = ' '
!!$  exp(:)%repotag(1:)       = ' '
!!$  exp(:)%start_fin(1:)     = ' '
!!$  exp(:)%mach(1:)          = ' '
!!$  exp(:)%dout(1:)          = ' '
!!$  exp(:)%forcing(1:)       = ' '
!!$  exp_file = 'experiments.txt'
!!$  !
!!$  call load_exp(exp_file)
!!$  !
!!$  do i = 1,num_exp
!!$!     write(*,'(a,$)') '<tr>'
!!$!     write(*,'(20(''<td>'',a,''</td>''),$)') &
!!$     write(*,'(20(a,'',''))') &
!!$          trim(adjustl(exp(i)%case)),&
!!$          trim(adjustl(exp(i)%model_id)),&
!!$          trim(adjustl(exp(i)%loc)),&
!!$          trim(adjustl(exp(i)%expt_id)),&
!!$          trim(adjustl(exp(i)%rip_code)),&
!!$          trim(adjustl(exp(i)%cmip)),&
!!$          trim(adjustl(exp(i)%run_refcase)),&
!!$          trim(adjustl(exp(i)%run_refdate)),&
!!$          trim(adjustl(exp(i)%begin_end)),&
!!$          trim(adjustl(exp(i)%grid)),&
!!$          trim(adjustl(exp(i)%compset)),&
!!$          trim(adjustl(exp(i)%repotag)),&
!!$          trim(adjustl(exp(i)%start_fin)),&
!!$          trim(adjustl(exp(i)%mach)),&
!!$          trim(adjustl(exp(i)%dout)),&
!!$          trim(adjustl(exp(i)%forcing)),&
!!$!     write(*,*) '</tr>'
!!$  enddo
!!$end program driver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine load_exp(exp_file)
  use exp_info
  implicit none
  !
  integer::iostat,i,j,iexp
  character(len=256)::exp_file
  character(len=384)::instring
  logical::does_exist
  !
  ! Get experiment information
  !
  ! Master CCSM experiment list
  !
  !  Columns   Field
  !   1 -  54  CCSM case name
  !  55 -  74  model_id (CCSM4, CESM-CAM5, CESM-BGC, CESM-CHEM, CESM-WACCM, CESM1-CAM5.1-FV2, etc)
  !  75 -  76  Location (NC = NCAR; NE = NERSC; OR = ORNL)
  !  80 -  94  Official MIP name, or very brief description (N/A if not applicable)
  !  95 - 104  RIP code (cm5) or realization number (cm3) (N/A if not applicable)
  ! 105 - 109  MIP (cm3 or cm5) experiment (N/A if not applicable)
  ! 110 - 149  RUN_REFCASE (parent case)
  ! 150 - 164  RUN_REFDATE (branch date, yyyy-mm-dd)
  ! 165 - 174  years of experiment (YYYY if unknown)
  ! 175 - 184  GRID (resolution)
  ! 185 - 204  COMPSET (N/A if not applicable)
  ! 205 - 224  REPOTAG (N/A if not applicable)
  ! 225 - 244  Calendar dates of simulation execution (yyyy/mm-yyyy/mm)
  ! 245 - 254  MACH (hardware)
  ! 255 - end  DOUT_L_MSROOT (history file location on archive)
  !
  !
  inquire(file=trim(exp_file),exist=does_exist)
  if (.not.(does_exist)) then
     write(*,*) 'Cannot find ',trim(exp_file),'. Dying.'
     stop
  endif
  !
  exp(:)%case(1:)          = ' '
  exp(:)%model_id(1:)      = ' '
  exp(:)%loc(1:)           = ' '
  exp(:)%expt_id(1:)       = ' '
  exp(:)%rip_code(1:)      = ' '
  exp(:)%cmip(1:)          = ' '
  exp(:)%run_refcase(1:)   = ' '
  exp(:)%run_refdate(1:)   = ' '
  exp(:)%begin_end(1:)     = ' '
  exp(:)%grid(1:)          = ' '
  exp(:)%compset(1:)       = ' '
  exp(:)%repotag(1:)       = ' '
  exp(:)%start_fin(1:)     = ' '
  exp(:)%mach(1:)          = ' '
  exp(:)%dout(1:)          = ' '
  exp(:)%forcing(1:)       = 'unknown forcings'
  !
  open(20,file=trim(adjustl(exp_file)),form='formatted')
  iostat = 0 ; iexp = 1
  do while (iostat == 0)
     instring(1:)  = ' '
     read(20,'(a)',iostat=iostat) instring
     if (iostat == 0) then
        if (instring(1:1) /= '#') then
           exp(iexp)%case(1:)        = adjustl(instring(  1: 54))
           exp(iexp)%model_id(1:)    = adjustl(instring( 55: 74))
           exp(iexp)%loc(1:)         = adjustl(instring( 75: 76))
           exp(iexp)%expt_id(1:)     = adjustl(instring( 80: 94))
           exp(iexp)%rip_code(1:)    = adjustl(instring( 95:104))
           if (instring(105:109) == 'cm5') exp(iexp)%cmip(1:) = 'CMIP5'
           if (instring(105:109) == 'gmp') exp(iexp)%cmip(1:) = 'GeoMIP'
           if (instring(105:109) == 'tmp') exp(iexp)%cmip(1:) = 'TAMIP'
           exp(iexp)%run_refcase(1:) = adjustl(instring(110:149))
           exp(iexp)%run_refdate(1:) = adjustl(instring(150:164))
           exp(iexp)%begin_end(1:)   = adjustl(instring(165:174))
           exp(iexp)%grid(1:)        = adjustl(instring(175:184))
           exp(iexp)%compset(1:)     = adjustl(instring(185:204))
           exp(iexp)%repotag(1:)     = adjustl(instring(205:224))
           exp(iexp)%start_fin(1:)   = adjustl(instring(225:244))
           exp(iexp)%mach(1:)        = adjustl(instring(245:254))
           exp(iexp)%dout(1:)        = adjustl(instring(255:len_trim(instring)))
           iexp = iexp + 1
        endif
     endif
  enddo
  close(20)
  num_exp = iexp - 1
  !
  ! Parse out beginning/ending years
  do iexp = 1,num_exp
     if (exp(iexp)%expt_id(1:5) == 'tamip') then
        read(exp(iexp)%expt_id(6:9),'(i4.4)') exp(iexp)%begyr
        exp(iexp)%endyr = exp(iexp)%begyr
        exp(iexp)%length = (exp(iexp)%endyr-exp(iexp)%begyr)+1
!        write(*,'(''load_expt: '',a,5x,i4)') trim(exp(iexp)%expt_id),exp(iexp)%begyr
     else
        if (exp(iexp)%begin_end /= 'YYYY-YYYY') then
           read(exp(iexp)%begin_end(1:4),'(i4.4)') exp(iexp)%begyr
           read(exp(iexp)%begin_end(6:9),'(i4.4)') exp(iexp)%endyr
           exp(iexp)%length = (exp(iexp)%endyr-exp(iexp)%begyr)+1
        else
           write(*,'(''load_expt: YEARS UNDEFINED FOR '',a,'' STOPPING.'')') trim(exp(iexp)%expt_id)
           stop
        endif
     endif
  enddo
  ! Possible forcings:
  ! N/A Nat Ant GHG SD SI SA TO SO Oz LU Sl Vl SS Ds BC MD OC AA
  !     Description 
  ! AA  anthropogenic aerosols (a mixture of aerosols, not explicitly defined here) 
  ! Ant anthropogenic forcing  (a mixture, not explicitly defined here)
  ! BC  black carbon 
  ! Ds  Dust 
  ! GHG well-mixed greenhouse gases (a mixture, not explicitly defined here) 
  ! LU  land-use change 
  ! MD  mineral dust 
  ! Nat natural forcing (a combination, not explicitly defined here)
  ! OC  organic carbon 
  ! Oz  (= TO + SO) ozone (= tropospheric and stratospheric ozone) 
  ! SA  (= SD + SI) anthropogenic sulfate aerosol direct and indirect effects 
  ! SD  anthropogenic sulfate aerosol, accounting only for direct effects 
  ! SI  anthropogenic sulfate aerosol, accounting only for indirect effects 
  ! Sl  solar irradiance
  ! SO  stratospheric ozone
  ! SS  sea salt 
  ! TO  tropospheric ozone 
  ! Vl  volcanic aerosol
  !
  do i = 1,num_exp
     select case (exp(i)%expt_id)
     case ('1pctCO2')
        exp(i)%forcing(1:)      = 'GHG (CO2 only)'
     case ('abrupt4xCO2')
        exp(i)%forcing(1:)      = 'GHG (CO2 only, instantaneously to 4X)'
     case ('amip')
        exp(i)%forcing(1:)      = 'Sl GHG Vl SS Ds SA BC MD OC Oz AA'
     case ('amip4K')
        exp(i)%forcing(1:)      = 'Sl GHG Vl SS Ds SA BC MD OC Oz AA (4K SST perturbation)'
     case ('amip4xCO2')
        exp(i)%forcing(1:)      = 'GHG (CO2 only, instantaneously to 4X) Sl Vl SS Ds SA BC MD OC Oz AA'
     case ('amipFuture')
        exp(i)%forcing(1:)      = 'Sl GHG Vl SS Ds SA BC MD OC Oz AA (and SST perturbation)'
     case ('aqua4K')
        exp(i)%forcing(1:)      = 'unknown forcings'
     case ('aqua4xCO2')
        exp(i)%forcing(1:)      = 'unknown forcings'
     case ('aquaControl')
        exp(i)%forcing(1:)      = 'unknown forcings'
     case ('decadal1961','decadal1966','decadal1971','decadal1975','decadal1976','decadal1980',&
           'decadal1981','decadal1985','decadal1986','decadal1990','decadal1991','decadal1995','decadal1996')
        exp(i)%forcing(1:)      = 'Sl GHG Vl SS Ds SD BC MD OC Oz AA LU'
     case ('decadal2000','decadal2001','decadal2002','decadal2003','decadal2004','decadal2005','decadal2006')
        exp(i)%forcing(1:)      = 'Sl GHG SS Ds SD BC MD OC Oz AA LU'
     case ('esmControl')
        exp(i)%forcing(1:)      = 'Sl GHG SS Ds SD BC MD OC Oz AA (all fixed at 1850 values)'
     case ('esmFdbk1')
        exp(i)%forcing(1:)      = 'GHG (CO2 only, 1% increase/year, radiation only) Sl SS Ds SD BC MD OC Oz AA (all fixed at 1850 values)'
     case ('esmFdbk2')
        exp(i)%forcing(1:)      = 'GHG (CO2 only, historical + RCP4.5, radiation only) Sl SS Ds SD BC MD OC Oz AA (all fixed at 1850 values)'
     case ('esmFixClim1')
        exp(i)%forcing(1:)      = 'GHG (CO2 only, 1% increase/year, carbon cycle only) Sl SS Ds SD BC MD OC Oz AA (all fixed at 1850 values)'
     case ('esmFixClim2')
        exp(i)%forcing(1:)      = 'GHG (CO2 only, historical + RCP4.5, carbon cycle only) Sl SS Ds SD BC MD OC Oz AA (all fixed at 1850 values)'
     case ('esmHistorical')
        exp(i)%forcing(1:)      = 'Sl GHG Vl SS Ds SD BC MD OC Oz AA LU'
     case ('esmrcp85')
        exp(i)%forcing(1:)      = 'Sl GHG SS Ds SA BC MD OC Oz AA'
     case ('historical')
        exp(i)%forcing(1:)      = 'Sl GHG Vl SS Ds SD BC MD OC Oz AA LU'
     case ('historicalExt')
        exp(i)%forcing(1:)      = 'unknown forcings'
     case ('historicalGHG')
        exp(i)%forcing(1:)      = 'GHG (time-varying over course of simulation), Sl Vl LU SS Ds SD BC MD OC Oz AA (all fixed at or cycled over 1850 values)'
     case ('historicalMisc')
        select case (trim(adjustl(exp(i)%case)))
        case ('b40.20th.aero.1deg.006','b40.20th.aero.1deg.008','b40.20th.aero.1deg.012')
           exp(i)%forcing(1:) = 'SS Ds SD BC MD OC AA (time-varying over course of simulation), Sl GHG Vl LU Oz (all fixed at or cycled over 1850 values)'
        case ('b40.20th.anthro.1deg.006','b40.20th.anthro.1deg.008','b40.20th.anthro.1deg.009','b40.20th.anthro.1deg.012')
           exp(i)%forcing(1:) = 'GHG Oz SS Ds SD BC MD OC AA LU (time-varying over course of simulation), Sl Vl (all fixed at or cycled over 1850 values)'
        case ('b40.20th.bcarb.1deg.008')
           exp(i)%forcing(1:) = 'BC OC (time-varying over course of simulation), GHG MD AA LU Oz SS Ds SD Sl Vl (all fixed at or cycled over 1850 values)'
        case ('b40.20th_SFland.1deg.001','b40.20th.land.1deg.008','b40.20th.land.1deg.012')
           exp(i)%forcing(1:) = 'LU (time-varying over course of simulation), GHG Oz SS Ds SD BC MD OC AA LU Sl Vl (all fixed at or cycled over 1850 values)'
        case ('b40.20th.oz.1deg.006','b40.20th.oz.1deg.008','b40.20th.oz.1deg.012')
           exp(i)%forcing(1:) = 'Oz (time-varying over course of simulation), GHG LU SS Ds SD BC MD OC AA LU Sl Vl (all fixed at or cycled over 1850 values)'
        case ('b40.20th.so4.1deg.008')
           exp(i)%forcing(1:) = 'SD (time-varying over course of simulation), GHG LU SS Ds Oz BC MD OC AA LU Sl Vl (all fixed at or cycled over 1850 values)'
        case ('b40.20th.solar.1deg.006','b40.20th.solar.1deg.008','b40.20th.solar.1deg.012')
           exp(i)%forcing(1:) = 'Sl (time-varying over course of simulation), Vl GHG LU SS Ds SD BC MD OC Oz AA (all fixed at or cycled over 1850 values)'
        case ('b40.20th.volc.1deg.008','b40.20th.volc.1deg.006','b40.20th.volc.1deg.012')
           exp(i)%forcing(1:) = 'Vl (time-varying over course of simulation), Sl GHG LU SS Ds SD BC MD OC Oz AA (all fixed at or cycled over 1850 values)'
        case default
           write(*,*) 'Unknown historicalMisc case: ',trim(adjustl(exp(i)%case))
           stop
        end select
     case ('historicalNat')
        exp(i)%forcing(1:)      = 'Sl Vl (time-varying over course of simulation), GHG LU SS Ds SD BC MD OC Oz AA (all fixed at or cycled over 1850 values)'
     case ('lgm')
        exp(i)%forcing(1:)      = 'Sl GHG Vl SS Ds SD BC MD OC Oz AA LU (all fixed at or cycled over 1850 values)'
     case ('midHolocene')
        exp(i)%forcing(1:)      = 'Sl GHG Vl SS Ds SD BC MD OC Oz AA LU (all fixed at or cycled over 1850 values)'
     case ('noVolcXXXX')
        exp(i)%forcing(1:)      = 'unknown forcings'
     case ('past1000')
        exp(i)%forcing(1:)      = 'Sl GHG Vl LU (time-varying over course of simulation), SS Ds SD BC MD OC Oz AA (all fixed at or cycled over 1850 values)'
     case ('piControl')
        exp(i)%forcing(1:)      = 'Sl GHG SS Ds SD BC MD OC Oz AA (all fixed at 1850 values)'
     case ('rcp26')
        exp(i)%forcing(1:)      = 'Sl GHG SS Ds SA BC MD OC Oz AA'
     case ('rcp45')
        exp(i)%forcing(1:)      = 'Sl GHG SS Ds SA BC MD OC Oz AA'
     case ('rcp60')
        exp(i)%forcing(1:)      = 'Sl GHG SS Ds SA BC MD OC Oz AA'
     case ('rcp85')
        exp(i)%forcing(1:)      = 'Sl GHG SS Ds SA BC MD OC Oz AA'
     case ('sst2030')
        exp(i)%forcing(1:)      = 'unknown forcings'
     case ('sstClim')
        exp(i)%forcing(1:)      = 'Sl GHG Vl SS Ds SA BC MD OC Oz AA (observed sea surface temps)'
     case ('sstClim4xCO2')
        exp(i)%forcing(1:)      = 'Sl GHG (CO2 to 4X) Vl SS Ds SA BC MD OC Oz AA (observed sea surface temps)'
     case ('sstClimAerosol')
        exp(i)%forcing(1:)      = 'unknown forcings'
     case ('sstClimSulfate')
        exp(i)%forcing(1:)      = 'unknown forcings'
     case ('volcIn2010')
        exp(i)%forcing(1:)      = 'unknown forcings'
     case ('G1')
        exp(i)%forcing(1:)      = 'GHG (CO2 at 4XCO2) Sl (reduced to balance TOA) SS Ds SD BC MD OC Oz AA LU'
     case ('G2')
        exp(i)%forcing(1:)      = 'GHG (CO2 at 1%/year increase) Sl (reduced to balance TOA) SS Ds SD BC MD OC Oz AA LU'
     case ('G3S')
        exp(i)%forcing(1:)      = 'GHG (RCP4.5) Sl (reduced to balance TOA) SS Ds SD BC MD OC Oz AA LU'
     case ('tamip200810','tamip200901','tamip200904','tamip200907')
        exp(i)%forcing(1:)      = 'Sl GHG Vl SS Ds SA BC MD OC Oz AA'
     end select
  enddo
  !
  ! Add grid long name to grid info
  !
  do i = 1,num_exp
     select case (trim(adjustl(exp(i)%grid)))
     case ('f09_f09')
        exp(i)%grid(1:) = trim(adjustl(exp(i)%grid))//' (0.9x1.25_0.9x1.25)'
     case ('f09_g16')
        exp(i)%grid(1:) = trim(adjustl(exp(i)%grid))//' (0.9x1.25_gx1v6)'
     case ('f19_f19')
        exp(i)%grid(1:) = trim(adjustl(exp(i)%grid))//' (1.9x2.5_1.9x2.5)'
     case ('f19_g16')
        exp(i)%grid(1:) = trim(adjustl(exp(i)%grid))//' (1.9x2.5_gx1v6)'
     case ('T31_g37')
        exp(i)%grid(1:) = trim(adjustl(exp(i)%grid))//' (T31_gx3v7)'
     case default
        exp(i)%grid(1:) = 'unknown resolution'
     end select
  enddo
  write(*,'(''Experiments loaded: '',i5,'' entries.'')') num_exp
end subroutine load_exp
