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
  character(len=256)::instring,exp_file
  logical::does_exist
  !
  ! Get experiment information
  !
  ! Master CCSM experiment list
  !
  !  Columns   Field
  !   1 -  40  CCSM case name
  !  45 -  59  model_id (CCSM4, CESM1, CCSM4-BGC, CCSM4-FSCHEM, CCSM4-WACCM, etc)
  !  60 -  61  Location (NC = NCAR; NE = NERSC; OR = ORNL)
  !  65 -  79  Official MIP name, or very brief description (N/A if not applicable
  !  80 -  89  RIP code (cm5) or realization number (cm3) (N/A if not applicable)
  !  90 -  94  MIP (cm3 or cm5) experiment (N/A if not applicable)
  !  95 - 134  RUN_REFCASE (parent case)
  ! 135 - 149  RUN_REFDATE (branch date, yyyy-mm-dd)
  ! 150 - 159  years of experiment (YYYY if unknown)
  ! 160 - 169  GRID (resolution)
  ! 170 - 189  COMPSET (N/A if not applicable)
  ! 190 - 209  REPOTAG (N/A if not applicable)
  ! 210 - 229  Calendar dates of simulation execution (yyyy/mm-yyyy/mm)
  ! 230 - 239  MACH (hardware)
  ! 240 - end  DOUT_L_MSROOT (history file location on archive)
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
  exp(:)%forcing(1:)       = ' '
  !
  open(20,file=trim(adjustl(exp_file)),form='formatted')
  iostat = 0 ; iexp = 1
  do while (iostat == 0)
     instring(1:)  = ' '
     read(20,'(a)',iostat=iostat) instring
     if (iostat == 0) then
        if (instring(1:1) /= '#') then
           exp(iexp)%case(1:)        = adjustl(instring(  1: 40))
           exp(iexp)%model_id(1:)    = adjustl(instring( 45: 59))
           exp(iexp)%loc(1:)         = adjustl(instring( 60: 61))
           exp(iexp)%expt_id(1:)     = adjustl(instring( 65: 79))
           exp(iexp)%rip_code(1:)    = adjustl(instring( 80: 89))
           exp(iexp)%cmip(1:)        = adjustl(instring( 90: 94))
           exp(iexp)%run_refcase(1:) = adjustl(instring( 95:134))
           exp(iexp)%run_refdate(1:) = adjustl(instring(135:149))
           exp(iexp)%begin_end(1:)   = adjustl(instring(150:159))
           exp(iexp)%grid(1:)        = adjustl(instring(160:169))
           exp(iexp)%compset(1:)     = adjustl(instring(170:189))
           exp(iexp)%repotag(1:)     = adjustl(instring(190:209))
           exp(iexp)%start_fin(1:)   = adjustl(instring(210:229))
           exp(iexp)%mach(1:)        = adjustl(instring(230:239))
           exp(iexp)%dout(1:)        = adjustl(instring(240:len_trim(instring)))
           iexp = iexp + 1
        endif
     endif
  enddo
  close(20)
  num_exp = iexp - 1
  !
  ! Parse out beginning/ending years
  do iexp = 1,num_exp
     if (exp(iexp)%begin_end /= 'YYYY-YYYY') then
        read(exp(iexp)%begin_end(1:4),'(i4.4)') exp(iexp)%begyr
        read(exp(iexp)%begin_end(6:9),'(i4.4)') exp(iexp)%endyr
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
     select case (trim(adjustl(exp(i)%expt_id)))
     case ('1pctCO2')
        exp(i)%forcing(1:)      = 'Sl GHG Vl SS Ds SA BC MD OC Oz AA'
     case ('abrupt4xCO2')
        exp(i)%forcing(1:)      = 'Sl GHG Vl SS Ds SA BC MD OC Oz AA'
     case ('amip')
        exp(i)%forcing(1:)      = 'Sl GHG Vl SS Ds SA BC MD OC Oz AA'
     case ('amip4K')
        exp(i)%forcing(1:)      = 'Sl GHG Vl SS Ds SA BC MD OC Oz AA'
     case ('amip4xCO2')
        exp(i)%forcing(1:)      = 'Sl GHG Vl SS Ds SA BC MD OC Oz AA'
     case ('amipFuture')
        exp(i)%forcing(1:)      = 'Sl GHG Vl SS Ds SA BC MD OC Oz AA'
     case ('aqua4K')
        exp(i)%forcing(1:)      = 'Sl GHG Vl SS Ds SA BC MD OC Oz AA'
     case ('aqua4xCO2')
        exp(i)%forcing(1:)      = 'Sl GHG Vl SS Ds SA BC MD OC Oz AA'
     case ('aquaControl')
        exp(i)%forcing(1:)      = 'Sl GHG Vl SS Ds SA BC MD OC Oz AA'
     case ('decadalXXXX')
        exp(i)%forcing(1:)      = 'Sl GHG Vl SS Ds SA BC MD OC Oz AA'
     case ('esmControl')
        exp(i)%forcing(1:)      = 'Sl GHG Vl SS Ds SA BC MD OC Oz AA'
     case ('esmFdbk1')
        exp(i)%forcing(1:)      = 'Sl GHG Vl SS Ds SA BC MD OC Oz AA'
     case ('esmFdbk2')
        exp(i)%forcing(1:)      = 'Sl GHG Vl SS Ds SA BC MD OC Oz AA'
     case ('esmFixClim1')
        exp(i)%forcing(1:)      = 'Sl GHG Vl SS Ds SA BC MD OC Oz AA'
     case ('esmFixClim2')
        exp(i)%forcing(1:)      = 'Sl GHG Vl SS Ds SA BC MD OC Oz AA'
     case ('esmHistorical')
        exp(i)%forcing(1:)      = 'Sl GHG Vl SS Ds SA BC MD OC Oz AA'
     case ('esmrcp85')
        exp(i)%forcing(1:)      = 'Sl GHG Vl SS Ds SA BC MD OC Oz AA'
     case ('historical')
        exp(i)%forcing(1:)      = 'Sl GHG Vl SS Ds SD BC MD OC Oz AA LU'
!        select case (trim(adjustl(exp(i)%case)))
!        case ('b40.20th.track1.deg.005','b40.20th.track1.deg.006','b40.20th.track1.deg.007','b40.20th.track1.deg.008','b40.20th.track1.deg.009','b40.20th.track1.deg.012')
!        case ('b40.20th.track1.1deg.008')
!!$           exp(i)%forcing_note(1:) = '\n'//&
!!$                                     'Sl  : SOLAR_TSI_Lean_1610-2007_annual_c090324.nc \n'//&
!!$                                     'GHG : ghg_hist_1765-2005_c091218.nc \n'//&
!!$                                     'Vl  : CCSM4_volcanic_1850-2008_prototype1.nc \n'//&
!!$                                     'SS  : ssam_camrt_c080918.nc, sscm_camrt_c080918.nc \n'//&
!!$                                     'Ds  :  \n'//&
!!$                                     'MD  :  \n'//&
!!$                                     'SD  : aero_1.9x2.5_L26_1850-2005_c091112.nc, sulfate_camrt_c080918.nc \n'//&
!!$                                     'BC  : bcpho_camrt_c080918.nc, bcphi_camrt_c080918.nc \n'//&
!!$                                     'OC  : ocpho_camrt_c080918.nc, ocphi_camrt_c080918.nc \n'//&
!!$                                     'Oz  : ozone_1.9x2.5_L26_1850-2005_c091112.nc \n'//&
!!$                                     'AA  :  \n'//&
!!$                                     'LU  :  \n'
!        case default
!           write(*,*) 'Unknown historical case: ',trim(adjustl(exp(i)%case))
!        end select
     case ('historicalExt')
        exp(i)%forcing(1:)      = 'Sl GHG Vl SS Ds SA BC MD OC Oz AA'
     case ('historicalGHG')
        exp(i)%forcing(1:)      = 'Sl GHG Vl SS Ds SA BC MD OC Oz AA'
     case ('historicalMisc')
        select case (trim(adjustl(exp(i)%case)))
        case ('b40.20th.aero.1deg.006')
           exp(i)%forcing(1:)      = 'Sl GHG Vl SS Ds SA BC MD OC Oz AA'
        case ('b40.20th.aero.1deg.008')
           exp(i)%forcing(1:)      = 'Sl GHG Vl SS Ds SA BC MD OC Oz AA'
        case ('b40.20th.aero.1deg.012')
           exp(i)%forcing(1:)      = 'Sl GHG Vl SS Ds SA BC MD OC Oz AA'
        case ('b40.20th.anthro.1deg.006')
           exp(i)%forcing(1:)      = 'Sl GHG Vl SS Ds SA BC MD OC Oz AA'
        case ('b40.20th.anthro.1deg.008')
           exp(i)%forcing(1:)      = 'Sl GHG Vl SS Ds SA BC MD OC Oz AA'
        case ('b40.20th.anthro.1deg.009')
           exp(i)%forcing(1:)      = 'Sl GHG Vl SS Ds SA BC MD OC Oz AA'
        case ('b40.20th.anthro.1deg.012')
           exp(i)%forcing(1:)      = 'Sl GHG Vl SS Ds SA BC MD OC Oz AA'
        case ('b40.20th.bcarb.1deg.008')
           exp(i)%forcing(1:)      = 'Sl GHG Vl SS Ds SA BC MD OC Oz AA'
        case ('b40.20th_SFland.1deg.001')
           exp(i)%forcing(1:)      = 'Sl GHG Vl SS Ds SA BC MD OC Oz AA'
        case ('b40.20th.land.1deg.008')
           exp(i)%forcing(1:)      = 'Sl GHG Vl SS Ds SA BC MD OC Oz AA'
        case ('b40.20th.land.1deg.012')
           exp(i)%forcing(1:)      = 'Sl GHG Vl SS Ds SA BC MD OC Oz AA'
        case ('b40.20th.oz.1deg.006')
           exp(i)%forcing(1:)      = 'Sl GHG Vl SS Ds SA BC MD OC Oz AA'
        case ('b40.20th.oz.1deg.008')
           exp(i)%forcing(1:)      = 'Sl GHG Vl SS Ds SA BC MD OC Oz AA'
        case ('b40.20th.oz.1deg.012')
           exp(i)%forcing(1:)      = 'Sl GHG Vl SS Ds SA BC MD OC Oz AA'
        case ('b40.20th.so4.1deg.008')
           exp(i)%forcing(1:)      = 'Sl GHG Vl SS Ds SA BC MD OC Oz AA'
        case ('b40.20th.solar.1deg.006')
           exp(i)%forcing(1:)      = 'Sl GHG Vl SS Ds SA BC MD OC Oz AA'
        case ('b40.20th.solar.1deg.008')
           exp(i)%forcing(1:)      = 'Sl GHG Vl SS Ds SA BC MD OC Oz AA'
        case ('b40.20th.solar.1deg.012')
           exp(i)%forcing(1:)      = 'Sl GHG Vl SS Ds SA BC MD OC Oz AA'
        case ('b40.20th.volc.1deg.008')
           exp(i)%forcing(1:)      = 'Sl GHG Vl SS Ds SA BC MD OC Oz AA'
        case ('b40.20th.volc.1deg.006')
           exp(i)%forcing(1:)      = 'Sl GHG Vl SS Ds SA BC MD OC Oz AA'
        case ('b40.20th.volc.1deg.012')
           exp(i)%forcing(1:)      = 'Sl GHG Vl SS Ds SA BC MD OC Oz AA'
        case default
           write(*,*) 'Unknown historicalMisc case: ',trim(adjustl(exp(i)%case))
           stop
        end select
     case ('historicalNat')
        exp(i)%forcing(1:)      = 'Sl GHG Vl SS Ds SA BC MD OC Oz AA'
     case ('lgm')
        exp(i)%forcing(1:)      = 'Sl GHG Vl SS Ds SA BC MD OC Oz AA'
     case ('midHolocene')
        exp(i)%forcing(1:)      = 'Sl GHG Vl SS Ds SA BC MD OC Oz AA'
     case ('noVolcXXXX')
        exp(i)%forcing(1:)      = 'Sl GHG Vl SS Ds SA BC MD OC Oz AA'
     case ('past1000')
        exp(i)%forcing(1:)      = 'Sl GHG Vl SS Ds SA BC MD OC Oz AA'
     case ('piControl')
        exp(i)%forcing(1:)      = 'Sl GHG Vl SS Ds SA BC MD OC Oz AA'
     case ('rcp26')
        exp(i)%forcing(1:)      = 'Sl GHG Vl SS Ds SA BC MD OC Oz AA'
     case ('rcp45')
        exp(i)%forcing(1:)      = 'Sl GHG Vl SS Ds SA BC MD OC Oz AA'
     case ('rcp60')
        exp(i)%forcing(1:)      = 'Sl GHG Vl SS Ds SA BC MD OC Oz AA'
     case ('rcp85')
        exp(i)%forcing(1:)      = 'Sl GHG Vl SS Ds SA BC MD OC Oz AA'
     case ('sst2030')
        exp(i)%forcing(1:)      = 'Sl GHG Vl SS Ds SA BC MD OC Oz AA'
     case ('sstClim')
        exp(i)%forcing(1:)      = 'Sl GHG Vl SS Ds SA BC MD OC Oz AA'
     case ('sstClim4xCO2')
        exp(i)%forcing(1:)      = 'Sl GHG Vl SS Ds SA BC MD OC Oz AA'
     case ('sstClimAerosol')
        exp(i)%forcing(1:)      = 'Sl GHG Vl SS Ds SA BC MD OC Oz AA'
     case ('sstClimSulfate')
        exp(i)%forcing(1:)      = 'Sl GHG Vl SS Ds SA BC MD OC Oz AA'
     case ('volcIn2010')
        exp(i)%forcing(1:)      = 'Sl GHG Vl SS Ds SA BC MD OC Oz AA'
     case default
        exp(i)%forcing(1:)      = 'unknown forcings'
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
