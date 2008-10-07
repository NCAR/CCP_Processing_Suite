PROGRAM main

  !-----------------------------------------------------------------------------
  !   This program regrids POP dipole grid data to a regular lat-lon grid
  !   using SCRIP-computed remapping weights stored in 'map_file', 
  !   and optionally, to a standard vertical grid.  Seperate SCRIP remapping 
  !   files should exist for each k-level of the native POP grid.
  ! 
  !   Usage:   
  !		src_file       == POP output on original grid
  !	        out_form       == Regridded output netcdf template
  !		dst_lclobber   == TRUE allows clobbering dst_file
  !		dst_z_standard == TRUE regrids to a standard z grid
  !				  FALSE preserves native z grid
  !		z_file         == data file which defines the standard z grid
  !		nk             == # of native z levels to process
  !		src_var(s)     == variable(s) to be processed
  !		map_file       == SCRIP horizontal remapping file template 
  !
  !   If dst_z_standard = .true. remaps vertically to a standard
  !   grid given by z_std, otherwise just performs horizontal
  !   regridding down to nk levels (must have nk corresponding
  !   mapping files).
  !
  !   (Based on K. Lindsay's /fs/cgd/data0/klindsay/POP_tools/regrid )
  !
  !   CVS:$Id: main_ccsm_ocn.F90,v 1.1 2005/03/11 21:16:42 strandwg Exp $
  !-----------------------------------------------------------------------------

  USE kinds_mod
  USE vars_mod
  USE popregrid

  !*****************************************************************************

  IMPLICIT NONE

  !*****************************************************************************

  NAMELIST /popregrid_nml/ &
       src_file, src_var, &
       nk, map_file, z_file, &
       out_form, dst_lclobber, dst_z_standard, dst_var

  INTEGER(KIND=int_kind), PARAMETER :: stdin = 5

  INTEGER(KIND=int_kind) :: ios,ix,jx,ivar

  !*****************************************************************************

  src_file       = 'unknown'
  src_var        = 'unknown'
  out_form       = 'unknown'
  dst_file       = 'unknown'
  dst_lclobber   = .TRUE.
  dst_z_standard = .TRUE.
  map_file       = 'unknown'
  z_file         = 'unknown'
  nk             = 1

  DO
     dst_var = 'unknown'
     READ(UNIT=stdin, NML=popregrid_nml, IOSTAT=ios, ERR=10, END=20)
     dst_var = src_var
     nvars = 1
     DO ivar = 1,max_vars
        IF ((INDEX(src_var(ivar), 'unknown')) /= 1) nvars = nvars + 1
     END DO
     nvars = nvars - 1
     WRITE(*,*) 'src_file        = ',TRIM(src_file)
     WRITE(*,*) 'nvars           = ',nvars
     ix = INDEX(out_form,'src_var')-1
     jx = ix + len_trim('src_var')+1
     DO ivar = 1,nvars
        dst_file(ivar) = out_form(1:ix)//trim(src_var(ivar))//out_form(jx:)
        WRITE(*,'('' src/dst_var      = '',a8,'' dst_file = '',a)') &
             src_var(ivar),TRIM(dst_file(ivar))
     END DO
     WRITE(*,*) 'map_file        = ', TRIM(map_file)
     WRITE(*,*) 'dst_z_standard  = ', dst_z_standard
     WRITE(*,*) 'z_file          = ', TRIM(z_file)
     CALL regrid_driver
     WRITE(*,*) ' '
     CYCLE
10   WRITE(*,*) 'error reading namelist, iostat = ', ios
20   EXIT
  END DO

END PROGRAM main
