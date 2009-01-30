MODULE ccsm_ice_rgd

  !-----------------------------------------------------------------------------
  !   CVS:$Id: ccsm_ice_rgd.F90,v 1.4 2005/04/26 17:52:16 strandwg Exp $
  !-----------------------------------------------------------------------------

  USE kinds_mod
  USE vars_mod
  USE msg_mod
  USE char_case
  USE const_mod
  USE nf_wrap
  USE nf_tools

  IMPLICIT NONE
  SAVE

  !*****************************************************************************

  INTEGER(KIND=int_kind),DIMENSION(max_vars)::src_varid, src_varid2
  INTEGER(KIND=int_kind)::map_ncid,src_ncid,src_unlimdimid
  INTEGER(KIND=int_kind)::grid_ncid, grid_varid
  INTEGER(KIND=int_kind)::src_imt, src_jmt, timelen
  INTEGER(KIND=int_kind)::src_grid_size,dst_grid_size,num_links

  CHARACTER(LEN=char_len):: src_var_units,src_var_lname, &
       src_depth_name, src_depth_edges_name, src_time_name, &
       src_var_coords,src_lon_name,src_lat_name,zunits

  LOGICAL(KIND=log_kind),DIMENSION(max_vars) :: src_var_tgrid,src_var_rotate

  INTEGER(KIND=int_kind), DIMENSION(2) ::  src_grid_dims,dst_grid_dims

  INTEGER(KIND=int_kind), DIMENSION(:), ALLOCATABLE :: &
       src_grid_imask,dst_grid_imask,src_address,dst_address

  REAL(KIND=dbl_kind), DIMENSION(:,:), ALLOCATABLE :: &
       src_grid_corner_lat,src_grid_corner_lon,dst_grid_corner_lat, &
       dst_grid_corner_lon

  REAL(KIND=dbl_kind), DIMENSION(:), ALLOCATABLE :: &
       src_grid_center_lat,dst_grid_center_lat,src_grid_center_lon, &
       dst_grid_center_lon,src_grid_area,dst_grid_area,src_grid_frac, &
       dst_grid_frac,remap_matrix,dtime

  REAL(KIND=real_kind), DIMENSION(:), ALLOCATABLE :: &
       src_lon,src_lat,src_depth, src_depth_edges,time

  REAL(KIND=real_kind), DIMENSION(:), ALLOCATABLE :: dst_depth

  INTEGER(KIND=int_kind)::dst_kmin, dst_kmax

  INTEGER(KIND=int_kind), DIMENSION(:,:), ALLOCATABLE :: IMASK

  REAL(KIND=real_kind), DIMENSION(:,:), ALLOCATABLE :: &
       TAREA, AREA, ANGLE, raw_mask, area_mask
  REAL(KIND=real_kind), DIMENSION(:,:), ALLOCATABLE :: &
       LON2D, LAT2D
  REAL(KIND=real_kind), DIMENSION(:), ALLOCATABLE :: &
       LON, LAT

  INTEGER(KIND=int_kind),DIMENSION(max_vars)::dst_ncid, dst_varid

  LOGICAL(KIND=log_kind) :: dst_needs_time, &
       dst_needs_depth, dst_needs_LAT, dst_needs_LON, &
       dst_needs_AREA, dst_needs_IMASK

  INTEGER(KIND=int_kind) :: dst_time_dimid, &
       dst_X_dimid, dst_Y_dimid, &
       dst_depth_dimid, dst_LAT_varid, dst_LON_varid, &
       dst_AREA_varid, dst_IMASK_varid

  CHARACTER(LEN=char_len) :: dst_time_name, &
       dst_X_name, dst_Y_name, &
       dst_depth_name, dst_LAT_name, dst_LON_name, &
       dst_AREA_name, dst_IMASK_name

  REAL(KIND=real_kind), PARAMETER :: default_msv = 1.0E+30

  REAL(KIND=dbl_kind), PARAMETER :: radius = 6371.22e5    ! radius of Earth (cm)

  !*****************************************************************************

CONTAINS

  !*****************************************************************************

  SUBROUTINE regrid_driver

    !---------------------------------------------------------------------------
    !   local variables
    !---------------------------------------------------------------------------

    CHARACTER(LEN=*), PARAMETER :: sub_name = 'regrid_driver'
    LOGICAL(KIND=log_kind) :: lstat_ok
    INTEGER(KIND=int_kind) :: k,l,ivar

!    WRITE(*,*) 'Entering ',trim(sub_name)
    CALL read_map_file(1)
    CALL read_src_file_info
    DO ivar = 1,nvars
       CALL dst_depth_init(ivar)
       CALL create_dst_file(ivar)
       CALL def_dst_file(ivar)
       CALL nf_enddef_wrap(dst_ncid(ivar))
       CALL put_aux_vars(ivar)
       DO k=1,max(src_km(ivar),1)
          CALL read_map_file(k)
          DO l=1,MAX(1,timelen)
             CALL regrid_wrap(k,l,ivar)
          END DO
       END DO
       CALL nf_close_wrap(dst_ncid(ivar))
    END DO
    CALL nf_close_wrap(src_ncid)
    CALL free_vars
!    WRITE(*,*) 'Exiting  ',trim(sub_name)

  END SUBROUTINE regrid_driver

  !*****************************************************************************

  SUBROUTINE read_src_file_info

    !---------------------------------------------------------------------------
    !   Assumed to be a POP history file
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    !   local variables
    !---------------------------------------------------------------------------

    CHARACTER(LEN=*), PARAMETER :: sub_name = 'read_src_file_info'
    INTEGER(KIND=int_kind)::grid_ncid, grid_varid
    INTEGER(KIND=int_kind) :: ndims,gndims,gvarid,avarid,varid,ivar,i
    INTEGER(KIND=int_kind), DIMENSION(4) :: dimids,gdimids

    !---------------------------------------------------------------------------
!    WRITE(*,*) 'Entering ',trim(sub_name)
    src_var_tgrid  = .false.
    src_var_rotate = .false.
    DO ivar = 1,nvars
       IF (src_var(ivar)=='u'.or.src_var(ivar)=='v') THEN
          src_var_rotate(ivar) = .true.
!          CALL msg_write(sub_name, 'vector rotating ',TRIM(src_var(ivar)))
          SELECT CASE (src_var(ivar))
          CASE ('u')
             src_var2(ivar) = 'v'
          CASE ('v')
             src_var2(ivar) = 'u'
          END SELECT
       END IF
    END DO

    CALL nf_open_wrap(src_file, 0, src_ncid)
    DO ivar = 1,nvars
       CALL nf_inq_varid_wrap(src_ncid, src_var(ivar), src_varid(ivar))
       CALL nf_inq_varndims_wrap(src_ncid, src_varid(ivar), ndims)
       
       IF (src_var_rotate(ivar)) &
            CALL nf_inq_varid_wrap(src_ncid,src_var2(ivar),src_varid2(ivar))
       
       IF (ndims > 4 .OR. ndims < 2) THEN
          CALL msg_write(sub_name, 'illegal ndims ', ndims)
          STOP
       END IF

       CALL nf_inq_vardimid_wrap(src_ncid, src_varid(ivar), dimids)
    
       CALL nf_inq_dim_wrap(src_ncid, dimids(1), src_lon_name, src_imt)
       IF ((INDEX(upper(src_lon_name), 'LON') == 0) .AND. &
            (INDEX(upper(src_lon_name), 'X') == 0)) THEN
          CALL msg_write(sub_name, 'unexpected src_lon_name ', &
               TRIM(src_lon_name))
          STOP
       END IF
       
       CALL nf_inq_dim_wrap(src_ncid, dimids(2), src_lat_name, src_jmt)
       IF ((INDEX(upper(src_lat_name), 'LAT') == 0) .AND. &
           (INDEX(upper(src_lat_name), 'Y') == 0)) THEN
          CALL msg_write(sub_name, 'unexpected src_lat_name ', &
               TRIM(src_lat_name))
          STOP
       END IF

       src_km(ivar) = 0
       timelen = 0

       IF (ndims > 2) THEN
          CALL nf_inq_dim_wrap(src_ncid, dimids(3), src_time_name, timelen)
          IF ((INDEX(upper(src_time_name), 'TIME') /= 1) .AND. &
              (INDEX(upper(src_time_name), 'DATE') /= 1) .AND. &
              (INDEX(upper(src_time_name), 'T') /= 1) .AND. &
              (INDEX(upper(src_time_name), 'MONTH') /= 1)) THEN
             CALL msg_write(sub_name, 'unexpected src_time_name ', &
                  TRIM(src_time_name))
             STOP
          END IF
       END IF
    END DO
    !
    !  get ANGLE
    IF (.NOT.(ALLOCATED(ANGLE))) ALLOCATE(ANGLE(src_imt,src_jmt))
    CALL nf_inq_varid_wrap(src_ncid, 'ANGLE', gvarid)
    CALL nf_get_var_wrap(src_ncid, gvarid, ANGLE)
    !
    ! Get and create mask using 'aice'
    IF (.NOT.(ALLOCATED(raw_mask)))   ALLOCATE(raw_mask(src_imt,src_jmt))
    IF (.NOT.(ALLOCATED(area_mask))) ALLOCATE(area_mask(dst_imt,dst_jmt))
    CALL nf_inq_varid_wrap(src_ncid, 'aice', avarid)
    CALL nf_get_var_wrap(src_ncid, avarid, raw_mask)
    CALL regrid_level(default_msv, raw_mask, area_mask)
    WHERE (area_mask <= 1.)
       area_mask = 0.
    ELSEWHERE
       area_mask = 1.
    ENDWHERE
    do ivar = 1,nvars
       IF (src_km(ivar) < nk) nk = MAX(1,src_km(ivar))
       IF (timelen > 0) THEN
          IF (.NOT.ALLOCATED(time)) ALLOCATE(time(timelen))
          IF (src_km(ivar) > 0) THEN
             CALL nf_read_dim(src_ncid, dimids(3), time)
          ELSE
             CALL nf_read_dim(src_ncid, dimids(3), time)
          END IF
       END IF
    END DO
!    WRITE(*,*) 'Exiting  ',trim(sub_name)

  END SUBROUTINE read_src_file_info

  !*****************************************************************************

  SUBROUTINE read_map_file(k)

    !---------------------------------------------------------------------------
    !   Assumes 'scrip' convention SCRIP remapping file
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    !   arguments
    !---------------------------------------------------------------------------

    INTEGER(KIND=int_kind), INTENT(IN) :: k

    !---------------------------------------------------------------------------
    !   local variables
    !---------------------------------------------------------------------------

    CHARACTER(LEN=*), PARAMETER :: sub_name = 'read_map_file'
    CHARACTER(LEN=char_len) :: tmpstr
    CHARACTER(LEN=2) :: kstr

    INTEGER(KIND=int_kind) :: ndims,ki,k0,varid,nx,ny,in
    INTEGER(KIND=int_kind), DIMENSION(4) :: dimids

    !---------------------------------------------------------------------------

    ! Overwrites map file with new k value
    ki = INDEX(map_file, 'k')
    k0 = IACHAR('0')
    kstr = CHAR(k0+k/10)//CHAR(k0+MOD(k,10))
    map_file(ki+1:ki+2) = kstr

!    WRITE(*,*) 'Entering ',trim(sub_name)
!    WRITE(*,*) 'Opening  ',trim(map_file)
    CALL nf_open_wrap(map_file, 0, map_ncid)
    CALL nf_inq_varid_wrap(map_ncid, 'src_grid_dims', varid)
    CALL nf_get_vara_wrap(map_ncid, varid, (/ 1 /), (/ 2 /), src_grid_dims)
    CALL nf_inq_varid_wrap(map_ncid, 'dst_grid_dims', varid)
    CALL nf_get_vara_wrap(map_ncid, varid, (/ 1 /), (/ 2 /), dst_grid_dims)
    CALL nf_inq_varid_wrap(map_ncid, 'src_grid_center_lat', varid)
    CALL nf_inq_vardimid_wrap(map_ncid, varid, dimids)
    CALL nf_inq_dim_wrap(map_ncid, dimids(1), tmpstr, src_grid_size)
    IF (ALLOCATED(src_grid_center_lat)) THEN
        DEALLOCATE(src_grid_center_lat)
        ALLOCATE(src_grid_center_lat(src_grid_size))
    ELSE
        ALLOCATE(src_grid_center_lat(src_grid_size))
    END IF
    IF (ALLOCATED(src_grid_center_lon)) THEN
        DEALLOCATE(src_grid_center_lon)
        ALLOCATE(src_grid_center_lon(src_grid_size))
    ELSE
        ALLOCATE(src_grid_center_lon(src_grid_size))
    END IF
    IF (ALLOCATED(src_grid_area)) THEN
        DEALLOCATE(src_grid_area)
        ALLOCATE(src_grid_area(src_grid_size))
    ELSE
        ALLOCATE(src_grid_area(src_grid_size))
    END IF
    IF (ALLOCATED(src_grid_imask)) THEN
        DEALLOCATE(src_grid_imask)
        ALLOCATE(src_grid_imask(src_grid_size))
    ELSE
        ALLOCATE(src_grid_imask(src_grid_size))
    END IF
    IF (ALLOCATED(src_grid_frac)) THEN
        DEALLOCATE(src_grid_frac)
        ALLOCATE(src_grid_frac(src_grid_size))
    ELSE
        ALLOCATE(src_grid_frac(src_grid_size))
    END IF
    IF (ALLOCATED(src_grid_corner_lat)) THEN
        DEALLOCATE(src_grid_corner_lat)
        ALLOCATE(src_grid_corner_lat(src_grid_size,4))
    ELSE
        ALLOCATE(src_grid_corner_lat(src_grid_size,4))
    END IF
    IF (ALLOCATED(src_grid_corner_lon)) THEN
        DEALLOCATE(src_grid_corner_lon)
        ALLOCATE(src_grid_corner_lon(src_grid_size,4))
    ELSE
        ALLOCATE(src_grid_corner_lon(src_grid_size,4))
    END IF
    CALL nf_get_vara_wrap(map_ncid, varid, (/ 1 /), (/ src_grid_size /), &
         src_grid_center_lat)
    CALL nf_inq_varid_wrap(map_ncid, 'src_grid_center_lon', varid)
    CALL nf_get_vara_wrap(map_ncid, varid, (/ 1 /), (/ src_grid_size /), &
         src_grid_center_lon)
    CALL nf_inq_varid_wrap(map_ncid, 'src_grid_area', varid)
    CALL nf_get_vara_wrap(map_ncid, varid, (/ 1 /), (/ src_grid_size /), &
         src_grid_area)
    CALL nf_inq_varid_wrap(map_ncid, 'src_grid_imask', varid)
    CALL nf_get_vara_wrap(map_ncid, varid, (/ 1 /), (/ src_grid_size /), &
         src_grid_imask)
    CALL nf_inq_varid_wrap(map_ncid, 'src_grid_frac', varid)
    CALL nf_get_vara_wrap(map_ncid, varid, (/ 1 /), (/ src_grid_size /), &
         src_grid_frac)
    CALL nf_inq_varid_wrap(map_ncid, 'src_grid_corner_lat', varid)
    CALL nf_get_vara_wrap(map_ncid, varid, (/ 1,1 /), (/ 4,src_grid_size /), &
         src_grid_corner_lat)
    CALL nf_inq_varid_wrap(map_ncid, 'src_grid_corner_lon', varid)
    CALL nf_get_vara_wrap(map_ncid, varid, (/ 1,1 /), (/ 4,src_grid_size /), &
         src_grid_corner_lon)
    CALL nf_inq_varid_wrap(map_ncid, 'dst_grid_center_lat', varid)
    CALL nf_inq_vardimid_wrap(map_ncid, varid, dimids)
    CALL nf_inq_dim_wrap(map_ncid, dimids(1), tmpstr, dst_grid_size)
    IF (ALLOCATED(dst_grid_center_lat)) THEN
       DEALLOCATE(dst_grid_center_lat)
       ALLOCATE(dst_grid_center_lat(dst_grid_size))
    ELSE
       ALLOCATE(dst_grid_center_lat(dst_grid_size))
    END IF
    IF (ALLOCATED(dst_grid_center_lon)) THEN
       DEALLOCATE(dst_grid_center_lon)
       ALLOCATE(dst_grid_center_lon(dst_grid_size))
    ELSE
       ALLOCATE(dst_grid_center_lon(dst_grid_size))
    END IF
    IF (ALLOCATED(dst_grid_area)) THEN
       DEALLOCATE(dst_grid_area)
       ALLOCATE(dst_grid_area(dst_grid_size))
    ELSE
       ALLOCATE(dst_grid_area(dst_grid_size))
    END IF
    IF (ALLOCATED(dst_grid_imask)) THEN
       DEALLOCATE(dst_grid_imask)
       ALLOCATE(dst_grid_imask(dst_grid_size))
    ELSE
       ALLOCATE(dst_grid_imask(dst_grid_size))
    END IF
    IF (ALLOCATED(dst_grid_frac)) THEN
       DEALLOCATE(dst_grid_frac)
       ALLOCATE(dst_grid_frac(dst_grid_size))
    ELSE
       ALLOCATE(dst_grid_frac(dst_grid_size))
    END IF
    IF (ALLOCATED(dst_grid_corner_lat)) THEN
       DEALLOCATE(dst_grid_corner_lat)
       ALLOCATE(dst_grid_corner_lat(dst_grid_size,4))
    ELSE
       ALLOCATE(dst_grid_corner_lat(dst_grid_size,4))
    END IF
    IF (ALLOCATED(dst_grid_corner_lon)) THEN
       DEALLOCATE(dst_grid_corner_lon)
       ALLOCATE(dst_grid_corner_lon(dst_grid_size,4))
    ELSE
       ALLOCATE(dst_grid_corner_lon(dst_grid_size,4))
    END IF
    CALL nf_get_vara_wrap(map_ncid, varid, (/ 1 /), (/ dst_grid_size /), &
         dst_grid_center_lat)
    CALL nf_inq_varid_wrap(map_ncid, 'dst_grid_center_lon', varid)
    CALL nf_get_vara_wrap(map_ncid, varid, (/ 1 /), (/ dst_grid_size /), &
         dst_grid_center_lon)
    CALL nf_inq_varid_wrap(map_ncid, 'dst_grid_area', varid)
    CALL nf_get_vara_wrap(map_ncid, varid, (/ 1 /), (/ dst_grid_size /), &
         dst_grid_area)
    CALL nf_inq_varid_wrap(map_ncid, 'dst_grid_imask', varid)
    CALL nf_get_vara_wrap(map_ncid, varid, (/ 1 /), (/ dst_grid_size /), &
         dst_grid_imask)
    CALL nf_inq_varid_wrap(map_ncid, 'dst_grid_frac', varid)
    CALL nf_get_vara_wrap(map_ncid, varid, (/ 1 /), (/ dst_grid_size /), &
         dst_grid_frac)
    CALL nf_inq_varid_wrap(map_ncid, 'dst_grid_corner_lat', varid)
    CALL nf_get_vara_wrap(map_ncid, varid, (/ 1,1 /), (/ 4,dst_grid_size /), &
         dst_grid_corner_lat)
    CALL nf_inq_varid_wrap(map_ncid, 'dst_grid_corner_lon', varid)
    CALL nf_get_vara_wrap(map_ncid, varid, (/ 1,1 /), (/ 4,dst_grid_size /), &
         dst_grid_corner_lon)
    dst_imt = dst_grid_dims(1)
    dst_jmt = dst_grid_dims(2)
    CALL nf_inq_varid_wrap(map_ncid, 'src_address', varid)
    CALL nf_inq_vardimid_wrap(map_ncid, varid, dimids)
    CALL nf_inq_dim_wrap(map_ncid, dimids(1), tmpstr, num_links)
    IF (ALLOCATED(src_address)) THEN
       DEALLOCATE(src_address)
       ALLOCATE(src_address(num_links))
    ELSE
       ALLOCATE(src_address(num_links))
    END IF
    IF (ALLOCATED(dst_address)) THEN
       DEALLOCATE(dst_address)
       ALLOCATE(dst_address(num_links))
    ELSE
       ALLOCATE(dst_address(num_links))
    END IF
    IF (ALLOCATED(remap_matrix)) THEN
       DEALLOCATE(remap_matrix)
       ALLOCATE(remap_matrix(num_links))
    ELSE
       ALLOCATE(remap_matrix(num_links))
    END IF
    CALL nf_get_vara_wrap(map_ncid, varid, (/ 1 /), (/ num_links /), &
         src_address)
    CALL nf_inq_varid_wrap(map_ncid, 'dst_address', varid)
    CALL nf_get_vara_wrap(map_ncid, varid, (/ 1 /), (/ num_links /), &
         dst_address)
    CALL nf_inq_varid_wrap(map_ncid, 'remap_matrix', varid)
    CALL nf_get_vara_wrap(map_ncid, varid, (/ 1,1 /), (/ 1,num_links /), &
         remap_matrix)
    CALL nf_close_wrap(map_ncid)
!    WRITE(*,*) 'Exiting  ',trim(sub_name)
  END SUBROUTINE read_map_file

  !*****************************************************************************

  SUBROUTINE nf_read_dim(ncid, dimid, dim_vals)

    !---------------------------------------------------------------------------
    !   arguments
    !---------------------------------------------------------------------------

    INTEGER(KIND=int_kind), INTENT(IN) :: ncid, dimid
    REAL(KIND=real_kind), DIMENSION(:), INTENT(OUT) :: dim_vals

    !---------------------------------------------------------------------------
    !   local variables
    !---------------------------------------------------------------------------

    CHARACTER(LEN=*), PARAMETER :: sub_name = 'nf_read_dim'

    INTEGER(KIND=int_kind) :: varid, ndims
    INTEGER(KIND=int_kind), DIMENSION(1) :: vardims
    CHARACTER(LEN=char_len) :: var_name

    !---------------------------------------------------------------------------

    CALL nf_inq_dimname_wrap(ncid, dimid, var_name)
    CALL nf_inq_varid_wrap(ncid, var_name, varid)

    CALL nf_inq_varndims_wrap(ncid, varid, ndims)
    IF (ndims /= 1) THEN
       CALL msg_write(sub_name, 'illegal ndims for ', var_name, ndims)
       STOP
    END IF

    CALL nf_inq_vardimid_wrap(ncid, varid, vardims)
    IF (vardims(1) /= dimid) THEN
       CALL msg_write(sub_name, 'illegal dimid for ', var_name, vardims(1))
       STOP
    END IF

    CALL nf_get_var_wrap(ncid, varid, dim_vals)

  END SUBROUTINE nf_read_dim

  !*****************************************************************************

  SUBROUTINE nf_read_dim_edges(ncid, dimid, dim_edge_name, dim_edge_vals)

    include 'netcdf.inc'

    !---------------------------------------------------------------------------
    !   arguments
    !---------------------------------------------------------------------------

    INTEGER(KIND=int_kind), INTENT(IN) :: ncid, dimid
    REAL(KIND=real_kind), DIMENSION(:), INTENT(OUT) :: dim_edge_vals
    CHARACTER(len=*), INTENT(OUT) :: dim_edge_name

    !---------------------------------------------------------------------------
    !   local variables
    !---------------------------------------------------------------------------

    CHARACTER(LEN=*), PARAMETER :: sub_name = 'nf_read_dim_edges'

    INTEGER(KIND=int_kind) :: stat, varid
    INTEGER(KIND=int_kind) :: ndims
    INTEGER(KIND=int_kind), DIMENSION(1) :: vardims
    INTEGER(KIND=int_kind) :: dimlen, dimlenp1, i
    REAL(KIND=real_kind), DIMENSION(:), ALLOCATABLE :: dim_vals
    REAL(KIND=real_kind) :: delta

    CHARACTER(LEN=char_len) :: var_name

    !---------------------------------------------------------------------------

    CALL nf_inq_dim_wrap(ncid, dimid, var_name, dimlen)
    CALL nf_inq_varid_wrap(ncid, var_name, varid)

    CALL nf_inq_varndims_wrap(ncid, varid, ndims)
    IF (ndims /= 1) THEN
       CALL msg_write(sub_name, 'illegal ndims for ', var_name, ndims)
       STOP
    END IF

    CALL nf_inq_vardimid_wrap(ncid, varid, vardims)
    IF (vardims(1) /= dimid) THEN
       CALL msg_write(sub_name, 'illegal dimid for ', var_name, vardims(1))
       STOP
    END IF

    dim_edge_name = ''
    CALL nf_get_att_wrap(ncid, varid, 'edges', dim_edge_name, &
         ALLOW=NF_ENOTATT, stat_out=stat)

    !---------------------------------------------------------------------------
    !   If no edge attribute exists for the variable corresponding to
    !   the dimension, generate edges.
    !---------------------------------------------------------------------------

    IF (stat == NF_ENOTATT) THEN
       dim_edge_name = TRIM(var_name) // '_edges'

       ALLOCATE(dim_vals(1:dimlen))
       CALL nf_get_var_wrap(ncid, varid, dim_vals)

       IF (dimlen == 1) THEN
          dim_edge_vals(1) = dim_vals(1)
          dim_edge_vals(2) = dim_vals(1)
       ELSE
          delta = dim_vals(2) - dim_vals(1)
          dim_edge_vals(1) = dim_vals(1) - p5 * delta
          DO i = 1,dimlen
             delta = dim_vals(i) - dim_edge_vals(i)
             dim_edge_vals(i+1) = dim_edge_vals(i) + c2 * delta
          END DO
       END IF
       DEALLOCATE(dim_vals)
    ELSE
       CALL nf_inq_varid_wrap(ncid, dim_edge_name, varid)

       CALL nf_inq_varndims_wrap(ncid, varid, ndims)
       IF (ndims /= 1) THEN
          CALL msg_write(sub_name, 'illegal ndims for ', dim_edge_name, ndims)
          STOP
       END IF

       CALL nf_inq_vardimid_wrap(ncid, varid, vardims)
       CALL nf_inq_dimlen_wrap(ncid, vardims(1), dimlenp1)
       IF (dimlenp1 /= dimlen + 1) THEN
          CALL msg_write(sub_name, 'invalid size for dimension ', &
               dim_edge_name, dimlenp1)
          STOP
       END IF

       CALL nf_get_var_wrap(ncid, varid, dim_edge_vals)
    END IF

  END SUBROUTINE nf_read_dim_edges

  !*****************************************************************************

  SUBROUTINE regrid_to_tgrid(msv,src_data)

    !---------------------------------------------------------------------------
    !   arguments
    !---------------------------------------------------------------------------
    REAL(KIND=real_kind), INTENT(IN) :: msv
    REAL(KIND=real_kind), DIMENSION(src_imt, src_jmt), INTENT(INOUT) :: src_data

    !---------------------------------------------------------------------------
    !   local variables
    !---------------------------------------------------------------------------
    REAL(KIND=real_kind), DIMENSION(src_imt, src_jmt) :: work
    REAL(KIND=real_kind), DIMENSION(4) :: corners
    REAL(KIND=real_kind), DIMENSION(4) :: wts
    REAL(KIND=real_kind) :: n
    INTEGER(KIND=int_kind) :: i,j

    work = msv
    DO j=2,src_jmt
      DO i=1,src_imt
         wts = (/ 1., 1., 1., 1. /)
         IF (i==1) THEN
           corners = (/ src_data(i,j),src_data(i,j-1),src_data(src_imt,j), &
                src_data(src_imt,j-1) /)
         ELSE
           corners = (/ src_data(i,j),src_data(i,j-1),src_data(i-1,j), &
                src_data(i-1,j-1) /)
         END IF
         where (corners == msv) wts = 0.
         n = sum(wts)
         IF (n .gt. 0) THEN
            wts = wts/n
            work(i,j) = corners(1)*wts(1) + corners(2)*wts(2) + &
                corners(3)*wts(3) + corners(4)*wts(4)
         END IF
      END DO
    END DO
    src_data = work

  END SUBROUTINE regrid_to_tgrid

  !*****************************************************************************

  SUBROUTINE dst_depth_init(ivar)

    INTEGER(KIND=int_kind), INTENT(IN) :: ivar
    !---------------------------------------------------------------------------
    !   local variables
    !---------------------------------------------------------------------------

    CHARACTER(LEN=*), PARAMETER :: sub_name = 'dst_depth_init'
    REAL(KIND=real_kind) :: ztmp
    INTEGER(KIND=int_kind) :: i

    !---------------------------------------------------------------------------
    !   If there is no depth information in src, do not put any in dst.
    !---------------------------------------------------------------------------

    IF (src_km(ivar) == 0) THEN
       dst_kmin = 0
       dst_kmax = 0
       RETURN
    END IF

    !---------------------------------------------------------------------------
    !   Set destination vertical grid
    !---------------------------------------------------------------------------

    IF (dst_z_standard) THEN
       open(10,file=z_file,err=100)
       read(10,*) dst_km
       IF (.NOT.(ALLOCATED(dst_depth))) ALLOCATE(dst_depth(dst_km))
       DO i=1,dst_km
         read(10,*) ztmp
         dst_depth(i) = ztmp
       END DO
       close(10)
       dst_kmin = 1
       dst_kmax = dst_km
    ELSE
       dst_km = src_km(ivar)
       IF (.NOT.(ALLOCATED(dst_depth))) ALLOCATE(dst_depth(dst_km))
       dst_depth = src_depth
       dst_kmin = 1
       dst_kmax = src_km(ivar)
    END IF
    RETURN
100 CALL msg_write(sub_name, 'Problem opening z_file ',z_file)
    STOP

  END SUBROUTINE dst_depth_init

  !*****************************************************************************

  SUBROUTINE create_dst_file(ivar)
    include 'netcdf.inc'
    !---------------------------------------------------------------------------
    !   arguments
    !---------------------------------------------------------------------------
    INTEGER(KIND=int_kind), INTENT(IN) :: ivar

    !---------------------------------------------------------------------------
    !   local variables
    !---------------------------------------------------------------------------

    LOGICAL(KIND=log_kind) :: dst_exists, lcreate
    INTEGER(KIND=int_kind) :: i, j

    !---------------------------------------------------------------------------
!    WRITE(*,*) 'Entering create_dst_file ',ivar
    dst_needs_time        = .true.
    dst_needs_depth       = dst_kmin > 0
    dst_needs_LAT         = .true.
    dst_needs_LON         = .true.
    dst_needs_AREA        = .true.
    dst_needs_IMASK       = .true.

    inquire(file=dst_file(ivar), exist=dst_exists)
    lcreate = dst_lclobber .OR. .NOT. dst_exists

    IF (lcreate) THEN
       CALL nf_create_wrap(dst_file(ivar), 0, dst_ncid(ivar))
    ELSE
       CALL nf_open_wrap(dst_file(ivar), NF_WRITE, dst_ncid(ivar))

       !------------------------------------------------------------------------
       !   Determine which dimensions & variables are already in the file.
       !------------------------------------------------------------------------

       IF (dst_needs_time) THEN
          CALL nf_dimlookup(dst_ncid(ivar), time, dst_time_name, dst_time_dimid)
          IF (dst_time_dimid > 0) dst_needs_time = .FALSE.
       END IF

       IF (dst_needs_depth) THEN
          CALL nf_dimlookup(dst_ncid(ivar), dst_depth(dst_kmin:dst_kmax), &
               dst_depth_name, dst_depth_dimid)
          IF (dst_depth_dimid > 0) dst_needs_depth = .FALSE.
       END IF

       CALL nf_varlookup(dst_ncid(ivar), LAT, dst_LAT_name, dst_LAT_varid)
       IF (dst_LAT_varid > 0) dst_needs_LAT = .FALSE.

       CALL nf_varlookup(dst_ncid(ivar), LON, dst_LON_name, dst_LON_varid)
       IF (dst_LON_varid > 0) dst_needs_LON = .FALSE.

       CALL nf_varlookup(dst_ncid(ivar), AREA, dst_AREA_name, dst_AREA_varid)
       
!       IF (dst_AREA_varid > 0) dst_needs_AREA = .FALSE.
       dst_needs_AREA = .FALSE.
       CALL nf_varlookup(dst_ncid(ivar), REAL(IMASK, real_kind), &
               dst_IMASK_name, dst_IMASK_varid)
!       IF (dst_IMASK_varid > 0) dst_needs_IMASK = .FALSE.
       dst_needs_IMASK = .FALSE.

       CALL nf_redef_wrap(dst_ncid(ivar))
    END IF
!    WRITE(*,*) 'Exiting  create_dst_file ',ivar

  END SUBROUTINE create_dst_file

  !*****************************************************************************

  SUBROUTINE def_dst_file(ivar)

    !---------------------------------------------------------------------------
    !   arguments
    !---------------------------------------------------------------------------
    INTEGER(KIND=int_kind), INTENT(IN) :: ivar

    !---------------------------------------------------------------------------
    !   define dimensions, variables and attributes
    !---------------------------------------------------------------------------

    INCLUDE 'netcdf.inc'

    !---------------------------------------------------------------------------
    !   local variables
    !---------------------------------------------------------------------------

    INTEGER(KIND=int_kind) :: tmp_varid, src_timeid, natts, attnum, tmp_dimid
    CHARACTER(LEN=char_len) :: att_name, att_string
    INTEGER(KIND=int_kind) :: stat
    REAL(KIND=real_kind) :: msv

    !---------------------------------------------------------------------------
!    WRITE(*,*) 'Entering def_dst_file ',ivar
    IF (dst_needs_time) THEN
       dst_time_name = src_time_name
       CALL nf_def_dim_wrap(dst_ncid(ivar), dst_time_name, NF_UNLIMITED, dst_time_dimid)
       CALL nf_def_var_wrap(dst_ncid(ivar), dst_time_name, NF_FLOAT, 1, &
            (/ dst_time_dimid /), tmp_varid)

       !------------------------------------------------------------------------
       !   copy attributes for src_time
       !------------------------------------------------------------------------

       CALL nf_inq_varid_wrap(src_ncid, src_time_name, src_timeid)
       CALL nf_inq_varnatts_wrap(src_ncid, src_timeid, natts)
       DO attnum = 1,natts
          CALL nf_inq_attname_wrap(src_ncid, src_timeid, attnum, att_name)
          if (trim(att_name) /= "bounds") then
             CALL nf_copy_att_wrap(src_ncid, src_timeid, att_name, dst_ncid(ivar), &
                  tmp_varid)
          END IF
       END DO
    END IF

    IF (dst_needs_depth) THEN
       dst_depth_name = src_depth_name
       CALL nf_def_dim_wrap(dst_ncid(ivar), dst_depth_name, dst_kmax-dst_kmin+1, &
            dst_depth_dimid)

       CALL nf_def_var_wrap(dst_ncid(ivar), dst_depth_name, NF_FLOAT, 1, &
            (/ dst_depth_dimid /), tmp_varid)

       att_name = 'units'
       att_string = 'meters'
       CALL nf_put_att_wrap(dst_ncid(ivar), tmp_varid, att_name, &
            LEN_TRIM(att_string), att_string)

       att_name = 'positive'
       att_string = 'down'
       CALL nf_put_att_wrap(dst_ncid(ivar), tmp_varid, att_name, &
            LEN_TRIM(att_string), att_string)
    END IF

    IF (dst_needs_LAT) THEN
       dst_LAT_name = 'lat'
       dst_X_name = 'lon'
       dst_Y_name = 'lat'
       CALL nf_def_dim_wrap(dst_ncid(ivar), dst_X_name, dst_imt, dst_X_dimid)
       CALL nf_def_dim_wrap(dst_ncid(ivar), dst_Y_name, dst_jmt, dst_Y_dimid)
!       CALL nf_def_var_wrap(dst_ncid(ivar), dst_LAT_name, NF_FLOAT, 2, &
!            (/ dst_X_dimid, dst_Y_dimid /), dst_LAT_varid)
       CALL nf_def_var_wrap(dst_ncid(ivar), dst_LAT_name, NF_FLOAT, 1, &
            (/ dst_Y_dimid /), dst_LAT_varid)

       att_name = 'units'
       att_string = 'degrees_north'
       CALL nf_put_att_wrap(dst_ncid(ivar), dst_LAT_varid, att_name, &
            LEN_TRIM(att_string), att_string)

       att_name = 'long_name'
       att_string = 'Latitude'
       CALL nf_put_att_wrap(dst_ncid(ivar), dst_LAT_varid, att_name, &
            LEN_TRIM(att_string), att_string)
    END IF

    IF (dst_needs_LON) THEN
       dst_LON_name = 'lon'
!       CALL nf_def_var_wrap(dst_ncid(ivar), dst_LON_name, NF_FLOAT, 2, &
!            (/ dst_X_dimid, dst_Y_dimid /), dst_LON_varid)
       CALL nf_def_var_wrap(dst_ncid(ivar), dst_LON_name, NF_FLOAT, 1, &
            (/ dst_X_dimid /), dst_LON_varid)

       att_name = 'units'
       att_string = 'degrees_east'
       CALL nf_put_att_wrap(dst_ncid(ivar), dst_LON_varid, att_name, &
            LEN_TRIM(att_string), att_string)

       att_name = 'long_name'
       att_string = 'Longitude'
       CALL nf_put_att_wrap(dst_ncid(ivar), dst_LON_varid, att_name, &
            LEN_TRIM(att_string), att_string)
    END IF

    !------------------------------------------------------------------------
    !   define dst_var
    !------------------------------------------------------------------------

    IF (dst_kmin == 0 .AND. timelen == 0) THEN
       CALL nf_def_var_wrap(dst_ncid(ivar), dst_var(ivar), NF_FLOAT, 2, &
            (/ dst_X_dimid, dst_Y_dimid /), dst_varid(ivar))
    ELSE IF (dst_kmin > 0 .AND. timelen == 0) THEN
       CALL nf_def_var_wrap(dst_ncid(ivar), dst_var(ivar), NF_FLOAT, 3, &
            (/ dst_X_dimid, dst_Y_dimid, dst_depth_dimid /), dst_varid(ivar))
    ELSE IF (dst_kmin == 0 .AND. timelen > 0) THEN
       CALL nf_def_var_wrap(dst_ncid(ivar), dst_var(ivar), NF_FLOAT, 3, &
            (/ dst_X_dimid, dst_Y_dimid, dst_time_dimid /), dst_varid(ivar))
    ELSE
       CALL nf_def_var_wrap(dst_ncid(ivar), dst_var(ivar), NF_FLOAT, 4, &
            (/ dst_X_dimid, dst_Y_dimid, dst_depth_dimid, dst_time_dimid /),&
            dst_varid(ivar))
    END IF

    !------------------------------------------------------------------------
    !   copy attributes for src_varid
    !------------------------------------------------------------------------

    CALL nf_inq_varnatts_wrap(src_ncid, src_varid(ivar), natts)
    DO attnum = 1,natts
       CALL nf_inq_attname_wrap(src_ncid, src_varid(ivar), attnum, att_name)
       CALL nf_copy_att_wrap(src_ncid, src_varid(ivar), att_name, dst_ncid(ivar), dst_varid(ivar))
    END DO

    CALL nf_get_att_wrap(dst_ncid(ivar), dst_varid(ivar), 'missing_value', msv, &
         ALLOW=NF_ENOTATT, stat_out=stat)
    IF (stat == NF_ENOTATT) THEN
       att_name = 'missing_value'
       CALL nf_put_att_wrap(dst_ncid(ivar), dst_varid(ivar), att_name, NF_FLOAT, 1, &
            default_msv)
    END IF
!    WRITE(*,*) 'Exiting  def_dst_file ',ivar

  END SUBROUTINE def_dst_file

  !*****************************************************************************

  SUBROUTINE put_aux_vars(ivar)
    !---------------------------------------------------------------------------
    !   arguments
    !---------------------------------------------------------------------------
    INTEGER(KIND=int_kind), INTENT(IN) :: ivar

    !---------------------------------------------------------------------------
    !   local variables
    !---------------------------------------------------------------------------

    REAL   (KIND=real_kind) :: time_var
    INTEGER(KIND=int_kind)  :: varid

    !---------------------------------------------------------------------------

    IF (dst_needs_time) THEN
       CALL nf_inq_varid_wrap(dst_ncid(ivar), dst_time_name, varid)
       CALL nf_put_vara_wrap(dst_ncid(ivar), varid, (/1/), (/timelen/), time )
!!       write(*,*) 'nf_inq_varid_wrap: ', trim(dst_time_name),' ',varid
!!       write(*,*) 'nf_put_var_wrap  : ', varid, timelen, time
    END IF

    IF (dst_kmin > 0 .AND. dst_needs_depth) THEN
       CALL nf_inq_varid_wrap(dst_ncid(ivar), dst_depth_name, varid)
       CALL nf_put_var_wrap(dst_ncid(ivar), varid, dst_depth(dst_kmin:dst_kmax))
    END IF

    IF (dst_needs_LAT) THEN
       IF (.NOT.ALLOCATED(LAT)) ALLOCATE(LAT(dst_jmt))
       LAT = dst_grid_center_lat(1:dst_imt*dst_jmt:dst_imt)
       CALL nf_put_var_wrap(dst_ncid(ivar), dst_LAT_varid, LAT)
    END IF

    IF (dst_needs_LON) THEN
       IF (.NOT.ALLOCATED(LON)) ALLOCATE(LON(dst_imt))
       LON = dst_grid_center_lon(1:dst_imt)
       CALL nf_put_var_wrap(dst_ncid(ivar), dst_LON_varid, LON)
    END IF

  END SUBROUTINE put_aux_vars

  !*****************************************************************************

  SUBROUTINE regrid_wrap(k,l,ivar)

    !---------------------------------------------------------------------------
    !   arguments
    !---------------------------------------------------------------------------

    INTEGER(KIND=int_kind), INTENT(IN) :: k    ! depth level
    INTEGER(KIND=int_kind), INTENT(IN) :: l    ! time level
    INTEGER(KIND=int_kind), INTENT(IN) :: ivar ! variable number

    !---------------------------------------------------------------------------
    !   local variables
    !---------------------------------------------------------------------------

    CHARACTER(LEN=*), PARAMETER :: sub_name = 'regrid_wrap'
    INTEGER(KIND=int_kind) :: stat
    REAL(KIND=real_kind) :: msv
    REAL(KIND=real_kind), DIMENSION(src_imt, src_jmt) :: src_data,src_data2
    REAL(KIND=real_kind), DIMENSION(dst_imt, dst_jmt) :: &
       src_data_above_on_dst, src_data_below_on_dst, dst_data, dst_data_prev
    INTEGER(KIND=int_kind) :: dst_lev, src_lev, i, j
    REAL(KIND=real_kind) :: dst_z, src_z, src_z_above, src_z_below, wz, &
       num, denom
    LOGICAL(KIND=log_kind) :: needs_fill

    !---------------------------------------------------------------------------
    !   use unlikely value if variable has no missing_value attribute
    !---------------------------------------------------------------------------
!    WRITE(*,*) 'Entering ',trim(sub_name),' ',trim(src_var(ivar)),' T: ',l,' Z: ',k
    !
    CALL nf_get_att_wrap(src_ncid, src_varid(ivar), 'missing_value', msv, &
         ALLOW=NF_ENOTATT, stat_out=stat)
    IF (stat == NF_ENOTATT) msv = default_msv
    CALL nf_get_vara_wrap(src_ncid, src_varid(ivar), (/ 1, 1, l /), &
         (/ src_imt, src_jmt, 1 /), src_data)
    IF (src_var_rotate(ivar)) THEN
       CALL nf_get_vara_wrap(src_ncid, src_varid2(ivar), (/ 1, 1, l /), &
            (/ src_imt, src_jmt, 1 /), src_data2)
       CALL rotate(msv,src_data,src_data2,ivar)
    END IF
    CALL regrid_level(msv, src_data, dst_data)
    !---------------------------------------------------------------------------
    !   Mask out areas where aice < 1%
    !---------------------------------------------------------------------------
    WHERE (area_mask == 0.)
       dst_data = 0.
    ENDWHERE
    ! Set zero velocities to msv
    IF ((src_var(ivar) == 'u').or.(src_var(ivar) == 'v')) THEN
       WHERE (dst_data == 0.)
          dst_data = msv
       ENDWHERE
    END IF
    CALL nf_put_vara_wrap(dst_ncid(ivar), dst_varid(ivar), (/ 1, 1, l /), &
         (/ dst_imt, dst_jmt, 1 /), dst_data)
!    WRITE(*,*) 'Exiting  ',trim(sub_name)
  END SUBROUTINE regrid_wrap

  !*****************************************************************************

  SUBROUTINE regrid_level(msv, src_data, dst_data)
    !---------------------------------------------------------------------------
    !   arguments
    REAL(KIND=real_kind), INTENT(IN) :: msv
    REAL(KIND=real_kind), DIMENSION(src_imt, src_jmt), INTENT(INOUT) :: src_data
    REAL(KIND=real_kind), DIMENSION(dst_imt, dst_jmt), INTENT(OUT)   :: dst_data
    !---------------------------------------------------------------------------
    !   local variables
    INTEGER(KIND=int_kind) :: in,jn,kn
    REAL(KIND=real_kind), DIMENSION(src_grid_size) :: tmp1
    REAL(KIND=real_kind), DIMENSION(dst_grid_size) :: tmp2
    !
!    WRITE(*,*) 'Entering regrid_level'
    tmp1 = reshape(src_data,(/src_grid_size/))
    tmp2 = msv
    DO in = 1, num_links
      IF (src_address(in) /= 0 .AND. tmp1(src_address(in)) /= msv) THEN
         IF (tmp2(dst_address(in)) == msv) THEN
            tmp2(dst_address(in)) = tmp1(src_address(in))*remap_matrix(in)
         ELSE
            tmp2(dst_address(in)) = tmp2(dst_address(in)) + &
                 tmp1(src_address(in))*remap_matrix(in)
         END IF
      END IF
    END DO
    dst_data = reshape(tmp2,(/ dst_imt,dst_jmt /))
!    WRITE(*,*) 'Exiting  regrid_level'
  END SUBROUTINE regrid_level

  !*****************************************************************************

  SUBROUTINE rotate(msv, data1, data2, ivar)

    !---------------------------------------------------------------------------
    !   arguments
    !---------------------------------------------------------------------------

    REAL(KIND=real_kind), INTENT(IN) :: msv
    REAL(KIND=real_kind), DIMENSION(src_imt, src_jmt), INTENT(INOUT) :: data1
    REAL(KIND=real_kind), DIMENSION(src_imt, src_jmt), INTENT(IN) :: data2
    INTEGER(KIND=int_kind), INTENT(IN) :: ivar

    !---------------------------------------------------------------------------
    !   local variables
    !---------------------------------------------------------------------------

    INTEGER(KIND=int_kind) :: in
    REAL(KIND=real_kind), DIMENSION(src_grid_size) :: tmp1
    REAL(KIND=real_kind), DIMENSION(dst_grid_size) :: tmp2

    SELECT CASE (src_var(ivar))
    CASE ('u')
      where(data1 /= msv) data1 = data1*cos(ANGLE) + data2*sin(-ANGLE)
    CASE ('v')
      where(data1 /= msv) data1 = data1*cos(ANGLE) - data2*sin(-ANGLE)
    CASE DEFAULT
      WRITE(*,*) trim(src_var(ivar)),' ROTATION NOT FOUND. ERR!'
      STOP
    END SELECT

  END SUBROUTINE rotate

  !*****************************************************************************

  SUBROUTINE compute_area(x,y,a)

    !---------------------------------------------------------------------------
    !   arguments
    !---------------------------------------------------------------------------
    REAL(KIND=real_kind), DIMENSION(dst_imt,dst_jmt), INTENT(IN) :: x,y
    REAL(KIND=real_kind), DIMENSION(dst_imt,dst_jmt), INTENT(OUT) :: a

    !---------------------------------------------------------------------------
    !   local variables
    !---------------------------------------------------------------------------
    INTEGER(KIND=int_kind) :: imt, jmt, i, j, ip1, im1, jp1, jm1
    REAL(KIND=real_kind), DIMENSION(dst_imt,dst_jmt) :: dxt,dyt

    ! For future development... add some way of computing destination
    ! grid area weights if they are not available from the scrip
    ! remapping file (scrip only includes areas for conservative
    ! mappings).

  END SUBROUTINE compute_area

  !*****************************************************************************

  INTEGER FUNCTION find_index(x, sorted_array)

    !---------------------------------------------------------------------------
    !   arguments
    !---------------------------------------------------------------------------

    REAL(KIND=real_kind) :: x
    REAL(KIND=real_kind), DIMENSION(:) :: sorted_array

    !---------------------------------------------------------------------------
    !   local variables
    !---------------------------------------------------------------------------

    INTEGER(KIND=int_kind) :: len, i

    len = SIZE(sorted_array)

    IF (x < sorted_array(1)) THEN
       find_index = 0
       RETURN
    END IF

    DO i = 1,len-1
       IF (x < sorted_array(i+1)) THEN
          find_index = i
          RETURN
       END IF
    END DO

    find_index = len

  END FUNCTION find_index

  !*****************************************************************************

  SUBROUTINE free_vars

    IF (ALLOCATED(LAT))             DEALLOCATE(LAT)
    IF (ALLOCATED(LON))             DEALLOCATE(LON)
    IF (ALLOCATED(AREA))            DEALLOCATE(AREA)
    IF (ALLOCATED(IMASK))           DEALLOCATE(IMASK)
    IF (ALLOCATED(ANGLE))           DEALLOCATE(ANGLE)
    IF (ALLOCATED(src_depth))       DEALLOCATE(src_depth)
    IF (ALLOCATED(dst_depth))       DEALLOCATE(dst_depth)
    IF (ALLOCATED(time))            DEALLOCATE(time)
    IF (ALLOCATED(src_grid_center_lon))  DEALLOCATE(src_grid_center_lon)
    IF (ALLOCATED(src_grid_center_lat))  DEALLOCATE(src_grid_center_lat)
    IF (ALLOCATED(src_grid_corner_lon))  DEALLOCATE(src_grid_corner_lon)
    IF (ALLOCATED(src_grid_corner_lat))  DEALLOCATE(src_grid_corner_lat)
    IF (ALLOCATED(src_grid_area))  DEALLOCATE(src_grid_area)
    IF (ALLOCATED(src_grid_imask))  DEALLOCATE(src_grid_imask)
    IF (ALLOCATED(src_grid_frac))  DEALLOCATE(src_grid_frac)
    IF (ALLOCATED(dst_grid_center_lon))  DEALLOCATE(dst_grid_center_lon)
    IF (ALLOCATED(dst_grid_center_lat))  DEALLOCATE(dst_grid_center_lat)
    IF (ALLOCATED(dst_grid_corner_lon))  DEALLOCATE(dst_grid_corner_lon)
    IF (ALLOCATED(dst_grid_corner_lat))  DEALLOCATE(dst_grid_corner_lat)
    IF (ALLOCATED(dst_grid_area))  DEALLOCATE(dst_grid_area)
    IF (ALLOCATED(dst_grid_imask))  DEALLOCATE(dst_grid_imask)
    IF (ALLOCATED(dst_grid_frac))  DEALLOCATE(dst_grid_frac)
    IF (ALLOCATED(src_address))  DEALLOCATE(src_address)
    IF (ALLOCATED(dst_address))  DEALLOCATE(dst_address)
    IF (ALLOCATED(remap_matrix))  DEALLOCATE(remap_matrix)

  END SUBROUTINE free_vars

  !*****************************************************************************

END MODULE ccsm_ice_rgd

