!
subroutine open_cdf(ncid,filename,read_only)
  !
  use counters_netcdf_jfl
  use definitions_netcdf_jfl
  !
  ! open existing netcdf file and return its ID number
  !
  implicit none
  !
  include 'netcdf.inc'
  !
  ! arguments
  !
  integer ncid
  character(*) filename
  logical, optional :: read_only
  !
  ! local
  !
  integer status
  !
  if(present(read_only)) then
     if(read_only) then
        status = nf_open(filename,NF_NOWRITE,ncid)
        if(status.ne.nf_noerr) call handle_err(status)
     else
        status = nf_open(filename,NF_WRITE,ncid)
        if(status.ne.nf_noerr) call handle_err(status)
     endif
  else
     status = nf_open(filename,NF_WRITE,ncid)
     if(status.ne.nf_noerr) call handle_err(status)
  endif
  !
  ! make sure the file ID is OK
  !
  if(ncid.lt.0) then
     print *,'Problem opening ',filename,' ncid = ',ncid
     call abort()
  else
     file_counter = file_counter + 1
  endif
  !
  nc_filename(file_counter) = filename
  !
  return
end subroutine open_cdf
!
subroutine create_cdf(ncid,filename)
  !
  use counters_netcdf_jfl
  use definitions_netcdf_jfl
  !
  ! create new netcdf file and return its ID number
  !
  implicit none
  !
  include 'netcdf.inc'
  !
  ! arguments
  !
  integer ncid
  character(*) filename
  !
  ! local
  !
  integer status
  !
  status = nf_create(filename,NF_CLOBBER,ncid)
  if(status.ne.nf_noerr) call handle_err(status)
  !
  ! make sure the file ID is OK
  !
  if(ncid.le.0) then
     print *,'Problem opening ',filename,' ncid = ',ncid
     call abort()
  else
     file_counter = file_counter + 1
  endif
  !
  nc_filename(file_counter) = filename
  !
  return
end subroutine create_cdf
!
subroutine end_def(ncid)
  !
  ! end definition mode for netcdf file
  !
  implicit none
  !
  include 'netcdf.inc'
  !
  ! argument
  !
  integer ncid
  !
  ! local
  !
  integer status
  !
  status = nf_enddef(ncid)
  if(status.ne.nf_noerr) call handle_err(status)
  !
  return
end subroutine end_def
!
subroutine close_cdf(ncid)
  !
  ! create new netcdf file and return its ID number
  !
  implicit none
  !
  include 'netcdf.inc'
  !
  ! arguments
  !
  integer ncid
  !
  ! local
  !
  integer status
  !
  status = nf_close(ncid)
  if(status.ne.nf_noerr) call handle_err(status)
  !
  return
end subroutine close_cdf
!
subroutine def_dim(ncid,dim_name,dim_length)
  !
  ! define dimension in netcdf file
  !
  use definitions_netcdf_jfl
  use counters_netcdf_jfl
  !
  implicit none
  !
  include 'netcdf.inc'
  !
  ! arguments
  !
  integer ncid,dim_length
  character(*) dim_name
  !
  ! local
  !
  integer status,dim_id
  !
  if(dim_length.gt.0) then
     status = nf_def_dim(ncid,dim_name,dim_length  ,dim_id)
  else
     status = nf_def_dim(ncid,dim_name,nf_unlimited,dim_id)
  endif
  if(status.ne.nf_noerr) call handle_err(status)
  !
  dim_counter = dim_counter + 1
  if(dim_counter.gt.dim_max) then
     print *,'maximum number of dimensions ',dim_max
     print *,'is reached, stopping'
     stop
  endif
  dim_info(dim_counter)%name      = dim_name
  dim_info(dim_counter)%id        = dim_id
  dim_info(dim_counter)%length    = dim_length
  dim_info(dim_counter)%file_name = nc_filename(file_counter)
  !
  return
end subroutine def_dim
!
subroutine def_var(ncid,var_name,format,type,units,long_name)
  !
  use definitions_netcdf_jfl
  use counters_netcdf_jfl
  !
  implicit none
  !
  include 'netcdf.inc'
  !
  ! arguments
  !
  integer ncid
  character(*)::var_name,format,type
  character(*),optional :: units,long_name
  !
  ! local
  !
  integer i,j,n
  integer status,xtype,nvdims,var_id
  integer, dimension(:), allocatable :: vdims ,length
  character*1  :: comma = ','
  character*32 :: work
  integer pos1,pos2,pos3,position(100),ncomma
  !
  ! remove possible comma at the end of format
  !
  if(format(len(format):len(format)).eq.comma) then
     format(len(format):len(format)) = ' '
  endif
  !
  ! special case for mif variables
  !
  if(format(:len_trim(format)).eq.'mif') then
     if(type.eq.'character') then
        nvdims = 3
        allocate(vdims (nvdims))
        allocate(length(nvdims))
        do i=1,dim_counter
           if(dim_info(i)%name.eq.'mif1') then
              vdims (2) = dim_info(i)%id
              length(2) = dim_info(i)%length
           elseif(dim_info(i)%name.eq.'mif2') then
              vdims (3) = dim_info(i)%id
              length(3) = dim_info(i)%length
           elseif(dim_info(i)%name.eq.'mif3') then
              vdims (1) = dim_info(i)%id
              length(1) = dim_info(i)%length
           endif
        end do
     else
        nvdims = 2
        allocate(vdims (nvdims))
        allocate(length(nvdims))
        do i=1,dim_counter
           if(dim_info(i)%name.eq.'mif1') then
              vdims (1) = dim_info(i)%id
              length(1) = dim_info(i)%length
           elseif(dim_info(i)%name.eq.'mif2') then
              vdims (2) = dim_info(i)%id
              length(2) = dim_info(i)%length
           endif
        end do
     endif
  else
     !
     ! default case
     !
     if(len_trim(format).eq.0) then
        print *,'problem with defining ',var_name
        print *,'no format defined'
        stop
     endif
     !
     ! find position of commas
     !
     ncomma = 0
     do n=1,len_trim(format)
        if(format(n:n).eq.comma) then
           ncomma = ncomma + 1
           position(ncomma) = n
        endif
     end do
     ncomma = ncomma + 1
     position(ncomma) = len_trim(format) + 1
     !
     nvdims = ncomma
     allocate(vdims (nvdims))
     allocate(length(nvdims))
     vdims = -1
     !
     pos1 = 1
     do i=1,nvdims
        pos2 = position(i) - 1
        pos3 = pos2 - pos1 + 1
        work = ' '
        work(:pos3) = format(pos1:pos2)
        pos1 = pos2 + 2
        do j=1,dim_counter
           if(dim_info(j)%name.eq.work(:pos3).and.&
                dim_info(j)%file_name.eq.nc_filename(file_counter)) then
              vdims (i) = dim_info(j)%id
              length(i) = dim_info(j)%length
              exit
           endif
        end do
     end do
     !
  endif
  !
  ! define external type
  !
  if(type.eq.'integer'  ) xtype = nf_int
  if(type.eq.'float'    ) xtype = nf_float
  if(type.eq.'double'   ) xtype = nf_double
  if(type.eq.'character') xtype = nf_char
  !
  ! define variable
  !
  status = nf_def_var(ncid,var_name,xtype,nvdims,vdims,var_id)
  if(status.ne.nf_noerr) call handle_err(status)
  !
  ! declare attributes
  !
  if(present(long_name)) then
     status = nf_put_att_text(ncid,var_id,'long_name'&
          ,len(long_name),long_name)
     if(status.ne.nf_noerr) call handle_err(status)
  endif
  !
  if(present(units)) then
     status = nf_put_att_text(ncid,var_id,'units'&
          ,len(units),units)
     if(status.ne.nf_noerr) call handle_err(status)
  endif
  !
  ! store information
  !
  var_counter = var_counter + 1
  if(var_counter.gt.var_max) then
     print *,'maximum number of variables ',var_max
     print *,'is reached, stopping'
     stop
  endif
  var_info(var_counter)%name      = var_name
  var_info(var_counter)%format    = format
  var_info(var_counter)%type      = type
  var_info(var_counter)%id        = var_id
  var_info(var_counter)%file_name = nc_filename(file_counter)
  var_info(var_counter)%nvdims    = nvdims
  do i=1,nvdims
     var_info(var_counter)%vdims(i) = length(i)
  end do
  !
  deallocate(vdims )
  deallocate(length)
  !
  return
end subroutine def_var
!
subroutine get_dims(ncid)
  !
  use definitions_netcdf_jfl
  use counters_netcdf_jfl
  !
  implicit none
  !
  include 'netcdf.inc'
  !
  ! arguments
  !
  integer ncid
  !
  ! local
  !
  integer n,num_dim
  integer length
  integer status
  character*32 name
  !
  status = nf_inq_ndims(ncid,num_dim)
  if(status.ne.nf_noerr) call handle_err(status)
  !
  if(num_dim.gt.dim_max) then
     print *,'number of dimensions found = ',num_dim
     print *,'maximum allowed            = ',dim_max
     stop
  endif
  !
  do n=1,num_dim
     name = ' '
     status = nf_inq_dim(ncid,n,name,length)
     if(status.ne.nf_noerr) call handle_err(status)
     dim_counter = dim_counter + 1
     dim_info(dim_counter)%id        = n
     dim_info(dim_counter)%name      = name
     dim_info(dim_counter)%length    = length
     dim_info(dim_counter)%file_name = nc_filename(file_counter)
     dim_info(dim_counter)%ncid      = ncid
  end do
  !
  return
end subroutine get_dims
!
subroutine get_vars(ncid)
  !
  use definitions_netcdf_jfl
  use counters_netcdf_jfl
  !
  implicit none
  !
  include 'netcdf.inc'
  !
  ! arguments
  !
  integer ncid
  !
  ! local
  !
  integer i,j,n,num_var
  integer status
  integer xtype,nvdims,length,natts
  integer pos1,pos2
  integer, dimension(:), allocatable :: dim_id
  real::missing_value,FillValue
  integer::int_missing_value,int_FillValue
  character(len=256)::name,aname,units
  character*1 :: comma = ','
  !
  status = nf_inq_nvars(ncid,num_var)
  if(status.ne.nf_noerr) call handle_err(status)
  !
  if(num_var.gt.var_max) then
     print *,'number of variables found = ',num_var
     print *,'maximum allowed           = ',var_max
     stop
  endif
  !
  do n=1,num_var
     !
     name          = ' '
     units         = ' '
     missing_value = -1.e32
     FillValue     = -1.e32
     int_missing_value = -1.e9
     int_FillValue     = -1.e9
     status = nf_inq_varname(ncid,n,name)
     if(status.ne.nf_noerr) call handle_err(status)
     !
     status = nf_inq_vartype(ncid,n,xtype)
     if(status.ne.nf_noerr) call handle_err(status)
     !
     status = nf_inq_varndims(ncid,n,nvdims)
     if(status.ne.nf_noerr) call handle_err(status)
     !
     status = nf_inq_varnatts(ncid,n,natts)
     if(status.ne.nf_noerr) call handle_err(status)
     if (natts /= 0) then
        do i=1,natts
           status=nf_inq_attname(ncid,n,i,aname)
           if(status.ne.nf_noerr) call handle_err(status)
           if (trim(aname) == 'units') then
              status = nf_get_att_text(ncid,n,aname,units)
              if(status.ne.nf_noerr) call handle_err(status)
           endif
           if (trim(aname) == 'missing_value') then
              if (xtype == nf_float) then
                 status = nf_get_att_real(ncid,n,aname,missing_value)
                 if(status.ne.nf_noerr) call handle_err(status)
              elseif (xtype == nf_int) then
                 status = nf_get_att_int(ncid,n,aname,int_missing_value)
                 if(status.ne.nf_noerr) call handle_err(status)
              endif
           endif
           if (trim(aname) == '_FillValue') then
              if (xtype == nf_float) then
                 status = nf_get_att_real(ncid,n,aname,FillValue)
                 if(status.ne.nf_noerr) call handle_err(status)
              elseif (xtype == nf_int) then
                 status = nf_get_att_int(ncid,n,aname,int_FillValue)
                 if(status.ne.nf_noerr) call handle_err(status)
              endif
           endif
        enddo
     endif
     !
     allocate(dim_id(nvdims))
     status = nf_inq_vardimid(ncid,n,dim_id)
     if(status.ne.nf_noerr) call handle_err(status)
     !
     var_counter = var_counter + 1
     var_info(var_counter)%ncid              = ncid
     var_info(var_counter)%id                = n
     var_info(var_counter)%name              = name
     var_info(var_counter)%units             = units
     var_info(var_counter)%missing_value     = missing_value
     var_info(var_counter)%FillValue         = FillValue
     var_info(var_counter)%int_missing_value = int_missing_value
     var_info(var_counter)%int_FillValue     = int_FillValue
     var_info(var_counter)%file_name         = nc_filename(file_counter)
     select case ( xtype )
     case ( nf_int )
        var_info(var_counter)%type = 'integer'
     case ( nf_float )
        var_info(var_counter)%type = 'float'
     case ( nf_double )
        var_info(var_counter)%type = 'double'
     case ( nf_char )
        var_info(var_counter)%type = 'character'
     end select
     var_info(var_counter)%nvdims = nvdims
     pos1 = 1
     pos2 = 1
     do i=1,nvdims
        do j=1,dim_counter
           if((dim_info(j)%id.eq.dim_id(i)).and.&
                (dim_info(j)%file_name.eq.nc_filename(file_counter))) then
              length = dim_info(j)%length
              name = ' '
              name = dim_info(j)%name
              exit
           endif
        end do
        var_info(var_counter)%vdims(i) = length
        pos2 = pos1 + len_trim(name) - 1
        var_info(var_counter)%format(pos1:pos2)&
             = name(:len_trim(name))
        pos2 = pos2 + 1
        var_info(var_counter)%format(pos2:pos2) = comma
        pos1 = pos2 + 1
     end do
     !
     ! remove last comma
     !
     var_info(var_counter)%format(pos2:pos2) = ' '
     !
     ! deallocate memory
     !
     deallocate(dim_id)
     !
  end do
  !
  return
end subroutine get_vars
!
subroutine write_var_i_0(ncid,var_name,array)
  !
  use definitions_netcdf_jfl
  use counters_netcdf_jfl
  !
  implicit none
  !
  include 'netcdf.inc'
  !
  ! arguments
  !
  integer ncid
  character(*) var_name
  integer array
  !
  include 'put_var_int.inc'
  !
  return
end subroutine write_var_i_0
!
subroutine write_var_i_1(ncid,var_name,array)
  !
  use definitions_netcdf_jfl
  use counters_netcdf_jfl
  !
  implicit none
  !
  include 'netcdf.inc'
  !
  ! arguments
  !
  integer ncid
  character(*) var_name
  integer array(:)
  !
  include 'put_var_int.inc'
  !
  return
end subroutine write_var_i_1
!
subroutine write_var_i_2(ncid,var_name,array)
  !
  use definitions_netcdf_jfl
  use counters_netcdf_jfl
  !
  implicit none
  !
  include 'netcdf.inc'
  !
  ! arguments
  !
  integer ncid
  character(*) var_name
  integer array(:,:)
  !
  include 'put_var_int.inc'
  !
  return
end subroutine write_var_i_2
!
subroutine write_var_i_3(ncid,var_name,array)
  !
  use definitions_netcdf_jfl
  use counters_netcdf_jfl
  !
  implicit none
  !
  include 'netcdf.inc'
  !
  ! arguments
  !
  integer ncid
  character(*) var_name
  integer array(:,:,:)
  !
  include 'put_var_int.inc'
  !
  return
end subroutine write_var_i_3
!
subroutine write_var_f_0(ncid,var_name,array)
  !
  use definitions_netcdf_jfl
  use counters_netcdf_jfl
  !
  implicit none
  !
  include 'netcdf.inc'
  !
  ! arguments
  !
  integer ncid
  character(*) var_name
  real array
  !
  include 'put_var_real.inc'
  !
  return
end subroutine write_var_f_0
!
subroutine write_var_f_1(ncid,var_name,array)
  !
  use definitions_netcdf_jfl
  use counters_netcdf_jfl
  !
  implicit none
  !
  include 'netcdf.inc'
  !
  ! arguments
  !
  integer ncid
  character(*) var_name
  real array(:)
  !
  include 'put_var_real.inc'
  !
  return
end subroutine write_var_f_1
!
subroutine write_var_f_2(ncid,var_name,array)
  !
  use definitions_netcdf_jfl
  use counters_netcdf_jfl
  !
  implicit none
  !
  include 'netcdf.inc'
  !
  ! arguments
  !
  integer ncid
  character(*) var_name
  real array(:,:)
  !
  include 'put_var_real.inc'
  !
  return
end subroutine write_var_f_2
!
subroutine write_var_f_3(ncid,var_name,array)
  !
  use definitions_netcdf_jfl
  use counters_netcdf_jfl
  !
  implicit none
  !
  include 'netcdf.inc'
  !
  ! arguments
  !
  integer ncid
  character(*) var_name
  real array(:,:,:)
  !
  include 'put_var_real.inc'
  !
  return
end subroutine write_var_f_3
!
subroutine write_var_f_4(ncid,var_name,array)
  !
  use definitions_netcdf_jfl
  use counters_netcdf_jfl
  !
  implicit none
  !
  include 'netcdf.inc'
  !
  ! arguments
  !
  integer ncid
  character(*) var_name
  real array(:,:,:,:)
  !
  include 'put_var_real.inc'
  !
  return
end subroutine write_var_f_4
!
subroutine write_var_d_0(ncid,var_name,array)
  !
  use definitions_netcdf_jfl
  use counters_netcdf_jfl
  !
  implicit none
  !
  include 'netcdf.inc'
  !
  ! arguments
  !
  integer ncid
  character(*) var_name
  real*8 array
  !
  include 'put_var_dble.inc'
  !
  return
end subroutine write_var_d_0
!
subroutine write_var_d_1(ncid,var_name,array)
  !
  use definitions_netcdf_jfl
  use counters_netcdf_jfl
  !
  implicit none
  !
  include 'netcdf.inc'
  !
  ! arguments
  !
  integer ncid
  character(*) var_name
  real*8 array(:)
  !
  include 'put_var_dble.inc'
  !
  return
end subroutine write_var_d_1
!
subroutine write_var_d_2(ncid,var_name,array)
  !
  use definitions_netcdf_jfl
  use counters_netcdf_jfl
  !
  implicit none
  !
  include 'netcdf.inc'
  !
  ! arguments
  !
  integer ncid
  character(*) var_name
  real*8 array(:,:)
  !
  include 'put_var_dble.inc'
  !
  return
end subroutine write_var_d_2
!
subroutine write_var_d_3(ncid,var_name,array)
  !
  use definitions_netcdf_jfl
  use counters_netcdf_jfl
  !
  implicit none
  !
  include 'netcdf.inc'
  !
  ! arguments
  !
  integer ncid
  character(*) var_name
  real*8 array(:,:,:)
  !
  include 'put_var_dble.inc'
  !
  return
end subroutine write_var_d_3
!
subroutine write_var_c_0(ncid,var_name,array)
  !
  use definitions_netcdf_jfl
  use counters_netcdf_jfl
  !
  implicit none
  !
  include 'netcdf.inc'
  !
  ! arguments
  !
  integer ncid
  character(*) var_name
  character(*) array
  !
  include 'put_var_text.inc'
  !
  return
end subroutine write_var_c_0
!
subroutine write_var_c_1(ncid,var_name,array)
  !
  use definitions_netcdf_jfl
  use counters_netcdf_jfl
  !
  implicit none
  !
  include 'netcdf.inc'
  !
  ! arguments
  !
  integer ncid
  character(*) var_name
  character(*) array(:)
  !
  include 'put_var_text.inc'
  !
  return
end subroutine write_var_c_1
!
subroutine write_var_c_2(ncid,var_name,array)
  !
  use definitions_netcdf_jfl
  use counters_netcdf_jfl
  !
  implicit none
  !
  include 'netcdf.inc'
  !
  ! arguments
  !
  integer ncid
  character(*) var_name
  character(*) array(:,:)
  !
  include 'put_var_text.inc'
  !
  return
end subroutine write_var_c_2
!
subroutine read_var_i_0(ncid,var_name,array)
  !
  use definitions_netcdf_jfl
  use counters_netcdf_jfl
  !
  implicit none
  !
  include 'netcdf.inc'
  !
  ! arguments
  !
  integer ncid
  character(*) var_name
  integer array
  !
  include 'get_var_int.inc'
  !
  return
end subroutine read_var_i_0
!
subroutine read_var_i_1(ncid,var_name,array)
  !
  use definitions_netcdf_jfl
  use counters_netcdf_jfl
  !
  implicit none
  !
  include 'netcdf.inc'
  !
  ! arguments
  !
  integer ncid
  character(*) var_name
  integer array(:)
  !
  include 'get_var_int.inc'
  !
  return
end subroutine read_var_i_1
!
subroutine read_var_i8_1(ncid,var_name,array)
  !
  use definitions_netcdf_jfl
  use counters_netcdf_jfl
  !
  implicit none
  !
  include 'netcdf.inc'
  !
  ! arguments
  !
  integer ncid
  character(*) var_name
  integer*8 array(:)
  !
  include 'get_var_int.inc'
  !
  return
end subroutine read_var_i8_1
!
subroutine read_var_i_2(ncid,var_name,array)
  !
  use definitions_netcdf_jfl
  use counters_netcdf_jfl
  !
  implicit none
  !
  include 'netcdf.inc'
  !
  ! arguments
  !
  integer ncid
  character(*) var_name
  integer array(:,:)
  !
  include 'get_var_int.inc'
  !
  return
end subroutine read_var_i_2
!
subroutine read_var_i_3(ncid,var_name,array)
  !
  use definitions_netcdf_jfl
  use counters_netcdf_jfl
  !
  implicit none
  !
  include 'netcdf.inc'
  !
  ! arguments
  !
  integer ncid
  character(*) var_name
  integer array(:,:,:)
  !
  include 'get_var_int.inc'
  !
  return
end subroutine read_var_i_3
!
subroutine read_var_f_0(ncid,var_name,array)
  !
  use definitions_netcdf_jfl
  use counters_netcdf_jfl
  !
  implicit none
  !
  include 'netcdf.inc'
  !
  ! arguments
  !
  integer ncid
  character(*) var_name
  real array
  !
  include 'get_var_real.inc'
  !
  return
end subroutine read_var_f_0
!
subroutine read_var_f_1(ncid,var_name,array)
  !
  use definitions_netcdf_jfl
  use counters_netcdf_jfl
  !
  implicit none
  !
  include 'netcdf.inc'
  !
  ! arguments
  !
  integer ncid
  character(*) var_name
  real array(:)
  !
  include 'get_var_real.inc'
  !
  return
end subroutine read_var_f_1
!
subroutine read_var_f_2(ncid,var_name,array)
  !
  use definitions_netcdf_jfl
  use counters_netcdf_jfl
  !
  implicit none
  !
  include 'netcdf.inc'
  !
  ! arguments
  !
  integer ncid
  character(*) var_name
  real array(:,:)
  !
  include 'get_var_real.inc'
  !
  return
end subroutine read_var_f_2
!
subroutine read_var_f_3(ncid,var_name,array)
  !
  use definitions_netcdf_jfl
  use counters_netcdf_jfl
  !
  implicit none
  !
  include 'netcdf.inc'
  !
  ! arguments
  !
  integer ncid
  character(*) var_name
  real array(:,:,:)
  !
  include 'get_var_real.inc'
  !
  return
end subroutine read_var_f_3
!
subroutine read_var_f_4(ncid,var_name,array)
  !
  use definitions_netcdf_jfl
  use counters_netcdf_jfl
  !
  implicit none
  !
  include 'netcdf.inc'
  !
  ! arguments
  !
  integer ncid
  character(*) var_name
  real array(:,:,:,:)
  !
  include 'get_var_real.inc'
  !
  return
end subroutine read_var_f_4
!
subroutine read_var_d_0(ncid,var_name,array)
  !
  use definitions_netcdf_jfl
  use counters_netcdf_jfl
  !
  implicit none
  !
  include 'netcdf.inc'
  !
  ! arguments
  !
  integer ncid
  character(*) var_name
  double precision array
  !
  include 'get_var_double.inc'
  !
  return
end subroutine read_var_d_0
!
subroutine read_var_d_1(ncid,var_name,array)
  !
  use definitions_netcdf_jfl
  use counters_netcdf_jfl
  !
  implicit none
  !
  include 'netcdf.inc'
  !
  ! arguments
  !
  integer ncid
  character(*) var_name
  double precision array(:)
  !
  include 'get_var_double.inc'
  !
  return
end subroutine read_var_d_1
!
subroutine read_var_d_2(ncid,var_name,array)
  !
  use definitions_netcdf_jfl
  use counters_netcdf_jfl
  !
  implicit none
  !
  include 'netcdf.inc'
  !
  ! arguments
  !
  integer ncid
  character(*) var_name
  double precision array(:,:)
  !
  include 'get_var_double.inc'
  !
  return
end subroutine read_var_d_2
!
subroutine read_var_d_3(ncid,var_name,array)
  !
  use definitions_netcdf_jfl
  use counters_netcdf_jfl
  !
  implicit none
  !
  include 'netcdf.inc'
  !
  ! arguments
  !
  integer ncid
  character(*) var_name
  double precision array(:,:,:)
  !
  include 'get_var_double.inc'
  !
  return
end subroutine read_var_d_3
!
subroutine read_var_d_4(ncid,var_name,array)
  !
  use definitions_netcdf_jfl
  use counters_netcdf_jfl
  !
  implicit none
  !
  include 'netcdf.inc'
  !
  ! arguments
  !
  integer ncid
  character(*) var_name
  double precision array(:,:,:,:)
  !
  include 'get_var_double.inc'
  !
  return
end subroutine read_var_d_4
!
subroutine read_var_c_0(ncid,var_name,array)
  !
  use definitions_netcdf_jfl
  use counters_netcdf_jfl
  !
  implicit none
  !
  include 'netcdf.inc'
  !
  ! arguments
  !
  integer ncid
  character(*) var_name
  character(*) array
  !
  include 'get_var_text.inc'
  !
  return
end subroutine read_var_c_0
!
subroutine read_var_c_1(ncid,var_name,array)
  !
  use definitions_netcdf_jfl
  use counters_netcdf_jfl
  !
  implicit none
  !
  include 'netcdf.inc'
  !
  ! arguments
  !
  integer ncid
  character(*) var_name
  character(*) array(:)
  !
  include 'get_var_text.inc'
  !
  return
end subroutine read_var_c_1
!
subroutine read_var_c_2(ncid,var_name,array)
  !
  use definitions_netcdf_jfl
  use counters_netcdf_jfl
  !
  implicit none
  !
  include 'netcdf.inc'
  !
  ! arguments
  !
  integer ncid
  character(*) var_name
  character(*) array(:,:)
  !
  include 'get_var_text.inc'
  !
  return
end subroutine read_var_c_2
subroutine write_att_text(ncid,var_name,att_name,att_val)
  !
  use definitions_netcdf_jfl
  use counters_netcdf_jfl
  !
  implicit none
  !
  include 'netcdf.inc'
  !
  ! arguments
  !
  integer ncid
  character(*) var_name
  character(*) att_name
  character(*) att_val
  !
  include 'put_att_text.inc'
  !
  return
end subroutine write_att_text
subroutine write_att_dble(ncid,var_name,att_name,att_val)
  !
  use definitions_netcdf_jfl
  use counters_netcdf_jfl
  !
  implicit none
  !
  include 'netcdf.inc'
  !
  ! arguments
  !
  integer ncid
  character(*) var_name
  character(*) att_name
  real*8       att_val
  !
  include 'put_att_dble.inc'
  !
  return
end subroutine write_att_dble
subroutine write_att_real(ncid,var_name,att_name,att_val)
  !
  use definitions_netcdf_jfl
  use counters_netcdf_jfl
  !
  implicit none
  !
  include 'netcdf.inc'
  !
  ! arguments
  !
  integer ncid
  character(*) var_name
  character(*) att_name
  real         att_val
  !
  include 'put_att_real.inc'
  !
  return
end subroutine write_att_real
subroutine write_att_int(ncid,var_name,att_name,att_val)
  !
  use definitions_netcdf_jfl
  use counters_netcdf_jfl
  !
  implicit none
  !
  include 'netcdf.inc'
  !
  ! arguments
  !
  integer ncid
  character(*) var_name
  character(*) att_name
  integer      att_val
  !
  include 'put_att_int.inc'
  !
  return
end subroutine write_att_int
!
subroutine read_att_text(ncid,var_name,att_name,att_val)
  !
  use definitions_netcdf_jfl
  use counters_netcdf_jfl
  !
  implicit none
  !
  include 'netcdf.inc'
  !
  ! arguments
  !
  integer ncid
  character(*) var_name
  character(*) att_name
  character(*) att_val
  !
  include 'get_att_text.inc'
  !
  return
end subroutine read_att_text
subroutine read_att_dble(ncid,var_name,att_name,att_val)
  !
  use definitions_netcdf_jfl
  use counters_netcdf_jfl
  !
  implicit none
  !
  include 'netcdf.inc'
  !
  ! arguments
  !
  integer ncid
  character(*) var_name
  character(*) att_name
  real*8       att_val
  !
  include 'get_att_dble.inc'
  !
  return
end subroutine read_att_dble
subroutine read_att_real(ncid,var_name,att_name,att_val)
  !
  use definitions_netcdf_jfl
  use counters_netcdf_jfl
  !
  implicit none
  !
  include 'netcdf.inc'
  !
  ! arguments
  !
  integer ncid
  character(*) var_name
  character(*) att_name
  real         att_val
  !
  include 'get_att_real.inc'
  !
  return
end subroutine read_att_real
subroutine read_att_int(ncid,var_name,att_name,att_val)
  !
  use definitions_netcdf_jfl
  use counters_netcdf_jfl
  !
  implicit none
  !
  include 'netcdf.inc'
  !
  ! arguments
  !
  integer ncid
  character(*) var_name
  character(*) att_name
  integer      att_val
  !
  include 'get_att_int.inc'
  !
  return
end subroutine read_att_int
