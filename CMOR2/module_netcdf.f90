!
! define module with counters and maximum dimensions
!
module counters_netcdf_jfl
  !
  integer :: dim_counter  = 0
  integer :: var_counter  = 0
  integer :: time_counter = 0
  integer :: file_counter = 0
  integer, parameter :: dim_max = 10000
  integer, parameter :: var_max = 10000
  integer, parameter :: fil_max = 10000
  !
end module counters_netcdf_jfl
!
! define module with structure definitions
!
module definitions_netcdf_jfl
  !
  use counters_netcdf_jfl
  !
  type DimInfo
     character(len=132)::name,file_name
     integer           ::id,length
  end type DimInfo
  !
  type VarInfo
     character(len=132)   ::name,format,type,file_name,units
     real                 ::missing_value,FillValue
     integer              ::id,nvdims
     integer,dimension(10)::vdims
  end type VarInfo
  !
  type(DimInfo), dimension(dim_max) :: dim_info
  type(VarInfo), dimension(var_max) :: var_info
  character(len=132),dimension(fil_max) :: nc_filename
  !
end module definitions_netcdf_jfl
!
! define module with explicit subroutine interface
! and definition of generic subroutines (read_var and
! write_var).
!
module interfaces_netcdf_jfl
  !
  interface
     subroutine def_dim(ncid,dim_name&
          ,dim_length)
       integer, intent(in) :: ncid,dim_length
       character(*), intent(in) :: dim_name
     end subroutine def_dim
  end interface
  !
  interface
     subroutine def_var(ncid,var_name&
          ,format&
          ,type&
          ,units&
          ,long_name)
       integer, intent(in) :: ncid
       character(*), intent(in)           :: var_name
       character(*), intent(in)           :: format
       character(*), intent(in)           :: type
       character(*), intent(in), optional :: units
       character(*), intent(in), optional :: long_name
     end subroutine def_var
  end interface
  !
  interface def_att
     subroutine write_att_text(ncid,var_name,att_name,att_val)
       integer     , intent(in) :: ncid
       character(*), intent(in) :: var_name
       character(*), intent(in) :: att_name
       character(*), intent(in) :: att_val
     end subroutine write_att_text
     subroutine write_att_real(ncid,var_name,att_name,att_val)
       integer     , intent(in) :: ncid
       character(*), intent(in) :: var_name
       character(*), intent(in) :: att_name
       real        , intent(in) :: att_val
     end subroutine write_att_real
     subroutine write_att_dble(ncid,var_name,att_name,att_val)
       integer     , intent(in) :: ncid
       character(*), intent(in) :: var_name
       character(*), intent(in) :: att_name
       real*8      , intent(in) :: att_val
     end subroutine write_att_dble
     subroutine write_att_int(ncid,var_name,att_name,att_val)
       integer     , intent(in) :: ncid
       character(*), intent(in) :: var_name
       character(*), intent(in) :: att_name
       integer     , intent(in) :: att_val
     end subroutine write_att_int
  end interface
  !
  interface get_att
     subroutine read_att_text(ncid,var_name,att_name,att_val)
       integer     , intent(in) :: ncid
       character(*), intent(in) :: var_name
       character(*), intent(in) :: att_name
       character(*), intent(in) :: att_val
     end subroutine read_att_text
     subroutine read_att_real(ncid,var_name,att_name,att_val)
       integer     , intent(in) :: ncid
       character(*), intent(in) :: var_name
       character(*), intent(in) :: att_name
       real        , intent(in) :: att_val
     end subroutine read_att_real
     subroutine read_att_dble(ncid,var_name,att_name,att_val)
       integer     , intent(in) :: ncid
       character(*), intent(in) :: var_name
       character(*), intent(in) :: att_name
       real*8      , intent(in) :: att_val
     end subroutine read_att_dble
     subroutine read_att_int(ncid,var_name,att_name,att_val)
       integer     , intent(in) :: ncid
       character(*), intent(in) :: var_name
       character(*), intent(in) :: att_name
       integer     , intent(in) :: att_val
     end subroutine read_att_int
  end interface
  !
  interface write_var
     subroutine write_var_c_0(ncid,var_name,array)
       integer     , intent(in) :: ncid
       character(*), intent(in) :: var_name
       character(*), intent(in) :: array
     end subroutine write_var_c_0
     subroutine write_var_c_1(ncid,var_name,array)
       integer     , intent(in) :: ncid
       character(*), intent(in) :: var_name
       character(*), intent(in) :: array(:)
     end subroutine write_var_c_1
     subroutine write_var_c_2(ncid,var_name,array)
       integer     , intent(in) :: ncid
       character(*), intent(in) :: var_name
       character(*), intent(in) :: array(:,:)
     end subroutine write_var_c_2
     subroutine write_var_i_0(ncid,var_name,array)
       integer     , intent(in) :: ncid
       character(*), intent(in) :: var_name
       integer     , intent(in) :: array
     end subroutine write_var_i_0
     subroutine write_var_i_1(ncid,var_name,array)
       integer     , intent(in) :: ncid
       character(*), intent(in) :: var_name
       integer     , intent(in) :: array(:)
     end subroutine write_var_i_1
     subroutine write_var_i_2(ncid,var_name,array)
       integer     , intent(in) :: ncid
       character(*), intent(in) :: var_name
       integer     , intent(in) :: array(:,:)
     end subroutine write_var_i_2
     subroutine write_var_i_3(ncid,var_name,array)
       integer     , intent(in) :: ncid
       character(*), intent(in) :: var_name
       integer     , intent(in) :: array(:,:,:)
     end subroutine write_var_i_3
     subroutine write_var_f_0(ncid,var_name,array)
       integer     , intent(in) :: ncid
       character(*), intent(in) :: var_name
       real        , intent(in) :: array
     end subroutine write_var_f_0
     subroutine write_var_f_1(ncid,var_name,array)
       integer     , intent(in) :: ncid
       character(*), intent(in) :: var_name
       real        , intent(in) :: array(:)
     end subroutine write_var_f_1
     subroutine write_var_f_2(ncid,var_name,array)
       integer     , intent(in) :: ncid
       character(*), intent(in) :: var_name
       real        , intent(in) :: array(:,:)
     end subroutine write_var_f_2
     subroutine write_var_f_3(ncid,var_name,array)
       integer     , intent(in) :: ncid
       character(*), intent(in) :: var_name
       real        , intent(in) :: array(:,:,:)
     end subroutine write_var_f_3
     subroutine write_var_f_4(ncid,var_name,array)
       integer     , intent(in) :: ncid
       character(*), intent(in) :: var_name
       real        , intent(in) :: array(:,:,:,:)
     end subroutine write_var_f_4
     subroutine write_var_d_0(ncid,var_name,array)
       integer     , intent(in) :: ncid
       character(*), intent(in) :: var_name
       real*8      , intent(in) :: array
     end subroutine write_var_d_0
     subroutine write_var_d_1(ncid,var_name,array)
       integer     , intent(in) :: ncid
       character(*), intent(in) :: var_name
       real*8      , intent(in) :: array(:)
     end subroutine write_var_d_1
     subroutine write_var_d_2(ncid,var_name,array)
       integer     , intent(in) :: ncid
       character(*), intent(in) :: var_name
       real*8      , intent(in) :: array(:,:)
     end subroutine write_var_d_2
     subroutine write_var_d_3(ncid,var_name,array)
       integer     , intent(in) :: ncid
       character(*), intent(in) :: var_name
       real*8      , intent(in) :: array(:,:,:)
     end subroutine write_var_d_3
  end interface
  !
  interface read_var
     subroutine read_var_c_0(ncid,var_name,array)
       integer     , intent(in) :: ncid
       character(*), intent(in) :: var_name
       character(*), intent(in) :: array
     end subroutine read_var_c_0
     subroutine read_var_c_1(ncid,var_name,array)
       integer     , intent(in) :: ncid
       character(*), intent(in) :: var_name
       character(*), intent(in) :: array(:)
     end subroutine read_var_c_1
     subroutine read_var_c_2(ncid,var_name,array)
       integer     , intent(in) :: ncid
       character(*), intent(in) :: var_name
       character(*), intent(in) :: array(:,:)
     end subroutine read_var_c_2
     subroutine read_var_i_0(ncid,var_name,array)
       integer     , intent(in) :: ncid
       character(*), intent(in) :: var_name
       integer     , intent(in) :: array
     end subroutine read_var_i_0
     subroutine read_var_i_1(ncid,var_name,array)
       integer     , intent(in) :: ncid
       character(*), intent(in) :: var_name
       integer     , intent(in) :: array(:)
     end subroutine read_var_i_1
     subroutine read_var_i8_1(ncid,var_name,array)
       integer     , intent(in) :: ncid
       character(*), intent(in) :: var_name
       integer*8   , intent(in) :: array(:)
     end subroutine read_var_i8_1
     subroutine read_var_i_2(ncid,var_name,array)
       integer     , intent(in) :: ncid
       character(*), intent(in) :: var_name
       integer     , intent(in) :: array(:,:)
     end subroutine read_var_i_2
     subroutine read_var_i_3(ncid,var_name,array)
       integer     , intent(in) :: ncid
       character(*), intent(in) :: var_name
       integer     , intent(in) :: array(:,:,:)
     end subroutine read_var_i_3
     subroutine read_var_f_0(ncid,var_name,array)
       integer     , intent(in) :: ncid
       character(*), intent(in) :: var_name
       real        , intent(in) :: array
     end subroutine read_var_f_0
     subroutine read_var_f_1(ncid,var_name,array)
       integer     , intent(in) :: ncid
       character(*), intent(in) :: var_name
       real        , intent(in) :: array(:)
     end subroutine read_var_f_1
     subroutine read_var_f_2(ncid,var_name,array)
       integer     , intent(in) :: ncid
       character(*), intent(in) :: var_name
       real        , intent(in) :: array(:,:)
     end subroutine read_var_f_2
     subroutine read_var_f_3(ncid,var_name,array)
       integer     , intent(in) :: ncid
       character(*), intent(in) :: var_name
       real        , intent(in) :: array(:,:,:)
     end subroutine read_var_f_3
     subroutine read_var_f_4(ncid,var_name,array)
       integer     , intent(in) :: ncid
       character(*), intent(in) :: var_name
       real        , intent(in) :: array(:,:,:,:)
     end subroutine read_var_f_4
     subroutine read_var_d_0(ncid,var_name,array)
       integer     , intent(in) :: ncid
       character(*), intent(in) :: var_name
       double precision, intent(in) :: array
     end subroutine read_var_d_0
     subroutine read_var_d_1(ncid,var_name,array)
       integer     , intent(in) :: ncid
       character(*), intent(in) :: var_name
       double precision, intent(in) :: array(:)
     end subroutine read_var_d_1
     subroutine read_var_d_2(ncid,var_name,array)
       integer     , intent(in) :: ncid
       character(*), intent(in) :: var_name
       double precision, intent(in) :: array(:,:)
     end subroutine read_var_d_2
  end interface
  !
end module interfaces_netcdf_jfl
