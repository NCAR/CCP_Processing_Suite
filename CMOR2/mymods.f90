!
! define module with maximum dims
!
module max_parms
  !
  integer,parameter::max_entries = 10000
  integer,parameter::max_exprmnt = 10000
  !
end module max_parms
!
! define module with table structure definitions
!
module table_info
  use max_parms
  type TableInfo
     character(len=512)::approx_interval,axis,axis_entry,baseURL,bounds_values,cell_measures,cell_methods,cf_version,climatology,cmor_version
     character(len=512)::comment,coords_attrib,dimensions,expt_id_ok,forcings,formula,frequency,generic_levels,index_only,long_name,missing_value
     character(len=512)::modeling_realm,must_call_cmor_grid,must_have_bounds,ok_max_mean_abs,ok_min_mean_abs,out_name,positive,product,project_id
     character(len=512)::required_global_attributes,standard_name,stored_direction,table_date,table_id,tolerance,type,units,valid_max,valid_min,value
     character(len=512)::requested,requested_bounds,variable_entry,z_bounds_factors,z_factors
  end type TableInfo
  integer::num_tab
  type(TableInfo),dimension(0:max_entries)::table
end module table_info
!
! define module with experiment structure definitions
!
module exp_info
  !
  use max_parms
  !
  type SimInfo
     character(len=512)::case,loc,model_id,expt_id,rip_code,cmip,run_refcase,run_refdate
     character(len=512)::begin_end,grid,compset,repotag,start_fin,mach,dout,forcing
  end type SimInfo
  integer::num_exp
  !
  type(SimInfo),dimension(max_entries)::exp
  !
end module exp_info
!
! define module with xwalk structure definitions
!
module xwalk_info
  !
  use max_parms
  !
  type XWInfo
     character(len=512)::varin2d,units2d,entry2d
  end type XWInfo
  integer::num_xw
  !
  type(XWInfo),dimension(max_entries)::xw
  !
end module xwalk_info
!
!
!
module atm_grid_info
  double precision,dimension(:), allocatable::alats,alons,slon,slat,plevs,zlevs,zlev_bnds
  double precision,dimension(:), allocatable::a_coeff,b_coeff,a_coeff_bnds,b_coeff_bnds
  double precision,dimension(:,:),allocatable::bnds_lat,bnds_lon
  double precision::p0
  integer::nlons,nlats,nlevs
end module
