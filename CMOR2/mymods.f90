!
! define module with maximum dims
!
module max_parms
  !
  integer,parameter::max_entries   = 10000
  integer,parameter::max_exprmnt   = 10000
  integer,parameter::max_cesm_vars =    10
  integer,parameter::max_ncfiles   =   200
  integer,parameter::max_nchunks   =  1000
  !
end module max_parms
!
! CMOR "table" information
!
module table_info
  use max_parms
  type TableInfo
     character(len=256)::approx_interval,axis,axis_entry,baseURL,bounds_values,cell_measures,cell_methods,cf_version,climatology,cmor_version
     character(len=256)::comment,coords_attrib,dimensions,expt_id_ok,forcings,formula,frequency,generic_levels,index_only,long_name,missing_value
     character(len=256)::modeling_realm,must_call_cmor_grid,must_have_bounds,ok_max_mean_abs,ok_min_mean_abs,out_name,positive,product,project_id
     character(len=256)::required_global_attributes,standard_name,stored_direction,table_date,table_id,tolerance,type,units,valid_max,valid_min,value
     character(len=256)::requested,requested_bounds,variable_entry,z_bounds_factors,z_factors
  end type TableInfo
  integer::num_tab,tab_found
  !
  ! table_ids(1) = Amon, day, OImon, Omon, etc.' table_ids(2) = CMIP5_grids. Always.
  !
  integer,dimension(10)::table_ids
  type(TableInfo),dimension(0:max_entries)::table
end module table_info
!
! CESM experiment information
!
module exp_info
  !
  use max_parms
  !
  type SimInfo
     character(len=256)::case,loc,model_id,expt_id,rip_code,cmip,run_refcase,run_refdate
     character(len=256)::begin_end,grid,compset,repotag,start_fin,mach,dout,forcing
     integer::begyr,endyr,length
  end type SimInfo
  integer::num_exp,exp_found,parent_found
  character(len=256)::case_read,comp_read
  !
  type(SimInfo),dimension(max_entries)::exp
  !
end module exp_info
!
! CMOR arguments information
!
module mycmor_info
  !
  use max_parms
  !
  type CMORInfo
     character(len=256)::table_file,outpath,experiment_id,institution,source,calendar,contact,history
     character(len=256)::comment,references,model_id,forcing,institute_id
     character(len=256)::parent_experiment_id,parent_experiment_rip,positive
     character(len=256)::forcing_note
     character(len=512)::ack_NC,ack_OR,ack_NE
     integer::realization,initialization_method,physics_version
     double precision::branch_time
  end type CMORInfo
  !
  type(CMORInfo)::mycmor
  !
end module mycmor_info
!
! Crosswalk (xwalk) information
!
module xwalk_info
  !
  use max_parms
  !
  type XWInfo
     character(len=256)::table,entry,realm,comment
     character(len=256),dimension(10)::cesm_vars
     integer::ncesm_vars
     real::priority
  end type XWInfo
  integer::num_xw,xw_found
  !
  type(XWInfo),dimension(max_entries)::xw
  !
end module xwalk_info
!
! Input netCDF files information
!
module files_info
  use max_parms
  !
  character(len=256),dimension(max_ncfiles,max_cesm_vars)::ncfile
  integer,dimension(max_cesm_vars)::nc_nfiles
  integer,dimension(max_ncfiles,max_cesm_vars)::ncid,var_found,ntimes
  logical,dimension(max_ncfiles,max_cesm_vars)::exists
  !
end module files_info
!
! Grid information
!
module grid_info
  real,dimension(:),    allocatable::atm_lats,atm_lons,atm_levs,atm_levs_bnds,atm_plevs,atm_ilevs,atm_ilevs_bnds
  real,dimension(:,:),  allocatable::atm_lats_bnds,atm_lons_bnds
  real,dimension(:),    allocatable::a_coeff,b_coeff,a_coeff_bnds,b_coeff_bnds
  real,dimension(:,:),  allocatable::ice_lats,ice_lons
  real,dimension(:,:,:),allocatable::ice_lats_bnds,ice_lons_bnds
  real,dimension(:),    allocatable::ocn_levs,ocn_trans_lats,ocn_trans_levs
  real,dimension(:,:),  allocatable::ocn_levs_bnds,ocn_trans_lats_bnds,ocn_trans_levs_bnds
  real,dimension(:,:,:),allocatable::ocn_lats_bnds,ocn_lons_bnds
  real,dimension(:,:),  allocatable::ocn_lats,ocn_lons
  integer,dimension(:,:),allocatable::kmt
  real,dimension(:),    allocatable::lnd_lats,lnd_lons,lnd_levs,lnd_levs_bnds
  real,dimension(:,:,:),allocatable::lnd_zsoi,lnd_dzsoi  ! CLM soil depth (m), CLM soil layer thickness (m)
  real,dimension(:,:),  allocatable::lnd_lats_bnds,lnd_lons_bnds
  real::p0
  integer::nlons,nlats,nlevs,nilevs,nplevs,naxes,zfactor_id
  integer::nlats_trans,nmoc_z,ntrans_reg,nmoc_comp,ntrans_comp
  integer,dimension(1)::grid_id
  integer,dimension(10)::axis_ids
  character(len=256),dimension(10)::dimnames,dimunits
  character(len=256)::time_units
end module grid_info
!
! Output time chunks information
!
module output_times_info
  use max_parms
  integer,dimension(max_nchunks)::nchunks,tidx1,tidx2
end module output_times_info
