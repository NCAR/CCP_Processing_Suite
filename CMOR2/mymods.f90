!
! define module with maximum dims
!
module max_parms
  !
  integer,parameter::max_entries   =   5000
  integer,parameter::max_exprmnt   =   1000
  integer,parameter::max_cesm_vars =     20
  integer,parameter::max_ncfiles   =   1000
  integer,parameter::max_nchunks   =   1000
!  integer,parameter::max_ntimes    = 400000
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
  integer,dimension(20)::table_ids
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
     character(len=256)::runbegend,mipbegend,grid,compset,repotag,start_fin,mach,dout,forcing
     integer::runbeg,runend,mipbeg,mipend,runlen,miplen
     character(len=256),dimension(10)::icase
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
     character(len=256)::references,model_id,forcing,institute_id
     character(len=256)::parent_experiment_id,parent_experiment_rip,positive
     character(len=256)::forcing_note
     character(len=512)::ack_NC,ack_OR,ack_NE
!     character(len=1024)::comment
     character(len=256)::comment
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
     character(len=256)::table,entry,realm,sname,dims,comment
     character(len=256),dimension(100)::cesm_vars
     integer::ncesm_vars
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
  character(len=256),dimension(max_ncfiles,max_cesm_vars)::ncfile_nh,ncfile_sh ! For OImon processing
  integer,dimension(max_cesm_vars)::nc_nfiles
  integer,dimension(max_cesm_vars)::nc_nfiles_nh,nc_nfiles_sh ! For OImon processing
  integer,dimension(max_ncfiles,max_cesm_vars)::myncid,var_found,ntimes
  integer,dimension(max_ncfiles,max_cesm_vars)::myncid_nh,myncid_sh ! For OImon processing
!  real   ,dimension(max_ncfiles,max_cesm_vars,max_ntimes)::timevals
  integer::ifile
  logical::all_continue
  !
end module files_info
!
! Grid information
!
module grid_info
  real,dimension(:),    allocatable::atm_lats,atm_lons,atm_levs,atm_levs_bnds,atm_ilevs,atm_ilevs_bnds
  real,dimension(:),    allocatable::cosp_tau,cosp_prs,cosp_ht,cosp_dbze,cosp_sza,gaussian_wts
  real,dimension(:),    allocatable::atm_plev31,atm_plev23,atm_plev17,atm_plev8,atm_plev7,atm_plev7_bnds,atm_plev3
  real,dimension(:,:),  allocatable::atm_lats_bnds,atm_lons_bnds,landfrac,phis,area
  real,dimension(:,:),  allocatable::cosp_tau_bnds,cosp_prs_bnds,cosp_ht_bnds,cosp_dbze_bnds
  real,dimension(:),    allocatable::a_coeff,b_coeff,a_coeff_bnds,b_coeff_bnds
  real,dimension(:),    allocatable::ah_coeff,bh_coeff
  real,dimension(:,:),  allocatable::ice_t_lats,ice_t_lons,ice_u_lats,ice_u_lons
  real,dimension(:,:,:),allocatable::ice_t_lats_bnds,ice_t_lons_bnds,ice_u_lats_bnds,ice_u_lons_bnds
  real,dimension(:),    allocatable::ocn_t_levs,ocn_trans_lats,ocn_trans_levs,ocn_t_dz
  real,dimension(:,:),  allocatable::ocn_t_levs_bnds,ocn_trans_lats_bnds,ocn_trans_levs_bnds
  real,dimension(:,:),  allocatable::ocn_t_lats,ocn_t_lons,ocn_t_area,ocn_t_hte,ocn_t_htn
  real,dimension(:,:,:),allocatable::ocn_t_lats_bnds,ocn_t_lons_bnds
  real,dimension(:,:),  allocatable::ocn_u_lats,ocn_u_lons
  real,dimension(:,:,:),allocatable::ocn_u_lats_bnds,ocn_u_lons_bnds
  integer,dimension(:,:),allocatable::kmt
  real,dimension(:),    allocatable::lnd_lats,lnd_lons,lnd_levs,lnd_levs_bnds,atm_sites
  real,dimension(:,:,:),allocatable::lnd_zsoi,lnd_dzsoi  ! CLM soil depth (m), CLM soil layer thickness (m)
  real,dimension(:,:),  allocatable::lnd_lats_bnds,lnd_lons_bnds
  real::p0Pa,p0mb
  integer::nlons,nlats,nlevs,nilevs,nplev23,nplev17,nplev8,nplev7,nplev3,nplev31
  integer::ncosp_tau,ncosp_prs,ncosp_ht,ncosp_dbze,nsites,ncosp_sza
  integer::naxes,zfactor_id
  integer::nlats_trans,nmoc_z,ntrans_reg,nmoc_comp,ntrans_comp
  integer,dimension(10)::grid_id
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
