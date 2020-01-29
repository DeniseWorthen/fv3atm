module module_cplfields

  !-----------------------------------------------------------------------------
  ! This module contains the fv3 Coupling Fields: export and import
  !
  !-----------------------------------------------------------------------------

  use ESMF
  use NUOPC

  implicit none

  private

  type :: CplFields
   character(len=64)      :: fieldname
   character(len=10)      :: fieldtype
   logical                :: fieldshare
   logical                :: fieldvalid
  end type CplFields

  integer, parameter :: maxExportFields =  80
  integer, parameter :: maxImportFields =  20

  type(CplFields), public :: CplExportFields(maxExportFields)
  type(CplFields), public :: CplImportFields(maxImportFields)

  integer, public :: NexportFields
  integer, public :: NimportFields

! Export Fields ----------------------------------------
  real(kind=8), allocatable, public :: exportData(:,:,:)
  type(ESMF_Field), target, public, allocatable :: exportFields(:)

! Import Fields ----------------------------------------
  type(ESMF_Field), target, public, allocatable :: importFields(:)

  ! Methods
  public fillExportFields
  public queryFieldList
  public cplFieldGet
  public cplfld_setup

!-----------------------------------------------------------------------------
  contains
!-----------------------------------------------------------------------------

  subroutine fillExportFields(data_a2oi, rc)
    ! Fill updated data into the export Fields.
    real(kind=8), target, intent(in)            :: data_a2oi(:,:,:)
    integer, intent(out), optional              :: rc

    integer                                     :: localrc
    integer                                     :: n,dimCount
    logical                                     :: isCreated
    type(ESMF_TypeKind_Flag)                    :: datatype
    real(kind=ESMF_KIND_R4), dimension(:,:), pointer   :: datar42d
    real(kind=ESMF_KIND_R8), dimension(:,:), pointer   :: datar82d
    
!
    if (present(rc)) rc=ESMF_SUCCESS

    do n=1, size(exportFields)
      isCreated = ESMF_FieldIsCreated(exportFields(n), rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
      if (isCreated) then
! set data 
        call ESMF_FieldGet(exportFields(n), dimCount=dimCount, typekind=datatype, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
        if ( datatype == ESMF_TYPEKIND_R8) then
           if ( dimCount == 2) then
             call ESMF_FieldGet(exportFields(n),farrayPtr=datar82d,localDE=0, rc=localrc)
             if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
             datar82d = data_a2oi(:,:,n)
           endif
        else if ( datatype == ESMF_TYPEKIND_R4) then
           if ( dimCount == 2) then
             call ESMF_FieldGet(exportFields(n),farrayPtr=datar82d,localDE=0, rc=localrc)
             if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
             datar42d = data_a2oi(:,:,n)
           endif
        endif
      endif
    enddo
  end subroutine fillExportFields
!
!------------------------------------------------------------------------------
!
  integer function queryFieldList(fieldlist, fieldname, abortflag, rc)
    ! returns integer index of first found fieldname in fieldlist
    ! by default, will abort if field not found, set abortflag to false
    ! to turn off the abort.
    ! return value of < 1 means the field was not found

    character(len=*),intent(in) :: fieldlist(:)
    character(len=*),intent(in) :: fieldname
    logical, optional           :: abortflag
    integer, optional           :: rc

    integer :: n
    logical :: labort

    labort = .true.
    if (present(abortflag)) then
      labort = abortflag
    endif

    queryFieldList = 0
    n = 1
    do while (queryFieldList < 1 .and. n <= size(fieldlist))
      if (trim(fieldlist(n)) == trim(fieldname)) then
        queryFieldList = n
      else
        n = n + 1
      endif
    enddo

    if (labort .and. queryFieldList < 1) then
      call ESMF_LogWrite('queryFieldList ABORT on fieldname '//trim(fieldname), &
                          ESMF_LOGMSG_INFO, line=__LINE__, file=__FILE__, rc=rc)
      CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    endif
  end function queryFieldList
!
!------------------------------------------------------------------------------
!
  subroutine cplStateGet(state, fieldList, fieldCount, rc)

    character(len=*), intent(in)            :: state
    type(ESMF_Field), pointer,     optional :: fieldList(:)
    integer,          intent(out), optional :: fieldCount
    integer,          intent(out), optional :: rc

    !--- begin
    if (present(rc)) rc = ESMF_SUCCESS

    select case (trim(state))
      case ('import','i')
        if (present(fieldList )) fieldList  => importFields
        if (present(fieldCount)) fieldCount =  size(importFields)
      case ('export','o')
        if (present(fieldList )) fieldList  => exportFields
        if (present(fieldCount)) fieldCount =  size(exportFields)
      case default
        call ESMF_LogSetError(ESMF_RC_ARG_OUTOFRANGE, &
          msg="state argument can only be import(i)/export(o).", &
          line=__LINE__, file=__FILE__, rcToReturn=rc)
        return
    end select

  end subroutine cplStateGet


  subroutine cplFieldGet(state, name, localDe, &
                         farrayPtr2d, farrayPtr3d, farrayPtr4d, rc)

    character(len=*),   intent(in)            :: state
    character(len=*),   intent(in)            :: name
    integer,            intent(in),  optional :: localDe
    real(ESMF_KIND_R8), pointer,     optional :: farrayPtr2d(:,:)
    real(ESMF_KIND_R8), pointer,     optional :: farrayPtr3d(:,:,:)
    real(ESMF_KIND_R8), pointer,     optional :: farrayPtr4d(:,:,:,:)
    integer,            intent(out), optional :: rc

    !--- local variables
    integer                    :: localrc
    integer                    :: de, item, fieldCount, rank
    logical                    :: isCreated
    type(ESMF_Field), pointer  :: fieldList(:)
    character(len=ESMF_MAXSTR) :: fieldName

    !--- begin
    if (present(rc)) rc = ESMF_SUCCESS

    if (present(farrayPtr2d)) nullify(farrayPtr2d)
    if (present(farrayPtr3d)) nullify(farrayPtr3d)
    if (present(farrayPtr4d)) nullify(farrayPtr4d)

    de = 0
    if (present(localDe)) de = localDe

    call cplStateGet(state, fieldList=fieldList, fieldCount=fieldCount, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return

    do item = 1, fieldCount
      isCreated = ESMF_FieldIsCreated(fieldList(item), rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
      if (isCreated) then
        call ESMF_FieldGet(fieldList(item), name=fieldName, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
        if (trim(fieldName) == trim(name)) then
          call ESMF_FieldGet(fieldList(item), rank=rank, rc=localrc)
          if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
          select case (rank)
            case (2)
              if (present(farrayPtr2d)) then
                call ESMF_FieldGet(fieldList(item), localDe=de, farrayPtr=farrayPtr2d, rc=localrc)
                if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
              end if
            case (3)
              if (present(farrayPtr3d)) then
                call ESMF_FieldGet(fieldList(item), localDe=de, farrayPtr=farrayPtr3d, rc=localrc)
                if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
              end if
            case (4)
              if (present(farrayPtr4d)) then
                call ESMF_FieldGet(fieldList(item), localDe=de, farrayPtr=farrayPtr4d, rc=localrc)
                if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__, rcToReturn=rc)) return
              end if
            case default
              call ESMF_LogSetError(ESMF_RC_NOT_IMPL, msg="field rank should be 2, 3, or 4.", &
                                    line=__LINE__, file=__FILE__, rcToReturn=rc)
              return
          end select
          exit
        end if
      end if
    end do

  end subroutine cplFieldGet
!
!------------------------------------------------------------------------------
!
  subroutine cplfld_setup

   ! Field types should be provided for proper handling
   ! according to the table below:
   !  g : soil levels (3D)
   !  i : interface (3D)
   !  l : model levels (3D)
   !  s : surface (2D)
   !  t : tracers (4D)

   ! Set CplExportFields%share to .true. if field is provided as memory reference
   ! to coupled components
   ! Set CplImportFields%share to .true. if field is provided as memory reference
   ! from coupled components

          integer :: ii, rc
   character(240) :: msgString

   NexportFields = 0
   NimportFields = 0

    ! set defaults
    ii = 0; rc = 0
    CplExportFields(:)%fieldname  = ' '
    CplExportFields(:)%fieldtype  = 's'
    CplExportFields(:)%fieldshare = .false.
    CplExportFields(:)%fieldvalid = .false.  !not used for ExportFields
    
    ii = ii + 1; call filltype(CplExportFields(ii),                    'inst_pres_interface', 'i',  .true.)
    ii = ii + 1; call filltype(CplExportFields(ii),                       'inst_pres_levels', 'l',  .true.)
    ii = ii + 1; call filltype(CplExportFields(ii),                    'inst_geop_interface', 'i',  .true.)
    ii = ii + 1; call filltype(CplExportFields(ii),                       'inst_geop_levels', 'l',  .true.)
    ii = ii + 1; call filltype(CplExportFields(ii),                       'inst_temp_levels', 'l',  .true.)
    ii = ii + 1; call filltype(CplExportFields(ii),                 'inst_zonal_wind_levels', 'l',  .true.)
    ii = ii + 1; call filltype(CplExportFields(ii),                 'inst_merid_wind_levels', 'l',  .true.)
    ii = ii + 1; call filltype(CplExportFields(ii),                      'inst_omega_levels', 'l',  .true.)
    ii = ii + 1; call filltype(CplExportFields(ii),                  'inst_tracer_mass_frac', 't',  .true.)
    ii = ii + 1; call filltype(CplExportFields(ii),                              'soil_type', 's',  .true.)
    ii = ii + 1; call filltype(CplExportFields(ii),                        'inst_pbl_height', 's',  .true.)
    ii = ii + 1; call filltype(CplExportFields(ii),                      'surface_cell_area', 's',  .true.)
    ii = ii + 1; call filltype(CplExportFields(ii),        'inst_convective_rainfall_amount', 's',  .true.)
    ii = ii + 1; call filltype(CplExportFields(ii),  'inst_exchange_coefficient_heat_levels', 'l',  .true.)
    ii = ii + 1; call filltype(CplExportFields(ii),   'inst_spec_humid_conv_tendency_levels', 'l',  .true.)
    ii = ii + 1; call filltype(CplExportFields(ii),                 'inst_friction_velocity', 's',  .true.)
    ii = ii + 1; call filltype(CplExportFields(ii),                   'inst_rainfall_amount', 's',  .true.)
    ii = ii + 1; call filltype(CplExportFields(ii),             'inst_soil_moisture_content', 'g',  .true.)
    ii = ii + 1; call filltype(CplExportFields(ii),                 'inst_up_sensi_heat_flx', 's',  .true.)
    ii = ii + 1; call filltype(CplExportFields(ii),                'inst_lwe_snow_thickness', 's',  .true.)
    ii = ii + 1; call filltype(CplExportFields(ii),                        'vegetation_type', 's',  .true.)
    ii = ii + 1; call filltype(CplExportFields(ii),              'inst_vegetation_area_frac', 's',  .true.)
    ii = ii + 1; call filltype(CplExportFields(ii),                 'inst_surface_roughness', 's',  .true.)
    ii = ii + 1; call filltype(CplExportFields(ii),                  'mean_zonal_moment_flx', 's', .false.)
    ii = ii + 1; call filltype(CplExportFields(ii),                  'mean_merid_moment_flx', 's', .false.)
    ii = ii + 1; call filltype(CplExportFields(ii),                    'mean_sensi_heat_flx', 's', .false.)
    ii = ii + 1; call filltype(CplExportFields(ii),                    'mean_laten_heat_flx', 's', .false.)
    ii = ii + 1; call filltype(CplExportFields(ii),                       'mean_down_lw_flx', 's', .false.)
    ii = ii + 1; call filltype(CplExportFields(ii),                       'mean_down_sw_flx', 's', .false.)
    ii = ii + 1; call filltype(CplExportFields(ii),                         'mean_prec_rate', 's', .false.)
    ii = ii + 1; call filltype(CplExportFields(ii),                  'inst_zonal_moment_flx', 's', .false.)
    ii = ii + 1; call filltype(CplExportFields(ii),                  'inst_merid_moment_flx', 's', .false.)
    ii = ii + 1; call filltype(CplExportFields(ii),                    'inst_sensi_heat_flx', 's', .false.)
    ii = ii + 1; call filltype(CplExportFields(ii),                    'inst_laten_heat_flx', 's', .false.)
    ii = ii + 1; call filltype(CplExportFields(ii),                       'inst_down_lw_flx', 's', .false.)
    ii = ii + 1; call filltype(CplExportFields(ii),                       'inst_down_sw_flx', 's',  .true.)
    ii = ii + 1; call filltype(CplExportFields(ii),                     'inst_temp_height2m', 's', .false.)
    ii = ii + 1; call filltype(CplExportFields(ii),               'inst_spec_humid_height2m', 's', .false.)
    ii = ii + 1; call filltype(CplExportFields(ii),              'inst_zonal_wind_height10m', 's', .false.)
    ii = ii + 1; call filltype(CplExportFields(ii),              'inst_merid_wind_height10m', 's', .false.)
    ii = ii + 1; call filltype(CplExportFields(ii),               'inst_temp_height_surface', 's',  .true.)
    ii = ii + 1; call filltype(CplExportFields(ii),               'inst_pres_height_surface', 's', .false.)
    ii = ii + 1; call filltype(CplExportFields(ii),                    'inst_surface_height', 's', .false.)
    ii = ii + 1; call filltype(CplExportFields(ii),                        'mean_net_lw_flx', 's', .false.)
    ii = ii + 1; call filltype(CplExportFields(ii),                        'mean_net_sw_flx', 's', .false.)
    ii = ii + 1; call filltype(CplExportFields(ii),                        'inst_net_lw_flx', 's', .false.)
    ii = ii + 1; call filltype(CplExportFields(ii),                        'inst_net_sw_flx', 's', .false.)
    ii = ii + 1; call filltype(CplExportFields(ii),                'mean_down_sw_ir_dir_flx', 's', .false.)
    ii = ii + 1; call filltype(CplExportFields(ii),                'mean_down_sw_ir_dif_flx', 's', .false.)
    ii = ii + 1; call filltype(CplExportFields(ii),               'mean_down_sw_vis_dir_flx', 's', .false.)
    ii = ii + 1; call filltype(CplExportFields(ii),               'mean_down_sw_vis_dif_flx', 's', .false.)
    ii = ii + 1; call filltype(CplExportFields(ii),                'inst_down_sw_ir_dir_flx', 's', .false.)
    ii = ii + 1; call filltype(CplExportFields(ii),                'inst_down_sw_ir_dif_flx', 's', .false.)
    ii = ii + 1; call filltype(CplExportFields(ii),               'inst_down_sw_vis_dir_flx', 's', .false.)
    ii = ii + 1; call filltype(CplExportFields(ii),               'inst_down_sw_vis_dif_flx', 's', .false.)
    ii = ii + 1; call filltype(CplExportFields(ii),                 'mean_net_sw_ir_dir_flx', 's', .false.)
    ii = ii + 1; call filltype(CplExportFields(ii),                 'mean_net_sw_ir_dif_flx', 's', .false.)
    ii = ii + 1; call filltype(CplExportFields(ii),                'mean_net_sw_vis_dir_flx', 's', .false.)
    ii = ii + 1; call filltype(CplExportFields(ii),                'mean_net_sw_vis_dif_flx', 's', .false.)
    ii = ii + 1; call filltype(CplExportFields(ii),                 'inst_net_sw_ir_dir_flx', 's', .false.)
    ii = ii + 1; call filltype(CplExportFields(ii),                 'inst_net_sw_ir_dif_flx', 's', .false.)
    ii = ii + 1; call filltype(CplExportFields(ii),                'inst_net_sw_vis_dir_flx', 's', .false.)
    ii = ii + 1; call filltype(CplExportFields(ii),                'inst_net_sw_vis_dif_flx', 's', .false.)
    ii = ii + 1; call filltype(CplExportFields(ii),                     'inst_land_sea_mask', 's',  .true.)
    ii = ii + 1; call filltype(CplExportFields(ii),                'inst_temp_height_lowest', 's', .false.)
    ii = ii + 1; call filltype(CplExportFields(ii),          'inst_spec_humid_height_lowest', 's', .false.)
    ii = ii + 1; call filltype(CplExportFields(ii),          'inst_zonal_wind_height_lowest', 's', .false.)
    ii = ii + 1; call filltype(CplExportFields(ii),          'inst_merid_wind_height_lowest', 's', .false.)
    ii = ii + 1; call filltype(CplExportFields(ii),                'inst_pres_height_lowest', 's', .false.)
    ii = ii + 1; call filltype(CplExportFields(ii),                     'inst_height_lowest', 's', .false.)
    ii = ii + 1; call filltype(CplExportFields(ii),                        'mean_fprec_rate', 's', .false.)
 
    NexportFields = ii
    if(NexportFields > maxExportFields)then
     msgString = 'NexportFields > MaxExportFields'
     call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO, rc=rc)
     if (ESMF_LogFoundError(rcToCheck=-1, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    endif
    allocate(exportFields(1:NexportFields))

!-----------------------------------------------------------------------------

    ! set defaults
    ii = 0; rc = 0
    CplImportFields(:)%fieldname  = ' '
    CplImportFields(:)%fieldtype  = 's'
    CplImportFields(:)%fieldshare = .false.
    CplImportFields(:)%fieldvalid = .false.  !set dynamically

    ii = ii + 1; call filltype(CplImportFields(ii),                  'inst_tracer_mass_frac', 't',  .true.)
    ii = ii + 1; call filltype(CplImportFields(ii),                              'land_mask', 's', .false.)
    ii = ii + 1; call filltype(CplImportFields(ii),            'sea_ice_surface_temperature', 's', .false.)
    ii = ii + 1; call filltype(CplImportFields(ii),                'sea_surface_temperature', 's', .false.)
    ii = ii + 1; call filltype(CplImportFields(ii),                           'ice_fraction', 's', .false.)
    ii = ii + 1; call filltype(CplImportFields(ii),                         'mean_up_lw_flx', 's', .false.)
    ii = ii + 1; call filltype(CplImportFields(ii),                    'mean_laten_heat_flx', 's', .false.)
    ii = ii + 1; call filltype(CplImportFields(ii),                    'mean_sensi_heat_flx', 's', .false.)
    ii = ii + 1; call filltype(CplImportFields(ii),                  'mean_zonal_moment_flx', 's', .false.)
    ii = ii + 1; call filltype(CplImportFields(ii),                  'mean_merid_moment_flx', 's', .false.)
    ii = ii + 1; call filltype(CplImportFields(ii),                        'mean_ice_volume', 's', .false.)
    ii = ii + 1; call filltype(CplImportFields(ii),                       'mean_snow_volume', 's', .false.)
    ii = ii + 1; call filltype(CplImportFields(ii),             'inst_tracer_up_surface_flx', 'u',  .true.)
    ii = ii + 1; call filltype(CplImportFields(ii),           'inst_tracer_down_surface_flx', 'd',  .true.)
    ii = ii + 1; call filltype(CplImportFields(ii),             'inst_tracer_clmn_mass_dens', 'c',  .true.)
    ii = ii + 1; call filltype(CplImportFields(ii),              'inst_tracer_anth_biom_flx', 'b',  .true.)

    NimportFields = ii
    if(NimportFields > MaxImportFields)then
     msgString = 'NimportFields > MaxImportFields'
     call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO, rc=rc)
     if (ESMF_LogFoundError(rcToCheck=-1, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=__FILE__)) return
    endif
    allocate(importFields(1:NimportFields))

  end subroutine cplfld_setup

  subroutine filltype(fldtyp,vname,vtype,vshare)

  character(len=*), intent(in) :: vname,vtype
  logical, intent(in) :: vshare

  type(CplFields), intent(out) :: fldtyp

  fldtyp%fieldname  = trim(vname)
  fldtyp%fieldtype  = trim(vtype)
  fldtyp%fieldshare = vshare

 end subroutine filltype

end module module_cplfields
