!     Subroutines for basic 3D linear elastic elements 



!==========================SUBROUTINE el_linelast_3dbasic ==============================
subroutine el_linelast_3dbasic(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
    n_properties, element_properties, element_coords, length_coord_array, &                      ! Input variables
    dof_increment, dof_total, length_dof_array, &                                                ! Input variables
    n_state_variables, initial_state_variables, &                                                ! Input variables
    updated_state_variables,element_stiffness,element_residual, fail)      ! Output variables                          ! Output variables
    use Types
    use ParamIO
    !  use Globals, only: TIME,DTIME  For a time dependent problem uncomment this line to access the time increment and total time
    use Mesh, only : node
    use Element_Utilities, only : N => shape_functions_3D
    use Element_Utilities, only : dNdxi => shape_function_derivatives_3D
    use Element_Utilities, only:  dNdx => shape_function_spatial_derivatives_3D
    use Element_Utilities, only:  dNbardx => vol_avg_shape_function_derivatives_3D
    use Element_Utilities, only : xi => integrationpoints_3D, w => integrationweights_3D
    use Element_Utilities, only : dxdxi => jacobian_3D
    use Element_Utilities, only : initialize_integration_points
    use Element_Utilities, only : calculate_shapefunctions
    use Element_Utilities, only : invert_small
    implicit none

    integer, intent( in )         :: lmn                                                    ! Element number
    integer, intent( in )         :: element_identifier                                     ! Flag identifying element type (specified in .in file)
    integer, intent( in )         :: n_nodes                                                ! # nodes on the element
    integer, intent( in )         :: n_properties                                           ! # properties for the element
    integer, intent( in )         :: length_coord_array                                     ! Total # coords
    integer, intent( in )         :: length_dof_array                                       ! Total # DOF
    integer, intent( in )         :: n_state_variables                                      ! # state variables for the element

    type (node), intent( in )     :: node_property_list(n_nodes)                  ! Data structure describing storage for nodal variables - see below
    !  type node
    !      sequence
    !      integer :: flag                          ! Integer identifier
    !      integer :: coord_index                   ! Index of first coordinate in coordinate array
    !      integer :: n_coords                      ! Total no. coordinates for the node
    !      integer :: dof_index                     ! Index of first DOF in dof array
    !      integer :: n_dof                         ! Total no. of DOF for node
    !   end type node
    !   Access these using node_property_list(k)%n_coords eg to find the number of coords for the kth node on the element

    real( prec ), intent( in )    :: element_coords(length_coord_array)                     ! Coordinates, stored as x1,(x2),(x3) for each node in turn
    real( prec ), intent( in )    :: dof_increment(length_dof_array)                        ! DOF increment, stored as du1,du2,du3,du4... for each node in turn
    real( prec ), intent( in )    :: dof_total(length_dof_array)                            ! accumulated DOF, same storage as for increment

    real( prec ), intent( in )    :: element_properties(n_properties)                       ! Element or material properties, stored in order listed in input file
    real( prec ), intent( in )    :: initial_state_variables(n_state_variables)             ! Element state variables.  Defined in this routine
  
    logical, intent( out )        :: fail                                                   ! Set to .true. to force a timestep cutback
    real( prec ), intent( inout ) :: updated_state_variables(n_state_variables)             ! State variables at end of time step
    real( prec ), intent( out )   :: element_stiffness(length_dof_array,length_dof_array)   ! Element stiffness (ROW,COLUMN)
    real( prec ), intent( out )   :: element_residual(length_dof_array)                     ! Element residual force (ROW)
          

    ! Local Variables
    integer      :: n_points,kint

    real (prec)  ::  strain(6), dstrain(6)           ! Strain vector contains [ep11, ep22, ep33, 2ep12, 2ep13, 2ep23]
    real (prec)  ::  stress(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  D(6,6)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)

    real (prec)  ::  B(6,length_dof_array)             ! strain = B*(dof_total+dof_increment)
    real (prec)  :: Bbar(6,length_dof_array)
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  :: elvol
    !     Subroutine to compute element stiffness matrix and residual force vector for 3D linear elastic elements
    !     El props are:

    !     element_properties(1)         sigma zero
    !     element_properties(2)         epsilon zero
    !     element_properties(3)         n
    !     element_properties(4)         K
    fail = .false.
    
    x = reshape(element_coords,(/3,length_coord_array/3/))
    !hypo_n = 10.d0
    !K_in = 1000.d0
    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

    call initialize_integration_points(n_points, n_nodes, xi, w)

    dNbardx = 0.d0
    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
        dNbardx(1:n_nodes,1:3) = dNdx*w(kint)*determinant
        elvol = w(kint)*determinant
    end do
    ! calculate the element volume
    dNbardx = dNbardx/elvol

    element_residual = 0.d0
    element_stiffness = 0.d0


    !     --  Loop over integration points
    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
        !write(6,*) 'dNdx = ', dNdx
        B = 0.d0
        B(1,1:3*n_nodes-2:3) = dNdx(1:n_nodes,1)
        B(2,2:3*n_nodes-1:3) = dNdx(1:n_nodes,2)
        B(3,3:3*n_nodes:3)   = dNdx(1:n_nodes,3)
        B(4,1:3*n_nodes-2:3) = dNdx(1:n_nodes,2)
        B(4,2:3*n_nodes-1:3) = dNdx(1:n_nodes,1)
        B(5,1:3*n_nodes-2:3) = dNdx(1:n_nodes,3)
        B(5,3:3*n_nodes:3)   = dNdx(1:n_nodes,1)
        B(6,2:3*n_nodes-1:3) = dNdx(1:n_nodes,3)
        B(6,3:3*n_nodes:3)   = dNdx(1:n_nodes,2)
        !write(6,*) 'B = ', B
        Bbar = 0.d0
        Bbar(1,1:3*n_nodes-2:3) = B(1,1:3*n_nodes-2:3) + 0.33D0*(dNbardx(1:n_nodes,1) -dNdx(1:n_nodes,1))
        Bbar(1,2:3*n_nodes-1:3) = B(1,2:3*n_nodes-1:3) + 0.33D0*(dNbardx(1:n_nodes,2)-dNdx(1:n_nodes,2))
        Bbar(1,3:3*n_nodes:3) = B(1,3:3*n_nodes:3) + 0.33D0*(dNbardx(1:n_nodes,3)-dNdx(1:n_nodes,3))
        Bbar(2,1:3*n_nodes-2:3) = B(2,1:3*n_nodes-2:3) + 0.33D0*(dNbardx(1:n_nodes,1) -dNdx(1:n_nodes,1))
        Bbar(2,2:3*n_nodes-1:3) = B(2,2:3*n_nodes-1:3) + 0.33D0*(dNbardx(1:n_nodes,2)-dNdx(1:n_nodes,2))
        Bbar(2,3:3*n_nodes:3) = B(2,3:3*n_nodes:3) + 0.33D0*(dNbardx(1:n_nodes,3)-dNdx(1:n_nodes,3))
        Bbar(3,1:3*n_nodes-2:3) = B(3,1:3*n_nodes-2:3) + 0.33D0*(dNbardx(1:n_nodes,1) -dNdx(1:n_nodes,1))
        Bbar(3,2:3*n_nodes-1:3) = B(3,2:3*n_nodes-1:3) + 0.33D0*(dNbardx(1:n_nodes,2)-dNdx(1:n_nodes,2))
        Bbar(3,3:3*n_nodes:3) = B(3,3:3*n_nodes:3) + 0.33D0*(dNbardx(1:n_nodes,3)-dNdx(1:n_nodes,3))
        Bbar(4,1:3*n_nodes-2:3) = B(4,1:3*n_nodes-2:3)
        Bbar(4,2:3*n_nodes-1:3) = B(4,2:3*n_nodes-1:3)
        Bbar(5,1:3*n_nodes-2:3) = B(5,1:3*n_nodes-2:3)
        Bbar(5,3:3*n_nodes:3)   = B(5,3:3*n_nodes:3)
        Bbar(6,2:3*n_nodes-1:3) = B(6,2:3*n_nodes-1:3)
        Bbar(6,3:3*n_nodes:3)   = B(6,3:3*n_nodes:3)
        !write(6,*) 'bbar = ', Bbar
        !write(6,*) 'dofs = ', dof_total
        !write(6,*) 'dof is = ', dof_increment
        strain = matmul(B,dof_total)
        dstrain = matmul(B,dof_increment)

        !stress = matmul(D,strain+dstrain)
        call hypoelasticmaterial(strain+dstrain,element_properties,n_properties,stress,D)

        element_residual(1:3*n_nodes) = element_residual(1:3*n_nodes) - matmul(transpose(B),stress)*w(kint)*determinant

        element_stiffness(1:3*n_nodes,1:3*n_nodes) = element_stiffness(1:3*n_nodes,1:3*n_nodes) &
            + matmul(transpose(B(1:6,1:3*n_nodes)),matmul(D,B(1:6,1:3*n_nodes)))*w(kint)*determinant

    end do

    return
end subroutine el_linelast_3dbasic


!==========================SUBROUTINE el_linelast_3dbasic_dynamic ==============================
subroutine el_linelast_3dbasic_dynamic(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
    n_properties, element_properties,element_coords, length_coord_array, &                               ! Input variables
    dof_increment, dof_total, length_dof_array,  &                                                       ! Input variables
    n_state_variables, initial_state_variables, &                                                        ! Input variables
    updated_state_variables,element_residual,element_deleted)                                            ! Output variables
    use Types
    use ParamIO
    use Mesh, only : node
    use Element_Utilities, only : N => shape_functions_3D
    use Element_Utilities, only:  dNdxi => shape_function_derivatives_3D
    use Element_Utilities, only:  dNdx => shape_function_spatial_derivatives_3D
    use Element_Utilities, only : xi => integrationpoints_3D, w => integrationweights_3D
    use Element_Utilities, only : dxdxi => jacobian_3D
    use Element_Utilities, only : initialize_integration_points
    use Element_Utilities, only : calculate_shapefunctions
    use Element_Utilities, only : invert_small
    implicit none

    integer, intent( in )         :: lmn                                                    ! Element number
    integer, intent( in )         :: element_identifier                                     ! Flag identifying element type (specified in .in file)
    integer, intent( in )         :: n_nodes                                                ! # nodes on the element
    integer, intent( in )         :: n_properties                                           ! # properties for the element
    integer, intent( in )         :: length_coord_array                                     ! Total # coords
    integer, intent( in )         :: length_dof_array                                       ! Total # DOF
    integer, intent( in )         :: n_state_variables                                      ! # state variables for the element

    type (node), intent( in )     :: node_property_list(n_nodes)                  ! Data structure describing storage for nodal variables - see below
    !  type node
    !      sequence
    !      integer :: flag                          ! Integer identifier
    !      integer :: coord_index                   ! Index of first coordinate in coordinate array
    !      integer :: n_coords                      ! Total no. coordinates for the node
    !      integer :: dof_index                     ! Index of first DOF in dof array
    !      integer :: n_dof                         ! Total no. of DOF for node
    !   end type node
    !   Access these using node_property_list(k)%n_coords eg to find the number of coords for the kth node on the element

    real( prec ), intent( in )    :: element_coords(length_coord_array)                     ! Coordinates, stored as x1,(x2),(x3) for each node in turn
    real( prec ), intent( in )    :: dof_increment(length_dof_array)                        ! DOF increment, stored as du1,du2,du3,du4... for each node in turn
    real( prec ), intent( in )    :: dof_total(length_dof_array)                            ! accumulated DOF, same storage as for increment

    real( prec ), intent( in )    :: element_properties(n_properties)                       ! Element or material properties, stored in order listed in input file
    real( prec ), intent( in )    :: initial_state_variables(n_state_variables)             ! Element state variables.  Defined in this routine
               
    real( prec ), intent( inout ) :: updated_state_variables(n_state_variables)             ! State variables at end of time step
    real( prec ), intent( out )   :: element_residual(length_dof_array)                     ! Element residual force (ROW)
          
    logical, intent( inout )      :: element_deleted                                        ! Set to .true. to delete element

    ! Local Variables
    integer      :: n_points,kint

    real (prec)  ::  strain(6), dstrain(6)             ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  stress(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  D(6,6)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(6,length_dof_array)             ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  :: E, xnu, D44, D11, D12              ! Material properties
    !
    !     Subroutine to compute element force vector for a linear elastodynamic problem
    !     El props are:

    !     element_properties(1)         Young's modulus
    !     element_properties(2)         Poisson's ratio
    
    x = reshape(element_coords,(/3,length_coord_array/3/))

    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

    call initialize_integration_points(n_points, n_nodes, xi, w)

    element_residual = 0.d0

    D = 0.d0
    E = element_properties(1)
    xnu = element_properties(2)
    d44 = 0.5D0*E/(1+xnu) 
    d11 = (1.D0-xnu)*E/( (1+xnu)*(1-2.D0*xnu) )
    d12 = xnu*E/( (1+xnu)*(1-2.D0*xnu) )
    D(1:3,1:3) = d12
    D(1,1) = d11
    D(2,2) = d11
    D(3,3) = d11
    D(4,4) = d44
    D(5,5) = d44
    D(6,6) = d44
  
    !     --  Loop over integration points
    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
        B = 0.d0
        B(1,1:3*n_nodes-2:3) = dNdx(1:n_nodes,1)
        B(2,2:3*n_nodes-1:3) = dNdx(1:n_nodes,2)
        B(3,3:3*n_nodes:3)   = dNdx(1:n_nodes,3)
        B(4,1:3*n_nodes-2:3) = dNdx(1:n_nodes,2)
        B(4,2:3*n_nodes-1:3) = dNdx(1:n_nodes,1)
        B(5,1:3*n_nodes-2:3) = dNdx(1:n_nodes,3)
        B(5,3:3*n_nodes:3)   = dNdx(1:n_nodes,1)
        B(6,2:3*n_nodes-1:3) = dNdx(1:n_nodes,3)
        B(6,3:3*n_nodes:3)   = dNdx(1:n_nodes,2)

        strain = matmul(B,dof_total)
        dstrain = matmul(B,dof_increment)

        stress = matmul(D,strain+dstrain)
        element_residual(1:3*n_nodes) = element_residual(1:3*n_nodes) - matmul(transpose(B),stress)*w(kint)*determinant

    end do
  
    return
end subroutine el_linelast_3dbasic_dynamic


!==========================SUBROUTINE fieldvars_linelast_3dbasic ==============================
subroutine fieldvars_linelast_3dbasic(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
    n_properties, element_properties,element_coords,length_coord_array, &                                ! Input variables
    dof_increment, dof_total, length_dof_array,  &                                                      ! Input variables
    n_state_variables, initial_state_variables,updated_state_variables, &                               ! Input variables
    n_field_variables,field_variable_names, &                                                           ! Field variable definition
    nodal_fieldvariables)      ! Output variables
    use Types
    use ParamIO
    use Mesh, only : node
    use Element_Utilities, only : N => shape_functions_3D
    use Element_Utilities, only: dNdxi => shape_function_derivatives_3D
    use Element_Utilities, only: dNdx => shape_function_spatial_derivatives_3D
    use Element_Utilities, only:  dNbardx => vol_avg_shape_function_derivatives_3D
    use Element_Utilities, only : xi => integrationpoints_3D, w => integrationweights_3D
    use Element_Utilities, only : dxdxi => jacobian_3D
    use Element_Utilities, only : initialize_integration_points
    use Element_Utilities, only : calculate_shapefunctions
    use Element_Utilities, only : invert_small
    implicit none

    integer, intent( in )         :: lmn                                                    ! Element number
    integer, intent( in )         :: element_identifier                                     ! Flag identifying element type (specified in .in file)
    integer, intent( in )         :: n_nodes                                                ! # nodes on the element
    integer, intent( in )         :: n_properties                                           ! # properties for the element
    integer, intent( in )         :: length_coord_array                                     ! Total # coords
    integer, intent( in )         :: length_dof_array                                       ! Total # DOF
    integer, intent( in )         :: n_state_variables                                      ! # state variables for the element
    integer, intent( in )         :: n_field_variables                                      ! # field variables

    type (node), intent( in )     :: node_property_list(n_nodes)                  ! Data structure describing storage for nodal variables - see below
    !  type node
    !      sequence
    !      integer :: flag                          ! Integer identifier
    !      integer :: coord_index                   ! Index of first coordinate in coordinate array
    !      integer :: n_coords                      ! Total no. coordinates for the node
    !      integer :: dof_index                     ! Index of first DOF in dof array
    !      integer :: n_dof                         ! Total no. of DOF for node
    !   end type node
    !   Access these using node_property_list(k)%n_coords eg to find the number of coords for the kth node on the element

    character (len=100), intent(in) :: field_variable_names(n_field_variables)

    real( prec ), intent( in )    :: element_coords(length_coord_array)                     ! Coordinates, stored as x1,x2,(x3) for each node in turn
    real( prec ), intent( in )    :: dof_increment(length_dof_array)                        ! DOF increment, stored as du1,du2,du3,du4... for each node in turn
    real( prec ), intent( in )    :: dof_total(length_dof_array)                            ! accumulated DOF, same storage as for increment

    real( prec ), intent( in )    :: initial_state_variables(n_state_variables)             ! Element state variables.  Defined in this routine
    real( prec ), intent( in )    :: updated_state_variables(n_state_variables)             ! State variables at end of time step
    real( prec ), intent( in )    :: element_properties(n_properties)                       ! Element or material properties, stored in order listed in input file
             
    real( prec ), intent( out )   :: nodal_fieldvariables(n_field_variables,n_nodes)        ! Nodal field variables
  
    ! Local Variables
    logical      :: strcmp
  
    integer      :: n_points,kint,k

    real (prec)  ::  strain(6), dstrain(6)             ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  stress(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  sdev(6)                           ! Deviatoric stress
    real (prec)  ::  D(6,6)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(6,length_dof_array)             ! strain = B*(dof_total+dof_increment)
    real (prec)  :: Bbar(6,length_dof_array)
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  :: elvol
    real (prec)  :: p, smises                          ! Pressure and Mises stress
    !
    !     Subroutine to compute element contribution to project element integration point data to nodes

    !     element_properties(1)         Young's modulus
    !     element_properties(2)         Poisson's ratio

    x = reshape(element_coords,(/3,length_coord_array/3/))

    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

    call initialize_integration_points(n_points, n_nodes, xi, w)

    dNbardx = 0.d0
    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
        dNbardx(1:n_nodes,1:3) = dNdx*w(kint)*determinant
        elvol = w(kint)*determinant
    end do
    ! calculate the element volume
    !write(6,*) 'dof_array = ', dof_array
    dNbardx = dNbardx/elvol

    nodal_fieldvariables = 0.d0

    !     --  Loop over integration points
    do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
        B = 0.d0
        B(1,1:3*n_nodes-2:3) = dNdx(1:n_nodes,1)
        B(2,2:3*n_nodes-1:3) = dNdx(1:n_nodes,2)
        B(3,3:3*n_nodes:3)   = dNdx(1:n_nodes,3)
        B(4,1:3*n_nodes-2:3) = dNdx(1:n_nodes,2)
        B(4,2:3*n_nodes-1:3) = dNdx(1:n_nodes,1)
        B(5,1:3*n_nodes-2:3) = dNdx(1:n_nodes,3)
        B(5,3:3*n_nodes:3)   = dNdx(1:n_nodes,1)
        B(6,2:3*n_nodes-1:3) = dNdx(1:n_nodes,3)
        B(6,3:3*n_nodes:3)   = dNdx(1:n_nodes,2)
        Bbar = 0.d0
        Bbar(1,1:3*n_nodes-2:3) = B(1,1:3*n_nodes-2:3) + 0.3333D0*(dNbardx(1:n_nodes,1) -dNdx(1:n_nodes,1))
        Bbar(1,2:3*n_nodes-1:3) = B(1,2:3*n_nodes-1:3) + 0.3333D0*(dNbardx(1:n_nodes,2)-dNdx(1:n_nodes,2))
        Bbar(1,3:3*n_nodes:3) = B(1,3:3*n_nodes:3) + 0.3333D0*(dNbardx(1:n_nodes,3)-dNdx(1:n_nodes,3))
        Bbar(2,1:3*n_nodes-2:3) = B(2,1:3*n_nodes-2:3) + 0.3333D0*(dNbardx(1:n_nodes,1) -dNdx(1:n_nodes,1))
        Bbar(2,2:3*n_nodes-1:3) = B(2,2:3*n_nodes-1:3) + 0.3333D0*(dNbardx(1:n_nodes,2)-dNdx(1:n_nodes,2))
        Bbar(2,3:3*n_nodes:3) = B(2,3:3*n_nodes:3) + 0.3333D0*(dNbardx(1:n_nodes,3)-dNdx(1:n_nodes,3))
        Bbar(3,1:3*n_nodes-2:3) = B(3,1:3*n_nodes-2:3) + 0.3333D0*(dNbardx(1:n_nodes,1) -dNdx(1:n_nodes,1))
        Bbar(3,2:3*n_nodes-1:3) = B(3,2:3*n_nodes-1:3) + 0.3333D0*(dNbardx(1:n_nodes,2)-dNdx(1:n_nodes,2))
        Bbar(3,3:3*n_nodes:3) = B(3,3:3*n_nodes:3) + 0.3333D0*(dNbardx(1:n_nodes,3)-dNdx(1:n_nodes,3))
        !Bbar(2,2:3*n_nodes-1:3) = B(2,2:3*n_nodes-1:3) + 0.33D0*(dNbardx(1:n_nodes,2)-dNdx(1:n_nodes,2))
        !Bbar(3,3:3*n_nodes:3)   = B(3,3:3*n_nodes:3) + 0.33D0*(dNbardx(1:n_nodes,3)-dNdx(1:n_nodes,3))
        Bbar(4,1:3*n_nodes-2:3) = B(4,1:3*n_nodes-2:3)
        Bbar(4,2:3*n_nodes-1:3) = B(4,2:3*n_nodes-1:3)
        Bbar(5,1:3*n_nodes-2:3) = B(5,1:3*n_nodes-2:3)
        Bbar(5,3:3*n_nodes:3)   = B(5,3:3*n_nodes:3)
        Bbar(6,2:3*n_nodes-1:3) = B(6,2:3*n_nodes-1:3)
        Bbar(6,3:3*n_nodes:3)   = B(6,3:3*n_nodes:3)

        strain = matmul(B,dof_total)
        dstrain = matmul(B,dof_increment)
        call hypoelasticmaterial(strain+dstrain,element_properties,n_properties,stress,D)
        !stress = matmul(D,strain+dstrain)
        p = sum(stress(1:3))/3.d0
        sdev = stress
        sdev(1:3) = sdev(1:3)-p
        smises = dsqrt( dot_product(sdev(1:3),sdev(1:3)) + 2.d0*dot_product(sdev(4:6),sdev(4:6)) )*dsqrt(1.5d0)
        ! In the code below the strcmp( string1, string2, nchar) function returns true if the first nchar characters in strings match
        do k = 1,n_field_variables
            if (strcmp(field_variable_names(k),'S11',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(1)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S22',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(2)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S33',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(3)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S12',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(4)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S13',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(5)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S23',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(6)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'SMISES',6) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + smises*N(1:n_nodes)*determinant*w(kint)
            endif
        end do
 
    end do
  
    return
end subroutine fieldvars_linelast_3dbasic

subroutine hypoelasticmaterial(strain,element_properties,n_properties,stress,D)

    use Types
    use ParamIO

    implicit none

    integer, intent (in) :: n_properties

    real (prec), intent (in) :: strain(6)
    real (prec), intent (in) :: element_properties(n_properties)

    real (prec), intent (out) :: stress(6)
    real (prec), intent (out) :: D(6,6)

    real (prec)  ::  e_dyadic_e(6,6)
    real (prec)  ::  e_ij(6)                           ! e, [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  :: sigzero, epzero, hypo_n, K_in              ! Material properties
    real (prec)  :: elvol, ep_kk, ep_e, sig_e
    real (prec)  :: E_t, E_s


    sigzero = element_properties(1)
    !write(6,*) 'dof_total = ', dof_total
    epzero = element_properties(2)
    !write(6,*) 'epzero = ', epzero
    hypo_n = element_properties(3)
    !write(6,*) 'n = ', hypo_n
    K_in = element_properties(4)
    !write(6,*) 'K = ', K_in

    D = 0.d0
    stress = 0.d0
        !write(6,*) 'total strain top = ', tot_strain(1)
        ep_kk = strain(1) + strain(2) + strain(3)

        !write(6,*) 'ep_kk = ', ep_kk
        e_ij(1) = strain(1)-(1.d0/3.d0)*ep_kk
        e_ij(2) = strain(2)-(1.d0/3.d0)*ep_kk
        e_ij(3) = strain(3)-(1.d0/3.d0)*ep_kk
        e_ij(4) = strain(4)
        e_ij(5) = strain(5)
        e_ij(6) = strain(6)
        !write(6,*) 'e_ij =',e_ij
        ep_e = sqrt((2.d0/3.d0)*(e_ij(1)**2.d0+e_ij(2)**2.d0+e_ij(3)**2.d0+e_ij(4)**2.d0+e_ij(5)**2.d0+e_ij(6)**2.d0))

        if (ep_e==0) then

            sig_e = sigzero*(sqrt(((1.d0+hypo_n**2.d0)/(hypo_n-1)**2.d0)-(hypo_n/(hypo_n-1.d0)&
                -ep_e/epzero)**2.d0)-1.d0/(hypo_n-1.d0))

                stress(1) = 0 !(2.d0/3.d0)*(sig_e/ep_e)*e_ij(1)+K_in*ep_kk
                stress(2) = 0 !(2.d0/3.d0)*(sig_e/ep_e)*e_ij(2)+K_in*ep_kk
                stress(3) = 0 !(2.d0/3.d0)*(sig_e/ep_e)*e_ij(3)+K_in*ep_kk
                stress(4) = 0 !(2.d0/3.d0)*(sig_e/ep_e)*e_ij(4)
                stress(5) = 0 !(2.d0/3.d0)*(sig_e/ep_e)*e_ij(5)
                stress(6) = 0 !(2.d0/3.d0)*(sig_e/ep_e)*e_ij(6)

                E_s = -sigzero*(hypo_n/(hypo_n-1.d0)-ep_e/epzero)/epzero*sqrt(((1.d0+hypo_n**2.d0)/(hypo_n-1.d0)**2.d0)&
                    -(hypo_n/(hypo_n-1.d0)-ep_e/epzero)**2.d0)

                D = 0.d0
                D(1,1) = (E_s/3.d0)*2.d0+(K_in-(2.d0*E_s/9.d0))
                D(1,2) = (K_in-(2.d0*E_s/9.d0))
                D(1,3) = (K_in-(2.d0*E_s/9.d0))
                D(2,1) = (K_in-(2.d0*E_s/9.d0))
                D(2,2) = (E_s/3.d0)*2.d0+(K_in-(2.d0*E_s/9.d0))
                D(2,3) = (K_in-(2.d0*E_s/9.d0))
                D(3,1) = (K_in-(2.d0*E_s/9.d0))
                D(3,2) = (K_in-(2.d0*E_s/9.d0))
                D(3,3) = (K_in-(2.d0*E_s/9.d0))
                D(4,4) = (E_s/3.d0)
                D(5,5) = (E_s/3.d0)
                D(6,6) = (E_s/3.d0)
                !write(6,*) 'epsilon was zero ',  hypo_n
        else

            if (ep_e < epzero)  sig_e = sigzero*(sqrt(((1.d0+hypo_n**2.d0)/(hypo_n-1)**2.d0)-(hypo_n/(hypo_n-1.d0)&
                -ep_e/epzero)**2.d0)-1.d0/(hypo_n-1.d0))
            if (ep_e > epzero)  sig_e = sigzero*(ep_e/epzero)**(1.d0/hypo_n)
            !write(6,*) 'sig_e = ', sig_e
            stress(1) = (2.d0/3.d0)*(sig_e/ep_e)*e_ij(1)+K_in*ep_kk
            stress(2) = (2.d0/3.d0)*(sig_e/ep_e)*e_ij(2)+K_in*ep_kk
            stress(3) = (2.d0/3.d0)*(sig_e/ep_e)*e_ij(3)+K_in*ep_kk
            stress(4) = (2.d0/3.d0)*(sig_e/ep_e)*e_ij(4)
            stress(5) = (2.d0/3.d0)*(sig_e/ep_e)*e_ij(5)
            stress(6) = (2.d0/3.d0)*(sig_e/ep_e)*e_ij(6)

            if (ep_e < epzero) then
                E_t = -sigzero*(hypo_n/(hypo_n-1.d0)-ep_e/epzero)/epzero*sqrt(((1.d0+hypo_n**2.d0)/(hypo_n-1.d0)**2.d0)&
                -(hypo_n/(hypo_n-1.d0)-ep_e/epzero)**2.d0)
            else
                E_t = (sigzero/hypo_n)*(ep_e/epzero)**((1.d0/hypo_n)-1.d0)
            endif
           ! write(6,*) 'E_t = ', E_t
            E_s = sig_e/ep_e
            !write(6,*) 'E_s = ', E_s
            D = 0.d0
            e_dyadic_e = 0.d0
            e_dyadic_e = spread(e_ij,dim=2,ncopies=6)*spread(e_ij,dim=1,ncopies=6)
            D(1,1) = (4.d0/(9.d0*ep_e**2.d0))*(E_t-E_s)*e_dyadic_e(1,1)+(E_s/3.d0)*2.d0+(K_in-(2.d0*E_s/9.d0))
            D(1,2) = (4.d0/(9.d0*ep_e**2.d0))*(E_t-E_s)*e_dyadic_e(1,2)+(K_in-(2.d0*E_s/9.d0))
            D(1,3) = (4.d0/(9.d0*ep_e**2.d0))*(E_t-E_s)*e_dyadic_e(1,3)+(K_in-(2.d0*E_s/9.d0))
            D(2,1) = (4.d0/(9.d0*ep_e**2.d0))*(E_t-E_s)*e_dyadic_e(2,1)+(K_in-(2.d0*E_s/9.d0))
            D(2,2) = (4.d0/(9.d0*ep_e**2.d0))*(E_t-E_s)*e_dyadic_e(2,2)+(E_s/3.d0)*2.d0+(K_in-(2.d0*E_s/9.d0))
            D(2,3) = (4.d0/(9.d0*ep_e**2.d0))*(E_t-E_s)*e_dyadic_e(2,3)+(K_in-(2.d0*E_s/9.d0))
            D(3,1) = (4.d0/(9.d0*ep_e**2.d0))*(E_t-E_s)*e_dyadic_e(3,1)+(K_in-(2.d0*E_s/9.d0))
            D(3,2) = (4.d0/(9.d0*ep_e**2.d0))*(E_t-E_s)*e_dyadic_e(3,2)+(K_in-(2.d0*E_s/9.d0))
            D(3,3) = (4.d0/(9.d0*ep_e**2.d0))*(E_t-E_s)*e_dyadic_e(3,3)+(E_s/3.d0)*2.d0+(K_in-(2.d0*E_s/9.d0))
            D(4,4) = (4.d0/(9.d0*ep_e**2.d0))*(E_t-E_s)*e_dyadic_e(4,4)+(E_s/3.d0)
            D(5,5) = (4.d0/(9.d0*ep_e**2.d0))*(E_t-E_s)*e_dyadic_e(5,5)+(E_s/3.d0)
            D(6,6) = (4.d0/(9.d0*ep_e**2.d0))*(E_t-E_s)*e_dyadic_e(6,6)+(E_s/3.d0)
            D(1,4) = (4.d0/(9.d0*ep_e**2.d0))*(E_t-E_s)*e_dyadic_e(1,4)
            D(1,5) = (4.d0/(9.d0*ep_e**2.d0))*(E_t-E_s)*e_dyadic_e(1,5)
            D(1,6) = (4.d0/(9.d0*ep_e**2.d0))*(E_t-E_s)*e_dyadic_e(1,6)
            D(2,4) = (4.d0/(9.d0*ep_e**2.d0))*(E_t-E_s)*e_dyadic_e(2,4)
            D(2,5) = (4.d0/(9.d0*ep_e**2.d0))*(E_t-E_s)*e_dyadic_e(2,5)
            D(2,6) = (4.d0/(9.d0*ep_e**2.d0))*(E_t-E_s)*e_dyadic_e(2,6)
            D(3,4) = (4.d0/(9.d0*ep_e**2.d0))*(E_t-E_s)*e_dyadic_e(3,4)
            D(3,5) = (4.d0/(9.d0*ep_e**2.d0))*(E_t-E_s)*e_dyadic_e(3,5)
            D(3,6) = (4.d0/(9.d0*ep_e**2.d0))*(E_t-E_s)*e_dyadic_e(3,6)
            D(4,1) = (4.d0/(9.d0*ep_e**2.d0))*(E_t-E_s)*e_dyadic_e(4,1)
            D(4,2) = (4.d0/(9.d0*ep_e**2.d0))*(E_t-E_s)*e_dyadic_e(4,2)
            D(4,3) = (4.d0/(9.d0*ep_e**2.d0))*(E_t-E_s)*e_dyadic_e(4,3)
            D(4,5) = (4.d0/(9.d0*ep_e**2.d0))*(E_t-E_s)*e_dyadic_e(4,5)
            D(4,6) = (4.d0/(9.d0*ep_e**2.d0))*(E_t-E_s)*e_dyadic_e(4,6)
            D(5,1) = (4.d0/(9.d0*ep_e**2.d0))*(E_t-E_s)*e_dyadic_e(5,1)
            D(5,2) = (4.d0/(9.d0*ep_e**2.d0))*(E_t-E_s)*e_dyadic_e(5,2)
            D(5,3) = (4.d0/(9.d0*ep_e**2.d0))*(E_t-E_s)*e_dyadic_e(5,3)
            D(5,4) = (4.d0/(9.d0*ep_e**2.d0))*(E_t-E_s)*e_dyadic_e(5,4)
            D(5,6) = (4.d0/(9.d0*ep_e**2.d0))*(E_t-E_s)*e_dyadic_e(5,6)
            D(6,1) = (4.d0/(9.d0*ep_e**2.d0))*(E_t-E_s)*e_dyadic_e(6,1)
            D(6,2) = (4.d0/(9.d0*ep_e**2.d0))*(E_t-E_s)*e_dyadic_e(6,2)
            D(6,3) = (4.d0/(9.d0*ep_e**2.d0))*(E_t-E_s)*e_dyadic_e(6,3)
            D(6,4) = (4.d0/(9.d0*ep_e**2.d0))*(E_t-E_s)*e_dyadic_e(6,4)
            D(6,5) = (4.d0/(9.d0*ep_e**2.d0))*(E_t-E_s)*e_dyadic_e(6,5)
        endif


    return
 end subroutine hypoelasticmaterial
