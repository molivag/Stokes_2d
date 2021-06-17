program Stokes
  use library

  implicit none


  integer :: i !Estas variables declaradas solo son de prueba para ir testeando la funcionalidad del codigom, se cambiaran por el bucle principal en compK
  real(8), allocatable, dimension(:,:) :: A_K
  real, allocatable, dimension(:,:) :: N, Nx, Ny

  ! real, dimension(2,2) :: Jaco

! - - - - - Aqui declaro funciones para ir probando el codigo
  ! real, dimension(3,4)                          :: MatH
  ! real, dimension(dim_prob,dim_prob)            :: Jaco, Jinv
  ! real, dimension(4,2*n_nodes_per_element)      :: B
  ! real, dimension(2*dim_prob, 2*dim_prob)       :: Jb
  ! real, dimension(3,dim_prob*dim_prob)          :: HJ
  ! real, dimension(3,2*n_nodes_per_element)      :: HJB
  ! real, dimension(3,2*n_nodes_per_element)     :: matmul3
! - - - - - Aqui declaro funciones para ir probando el codigo

  call ReadRealFile(10,"nodes.dat", 341,3, nodes)
  call ReadIntegerFile(20,"elements.dat", 100,9, elements)
  call ReadReal(30,"materials.dat", materials)
  call ReadIntegerFile(40,"pnodes.dat", 341,2, pnodes)
  call ReadIntegerFile(50,"pelements.dat", 100,5, pelements)
  call GetQuadGauss(2,2,gauss_points, gauss_weights)

  ! allocate( N(Nne,size(gauss_points,1)),Nx(Nne,size(gauss_points,1)) ,Ny(Nne,size(gauss_points,1)) )
  call CompNDNatPointsQuad8(gauss_points, N, Nx, Ny)

  allocate(A_K(2*n_nodes+n_pnodes, 2*n_nodes+n_pnodes))
  
  call GlobalK( A_K, Nx, Ny)

  do i = 445,524
    print*,A_K(i,453)
    end do

end program Stokes
