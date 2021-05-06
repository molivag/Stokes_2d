program Stokes
  use library

  implicit none


  integer :: i,j,  el_num, Gp !Estas variables declaradas solo son de prueba para ir testeando la funcionalidad del codigom, se cambiaran por el bucle principal en compK
  real, allocatable, dimension(:,:) :: A_K
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
  call ReadRealFile(40,"pnodes.dat", 341,2, pnodes)
  call ReadIntegerFile(50,"pelements.dat", 100,5, pelements)
  call GetQuadGauss(2,2,gauss_points, gauss_weights)

  allocate( N(Nne,size(gauss_points,1)),Nx(Nne,size(gauss_points,1)) ,Ny(Nne,size(gauss_points,1)) )

  call CompNDNatPointsQuad8(gauss_points, N, Nx, Ny)

  allocate(A_K(2*n_nodes+n_pnodes, 2*n_nodes+n_pnodes))
  

  el_num = 1
  call SetElementNodes(el_num, elements, nodes, element_nodes, node_id_map)

  Gp = 1
  ! Jaco = J2D( Nx, Ny, Gp)
  ! Jinv = inv2x2(Jaco)

  ! MatH = compH()

  ! B = compBmat(Nx, Ny, Gp)

  ! Jb = buildJb(Jinv)

  ! HJ = matmul(MatH,Jb)
  ! HJB = matmul(HJ,B)

  call GlobalK( A_K, Nx, Ny)



  do i =1,100
    print*, A_K(i,1)
  end do

  




  ! print*,' '
  ! print*,'shape de Jaco',shape(Jaco)
  ! print*,' '
  ! do i =1 ,dim_prob
  !   print*, Jaco(i,:)
  ! end do

  ! print*,' '
  ! print*,'shape de Jinv',shape(Jinv)
  ! print*,' '
  ! do i =1 ,dim_prob
  !  print*, Jinv(i,:)
  ! end do


  ! print*,' '
  ! print*,'shape de H',shape(MatH)
  ! print*,' '
  ! do i =1 ,3
  !   print*, MatH(i,:)
  ! end do


  ! print*,' '
  ! print*,'shape de Jb',shape(Jb)
  ! print*,' '
  ! do i =1 , 4
  !   print*, Jb(i,:)
  ! end do
  ! print*,' '

  ! print*,' '
  ! print*,'shape de B',shape(B)
  ! print*,' '
  ! do i =1 ,dim_prob*2
  !   print*, B(i,:)
  ! end do

  ! print*,' '
  ! print*,'shape de HJB',shape(HJB)
  ! print*,' '
  ! print*,' '

  ! do i =1 ,3
  !    print*, HJB(i,:)
  ! end do




end program Stokes
