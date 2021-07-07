program Stokes
  use library

  implicit none
  ! - - - - - - - - - - * * * Variables que se usan aqui en main * * * * * * * - - - - - - - - - -
  real, allocatable, dimension(:,:) :: Fbcsvp
  integer           :: NoBV, NoBVcol, i, j, mrow, ncol
  !========== S O L V E R
  external :: mkl_dgetrfnp, dgetrf, dgetrs
  integer                                :: S_m, S_n, S_lda, S_ldb, S_infoLU, S_info, S_nrhs
  integer, allocatable, dimension(:,:)   :: S_ipiv
  character*1                            :: S_trans
  ! - - - - - - - - - - - - - - - * * * Fin * * * * * * * - - - - - - - - - - - - - - - - 
  
  ! integer :: i, j !Estas variables declaradas solo son de prueba para ir testeando la funcionalidad del codigom, se cambiaran por el bucle principal en compK
  double precision, allocatable, dimension(:,:) :: A_K
  double precision, allocatable, dimension(:,:) :: N, Nx, Ny
  double precision, allocatable, dimension(:,:) :: Sv

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
  
  call SetBounCond( NoBV, NoBVcol) !Esta funcion crea el archivo bcsVP.dat
  
  allocate( Fbcsvp(NoBV, NoBVcol) ) !Designo la memoria para la matriz de nodos con valor en la frontera
  call ReadMixFile(60,"Fbcsvp.dat", NoBV, NoBVcol, Fbcsvp)!Llamo el archivo de valores en la frontera y lo guardo en Fbcsvp
  print*,''

  call GlobalK( A_K, Nx, Ny)

  allocate(Sv(2*n_nodes+n_pnodes, 1))
  Sv = 0 !initializing source vector (Sv) 
  !========== Una vez calculada la matriz global y el vector de fuente (Sv), les aplicamos las condiciones
  ! de frontera esta subrutina anterior usa como input Sv y A_K y los entrega de nuevo con las BCS aplicadas 
  call ApplyBoundCond(NoBV, Fbcsvp, A_K, Sv )
  !Despues de este ultimo call, obtenemos la matriz y vector global con condiciones de frontera
  
  
  print*,'!==================== S O L V E R (L A P A C K) ====================!'
  print*,'!========== LU FACTORIZATION A = P*L*U'
  S_m   = size(A_K,1)
  S_n   = size(A_K,2)
  S_lda = max(1,size(A_K,1))
  allocate( S_ipiv( 1, max(1,min(S_m, S_n)) ) )
  S_trans = 'N'
  S_nrhs  = size(SV,1)
  S_ldb = max(1,size(Sv,1))
  
  print*,''
  print*,'= = = = = Solver Parameters = = = = = = = = = ='
  print*,'= shape of S_ipiv  ',shape(S_ipiv)
  print*,'= Leading dimension of A_K     ',S_lda
  print*,'= Leading dimension of Sv      ',S_ldb
  print*,'= Number of right-hand sides   ',S_nrhs
  print*,'= = = = = = = = = = = = = = = = = = = = = = = ='
  call dgetrf( S_m, S_n, A_K, S_lda, S_ipiv, S_infoLU )
  ! call mkl_dgetrfnp( S_m, S_n, A_K, S_lda, S_info )
  print*,''
  PRINT *, "!========== Factorization completed "
  PRINT*, 'S_info for LU Factorization', S_infoLU
  PRINT*, 'SHAPE OF A_K_LU',SHAPE(A_K)
  PRINT*, ' '  
  print*,'!========== SOLVING SYSTEM OF EQUATIONS '

  call dgetrs( S_trans, S_n, S_nrhs, A_K, S_lda, S_ipiv, Sv, S_ldb, S_info )

  ! PRINT*, 'S_info for sol. of  Ax=b   ', S_info
  ! PRINT*, 'SHAPE OF A_K_LU',SHAPE(Sv)
  ! PRINT*, ' '  

  ! DEALLOCATE(A_K, Sv, S_ipiv)
  
  ! do i = 1, 803
  !   print*, S_ipiv(1,i)
  ! end do
  
  ! ========== Escribir en archivo la matriz global
  
  ! mrow = 2*n_nodes+n_pnodes 
  ! ncol = 2*n_nodes+n_pnodes
  ! open(unit=70, file='A_K_LU.dat', ACTION="write", STATUS="replace ")
  ! do i=1,2*n_nodes+n_pnodes 
  !   write(70, '(1000F20.9)')( A_K(i,j) ,j=1,2*n_nodes+n_pnodes)
  ! end do
  ! close(70)

  mrow = 2*n_nodes+n_pnodes 
  ncol = 2*n_nodes+n_pnodes
  open(unit=71, file='Sv.dat', ACTION="write", STATUS="replace ")
  do i=1,2*n_nodes+n_pnodes 
    write(71, '(1000F10.7)') Sv(i,1)
  end do
  close(71)

end program Stokes
