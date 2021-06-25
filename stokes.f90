program Stokes
  use library

  implicit none
  ! - - - - - - - - - - * * * Variables que se usan aqui en main * * * * * * * - - - - - - - - - -
  real, allocatable, dimension(:,:) :: Fbcsvp
  integer           :: NoBV, NoBVcol, mrow, ncol, i, j 
  !========== S O L V E R
  integer(8)        :: S_lda, S_info, S_m, S_n, S_ipiv
  ! - - - - - - - - - - - - - - - * * * Fin * * * * * * * - - - - - - - - - - - - - - - - 
  
  ! integer :: i, j !Estas variables declaradas solo son de prueba para ir testeando la funcionalidad del codigom, se cambiaran por el bucle principal en compK
  double precision, allocatable, dimension(:,:) :: A_K
  double precision, allocatable, dimension(:,:) :: N, Nx, Ny
  double precision, dimension(2*n_nodes+n_pnodes, 1) :: Sv

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
  
  call GlobalK( A_K, Nx, Ny)

  Sv = 0 !initializing source vector (Sv) 
  !========== Una vez calculada la matriz global y el vector de fuente (Sv), les aplicamos las condiciones
  ! de frontera esta subrutina anterior usa como input Sv y A_K y los entrega de nuevo con las BCS aplicadas 
  call ApplyBoundCond(NoBV, Fbcsvp, A_K, Sv )
  !Despues de este ultimo call, obtenemos la matriz y vector global con condiciones de frontera

  ! !========== S O L V E R (L A P A C K) ==========
  ! !========== LU FACTORIZATION 
  ! S_m   = size(A_K,1)
  ! S_n   = size(A_K,2)
  ! S_lda = max(1,size(A_K,1))

  ! print*,' '
  ! print*,'Leading dimension of A_K',S_lda
  ! print*,' '

  ! call dgetrf( S_m, S_n, A_K, S_lda, S_ipiv, S_info )
  
  !========== Escribir en archivo la matriz global
  
  mrow = 2*n_nodes+n_pnodes 
  ncol = 2*n_nodes+n_pnodes
  open(unit=2, file='A_K.dat', ACTION="write", STATUS="new")
  do i=1,2*n_nodes+n_pnodes 
    write(2, '(1000F20.9)')( A_K(i,j) ,j=1,2*n_nodes+n_pnodes)
  end do
  close(2)

end program Stokes
