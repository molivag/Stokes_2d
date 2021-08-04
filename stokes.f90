program Stokes
  use library
  use Parameters
  use Isoparametric


  implicit none
  ! - - - - - - - - - - * * * Variables que se usan aqui en main * * * * * * * - - - - - - - - - -
  double precision, allocatable, dimension(:,:) :: A_K, Sv, N, dN_dxi, dN_deta
  real, allocatable, dimension(:,:)             :: Fbcsvp
  integer                                       :: NoBV, NoBVcol
  external                                      :: SLEEP
  
  !========== S O L V E R ==========
  external                               :: mkl_dgetrfnp, dgetrf, dgetrs
  integer(4)                             :: S_m, S_n, S_lda, S_ldb, S_infoLU, S_nrhs , S_infoSOL
  integer, allocatable, dimension(:,:)   :: S_ipiv
  character*1                            :: S_trans
  ! - - - - - - - - - - - - - - - * * * Fin * * * * * * * - - - - - - - - - - - - - - - - 
  
  call GeneralInfo( )
  call ReadIntegerFile(20,"elements.dat", Nelem, nUne + 1, elements)  
  call ReadRealFile(10,"nodes.dat", n_nodes,3, nodes) !Para dreducir el numero de subrutinas, usar la sentencia option para 
  call ReadReal(30,"materials.dat", materials)    !Para dreducir el numero de subrutinas, usar la sentencia option para      
  call ReadIntegerFile(40,"pnodes.dat", n_nodes,2, pnodes)
  call ReadIntegerFile(50,"pelements.dat", Nelem,nPne + 1, pelements)
  call GetQuadGauss(ngp,ngp,gauss_points, gauss_weights)
  call ShapeFunctions(gauss_points, nUne, N, dN_dxi, dN_deta)

  allocate(A_K(2*n_nodes+n_pnodes, 2*n_nodes+n_pnodes))
  call SetBounCond( NoBV, NoBVcol) !Esta funcion crea el archivo bcsVP.dat
  allocate( Fbcsvp(NoBV, NoBVcol) ) !Designo la memoria para la matriz de nodos con valor en la frontera
  call ReadMixFile(60,"Fbcsvp.dat", NoBV, NoBVcol, Fbcsvp)!Llamo el archivo de valores en la frontera y lo guardo en Fbcsvp

  call GlobalK( A_K, totGp, dN_dxi, dN_deta)

  DEALLOCATE( N)
  DEALLOCATE( dN_dxi)
  DEALLOCATE( dN_deta)
  
  allocate(Sv(2*n_nodes+n_pnodes, 1))
  Sv = 0 !initializing source vector (Sv) 
  
  
  !========== Una vez calculada la matriz global y el vector de fuente (Sv), les aplicamos las condiciones
  ! de frontera esta subrutina anterior usa como input Sv y A_K y los entrega de nuevo con las BCS aplicadas 
  call ApplyBoundCond(NoBV, Fbcsvp, A_K, Sv )
  !Despues de este ultimo call, obtenemos la matriz y vector global con condiciones de frontera
  DEALLOCATE( Fbcsvp)

  call writeMatrix(A_K, Sv)
  stop


  print*,' '
  print*,'!==================== S O L V E R (L A P A C K) ====================!'
  S_m   = size(A_K,1)
  S_n   = size(A_K,2)
  S_lda = max(1,size(A_K,1))
  S_trans = 'N'
  S_nrhs  = size(Sv,1)
  S_ldb = max(1,size(Sv,1))
  allocate( S_ipiv( 1, max(1,min(S_m, S_n)) ) )

  print*,''
  print*,'= = = = = Solver Parameters = = = = = = = = ='
  print*,'* Shape of S_ipiv  ',shape(S_ipiv), '|'
  print*,'* Leading dimension of A_K     ',S_lda, '|'
  print*,'* Leading dimension of Sv      ',S_ldb, '|'
  print*,'* Number of right-hand sides   ',S_nrhs, '|'
  print*,'= = = = = = = = = = = = = = = = = = = = = = ='
  print*,''
  
  print*,'!=== INITIALIZING LU FACTORIZATION A = P*L*U'
  ! call sleep(2)
  call dgetrf( S_m, S_n, A_K, S_lda, S_ipiv, S_infoLU )
  ! call mkl_dgetrfnp( S_m, S_n, A_K, S_lda, S_infoLU )
  call MKLfactoResult(S_infoLU)
  ! stop
  print*, '.'
  print*,'!========== SOLVING SYSTEM OF EQUATIONS '
  ! call sleep(2)
  call dgetrs( S_trans, S_n, S_nrhs, A_K, S_lda, S_ipiv, Sv, S_ldb, S_infoSOL )
  call MKLsolverResult (S_infoSOL)
  
  call writeMatrix(A_K, Sv)

    ! ========== Escribir en archivo la matriz global
  ! mrow = 2*n_nodes+n_pnodes 
  ! ncol = 2*n_nodes+n_pnodes
  ! open(unit=70, file='A_K_LU.dat', ACTION="write", STATUS="replace ")
  ! do i=1,2*n_nodes+n_pnodes 
  !   write(70, '(1000F20.7)')( A_K(i,j) ,j=1,2*n_nodes+n_pnodes)
  ! end do
  ! close(70)

  ! mrow = 2*n_nodes+n_pnodes 
  ! ncol = 2*n_nodes+n_pnodes
  ! open(unit=71, file='Sv.dat', ACTION="write", STATUS="replace ")
  ! do i=1,2*n_nodes+n_pnodes 
  !   write(71, '(1000F20.7)') Sv(i,1)
  ! end do
  ! close(71)


  ! DEALLOCATE( N)
  ! DEALLOCATE( dN_dxi)
  ! DEALLOCATE( dN_deta)
  ! DEALLOCATE( Fbcsvp)
  ! DEALLOCATE( S_ipiv)
  ! DEALLOCATE( A_K)
  ! DEALLOCATE( Sv )
  

end program Stokes
