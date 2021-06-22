program Stokes
  use library

  implicit none
  ! - - - - - - - - - - * * * Variables que se usan aqui en main * * * * * * * - - - - - - - - - -
  integer :: NoBV, NoBVcol, mrow, ncol
  real, allocatable, dimension(:,:) :: Fbcsvp
  ! - - - - - - - - - - - - - - - * * * Fin * * * * * * * - - - - - - - - - - - - - - - - 
  
  integer :: i, j !Estas variables declaradas solo son de prueba para ir testeando la funcionalidad del codigom, se cambiaran por el bucle principal en compK
  real(8), allocatable, dimension(:,:) :: A_K
  real, allocatable, dimension(:,:) :: N, Nx, Ny
  real(8), dimension(2*n_nodes+n_pnodes, 1) :: Sv

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
  
  ! do i = 1,NoBV
  !   print*,Fbcsvp(i,3)
  ! end do

  call GlobalK( A_K, Nx, Ny)

  ! print*, ' '
  ! print*, 'A_K(640,802)',A_K(640,802)
  ! print*, ' '
  ! print*, 'A_K(640,803)',A_K(640,803)
  ! print*, ' '
  ! print*, 'n_nodes',n_nodes
  ! print*, 'n_pnodes',n_pnodes
  ! print*, ' '

  ! do i = 1,803
  !   print*,A_K(i,803)
  ! end do

  Sv = 0 !initializing source vector (Sv) 

  ! una vez calculada la matriz global y el vector de fuente (Sv), les aplicamos las condiciones de frontera
  ! call ApplyBoundCond(A_K, Sv, NoBV, Fbcsvp )

  !Despues de este ultimo call, obtenemos la matriz y vector global con condiciones de frontera
  ! Aqui entraria el solver, este deberia estar en un modulo distinto
  

  ! do j=1,2*n_nodes+n_pnodes
  !   do i = 1,2*n_nodes+n_pnodes
  !     print*,A_K(i,j)
  !   end do
  ! end do

  ! mrow = 2*n_nodes+n_pnodes
  ! ncol = 2*n_nodes+n_pnodes
  
  ! open(unit=2, file='globalK.dat', ACTION="write", STATUS="replace")
  ! do i=1,ncol
  !   write(2, '(1000F14.7)')( A_K(i,j) ,j=1,mrow)
  ! end do

  ! close(2)

end program Stokes
