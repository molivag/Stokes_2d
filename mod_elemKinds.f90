module Isoparametric
  implicit none

  ! ! ! Global variables related with Isoparametric module ! ! !

  integer :: totGp
  
  ! ! ! End of GLobal Variables ! ! !
  
  contains
  
    subroutine GetQuadGauss(fila, columna, gauss_points, gauss_weights, totGp)!, xi, eta)
      implicit none

      integer               :: i,j,k
      integer, intent(out)  :: totGp
      integer, intent(in)   :: fila, columna
      double precision, allocatable, dimension(:,:),intent(out) :: gauss_points, gauss_weights
      double precision, allocatable, dimension(:,:) :: w1, w2, w, x
      double precision, allocatable, dimension(:,:) :: xi, eta

      allocate(gauss_points(fila*columna,2),gauss_weights(fila*columna,1))
      allocate(w1(columna,1),w2(1,columna), x(columna,1))
      allocate(w(columna,columna))
      allocate(xi(size(gauss_points,1),1),eta(size(gauss_points,1),1))

      totGp = size(gauss_points,1) ! y declarar totGp con SAVE para tener siempre el valor de la variable totGp 
      gauss_points = 0
      gauss_weights = 0

      if ( fila == 1 .and. columna == 1 ) then
        x = 0.0
        w = 4.0
        else if (fila == 2 .and. columna == 2 ) then
          x = reshape([-1.0/sqrt(3.0), 1.0/sqrt(3.0) ], [columna,1])
          w1 = reshape([1.0, 1.0],[columna,1])
          w2 = reshape([1.0, 1.0],[1,columna])
          w = matmul(w1,w2)
        else if (fila == 3 .and. columna == 3 ) then
          x = reshape([-sqrt(3.0/5), 0.0, sqrt(3.0/5)], [columna,1])
          w1 = reshape([5.0/9, 8/9.0, 5.0/9],[columna,1])
          w2 = reshape([5.0/9, 8/9.0, 5.0/9],[1,columna])
          w = matmul(w1,w2)
        else
        print*, 'Error calling GetQuadGauss\n'
      end if

      !for-loop block: set up the 2-D Gauss points and weights
      k=1
      do i=1,fila
        do j=1, columna
          gauss_points(k,1:2)= [x(i,1), x(j,1)]
          gauss_weights(k,1) = w(i,j)     !En fortran un vector debe especificarse
          k=k+1                           !como un arreglo de rango 2 es decir (1,n) o (n,1)
        end do                            !mientras que en matlab solo con escribir (n)
      end do

      ! xi  = gauss_points(:,1)
      ! eta = gauss_points(:,2)
      xi(:,1)  = gauss_points(:,1)     ! xi-coordinate of point j
      eta(:,1) = gauss_points(:,2)

      DEALLOCATE(w1, w2, x, w)
      !Esta funcion no afecta al resultado pues se ha liberado la memoria para calcular pesos y puntos
      !de Gauss mas no las variables que contienen pesos y puntos de gauus.
    end subroutine GetQuadGauss



  
    subroutine Quad8Nodes(gauss_points,  N, dN_dxi, dN_deta)
      implicit None

      ! Shape functions for square (quadrilaters) linear elements
			!
			!  |
			!  |
			!  | o- - o - -o
			!  Y |         |
			!  | o         o
			!  | |         |
			!  | o- - o - -o
			!  |
			!  +--------X-------->
      

      integer, parameter  :: Nne = 8          !Esto se puede quitar poniendo un modulo con todos las variables globales incluyendo las de mod_library
      integer, parameter  :: dim_prob = 2

      ! integer,                          intent(in) :: totGp
      double precision, dimension(:,:), intent(in) :: gauss_points
      double precision, allocatable, dimension(:,:), intent(out) :: N, dN_dxi, dN_deta
      double precision, dimension(size(gauss_points,1)) :: xi_vector, eta_vector
      integer, dimension(Nne,dim_prob) :: master_nodes
      double precision    :: xi, eta, mn_xi, mn_eta
      integer             :: i, j, jj, k

      ! !number of gauss points
      ! totGp = size(gauss_points,1) ! esta puede quedar como variable global si se usa en alguna otra subrutina
      !                           ! si solo se usa aqui, entonces variable como local-----> Si se usa en otra rutina, en compK
      allocate( N(Nne,totGp),dN_dxi(Nne,totGp),dN_deta(Nne,totGp) )

      N = 0.0
      dN_dxi  = 0.0
      dN_deta = 0.0

      xi_vector  = gauss_points(:,1)     ! xi-coordinate of point j
      eta_vector = gauss_points(:,2)

      !coordinates of the nodes of the master element
      master_nodes = reshape([1, -1, -1, 1, 0, -1, 0, 1, 1, 1, -1, -1, 1, 0, -1, 0], [Nne,dim_prob])
      !NOTA ** Para que el reshape funcione correctamente, o que produzca el par de valores deseado, primero se deben
      !colocar todos los valores en x, luego todos los de y y luego, si hubiera todos los de z para que al acomodarse salga el par
      !suponiendo un reshape de 3,2 debe acomodarse x1, x2, x3, y1, y2, y3 DUDA *Siempre es asi*

      ! dN(xi,eta)/dx = dN/dxi(dxi\dx) + dN/deta(deta/dx)
      ! dN(xi,eta)/dy = dN/dxi(dxi\dy) + dN/deta(deta/dy)
      ! Aqui se calculan as funciones de forma N y parte de las derivadas dN/dxi and dN_deta
      ! mas no las derivadas dN/dx and dN/dy completas

      do j = 1, totGp
        xi  = xi_vector(j)      ! xi-coordinate of point j
        eta = eta_vector(j)     ! eta-coordinate of point j

        N(5,j)=1.0/2*(1-xi**2)*(1+eta)
        N(6,j)=1.0/2*(1-xi)*(1-eta**2)
        N(7,j)=1.0/2*(1-xi**2)*(1-eta)
        N(8,j)=1.0/2*(1+xi)*(1-eta**2)
        dN_dxi(5,j)=-xi*(1+eta)
        dN_dxi(6,j)=-1.0/2*(1-eta**2)
        dN_dxi(7,j)=-xi*(1-eta)
        dN_dxi(8,j)=1.0/2*(1-eta**2)
        dN_deta(5,j)=1.0/2*(1-xi**2)
        dN_deta(6,j)=(1-xi)*(-eta)
        dN_deta(7,j)=-1.0/2*(1-xi**2)
        dN_deta(8,j)=(1+xi)*(-eta)

        do i = 1, 4
          mn_xi = master_nodes(i,1)
          mn_eta= master_nodes(i,2)
          if (i==1) then
            jj=8
          else
            jj=i+3
          end if
          k=i+4
          N(i,j)=(1.0 + mn_xi*xi)*(1.0 + mn_eta*eta)/4.0 - 1.0/2*(N(jj,j)+N(k,j))
          dN_dxi(i,j)= mn_xi*(1.0 + mn_eta*eta)/4.0 - 1.0/2*(dN_dxi(jj,j)+dN_dxi(k,j))
          dN_deta(i,j)= mn_eta*(1.0 + mn_xi*xi)/4.0 - 1.0/2*(dN_deta(jj,j)+dN_deta(k,j))

        end do

      end do

    end subroutine Quad8Nodes


    subroutine Quad4Nodes (gauss_points,  N, dN_dxi, dN_deta)
      implicit None

      !CompNDNatPointsQuad8
      ! Shape functions for square (quadrilaters) linear elements
			!
			!  |
			!  |
			!  | o- - - -o
			!  Y |       |
			!  | |       |
			!  | o- - - -o
			!  |
			!  +--------X-------->
      

      integer, parameter  :: Nne = 4        !Esto se puede quitar poniendo un modulo con todos las variables globales incluyendo las de mod_library
      integer, parameter  :: dim_prob = 2

      ! integer,                          intent(in) :: totGp
      double precision, dimension(:,:), intent(in) :: gauss_points
      double precision, allocatable, dimension(:,:), intent(out) :: N
      double precision, allocatable, dimension(:,:), intent(out), optional :: dN_dxi, dN_deta
      double precision, dimension(size(gauss_points,1)) :: xi_vector, eta_vector
      integer, dimension(Nne,dim_prob) :: master_nodes
      double precision    :: xi, eta, mn_xi, mn_eta
      integer             :: i, j

      !si solo se usa aqui, entonces variable como local-----> Si se usa en otra rutina, en compK
      !Al parecer este alojamiento de memoria evita el error despues de ejecutar LAPACK
      


      xi_vector  = gauss_points(:,1)     ! xi-coordinate of point j
      eta_vector = gauss_points(:,2)

      !coordinates of the nodes of the master element
      master_nodes = reshape([1, -1, -1, 1, 1, 1, -1, -1], [Nne,dim_prob])
      !NOTA ** Para que el reshape funcione correctamente, o que produzca el par de valores deseado, primero se deben
      !colocar todos los valores en x, luego todos los de y y luego, si hubiera todos los de z para que al acomodarse salga el par
      !suponiendo un reshape de 3,2 debe acomodarse x1, x2, x3, y1, y2, y3 DUDA *Siempre es asi*

      ! dN(xi,eta)/dx = dN/dxi(dxi\dx) + dN/deta(deta/dx)
      ! dN(xi,eta)/dy = dN/dxi(dxi\dy) + dN/deta(deta/dy)
      ! Aqui se calculan as funciones de forma N y parte de las derivadas dN/dxi and dN_deta
      ! mas no las derivadas dN/dx and dN/dy completas
           
      !do loop: compute N, dN_dxi, dN_deta
      if (present(dN_dxi) .and. present(dN_deta))then
        allocate(dN_dxi(Nne,totGp) )
        allocate(dN_deta(Nne,totGp) )
        dN_dxi  = 0.0
        dN_deta = 0.0
        
        do j=1,totGp                              ! columns for point 1,2 ...
          xi=xi_vector(j);                      ! xi-coordinate of point j 
          eta=eta_vector(j);                    ! eta-coordinate of point j 
          do i=1,4                              ! rows for N1, N2, ...
            mn_xi = master_nodes(i,1)
            mn_eta= master_nodes(i,2)
            dN_dxi(i,j)= mn_xi*(1.0 + mn_eta*eta)/4.0             ! dNi/dxi(xi,eta)
            dN_deta(i,j)= mn_eta*(1.0 + mn_xi*xi )/4.0            ! dNi/deta(xi,eta
          end do
        end do
      else
        allocate(N(Nne,totGp) )
        N     = 0.0
        do j=1,totGp                              ! columns for point 1,2 ...
          xi=xi_vector(j);                      ! xi-coordinate of point j 
          eta=eta_vector(j);                    ! eta-coordinate of point j 
          do i=1,4                              ! rows for N1, N2, ...
            mn_xi = master_nodes(i,1)
            mn_eta= master_nodes(i,2)
            N(i,j)=(1.0 + mn_xi*xi)*(1.0 + mn_eta*eta)/4.0        ! Ni(xi,eta)
          end do
        end do

      endif


      

    end subroutine Quad4Nodes



  	! != = = = = = = = = = = = = = = = = = = = = = =
    subroutine Triang3pElem( gauss_points, N, dN_dxi, dN_deta)
      implicit none

			!  
			!  |
			!  |        o
			!  |       / \
			!  |      /   \
			!  Y     /     \
			!  |    /       \
			!  |   /         \
			!  |  o-----------o
			!  |
			!  +--------X-------->



      integer, parameter  :: Nne = 4
      integer, parameter  :: dim_prob = 2   !Esto se puede quitar poniendo un modulo con todos las variables globales incluyendo las de mod_library

      double precision, dimension(:,:), intent(in) :: gauss_points
      double precision, allocatable, dimension(:,:), intent(out) :: N
      double precision, allocatable, dimension(:,:), intent(out), optional :: dN_dxi, dN_deta
      double precision, dimension(size(gauss_points,1)) :: xi_vector, eta_vector
      double precision    :: xi, eta
      integer             :: j


      allocate(dN_dxi(Nne,totGp) )
      allocate(dN_deta(Nne,totGp) )
      N = 0.0
      dN_dxi  = 0.0
      dN_deta = 0.0

      xi_vector  = gauss_points(:,1)     ! xi-coordinate of point j
      eta_vector = gauss_points(:,2)

      do j=1,totGp
        xi=xi_vector(j);                      ! xi-coordinate of point j 
        eta=eta_vector(j); 

        N(1,j)      =  2*(1-xi-eta)*(1-xi-eta-0.5)
        dN_dxi(1,j) = -2*(1-xi-eta-0.5)-2*(1-xi-eta) 
        dN_deta(1,j)= -2*(1-xi-eta-0.5)-2*(1-xi-eta) 

        N(2,j)      = 2*xi*(xi-0.5)
        dN_dxi(2,j) = 2*(xi-0.5) + 2*xi 
        dN_deta(2,j)= 0

        N(3,j)      =  2*eta*(eta-0.5)
        dN_dxi(3,j) = 0
        dN_deta(3,j)= 2*(eta-0.5) + 2*eta 

      end do
      
      

    end subroutine Triang3pElem



  	



  !end contains
    
end module Isoparametric


