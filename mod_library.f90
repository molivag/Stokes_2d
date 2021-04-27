
module library
  Implicit None

  ! ! ! Aqui se declaran las variables globales ! ! !

  real                      :: materials
  real, dimension(341,3)    :: nodes
  integer, dimension(100,9) :: elements
  real, dimension(341,2)    :: pnodes
  integer, dimension(100,5) :: pelements
  integer, parameter        :: n_nodes = size(nodes,1)
  integer, parameter        :: n_pnodes = 121 !Duda, como lo tomo desde el txt es decir   n_pnodes = maxval(pnodes(:,2))
  integer, parameter        :: n_elements = size(elements,1)
  integer, parameter        :: Nne = size(elements,2)-1
  integer, parameter        :: dim_prob = size(nodes,2)-1 !dimension del problema

  real, allocatable, dimension(:,:) :: gauss_points, gauss_weights !Verificar si debe ser global---> Si, se usa en la funcion ComputeK
  real, allocatable, dimension(:,:) :: N, Nx, Ny
  real, dimension(Nne,dim_prob)     :: element_nodes
  integer, dimension(Nne,1)         :: node_id_map

  ! ! ! Fin de variables globales ! ! !


  contains

    subroutine ReadRealFile(UnitNum, FileName, NumRows, NumCols, Real_Array)

      integer :: i, j, status, UnitNum, NumRows, NumCols
      character (len=*), intent (in) :: FileName
      real, dimension (1:NumRows, 1:NumCols), intent (out) :: Real_Array


      open (unit = UnitNum, file =FileName, status='old', action='read' , iostat = status)

      ! read in values
      read(UnitNum,*) ((Real_Array(i,j), j=1,NumCols), i=1,NumRows)
      print *, "Status_Real_File ", status

      close (UnitNum)

    end subroutine

    subroutine ReadIntegerFile(UnitNum, FileName, NumRows, NumCols, IntegerArray)

      integer :: i, j, status, UnitNum, NumRows, NumCols
      character (len=*), intent (in) :: FileName
      integer, dimension (1:NumRows, 1:NumCols), intent (out) :: IntegerArray


      open (unit = UnitNum, file =FileName, status='old', action='read' , iostat = status)

      ! read in values
      read(UnitNum,*) ((IntegerArray(i,j), j=1,NumCols), i=1,NumRows)
      print *, "Status_Int_File  ", status

      close (UnitNum)

    end subroutine

    subroutine ReadReal(UnitNum, FileName, value)

      integer :: status, UnitNum
      character (len=*), intent (in) :: FileName
      real :: value


      open (unit = UnitNum, file =FileName, status='old', action='read' , iostat = status)

      ! read in values
      read(UnitNum,*) value
      print *, "Status_Single_Val", status

      close (UnitNum)

    end subroutine

    subroutine GetQuadGauss(fila, columna, gauss_points, gauss_weights)
      implicit none

      integer,intent(in) :: fila, columna
      real, allocatable, dimension(:,:),intent(out) :: gauss_points, gauss_weights
      integer :: i,j,k
      real, allocatable, dimension(:,:) :: w1, w2, w, x

      allocate(gauss_points(fila*columna,2),gauss_weights(fila*columna,1))
      allocate(w1(columna,1),w2(1,columna), x(columna,1))
      allocate(w(columna,columna))

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

      DEALLOCATE(w1, w2, x, w)
      !Esta funcion no afecta al resultado pues se ha liberado la memoria para calcular pesos y puntos
      !de Gauss mas no las variables que contienen pesos y puntos de gauus.
    end subroutine GetQuadGauss

    subroutine CompNDNatPointsQuad8(gauss_points, N, Nx, Ny)
      implicit none

      real, dimension(:,:), intent(in) :: gauss_points
      real, allocatable, dimension(:,:), intent(out) :: N, Nx, Ny
      real, dimension(size(gauss_points,1)) :: xi_vector, eta_vector
      integer, dimension(Nne,dim_prob) :: master_nodes
      real    :: xi, eta, mn_xi, mn_eta
      integer :: ngp, i, j, jj, k

      !number of gauss points
      ngp = size(gauss_points,1) ! esta puede quedar como variable global si se usa en alguna otra subrutina
                                ! si solo se usa aqui, entonces variable como local-----> Si se usa en otra rutina, en compK
      allocate( N(Nne,ngp),Nx(Nne,ngp),Ny(Nne,ngp) )

      N  = 0.0
      Nx = 0.0
      Ny = 0.0

      xi_vector  = gauss_points(:,1)     ! xi-coordinate of point j
      eta_vector = gauss_points(:,2)

      !coordinates of the nodes of the master element
      master_nodes = reshape([1, -1, -1, 1, 0, -1, 0, 1, 1, 1, -1, -1, 1, 0, -1, 0], [Nne,dim_prob])
      !NOTA ** Para que el reshape funcione correctamente, o que produzca el par de valores deseado, primero se deben
      !colocar todos los valores en x, luego todos los de y y luego, si hubiera todos los de z para que al acomodarse salga el par
      !suponiendo un reshape de 3,2 debe acomodarse x1, x2, x3, y1, y2, y3 DUDA *Siempre es asi*

      do j = 1, ngp
        xi  = xi_vector(j)      ! xi-coordinate of point j
        eta = eta_vector(j)     ! eta-coordinate of point j

        N(5,j)=1.0/2*(1-xi**2)*(1+eta)
        N(6,j)=1.0/2*(1-xi)*(1-eta**2)
        N(7,j)=1.0/2*(1-xi**2)*(1-eta)
        N(8,j)=1.0/2*(1+xi)*(1-eta**2)
        Nx(5,j)=-xi*(1+eta)
        Nx(6,j)=-1.0/2*(1-eta**2)
        Nx(7,j)=-xi*(1-eta)
        Nx(8,j)=1.0/2*(1-eta**2)
        Ny(5,j)=1.0/2*(1-xi**2)
        Ny(6,j)=(1-xi)*(-eta)
        Ny(7,j)=-1.0/2*(1-xi**2)
        Ny(8,j)=(1+xi)*(-eta)

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
          Nx(i,j)= mn_xi*(1.0 + mn_eta*eta)/4.0 - 1.0/2*(Nx(jj,j)+Nx(k,j))
          Ny(i,j)= mn_eta*(1.0 + mn_xi*xi)/4.0 - 1.0/2*(Ny(jj,j)+Ny(k,j))

        end do

      end do

    end subroutine CompNDNatPointsQuad8

    subroutine SetElementNodes(elm_num, element_nodes, node_id_map)

      implicit none

      ! integer, dimension(100,9), intent(in):: elements
      ! real, dimension(341,3), intent(in):: nodes
      integer,intent(in) :: elm_num ! number of element for each elemental integral in the main do ok K global
      real, dimension(Nne,dim_prob), intent(out) ::  element_nodes
      integer, dimension(Nne,1), intent(out)     :: node_id_map

      integer :: i,j, global_node_id


      element_nodes = 0
      node_id_map = 0

      do i = 1, Nne
        global_node_id = elements(elm_num,i+1)
        do j=1 ,dim_prob
          element_nodes(i,j) = nodes(global_node_id,j+1)
        end do
        node_id_map(i,1) = global_node_id;
      end do

    end subroutine SetElementNodes

    function J2D( Nx, Ny, Gp) result(Jacobian)
      implicit none

      ! real, dimension(Nne,dim_prob), intent(in)            :: element_nodes
      real, dimension(Nne,size(gauss_points) ), intent(in) :: Nx, Ny
      integer, intent(in)                                                  :: Gp !esta variable se usara en el lazo principal con el numero de punto de gauss para evaluar las integrales elementales
      real, dimension(dim_prob,Nne)                        :: Basis2D
      real, dimension(1,Nne)                               :: Nxi, Neta
      ! real, dimension(dim_prob,dim_prob)                                   :: J2D
      real, dimension(dim_prob,dim_prob)                                   :: Jacobian

      !con estas instrucciones extraigo la columna de Nx como renglon y lo guardo en Nxi, Gp se
      !ira moviendo conforme la funcion J2D sea llamada en el lazo principal para cada elemento lo mismo para Neta con Ny
      Nxi  = spread(Nx(:,Gp),dim = 1, ncopies= 1)
      Neta = spread(Ny(:,Gp),dim = 1, ncopies= 1)

      Basis2D(1,:) = Nxi(1,:)
      Basis2D(2,:) = Neta(1,:)

      Jacobian = matmul(Basis2D,element_nodes)

      ! J = J2D

      ! - - - * * * D U D A * * * - - -
        !Si el nombre de la funcion es el mismo que la variable donde se guarda, entonces no puedo declararla como
        ! variable global, ¿Como debo hacerlo?

        !Si lo dejo como

        !J = Matmul(Basis2D,element_nodes)

        !Me marca un warning
      ! - - - * * * D U D A * * * - - -



      return
    end function J2D

    function inv2x2(A)  result(Jinv)

      implicit none

      real, dimension(dim_prob,dim_prob), intent(in)       :: A
      real, dimension(dim_prob,dim_prob)     :: Jinv

      double precision, parameter :: EPS = 1.0E-10
      real :: det
      real, dimension(2,2) :: cofactor


      det =   A(1,1)*A(2,2) - A(1,2)*A(2,1)

      if (abs(det) .le. EPS) then
        Jinv = 0.0D0
        return
      end IF

      cofactor(1,1) = +A(2,2)
      cofactor(1,2) = -A(2,1)
      cofactor(2,1) = -A(1,2)
      cofactor(2,2) = +A(1,1)

      Jinv = transpose(cofactor) / det

      return

    end function inv2x2

    function buildJb(A)
      !Funcion que construye una matriz de 4 x 4 en bloques de 2 para el caso 2D

      implicit none

      real, dimension(dim_prob,dim_prob), intent (in)   :: A
      real, dimension(2*dim_prob, 2*dim_prob)           :: buildJb

      buildJb(1:2,1:2) = A
      buildJb(3:4,3:4) = A

    end function

    function m22det(A)

      implicit none
      real :: m22det
      real, dimension(2,2), intent(in)  :: A



      m22det =   A(1,1)*A(2,2) - A(1,2)*A(2,1)

      return

    end function m22det

    function CompH()
      implicit None

      ! integer :: CompH
      ! integer, dimension(3,4) :: H
      integer, dimension(3,4) :: CompH
      CompH = 0
      CompH(1,1)=1;
      CompH(2,4)=1;
      CompH(3,2)=1;
      CompH(3,3)=1;


      ! CompH = H

      ! - - - * * * D U D A * * *
        !no puedo colocar direwctamente el nombre d ela funcion (la funcion misma) como variable global y debo pasarselo a otra variable y esa si ponerla como
        !vbariable global por eso hago el cambio de CompH = H y H esta como variable global. Es Asi?
      ! - - - * * * D U D A * * *

      return

    end function CompH

    function compBmat(Nx, Ny, Gp)

      implicit none

      real, dimension(Nne,size(gauss_points) ) :: Nx, Ny
      integer, intent (in) :: Gp

      ! real, dimension(4, 2*Nne)  :: B
      real, dimension(4, 2*Nne)  :: compBmat
      real, dimension(1, Nne)    :: Nxi, Neta

      integer::  i

      ! B = 0
      Nxi  = spread(Nx(:,Gp),dim = 1, ncopies= 1)
      Neta = spread(Ny(:,Gp),dim = 1, ncopies= 1)


      do i=1, Nne
        compBmat(1,2*i-1)= Nxi(1,i)
        compBmat(3,2*i)  = Nxi(1,i)
        compBmat(2,2*i-1)= Neta(1,i)
        compBmat(4,2*i)  = Neta(1,i)
      end do

      ! compBmat = B
      return
      ! - - - * * * D U D A * * * - - -
        !En matlab basta con  Nxi(i) aqui no es posible indicar un vector solo con una dimension?
        !Siempore se debe indicar matriz como un vector fila o vector columna?
      ! - - - * * * D U D A * * * - - -

    end function compBmat

    function AssembleK(K, ke, ndDOF)

      implicit none
      real, dimension(2*n_nodes+n_pnodes, 2*n_nodes+n_pnodes)  :: K !Global Stiffnes matrix
      real, dimension(2*Nne, 2*Nne)                            :: ke
      real, dimension(2*n_nodes+n_pnodes, 2*n_nodes+n_pnodes)  :: AssembleK
      integer :: ndDOF, i, j, row_node, row, col_node, col !nodal Degrees of Freedom

      !K debe llevar inout por que entra como variable (IN) pero en esta funcion se modifica (out)

      do i = 1, Nne
        row_node = node_id_map(i,1)
        row = ndDOF*row_node - (ndDOF-1)

        do j = 1, Nne
          col_node = node_id_map(j,1)
          col = ndDOF*col_node - (ndDOF-1)
          AssembleK(row:row+ndDOF-1, col:col+ndDOF-1) =  K(row:row+ndDOF-1, col:col+ndDOF-1) + &
          ke((i-1)*ndDOF+1:i*ndDOF,(j-1)*ndDOF+1:j*ndDOF)
        enddo

      enddo


      return

    end function AssembleK

    ! subroutine GlobalK( A_K ) !Al tener un solo parametro de salida puedo declararla como funcion

    !   implicit none

    !   !- - - * * * DUDA * * * - - -

    !     !Esto ya estaba como variable global, por que edebo declararlo de nuevo aqui
    !     ! real                      :: materials
    !     ! real, dimension(341,3)    :: nodes
    !     ! integer, dimension(100,9) :: elements
    !     ! real, dimension(341,2)    :: pnodes
    !     ! integer, dimension(100,5) :: pelements

    !     ! ----> Pide declararlo porque lo estoy ponioendo como argumento de entrada en la subrutina, sino lo pongo entonces dentro de la subrutina
    !     ! lo toma de la parte global donde ha sido declarado



    !     ! n_nodes               Ya declarado como variable global
    !     ! n_pnodes              Ya declarado como variable global
    !     ! n_elements            Ya declarado como variable global
    !     ! Nne   Ya declarado como variable global
    !   !- - - * * * * * * * * * - - -

    !   real, allocatable, dimension(:,:)       :: K
    !   real, dimension(2*n_nodes+n_pnodes, 2*n_nodes+n_pnodes),intent(out) :: A_K!Global Stiffnes matrix
    !   ! real, dimension(2*n_nodes+n_pnodes, 2*n_nodes+n_pnodes) :: AssembleK

    !   real, dimension(2*Nne, 2*Nne)           :: ke
    !   real, dimension(dim_prob, dim_prob)     :: Jaco, Jinv
    !   real, dimension(3,3)                    :: cc, C
    !   real, dimension(2*dim_prob, 2*dim_prob) :: Jb
    !   real, dimension(4,2*Nne)                :: B
    !   real, dimension(3,dim_prob*dim_prob)    :: HJ
    !   real, dimension(3,2*Nne)                :: HJB
    !   real, dimension(2*Nne,3)                :: HJB_T
    !   real, dimension(3,4)                    :: H
    !   real                                    :: detJ
    !   integer                                 :: gp, ngp, e

    !   allocate(K(2*n_nodes+n_pnodes, 2*n_nodes+n_pnodes) )

    !   K  = 0.0
    !   cc = reshape([2, 0, 0, 0, 2, 0, 0, 0, 1],[3,3])
    !   C  = materials * cc
    !   H  = CompH()
    !   ngp= size(gauss_points,1) !TMB PODRIA SER VARIABLE PERMANENTE CON SAVE

    !   !elements loop
    !   do e = 1, n_nodes
    !     ke = 0
    !     Jb = 0
    !     call SetElementNodes(e, elements, nodes, element_nodes, node_id_map)
    !     !do-loop: compute element stiffness matrix ke
    !     do gp  = 1, ngp
    !       Jaco = J2D(element_nodes, Nx, Ny, gp)
    !       detJ = m22det(Jaco)
    !       Jinv = inv2x2(Jaco)
    !       Jb   = buildJb (Jinv)
    !       B    = compBmat( Nx, Ny, gp)
    !       HJ   = matmul(H,Jb)
    !       HJB  = matmul(HJ,B)
    !      HJB_T = transpose(HJB)
    !      print*, shape(HJB_T)
    !      !aqui marcaba error por que gauss_weights es un vector columna y debe indicarse con dos indices
    !      ke = 1 !ke + HJB_T * C * HJB * detJ * gauss_weights(gp,1)
    !     end do

    !     A_K = AssembleK(K, ke, 2) ! assemble global K

    !   end do

    ! end subroutine GlobalK










end module library
