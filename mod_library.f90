module library
  use Isoparametric
  Implicit None

  ! ! ! Aqui se declaran las variables globales ! ! !

  real                      :: materials
  real,    dimension(341,3) :: nodes
  integer, dimension(100,9) :: elements
  integer, dimension(341,2) :: pnodes
  integer, dimension(100,5) :: pelements
  integer, parameter        :: n_nodes = size(nodes,1)
  integer, parameter        :: n_pnodes = 121 !Duda, como lo tomo desde el txt es decir   n_pnodes = maxval(pnodes(:,2))
  integer, parameter        :: n_elements = size(elements,1)
  integer, parameter        :: Nne = size(elements,2)-1     !Number of nodes in the element
  integer, parameter        :: Npne = size(pelements,2)-1   !Number of pnodes in the eleement
  integer, parameter        :: dim_prob = size(nodes,2)-1 !dimension del problema

  double precision, allocatable, dimension(:,:) :: gauss_points, gauss_weights !Verificar si debe ser global---> Si, se usa en la funcion ComputeK

  ! ! ! Fin de variables globales ! ! !


  contains

    subroutine ReadRealFile(UnitNum, FileName, NumRows, NumCols, Real_Array)
      implicit none

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

      integer :: i, j, status
      integer, intent(in)            :: UnitNum, NumRows, NumCols
      character (len=*), intent (in) :: FileName
      integer, dimension (1:NumRows, 1:NumCols), intent (out) :: IntegerArray


      open (unit = UnitNum, file =FileName, status='old', action='read' , iostat = status)

      ! read in values
      read(UnitNum,*) ((IntegerArray(i,j), j=1,NumCols), i=1,NumRows)
      print *, "Status_Int_File  ", status
      print*, "Shape of ",FileName," is ", shape(IntegerArray)
      close (UnitNum)

    end subroutine ReadIntegerFile

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
<<<<<<< HEAD
    
=======

>>>>>>> Se cambiaron nombres de Nx por dN_dxi y Ny por dN_deta en todo los files .f90
    subroutine ReadMixFile(UnitNum, FileName, NumRows, NumCols, Real_Array)
      implicit none

      ! - - - - - - - - - - * * * * * * * * * * - - - - - - - - - -
      ! Rutina que lee un conjunto de datos en el formato indicado
      ! en la etiqueta 22
      !- - - - - - - - - - * * * * * * * * * * - - - - - - - - - -

      integer :: i, j, status, UnitNum, NumRows, NumCols
      character (len=*), intent (in) :: FileName
      real, dimension (1:NumRows, 1:NumCols), intent (out) :: Real_Array


      open (unit = UnitNum, file =FileName, status='old', action='read' , iostat = status)

      ! read in values
      read(UnitNum,22) ((Real_Array(i,j), j=1,NumCols), i=1,NumRows)
      print *, "Status_Mix_File  ", status

      22 format(3F13.10)

      close (UnitNum)

    end subroutine
<<<<<<< HEAD

    subroutine GetQuadGauss(fila, columna, gauss_points, gauss_weights)
      implicit none

      integer,intent(in) :: fila, columna
      double precision, allocatable, dimension(:,:),intent(out) :: gauss_points, gauss_weights
      integer :: i,j,k
      double precision, allocatable, dimension(:,:) :: w1, w2, w, x

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


      ! Aqui puedo poner un ngp = size(gauss_points,1) ! y declarar ngp con SAVE para tener siempre el valor de la variable ngp 

      ! DEALLOCATE(w1, w2, x, w)
      !Esta funcion no afecta al resultado pues se ha liberado la memoria para calcular pesos y puntos
      !de Gauss mas no las variables que contienen pesos y puntos de gauus.
    end subroutine GetQuadGauss

    subroutine CompNDNatPointsQuad8(gauss_points, N, Nx, Ny)
      implicit none

      double precision, dimension(:,:), intent(in) :: gauss_points
      double precision, allocatable, dimension(:,:), intent(out) :: N, Nx, Ny
      double precision, dimension(size(gauss_points,1)) :: xi_vector, eta_vector
      integer, dimension(Nne,dim_prob) :: master_nodes
      double precision    :: xi, eta, mn_xi, mn_eta
      integer :: ngp, i, j, jj, k

      !number of gauss points
      ngp = size(gauss_points,1) ! esta puede quedar como variable global si se usa en alguna otra subrutina
                                ! si solo se usa aqui, entonces variable como local-----> Si se usa en otra rutina, en compK
      allocate( N(Nne,ngp), Nx(Nne,ngp), Ny(Nne,ngp) )

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

    subroutine CompNDNatPointsQuad4(gauss_points, Np)
      implicit none

      double precision, dimension(:,:), intent(in) :: gauss_points
      double precision, allocatable, dimension(:,:), intent(out) :: Np
      double precision, dimension(size(gauss_points,1)) :: xi_vector, eta_vector
      integer, dimension(Npne,dim_prob) :: master_nodes
      double precision    :: xi, eta, mn_xi, mn_eta
      integer :: ngp, i, j

      ngp = size(gauss_points,1) 
      !number of gauss points
      ! esta puede quedar como variable global si se usa en alguna otra subrutina
                                ! si solo se usa aqui, entonces variable como local-----> Si se usa en otra rutina, en compK
      allocate( Np(Npne,ngp))

      Np  = 0.0
      xi_vector  = gauss_points(:,1)     ! xi-coordinate of point j
      eta_vector = gauss_points(:,2)

      !coordinates of the nodes of the master element
      master_nodes = reshape([1, -1, -1, 1, 1, 1, -1, -1], [Npne,dim_prob])
      !NOTA ** Para que el reshape funcione correctamente, o que produzca el par de valores deseado, primero se deben
      !colocar todos los valores en x, luego todos los de y y luego, si hubiera, todos los de z para que al acomodarse salga el par
      !suponiendo un reshape de 3,2 debe acomodarse x1, x2, x3, y1, y2, y3 DUDA *Siempre debe ser es asi*

      do j = 1, ngp
        xi  = xi_vector(j)      ! xi-coordinate of point j
        eta = eta_vector(j)     ! eta-coordinate of point j
        do i = 1, 4
          mn_xi = master_nodes(i,1)
          mn_eta= master_nodes(i,2)
          Np(i,j)=(1.0 + mn_xi*xi)*(1.0 + mn_eta*eta)/4.0 
        end do
      end do

    end subroutine CompNDNatPointsQuad4
                              
=======
                             
>>>>>>> Se cambiaron nombres de Nx por dN_dxi y Ny por dN_deta en todo los files .f90
    subroutine SetElementNodes(elm_num, elements, nodes, element_nodes, node_id_map)
      implicit none

      integer, dimension(100,9),  intent(in)::  elements
      real, dimension(341,3), intent(in)    ::  nodes
      integer,intent(in)                    :: elm_num ! number of element for each elemental integral in do of K global
      real, dimension(Nne,dim_prob), intent(out) :: element_nodes
      integer, dimension(Nne,1), intent(out)     :: node_id_map
      integer                               :: i,j, global_node_id


      element_nodes = 0.0
      node_id_map = 0.0

      do i = 1, Nne
        global_node_id = elements(elm_num,i+1)
        do j=1 ,dim_prob
          element_nodes(i,j) = nodes(global_node_id,j+1)
        end do
        node_id_map(i,1) = global_node_id
      end do

    end subroutine SetElementNodes

    subroutine PreassureElemNods(elm_num, pelements, nodes, pelement_nodes, pnode_id_map)
      implicit none

      integer, dimension(100,5),  intent(in)::  pelements
      real, dimension(341,3), intent(in)    ::  nodes
      integer,intent(in)                    :: elm_num ! number of element for each elemental integral in do of K global
      real, dimension(Npne,dim_prob), intent(out) :: pelement_nodes
      integer, dimension(Npne,1), intent(out)     :: pnode_id_map
      integer                               :: i,j, global_node_id


      pelement_nodes = 0.0
      pnode_id_map = 0.0

      do i = 1, Npne
        global_node_id = pelements(elm_num,i+1)
        do j=1 ,dim_prob
          pelement_nodes(i,j) = nodes(global_node_id,j+1)
        end do
        pnode_id_map(i,1) = global_node_id
      end do

    end subroutine PreassureElemNods

    function J2D( element_nodes, dN_dxi, dN_deta, Gp)
      implicit none

      real, dimension(Nne,dim_prob), intent(in)            :: element_nodes
      double precision, dimension(Nne,size(gauss_points) ), intent(in) :: dN_dxi, dN_deta
      integer, intent(in)                                  :: Gp !esta variable se usara en el lazo principal con el numero de punto de gauss para evaluar las integrales elementales
      double precision, dimension(dim_prob,Nne)                        :: Basis2D
      double precision, dimension(1,Nne)                               :: Nxi, Neta
      double precision, dimension(dim_prob,dim_prob)                   :: J2D


      !con estas instrucciones extraigo la columna de Nx como renglon y lo guardo en Nxi, Gp se
      !ira moviendo conforme la funcion J2D sea llamada en el lazo principal para cada elemento lo mismo para Neta con dN_deta
      Nxi  = spread(dN_dxi(:,Gp),dim = 1, ncopies= 1)
      Neta = spread(dN_deta(:,Gp),dim = 1, ncopies= 1)

      !Las siguientes tres lineas realizan de forma implicita el calculo de las derivadas
      !espaciales es decir dN/dx and dN/dy (eq. 5.114 - 5.117). Las derivadas espaciales 
      !no se calcula explicitamente, en su lugar se usa:
        
      !            d/dy = (d/deta)J^-1      (ver eq. 5.77 y 5.109)

      ! Esta forma de 
      Basis2D(1,:) = Nxi(1,:)
      Basis2D(2,:) = Neta(1,:)
      J2D = matmul(Basis2D,element_nodes) !Aqui se usan directamente las derivadas (eqs 5.114-5.117) de las coordenadas fisicas
                                          !respecto de las coordenadas del master element (isoparametric domain) para llenar la matriz Jacobiano.
      ! De la subroutina Quad4Nodes or Quad8Nodes  ya se tienen las derivadas de las funciones de forma respecto a las coordenadas del master element
      ! es decir dN/dxi and dN/deta contenidas en Basis 2D. Luego, se multiplican por element_nodes para completar el Jacobiano.
      

      ! - - - * * * D U D A * * * - - -
        !Si el nombre de la funcion es el mismo que la variable donde se guarda, entonces no puedo declararla como
        ! variable global, ¿Como debo hacerlo?

        !Si lo dejo como

        !J = Matmul(Basis2D,element_nodes)

        !Me marca un warning
      ! - - - * * * D U D A * * * - - -

      return
    end function J2D

    function inv2x2(A)

      implicit none

      double precision, dimension(dim_prob,dim_prob), intent(in)  :: A
      double precision, dimension(dim_prob,dim_prob)              :: inv2x2

      double precision, parameter :: EPS = 1.0E-10
      double precision :: det
      double precision, dimension(2,2) :: cofactor


      det =   A(1,1)*A(2,2) - A(1,2)*A(2,1)

      if (abs(det) .le. EPS) then
        inv2x2 = 0.0D0
        return
      end IF

      cofactor(1,1) = +A(2,2)
      cofactor(1,2) = -A(2,1)
      cofactor(2,1) = -A(1,2)
      cofactor(2,2) = +A(1,1)

      inv2x2 = transpose(cofactor) / det

      return

    end function inv2x2

    function buildJb(A)
      !Funcion que construye una matriz de 4 x 4 en bloques de 2 para el caso 2D

      implicit none

      double precision, dimension(dim_prob,dim_prob), intent (in)   :: A
      double precision, dimension(2*dim_prob, 2*dim_prob)           :: buildJb

      buildJb(1:2,1:2) = A
      buildJb(3:4,3:4) = A

    end function

    function m22det(A)

      implicit none
      double precision :: m22det
      double precision, dimension(2,2), intent(in)  :: A



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

    function compBmat(dN_dxi, dN_deta, Gp)

      implicit none

      double precision, dimension(Nne,size(gauss_points) ) :: dN_dxi, dN_deta
      integer, intent (in) :: Gp

      ! real, dimension(4, 2*Nne)  :: B
      double precision, dimension(4, 2*Nne)  :: compBmat
      double precision, dimension(1, Nne)    :: Nxi, Neta

      integer::  i

      ! B = 0
      Nxi  = spread(dN_dxi(:,Gp),dim = 1, ncopies= 1)
      Neta = spread(dN_deta(:,Gp),dim = 1, ncopies= 1)


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

    subroutine AssembleK(K, ke, node_id_map, ndDOF)

      implicit none
      real(8), dimension(2*n_nodes+n_pnodes, 2*n_nodes+n_pnodes),intent(in out)  :: K !Global Stiffnes matrix
      real(8), dimension(2*Nne, 2*Nne), intent(in)    :: ke
      integer, dimension(Nne,1), intent(in)           :: node_id_map
      integer, intent(in)                             :: ndDOF 
      integer :: i, j, row_node, row, col_node, col !nodal Degrees of Freedom

      !K debe llevar inout por que entra como variable (IN) pero en esta funcion se modifica (out)

      do i = 1, Nne
        row_node = node_id_map(i,1)
        row = ndDOF*row_node - (ndDOF-1)

        do j = 1, Nne
          col_node = node_id_map(j,1)
          col = ndDOF*col_node - (ndDOF-1)
          K(row:row+ndDOF-1, col:col+ndDOF-1) =  K(row:row+ndDOF-1, col:col+ndDOF-1) + &
          ke((i-1)*ndDOF+1:i*ndDOF,(j-1)*ndDOF+1:j*ndDOF)
        enddo

      enddo


      return

    end subroutine AssembleK

    subroutine GlobalK( A_K, dN_dxi, dN_deta) !Al tener un solo parametro de salida puedo declararla como funcion

      implicit none

      !- - - * * * DUDA * * * - - -

        !Esto ya estaba como variable global, por que edebo declararlo de nuevo aqui
        ! real                      :: materials
        ! real, dimension(341,3)    :: nodes
        ! integer, dimension(100,9) :: elements
        ! real, dimension(341,2)    :: pnodes
        ! integer, dimension(100,5) :: pelements

        ! ----> Pide declararlo porque lo estoy ponioendo como argumento de entrada en la subrutina, sino lo pongo entonces dentro de la subrutina
        ! lo toma de la parte global donde ha sido declarado



        ! n_nodes               Ya declarado como variable global
        ! n_pnodes              Ya declarado como variable global
        ! n_elements            Ya declarado como variable global
        ! Nne   Ya declarado como variable global
      !- - - * * * * * * * * * - - -

      double precision, dimension(2*n_nodes+n_pnodes, 2*n_nodes+n_pnodes),intent(out) :: A_K  !Global Stiffnes matrix
      double precision, dimension(Nne,size(gauss_points,1)), intent(in)               :: dN_dxi, dN_deta
      double precision, allocatable, dimension(:,:)       :: Np, dd, ddd
      double precision, dimension(2*Nne, 2*Nne)        :: ke
      double precision, dimension(dim_prob, dim_prob)     :: Jaco, Jinv
      double precision                                    :: detJ
      real,  dimension(3,3)                               :: cc, C
      double precision, dimension(2*dim_prob, 2*dim_prob) :: Jb
      double precision, dimension(4,2*Nne)                :: B
      double precision, dimension(3,dim_prob*dim_prob)    :: HJ
      double precision, dimension(3,2*Nne)                :: HJB
      double precision, dimension(2*Nne,3)                :: HJB_T
      real, dimension(3,4)                    :: H
      double precision, dimension(16,3)                   :: part1
      double precision, dimension(16,16)                  :: part2
      double precision, dimension(16,16)                  :: part3
      real(8), allocatable, dimension(:,:)    :: K12, K12_T!Lo puse allocatable por que marca error en la memoria 
      ! Array 'k12' at (1) is larger than limit set by '-fmax-stack-var-size=', moved from stack to static storage. This makes the procedure unsafe when called recursively, 
      !or concurrently from multiple threads. Consider using '-frecursive', or increase the '-fmax-stack-var-size=' limit, or change the code to use an ALLOCATABLE array. [-Wsurprising]
      
      double precision, dimension(16,4)                   :: kep
      double precision, dimension(8,2)                    :: part4
      double precision, dimension(2,1)                    :: part5
      double precision, dimension(2,1)                    :: A
      double precision, dimension(4,1)                    :: part6
      double precision, dimension(1,4)                    :: part7
      double precision, dimension(16,4)                   :: part8
      double precision, dimension(4*Npne,1)               :: dN
      integer, dimension(Nne,1)               :: node_id_map
      real, dimension(Nne,dim_prob)           :: element_nodes
      integer, dimension(Npne,1)              :: pnode_id_map
      real, dimension(Npne,dim_prob)          :: pelement_nodes

      integer                                 :: gp, ngp, e, i,j, row_node, row 
      integer                                 :: col_node, pnode_id, col!, mrow, ncol
      

      ngp = size(gauss_points,1) !TMB PODRIA SER VARIABLE PERMANENTE CON SAVE

      A_K  = 0.0
      cc = reshape([2, 0, 0, 0, 2, 0, 0, 0, 1],[3,3])
      C  = materials * cc
      H  = CompH()
      
      !elements loop for K1-1 block Global K
      do e = 1, n_elements
        ke = 0
        Jb = 0
        call SetElementNodes(e, elements, nodes, element_nodes, node_id_map)
        !do-loop: compute element stiffness matrix ke
        do gp  = 1, ngp
          Jaco = J2D(element_nodes, dN_dxi, dN_deta, gp)
          detJ = m22det(Jaco)
          Jinv = inv2x2(Jaco)
          Jb   = buildJb (Jinv)
          B    = compBmat( dN_dxi, dN_deta, gp)
          HJ   = matmul(H,Jb)
          HJB  = matmul(HJ,B)
         HJB_T = transpose(HJB)
         part1 = matmul(HJB_T,C)
         part2 = matmul(part1,HJB)
         part3 = part2 * detJ
         !aqui marcaba error por que gauss_weights es un vector columna y debe indicarse con dos indices
         ke = ke + part3 * gauss_weights(gp,1)  !
        end do

        call AssembleK(A_K, ke, node_id_map, 2) ! assemble global K
      
      end do

      !Setup for K1-2 block
      ! Npne is declared at the top as parameter
      allocate (K12(n_nodes*2,n_pnodes),K12_T(n_pnodes,n_nodes*2))
      K12 = 0
      
      call Quad4Nodes(gauss_points, ngp, Np, dd,ddd)

      do e = 1, n_elements
        kep = 0.0
        call SetElementNodes(e, elements, nodes, element_nodes, node_id_map)
        call PreassureElemNods(e, pelements, nodes, pelement_nodes, pnode_id_map) !--Arreglar esto para que sea con p en todos los arguments
        ! for-loop: compute element stiffness matrix ke     
        do gp  = 1, ngp
          Jaco = J2D(element_nodes, dN_dxi, dN_deta, gp)
          detJ = m22det(Jaco)
          Jinv = inv2x2(Jaco)
          dN = 0.0
          do j = 1, Nne
            part4(j,:) = [ dN_dxi(j,gp), dN_deta(j,gp) ]  
            part5 = reshape([part4(j,:)],[2,1])           !--Revisar por que en la linea 514 y 515 si se puede hacer
            A =  matmul(Jinv,part5)           !Tuve que separar todas las multiplicaciones para que funcione 
                                              ! pues el resultado de matmul debe guardarse en otra variable, A sino marca error
            dN(2*j-1:2*j ,1)= A(:,1)          !quiza si necesite dividri la operacion de este matmul
          end do
          part6(:,1) = Np(:,gp)
          part7 = transpose(part6)
          part8 = matmul(dn,part7)
          kep = kep + part8 * (detJ*gauss_weights(gp,1)) !----> Verificar si hace falta la parte9 separando part8 * detJ
        end do  
      
        ! for-loop: assemble ke into global Kp
        do i = 1, Nne
          
          row_node = node_id_map(i,1)
          row = 2*row_node - 1
          do j = 1, Npne
            col_node = pnode_id_map(j,1)
            pnode_id = pnodes(col_node,2)
            col = pnode_id
            K12(row:row+1, col) = K12(row:row+1, col) + kep(2*i-1:i*2, j)
          enddo
        
        enddo 

        !Desde aqui
          ! print*, 'element number', e

          ! if ( e <= 15 ) then
          !   do i =1,10
          !     print*, K12(i,e)
          !   end do    
          ! else if (e >= 16 .or. e <= 35 )then
          !     do i = 10,25
          !     print*,K12(i,e)
          !     end do
          ! else if ( e >= 36 .or. e<= 65 )then
          !   do i = 55,85
          !     print*,K12(i,e)
          !     end do
          !   end if
                
          !   print*, ' '
          !   print*, ' '
          !   print*, 'pelements_nodes'
          !   print*, ' '
          !   do i = 1,Npne
          !     print*, pelement_nodes(i,:)
          !   end do
          !   print*, ' '
          !   print*, 'pnode_id_map'
          !   print*, ' '
          !   do i = 1,Npne
          !     print*, pnode_id_map(i,:)
          !   end do
        !Y hasta aqui para comprobar la matriz K12
      end do

      !========== Filling the symetric (upper and lower) part of K ==========
      !========== Upper
      do i = 1, 2*n_nodes
        do j = 2*n_nodes+1, (2*n_nodes+n_pnodes)
          A_K(i, j) = -K12(i,j-682)
        end do
      end do
      !========== Lower
      K12_T = transpose(-K12)
      do i = 2*n_nodes+1, (2*n_nodes+n_pnodes)
        do j = 1, 2*n_nodes 
          A_K(i, j) = K12_T(i-682,j)
        end do
      end do

      ! !========== Escribe en archivo la parte simetrica de la matriz global 
      ! mrow = 2*n_nodes
      ! ncol = n_pnodes
      ! open(unit=1, file='K12_T.dat', ACTION="write", STATUS="replace")
      ! do i=1,ncol
      !   write(1, '(1000F14.7)')( K12_T(i,j) ,j=1,mrow)
      ! end do
      ! close(1)
      ! !========== Escribir en archivo la matriz global
      ! mrow = 2*n_nodes+n_pnodes 
      ! ncol = 2*n_nodes+n_pnodes
      ! open(unit=2, file='completeA_K.dat', ACTION="write", STATUS="replace")
      ! do i=1,mrow
      !   write(2, '(1000F14.7)')( A_K(i,j) ,j=1,ncol)
      ! end do
      ! close(2)
      
    end subroutine GlobalK

    subroutine SetBounCond( NoBV, NoBVcol )
      !========================================================================
      !Esta subroutina revisa todos los nodos de la malla y define el tipo de
      !nodo en la frontera. Abre un archivo en donde comenzara a escribir, 
      ! en la primer columna: el numero de nodo. 
      ! La segunda columna tendra el tipo de nodo
      ! 1 = ux (componente x de la velocidad) 
      ! 2 = uy (componente y de la velocidad) 
      ! 3 = para la presion 
      !La tercera columna asigna el valor correspondiente de la condicion de forntera
      !=========================================================================
      implicit none

      ! real, dimension(341,3)    :: nodes
      ! integer, parameter        :: n_nodes = size(nodes,1)
      
      integer, intent(out) :: NoBV, NoBVcol
      integer :: ierror, a ,b, c, i
      real    :: x, y

      ! call ReadRealFile(10,"nodes.dat", 341,3, nodes)
      ! inicializamos los contadores
      ! Los contadores son para que cada vez
      ! que un if se cumpla, se sume el numero
      ! equivalente a los renglones escritos en
      ! el archivo de texto que se esta creando
      ! y asi se tenga el numero total de nodos en
      ! las fronteras
      a = 0
      b = 0
      c = 0 
      
      open(unit=100, file='Fbcsvp.dat',Status= 'replace', action= 'write',iostat=ierror)
      NoBVcol = size(nodes,2)     
     
      do i =1, n_nodes
        x=nodes(i,2)
        y=nodes(i,3)
        if(y.eq.1.0) then
         write(100,50) real(i), 1.0, 1.0
         write(100,50) real(i), 2.0, 0.0
         a=a+2
         if(x.eq.0.5)then
            write(100,50) real(i),3.0,0.0
            b=b+1
         end if

        else if (x.eq.0.0 .or. y.eq.0.0 .or. x.eq.1.0)then
         write(100,50) real(i), 1.0, 0.0
         write(100,50) real(i), 2.0, 0.0
         c=c+2

        end if
        NoBV = a+b+c
      end do
    
      close(100)
    
      50 format(3F15.10)

   
    end subroutine SetBounCond  

    subroutine ApplyBoundCond( NoBV, Fbcsvp, A_K, Sv )
      ! - - - - - - - - - - * * * * * * * * * * - - - - - - - 
      ! Set velocity (u) and preasure (p) boundary condition by penalty method
      ! - - - - - - - - - - * * * * * * * * * * - - - - - - - - - -
      implicit none
      
      real, dimension(NoBV,3), intent(in) :: Fbcsvp
      double precision :: param, coeff
      double precision, dimension(2*n_nodes+n_pnodes, 2*n_nodes+n_pnodes),intent(in out) :: A_K  !Global Stiffnes matrix
      double precision, dimension(2*n_nodes+n_pnodes, 1), intent(in out) :: Sv
      integer :: preasure_row, NoBV, i, component, node_id, pnode_id

      !Esencialmente la siguiente instruccion hace: A_K(1*2-1,:) = A_K(1,:) Es decir, obtene el valor maximo de
      !la primera fila de la matriz global K (A_K). No le veo el caso pero lo dejamos asi.
      param = maxval(A_K(int(Fbcsvp(1,1))*2-1,:))
      coeff = abs(param) * 1.0E7
      
      preasure_row = 2*n_nodes

      do i =1, NoBV
        node_id = int(Fbcsvp(i,1)) !se pone este int() pq la 1a y 2a col de Fbcsvp esta leida como integer pero 
        component = int(Fbcsvp(i,2))!la matriz completa esta declarada como real en esta subroutina y en el main.
        if ( component .le. 2 ) then
          A_K(2*node_id-2+component, 2*node_id-2 +component) = coeff
          Sv(2*node_id-2+component, 1) = Fbcsvp(i,3)*coeff 
        else                                                     
          pnode_id = pnodes(node_id,2)
          A_K(preasure_row+pnode_id, preasure_row + pnode_id) = coeff
          Sv(preasure_row+pnode_id,1) = Fbcsvp(i,3)*coeff 
        end if
      end do

    end subroutine ApplyBoundCond
  


  !Fin de contains






end module library
