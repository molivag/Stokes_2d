module Parameters
  implicit none
  
  integer, parameter :: DimPr     = 2     !Dimension del problema 
  integer, parameter :: Nelem     = 100   !Number of elements
  integer, parameter :: n_nodes   = 341   !Total number of velocity nodes
  integer, parameter :: n_pnodes  = 121   !Total number of preasure nodes MAXVAL(pnodes,2)
  integer, parameter :: nUne      = 8     !Number of velocity nodes in the element
  integer, parameter :: nPne      = 4     !Number of preasure nodes in the element
  integer, parameter :: Dof       = 3     !Degrees of fredoom: 2 for velocity + 1 one for preasure
  integer, parameter :: ngp       = 2     !Number of Gauss points for quadrature  

  integer, dimension(Nelem, nUne + 1)   :: elements
  integer, dimension(Nelem, nPne + 1)   :: pelements  
  real,    dimension(n_nodes, DimPr + 1):: nodes
  integer, dimension(n_nodes, 2)        :: pnodes
  real                                  :: materials

  ! integer, dimension(100,9) :: elements
  ! integer, dimension(100,5) :: pelements  
  ! real,    dimension(341,3) :: nodes
  ! integer, dimension(341,2) :: pnodes
  ! real                      :: materials


  ! integer, parameter        :: Nelem = size(elements,1)
  ! integer, parameter        :: n_nodes = size(nodes,1)
  ! integer, parameter        :: n_pnodes = 121 !Duda, como lo tomo desde el txt es decir   n_pnodes = maxval(pnodes(:,2))
  ! integer, parameter        :: nUne = size(elements,2)-1    !Number of velocity nodes in the element
  ! integer, parameter        :: nPne = size(pelements,2)-1   !Number of preasure nodes in the eleement
  ! integer, parameter        :: DimPr = size(nodes,2)-1   !Dimension del problema

  != = = = = =  About Isoparametric Mapping  = = = = =  
  
  double precision, allocatable, dimension(:,:) :: gauss_points, gauss_weights !Verificar si debe ser global---> Si, se usa en la funcion ComputeK
  character(len=6)  :: ElemType !could be: Trian, Quad
  integer           :: totGp

  
end module Parameters

