      program pgrm_03_01
!
!     This program computes the inverse square-root of a matrix.
!
!     At execution time, the program expects 2 command line arguments: (1) nDim;
!     and (2) the matrix being raised to the (1/2) power.
!
!
      implicit none
      integer,parameter::unitIn=10
      integer::i,iError,nDim,lenSym
      real,dimension(:),allocatable::inputSymMatrix
      real,dimension(:,:),allocatable::inputSqMatrix,invSqrtInputMatrix
      character(len=256)::cmdlineArg
!
!
!     Begin by reading the leading dimension of the matrix and the input file
!     name from the command line. Then, open the file and read the input matrix,
!     inputSymMatrix
!
      call Get_Command_Argument(1,cmdlineArg)
      read(cmdlineArg,'(I)') nDim
      lenSym = (nDim*(nDim+1))/2
      allocate(inputSymMatrix(lenSym),inputSqMatrix(nDim,nDim),  &
        invSqrtInputMatrix(nDim,nDim))
      call Get_Command_Argument(2,cmdlineArg)
      open(Unit=unitIn,File=TRIM(cmdlineArg),Status='OLD',IOStat=iError)
      if(iError.ne.0) then
        write(*,*)' Error opening input file.'
        STOP
      endIf
      do i = 1,lenSym
        read(unitIn,*) inputSymMatrix(i)
      endDo
      close(Unit=unitIn)
!
!     Form the square-root of inputSymMatrix. The result is loaded into square
!     (full storage) matrix invSqrtInputMatrix.
!
      write(*,*)' The matrix loaded (column) upper-triangle packed:'
      call SymmetricPacked2Matrix_UpperPacked(nDim,inputSymMatrix,
        inputSqMatrix)
      write(*,*)' Input Matrix:'
      call Print_Matrix_Full_Real(inputSqMatrix,nDim,nDim)
      Call invsqrt_symmatrix(nDim,inputSqMatrix,invSqrtInputMatrix)
      write(*,*)' Inverse sqrt Matrix:'
      Call Print_Matrix_Full_Real(invSqrtInputMatrix,nDim,nDim)
      write(*,*)' Matrix Product that should be the identity.'
      call Print_Matrix_Full_Real(MatMul(MatMul(invSqrtInputMatrix, &
        invSqrtInputMatrix),inputSqMatrix),nDim,nDim)

!
      end program pgrm_03_01


      Subroutine SymmetricPacked2Matrix_UpperPac(N,ArrayIn,AMatOut)
      Implicit None
      Real,Dimension(N*N)::ArrayIn
      Real,Dimension(N,N)::AMatOut
      Integer::i,j,k,N


      
      Do j=1, N
        Do i=1, j
        AMatOut(i,j)=ArrayIn(k)
        AMatOut(j,i)=AMatOut(i,j)
        k=k+1
        EndDo
        EndDo

      Return
      End Subroutine SymmetricPacked2Matrix_UpperPac
      Subroutine Print_Matrix_Full_Real(AMat,M,N)
      
       implicit none
       integer,intent(in)::M,N
       real,dimension(M,N),intent(in)::AMat
       integer,parameter::IOut=6,NColumns=5
       integer::i,j,IFirst,ILast

      1000 Format(1x,A)
      2000 Format(5x,5(7x,I7))
      2010 Format(1x,I7,5F14.6)
        Do IFirst = 1,N,NColumns-1,N
          ILast = Min(IFIrst+NColumn-1,N)
          write(IOut,2000_ (i,i=IFirst,ILast)
          Do i = 1,M
            write(IOut,2010) i,AMat(i,j),j=IFIrst,ILast_
          endDo
        endDo
        Return
        End Subroutine Print_Matrix_Full_Real

        Subroutine invsqrt_symmatrix(N,inputSymMatrix ,Mat)

        implicit none
        integer,intent(In)::N
        real,dimension((N*(N+1))/2)),intent(In)::inputSymMatrix
        real,dimension(N,N),intent(Out)::Mat

        Integer::i,j,k,IError,s
        Real,Dimension(:,:),Allocatable::EVecs,Temp_Matrix,Evalmat
        Real,Dimension(:),Allocatable::EVals,Temp_Vector,Temp2

        Allocate(EVals(N),EVecs(N,N),Temp_Vector(3*N),Evalmat(N,N))
        Allocate(Temp_Matrix(N,N),Temp2((N*(N+1))/2))
        Do i = 1,(N*(N+1))/2
          Temp2(i)= inputSymMatrix(i)
        EndDo
        Call SSPEV("V","U",N,Temp2,EVals,EVecs,N,  &
        Temp_Vector,IError)        
       
         Do i= 1,N
           Evalmat(i,i) =1/(sqrt(Evals(i)))
         EndDo
        Mat = Matmul(Matmul(EVecs,Evalmat),transpose(EVecs))

        Return
      End Subroutine
