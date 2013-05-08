!> @ingroup Library
!> @{
!> @defgroup Lib_VTK_IOLibrary Lib_VTK_IO
!> @}

!> @ingroup Interface
!> @{
!> @defgroup Lib_VTK_IOInterface Lib_VTK_IO
!> @}

!> @ingroup PrivateVarPar
!> @{
!> @defgroup Lib_VTK_IOPrivateVarPar Lib_VTK_IO
!> @}

!> @ingroup PublicProcedure
!> @{
!> @defgroup Lib_VTK_IOPublicProcedure Lib_VTK_IO
!> @}

!> @ingroup PrivateProcedure
!> @{
!> @defgroup Lib_VTK_IOPrivateProcedure Lib_VTK_IO
!> @}

!> @brief   This is a library of functions for Input and Output pure Fortran data in VTK format.
!> @details It is useful for Paraview visualization tool. Even though there are many wrappers/porting of the VTK source
!>          code (C++ code), there is not a Fortran one. This library is not a porting or a wrapper of the VTK code,
!>          but it only an exporter/importer of the VTK data format written in pure Fortran language (standard Fortran 2003)
!>          that can be used by Fortran coders (yes, there are still a lot of these brave coders...) without mixing Fortran with
!>          C++ language. Fortran is still the best language for high performance computing for scientific purpose, like CFD
!>          computing. It is necessary a tool to deal with VTK standard directly by Fortran code. The library was made to fill
!>          this empty: it is a simple Fortran module able to export native Fortran data into VTK data format and to import VTK
!>          data into a Fortran code, both in ascii and binary file format.
!>
!>          The library provides an automatic way to deal with VTK data format: all the formatting processes is nested into the
!>          library and users communicate with it by a simple API passing only native Fortran data (native Fortran scalar, vector
!>          and matrix).
!>
!>          The library is still in developing and testing, this is first usable release, but there are not all the features of
!>          the stable release (the importer is totally absent and the exporter is not complete). Surely there are a lot of bugs
!>          and the programming style is not the best, but the exporters are far-complete.
!>
!>          The supported VTK features are:
!>          - Exporters:
!>            - Legacy standard:
!>              - Structured Points;
!>              - Structured Grid;
!>              - Unstructured Grid;
!>              - Polydata (\b missing);
!>              - Rectilinear Grid;
!>              - Field (\b missing);
!>            - XML standard:
!>              - serial dataset:
!>                - Image Data (\b missing);
!>                - Polydata (\b missing);
!>                - Rectilinear Grid;
!>                - Structured Grid;
!>                - Unstructured Grid;
!>              - parallel (partitioned) dataset:
!>                - Image Data (\b missing);
!>                - Polydata (\b missing);
!>                - Rectilinear Grid;
!>                - Structured Grid;
!>                - Unstructured Grid;
!>              - composite dataset:
!>                - vtkMultiBlockDataSet;
!>          - Importers are \b missing.
!>
!>          @libvtk can handle multiple concurrent files and it is \b thread/processor-safe (meaning that can be safely used into
!>          parallel frameworks as OpenMP or MPI).
!>
!>          The library is an open source project, it is distributed under the GPL v3. Anyone is interest to use, to develop or
!>          to contribute to @libvtk is welcome.
!>
!>          It can be found at: https://github.com/szaghi/Lib_VTK_IO
!>
!> @par VTK_Standard
!>      VTK, Visualization Toolkit, is an open source software that provides a powerful framework for the computer graphic, for
!>      the images processing and for 3D rendering. It is widely used in the world and so it has a very large community of users,
!>      besides the Kitware (The Kitware homepage can be found here: http://public.kitware.com) company provides professional
!>      support. The toolkit is written in C++ and a lot of porting/wrappers for Tcl/Tk, Java and Python are provided, unlucky
!>      there aren't wrappers for Fortran.
!>
!>      Because of its good features the VTK toolkit has been used to develop a large set of open source programs. For my work
!>      the most important family of programs is the scientific visualization programs. A lot of high-quality scientific
!>      visualization tool are available on the web but for me the best is ParaView: I think that it is one of the best
!>      scientific visualization program in the world and it is open source! Paraview is based on VTK.
!> @par Paraview
!>      ParaView (The ParaView homepage can be found here: http://www.paraview.org) is an open source software voted to scientific
!>      visualization and able to use the power of parallel architectures. It has an architecture client-server in order to make
!>      easy the remote visualization of very large set of data. Because it is based on VTK it inherits all VTK features. ParaView
!>      is very useful for Computational Fluid Dynamics visualizations because it provides powerful post-processing tools, it
!>      provides a very large set of importers for the most used format like Plot3D and HDF (the list is very large). It is easy to
!>      extend ParaView because it supports all the scripting language supported by VTK.
!> @note All the @libvtk functions are <b>I4P integer functions</b>: the returned integer output is 0 if the function calling has
!> been completed right while it is >0  if some errors occur (the error handling is only at its embryonic phase). Therefore the
!> functions calling must be done in the following way: \n
!> @code
!> ...
!> integer(I4P):: E_IO
!> ...
!> E_IO = VTK_INI(....
!> ... @endcode
!> @libvtk programming style is based on two main principles: <em>portable kind-precision</em> of reals and integers
!> variables and <em>dynamic dispatching</em>. Using <em>dynamic dispatching</em> @libvtk has a simple API. The user calls
!> a generic procedure (VTK_INI, VTK_GEO,...) and the library, depending on the type and number of the inputs passed, calls the
!> correct internal function (i.e. VTK_GEO for R8P real type if the input passed is R8P real type). By this interface only few
!> functions are used without the necessity of calling a different function for every different inputs type.
!> Dynamic dispatching is based on the internal kind-precision selecting convention: Fortran 90/95 standard has introduced some
!> useful functions to achieve the portability of reals and integers precision and @libvtk uses these functions to define portable
!> kind-precision; to this aim @libvtk uses IR_Precision module.
!> @author    Stefano Zaghi
!> @version   1.1
!> @date      2013-04-26
!> @par News
!>      - Correct bug affecting binary output;
!>      - implement concurrent multiple files IO capability;
!>      - implement FieldData tag for XML files, useful for tagging dataset with global auxiliary data, e.g. time, time step, ecc;
!>      - implement Parallel (Partitioned) XML files support (.pvtu,.pvts,.pvtr);
!>      - implement Driver testing program for providing practical examples of @libvtk usage;
!>      - added support for parallel framework, namely OpenMP (thread-safe) and MPI (process-safe).
!> @copyright GNU Public License version 3.
!> @note The supported compilers are GNU gfortran 4.7.x (or higher) and Intel Fortran 12.x (or higher). @libvtk needs a modern
!> compiler providing support for some Fortran standard 2003 features.
!> @todo \b CompleteExporter: Complete the exporters
!> @todo \b CompleteImporter: Complete the importers
!> @todo \b DocExamples: Complete the documentation of examples
!> @todo \b g95_test: Test g95 compiler
!> @bug <b>Array-Reshape</b>: \n Fortran allows automatic reshape of arrays, e.g. 2D array can be automatically (in the
!>                            function calling) transformed  to a 1D array with the same number of element of 2D array. The use of
!>                            dynamic dispatching for @libvtk functions by means of generic interfaces had disable this feature:
!>                            dynamic dispatching use the array-shape information to detect, at compile-time,
!>                            the correct function to be called inside the generic interface functions. Thus automatic reshaping
!>                            of arrays at calling function phase is not allowed. \n
!>                            Instead an explicit reshape can be used by means of the Fortran built-in function \em reshape.
!>                            As an example considering a call to the generic function \em VTK_VAR_XML an explicit array reshape
!>                            could be: \n \n
!>                            E_IO = VTK_VAR_XML(NC_NN=nn,varname='u',var=\b reshape(u(ni1:ni2,nj1:nj2,nk1:nk2),(/nn/))) \n \n
!>                            where built in function \em reshape has explicitly being used in the calling to VTK_VAR_XML.
!> @bug <b>XML-Efficiency</b>: \n This is not properly a bug. There is an inefficiency when saving XML binary file. To write XML
!>                             binary @libvtk uses a temporary scratch file to save binary data while saving all formatting data to
!>                             the final XML file. Only when all XML formatting data have been written the scratch file is rewind
!>                             and the binary data is saved in the final tag of XML file as \b raw data. This approach is not
!>                             efficient.
!> @ingroup Lib_VTK_IOLibrary
module Lib_VTK_IO
!-----------------------------------------------------------------------------------------------------------------------------------
USE IR_Precision                                                                ! Integers and reals precision definition.
USE, intrinsic:: ISO_FORTRAN_ENV, only: stdout=>OUTPUT_UNIT, stderr=>ERROR_UNIT ! Standard output/error logical units.
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
save
! functions for VTK XML
public:: VTK_INI_XML
public:: VTK_FLD_XML
public:: VTK_GEO_XML
public:: VTK_CON_XML
public:: VTK_DAT_XML
public:: VTK_VAR_XML
public:: VTK_END_XML
! functions for VTM XML
public:: VTM_INI_XML
public:: VTM_BLK_XML
public:: VTM_WRF_XML
public:: VTM_END_XML
! functions for PVTK XML
public:: PVTK_INI_XML
public:: PVTK_GEO_XML
public:: PVTK_DAT_XML
public:: PVTK_VAR_XML
public:: PVTK_END_XML
! functions for VTK LEGACY
public:: VTK_INI
public:: VTK_GEO
public:: VTK_CON
public:: VTK_DAT
public:: VTK_VAR
public:: VTK_END
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @brief Function for saving mesh with different topologies in VTK-legacy standard.
!> VTK_GEO is an interface to 8 different functions, there are 2 functions for each of 4 different topologies actually supported:
!> one function for mesh coordinates with R8P precision and one for mesh coordinates with R4P precision.
!> @remark This function must be called after VTK_INI. It saves the mesh geometry. The inputs that must be passed change depending
!> on the topologies chosen. Not all VTK topologies have been implemented (\em polydata topologies are absent).
!> @note Examples of usage are: \n
!> \b Structured points calling: \n
!> @code ...
!> integer(I4P):: Nx,Ny,Nz
!> real(I8P)::    X0,Y0,Z0,Dx,Dy,Dz
!> ...
!> E_IO=VTK_GEO(Nx,Ny,Nz,X0,Y0,Z0,Dx,Dy,Dz)
!> ... @endcode
!> \b Structured grid calling: \n
!> @code ...
!> integer(I4P):: Nx,Ny,Nz,Nnodes
!> real(R8P)::    X(1:Nnodes),Y(1:Nnodes),Z(1:Nnodes)
!> ...
!> E_IO=VTK_GEO(Nx,Ny,Nz,Nnodes,X,Y,Z)
!> ... @endcode
!> \b Rectilinear grid calling: \n
!> @code ...
!> integer(I4P):: Nx,Ny,Nz
!> real(R8P)::    X(1:Nx),Y(1:Ny),Z(1:Nz)
!> ...
!> E_IO=VTK_GEO(Nx,Ny,Nz,X,Y,Z)
!> ... @endcode
!> \b Unstructured grid calling: \n
!> @code ...
!> integer(I4P):: NN
!> real(R4P)::    X(1:NN),Y(1:NN),Z(1:NN)
!> ...
!> E_IO=VTK_GEO(NN,X,Y,Z)
!> ... @endcode
!> @ingroup Lib_VTK_IOInterface
interface VTK_GEO
  module procedure VTK_GEO_UNST_R8, & ! real(R8P) UNSTRUCTURED_GRID
                   VTK_GEO_UNST_R4, & ! real(R4P) UNSTRUCTURED_GRID
                   VTK_GEO_STRP_R8, & ! real(R8P) STRUCTURED_POINTS
                   VTK_GEO_STRP_R4, & ! real(R4P) STRUCTURED_POINTS
                   VTK_GEO_STRG_R8, & ! real(R8P) STRUCTURED_GRID
                   VTK_GEO_STRG_R4, & ! real(R4P) STRUCTURED_GRID
                   VTK_GEO_RECT_R8, & ! real(R8P) RECTILINEAR_GRID
                   VTK_GEO_RECT_R4    ! real(R4P) RECTILINEAR_GRID
endinterface
!> @brief Function for saving data variable(s) in VTK-legacy standard.
!> VTK_VAR is an interface to 8 different functions, there are 3 functions for scalar variables, 3 functions for vectorial
!> variables and 2 functions texture variables: scalar and vectorial data can be R8P, R4P and I4P data while texture variables can
!> be only R8P or R4P.
!> This function saves the data variables related to geometric mesh.
!> @remark The inputs that must be passed change depending on the data
!> variables type.
!> @note Examples of usage are: \n
!> \b Scalar data calling: \n
!> @code ...
!> integer(I4P):: NN
!> real(R4P)::    var(1:NN)
!> ...
!> E_IO=VTK_VAR(NN,'Sca',var)
!> ... @endcode
!> \b Vectorial data calling: \n
!> @code ...
!> integer(I4P):: NN
!> real(R4P)::    varX(1:NN),varY(1:NN),varZ(1:NN)
!> ...
!> E_IO=VTK_VAR('vect',NN,'Vec',varX,varY,varZ)
!> ... @endcode
!> @ingroup Lib_VTK_IOInterface
interface VTK_VAR
  module procedure VTK_VAR_SCAL_R8, & ! real(R8P)    scalar
                   VTK_VAR_SCAL_R4, & ! real(R4P)    scalar
                   VTK_VAR_SCAL_I4, & ! integer(I4P) scalar
                   VTK_VAR_VECT_R8, & ! real(R8P)    vectorial
                   VTK_VAR_VECT_R4, & ! real(R4P)    vectorial
                   VTK_VAR_VECT_I4, & ! integer(I4P) vectorial
                   VTK_VAR_TEXT_R8, & ! real(R8P)    vectorial (texture)
                   VTK_VAR_TEXT_R4    ! real(R4P)    vectorial (texture)
endinterface
!> @brief Function for saving field data (global auxiliary data, eg time, step number, dataset name, etc).
!> VTK_FLD_XML is an interface to 7 different functions, there are 2 functions for real field data, 4 functions for integer one
!> and one function for open and close field data tag.
!> @remark VTK_FLD_XML must be called after VTK_INI_XML and befor VTK_GEO_XML. It must always called three times at least: 1) for
!> opening the FieldData tag, 2) for saving at least one FieldData entry and 3) for closing the FieldData tag.
!> Examples of usage are: \n
!> \b saving the time and step cicle counter of current dataset: \n
!> @code ...
!> real(R8P)::    time
!> integer(I4P):: step
!> ...
!> E_IO=VTK_FLD_XML(fld_action='open')
!> E_IO=VTK_FLD_XML(fld=time,fname='TIME')
!> E_IO=VTK_FLD_XML(fld=step,fname='CYCLE')
!> E_IO=VTK_FLD_XML(fld_action='close')
!> ... @endcode
!> @ingroup Lib_VTK_IOInterface
interface VTK_FLD_XML
  module procedure VTK_FLD_XML_OC, & ! open/close field data tag
                   VTK_FLD_XML_R8, & ! real(R8P)    scalar
                   VTK_FLD_XML_R4, & ! real(R4P)    scalar
                   VTK_FLD_XML_I8, & ! integer(I8P) scalar
                   VTK_FLD_XML_I4, & ! integer(I4P) scalar
                   VTK_FLD_XML_I2, & ! integer(I2P) scalar
                   VTK_FLD_XML_I1    ! integer(I1P) scalar
endinterface
!> @brief Function for saving mesh with different topologies in VTK-XML standard.
!> VTK_GEO_XML is an interface to 7 different functions, there are 2 functions for each of 3 topologies supported and a function
!> for closing XML pieces: one function for mesh coordinates with R8P precision and one for mesh coordinates with R4P precision.
!> @remark VTK_GEO_XML must be called after VTK_INI_XML. It saves the mesh geometry. The inputs that must be passed
!> change depending on the topologies chosen. Not all VTK topologies have been implemented (\em polydata topologies are absent).
!> @note The XML standard is more powerful than legacy. XML file can contain more than 1 mesh with its
!> associated variables. Thus there is the necessity to close each \em pieces that compose the data-set saved in the
!> XML file. The VTK_GEO_XML called in the <em>close piece</em> format is used just to close the
!> current piece before saving another piece or closing the file. \n
!> Examples of usage are: \n
!> \b Structured grid calling: \n
!> @code ...
!> integer(I4P):: nx1,nx2,ny1,ny2,nz1,nz2,NN
!> real(R8P)::    X(1:NN),Y(1:NN),Z(1:NN)
!> ...
!> E_IO=VTK_GEO_XML(nx1,nx2,ny1,ny2,nz1,nz2,Nn,X,Y,Z)
!> ... @endcode
!> \b Rectilinear grid calling: \n
!> @code ...
!> integer(I4P):: nx1,nx2,ny1,ny2,nz1,nz2
!> real(R8P)::    X(nx1:nx2),Y(ny1:ny2),Z(nz1:nz2)
!> ...
!> E_IO=VTK_GEO_XML(nx1,nx2,ny1,ny2,nz1,nz2,X,Y,Z)
!> ... @endcode
!> \b Unstructured grid calling: \n
!> @code ...
!> integer(I4P):: Nn,Nc
!> real(R8P)::    X(1:Nn),Y(1:Nn),Z(1:Nn)
!> ...
!> E_IO=VTK_GEO_XML(Nn,Nc,X,Y,Z)
!> ... @endcode
!> \b Closing piece calling: 1n
!> @code ...
!> E_IO=VTK_GEO_XML()
!> ... @endcode
!> @ingroup Lib_VTK_IOInterface
interface VTK_GEO_XML
  module procedure VTK_GEO_XML_STRG_R4, & ! real(R4P) StructuredGrid
                   VTK_GEO_XML_STRG_R8, & ! real(R8P) StructuredGrid
                   VTK_GEO_XML_RECT_R8, & ! real(R8P) RectilinearGrid
                   VTK_GEO_XML_RECT_R4, & ! real(R4P) RectilinearGrid
                   VTK_GEO_XML_UNST_R8, & ! real(R8P) UnstructuredGrid
                   VTK_GEO_XML_UNST_R4, & ! real(R4P) UnstructuredGrid
                   VTK_GEO_XML_CLOSEP     ! closing tag "Piece" function
endinterface
!> @brief Function for saving data variable(s) in VTK-XML standard.
!> VTK_VAR_XML is an interface to 18 different functions, there are 6 functions for scalar variables, 6 functions for vectorial
!> variables and 6 functions for list variables: for all of 3 types of data the precision can be R8P, R4P, I8P, I4P, I2P and I1P.
!> This function saves the data variables related to geometric mesh.
!> @remark The inputs that must be passed change depending on the data variables type.
!> @note Examples of usage are: \n
!> \b Scalar data calling: \n
!> @code ...
!> integer(I4P):: NN
!> real(R8P)::    var(1:NN)
!> ...
!> E_IO=VTK_VAR_XML(NN,'Sca',var)
!> ... @endcode
!> \b Vectorial data calling: \n
!> @code ...
!> integer(I4P):: NN
!> real(R8P)::    varX(1:NN),varY(1:NN),varZ(1:NN),
!> ...
!> E_IO=VTK_VAR_XML(NN,'Vec',varX,varY,varZ)
!> ... @endcode
!> @ingroup Lib_VTK_IOInterface
interface VTK_VAR_XML
  module procedure VTK_VAR_XML_SCAL_R8, & ! real(R8P)    scalar
                   VTK_VAR_XML_SCAL_R4, & ! real(R4P)    scalar
                   VTK_VAR_XML_SCAL_I8, & ! integer(I8P) scalar
                   VTK_VAR_XML_SCAL_I4, & ! integer(I4P) scalar
                   VTK_VAR_XML_SCAL_I2, & ! integer(I2P) scalar
                   VTK_VAR_XML_SCAL_I1, & ! integer(I1P) scalar
                   VTK_VAR_XML_VECT_R8, & ! real(R8P)    vectorial
                   VTK_VAR_XML_VECT_R4, & ! real(R4P)    vectorial
                   VTK_VAR_XML_VECT_I8, & ! integer(I4P) vectorial
                   VTK_VAR_XML_VECT_I4, & ! integer(I4P) vectorial
                   VTK_VAR_XML_VECT_I2, & ! integer(I4P) vectorial
                   VTK_VAR_XML_VECT_I1, & ! integer(I4P) vectorial
                   VTK_VAR_XML_LIST_R8, & ! real(R8P)    list
                   VTK_VAR_XML_LIST_R4, & ! real(R4P)    list
                   VTK_VAR_XML_LIST_I8, & ! integer(I4P) list
                   VTK_VAR_XML_LIST_I4, & ! integer(I4P) list
                   VTK_VAR_XML_LIST_I2, & ! integer(I2P) list
                   VTK_VAR_XML_LIST_I1    ! integer(I1P) list
endinterface
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!> @ingroup Lib_VTK_IOPrivateVarPar
!> @{
! The library uses a small set of internal variables that are private (not accessible from the outside). The following are
! private variables.
! Parameters:
integer(I4P), parameter:: maxlen  = 500      !< Max number of characters of static string.
character(1), parameter:: end_rec = char(10) !< End-character for binary-record finalize.
integer(I4P), parameter:: ascii   = 0        !< Ascii-output-format parameter identifier.
integer(I4P), parameter:: binary  = 1        !< Binary-output-format parameter identifier.
! VTK file data:
type:: Type_VTK_File
  integer(I4P)::          f        = ascii !< Current output-format (initialized to ascii format).
  character(len=maxlen):: topology = ''    !< Mesh topology.
  integer(I4P)::          u        = 0_I4P !< Logical unit.
  integer(I4P)::          ua       = 0_I4P !< Logical unit for raw binary XML append file.
#ifdef HUGE
  integer(I8P)::          N_Byte   = 0_I8P !< Number of byte to be written/read.
#else
  integer(I4P)::          N_Byte   = 0_I4P !< Number of byte to be written/read.
#endif
  integer(I8P)::          ioffset  = 0_I8P !< Offset pointer.
  integer(I4P)::          indent   = 0_I4P !< Indent pointer.
  contains
    procedure, non_overridable:: byte_update ! Procedure for updating N_Byte and ioffset pointer.
endtype Type_VTK_File
type(Type_VTK_File), allocatable:: vtk(:)       !< Global data of VTK files [1:Nvtk].
integer(I4P)::                     Nvtk = 0_I4P !< Number of (concurrent) VTK files.
integer(I4P)::                     f    = 0_I4P !< Current VTK file index.
! VTM file data:
type:: Type_VTM_File
  integer(I4P):: u      = 0_I4P !< Logical unit.
  integer(I4P):: blk    = 0_I4P !< Block index.
  integer(I4P):: indent = 0_I4P !< Indent pointer.
endtype Type_VTM_File
type(Type_VTM_File):: vtm !< Global data of VTM files.
!> @}
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  ! The library uses four auxiliary procedures that are private thus they cannot be called outside the library.
  !> @ingroup Lib_VTK_IOPrivateProcedure
  !> @{
  !> @brief Function for getting a free logic unit. The users of @libvtk does not know which is the logical
  !>        unit: @libvtk uses this information without boring the users. The logical unit used is safe-free: if the program
  !>        calling @libvtk has others logical units used @libvtk will never use these units, but will choice one that is free.
  !>@return Free_Unit
  integer function Get_Unit(Free_Unit)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer, intent(OUT), optional:: Free_Unit !< Free logic unit.
  integer::                        n1        !< Counter.
  integer::                        ios       !< Inquiring flag.
  logical::                        lopen     !< Inquiring flag.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  Get_Unit = -1
  n1=1
  do
    if ((n1/=stdout).AND.(n1/=stderr)) then
      inquire(unit=n1,opened=lopen,iostat=ios)
      if (ios==0) then
        if (.NOT.lopen) then
          Get_Unit = n1 ; if (present(Free_Unit)) Free_Unit = Get_Unit
          return
        endif
      endif
    endif
    n1=n1+1
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Get_Unit

  !> @brief Function for converting lower case characters of a string to upper case ones. @libvtk uses this function in
  !>        order to achieve case-insensitive: all character variables used within @libvtk functions are pre-processed by
  !>        Uppper_Case function before these variables are used. So the users can call @libvtk functions without pay attention of
  !>        the case of the keywords passed to the functions: calling the function VTK_INI with the string
  !>        <em>E_IO = VTK_INI('Ascii',...)</em> is equivalent to <em>E_IO = VTK_INI('ASCII',...)</em>.
  !>@return Upper_Case
  elemental function Upper_Case(string)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  character(len=*), intent(IN):: string     !< String to be converted.
  character(len=len(string))::   Upper_Case !< Converted string.
  integer::                      n1         !< Characters counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  Upper_Case = string
  do n1=1,len(string)
    select case(ichar(string(n1:n1)))
    case(97:122)
      Upper_Case(n1:n1)=char(ichar(string(n1:n1))-32) ! Upper case conversion
    endselect
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction Upper_Case

  !> @brief Subroutine for updating N_Byte and ioffset pointer.
  elemental subroutine byte_update(vtk,N_Byte)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  class(Type_VTK_File), intent(INOUT):: vtk    !< Global data of VTK file.
#ifdef HUGE
  integer(I8P),         intent(IN)::    N_Byte !< Number of bytes saved.
#else
  integer(I4P),         intent(IN)::    N_Byte !< Number of bytes saved.
#endif
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  vtk%N_Byte = N_Byte
#ifdef HUGE
  vtk%ioffset = vtk%ioffset + BYI8P + N_Byte
#else
  vtk%ioffset = vtk%ioffset + BYI4P + N_Byte
#endif
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine byte_update

  !> @brief Subroutine for updating (adding and removing elements into) vtk array.
  pure subroutine vtk_update(act,cf,Nvtk,vtk)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  character(*),                     intent(IN)::    act        !< Action: 'ADD' one more element, 'REMOVE' current element file.
  integer(I4P),                     intent(INOUT):: cf         !< Current file index (for concurrent files IO).
  integer(I4P),                     intent(INOUT):: Nvtk       !< Number of (concurrent) VTK files.
  type(Type_VTK_File), allocatable, intent(INOUT):: vtk(:)     !< VTK files data.
  type(Type_VTK_File), allocatable::                vtk_tmp(:) !< Temporary array of VTK files data.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  select case(Upper_Case(trim(act)))
  case('ADD')
    if (Nvtk>0_I4P) then
      allocate(vtk_tmp(1:Nvtk))
      vtk_tmp = vtk
      deallocate(vtk)
      Nvtk = Nvtk + 1
      allocate(vtk(1:Nvtk))
      vtk(1:Nvtk-1) = vtk_tmp
      deallocate(vtk_tmp)
      cf = Nvtk
    else
      Nvtk = 1_I4P
      allocate(vtk(1:Nvtk))
      cf = Nvtk
    endif
  case('REMOVE')
    if (Nvtk>1_I4P) then
      allocate(vtk_tmp(1:Nvtk-1))
      if (cf==Nvtk) then
        vtk_tmp = vtk(1:Nvtk-1)
      else
        vtk_tmp(1 :cf-1) = vtk(1   :cf-1)
        vtk_tmp(cf:    ) = vtk(cf+1:    )
      endif
      deallocate(vtk)
      Nvtk = Nvtk - 1
      allocate(vtk(1:Nvtk))
      vtk = vtk_tmp
      deallocate(vtk_tmp)
      cf = 1_I4P
    else
      Nvtk = 0_I4P
      if (allocated(vtk)) deallocate(vtk)
      cf = Nvtk
    endif
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine vtk_update
  !> @}

  !> @brief Function for initializing VTK-XML file.
  !> The XML standard is more powerful than legacy one. It is flexible but on the other hand is (but not so more using a library
  !> like @libvtk...) complex than legacy standard. The output of XML functions is a well-formated XML file at least for the ascii
  !> format (in the binary format @libvtk uses raw-data format that does not produce a well formated XML file).
  !> Note that the XML functions have the same name of legacy functions with the suffix \em XML.
  !> @remark This function must be the first to be called.
  !> @note An example of usage is: \n
  !> @code ...
  !> integer(I4P):: nx1,nx2,ny1,ny2,nz1,nz2
  !> ...
  !> E_IO = VTK_INI_XML('BINARY','XML_RECT_BINARY.vtr','RectilinearGrid',nx1=nx1,nx2=nx2,ny1=ny1,ny2=ny2,nz1=nz1,nz2=nz2)
  !> ... @endcode
  !> Note that the file extension is necessary in the file name. The XML standard has different extensions for each
  !> different topologies (e.g. \em vtr for rectilinear topology). See the VTK-standard file for more information.
  !> @return E_IO: integer(I4P) error flag
  function VTK_INI_XML(cf,nx1,nx2,ny1,ny2,nz1,nz2,output_format,filename,mesh_topology) result(E_IO)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(OUT), optional:: cf            !< Current file index (for concurrent files IO).
  integer(I4P), intent(IN),  optional:: nx1           !< Initial node of x axis.
  integer(I4P), intent(IN),  optional:: nx2           !< Final node of x axis.
  integer(I4P), intent(IN),  optional:: ny1           !< Initial node of y axis.
  integer(I4P), intent(IN),  optional:: ny2           !< Final node of y axis.
  integer(I4P), intent(IN),  optional:: nz1           !< Initial node of z axis.
  integer(I4P), intent(IN),  optional:: nz2           !< Final node of z axis.
  character(*), intent(IN)::            output_format !< Output format: ASCII or BINARY.
  character(*), intent(IN)::            filename      !< File name.
  character(*), intent(IN)::            mesh_topology !< Mesh topology.
  integer(I4P)::                        E_IO          !< Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done.
  character(len=maxlen)::               s_buffer      !< Buffer string.
  integer(I4P)::                        rf            !< Real file index.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.ir_initialized) call IR_Init
  call vtk_update(act='add',cf=rf,Nvtk=Nvtk,vtk=vtk)
  f = rf
  if (present(cf)) cf = rf
  vtk(rf)%topology = trim(mesh_topology)
  select case(trim(Upper_Case(output_format)))
  case('ASCII')
    vtk(rf)%f = ascii
    open(unit=Get_Unit(vtk(rf)%u),file=trim(filename),form='FORMATTED',&
         access='SEQUENTIAL',action='WRITE',status='REPLACE',iostat=E_IO)
    ! writing header of file
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)'<?xml version="1.0"?>'
    if (endian==endianL) then
      s_buffer = '<VTKFile type="'//trim(vtk(rf)%topology)//'" version="0.1" byte_order="LittleEndian">'
    else
      s_buffer = '<VTKFile type="'//trim(vtk(rf)%topology)//'" version="0.1" byte_order="BigEndian">'
    endif
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)trim(s_buffer) ; vtk(rf)%indent = 2
    select case(trim(vtk(rf)%topology))
    case('RectilinearGrid','StructuredGrid')
      s_buffer = repeat(' ',vtk(rf)%indent)//'<'//trim(vtk(rf)%topology)//' WholeExtent="'//&
                 trim(str(n=nx1))//' '//trim(str(n=nx2))//' '//                             &
                 trim(str(n=ny1))//' '//trim(str(n=ny2))//' '//                             &
                 trim(str(n=nz1))//' '//trim(str(n=nz2))//'">'
    case('UnstructuredGrid')
      s_buffer = repeat(' ',vtk(rf)%indent)//'<'//trim(vtk(rf)%topology)//'>'
    endselect
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)trim(s_buffer) ; vtk(rf)%indent = vtk(rf)%indent + 2
  case('BINARY')
    vtk(rf)%f = binary
    open(unit=Get_Unit(vtk(rf)%u),file=trim(filename),&
         form='UNFORMATTED',access='STREAM',action='WRITE',status='REPLACE',iostat=E_IO)
    ! writing header of file
    write(unit=vtk(rf)%u,iostat=E_IO)'<?xml version="1.0"?>'//end_rec
    if (endian==endianL) then
      s_buffer = '<VTKFile type="'//trim(vtk(rf)%topology)//'" version="0.1" byte_order="LittleEndian">'
    else
      s_buffer = '<VTKFile type="'//trim(vtk(rf)%topology)//'" version="0.1" byte_order="BigEndian">'
    endif
    write(unit=vtk(rf)%u,iostat=E_IO)trim(s_buffer)//end_rec ; vtk(rf)%indent = 2
    select case(trim(vtk(rf)%topology))
    case('RectilinearGrid','StructuredGrid')
      s_buffer = repeat(' ',vtk(rf)%indent)//'<'//trim(vtk(rf)%topology)//' WholeExtent="'//&
                 trim(str(n=nx1))//' '//trim(str(n=nx2))//' '//                             &
                 trim(str(n=ny1))//' '//trim(str(n=ny2))//' '//                             &
                 trim(str(n=nz1))//' '//trim(str(n=nz2))//'">'
    case('UnstructuredGrid')
      s_buffer = repeat(' ',vtk(rf)%indent)//'<'//trim(vtk(rf)%topology)//'>'
    endselect
    write(unit=vtk(rf)%u,iostat=E_IO)trim(s_buffer)//end_rec ; vtk(rf)%indent = vtk(rf)%indent + 2
    ! opening the SCRATCH file used for appending raw binary data
    open(unit=Get_Unit(vtk(rf)%ua), form='UNFORMATTED', access='STREAM', action='READWRITE', status='SCRATCH', iostat=E_IO)
    vtk(rf)%ioffset = 0 ! initializing offset pointer
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_INI_XML

  !> @ingroup Lib_VTK_IOPrivateProcedure
  !> @{
  !> Function for open/close field data tag.
  !> @return E_IO: integer(I4P) error flag
  function VTK_FLD_XML_OC(cf,fld_action) result(E_IO)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN), optional:: cf         !< Current file index (for concurrent files IO).
  character(*), intent(IN)::           fld_action !< Field data tag action: OPEN or CLOSE tag.
  integer(I4P)::                       E_IO       !< Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done.
  integer(I4P)::                       rf         !< Real file index.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  rf = f
  if (present(cf)) then
    rf = cf ; f = cf
  endif
  select case(trim(Upper_Case(fld_action)))
  case('OPEN')
    select case(vtk(rf)%f)
    case(ascii)
      write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)repeat(' ',vtk(rf)%indent)//'<FieldData>' ; vtk(rf)%indent = vtk(rf)%indent + 2
    case(binary)
      write(unit=vtk(rf)%u,iostat=E_IO)repeat(' ',vtk(rf)%indent)//'<FieldData>'//end_rec ; vtk(rf)%indent = vtk(rf)%indent + 2
    endselect
  case('CLOSE')
    select case(vtk(rf)%f)
    case(ascii)
      vtk(rf)%indent = vtk(rf)%indent - 2 ; write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)repeat(' ',vtk(rf)%indent)//'</FieldData>'
    case(binary)
      vtk(rf)%indent = vtk(rf)%indent - 2 ; write(unit=vtk(rf)%u,iostat=E_IO)repeat(' ',vtk(rf)%indent)//'</FieldData>'//end_rec
    endselect
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_FLD_XML_OC

  !> Function for saving field data (global auxiliary data, e.g. time, step number, data set name...) (R8P).
  !> @return E_IO: integer(I4P) error flag
  function VTK_FLD_XML_R8(cf,fld,fname) result(E_IO)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN), optional:: cf       !< Current file index (for concurrent files IO).
  real(R8P),    intent(IN)::           fld      !< Field data value.
  character(*), intent(IN)::           fname    !< Field data name.
  integer(I4P)::                       E_IO     !< Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done.
  character(len=maxlen)::              s_buffer !< Buffer string.
  integer(I4P)::                       rf       !< Real file index.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  rf = f
  if (present(cf)) then
    rf = cf ; f = cf
  endif
  select case(vtk(rf)%f)
  case(ascii)
    s_buffer=repeat(' ',vtk(rf)%indent)//'<DataArray type="Float64" NumberOfTuples="1" Name="'//trim(fname)//'" format="ascii">'//&
             trim(str(n=fld))//'</DataArray>'
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)trim(s_buffer)
  case(binary)
    s_buffer=repeat(' ',vtk(rf)%indent)//'<DataArray type="Float64" NumberOfTuples="1" Name="'//trim(fname)// &
             '" format="appended" offset="'//trim(str(.true.,vtk(rf)%ioffset))//'"/>'
    write(unit=vtk(rf)%u,iostat=E_IO)trim(s_buffer)//end_rec
    call vtk(rf)%byte_update(N_Byte = int(BYR8P,I4P))
    write(unit=vtk(rf)%ua,iostat=E_IO)vtk(rf)%N_Byte,'R8',1_I4P
    write(unit=vtk(rf)%ua,iostat=E_IO)fld
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_FLD_XML_R8

  !> Function for saving field data (global auxiliary data, e.g. time, step number, data set name...) (R4P).
  !> @return E_IO: integer(I4P) error flag
  function VTK_FLD_XML_R4(cf,fld,fname) result(E_IO)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN), optional:: cf       !< Current file index (for concurrent files IO).
  real(R4P),    intent(IN)::           fld      !< Field data value.
  character(*), intent(IN)::           fname    !< Field data name.
  integer(I4P)::                       E_IO     !< Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done.
  character(len=maxlen)::              s_buffer !< Buffer string.
  integer(I4P)::                       rf       !< Real file index.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  rf = f
  if (present(cf)) then
    rf = cf ; f = cf
  endif
  select case(vtk(rf)%f)
  case(ascii)
    s_buffer=repeat(' ',vtk(rf)%indent)//'<DataArray type="Float32" NumberOfTuples="1" Name="'//trim(fname)//'" format="ascii">'//&
             trim(str(n=fld))//'</DataArray>'
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)trim(s_buffer)
  case(binary)
    s_buffer=repeat(' ',vtk(rf)%indent)//'<DataArray type="Float32" NumberOfTuples="1" Name="'//trim(fname)// &
             '" format="appended" offset="'//trim(str(.true.,vtk(rf)%ioffset))//'"/>'
    write(unit=vtk(rf)%u,iostat=E_IO)trim(s_buffer)//end_rec
    call vtk(rf)%byte_update(N_Byte = int(BYR4P,I4P))
    write(unit=vtk(rf)%ua,iostat=E_IO)vtk(rf)%N_Byte,'R4',1_I4P
    write(unit=vtk(rf)%ua,iostat=E_IO)fld
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_FLD_XML_R4

  !> Function for saving field data (global auxiliary data, e.g. time, step number, data set name...) (I8P).
  !> @return E_IO: integer(I4P) error flag
  function VTK_FLD_XML_I8(cf,fld,fname) result(E_IO)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN), optional:: cf       !< Current file index (for concurrent files IO).
  integer(I8P), intent(IN)::           fld      !< Field data value.
  character(*), intent(IN)::           fname    !< Field data name.
  integer(I4P)::                       E_IO     !< Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done.
  character(len=maxlen)::              s_buffer !< Buffer string.
  integer(I4P)::                       rf       !< Real file index.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  rf = f
  if (present(cf)) then
    rf = cf ; f = cf
  endif
  select case(vtk(rf)%f)
  case(ascii)
    s_buffer = repeat(' ',vtk(rf)%indent)//'<DataArray type="Int64" NumberOfTuples="1" Name="'//trim(fname)//'" format="ascii">'// &
               trim(str(n=fld))//'</DataArray>'
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)trim(s_buffer)
  case(binary)
    s_buffer = repeat(' ',vtk(rf)%indent)//'<DataArray type="Int64" NumberOfTuples="1" Name="'//trim(fname)// &
               '" format="appended" offset="'//trim(str(.true.,vtk(rf)%ioffset))//'"/>'
    write(unit=vtk(rf)%u,iostat=E_IO)trim(s_buffer)//end_rec
    call vtk(rf)%byte_update(N_Byte = int(BYI8P,I4P))
    write(unit=vtk(rf)%ua,iostat=E_IO)vtk(rf)%N_Byte,'I8',1_I4P
    write(unit=vtk(rf)%ua,iostat=E_IO)fld
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_FLD_XML_I8

  !> Function for saving field data (global auxiliary data, e.g. time, step number, data set name...) (I4P).
  !> @return E_IO: integer(I4P) error flag
  function VTK_FLD_XML_I4(cf,fld,fname) result(E_IO)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN), optional:: cf       !< Current file index (for concurrent files IO).
  integer(I4P), intent(IN)::           fld      !< Field data value.
  character(*), intent(IN)::           fname    !< Field data name.
  integer(I4P)::                       E_IO     !< Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done.
  character(len=maxlen)::              s_buffer !< Buffer string.
  integer(I4P)::                       rf       !< Real file index.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  rf = f
  if (present(cf)) then
    rf = cf ; f = cf
  endif
  select case(vtk(rf)%f)
  case(ascii)
    s_buffer = repeat(' ',vtk(rf)%indent)//'<DataArray type="Int32" NumberOfTuples="1" Name="'//trim(fname)//'" format="ascii">'// &
               trim(str(n=fld))//'</DataArray>'
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)trim(s_buffer)
  case(binary)
    s_buffer = repeat(' ',vtk(rf)%indent)//'<DataArray type="Int32" NumberOfTuples="1" Name="'//trim(fname)// &
               '" format="appended" offset="'//trim(str(.true.,vtk(rf)%ioffset))//'"/>'
    write(unit=vtk(rf)%u,iostat=E_IO)trim(s_buffer)//end_rec
    call vtk(rf)%byte_update(N_Byte = int(BYI4P,I4P))
    write(unit=vtk(rf)%ua,iostat=E_IO)vtk(rf)%N_Byte,'I4',1_I4P
    write(unit=vtk(rf)%ua,iostat=E_IO)fld
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_FLD_XML_I4

  !> Function for saving field data (global auxiliary data, e.g. time, step number, data set name...) (I2P).
  !> @return E_IO: integer(I4P) error flag
  function VTK_FLD_XML_I2(cf,fld,fname) result(E_IO)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN), optional:: cf       !< Current file index (for concurrent files IO).
  integer(I2P), intent(IN)::           fld      !< Field data value.
  character(*), intent(IN)::           fname    !< Field data name.
  integer(I4P)::                       E_IO     !< Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done.
  character(len=maxlen)::              s_buffer !< Buffer string.
  integer(I4P)::                       rf       !< Real file index.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  rf = f
  if (present(cf)) then
    rf = cf ; f = cf
  endif
  select case(vtk(rf)%f)
  case(ascii)
    s_buffer = repeat(' ',vtk(rf)%indent)//'<DataArray type="Int16" NumberOfTuples="1" Name="'//trim(fname)//'" format="ascii">'// &
               trim(str(n=fld))//'</DataArray>'
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)trim(s_buffer)
  case(binary)
    s_buffer = repeat(' ',vtk(rf)%indent)//'<DataArray type="Int16" NumberOfTuples="1" Name="'//trim(fname)// &
               '" format="appended" offset="'//trim(str(.true.,vtk(rf)%ioffset))//'"/>'
    write(unit=vtk(rf)%u,iostat=E_IO)trim(s_buffer)//end_rec
    call vtk(rf)%byte_update(N_Byte = int(BYI2P,I4P))
    write(unit=vtk(rf)%ua,iostat=E_IO)vtk(rf)%N_Byte,'I2',1_I4P
    write(unit=vtk(rf)%ua,iostat=E_IO)fld
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_FLD_XML_I2

  !> Function for saving field data (global auxiliary data, e.g. time, step number, data set name...) (I1P).
  !> @return E_IO: integer(I4P) error flag
  function VTK_FLD_XML_I1(cf,fld,fname) result(E_IO)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN), optional:: cf       !< Current file index (for concurrent files IO).
  integer(I1P), intent(IN)::           fld      !< Field data value.
  character(*), intent(IN)::           fname    !< Field data name.
  integer(I4P)::                       E_IO     !< Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done.
  character(len=maxlen)::              s_buffer !< Buffer string.
  integer(I4P)::                       rf       !< Real file index.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  rf = f
  if (present(cf)) then
    rf = cf ; f = cf
  endif
  select case(vtk(rf)%f)
  case(ascii)
    s_buffer = repeat(' ',vtk(rf)%indent)//'<DataArray type="Int8" NumberOfTuples="1" Name="'//trim(fname)//'" format="ascii">'// &
               trim(str(n=fld))//'</DataArray>'
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)trim(s_buffer)
  case(binary)
    s_buffer = repeat(' ',vtk(rf)%indent)//'<DataArray type="Int8" NumberOfTuples="1" Name="'//trim(fname)// &
               '" format="appended" offset="'//trim(str(.true.,vtk(rf)%ioffset))//'"/>'
    write(unit=vtk(rf)%u,iostat=E_IO)trim(s_buffer)//end_rec
    call vtk(rf)%byte_update(N_Byte = int(BYI1P,I4P))
    write(unit=vtk(rf)%ua,iostat=E_IO)vtk(rf)%N_Byte,'I1',1_I4P
    write(unit=vtk(rf)%ua,iostat=E_IO)fld
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_FLD_XML_I1

  !> Function for saving mesh with \b StructuredGrid topology (R8P).
  !> @return E_IO: integer(I4P) error flag
  function VTK_GEO_XML_STRG_R8(cf,nx1,nx2,ny1,ny2,nz1,nz2,NN,X,Y,Z) result(E_IO)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN), optional:: cf       !< Current file index (for concurrent files IO).
  integer(I4P), intent(IN)::           nx1      !< Initial node of x axis.
  integer(I4P), intent(IN)::           nx2      !< Final node of x axis.
  integer(I4P), intent(IN)::           ny1      !< Initial node of y axis.
  integer(I4P), intent(IN)::           ny2      !< Final node of y axis.
  integer(I4P), intent(IN)::           nz1      !< Initial node of z axis.
  integer(I4P), intent(IN)::           nz2      !< Final node of z axis.
  integer(I4P), intent(IN)::           NN       !< Number of all nodes.
  real(R8P),    intent(IN)::           X(1:NN)  !< X coordinates.
  real(R8P),    intent(IN)::           Y(1:NN)  !< Y coordinates.
  real(R8P),    intent(IN)::           Z(1:NN)  !< Z coordinates.
  integer(I4P)::                       E_IO     !< Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done.
  character(len=maxlen)::              s_buffer !< Buffer string.
  integer(I4P)::                       rf       !< Real file index.
  integer(I4P)::                       n1       !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  rf = f
  if (present(cf)) then
    rf = cf ; f = cf
  endif
  select case(vtk(rf)%f)
  case(ascii)
    s_buffer = repeat(' ',vtk(rf)%indent)//'<Piece Extent="'//trim(str(n=nx1))//' '//trim(str(n=nx2))//' '// &
                                                      trim(str(n=ny1))//' '//trim(str(n=ny2))//' '// &
                                                      trim(str(n=nz1))//' '//trim(str(n=nz2))//'">'
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)trim(s_buffer) ; vtk(rf)%indent = vtk(rf)%indent + 2
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)repeat(' ',vtk(rf)%indent)//'<Points>' ; vtk(rf)%indent = vtk(rf)%indent + 2
    s_buffer = repeat(' ',vtk(rf)%indent)//'<DataArray type="Float64" NumberOfComponents="3" Name="Points" format="ascii">'
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)trim(s_buffer)
    write(unit=vtk(rf)%u,fmt='(3('//FR8P//',1X))',iostat=E_IO)(X(n1),Y(n1),Z(n1),n1=1,NN)
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)repeat(' ',vtk(rf)%indent)//'</DataArray>' ; vtk(rf)%indent = vtk(rf)%indent - 2
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)repeat(' ',vtk(rf)%indent)//'</Points>'
  case(binary)
    s_buffer = repeat(' ',vtk(rf)%indent)//'<Piece Extent="'//trim(str(n=nx1))//' '//trim(str(n=nx2))//' '// &
                                                      trim(str(n=ny1))//' '//trim(str(n=ny2))//' '// &
                                                      trim(str(n=nz1))//' '//trim(str(n=nz2))//'">'
    write(unit=vtk(rf)%u,iostat=E_IO)trim(s_buffer)//end_rec ; vtk(rf)%indent = vtk(rf)%indent + 2
    write(unit=vtk(rf)%u,iostat=E_IO)repeat(' ',vtk(rf)%indent)//'<Points>'//end_rec ; vtk(rf)%indent = vtk(rf)%indent + 2
    s_buffer = repeat(' ',vtk(rf)%indent)//                                                                  &
               '<DataArray type="Float64" NumberOfComponents="3" Name="Points" format="appended" offset="'// &
               trim(str(.true.,vtk(rf)%ioffset))//'"/>'
    write(unit=vtk(rf)%u,iostat=E_IO)trim(s_buffer)//end_rec
    call vtk(rf)%byte_update(N_Byte = 3*NN*BYR8P)
    write(unit=vtk(rf)%ua,iostat=E_IO)vtk(rf)%N_Byte,'R8',3*NN
    write(unit=vtk(rf)%ua,iostat=E_IO)(X(n1),Y(n1),Z(n1),n1=1,NN)
    vtk(rf)%indent = vtk(rf)%indent - 2 ; write(unit=vtk(rf)%u,iostat=E_IO)repeat(' ',vtk(rf)%indent)//'</Points>'//end_rec
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_GEO_XML_STRG_R8

  !> Function for saving mesh with \b StructuredGrid topology (R4P).
  !> @return E_IO: integer(I4P) error flag
  function VTK_GEO_XML_STRG_R4(cf,nx1,nx2,ny1,ny2,nz1,nz2,NN,X,Y,Z) result(E_IO)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN), optional:: cf       !< Current file index (for concurrent files IO).
  integer(I4P), intent(IN)::           nx1      !< Initial node of x axis.
  integer(I4P), intent(IN)::           nx2      !< Final node of x axis.
  integer(I4P), intent(IN)::           ny1      !< Initial node of y axis.
  integer(I4P), intent(IN)::           ny2      !< Final node of y axis.
  integer(I4P), intent(IN)::           nz1      !< Initial node of z axis.
  integer(I4P), intent(IN)::           nz2      !< Final node of z axis.
  integer(I4P), intent(IN)::           NN       !< Number of all nodes.
  real(R4P),    intent(IN)::           X(1:NN)  !< X coordinates.
  real(R4P),    intent(IN)::           Y(1:NN)  !< Y coordinates.
  real(R4P),    intent(IN)::           Z(1:NN)  !< Z coordinates.
  integer(I4P)::                       E_IO     !< Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done.
  character(len=maxlen)::              s_buffer !< Buffer string.
  integer(I4P)::                       rf       !< Real file index.
  integer(I4P)::                       n1       !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  rf = f
  if (present(cf)) then
    rf = cf ; f = cf
  endif
  select case(vtk(rf)%f)
  case(ascii)
    s_buffer = repeat(' ',vtk(rf)%indent)//'<Piece Extent="'//trim(str(n=nx1))//' '//trim(str(n=nx2))//' '// &
                                                      trim(str(n=ny1))//' '//trim(str(n=ny2))//' '// &
                                                      trim(str(n=nz1))//' '//trim(str(n=nz2))//'">'
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)trim(s_buffer) ; vtk(rf)%indent = vtk(rf)%indent + 2
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)repeat(' ',vtk(rf)%indent)//'<Points>' ; vtk(rf)%indent = vtk(rf)%indent + 2
    s_buffer = repeat(' ',vtk(rf)%indent)//'<DataArray type="Float32" NumberOfComponents="3" Name="Points" format="ascii">'
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)trim(s_buffer)
    write(unit=vtk(rf)%u,fmt='(3('//FR4P//',1X))',iostat=E_IO)(X(n1),Y(n1),Z(n1),n1=1,NN)
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)repeat(' ',vtk(rf)%indent)//'</DataArray>' ; vtk(rf)%indent = vtk(rf)%indent - 2
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)repeat(' ',vtk(rf)%indent)//'</Points>'
  case(binary)
    s_buffer = repeat(' ',vtk(rf)%indent)//'<Piece Extent="'//trim(str(n=nx1))//' '//trim(str(n=nx2))//' '// &
                                                      trim(str(n=ny1))//' '//trim(str(n=ny2))//' '// &
                                                      trim(str(n=nz1))//' '//trim(str(n=nz2))//'">'
    write(unit=vtk(rf)%u,iostat=E_IO)trim(s_buffer)//end_rec ; vtk(rf)%indent = vtk(rf)%indent + 2
    write(unit=vtk(rf)%u,iostat=E_IO)repeat(' ',vtk(rf)%indent)//'<Points>'//end_rec ; vtk(rf)%indent = vtk(rf)%indent + 2
    s_buffer = repeat(' ',vtk(rf)%indent)//                                                                  &
               '<DataArray type="Float32" NumberOfComponents="3" Name="Points" format="appended" offset="'// &
               trim(str(.true.,vtk(rf)%ioffset))//'"/>'
    write(unit=vtk(rf)%u,iostat=E_IO)trim(s_buffer)//end_rec
    call vtk(rf)%byte_update(N_Byte = 3*NN*BYR4P)
    write(unit=vtk(rf)%ua,iostat=E_IO)vtk(rf)%N_Byte,'R4',3*NN
    write(unit=vtk(rf)%ua,iostat=E_IO)(X(n1),Y(n1),Z(n1),n1=1,NN)
    vtk(rf)%indent = vtk(rf)%indent - 2 ; write(unit=vtk(rf)%u,iostat=E_IO)repeat(' ',vtk(rf)%indent)//'</Points>'//end_rec
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_GEO_XML_STRG_R4

  !> Function for saving mesh with \b RectilinearGrid topology (R8P).
  !> @return E_IO: integer(I4P) error flag
  function VTK_GEO_XML_RECT_R8(cf,nx1,nx2,ny1,ny2,nz1,nz2,X,Y,Z) result(E_IO)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN), optional:: cf         !< Current file index (for concurrent files IO).
  integer(I4P), intent(IN)::           nx1        !< Initial node of x axis.
  integer(I4P), intent(IN)::           nx2        !< Final node of x axis.
  integer(I4P), intent(IN)::           ny1        !< Initial node of y axis.
  integer(I4P), intent(IN)::           ny2        !< Final node of y axis.
  integer(I4P), intent(IN)::           nz1        !< Initial node of z axis.
  integer(I4P), intent(IN)::           nz2        !< Final node of z axis.
  real(R8P),    intent(IN)::           X(nx1:nx2) !< X coordinates.
  real(R8P),    intent(IN)::           Y(ny1:ny2) !< Y coordinates.
  real(R8P),    intent(IN)::           Z(nz1:nz2) !< Z coordinates.
  integer(I4P)::                       E_IO       !< Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done.
  character(len=maxlen)::              s_buffer   !< Buffer string.
  integer(I4P)::                       rf         !< Real file index.
  integer(I4P)::                       n1         !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  rf = f
  if (present(cf)) then
    rf = cf ; f = cf
  endif
  select case(vtk(rf)%f)
  case(ascii)
    s_buffer = repeat(' ',vtk(rf)%indent)//'<Piece Extent="'//trim(str(n=nx1))//' '//trim(str(n=nx2))//' '// &
                                                      trim(str(n=ny1))//' '//trim(str(n=ny2))//' '// &
                                                      trim(str(n=nz1))//' '//trim(str(n=nz2))//'">'
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)trim(s_buffer) ; vtk(rf)%indent = vtk(rf)%indent + 2
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)repeat(' ',vtk(rf)%indent)//'<Coordinates>' ; vtk(rf)%indent = vtk(rf)%indent + 2
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)repeat(' ',vtk(rf)%indent)//'<DataArray type="Float64" Name="X" format="ascii">'
    write(unit=vtk(rf)%u,fmt=FR8P, iostat=E_IO)(X(n1),n1=nx1,nx2)
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)repeat(' ',vtk(rf)%indent)//'</DataArray>'
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)repeat(' ',vtk(rf)%indent)//'<DataArray type="Float64" Name="Y" format="ascii">'
    write(unit=vtk(rf)%u,fmt=FR8P, iostat=E_IO)(Y(n1),n1=ny1,ny2)
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)repeat(' ',vtk(rf)%indent)//'</DataArray>'
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)repeat(' ',vtk(rf)%indent)//'<DataArray type="Float64" Name="Z" format="ascii">'
    write(unit=vtk(rf)%u,fmt=FR8P, iostat=E_IO)(Z(n1),n1=nz1,nz2)
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)repeat(' ',vtk(rf)%indent)//'</DataArray>' ; vtk(rf)%indent = vtk(rf)%indent - 2
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)repeat(' ',vtk(rf)%indent)//'</Coordinates>'
  case(binary)
    s_buffer = repeat(' ',vtk(rf)%indent)//'<Piece Extent="'//trim(str(n=nx1))//' '//trim(str(n=nx2))//' '// &
                                                      trim(str(n=ny1))//' '//trim(str(n=ny2))//' '// &
                                                      trim(str(n=nz1))//' '//trim(str(n=nz2))//'">'
    write(unit=vtk(rf)%u,iostat=E_IO)trim(s_buffer)//end_rec ; vtk(rf)%indent = vtk(rf)%indent + 2
    write(unit=vtk(rf)%u,iostat=E_IO)repeat(' ',vtk(rf)%indent)//'<Coordinates>'//end_rec ; vtk(rf)%indent = vtk(rf)%indent + 2
    s_buffer = repeat(' ',vtk(rf)%indent)//'<DataArray type="Float64" Name="X" format="appended" offset="'//&
               trim(str(.true.,vtk(rf)%ioffset))//'"/>'
    write(unit=vtk(rf)%u,iostat=E_IO)trim(s_buffer)//end_rec
    call vtk(rf)%byte_update(N_Byte = (nx2-nx1+1)*BYR8P)
    write(unit=vtk(rf)%ua,iostat=E_IO)vtk(rf)%N_Byte,'R8',(nx2-nx1+1)
    write(unit=vtk(rf)%ua,iostat=E_IO)(X(n1),n1=nx1,nx2)
    s_buffer = repeat(' ',vtk(rf)%indent)//'<DataArray type="Float64" Name="Y" format="appended" offset="'//&
               trim(str(.true.,vtk(rf)%ioffset))//'"/>'
    write(unit=vtk(rf)%u,iostat=E_IO)trim(s_buffer)//end_rec
    call vtk(rf)%byte_update(N_Byte = (ny2-ny1+1)*BYR8P)
    write(unit=vtk(rf)%ua,iostat=E_IO)vtk(rf)%N_Byte,'R8',(ny2-ny1+1)
    write(unit=vtk(rf)%ua,iostat=E_IO)(Y(n1),n1=ny1,ny2)
    s_buffer = repeat(' ',vtk(rf)%indent)//'<DataArray type="Float64" Name="Z" format="appended" offset="'//&
               trim(str(.true.,vtk(rf)%ioffset))//'"/>'
    write(unit=vtk(rf)%u,iostat=E_IO)trim(s_buffer)//end_rec
    call vtk(rf)%byte_update(N_Byte = (nz2-nz1+1)*BYR8P)
    write(unit=vtk(rf)%ua,iostat=E_IO)vtk(rf)%N_Byte,'R8',(nz2-nz1+1)
    write(unit=vtk(rf)%ua,iostat=E_IO)(Z(n1),n1=nz1,nz2)
    vtk(rf)%indent = vtk(rf)%indent - 2 ; write(unit=vtk(rf)%u,iostat=E_IO)repeat(' ',vtk(rf)%indent)//'</Coordinates>'//end_rec
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_GEO_XML_RECT_R8

  !> Function for saving mesh with \b RectilinearGrid topology (R4P).
  !> @return E_IO: integer(I4P) error flag
  function VTK_GEO_XML_RECT_R4(cf,nx1,nx2,ny1,ny2,nz1,nz2,X,Y,Z) result(E_IO)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN), optional:: cf         !< Current file index (for concurrent files IO).
  integer(I4P), intent(IN)::           nx1        !< Initial node of x axis.
  integer(I4P), intent(IN)::           nx2        !< Final node of x axis.
  integer(I4P), intent(IN)::           ny1        !< Initial node of y axis.
  integer(I4P), intent(IN)::           ny2        !< Final node of y axis.
  integer(I4P), intent(IN)::           nz1        !< Initial node of z axis.
  integer(I4P), intent(IN)::           nz2        !< Final node of z axis.
  real(R4P),    intent(IN)::           X(nx1:nx2) !< X coordinates.
  real(R4P),    intent(IN)::           Y(ny1:ny2) !< Y coordinates.
  real(R4P),    intent(IN)::           Z(nz1:nz2) !< Z coordinates.
  integer(I4P)::                       E_IO       !< Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done.
  character(len=maxlen)::              s_buffer   !< Buffer string.
  integer(I4P)::                       rf         !< Real file index.
  integer(I4P)::                       n1         !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  rf = f
  if (present(cf)) then
    rf = cf ; f = cf
  endif
  select case(vtk(rf)%f)
  case(ascii)
    s_buffer = repeat(' ',vtk(rf)%indent)//'<Piece Extent="'//trim(str(n=nx1))//' '//trim(str(n=nx2))//' '// &
                                                      trim(str(n=ny1))//' '//trim(str(n=ny2))//' '// &
                                                      trim(str(n=nz1))//' '//trim(str(n=nz2))//'">'
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)trim(s_buffer) ; vtk(rf)%indent = vtk(rf)%indent + 2
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)repeat(' ',vtk(rf)%indent)//'<Coordinates>' ; vtk(rf)%indent = vtk(rf)%indent + 2
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)repeat(' ',vtk(rf)%indent)//'<DataArray type="Float32" Name="X" format="ascii">'
    write(unit=vtk(rf)%u,fmt=FR4P, iostat=E_IO)(X(n1),n1=nx1,nx2)
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)repeat(' ',vtk(rf)%indent)//'</DataArray>'
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)repeat(' ',vtk(rf)%indent)//'<DataArray type="Float32" Name="Y" format="ascii">'
    write(unit=vtk(rf)%u,fmt=FR4P, iostat=E_IO)(Y(n1),n1=ny1,ny2)
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)repeat(' ',vtk(rf)%indent)//'</DataArray>'
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)repeat(' ',vtk(rf)%indent)//'<DataArray type="Float32" Name="Z" format="ascii">'
    write(unit=vtk(rf)%u,fmt=FR4P, iostat=E_IO)(Z(n1),n1=nz1,nz2)
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)repeat(' ',vtk(rf)%indent)//'</DataArray>' ; vtk(rf)%indent = vtk(rf)%indent - 2
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)repeat(' ',vtk(rf)%indent)//'</Coordinates>'
  case(binary)
    s_buffer = repeat(' ',vtk(rf)%indent)//'<Piece Extent="'//trim(str(n=nx1))//' '//trim(str(n=nx2))//' '// &
                                                      trim(str(n=ny1))//' '//trim(str(n=ny2))//' '// &
                                                      trim(str(n=nz1))//' '//trim(str(n=nz2))//'">'
    write(unit=vtk(rf)%u,iostat=E_IO)trim(s_buffer)//end_rec ; vtk(rf)%indent = vtk(rf)%indent + 2
    write(unit=vtk(rf)%u,iostat=E_IO)repeat(' ',vtk(rf)%indent)//'<Coordinates>'//end_rec ; vtk(rf)%indent = vtk(rf)%indent + 2
    s_buffer = repeat(' ',vtk(rf)%indent)//'<DataArray type="Float32" Name="X" format="appended" offset="'//&
               trim(str(.true.,vtk(rf)%ioffset))//'"/>'
    write(unit=vtk(rf)%u,iostat=E_IO)trim(s_buffer)//end_rec
    call vtk(rf)%byte_update(N_Byte = (nx2-nx1+1)*BYR4P)
    write(unit=vtk(rf)%ua,iostat=E_IO)vtk(rf)%N_Byte,'R4',(nx2-nx1+1)
    write(unit=vtk(rf)%ua,iostat=E_IO)(X(n1),n1=nx1,nx2)
    s_buffer = repeat(' ',vtk(rf)%indent)//'<DataArray type="Float32" Name="Y" format="appended" offset="'//&
               trim(str(.true.,vtk(rf)%ioffset))//'"/>'
    write(unit=vtk(rf)%u,iostat=E_IO)trim(s_buffer)//end_rec
    call vtk(rf)%byte_update(N_Byte = (ny2-ny1+1)*BYR4P)
    write(unit=vtk(rf)%ua,iostat=E_IO)vtk(rf)%N_Byte,'R4',(ny2-ny1+1)
    write(unit=vtk(rf)%ua,iostat=E_IO)(Y(n1),n1=ny1,ny2)
    s_buffer = repeat(' ',vtk(rf)%indent)//'<DataArray type="Float32" Name="Z" format="appended" offset="'//&
               trim(str(.true.,vtk(rf)%ioffset))//'"/>'
    write(unit=vtk(rf)%u,iostat=E_IO)trim(s_buffer)//end_rec
    call vtk(rf)%byte_update(N_Byte = (nz2-nz1+1)*BYR4P)
    write(unit=vtk(rf)%ua,iostat=E_IO)vtk(rf)%N_Byte,'R4',(nz2-nz1+1)
    write(unit=vtk(rf)%ua,iostat=E_IO)(Z(n1),n1=nz1,nz2)
    vtk(rf)%indent = vtk(rf)%indent - 2  ; write(unit=vtk(rf)%u,iostat=E_IO)repeat(' ',vtk(rf)%indent)//'</Coordinates>'//end_rec
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_GEO_XML_RECT_R4

  !> Function for saving mesh with \b UnstructuredGrid topology (R8P).
  !> @return E_IO: integer(I4P) error flag
  function VTK_GEO_XML_UNST_R8(cf,NN,NC,X,Y,Z) result(E_IO)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN), optional:: cf       !< Current file index (for concurrent files IO).
  integer(I4P), intent(IN)::           NN       !< Number of nodes.
  integer(I4P), intent(IN)::           NC       !< Number of cells.
  real(R8P),    intent(IN)::           X(1:NN)  !< X coordinates.
  real(R8P),    intent(IN)::           Y(1:NN)  !< Y coordinates.
  real(R8P),    intent(IN)::           Z(1:NN)  !< Z coordinates.
  integer(I4P)::                       E_IO     !< Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done.
  character(len=maxlen)::              s_buffer !< Buffer string.
  integer(I4P)::                       rf       !< Real file index.
  integer(I4P)::                       n1       !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  rf = f
  if (present(cf)) then
    rf = cf ; f = cf
  endif
  select case(vtk(rf)%f)
  case(ascii)
    s_buffer = repeat(' ',vtk(rf)%indent)//'<Piece NumberOfPoints="'//trim(str(n=NN))//'" NumberOfCells="'//trim(str(n=NC))//'">'
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)trim(s_buffer) ; vtk(rf)%indent = vtk(rf)%indent + 2
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)repeat(' ',vtk(rf)%indent)//'<Points>' ; vtk(rf)%indent = vtk(rf)%indent + 2
    s_buffer = repeat(' ',vtk(rf)%indent)//'<DataArray type="Float64" NumberOfComponents="3" Name="Points" format="ascii">'
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)trim(s_buffer)
    write(unit=vtk(rf)%u,fmt='(3('//FR8P//',1X))',iostat=E_IO)(X(n1),Y(n1),Z(n1),n1=1,NN)
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)repeat(' ',vtk(rf)%indent)//'</DataArray>' ; vtk(rf)%indent = vtk(rf)%indent - 2
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)repeat(' ',vtk(rf)%indent)//'</Points>'
  case(binary)
    s_buffer = repeat(' ',vtk(rf)%indent)//'<Piece NumberOfPoints="'//trim(str(n=NN))//'" NumberOfCells="'//trim(str(n=NC))//'">'
    write(unit=vtk(rf)%u,iostat=E_IO)trim(s_buffer)//end_rec ; vtk(rf)%indent = vtk(rf)%indent + 2
    write(unit=vtk(rf)%u,iostat=E_IO)repeat(' ',vtk(rf)%indent)//'<Points>'//end_rec ; vtk(rf)%indent = vtk(rf)%indent + 2
    s_buffer = repeat(' ',vtk(rf)%indent)//                                                                  &
               '<DataArray type="Float64" NumberOfComponents="3" Name="Points" format="appended" offset="'// &
               trim(str(.true.,vtk(rf)%ioffset))//'"/>'
    write(unit=vtk(rf)%u,iostat=E_IO)trim(s_buffer)//end_rec
    call vtk(rf)%byte_update(N_Byte = 3*NN*BYR8P)
    write(unit=vtk(rf)%ua,iostat=E_IO)vtk(rf)%N_Byte,'R8',3*NN
    write(unit=vtk(rf)%ua,iostat=E_IO)(X(n1),Y(n1),Z(n1),n1=1,NN)
    vtk(rf)%indent = vtk(rf)%indent - 2 ; write(unit=vtk(rf)%u,iostat=E_IO)repeat(' ',vtk(rf)%indent)//'</Points>'//end_rec
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_GEO_XML_UNST_R8

  !> Function for saving mesh with \b UnstructuredGrid topology (R4P).
  !> @return E_IO: integer(I4P) error flag
  function VTK_GEO_XML_UNST_R4(cf,NN,NC,X,Y,Z) result(E_IO)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN), optional:: cf       !< Current file index (for concurrent files IO).
  integer(I4P), intent(IN)::           NN       !< Number of nodes.
  integer(I4P), intent(IN)::           NC       !< Number of cells.
  real(R4P),    intent(IN)::           X(1:NN)  !< X coordinates.
  real(R4P),    intent(IN)::           Y(1:NN)  !< Y coordinates.
  real(R4P),    intent(IN)::           Z(1:NN)  !< Z coordinates.
  integer(I4P)::                       E_IO     !< Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done.
  character(len=maxlen)::              s_buffer !< Buffer string.
  integer(I4P)::                       rf       !< Real file index.
  integer(I4P)::                       n1       !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  rf = f
  if (present(cf)) then
    rf = cf ; f = cf
  endif
  select case(vtk(rf)%f)
  case(ascii)
    s_buffer = repeat(' ',vtk(rf)%indent)//'<Piece NumberOfPoints="'//trim(str(n=NN))//'" NumberOfCells="'//trim(str(n=NC))//'">'
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)trim(s_buffer) ; vtk(rf)%indent = vtk(rf)%indent + 2
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)repeat(' ',vtk(rf)%indent)//'<Points>' ; vtk(rf)%indent = vtk(rf)%indent + 2
    s_buffer = repeat(' ',vtk(rf)%indent)//'<DataArray type="Float32" NumberOfComponents="3" Name="Points" format="ascii">'
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)trim(s_buffer)
    write(unit=vtk(rf)%u,fmt='(3('//FR4P//',1X))',iostat=E_IO)(X(n1),Y(n1),Z(n1),n1=1,NN)
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)repeat(' ',vtk(rf)%indent)//'</DataArray>' ; vtk(rf)%indent = vtk(rf)%indent - 2
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)repeat(' ',vtk(rf)%indent)//'</Points>'
  case(binary)
    s_buffer = repeat(' ',vtk(rf)%indent)//'<Piece NumberOfPoints="'//trim(str(n=NN))//'" NumberOfCells="'//trim(str(n=NC))//'">'
    write(unit=vtk(rf)%u,iostat=E_IO)trim(s_buffer)//end_rec ; vtk(rf)%indent = vtk(rf)%indent + 2
    write(unit=vtk(rf)%u,iostat=E_IO)repeat(' ',vtk(rf)%indent)//'<Points>'//end_rec ; vtk(rf)%indent = vtk(rf)%indent + 2
    s_buffer = repeat(' ',vtk(rf)%indent)//                                                                  &
               '<DataArray type="Float32" NumberOfComponents="3" Name="Points" format="appended" offset="'// &
               trim(str(.true.,vtk(rf)%ioffset))//'"/>'
    write(unit=vtk(rf)%u,iostat=E_IO)trim(s_buffer)//end_rec
    call vtk(rf)%byte_update(N_Byte = 3*NN*BYR4P)
    write(unit=vtk(rf)%ua,iostat=E_IO)vtk(rf)%N_Byte,'R4',3*NN
    write(unit=vtk(rf)%ua,iostat=E_IO)(X(n1),Y(n1),Z(n1),n1=1,NN)
    vtk(rf)%indent = vtk(rf)%indent - 2 ; write(unit=vtk(rf)%u,iostat=E_IO)repeat(' ',vtk(rf)%indent)//'</Points>'//end_rec
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_GEO_XML_UNST_R4

  !> @brief Function for closing mesh block data.
  !> @return E_IO: integer(I4P) error flag
  function VTK_GEO_XML_CLOSEP(cf) result(E_IO)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN), optional:: cf   !< Current file index (for concurrent files IO).
  integer(I4P)::                       E_IO !< Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done.
  integer(I4P)::                       rf   !< Real file index.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  rf = f
  if (present(cf)) then
    rf = cf ; f = cf
  endif
  vtk(rf)%indent = vtk(rf)%indent - 2
  select case(vtk(rf)%f)
  case(ascii)
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)repeat(' ',vtk(rf)%indent)//'</Piece>'
  case(binary)
    write(unit=vtk(rf)%u,iostat=E_IO)repeat(' ',vtk(rf)%indent)//'</Piece>'//end_rec
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_GEO_XML_CLOSEP
  !> @}

  !> Function that \b must be used when unstructured grid is used, it saves the connectivity of the unstructured gird.
  !> @note The vector \b connect must follow the VTK-legacy standard. It is passed as \em assumed-shape array
  !> because its dimensions is related to the mesh dimensions in a complex way. Its dimensions can be calculated by the following
  !> equation: \f$dc = dc = \sum\limits_{i = 1}^{NC} {nvertex_i }\f$.
  !> Note that this equation is different from the legacy one. The XML connectivity convention is quite different from the
  !> legacy standard. As an example suppose we have a mesh composed by 2 cells, one hexahedron (8 vertices) and one pyramid with
  !> square basis (5 vertices) and suppose that the basis of pyramid is constitute by a face of the hexahedron and so the two cells
  !> share 4 vertices. The above equation gives \f$dc=8+5=13\f$. The connectivity vector for this mesh can be: \n
  !> first cell \n
  !> connect(1)  = 0 identification flag of \f$1^\circ\f$ vertex of 1° cell \n
  !> connect(2)  = 1 identification flag of \f$2^\circ\f$ vertex of 1° cell \n
  !> connect(3)  = 2 identification flag of \f$3^\circ\f$ vertex of 1° cell \n
  !> connect(4)  = 3 identification flag of \f$4^\circ\f$ vertex of 1° cell \n
  !> connect(5)  = 4 identification flag of \f$5^\circ\f$ vertex of 1° cell \n
  !> connect(6)  = 5 identification flag of \f$6^\circ\f$ vertex of 1° cell \n
  !> connect(7)  = 6 identification flag of \f$7^\circ\f$ vertex of 1° cell \n
  !> connect(8)  = 7 identification flag of \f$8^\circ\f$ vertex of 1° cell \n
  !> second cell \n
  !> connect(9 ) = 0 identification flag of \f$1^\circ\f$ vertex of 2° cell \n
  !> connect(10) = 1 identification flag of \f$2^\circ\f$ vertex of 2° cell \n
  !> connect(11) = 2 identification flag of \f$3^\circ\f$ vertex of 2° cell \n
  !> connect(12) = 3 identification flag of \f$4^\circ\f$ vertex of 2° cell \n
  !> connect(13) = 8 identification flag of \f$5^\circ\f$ vertex of 2° cell \n
  !> Therefore this connectivity vector convention is more simple than the legacy convention, now we must create also the
  !> \em offset vector that contains the data now missing in the \em connect vector. The offset
  !> vector for this mesh can be: \n
  !> first cell \n
  !> offset(1) = 8  => summ of nodes of \f$1^\circ\f$ cell \n
  !> second cell \n
  !> offset(2) = 13 => summ of nodes of \f$1^\circ\f$ and \f$2^\circ\f$ cells \n
  !> The value of every cell-offset can be calculated by the following equation: \f$offset_c=\sum\limits_{i=1}^{c}{nvertex_i}\f$
  !> where \f$offset_c\f$ is the value of \f$c^{th}\f$ cell and \f$nvertex_i\f$ is the number of vertices of \f$i^{th}\f$ cell.
  !> The function VTK_CON_XML does not calculate the connectivity and offset vectors: it writes the connectivity and offset
  !> vectors conforming the VTK-XML standard, but does not calculate them.
  !> The vector variable \em cell_type must conform the VTK-XML standard (see the file VTK-Standard at the
  !> Kitware homepage) that is the same of the legacy standard. It contains the
  !> \em type of each cells. For the above example this vector is: \n
  !> first cell \n
  !> cell_type(1) = 12 hexahedron type of \f$1^\circ\f$ cell \n
  !> second cell \n
  !> cell_type(2) = 14 pyramid type of \f$2^\circ\f$ cell \n
  !> @return E_IO: integer(I4P) error flag
  function VTK_CON_XML(cf,NC,connect,offset,cell_type) result(E_IO)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN), optional:: cf              !< Current file index (for concurrent files IO).
  integer(I4P), intent(IN)::           NC              !< Number of cells.
  integer(I4P), intent(IN)::           connect(:)      !< Mesh connectivity.
  integer(I4P), intent(IN)::           offset(1:NC)    !< Cell offset.
  integer(I1P), intent(IN)::           cell_type(1:NC) !< VTK cell type.
  integer(I4P)::                       E_IO            !< Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done.
  character(len=maxlen)::              s_buffer        !< Buffer string.
  integer(I4P)::                       rf              !< Real file index.
  integer(I4P)::                       n1              !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  rf = f
  if (present(cf)) then
    rf = cf ; f = cf
  endif
  select case(vtk(rf)%f)
  case(ascii)
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)repeat(' ',vtk(rf)%indent)//'<Cells>' ; vtk(rf)%indent = vtk(rf)%indent + 2
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)repeat(' ',vtk(rf)%indent)//&
                                              '<DataArray type="Int32" Name="connectivity" format="ascii">'
    write(unit=vtk(rf)%u,fmt=FI4P, iostat=E_IO)(connect(n1),n1=1,size(connect))
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)repeat(' ',vtk(rf)%indent)//'</DataArray>'
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)repeat(' ',vtk(rf)%indent)//'<DataArray type="Int32" Name="offsets" format="ascii">'
    write(unit=vtk(rf)%u,fmt=FI4P, iostat=E_IO)(offset(n1),n1=1,NC)
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)repeat(' ',vtk(rf)%indent)//'</DataArray>'
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)repeat(' ',vtk(rf)%indent)//'<DataArray type="Int8" Name="types" format="ascii">'
    write(unit=vtk(rf)%u,fmt=FI1P, iostat=E_IO)(cell_type(n1),n1=1,NC)
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)repeat(' ',vtk(rf)%indent)//'</DataArray>' ; vtk(rf)%indent = vtk(rf)%indent - 2
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)repeat(' ',vtk(rf)%indent)//'</Cells>'
  case(binary)
    write(unit=vtk(rf)%u,iostat=E_IO)repeat(' ',vtk(rf)%indent)//'<Cells>'//end_rec ; vtk(rf)%indent = vtk(rf)%indent + 2
    s_buffer = repeat(' ',vtk(rf)%indent)//'<DataArray type="Int32" Name="connectivity" format="appended" offset="'// &
               trim(str(.true.,vtk(rf)%ioffset))//'"/>'
    write(unit=vtk(rf)%u,iostat=E_IO)trim(s_buffer)//end_rec
    call vtk(rf)%byte_update(N_Byte = size(connect)*BYI4P)
    write(unit=vtk(rf)%ua,iostat=E_IO)vtk(rf)%N_Byte,'I4',size(connect)
    write(unit=vtk(rf)%ua,iostat=E_IO)(connect(n1),n1=1,size(connect))
    s_buffer = repeat(' ',vtk(rf)%indent)//'<DataArray type="Int32" Name="offsets" format="appended" offset="'// &
               trim(str(.true.,vtk(rf)%ioffset))//'"/>'
    write(unit=vtk(rf)%u,iostat=E_IO)trim(s_buffer)//end_rec
    call vtk(rf)%byte_update(N_Byte = NC*BYI4P)
    write(unit=vtk(rf)%ua,iostat=E_IO)vtk(rf)%N_Byte,'I4',NC
    write(unit=vtk(rf)%ua,iostat=E_IO)(offset(n1),n1=1,NC)
    s_buffer = repeat(' ',vtk(rf)%indent)//'<DataArray type="Int8" Name="types" format="appended" offset="'// &
               trim(str(.true.,vtk(rf)%ioffset))//'"/>'
    write(unit=vtk(rf)%u,iostat=E_IO)trim(s_buffer)//end_rec
    call vtk(rf)%byte_update(N_Byte = NC*BYI1P)
    write(unit=vtk(rf)%ua,iostat=E_IO)vtk(rf)%N_Byte,'I1',NC
    write(unit=vtk(rf)%ua,iostat=E_IO)(cell_type(n1),n1=1,NC)
    vtk(rf)%indent = vtk(rf)%indent - 2 ; write(unit=vtk(rf)%u,iostat=E_IO)repeat(' ',vtk(rf)%indent)//'</Cells>'//end_rec
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_CON_XML

  !> Function that \b must be called before saving the data related to geometric mesh, this function initializes the
  !> saving of data variables indicating the \em type (node or cell centered) of variables that will be saved.
  !> @note A single file can contain both cell and node centered variables. In this case the VTK_DAT_XML function must be
  !> called two times, before saving cell-centered variables and before saving node-centered variables.
  !> Examples of usage are: \n
  !> \b Opening node piece: \n
  !> @code ...
  !> E_IO=VTK_DAT_XML('node','OPEN')
  !> ... @endcode
  !> \b Closing node piece: \n
  !> @code ...
  !> E_IO=VTK_DAT_XML('node','CLOSE')
  !> ... @endcode
  !> @return E_IO: integer(I4P) error flag
  !> @ingroup Lib_VTK_IOPublicProcedure
  function VTK_DAT_XML(cf,var_location,var_block_action) result(E_IO)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN), optional:: cf               !< Current file index (for concurrent files IO).
  character(*), intent(IN)::           var_location     !< Location of saving variables: CELL or NODE centered.
  character(*), intent(IN)::           var_block_action !< Variables block action: OPEN or CLOSE block.
  integer(I4P)::                       E_IO             !< Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done.
  integer(I4P)::                       rf               !< Real file index.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  rf = f
  if (present(cf)) then
    rf = cf ; f = cf
  endif
  select case(vtk(rf)%f)
  case(ascii)
    select case(trim(Upper_Case(var_location)))
    case('CELL')
      select case(trim(Upper_Case(var_block_action)))
      case('OPEN')
        write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)repeat(' ',vtk(rf)%indent)//'<CellData>' ; vtk(rf)%indent = vtk(rf)%indent + 2
      case('CLOSE')
        vtk(rf)%indent = vtk(rf)%indent - 2 ; write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)repeat(' ',vtk(rf)%indent)//'</CellData>'
      endselect
    case('NODE')
      select case(trim(Upper_Case(var_block_action)))
      case('OPEN')
        write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)repeat(' ',vtk(rf)%indent)//'<PointData>' ; vtk(rf)%indent = vtk(rf)%indent + 2
      case('CLOSE')
        vtk(rf)%indent = vtk(rf)%indent - 2 ; write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)repeat(' ',vtk(rf)%indent)//'</PointData>'
      endselect
    endselect
  case(binary)
    select case(trim(Upper_Case(var_location)))
    case('CELL')
      select case(trim(Upper_Case(var_block_action)))
      case('OPEN')
        write(unit=vtk(rf)%u,iostat=E_IO)repeat(' ',vtk(rf)%indent)//'<CellData>'//end_rec ; vtk(rf)%indent = vtk(rf)%indent + 2
      case('CLOSE')
        vtk(rf)%indent = vtk(rf)%indent - 2 ; write(unit=vtk(rf)%u,iostat=E_IO)repeat(' ',vtk(rf)%indent)//'</CellData>'//end_rec
      endselect
    case('NODE')
      select case(trim(Upper_Case(var_block_action)))
      case('OPEN')
        write(unit=vtk(rf)%u,iostat=E_IO)repeat(' ',vtk(rf)%indent)//'<PointData>'//end_rec ; vtk(rf)%indent = vtk(rf)%indent + 2
      case('CLOSE')
        vtk(rf)%indent = vtk(rf)%indent - 2 ; write(unit=vtk(rf)%u,iostat=E_IO)repeat(' ',vtk(rf)%indent)//'</PointData>'//end_rec
      endselect
    endselect
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_DAT_XML

  !> @ingroup Lib_VTK_IOPrivateProcedure
  !> @{

  !> Function for saving field of scalar variable (R8P).
  !> @return E_IO: integer(I4P) error flag
  function VTK_VAR_XML_SCAL_R8(cf,NC_NN,varname,var) result(E_IO)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN), optional:: cf           !< Current file index (for concurrent files IO).
  integer(I4P), intent(IN)::           NC_NN        !< Number of cells or nodes.
  character(*), intent(IN)::           varname      !< Variable name.
  real(R8P),    intent(IN)::           var(1:NC_NN) !< Variable to be saved.
  integer(I4P)::                       E_IO         !< Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done.
  character(len=maxlen)::              s_buffer     !< Buffer string.
  integer(I4P)::                       rf           !< Real file index.
  integer(I4P)::                       n1           !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  rf = f
  if (present(cf)) then
    rf = cf ; f = cf
  endif
  select case(vtk(rf)%f)
  case(ascii)
    s_buffer = repeat(' ',vtk(rf)%indent)//'<DataArray type="Float64" Name="'//trim(varname)//&
               '" NumberOfComponents="1" format="ascii">'
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)trim(s_buffer)
    write(unit=vtk(rf)%u,fmt=FR8P,iostat=E_IO)(var(n1),n1=1,NC_NN)
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)repeat(' ',vtk(rf)%indent)//'</DataArray>'
  case(binary)
    s_buffer = repeat(' ',vtk(rf)%indent)//'<DataArray type="Float64" Name="'//trim(varname)// &
               '" NumberOfComponents="1" format="appended" offset="'//trim(str(.true.,vtk(rf)%ioffset))//'"/>'
    write(unit=vtk(rf)%u,iostat=E_IO)trim(s_buffer)//end_rec
    call vtk(rf)%byte_update(N_Byte = NC_NN*BYR8P)
    write(unit=vtk(rf)%ua,iostat=E_IO)vtk(rf)%N_Byte,'R8',NC_NN
    write(unit=vtk(rf)%ua,iostat=E_IO)(var(n1),n1=1,NC_NN)
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_VAR_XML_SCAL_R8

  !> Function for saving field of scalar variable (R4P).
  !> @return E_IO: integer(I4P) error flag
  function VTK_VAR_XML_SCAL_R4(cf,NC_NN,varname,var) result(E_IO)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN), optional:: cf           !< Current file index (for concurrent files IO).
  integer(I4P), intent(IN)::           NC_NN        !< Number of cells or nodes.
  character(*), intent(IN)::           varname      !< Variable name.
  real(R4P),    intent(IN)::           var(1:NC_NN) !< Variable to be saved.
  integer(I4P)::                       E_IO         !< Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done.
  character(len=maxlen)::              s_buffer     !< Buffer string.
  integer(I4P)::                       rf           !< Real file index.
  integer(I4P)::                       n1           !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  rf = f
  if (present(cf)) then
    rf = cf ; f = cf
  endif
  select case(vtk(rf)%f)
  case(ascii)
    s_buffer = repeat(' ',vtk(rf)%indent)//'<DataArray type="Float32" Name="'//trim(varname)//&
               '" NumberOfComponents="1" format="ascii">'
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)trim(s_buffer)
    write(unit=vtk(rf)%u,fmt=FR4P,iostat=E_IO)(var(n1),n1=1,NC_NN)
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)repeat(' ',vtk(rf)%indent)//'</DataArray>'
  case(binary)
    s_buffer = repeat(' ',vtk(rf)%indent)//'<DataArray type="Float32" Name="'//trim(varname)// &
               '" NumberOfComponents="1" format="appended" offset="'//trim(str(.true.,vtk(rf)%ioffset))//'"/>'
    write(unit=vtk(rf)%u,iostat=E_IO)trim(s_buffer)//end_rec
    call vtk(rf)%byte_update(N_Byte = NC_NN*BYR4P)
    write(unit=vtk(rf)%ua,iostat=E_IO)vtk(rf)%N_Byte,'R4',NC_NN
    write(unit=vtk(rf)%ua,iostat=E_IO)(var(n1),n1=1,NC_NN)
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_VAR_XML_SCAL_R4

  !> Function for saving field of scalar variable (I8P).
  !> @return E_IO: integer(I4P) error flag
  function VTK_VAR_XML_SCAL_I8(cf,NC_NN,varname,var) result(E_IO)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN), optional:: cf           !< Current file index (for concurrent files IO).
  integer(I4P), intent(IN)::           NC_NN        !< Number of cells or nodes.
  character(*), intent(IN)::           varname      !< Variable name.
  integer(I8P), intent(IN)::           var(1:NC_NN) !< Variable to be saved.
  integer(I4P)::                       E_IO         !< Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done.
  character(len=maxlen)::              s_buffer     !< Buffer string.
  integer(I4P)::                       rf           !< Real file index.
  integer(I4P)::                       n1           !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  rf = f
  if (present(cf)) then
    rf = cf ; f = cf
  endif
  select case(vtk(rf)%f)
  case(ascii)
    s_buffer = repeat(' ',vtk(rf)%indent)//'<DataArray type="Int64" Name="'//trim(varname)//&
               '" NumberOfComponents="1" format="ascii">'
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)trim(s_buffer)
    write(unit=vtk(rf)%u,fmt=FI8P,iostat=E_IO)(var(n1),n1=1,NC_NN)
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)'</DataArray>'
  case(binary)
    s_buffer = repeat(' ',vtk(rf)%indent)//'<DataArray type="Int64" Name="'//trim(varname)// &
               '" NumberOfComponents="1" format="appended" offset="'//trim(str(.true.,vtk(rf)%ioffset))//'"/>'
    write(unit=vtk(rf)%u,iostat=E_IO)trim(s_buffer)//end_rec
    call vtk(rf)%byte_update(N_Byte = int(NC_NN*BYI8P,I4P))
    write(unit=vtk(rf)%ua,iostat=E_IO)vtk(rf)%N_Byte,'I8',NC_NN
    write(unit=vtk(rf)%ua,iostat=E_IO)(var(n1),n1=1,NC_NN)
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_VAR_XML_SCAL_I8

  !> Function for saving field of scalar variable (I4P).
  !> @return E_IO: integer(I4P) error flag
  function VTK_VAR_XML_SCAL_I4(cf,NC_NN,varname,var) result(E_IO)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN), optional:: cf           !< Current file index (for concurrent files IO).
  integer(I4P), intent(IN)::           NC_NN        !< Number of cells or nodes.
  character(*), intent(IN)::           varname      !< Variable name.
  integer(I4P), intent(IN)::           var(1:NC_NN) !< Variable to be saved.
  integer(I4P)::                       E_IO         !< Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done.
  character(len=maxlen)::              s_buffer     !< Buffer string.
  integer(I4P)::                       rf           !< Real file index.
  integer(I4P)::                       n1           !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  rf = f
  if (present(cf)) then
    rf = cf ; f = cf
  endif
  select case(vtk(rf)%f)
  case(ascii)
    s_buffer = repeat(' ',vtk(rf)%indent)//'<DataArray type="Int32" Name="'//trim(varname)//&
               '" NumberOfComponents="1" format="ascii">'
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)trim(s_buffer)
    write(unit=vtk(rf)%u,fmt=FI4P,iostat=E_IO)(var(n1),n1=1,NC_NN)
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)repeat(' ',vtk(rf)%indent)//'</DataArray>'
  case(binary)
    s_buffer = repeat(' ',vtk(rf)%indent)//'<DataArray type="Int32" Name="'//trim(varname)// &
               '" NumberOfComponents="1" format="appended" offset="'//trim(str(.true.,vtk(rf)%ioffset))//'"/>'
    write(unit=vtk(rf)%u,iostat=E_IO)trim(s_buffer)//end_rec
    call vtk(rf)%byte_update(N_Byte = NC_NN*BYI4P)
    write(unit=vtk(rf)%ua,iostat=E_IO)vtk(rf)%N_Byte,'I4',NC_NN
    write(unit=vtk(rf)%ua,iostat=E_IO)(var(n1),n1=1,NC_NN)
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_VAR_XML_SCAL_I4

  !> Function for saving field of scalar variable (I2P).
  !> @return E_IO: integer(I4P) error flag
  function VTK_VAR_XML_SCAL_I2(cf,NC_NN,varname,var) result(E_IO)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN), optional:: cf           !< Current file index (for concurrent files IO).
  integer(I4P), intent(IN)::           NC_NN        !< Number of cells or nodes.
  character(*), intent(IN)::           varname      !< Variable name.
  integer(I2P), intent(IN)::           var(1:NC_NN) !< Variable to be saved.
  integer(I4P)::                       E_IO         !< Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done.
  character(len=maxlen)::              s_buffer     !< Buffer string.
  integer(I4P)::                       rf           !< Real file index.
  integer(I4P)::                       n1           !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  rf = f
  if (present(cf)) then
    rf = cf ; f = cf
  endif
  select case(vtk(rf)%f)
  case(ascii)
    s_buffer = repeat(' ',vtk(rf)%indent)//'<DataArray type="Int16" Name="'//trim(varname)//&
               '" NumberOfComponents="1" format="ascii">'
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)trim(s_buffer)
    write(unit=vtk(rf)%u,fmt=FI2P, iostat=E_IO)(var(n1),n1=1,NC_NN)
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)repeat(' ',vtk(rf)%indent)//'</DataArray>'
  case(binary)
    s_buffer = repeat(' ',vtk(rf)%indent)//'<DataArray type="Int16" Name="'//trim(varname)// &
               '" NumberOfComponents="1" format="appended" offset="'//trim(str(.true.,vtk(rf)%ioffset))//'"/>'
    write(unit=vtk(rf)%u,iostat=E_IO)trim(s_buffer)//end_rec
    call vtk(rf)%byte_update(N_Byte = NC_NN*BYI2P)
    write(unit=vtk(rf)%ua,iostat=E_IO)vtk(rf)%N_Byte,'I2',NC_NN
    write(unit=vtk(rf)%ua,iostat=E_IO)(var(n1),n1=1,NC_NN)
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_VAR_XML_SCAL_I2

  !> Function for saving field of scalar variable (I1P).
  !> @return E_IO: integer(I4P) error flag
  function VTK_VAR_XML_SCAL_I1(cf,NC_NN,varname,var) result(E_IO)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN), optional:: cf           !< Current file index (for concurrent files IO).
  integer(I4P), intent(IN)::           NC_NN        !< Number of cells or nodes.
  character(*), intent(IN)::           varname      !< Variable name.
  integer(I1P), intent(IN)::           var(1:NC_NN) !< Variable to be saved.
  integer(I4P)::                       E_IO         !< Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done.
  character(len=maxlen)::              s_buffer     !< Buffer string.
  integer(I4P)::                       rf           !< Real file index.
  integer(I4P)::                       n1           !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  rf = f
  if (present(cf)) then
    rf = cf ; f = cf
  endif
  select case(vtk(rf)%f)
  case(ascii)
    s_buffer=repeat(' ',vtk(rf)%indent)//'<DataArray type="Int8" Name="'//trim(varname)//'" NumberOfComponents="1" format="ascii">'
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)trim(s_buffer)
    write(unit=vtk(rf)%u,fmt=FI1P, iostat=E_IO)(var(n1),n1=1,NC_NN)
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)repeat(' ',vtk(rf)%indent)//'</DataArray>'
  case(binary)
    s_buffer=repeat(' ',vtk(rf)%indent)//'<DataArray type="Int8" Name="'//trim(varname)// &
             '" NumberOfComponents="1" format="appended" offset="'//trim(str(.true.,vtk(rf)%ioffset))//'"/>'
    write(unit=vtk(rf)%u,iostat=E_IO)trim(s_buffer)//end_rec
    call vtk(rf)%byte_update(N_Byte = NC_NN*BYI1P)
    write(unit=vtk(rf)%ua,iostat=E_IO)vtk(rf)%N_Byte,'I1',NC_NN
    write(unit=vtk(rf)%ua,iostat=E_IO)(var(n1),n1=1,NC_NN)
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_VAR_XML_SCAL_I1

  !> Function for saving field of vectorial variable (R8P).
  !> @return E_IO: integer(I4P) error flag
  function VTK_VAR_XML_VECT_R8(cf,NC_NN,varname,varX,varY,varZ) result(E_IO)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN), optional:: cf            !< Current file index (for concurrent files IO).
  integer(I4P), intent(IN)::           NC_NN         !< Number of cells or nodes.
  character(*), intent(IN)::           varname       !< Variable name.
  real(R8P),    intent(IN)::           varX(1:NC_NN) !< X component.
  real(R8P),    intent(IN)::           varY(1:NC_NN) !< Y component.
  real(R8P),    intent(IN)::           varZ(1:NC_NN) !< Z component.
  integer(I4P)::                       E_IO          !< Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done.
  character(len=maxlen)::              s_buffer      !< Buffer string.
  integer(I4P)::                       rf            !< Real file index.
  integer(I4P)::                       n1            !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  rf = f
  if (present(cf)) then
    rf = cf ; f = cf
  endif
  select case(vtk(rf)%f)
  case(ascii)
    s_buffer = repeat(' ',vtk(rf)%indent)//'<DataArray type="Float64" Name="'//trim(varname)//&
               '" NumberOfComponents="3" format="ascii">'
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)trim(s_buffer)
    write(unit=vtk(rf)%u,fmt='(3('//FR8P//',1X))',iostat=E_IO)(varX(n1),varY(n1),varZ(n1),n1=1,NC_NN)
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)repeat(' ',vtk(rf)%indent)//'</DataArray>'
  case(binary)
    s_buffer = repeat(' ',vtk(rf)%indent)//'<DataArray type="Float64" Name="'//trim(varname)// &
               '" NumberOfComponents="3" format="appended" offset="'//trim(str(.true.,vtk(rf)%ioffset))//'"/>'
    write(unit=vtk(rf)%u,iostat=E_IO)trim(s_buffer)//end_rec
    call vtk(rf)%byte_update(N_Byte = 3*NC_NN*BYR8P)
    write(unit=vtk(rf)%ua,iostat=E_IO)vtk(rf)%N_Byte,'R8',3*NC_NN
    write(unit=vtk(rf)%ua,iostat=E_IO)(varX(n1),varY(n1),varZ(n1),n1=1,NC_NN)
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_VAR_XML_VECT_R8

  !> Function for saving field of vectorial variable (R4P).
  !> @return E_IO: integer(I4P) error flag
  function VTK_VAR_XML_VECT_R4(cf,NC_NN,varname,varX,varY,varZ) result(E_IO)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN), optional:: cf            !< Current file index (for concurrent files IO).
  integer(I4P), intent(IN)::           NC_NN         !< Number of cells or nodes.
  character(*), intent(IN)::           varname       !< Variable name.
  real(R4P),    intent(IN)::           varX(1:NC_NN) !< X component.
  real(R4P),    intent(IN)::           varY(1:NC_NN) !< Y component.
  real(R4P),    intent(IN)::           varZ(1:NC_NN) !< Z component.
  integer(I4P)::                       E_IO          !< Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done.
  character(len=maxlen)::              s_buffer      !< Buffer string.
  integer(I4P)::                       rf            !< Real file index.
  integer(I4P)::                       n1            !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  rf = f
  if (present(cf)) then
    rf = cf ; f = cf
  endif
  select case(vtk(rf)%f)
  case(ascii)
    s_buffer = repeat(' ',vtk(rf)%indent)//'<DataArray type="Float32" Name="'//trim(varname)//&
               '" NumberOfComponents="3" format="ascii">'
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)trim(s_buffer)
    write(unit=vtk(rf)%u,fmt='(3('//FR4P//',1X))',iostat=E_IO)(varX(n1),varY(n1),varZ(n1),n1=1,NC_NN)
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)repeat(' ',vtk(rf)%indent)//'</DataArray>'
  case(binary)
    s_buffer = repeat(' ',vtk(rf)%indent)//'<DataArray type="Float32" Name="'//trim(varname)// &
               '" NumberOfComponents="3" format="appended" offset="'//trim(str(.true.,vtk(rf)%ioffset))//'"/>'
    write(unit=vtk(rf)%u,iostat=E_IO)trim(s_buffer)//end_rec
    call vtk(rf)%byte_update(N_Byte = 3*NC_NN*BYR4P)
    write(unit=vtk(rf)%ua,iostat=E_IO)vtk(rf)%N_Byte,'R4',3*NC_NN
    write(unit=vtk(rf)%ua,iostat=E_IO)(varX(n1),varY(n1),varZ(n1),n1=1,NC_NN)
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_VAR_XML_VECT_R4

  !> Function for saving field of vectorial variable (I8P).
  !> @return E_IO: integer(I4P) error flag
  function VTK_VAR_XML_VECT_I8(cf,NC_NN,varname,varX,varY,varZ) result(E_IO)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN), optional:: cf            !< Current file index (for concurrent files IO).
  integer(I4P), intent(IN)::           NC_NN         !< Number of cells or nodes.
  character(*), intent(IN)::           varname       !< Variable name.
  integer(I8P), intent(IN)::           varX(1:NC_NN) !< X component.
  integer(I8P), intent(IN)::           varY(1:NC_NN) !< Y component.
  integer(I8P), intent(IN)::           varZ(1:NC_NN) !< Z component.
  integer(I4P)::                       E_IO          !< Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done.
  character(len=maxlen)::              s_buffer      !< Buffer string.
  integer(I4P)::                       rf            !< Real file index.
  integer(I4P)::                       n1            !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  rf = f
  if (present(cf)) then
    rf = cf ; f = cf
  endif
  select case(vtk(rf)%f)
  case(ascii)
    s_buffer = repeat(' ',vtk(rf)%indent)//'<DataArray type="Int64" Name="'//trim(varname)//&
               '" NumberOfComponents="3" format="ascii">'
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)trim(s_buffer)
    write(unit=vtk(rf)%u,fmt='(3('//FI8P//',1X))',iostat=E_IO)(varX(n1),varY(n1),varZ(n1),n1=1,NC_NN)
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)repeat(' ',vtk(rf)%indent)//'</DataArray>'
  case(binary)
    s_buffer = repeat(' ',vtk(rf)%indent)//'<DataArray type="Int64" Name="'//trim(varname)// &
               '" NumberOfComponents="3" format="appended" offset="'//trim(str(.true.,vtk(rf)%ioffset))//'"/>'
    write(unit=vtk(rf)%u,iostat=E_IO)trim(s_buffer)//end_rec
    call vtk(rf)%byte_update(N_Byte = int(3*NC_NN*BYI8P,I4P))
    write(unit=vtk(rf)%ua,iostat=E_IO)vtk(rf)%N_Byte,'I8',3*NC_NN
    write(unit=vtk(rf)%ua,iostat=E_IO)(varX(n1),varY(n1),varZ(n1),n1=1,NC_NN)
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_VAR_XML_VECT_I8

  !> Function for saving field of vectorial variable (I4P).
  !> @return E_IO: integer(I4P) error flag
  function VTK_VAR_XML_VECT_I4(cf,NC_NN,varname,varX,varY,varZ) result(E_IO)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN), optional:: cf            !< Current file index (for concurrent files IO).
  integer(I4P), intent(IN)::           NC_NN         !< Number of cells or nodes.
  character(*), intent(IN)::           varname       !< Variable name.
  integer(I4P), intent(IN)::           varX(1:NC_NN) !< X component.
  integer(I4P), intent(IN)::           varY(1:NC_NN) !< Y component.
  integer(I4P), intent(IN)::           varZ(1:NC_NN) !< Z component.
  integer(I4P)::                       E_IO          !< Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done.
  character(len=maxlen)::              s_buffer      !< Buffer string.
  integer(I4P)::                       rf            !< Real file index.
  integer(I4P)::                       n1            !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  rf = f
  if (present(cf)) then
    rf = cf ; f = cf
  endif
  select case(vtk(rf)%f)
  case(ascii)
    s_buffer = repeat(' ',vtk(rf)%indent)//'<DataArray type="Int32" Name="'//trim(varname)//&
               '" NumberOfComponents="3" format="ascii">'
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)trim(s_buffer)
    write(unit=vtk(rf)%u,fmt='(3('//FI4P//',1X))',iostat=E_IO)(varX(n1),varY(n1),varZ(n1),n1=1,NC_NN)
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)repeat(' ',vtk(rf)%indent)//'</DataArray>'
  case(binary)
    s_buffer = repeat(' ',vtk(rf)%indent)//'<DataArray type="Int32" Name="'//trim(varname)// &
               '" NumberOfComponents="3" format="appended" offset="'//trim(str(.true.,vtk(rf)%ioffset))//'"/>'
    write(unit=vtk(rf)%u,iostat=E_IO)trim(s_buffer)//end_rec
    call vtk(rf)%byte_update(N_Byte = 3*NC_NN*BYI4P)
    write(unit=vtk(rf)%ua,iostat=E_IO)vtk(rf)%N_Byte,'I4',3*NC_NN
    write(unit=vtk(rf)%ua,iostat=E_IO)(varX(n1),varY(n1),varZ(n1),n1=1,NC_NN)
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_VAR_XML_VECT_I4

  !> Function for saving field of vectorial variable (I2P).
  !> @return E_IO: integer(I4P) error flag
  function VTK_VAR_XML_VECT_I2(cf,NC_NN,varname,varX,varY,varZ) result(E_IO)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN), optional:: cf            !< Current file index (for concurrent files IO).
  integer(I4P), intent(IN)::           NC_NN         !< Number of cells or nodes.
  character(*), intent(IN)::           varname       !< Variable name.
  integer(I2P), intent(IN)::           varX(1:NC_NN) !< X component.
  integer(I2P), intent(IN)::           varY(1:NC_NN) !< Y component.
  integer(I2P), intent(IN)::           varZ(1:NC_NN) !< Z component.
  integer(I4P)::                       E_IO          !< Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done.
  character(len=maxlen)::              s_buffer      !< Buffer string.
  integer(I4P)::                       rf            !< Real file index.
  integer(I4P)::                       n1            !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  rf = f
  if (present(cf)) then
    rf = cf ; f = cf
  endif
  select case(vtk(rf)%f)
  case(ascii)
    s_buffer = repeat(' ',vtk(rf)%indent)//'<DataArray type="Int16" Name="'//trim(varname)//&
               '" NumberOfComponents="3" format="ascii">'
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)trim(s_buffer)
    write(unit=vtk(rf)%u,fmt='(3('//FI2P//',1X))',iostat=E_IO)(varX(n1),varY(n1),varZ(n1),n1=1,NC_NN)
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)repeat(' ',vtk(rf)%indent)//'</DataArray>'
  case(binary)
    s_buffer = repeat(' ',vtk(rf)%indent)//'<DataArray type="Int16" Name="'//trim(varname)// &
               '" NumberOfComponents="3" format="appended" offset="'//trim(str(.true.,vtk(rf)%ioffset))//'"/>'
    write(unit=vtk(rf)%u,iostat=E_IO)trim(s_buffer)//end_rec
    call vtk(rf)%byte_update(N_Byte = 3*NC_NN*BYI2P)
    write(unit=vtk(rf)%ua,iostat=E_IO)vtk(rf)%N_Byte,'I2',3*NC_NN
    write(unit=vtk(rf)%ua,iostat=E_IO)(varX(n1),varY(n1),varZ(n1),n1=1,NC_NN)
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_VAR_XML_VECT_I2

  !> Function for saving field of vectorial variable (I1P).
  !> @return E_IO: integer(I4P) error flag
  function VTK_VAR_XML_VECT_I1(cf,NC_NN,varname,varX,varY,varZ) result(E_IO)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN), optional:: cf            !< Current file index (for concurrent files IO).
  integer(I4P), intent(IN)::           NC_NN         !< Number of cells or nodes.
  character(*), intent(IN)::           varname       !< Variable name.
  integer(I1P), intent(IN)::           varX(1:NC_NN) !< X component.
  integer(I1P), intent(IN)::           varY(1:NC_NN) !< Y component.
  integer(I1P), intent(IN)::           varZ(1:NC_NN) !< Z component.
  integer(I4P)::                       E_IO          !< Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done.
  character(len=maxlen)::              s_buffer      !< Buffer string.
  integer(I4P)::                       rf            !< Real file index.
  integer(I4P)::                       n1            !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  rf = f
  if (present(cf)) then
    rf = cf ; f = cf
  endif
  select case(vtk(rf)%f)
  case(ascii)
    s_buffer=repeat(' ',vtk(rf)%indent)//'<DataArray type="Int8" Name="'//trim(varname)//'" NumberOfComponents="3" format="ascii">'
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)trim(s_buffer)
    write(unit=vtk(rf)%u,fmt='(3('//FI1P//',1X))',iostat=E_IO)(varX(n1),varY(n1),varZ(n1),n1=1,NC_NN)
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)repeat(' ',vtk(rf)%indent)//'</DataArray>'
  case(binary)
    s_buffer=repeat(' ',vtk(rf)%indent)//'<DataArray type="Int8" Name="'//trim(varname)// &
             '" NumberOfComponents="3" format="appended" offset="'//trim(str(.true.,vtk(rf)%ioffset))//'"/>'
    write(unit=vtk(rf)%u,iostat=E_IO)trim(s_buffer)//end_rec
    vtk(rf)%N_Byte = 3*NC_NN*BYI1P
    call vtk(rf)%byte_update(N_Byte = 3*NC_NN*BYI1P)
    write(unit=vtk(rf)%ua,iostat=E_IO)vtk(rf)%N_Byte,'I1',3*NC_NN
    write(unit=vtk(rf)%ua,iostat=E_IO)(varX(n1),varY(n1),varZ(n1),n1=1,NC_NN)
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_VAR_XML_VECT_I1

  !> Function for saving field of list variable (R8P).
  !> @return E_IO: integer(I4P) error flag
  function VTK_VAR_XML_LIST_R8(cf,NC_NN,N_COL,varname,var) result(E_IO)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN), optional:: cf         !< Current file index (for concurrent files IO).
  integer(I4P), intent(IN)::           NC_NN      !< Number of cells or nodes.
  integer(I4P), intent(IN)::           N_COL      !< Number of columns.
  character(*), intent(IN)::           varname    !< Variable name.
  real(R8P),    intent(IN)::           var(1:,1:) !< Components.
  integer(I4P)::                       E_IO       !< Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done.
  character(len=maxlen)::              s_buffer   !< Buffer string.
  integer(I4P)::                       rf         !< Real file index.
  integer(I4P)::                       n1,n2      !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  rf = f
  if (present(cf)) then
    rf = cf ; f = cf
  endif
  select case(vtk(rf)%f)
  case(ascii)
    s_buffer = repeat(' ',vtk(rf)%indent)//'<DataArray type="Float64" Name="'//trim(varname)//'" NumberOfComponents="'// &
               trim(str(.true.,N_COL))//'" format="ascii">'
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)trim(s_buffer)
    do n1=1,NC_NN
      write(unit=vtk(rf)%u,fmt=FR8P,iostat=E_IO)(var(n1,n2),n2=1,N_COL)
    enddo
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)repeat(' ',vtk(rf)%indent)//'</DataArray>'
  case(binary)
    s_buffer = repeat(' ',vtk(rf)%indent)//'<DataArray type="Float64" Name="'//trim(varname)//'" NumberOfComponents="'// &
               trim(str(.true.,N_COL))//'" format="appended" offset="'//trim(str(.true.,vtk(rf)%ioffset))//'"/>'
    write(unit=vtk(rf)%u,iostat=E_IO)trim(s_buffer)//end_rec
    call vtk(rf)%byte_update(N_Byte = N_COL*NC_NN*BYR8P)
    write(unit=vtk(rf)%ua,iostat=E_IO)vtk(rf)%N_Byte,'R8',N_COL*NC_NN
    do n1=1,NC_NN
      write(unit=vtk(rf)%ua,iostat=E_IO)var(n1,:)
    enddo
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_VAR_XML_LIST_R8

  !> Function for saving field of list variable (R4P).
  !> @return E_IO: integer(I4P) error flag
  function VTK_VAR_XML_LIST_R4(cf,NC_NN,N_COL,varname,var) result(E_IO)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN), optional:: cf         !< Current file index (for concurrent files IO).
  integer(I4P), intent(IN)::           NC_NN      !< Number of cells or nodes.
  integer(I4P), intent(IN)::           N_COL      !< Number of columns.
  character(*), intent(IN)::           varname    !< Variable name.
  real(R4P),    intent(IN)::           var(1:,1:) !< Components.
  integer(I4P)::                       E_IO       !< Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done.
  character(len=maxlen)::              s_buffer   !< Buffer string.
  integer(I4P)::                       rf         !< Real file index.
  integer(I4P)::                       n1,n2      !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  rf = f
  if (present(cf)) then
    rf = cf ; f = cf
  endif
  select case(vtk(rf)%f)
  case(ascii)
    s_buffer = repeat(' ',vtk(rf)%indent)//'<DataArray type="Float32" Name="'//trim(varname)//'" NumberOfComponents="'// &
               trim(str(.true.,N_COL))//'" format="ascii">'
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)trim(s_buffer)
    do n1=1,NC_NN
      write(unit=vtk(rf)%u,fmt=FR4P,iostat=E_IO)(var(n1,n2),n2=1,N_COL)
    enddo
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)repeat(' ',vtk(rf)%indent)//'</DataArray>'
  case(binary)
    s_buffer = repeat(' ',vtk(rf)%indent)//'<DataArray type="Float32" Name="'//trim(varname)//'" NumberOfComponents="'// &
               trim(str(.true.,N_COL))//'" format="appended" offset="'//trim(str(.true.,vtk(rf)%ioffset))//'"/>'
    write(unit=vtk(rf)%u,iostat=E_IO)trim(s_buffer)//end_rec
    call vtk(rf)%byte_update(N_Byte = N_COL*NC_NN*BYR4P)
    write(unit=vtk(rf)%ua,iostat=E_IO)vtk(rf)%N_Byte,'R4',N_COL*NC_NN
    do n1=1,NC_NN
      write(unit=vtk(rf)%ua,iostat=E_IO)var(n1,:)
    enddo
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_VAR_XML_LIST_R4

  !> Function for saving field of list variable (I8P).
  !> @return E_IO: integer(I4P) error flag
  function VTK_VAR_XML_LIST_I8(cf,NC_NN,N_COL,varname,var) result(E_IO)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN), optional:: cf         !< Current file index (for concurrent files IO).
  integer(I4P), intent(IN)::           NC_NN      !< Number of cells or nodes.
  integer(I4P), intent(IN)::           N_COL      !< Number of columns.
  character(*), intent(IN)::           varname    !< Variable name.
  integer(I8P), intent(IN)::           var(1:,1:) !< Components.
  integer(I4P)::                       E_IO       !< Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done.
  character(len=maxlen)::              s_buffer   !< Buffer string.
  integer(I4P)::                       rf         !< Real file index.
  integer(I4P)::                       n1,n2      !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  rf = f
  if (present(cf)) then
    rf = cf ; f = cf
  endif
  select case(vtk(rf)%f)
  case(ascii)
    s_buffer = repeat(' ',vtk(rf)%indent)//'<DataArray type="Int64" Name="'//trim(varname)//'" NumberOfComponents="'// &
               trim(str(.true.,N_COL))//'" format="ascii">'
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)trim(s_buffer)
    do n1=1,NC_NN
      write(unit=vtk(rf)%u,fmt=FI8P,iostat=E_IO)(var(n1,n2),n2=1,N_COL)
    enddo
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)repeat(' ',vtk(rf)%indent)//'</DataArray>'
  case(binary)
    s_buffer = repeat(' ',vtk(rf)%indent)//'<DataArray type="Int64" Name="'//trim(varname)//'" NumberOfComponents="'// &
               trim(str(.true.,N_COL))//'" format="appended" offset="'//trim(str(.true.,vtk(rf)%ioffset))//'"/>'
    write(unit=vtk(rf)%u,iostat=E_IO)trim(s_buffer)//end_rec
    call vtk(rf)%byte_update(N_Byte = int(N_COL*NC_NN*BYI8P,I4P))
    write(unit=vtk(rf)%ua,iostat=E_IO)vtk(rf)%N_Byte,'I8',N_COL*NC_NN
    do n1=1,NC_NN
      write(unit=vtk(rf)%ua,iostat=E_IO)var(n1,:)
    enddo
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_VAR_XML_LIST_I8

  !> Function for saving field of list variable (I4P).
  !> @return E_IO: integer(I4P) error flag
  function VTK_VAR_XML_LIST_I4(cf,NC_NN,N_COL,varname,var) result(E_IO)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN), optional:: cf         !< Current file index (for concurrent files IO).
  integer(I4P), intent(IN)::           NC_NN      !< Number of cells or nodes.
  integer(I4P), intent(IN)::           N_COL      !< Number of columns.
  character(*), intent(IN)::           varname    !< Variable name.
  integer(I4P), intent(IN)::           var(1:,1:) !< Components.
  integer(I4P)::                       E_IO       !< Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done.
  character(len=maxlen)::              s_buffer   !< Buffer string.
  integer(I4P)::                       rf         !< Real file index.
  integer(I4P)::                       n1,n2      !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  rf = f
  if (present(cf)) then
    rf = cf ; f = cf
  endif
  select case(vtk(rf)%f)
  case(ascii)
    s_buffer = repeat(' ',vtk(rf)%indent)//'<DataArray type="Int32" Name="'//trim(varname)//'" NumberOfComponents="'// &
               trim(str(.true.,N_COL))//'" format="ascii">'
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)trim(s_buffer)
    do n1=1,NC_NN
      write(unit=vtk(rf)%u,fmt=FI4P,iostat=E_IO)(var(n1,n2),n2=1,N_COL)
    enddo
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)repeat(' ',vtk(rf)%indent)//'</DataArray>'
  case(binary)
    s_buffer = repeat(' ',vtk(rf)%indent)//'<DataArray type="Int32" Name="'//trim(varname)//'" NumberOfComponents="'// &
               trim(str(.true.,N_COL))//'" format="appended" offset="'//trim(str(.true.,vtk(rf)%ioffset))//'"/>'
    write(unit=vtk(rf)%u,iostat=E_IO)trim(s_buffer)//end_rec
    call vtk(rf)%byte_update(N_Byte = N_COL*NC_NN*BYI4P)
    write(unit=vtk(rf)%ua,iostat=E_IO)vtk(rf)%N_Byte,'I4',N_COL*NC_NN
    do n1=1,NC_NN
      write(unit=vtk(rf)%ua,iostat=E_IO)var(n1,:)
    enddo
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_VAR_XML_LIST_I4

  !> Function for saving field of list variable (I2P).
  !> @return E_IO: integer(I4P) error flag
  function VTK_VAR_XML_LIST_I2(cf,NC_NN,N_COL,varname,var) result(E_IO)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN), optional:: cf         !< Current file index (for concurrent files IO).
  integer(I4P), intent(IN)::           NC_NN      !< Number of cells or nodes.
  integer(I4P), intent(IN)::           N_COL      !< Number of columns.
  character(*), intent(IN)::           varname    !< Variable name.
  integer(I2P), intent(IN)::           var(1:,1:) !< Components.
  integer(I4P)::                       E_IO       !< Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done.
  character(len=maxlen)::              s_buffer   !< Buffer string.
  integer(I4P)::                       rf         !< Real file index.
  integer(I4P)::                       n1,n2      !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  rf = f
  if (present(cf)) then
    rf = cf ; f = cf
  endif
  select case(vtk(rf)%f)
  case(ascii)
    s_buffer = repeat(' ',vtk(rf)%indent)//'<DataArray type="Int16" Name="'//trim(varname)//'" NumberOfComponents="'// &
               trim(str(.true.,N_COL))//'" format="ascii">'
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)trim(s_buffer)
    do n1=1,NC_NN
      write(unit=vtk(rf)%u,fmt=FI2P,iostat=E_IO)(var(n1,n2),n2=1,N_COL)
    enddo
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)repeat(' ',vtk(rf)%indent)//'</DataArray>'
  case(binary)
    s_buffer = repeat(' ',vtk(rf)%indent)//'<DataArray type="Int16" Name="'//trim(varname)//'" NumberOfComponents="'// &
               trim(str(.true.,N_COL))//'" format="appended" offset="'//trim(str(.true.,vtk(rf)%ioffset))//'"/>'
    write(unit=vtk(rf)%u,iostat=E_IO)trim(s_buffer)//end_rec
    call vtk(rf)%byte_update(N_Byte = N_COL*NC_NN*BYI2P)
    write(unit=vtk(rf)%ua,iostat=E_IO)vtk(rf)%N_Byte,'I2',N_COL*NC_NN
    do n1=1,NC_NN
      write(unit=vtk(rf)%ua,iostat=E_IO)var(n1,:)
    enddo
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_VAR_XML_LIST_I2

  !> Function for saving field of list variable (I1P).
  !> @return E_IO: integer(I4P) error flag
  function VTK_VAR_XML_LIST_I1(cf,NC_NN,N_COL,varname,var) result(E_IO)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN), optional:: cf         !< Current file index (for concurrent files IO).
  integer(I4P), intent(IN)::           NC_NN      !< Number of cells or nodes.
  integer(I4P), intent(IN)::           N_COL      !< Number of columns.
  character(*), intent(IN)::           varname    !< Variable name.
  integer(I1P), intent(IN)::           var(1:,1:) !< Components.
  integer(I4P)::                       E_IO       !< Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done.
  character(len=maxlen)::              s_buffer   !< Buffer string.
  integer(I4P)::                       rf         !< Real file index.
  integer(I4P)::                       n1,n2      !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  rf = f
  if (present(cf)) then
    rf = cf ; f = cf
  endif
  select case(vtk(rf)%f)
  case(ascii)
    s_buffer = repeat(' ',vtk(rf)%indent)//'<DataArray type="Int8" Name="'//trim(varname)//'" NumberOfComponents="'// &
               trim(str(.true.,N_COL))//'" format="ascii">'
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)trim(s_buffer)
    do n1=1,NC_NN
      write(unit=vtk(rf)%u,fmt=FI1P,iostat=E_IO)(var(n1,n2),n2=1,N_COL)
    enddo
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)repeat(' ',vtk(rf)%indent)//'</DataArray>'
  case(binary)
    s_buffer = repeat(' ',vtk(rf)%indent)//'<DataArray type="Int8" Name="'//trim(varname)//'" NumberOfComponents="'// &
               trim(str(.true.,N_COL))//'" format="appended" offset="'//trim(str(.true.,vtk(rf)%ioffset))//'"/>'
    write(unit=vtk(rf)%u,iostat=E_IO)trim(s_buffer)//end_rec
    call vtk(rf)%byte_update(N_Byte = N_COL*NC_NN*BYI1P)
    write(unit=vtk(rf)%ua,iostat=E_IO)vtk(rf)%N_Byte,'I1',N_COL*NC_NN
    do n1=1,NC_NN
      write(unit=vtk(rf)%ua,iostat=E_IO)var(n1,:)
    enddo
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_VAR_XML_LIST_I1
  !> @}

  !> @brief Function for finalizing the VTK-XML file.
  !> @note An example of usage is: \n
  !> @code ...
  !> E_IO = VTK_END_XML()
  !> ... @endcode
  !> @return E_IO: integer(I4P) error flag
  !> @ingroup Lib_VTK_IOPublicProcedure
  function VTK_END_XML(cf) result(E_IO)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(INOUT), optional:: cf       !< Current file index (for concurrent files IO).
  integer(I4P)::                          E_IO     !< Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done.
  character(2)::                          var_type !< Varable type = R8,R4,I8,I4,I2,I1.
  real(R8P),    allocatable::             v_R8(:)  !< R8 vector for IO in AppendData.
  real(R4P),    allocatable::             v_R4(:)  !< R4 vector for IO in AppendData.
  integer(I8P), allocatable::             v_I8(:)  !< I8 vector for IO in AppendData.
  integer(I4P), allocatable::             v_I4(:)  !< I4 vector for IO in AppendData.
  integer(I2P), allocatable::             v_I2(:)  !< I2 vector for IO in AppendData.
  integer(I1P), allocatable::             v_I1(:)  !< I1 vector for IO in AppendData.
  integer(I4P)::                          rf       !< Real file index.
#ifdef HUGE
  integer(I8P)::                          N_v      !< Vector dimension.
  integer(I8P)::                          n1       !< Counter.
#else
  integer(I4P)::                          N_v      !< Vector dimension.
  integer(I4P)::                          n1       !< Counter.
#endif
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  rf = f
  if (present(cf)) then
    rf = cf ; f = cf
  endif
  select case(vtk(rf)%f)
  case(ascii)
    vtk(rf)%indent = vtk(rf)%indent - 2
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)repeat(' ',vtk(rf)%indent)//'</'//trim(vtk(rf)%topology)//'>'
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)'</VTKFile>'
  case(binary)
    vtk(rf)%indent = vtk(rf)%indent - 2
    write(unit  =vtk(rf)%u, iostat=E_IO)repeat(' ',vtk(rf)%indent)//'</'//trim(vtk(rf)%topology)//'>'//end_rec
    write(unit  =vtk(rf)%u, iostat=E_IO)repeat(' ',vtk(rf)%indent)//'<AppendedData encoding="raw">'//end_rec
    write(unit  =vtk(rf)%u, iostat=E_IO)'_'
    endfile(unit=vtk(rf)%ua,iostat=E_IO)
    rewind(unit =vtk(rf)%ua,iostat=E_IO)
    do
      read(unit=vtk(rf)%ua,iostat=E_IO,end=100)vtk(rf)%N_Byte,var_type,N_v
      select case(var_type)
      case('R8')
        allocate(v_R8(1:N_v))
        read(unit =vtk(rf)%ua,iostat=E_IO)(v_R8(n1),n1=1,N_v)
        write(unit=vtk(rf)%u, iostat=E_IO)int(vtk(rf)%N_Byte,I4P),(v_R8(n1),n1=1,N_v)
        deallocate(v_R8)
      case('R4')
        allocate(v_R4(1:N_v))
        read(unit =vtk(rf)%ua,iostat=E_IO)(v_R4(n1),n1=1,N_v)
        write(unit=vtk(rf)%u, iostat=E_IO)int(vtk(rf)%N_Byte,I4P),(v_R4(n1),n1=1,N_v)
        deallocate(v_R4)
      case('I8')
        allocate(v_I8(1:N_v))
        read(unit =vtk(rf)%ua,iostat=E_IO)(v_I8(n1),n1=1,N_v)
        write(unit=vtk(rf)%u, iostat=E_IO)int(vtk(rf)%N_Byte,I4P),(v_I8(n1),n1=1,N_v)
        deallocate(v_I8)
      case('I4')
        allocate(v_I4(1:N_v))
        read(unit =vtk(rf)%ua,iostat=E_IO)(v_I4(n1),n1=1,N_v)
        write(unit=vtk(rf)%u, iostat=E_IO)int(vtk(rf)%N_Byte,I4P),(v_I4(n1),n1=1,N_v)
        deallocate(v_I4)
      case('I2')
        allocate(v_I2(1:N_v))
        read(unit =vtk(rf)%ua,iostat=E_IO)(v_I2(n1),n1=1,N_v)
        write(unit=vtk(rf)%u, iostat=E_IO)int(vtk(rf)%N_Byte,I4P),(v_I2(n1),n1=1,N_v)
        deallocate(v_I2)
      case('I1')
        allocate(v_I1(1:N_v))
        read(unit =vtk(rf)%ua,iostat=E_IO)(v_I1(n1),n1=1,N_v)
        write(unit=vtk(rf)%u, iostat=E_IO)int(vtk(rf)%N_Byte,I4P),(v_I1(n1),n1=1,N_v)
        deallocate(v_I1)
      case default
        E_IO = 1
        write (stderr,'(A)')' bad var_type = '//var_type
        write (stderr,'(A)')' N_Byte = '//trim(str(n=vtk(rf)%N_Byte))//' N_v = '//trim(str(n=N_v))
        return
      endselect
    enddo
    100 continue
    write(unit=vtk(rf)%u,iostat=E_IO)end_rec
    write(unit=vtk(rf)%u,iostat=E_IO)repeat(' ',vtk(rf)%indent)//'</AppendedData>'//end_rec
    write(unit=vtk(rf)%u,iostat=E_IO)'</VTKFile>'//end_rec
    close(unit=vtk(rf)%ua,iostat=E_IO)
  endselect
  close(unit=vtk(rf)%u,iostat=E_IO)
  call vtk_update(act='remove',cf=rf,Nvtk=Nvtk,vtk=vtk)
  f = rf
  if (present(cf)) cf = rf
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_END_XML

  !> The VTK_VTM_XML function is used for initializing a VTM (VTK Multiblocks) XML file that is a wrapper to a set of VTK-XML files.
  !> @return E_IO: integer(I4P) error flag
  !> @ingroup Lib_VTK_IOPublicProcedure
  function VTM_INI_XML(filename) result(E_IO)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  character(*), intent(IN):: filename !< File name of output VTM file.
  integer(I4P)::             E_IO     !< Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done.
  character(len=maxlen)::    s_buffer !< Buffer string.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.ir_initialized) call IR_Init
  if (endian==endianL) then
    s_buffer='<VTKFile type="vtkMultiBlockDataSet" version="1.0" byte_order="LittleEndian">'
  else
    s_buffer='<VTKFile type="vtkMultiBlockDataSet" version="1.0" byte_order="BigEndian">'
  endif
  open(unit=Get_Unit(vtm%u),file=trim(filename),form='FORMATTED',access='SEQUENTIAL',action='WRITE',status='REPLACE',iostat=E_IO)
  write(unit=vtm%u,fmt='(A)',iostat=E_IO)'<?xml version="1.0"?>'
  write(unit=vtm%u,fmt='(A)',iostat=E_IO)trim(s_buffer) ; vtm%indent = 2
  write(unit=vtm%u,fmt='(A)',iostat=E_IO)repeat(' ',vtm%indent)//'<vtkMultiBlockDataSet>' ; vtm%indent = vtm%indent + 2
  vtm%blk = -1
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction VTM_INI_XML

  !> The VTM_BLK_XML function is used for opening or closing a block level of a VTM file.
  !> @return E_IO: integer(I4P) error flag
  !> @ingroup Lib_VTK_IOPublicProcedure
  function VTM_BLK_XML(block_action) result(E_IO)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  character(*), intent(IN):: block_action !< Block action: OPEN or CLOSE block.
  integer(I4P)::             E_IO         !< Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  select case(trim(Upper_Case(block_action)))
  case('OPEN')
    vtm%blk = vtm%blk + 1
    write(unit=vtm%u,fmt='(A)',iostat=E_IO)repeat(' ',vtm%indent)//'<Block index="'//trim(str(.true.,vtm%blk))//'">'
    vtm%indent = vtm%indent + 2
  case('CLOSE')
    vtm%indent = vtm%indent - 2 ; write(unit=vtm%u,fmt='(A)',iostat=E_IO)repeat(' ',vtm%indent)//'</Block>'
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction VTM_BLK_XML

  !> The VTM_WRF_XML function is used for saving the list of VTK-XML wrapped files by the actual block of the mutliblock VTM file.
  !> @return E_IO: integer(I4P) error flag
  !> @ingroup Lib_VTK_IOPublicProcedure
  function VTM_WRF_XML(flist) result(E_IO)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  character(*), intent(IN):: flist(:) !< List of VTK-XML wrapped files.
  integer(I4P)::             E_IO     !< Input/Output inquiring flag: 0 if IO is done, > 0 if IO is not done.
  integer(I4P)::             f        !< File counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  do f=1,size(flist)
    write(unit=vtm%u,fmt='(A)',iostat=E_IO)repeat(' ',vtm%indent)//'<DataSet index="'//trim(str(.true.,f-1))//'" file="'// &
                                           adjustl(trim(flist(f)))//'"/>'
  enddo
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction VTM_WRF_XML

  !> Function for finalizing open file, it has not inputs, @libvtk manages the file unit without the
  !> user's action.
  !> @return E_IO: integer(I4P) error flag
  !> @ingroup Lib_VTK_IOPublicProcedure
  function VTM_END_XML() result(E_IO)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P):: E_IO !< Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  vtm%indent = vtm%indent - 2
  write(unit=vtm%u,fmt='(A)',iostat=E_IO)repeat(' ',vtm%indent)//'</vtkMultiBlockDataSet>'
  write(unit=vtm%u,fmt='(A)',iostat=E_IO)'</VTKFile>'
  close(unit=vtm%u)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction VTM_END_XML

  !> @brief Function for initializing parallel (partitioned) VTK-XML file.
  !> @return E_IO: integer(I4P) error flag
  !> @ingroup Lib_VTK_IOPublicProcedure
  function PVTK_INI_XML(cf,nx1,nx2,ny1,ny2,nz1,nz2,filename,mesh_topology,tp) result(E_IO)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(OUT), optional:: cf            !< Current file index (for concurrent files IO).
  integer(I4P), intent(IN),  optional:: nx1           !< Initial node of x axis.
  integer(I4P), intent(IN),  optional:: nx2           !< Final node of x axis.
  integer(I4P), intent(IN),  optional:: ny1           !< Initial node of y axis.
  integer(I4P), intent(IN),  optional:: ny2           !< Final node of y axis.
  integer(I4P), intent(IN),  optional:: nz1           !< Initial node of z axis.
  integer(I4P), intent(IN),  optional:: nz2           !< Final node of z axis.
  character(*), intent(IN)::            filename      !< File name.
  character(*), intent(IN)::            mesh_topology !< Mesh topology.
  character(*), intent(IN)::            tp            !< Type of geometry representation (Float32, Float64, ecc).
  integer(I4P)::                        E_IO          !< Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done.
  character(len=maxlen)::               s_buffer      !< Buffer string.
  integer(I4P)::                        rf            !< Real file index.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.ir_initialized) call IR_Init
  call vtk_update(act='add',cf=rf,Nvtk=Nvtk,vtk=vtk)
  f = rf
  if (present(cf)) cf = rf
  vtk(rf)%topology = trim(mesh_topology)
  open(unit=Get_Unit(vtk(rf)%u),file=trim(filename),&
       form='FORMATTED',access='SEQUENTIAL',action='WRITE',status='REPLACE',iostat=E_IO)
  write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)'<?xml version="1.0"?>'
  if (endian==endianL) then
    s_buffer = '<VTKFile type="'//trim(vtk(rf)%topology)//'" version="0.1" byte_order="LittleEndian">'
  else
    s_buffer = '<VTKFile type="'//trim(vtk(rf)%topology)//'" version="0.1" byte_order="BigEndian">'
  endif
  write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)trim(s_buffer) ; vtk(rf)%indent = 2
  select case(trim(vtk(rf)%topology))
  case('PRectilinearGrid')
    s_buffer = repeat(' ',vtk(rf)%indent)//'<'//trim(vtk(rf)%topology)//' WholeExtent="'//&
               trim(str(n=nx1))//' '//trim(str(n=nx2))//' '//                             &
               trim(str(n=ny1))//' '//trim(str(n=ny2))//' '//                             &
               trim(str(n=nz1))//' '//trim(str(n=nz2))//'" GhostLevel="#">'
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)trim(s_buffer) ; vtk(rf)%indent = vtk(rf)%indent + 2
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)repeat(' ',vtk(rf)%indent)//'<PCoordinates>' ; vtk(rf)%indent = vtk(rf)%indent + 2
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)repeat(' ',vtk(rf)%indent)//'<PDataArray type="'//trim(tp)//'"/>'
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)repeat(' ',vtk(rf)%indent)//'<PDataArray type="'//trim(tp)//'"/>'
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)repeat(' ',vtk(rf)%indent)//'<PDataArray type="'//trim(tp)//'"/>'
    vtk(rf)%indent = vtk(rf)%indent - 2 ; write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)repeat(' ',vtk(rf)%indent)//'</PCoordinates>'
  case('PStructuredGrid')
    s_buffer = repeat(' ',vtk(rf)%indent)//'<'//trim(vtk(rf)%topology)//' WholeExtent="'//&
               trim(str(n=nx1))//' '//trim(str(n=nx2))//' '//                             &
               trim(str(n=ny1))//' '//trim(str(n=ny2))//' '//                             &
               trim(str(n=nz1))//' '//trim(str(n=nz2))//'" GhostLevel="#">'
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)trim(s_buffer) ; vtk(rf)%indent = vtk(rf)%indent + 2
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)repeat(' ',vtk(rf)%indent)//'<PPoints>' ; vtk(rf)%indent = vtk(rf)%indent + 2
    s_buffer = repeat(' ',vtk(rf)%indent)//'<PDataArray type="'//trim(tp)//'" NumberOfComponents="3" Name="Points"/>'
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)trim(s_buffer)
    vtk(rf)%indent = vtk(rf)%indent - 2 ; write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)repeat(' ',vtk(rf)%indent)//'</PPoints>'
  case('PUnstructuredGrid')
    s_buffer = repeat(' ',vtk(rf)%indent)//'<'//trim(vtk(rf)%topology)//' GhostLevel="0">'
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)trim(s_buffer) ; vtk(rf)%indent = vtk(rf)%indent + 2
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)repeat(' ',vtk(rf)%indent)//'<PPoints>' ; vtk(rf)%indent = vtk(rf)%indent + 2
    s_buffer = repeat(' ',vtk(rf)%indent)//'<PDataArray type="'//trim(tp)//'" NumberOfComponents="3" Name="Points"/>'
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)trim(s_buffer)
    vtk(rf)%indent = vtk(rf)%indent - 2 ; write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)repeat(' ',vtk(rf)%indent)//'</PPoints>'
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction PVTK_INI_XML

  !> Function for saving piece geometry source for parallel (partitioned) VTK-XML file.
  !> @return E_IO: integer(I4P) error flag
  !> @ingroup Lib_VTK_IOPublicProcedure
  function PVTK_GEO_XML(cf,nx1,nx2,ny1,ny2,nz1,nz2,source) result(E_IO)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN), optional:: cf       !< Current file index (for concurrent files IO).
  integer(I4P), intent(IN), optional:: nx1      !< Initial node of x axis.
  integer(I4P), intent(IN), optional:: nx2      !< Final node of x axis.
  integer(I4P), intent(IN), optional:: ny1      !< Initial node of y axis.
  integer(I4P), intent(IN), optional:: ny2      !< Final node of y axis.
  integer(I4P), intent(IN), optional:: nz1      !< Initial node of z axis.
  integer(I4P), intent(IN), optional:: nz2      !< Final node of z axis.
  character(*), intent(IN)::           source   !< Source file name containing the piece data.
  integer(I4P)::                       E_IO     !< Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done.
  character(len=maxlen)::              s_buffer !< Buffer string.
  integer(I4P)::                       rf       !< Real file index.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  rf = f
  if (present(cf)) then
    rf = cf ; f = cf
  endif
  select case (vtk(rf)%topology)
  case('PRectilinearGrid','PStructuredGrid')
    s_buffer = repeat(' ',vtk(rf)%indent)//'<Piece Extent="'// &
               trim(str(n=nx1))//' '//trim(str(n=nx2))//' '// &
               trim(str(n=ny1))//' '//trim(str(n=ny2))//' '// &
               trim(str(n=nz1))//' '//trim(str(n=nz2))//'" Source="'//trim(source)//'"/>'
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)trim(s_buffer)
  case('PUnstructuredGrid')
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)repeat(' ',vtk(rf)%indent)//'<Piece Source="'//trim(source)//'"/>'
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction PVTK_GEO_XML

  !> Function that \b must be called before saving the data related to geometric mesh, this function initializes the
  !> saving of data variables indicating the \em type (node or cell centered) of variables that will be saved.
  !> @ingroup Lib_VTK_IOPublicProcedure
  function PVTK_DAT_XML(cf,var_location,var_block_action) result(E_IO)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN), optional:: cf               !< Current file index (for concurrent files IO).
  character(*), intent(IN)::           var_location     !< Location of saving variables: CELL or NODE centered.
  character(*), intent(IN)::           var_block_action !< Variables block action: OPEN or CLOSE block.
  integer(I4P)::                       E_IO             !< Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done.
  integer(I4P)::                       rf               !< Real file index.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  rf = f
  if (present(cf)) then
    rf = cf ; f = cf
  endif
  select case(trim(Upper_Case(var_location)))
  case('CELL')
    select case(trim(Upper_Case(var_block_action)))
    case('OPEN')
      write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)repeat(' ',vtk(rf)%indent)//'<PCellData>' ; vtk(rf)%indent = vtk(rf)%indent + 2
    case('CLOSE')
      vtk(rf)%indent = vtk(rf)%indent - 2 ; write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)repeat(' ',vtk(rf)%indent)//'</PCellData>'
    endselect
  case('NODE')
    select case(trim(Upper_Case(var_block_action)))
    case('OPEN')
      write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)repeat(' ',vtk(rf)%indent)//'<PPointData>' ; vtk(rf)%indent = vtk(rf)%indent + 2
    case('CLOSE')
      vtk(rf)%indent = vtk(rf)%indent - 2 ; write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)repeat(' ',vtk(rf)%indent)//'</PPointData>'
    endselect
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction PVTK_DAT_XML

  !> Function for saving variable associated to nodes or cells geometry.
  !> @return E_IO: integer(I4P) error flag
  !> @ingroup Lib_VTK_IOPublicProcedure
  function PVTK_VAR_XML(cf,Nc,varname,tp) result(E_IO)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN), optional:: cf       !< Current file index (for concurrent files IO).
  integer(I4P), intent(IN), optional:: Nc       !< Number of components of variable.
  character(*), intent(IN)::           varname  !< Variable name.
  character(*), intent(IN)::           tp       !< Type of data representation (Float32, Float64, ecc).
  integer(I4P)::                       E_IO     !< Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done.
  character(len=maxlen)::              s_buffer !< Buffer string.
  integer(I4P)::                       rf       !< Real file index.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  rf = f
  if (present(cf)) then
    rf = cf ; f = cf
  endif
  if (present(Nc)) then
    s_buffer = repeat(' ',vtk(rf)%indent)//'<PDataArray type="'//trim(tp)//'" Name="'//trim(varname)//&
               '" NumberOfComponents="'//trim(str(.true.,Nc))//'"/>'
  else
    s_buffer = repeat(' ',vtk(rf)%indent)//'<PDataArray type="'//trim(tp)//'" Name="'//trim(varname)//'"/>'
  endif
  write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)trim(s_buffer)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction PVTK_VAR_XML

  !> @brief Function for finalizing the parallel (partitioned) VTK-XML file.
  !> @return E_IO: integer(I4P) error flag
  !> @ingroup Lib_VTK_IOPublicProcedure
  function PVTK_END_XML(cf) result(E_IO)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(INOUT), optional:: cf   !< Current file index (for concurrent files IO).
  integer(I4P)::                          E_IO !< Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done.
  integer(I4P)::                          rf   !< Real file index.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  rf = f
  if (present(cf)) then
    rf = cf ; f = cf
  endif
  vtk(rf)%indent = vtk(rf)%indent - 2
  write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)repeat(' ',vtk(rf)%indent)//'</'//trim(vtk(rf)%topology)//'>'
  write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)'</VTKFile>'
  close(unit=vtk(rf)%u,iostat=E_IO)
  call vtk_update(act='remove',cf=rf,Nvtk=Nvtk,vtk=vtk)
  f = rf
  if (present(cf)) cf = rf
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction PVTK_END_XML

  !> @brief Function for initializing VTK-legacy file.
  !> @remark This function must be the first to be called.
  !> @note An example of usage is: \n
  !> @code ...
  !> E_IO=VTK_INI('Binary','example.vtk','VTK legacy file','UNSTRUCTURED_GRID')
  !> ... @endcode
  !> @return E_IO: integer(I4P) error flag
  !> @ingroup Lib_VTK_IOPublicProcedure
  function VTK_INI(cf,output_format,filename,title,mesh_topology) result(E_IO)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(OUT), optional:: cf            !< Current file index (for concurrent files IO).
  character(*), intent(IN)::            output_format !< Output format: ASCII or BINARY.
  character(*), intent(IN)::            filename      !< Name of file.
  character(*), intent(IN)::            title         !< Title.
  character(*), intent(IN)::            mesh_topology !< Mesh topology.
  integer(I4P)::                        E_IO          !< Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done.
  integer(I4P)::                        rf            !< Real file index.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (.not.ir_initialized) call IR_Init
  call vtk_update(act='add',cf=rf,Nvtk=Nvtk,vtk=vtk)
  f = rf
  if (present(cf)) cf = rf
  vtk(rf)%topology = trim(mesh_topology)
  select case(trim(Upper_Case(output_format)))
  case('ASCII')
    vtk(rf)%f = ascii
    open(unit=Get_Unit(vtk(rf)%u),file=trim(filename),form='FORMATTED',&
         access='SEQUENTIAL',action='WRITE',status='REPLACE',iostat=E_IO)
    ! writing header of file
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)'# vtk DataFile Version 3.0'
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)trim(title)
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)trim(Upper_Case(output_format))
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)'DATASET '//trim(vtk(rf)%topology)
  case('BINARY')
    vtk(rf)%f = binary
    open(unit=Get_Unit(vtk(rf)%u),file=trim(filename),&
         form='UNFORMATTED',access='STREAM',action='WRITE',status='REPLACE',iostat=E_IO)
    ! writing header of file
    write(unit=vtk(rf)%u,iostat=E_IO)'# vtk DataFile Version 3.0'//end_rec
    write(unit=vtk(rf)%u,iostat=E_IO)trim(title)//end_rec
    write(unit=vtk(rf)%u,iostat=E_IO)trim(Upper_Case(output_format))//end_rec
    write(unit=vtk(rf)%u,iostat=E_IO)'DATASET '//trim(vtk(rf)%topology)//end_rec
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_INI

  !> @ingroup Lib_VTK_IOPrivateProcedure
  !> @{
  !> Function for saving mesh with \b STRUCTURED_POINTS topology (R8P).
  !> @return E_IO: integer(I4P) error flag
  function VTK_GEO_STRP_R8(cf,Nx,Ny,Nz,X0,Y0,Z0,Dx,Dy,Dz) result(E_IO)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN), optional:: cf       !< Current file index (for concurrent files IO).
  integer(I4P), intent(IN)::           Nx       !< Number of nodes in x direction.
  integer(I4P), intent(IN)::           Ny       !< Number of nodes in y direction.
  integer(I4P), intent(IN)::           Nz       !< Number of nodes in z direction.
  real(R8P),    intent(IN)::           X0       !< X coordinate of origin.
  real(R8P),    intent(IN)::           Y0       !< Y coordinate of origin.
  real(R8P),    intent(IN)::           Z0       !< Z coordinate of origin.
  real(R8P),    intent(IN)::           Dx       !< Space step in x direction.
  real(R8P),    intent(IN)::           Dy       !< Space step in y direction.
  real(R8P),    intent(IN)::           Dz       !< Space step in z direction.
  integer(I4P)::                       E_IO     !< Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done.
  character(len=maxlen)::              s_buffer !< Buffer string.
  integer(I4P)::                       rf       !< Real file index.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  rf = f
  if (present(cf)) then
    rf = cf ; f = cf
  endif
  select case(vtk(rf)%f)
  case(ascii)
    write(unit=vtk(rf)%u,fmt='(A,3'//FI4P//')',iostat=E_IO)'DIMENSIONS ',Nx,Ny,Nz
    write(unit=vtk(rf)%u,fmt='(A,3'//FR8P//')',iostat=E_IO)'ORIGIN ',X0,Y0,Z0
    write(unit=vtk(rf)%u,fmt='(A,3'//FR8P//')',iostat=E_IO)'SPACING ',Dx,Dy,Dz
  case(binary)
    write(s_buffer,      fmt='(A,3'//FI4P//')',iostat=E_IO)'DIMENSIONS ',Nx,Ny,Nz
    write(unit=vtk(rf)%u,                      iostat=E_IO)trim(s_buffer)//end_rec
    write(s_buffer,      fmt='(A,3'//FR8P//')',iostat=E_IO)'ORIGIN ',X0,Y0,Z0
    write(unit=vtk(rf)%u,                      iostat=E_IO)trim(s_buffer)//end_rec
    write(s_buffer,      fmt='(A,3'//FR8P//')',iostat=E_IO)'SPACING ',Dx,Dy,Dz
    write(unit=vtk(rf)%u,                      iostat=E_IO)trim(s_buffer)//end_rec
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_GEO_STRP_R8

  !> Function for saving mesh with \b STRUCTURED_POINTS topology (R4P).
  !> @return E_IO: integer(I4P) error flag
  function VTK_GEO_STRP_R4(cf,Nx,Ny,Nz,X0,Y0,Z0,Dx,Dy,Dz) result(E_IO)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN), optional:: cf       !< Current file index (for concurrent files IO).
  integer(I4P), intent(IN)::           Nx       !< Number of nodes in x direction.
  integer(I4P), intent(IN)::           Ny       !< Number of nodes in y direction.
  integer(I4P), intent(IN)::           Nz       !< Number of nodes in z direction.
  real(R4P),    intent(IN)::           X0       !< X coordinate of origin.
  real(R4P),    intent(IN)::           Y0       !< Y coordinate of origin.
  real(R4P),    intent(IN)::           Z0       !< Z coordinate of origin.
  real(R4P),    intent(IN)::           Dx       !< Space step in x direction.
  real(R4P),    intent(IN)::           Dy       !< Space step in y direction.
  real(R4P),    intent(IN)::           Dz       !< Space step in z direction.
  integer(I4P)::                       E_IO     !< Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done.
  character(len=maxlen)::              s_buffer !< Buffer string.
  integer(I4P)::                       rf       !< Real file index.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  rf = f
  if (present(cf)) then
    rf = cf ; f = cf
  endif
  select case(vtk(rf)%f)
  case(ascii)
    write(unit=vtk(rf)%u,fmt='(A,3'//FI4P//')',iostat=E_IO)'DIMENSIONS ',Nx,Ny,Nz
    write(unit=vtk(rf)%u,fmt='(A,3'//FR4P//')',iostat=E_IO)'ORIGIN ',X0,Y0,Z0
    write(unit=vtk(rf)%u,fmt='(A,3'//FR4P//')',iostat=E_IO)'SPACING ',Dx,Dy,Dz
  case(binary)
    write(s_buffer,      fmt='(A,3'//FI4P//')',iostat=E_IO)'DIMENSIONS ',Nx,Ny,Nz
    write(unit=vtk(rf)%u,                      iostat=E_IO)trim(s_buffer)//end_rec
    write(s_buffer,      fmt='(A,3'//FR4P//')',iostat=E_IO)'ORIGIN ',X0,Y0,Z0
    write(unit=vtk(rf)%u,                      iostat=E_IO)trim(s_buffer)//end_rec
    write(s_buffer,      fmt='(A,3'//FR4P//')',iostat=E_IO)'SPACING ',Dx,Dy,Dz
    write(unit=vtk(rf)%u,                      iostat=E_IO)trim(s_buffer)//end_rec
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_GEO_STRP_R4

  !> Function for saving mesh with \b STRUCTURED_GRID topology (R8P).
  !> @return E_IO: integer(I4P) error flag
  function VTK_GEO_STRG_R8(cf,Nx,Ny,Nz,NN,X,Y,Z) result(E_IO)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN), optional:: cf       !< Current file index (for concurrent files IO).
  integer(I4P), intent(IN)::           Nx       !< Number of nodes in x direction.
  integer(I4P), intent(IN)::           Ny       !< Number of nodes in y direction.
  integer(I4P), intent(IN)::           Nz       !< Number of nodes in z direction.
  integer(I4P), intent(IN)::           NN       !< Number of all nodes.
  real(R8P),    intent(IN)::           X(1:NN)  !< X coordinates.
  real(R8P),    intent(IN)::           Y(1:NN)  !< Y coordinates.
  real(R8P),    intent(IN)::           Z(1:NN)  !< Z coordinates.
  integer(I4P)::                       E_IO     !< Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done.
  character(len=maxlen)::              s_buffer !< Buffer string.
  integer(I4P)::                       rf       !< Real file index.
  integer(I4P)::                       n1       !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  rf = f
  if (present(cf)) then
    rf = cf ; f = cf
  endif
  select case(vtk(rf)%f)
  case(ascii)
    write(unit=vtk(rf)%u,fmt='(A,3'//FI4P//')', iostat=E_IO)'DIMENSIONS ',Nx,Ny,Nz
    write(unit=vtk(rf)%u,fmt='(A,'//FI4P//',A)',iostat=E_IO)'POINTS ',NN,' double'
    write(unit=vtk(rf)%u,fmt='(3'//FR8P//')',   iostat=E_IO)(X(n1),Y(n1),Z(n1),n1=1,NN)
  case(binary)
    write(s_buffer,      fmt='(A,3'//FI4P//')', iostat=E_IO)'DIMENSIONS ',Nx,Ny,Nz
    write(unit=vtk(rf)%u,                       iostat=E_IO)trim(s_buffer)//end_rec
    write(s_buffer,      fmt='(A,'//FI4P//',A)',iostat=E_IO)'POINTS ',NN,' double'
    write(unit=vtk(rf)%u,                       iostat=E_IO)trim(s_buffer)//end_rec
    write(unit=vtk(rf)%u,                       iostat=E_IO)(X(n1),Y(n1),Z(n1),n1=1,NN)
    write(unit=vtk(rf)%u,                       iostat=E_IO)end_rec
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_GEO_STRG_R8

  !> Function for saving mesh with \b STRUCTURED_GRID topology (R4P).
  !> @return E_IO: integer(I4P) error flag
  function VTK_GEO_STRG_R4(cf,Nx,Ny,Nz,NN,X,Y,Z) result(E_IO)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN), optional:: cf       !< Current file index (for concurrent files IO).
  integer(I4P), intent(IN)::           Nx       !< Number of nodes in x direction.
  integer(I4P), intent(IN)::           Ny       !< Number of nodes in y direction.
  integer(I4P), intent(IN)::           Nz       !< Number of nodes in z direction.
  integer(I4P), intent(IN)::           NN       !< Number of all nodes.
  real(R4P),    intent(IN)::           X(1:NN)  !< X coordinates.
  real(R4P),    intent(IN)::           Y(1:NN)  !< Y coordinates.
  real(R4P),    intent(IN)::           Z(1:NN)  !< Z coordinates.
  integer(I4P)::                       E_IO     !< Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done.
  character(len=maxlen)::              s_buffer !< Buffer string.
  integer(I4P)::                       rf       !< Real file index.
  integer(I4P)::                       n1       !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  rf = f
  if (present(cf)) then
    rf = cf ; f = cf
  endif
  select case(vtk(rf)%f)
  case(ascii)
    write(unit=vtk(rf)%u,fmt='(A,3'//FI4P//')', iostat=E_IO)'DIMENSIONS ',Nx,Ny,Nz
    write(unit=vtk(rf)%u,fmt='(A,'//FI4P//',A)',iostat=E_IO)'POINTS ',NN,' float'
    write(unit=vtk(rf)%u,fmt='(3'//FR4P//')',   iostat=E_IO)(X(n1),Y(n1),Z(n1),n1=1,NN)
  case(binary)
    write(s_buffer,      fmt='(A,3'//FI4P//')', iostat=E_IO)'DIMENSIONS ',Nx,Ny,Nz
    write(unit=vtk(rf)%u,                       iostat=E_IO)trim(s_buffer)//end_rec
    write(s_buffer,      fmt='(A,'//FI4P//',A)',iostat=E_IO)'POINTS ',NN,' float'
    write(unit=vtk(rf)%u,                       iostat=E_IO)trim(s_buffer)//end_rec
    write(unit=vtk(rf)%u,                       iostat=E_IO)(X(n1),Y(n1),Z(n1),n1=1,NN)
    write(unit=vtk(rf)%u,                       iostat=E_IO)end_rec
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_GEO_STRG_R4

  !> Function for saving mesh with \b RECTILINEAR_GRID topology (R8P).
  !> @return E_IO: integer(I4P) error flag
  function VTK_GEO_RECT_R8(cf,Nx,Ny,Nz,X,Y,Z) result(E_IO)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN), optional:: cf       !< Current file index (for concurrent files IO).
  integer(I4P), intent(IN)::           Nx       !< Number of nodes in x direction.
  integer(I4P), intent(IN)::           Ny       !< Number of nodes in y direction.
  integer(I4P), intent(IN)::           Nz       !< Number of nodes in z direction.
  real(R8P),    intent(IN)::           X(1:Nx)  !< X coordinates.
  real(R8P),    intent(IN)::           Y(1:Ny)  !< Y coordinates.
  real(R8P),    intent(IN)::           Z(1:Nz)  !< Z coordinates.
  integer(I4P)::                       E_IO     !< Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done.
  character(len=maxlen)::              s_buffer !< Buffer string.
  integer(I4P)::                       rf       !< Real file index.
  integer(I4P)::                       n1       !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  rf = f
  if (present(cf)) then
    rf = cf ; f = cf
  endif
  select case(vtk(rf)%f)
  case(ascii)
    write(unit=vtk(rf)%u,fmt='(A,3'//FI4P//')', iostat=E_IO)'DIMENSIONS ',Nx,Ny,Nz
    write(unit=vtk(rf)%u,fmt='(A,'//FI4P//',A)',iostat=E_IO)'X_COORDINATES ',Nx,' double'
    write(unit=vtk(rf)%u,fmt=FR8P,              iostat=E_IO)(X(n1),n1=1,Nx)
    write(unit=vtk(rf)%u,fmt='(A,'//FI4P//',A)',iostat=E_IO)'Y_COORDINATES ',Ny,' double'
    write(unit=vtk(rf)%u,fmt=FR8P,              iostat=E_IO)(Y(n1),n1=1,Ny)
    write(unit=vtk(rf)%u,fmt='(A,'//FI4P//',A)',iostat=E_IO)'Z_COORDINATES ',Nz,' double'
    write(unit=vtk(rf)%u,fmt=FR8P,              iostat=E_IO)(Z(n1),n1=1,Nz)
  case(binary)
    write(s_buffer,      fmt='(A,3'//FI4P//')', iostat=E_IO)'DIMENSIONS ',Nx,Ny,Nz
    write(unit=vtk(rf)%u,                       iostat=E_IO)trim(s_buffer)//end_rec
    write(s_buffer,      fmt='(A,'//FI4P//',A)',iostat=E_IO)'X_COORDINATES ',Nx,' double'
    write(unit=vtk(rf)%u,                       iostat=E_IO)trim(s_buffer)//end_rec
    write(unit=vtk(rf)%u,                       iostat=E_IO)(X(n1),n1=1,Nx)
    write(unit=vtk(rf)%u,                       iostat=E_IO)end_rec
    write(s_buffer,      fmt='(A,'//FI4P//',A)',iostat=E_IO)'Y_COORDINATES ',Ny,' double'
    write(unit=vtk(rf)%u,                       iostat=E_IO)trim(s_buffer)//end_rec
    write(unit=vtk(rf)%u,                       iostat=E_IO)(Y(n1),n1=1,Ny)
    write(unit=vtk(rf)%u,                       iostat=E_IO)end_rec
    write(s_buffer,      fmt='(A,'//FI4P//',A)',iostat=E_IO)'Z_COORDINATES ',Nz,' double'
    write(unit=vtk(rf)%u,                       iostat=E_IO)trim(s_buffer)//end_rec
    write(unit=vtk(rf)%u,                       iostat=E_IO)(Z(n1),n1=1,Nz)
    write(unit=vtk(rf)%u,                       iostat=E_IO)end_rec
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_GEO_RECT_R8

  !> Function for saving mesh with \b RECTILINEAR_GRID topology (R4P).
  !> @return E_IO: integer(I4P) error flag
  function VTK_GEO_RECT_R4(cf,Nx,Ny,Nz,X,Y,Z) result(E_IO)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN), optional:: cf       !< Current file index (for concurrent files IO).
  integer(I4P), intent(IN)::           Nx       !< Number of nodes in x direction.
  integer(I4P), intent(IN)::           Ny       !< Number of nodes in y direction.
  integer(I4P), intent(IN)::           Nz       !< Number of nodes in z direction.
  real(R4P),    intent(IN)::           X(1:Nx)  !< X coordinates.
  real(R4P),    intent(IN)::           Y(1:Ny)  !< Y coordinates.
  real(R4P),    intent(IN)::           Z(1:Nz)  !< Z coordinates.
  integer(I4P)::                       E_IO     !< Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done.
  character(len=maxlen)::              s_buffer !< Buffer string.
  integer(I4P)::                       rf       !< Real file index.
  integer(I4P)::                       n1       !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  rf = f
  if (present(cf)) then
    rf = cf ; f = cf
  endif
  select case(vtk(rf)%f)
  case(ascii)
    write(unit=vtk(rf)%u,fmt='(A,3'//FI4P//')', iostat=E_IO)'DIMENSIONS ',Nx,Ny,Nz
    write(unit=vtk(rf)%u,fmt='(A,'//FI4P//',A)',iostat=E_IO)'X_COORDINATES ',Nx,' float'
    write(unit=vtk(rf)%u,fmt=FR4P,              iostat=E_IO)(X(n1),n1=1,Nx)
    write(unit=vtk(rf)%u,fmt='(A,'//FI4P//',A)',iostat=E_IO)'Y_COORDINATES ',Ny,' float'
    write(unit=vtk(rf)%u,fmt=FR4P,              iostat=E_IO)(Y(n1),n1=1,Ny)
    write(unit=vtk(rf)%u,fmt='(A,'//FI4P//',A)',iostat=E_IO)'Z_COORDINATES ',Nz,' float'
    write(unit=vtk(rf)%u,fmt=FR4P,              iostat=E_IO)(Z(n1),n1=1,Nz)
  case(binary)
    write(s_buffer,      fmt='(A,3'//FI4P//')', iostat=E_IO)'DIMENSIONS ',Nx,Ny,Nz
    write(unit=vtk(rf)%u,                       iostat=E_IO)trim(s_buffer)//end_rec
    write(s_buffer,      fmt='(A,'//FI4P//',A)',iostat=E_IO)'X_COORDINATES ',Nx,' float'
    write(unit=vtk(rf)%u,                       iostat=E_IO)trim(s_buffer)//end_rec
    write(unit=vtk(rf)%u,                       iostat=E_IO)(X(n1),n1=1,Nx)
    write(unit=vtk(rf)%u,                       iostat=E_IO)end_rec
    write(s_buffer,      fmt='(A,'//FI4P//',A)',iostat=E_IO)'Y_COORDINATES ',Ny,' float'
    write(unit=vtk(rf)%u,                       iostat=E_IO)trim(s_buffer)//end_rec
    write(unit=vtk(rf)%u,                       iostat=E_IO)(Y(n1),n1=1,Ny)
    write(unit=vtk(rf)%u,                       iostat=E_IO)end_rec
    write(s_buffer,      fmt='(A,'//FI4P//',A)',iostat=E_IO)'Z_COORDINATES ',Nz,' float'
    write(unit=vtk(rf)%u,                       iostat=E_IO)trim(s_buffer)//end_rec
    write(unit=vtk(rf)%u,                       iostat=E_IO)(Z(n1),n1=1,Nz)
    write(unit=vtk(rf)%u,                       iostat=E_IO)end_rec
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_GEO_RECT_R4

  !> Function for saving mesh with \b UNSTRUCTURED_GRID topology (R8P).
  !> @return E_IO: integer(I4P) error flag
  function VTK_GEO_UNST_R8(cf,NN,X,Y,Z) result(E_IO)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN), optional:: cf       !< Current file index (for concurrent files IO).
  integer(I4P), intent(IN)::           NN       !< Number of nodes.
  real(R8P),    intent(IN)::           X(1:NN)  !< X coordinates of all nodes.
  real(R8P),    intent(IN)::           Y(1:NN)  !< Y coordinates of all nodes.
  real(R8P),    intent(IN)::           Z(1:NN)  !< Z coordinates of all nodes.
  integer(I4P)::                       E_IO     !< Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done.
  character(len=maxlen)::              s_buffer !< Buffer string.
  integer(I4P)::                       rf       !< Real file index.
  integer(I4P)::                       n1       !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  rf = f
  if (present(cf)) then
    rf = cf ; f = cf
  endif
  select case(vtk(rf)%f)
  case(ascii)
    write(unit=vtk(rf)%u,fmt='(A,'//FI4P//',A)',iostat=E_IO)'POINTS ',NN,' double'
    write(unit=vtk(rf)%u,fmt='(3'//FR8P//')',   iostat=E_IO)(X(n1),Y(n1),Z(n1),n1=1,NN)
  case(binary)
    write(s_buffer,      fmt='(A,'//FI4P//',A)',iostat=E_IO)'POINTS ',NN,' double'
    write(unit=vtk(rf)%u,                       iostat=E_IO)trim(s_buffer)//end_rec
    write(unit=vtk(rf)%u,                       iostat=E_IO)(X(n1),Y(n1),Z(n1),n1=1,NN)
    write(unit=vtk(rf)%u,                       iostat=E_IO)end_rec
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_GEO_UNST_R8

  !> Function for saving mesh with \b UNSTRUCTURED_GRID topology (R4P).
  !> @return E_IO: integer(I4P) error flag
  function VTK_GEO_UNST_R4(cf,NN,X,Y,Z) result(E_IO)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN), optional:: cf       !< Current file index (for concurrent files IO).
  integer(I4P), intent(IN)::           NN       !< number of nodes.
  real(R4P),    intent(IN)::           X(1:NN)  !< x coordinates of all nodes.
  real(R4P),    intent(IN)::           Y(1:NN)  !< y coordinates of all nodes.
  real(R4P),    intent(IN)::           Z(1:NN)  !< z coordinates of all nodes.
  integer(I4P)::                       E_IO     !< Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done.
  character(len=maxlen)::              s_buffer !< buffer string.
  integer(I4P)::                       rf       !< Real file index.
  integer(I4P)::                       n1       !< counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  rf = f
  if (present(cf)) then
    rf = cf ; f = cf
  endif
  select case(vtk(rf)%f)
  case(ascii)
    write(unit=vtk(rf)%u,fmt='(A,'//FI4P//',A)',iostat=E_IO)'POINTS ',NN,' float'
    write(unit=vtk(rf)%u,fmt='(3'//FR4P//')',   iostat=E_IO)(X(n1),Y(n1),Z(n1),n1=1,NN)
  case(binary)
    write(s_buffer,      fmt='(A,'//FI4P//',A)',iostat=E_IO)'POINTS ',NN,' float'
    write(unit=vtk(rf)%u,                       iostat=E_IO)trim(s_buffer)//end_rec
    write(unit=vtk(rf)%u,                       iostat=E_IO)(X(n1),Y(n1),Z(n1),n1=1,NN)
    write(unit=vtk(rf)%u,                       iostat=E_IO)end_rec
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_GEO_UNST_R4
  !> @}

  !> Function that \b must be used when unstructured grid is used, it saves the connectivity of the unstructured gird.
  !> @note The vector \b connect must follow the VTK-legacy standard. It is passed as \em assumed-shape array
  !> because its dimensions is related to the mesh dimensions in a complex way. Its dimensions can be calculated by the following
  !> equation: \f$dc = NC + \sum\limits_{i = 1}^{NC} {nvertex_i }\f$
  !> where \f$dc\f$ is connectivity vector dimension and \f$nvertex_i\f$ is the number of vertices of \f$i^{th}\f$ cell. The VTK-
  !> legacy standard for the mesh connectivity is quite obscure at least at first sight. It is more simple analyzing an example.
  !> Suppose we have a mesh composed by 2 cells, one hexahedron (8 vertices) and one pyramid with square basis (5 vertices) and
  !> suppose that the basis of pyramid is constitute by a face of the hexahedron and so the two cells share 4 vertices.
  !> The above equation !> gives \f$dc=2+8+5=15\f$. The connectivity vector for this mesh can be: \n
  !> first cell \n
  !> connect(1)  = 8  number of vertices of \f$1^\circ\f$ cell \n
  !> connect(2)  = 0  identification flag of \f$1^\circ\f$ vertex of 1° cell \n
  !> connect(3)  = 1  identification flag of \f$2^\circ\f$ vertex of 1° cell \n
  !> connect(4)  = 2  identification flag of \f$3^\circ\f$ vertex of 1° cell \n
  !> connect(5)  = 3  identification flag of \f$4^\circ\f$ vertex of 1° cell \n
  !> connect(6)  = 4  identification flag of \f$5^\circ\f$ vertex of 1° cell \n
  !> connect(7)  = 5  identification flag of \f$6^\circ\f$ vertex of 1° cell \n
  !> connect(8)  = 6  identification flag of \f$7^\circ\f$ vertex of 1° cell \n
  !> connect(9)  = 7  identification flag of \f$8^\circ\f$ vertex of 1° cell \n
  !> second cell \n
  !> connect(10) = 5  number of vertices of \f$2^\circ \f$cell \n
  !> connect(11) = 0  identification flag of \f$1^\circ\f$ vertex of 2° cell \n
  !> connect(12) = 1  identification flag of \f$2^\circ\f$ vertex of 2° cell \n
  !> connect(13) = 2  identification flag of \f$3^\circ\f$ vertex of 2° cell \n
  !> connect(14) = 3  identification flag of \f$4^\circ\f$ vertex of 2° cell \n
  !> connect(15) = 8  identification flag of \f$5^\circ\f$ vertex of 2° cell \n
  !> Note that the first 4 identification flags of pyramid vertices as the same of the first 4 identification flags of
  !> the hexahedron because the two cells share this face. It is also important to note that the identification flags start
  !> form $0$ value: this is impose to the VTK standard. The function VTK_CON does not calculate the connectivity vector: it
  !> writes the connectivity vector conforming the VTK standard, but does not calculate it.
  !> The vector variable \em cell_type must conform the VTK-legacy standard (see the file VTK-Standard at the
  !> Kitware homepage). It contains the
  !> \em type of each cells. For the above example this vector is: \n
  !> first cell \n
  !> cell_type(1) = 12 hexahedron type of \f$1^\circ\f$ cell \n
  !> second cell \n
  !> cell_type(2) = 14 pyramid type of \f$2^\circ\f$ cell \n
  !> @return E_IO: integer(I4P) error flag
  !> @ingroup Lib_VTK_IOPublicProcedure
  function VTK_CON(cf,NC,connect,cell_type) result(E_IO)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN), optional:: cf              !< Current file index (for concurrent files IO).
  integer(I4P), intent(IN)::           NC              !< Number of cells.
  integer(I4P), intent(IN)::           connect(:)      !< Mesh connectivity.
  integer(I4P), intent(IN)::           cell_type(1:NC) !< VTK cell type.
  integer(I4P)::                       E_IO            !< Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done.
  character(len=maxlen)::              s_buffer        !< Buffer string.
  integer(I4P)::                       ncon            !< Dimension of connectivity vector.
  integer(I4P)::                       rf              !< Real file index.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  rf = f
  if (present(cf)) then
    rf = cf ; f = cf
  endif
  ncon = size(connect,1)
  select case(vtk(rf)%f)
  case(ascii)
    write(unit=vtk(rf)%u,fmt='(A,2'//FI4P//')',iostat=E_IO)'CELLS ',NC,ncon
    write(unit=vtk(rf)%u,fmt=FI4P,             iostat=E_IO)connect
    write(unit=vtk(rf)%u,fmt='(A,'//FI4P//')', iostat=E_IO)'CELL_TYPES ',NC
    write(unit=vtk(rf)%u,fmt=FI4P,             iostat=E_IO)cell_type
  case(binary)
    write(s_buffer,      fmt='(A,2'//FI4P//')',iostat=E_IO)'CELLS ',NC,ncon
    write(unit=vtk(rf)%u,                      iostat=E_IO)trim(s_buffer)//end_rec
    write(unit=vtk(rf)%u,                      iostat=E_IO)connect
    write(unit=vtk(rf)%u,                      iostat=E_IO)end_rec
    write(s_buffer,      fmt='(A,'//FI4P//')', iostat=E_IO)'CELL_TYPES ',NC
    write(unit=vtk(rf)%u,                      iostat=E_IO)trim(s_buffer)//end_rec
    write(unit=vtk(rf)%u,                      iostat=E_IO)cell_type
    write(unit=vtk(rf)%u,                      iostat=E_IO)end_rec
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_CON

  !> Function that \b must be called before saving the data related to geometric mesh, this function initializes the
  !> saving of data variables indicating the \em type (node or cell centered) of variables that will be saved.
  !> @note A single file can contain both cell and node centered variables. In this case the VTK_DAT function must be
  !> called two times, before saving cell-centered variables and before saving node-centered variables.
  !> Examples of usage are: \n
  !> \b Node piece: \n
  !> @code ...
  !> E_IO=VTK_DAT(50,'node')
  !> ... @endcode
  !> \b Cell piece: \n
  !> @code ...
  !> E_IO=VTK_DAT(50,'cell')
  !> ... @endcode
  !> @return E_IO: integer(I4P) error flag
  !> @ingroup Lib_VTK_IOPublicProcedure
  function VTK_DAT(cf,NC_NN,var_location) result(E_IO)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN), optional:: cf           !< Current file index (for concurrent files IO).
  integer(I4P), intent(IN)::           NC_NN        !< Number of cells or nodes of field.
  character(*), intent(IN)::           var_location !< Location of saving variables: cell for cell-centered, node for node-centered.
  integer(I4P)::                       E_IO         !< Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done.
  character(len=maxlen)::              s_buffer     !< Buffer string.
  integer(I4P)::                       rf           !< Real file index.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  rf = f
  if (present(cf)) then
    rf = cf ; f = cf
  endif
  select case(vtk(rf)%f)
  case(ascii)
    select case(trim(Upper_Case(var_location)))
    case('CELL')
      write(unit=vtk(rf)%u,fmt='(A,'//FI4P//')',iostat=E_IO)'CELL_DATA ',NC_NN
    case('NODE')
      write(unit=vtk(rf)%u,fmt='(A,'//FI4P//')',iostat=E_IO)'POINT_DATA ',NC_NN
    endselect
  case(binary)
    select case(trim(Upper_Case(var_location)))
    case('CELL')
      write(s_buffer,fmt='(A,'//FI4P//')',iostat=E_IO)'CELL_DATA ',NC_NN
      write(unit=vtk(rf)%u,iostat=E_IO)trim(s_buffer)//end_rec
    case('NODE')
      write(s_buffer,fmt='(A,'//FI4P//')',iostat=E_IO)'POINT_DATA ',NC_NN
      write(unit=vtk(rf)%u,iostat=E_IO)trim(s_buffer)//end_rec
    endselect
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_DAT

  !> @ingroup Lib_VTK_IOPrivateProcedure
  !> @{
  !> Function for saving field of scalar variable (R8P).
  !> @return E_IO: integer(I4P) error flag
  function VTK_VAR_SCAL_R8(cf,NC_NN,varname,var) result(E_IO)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN), optional:: cf           !< Current file index (for concurrent files IO).
  integer(I4P), intent(IN)::           NC_NN        !< Number of nodes or cells.
  character(*), intent(IN)::           varname      !< Variable name.
  real(R8P),    intent(IN)::           var(1:NC_NN) !< Variable to be saved.
  integer(I4P)::                       E_IO         !< Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done.
  integer(I4P)::                       rf           !< Real file index.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  rf = f
  if (present(cf)) then
    rf = cf ; f = cf
  endif
  select case(vtk(rf)%f)
  case(ascii)
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)'SCALARS '//trim(varname)//' double 1'
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)'LOOKUP_TABLE default'
    write(unit=vtk(rf)%u,fmt=FR8P, iostat=E_IO)var
  case(binary)
    write(unit=vtk(rf)%u,iostat=E_IO)'SCALARS '//trim(varname)//' double 1'//end_rec
    write(unit=vtk(rf)%u,iostat=E_IO)'LOOKUP_TABLE default'//end_rec
    write(unit=vtk(rf)%u,iostat=E_IO)var
    write(unit=vtk(rf)%u,iostat=E_IO)end_rec
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_VAR_SCAL_R8

  !> Function for saving field of scalar variable (R4P).
  !> @return E_IO: integer(I4P) error flag
  function VTK_VAR_SCAL_R4(cf,NC_NN,varname,var) result(E_IO)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN), optional:: cf           !< Current file index (for concurrent files IO).
  integer(I4P), intent(IN)::           NC_NN        !< Number of nodes or cells.
  character(*), intent(IN)::           varname      !< Variable name.
  real(R4P),    intent(IN)::           var(1:NC_NN) !< Variable to be saved.
  integer(I4P)::                       E_IO         !< Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done.
  integer(I4P)::                       rf           !< Real file index.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  rf = f
  if (present(cf)) then
    rf = cf ; f = cf
  endif
  select case(vtk(rf)%f)
  case(ascii)
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)'SCALARS '//trim(varname)//' float 1'
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)'LOOKUP_TABLE default'
    write(unit=vtk(rf)%u,fmt=FR4P, iostat=E_IO)var
  case(binary)
    write(unit=vtk(rf)%u,iostat=E_IO)'SCALARS '//trim(varname)//' float 1'//end_rec
    write(unit=vtk(rf)%u,iostat=E_IO)'LOOKUP_TABLE default'//end_rec
    write(unit=vtk(rf)%u,iostat=E_IO)var
    write(unit=vtk(rf)%u,iostat=E_IO)end_rec
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_VAR_SCAL_R4

  !> Function for saving field of scalar variable (I4P).
  !> @return E_IO: integer(I4P) error flag
  function VTK_VAR_SCAL_I4(cf,NC_NN,varname,var) result(E_IO)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN), optional:: cf           !< Current file index (for concurrent files IO).
  integer(I4P), intent(IN)::           NC_NN        !< Number of nodes or cells.
  character(*), intent(IN)::           varname      !< Variable name.
  integer(I4P), intent(IN)::           var(1:NC_NN) !< Variable to be saved.
  integer(I4P)::                       E_IO         !< Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done.
  integer(I4P)::                       rf           !< Real file index.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  rf = f
  if (present(cf)) then
    rf = cf ; f = cf
  endif
  select case(vtk(rf)%f)
  case(ascii)
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)'SCALARS '//trim(varname)//' int 1'
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)'LOOKUP_TABLE default'
    write(unit=vtk(rf)%u,fmt=FI4P, iostat=E_IO)var
  case(binary)
    write(unit=vtk(rf)%u,iostat=E_IO)'SCALARS '//trim(varname)//' int 1'//end_rec
    write(unit=vtk(rf)%u,iostat=E_IO)'LOOKUP_TABLE default'//end_rec
    write(unit=vtk(rf)%u,iostat=E_IO)var
    write(unit=vtk(rf)%u,iostat=E_IO)end_rec
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_VAR_SCAL_I4

  !> Function for saving field of vectorial variable (R8P).
  !> @return E_IO: integer(I4P) error flag
  function VTK_VAR_VECT_R8(cf,vec_type,NC_NN,varname,varX,varY,varZ) result(E_IO)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN), optional:: cf            !< Current file index (for concurrent files IO).
  character(*), intent(IN)::           vec_type      !< Vector type: vect = generic vector , norm = normal vector.
  integer(I4P), intent(IN)::           NC_NN         !< Number of nodes or cells.
  character(*), intent(IN)::           varname       !< Variable name.
  real(R8P),    intent(IN)::           varX(1:NC_NN) !< X component of vector.
  real(R8P),    intent(IN)::           varY(1:NC_NN) !< Y component of vector.
  real(R8P),    intent(IN)::           varZ(1:NC_NN) !< Z component of vector.
  integer(I4P)::                       E_IO          !< Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done.
  integer(I4P)::                       rf            !< Real file index.
  integer(I8P)::                       n1            !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  rf = f
  if (present(cf)) then
    rf = cf ; f = cf
  endif
  select case(vtk(rf)%f)
  case(ascii)
    select case(Upper_Case(trim(vec_type)))
    case('VECT')
      write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)'VECTORS '//trim(varname)//' double'
    case('NORM')
      write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)'NORMALS '//trim(varname)//' double'
    endselect
    write(unit=vtk(rf)%u,fmt='(3'//FR8P//')',iostat=E_IO)(varX(n1),varY(n1),varZ(n1),n1=1,NC_NN)
  case(binary)
    select case(Upper_Case(trim(vec_type)))
    case('VECT')
      write(unit=vtk(rf)%u,iostat=E_IO)'VECTORS '//trim(varname)//' double'//end_rec
    case('NORM')
      write(unit=vtk(rf)%u,iostat=E_IO)'NORMALS '//trim(varname)//' double'//end_rec
    endselect
    write(unit=vtk(rf)%u,iostat=E_IO)(varX(n1),varY(n1),varZ(n1),n1=1,NC_NN)
    write(unit=vtk(rf)%u,iostat=E_IO)end_rec
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_VAR_VECT_R8

  !> Function for saving field of vectorial variable (R4P).
  !> @return E_IO: integer(I4P) error flag
  function VTK_VAR_VECT_R4(cf,vec_type,NC_NN,varname,varX,varY,varZ) result(E_IO)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN), optional:: cf            !< Current file index (for concurrent files IO).
  character(*), intent(IN)::           vec_type      !< Vector type: vect = generic vector , norm = normal vector.
  integer(I4P), intent(IN)::           NC_NN         !< Number of nodes or cells.
  character(*), intent(IN)::           varname       !< Variable name.
  real(R4P),    intent(IN)::           varX(1:NC_NN) !< X component of vector.
  real(R4P),    intent(IN)::           varY(1:NC_NN) !< Y component of vector.
  real(R4P),    intent(IN)::           varZ(1:NC_NN) !< Z component of vector.
  integer(I4P)::                       E_IO          !< Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done.
  integer(I4P)::                       rf            !< Real file index.
  integer(I8P)::                       n1            !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  rf = f
  if (present(cf)) then
    rf = cf ; f = cf
  endif
  select case(vtk(rf)%f)
  case(ascii)
    select case(Upper_Case(trim(vec_type)))
    case('vect')
      write(unit=vtk(rf)%u,fmt='(A)',          iostat=E_IO)'VECTORS '//trim(varname)//' float'
    case('norm')
      write(unit=vtk(rf)%u,fmt='(A)',          iostat=E_IO)'NORMALS '//trim(varname)//' float'
    endselect
    write(unit=vtk(rf)%u,fmt='(3'//FR4P//')',iostat=E_IO)(varX(n1),varY(n1),varZ(n1),n1=1,NC_NN)
  case(binary)
    select case(Upper_Case(trim(vec_type)))
    case('vect')
      write(unit=vtk(rf)%u,iostat=E_IO)'VECTORS '//trim(varname)//' float'//end_rec
    case('norm')
      write(unit=vtk(rf)%u,iostat=E_IO)'NORMALS '//trim(varname)//' float'//end_rec
    endselect
    write(unit=vtk(rf)%u,iostat=E_IO)(varX(n1),varY(n1),varZ(n1),n1=1,NC_NN)
    write(unit=vtk(rf)%u,iostat=E_IO)end_rec
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_VAR_VECT_R4

  !> Function for saving field of vectorial variable (I4P).
  !> @return E_IO: integer(I4P) error flag
  function VTK_VAR_VECT_I4(cf,NC_NN,varname,varX,varY,varZ) result(E_IO)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN), optional:: cf            !< Current file index (for concurrent files IO).
  integer(I4P), intent(IN)::           NC_NN         !< Number of nodes or cells.
  character(*), intent(IN)::           varname       !< Variable name.
  integer(I4P), intent(IN)::           varX(1:NC_NN) !< X component of vector.
  integer(I4P), intent(IN)::           varY(1:NC_NN) !< Y component of vector.
  integer(I4P), intent(IN)::           varZ(1:NC_NN) !< Z component of vector.
  integer(I4P)::                       E_IO          !< Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done.
  integer(I4P)::                       rf            !< Real file index.
  integer(I8P)::                       n1            !< Counter.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  rf = f
  if (present(cf)) then
    rf = cf ; f = cf
  endif
  select case(vtk(rf)%f)
  case(ascii)
    write(unit=vtk(rf)%u,fmt='(A)',iostat=E_IO)'VECTORS '//trim(varname)//' int'
    write(unit=vtk(rf)%u,fmt='(3'//FI4P//')',iostat=E_IO)(varX(n1),varY(n1),varZ(n1),n1=1,NC_NN)
  case(binary)
    write(unit=vtk(rf)%u,iostat=E_IO)'VECTORS '//trim(varname)//' int'//end_rec
    write(unit=vtk(rf)%u,iostat=E_IO)(varX(n1),varY(n1),varZ(n1),n1=1,NC_NN)
    write(unit=vtk(rf)%u,iostat=E_IO)end_rec
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_VAR_VECT_I4

  !> Function for saving texture variable (R8P).
  !> @return E_IO: integer(I4P) error flag
  function VTK_VAR_TEXT_R8(cf,NC_NN,dimm,varname,textCoo) result(E_IO)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN), optional:: cf                      !< Current file index (for concurrent files IO).
  integer(I4P), intent(IN)::           NC_NN                   !< Number of nodes or cells.
  integer(I4P), intent(IN)::           dimm                    !< Texture dimensions.
  character(*), intent(IN)::           varname                 !< Variable name.
  real(R8P),    intent(IN)::           textCoo(1:NC_NN,1:dimm) !< Texture.
  integer(I4P)::                       E_IO              !< Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done.
  character(len=maxlen)::              s_buffer                !< Buffer string.
  integer(I4P)::                       rf                      !< Real file index.
  integer(I8P)::                       n1,n2                   !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  rf = f
  if (present(cf)) then
    rf = cf ; f = cf
  endif
  select case(vtk(rf)%f)
  case(ascii)
    write(unit=vtk(rf)%u,fmt='(A,1X,'//FI4P//',1X,A)',iostat=E_IO)'TEXTURE_COORDINATES '//trim(varname),dimm,' double'
    write(s_buffer,fmt='(I1)',iostat=E_IO)dimm
    s_buffer='('//trim(s_buffer)//FR4P//')'
    write(unit=vtk(rf)%u,fmt=trim(s_buffer),iostat=E_IO)((textCoo(n1,n2),n2=1,dimm),n1=1,NC_NN)
  case(binary)
    write(s_buffer,fmt='(A,1X,'//FI4P//',1X,A)',iostat=E_IO)'TEXTURE_COORDINATES '//trim(varname),dimm,' double'
    write(unit=vtk(rf)%u,iostat=E_IO)trim(s_buffer)//end_rec
    write(unit=vtk(rf)%u,iostat=E_IO)((textCoo(n1,n2),n2=1,dimm),n1=1,NC_NN)
    write(unit=vtk(rf)%u,iostat=E_IO)end_rec
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_VAR_TEXT_R8

  !> Function for saving texture variable (R4P).
  !> @return E_IO: integer(I4P) error flag
  function VTK_VAR_TEXT_R4(cf,NC_NN,dimm,varname,textCoo) result(E_IO)
  !---------------------------------------------------------------------------------------------------------------------------------
  !! Function for saving texture variable (R4P).
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(IN), optional:: cf                      !< Current file index (for concurrent files IO).
  integer(I4P), intent(IN)::           NC_NN                   !< Number of nodes or cells.
  integer(I4P), intent(IN)::           dimm                    !< Texture dimensions.
  character(*), intent(IN)::           varname                 !< Variable name.
  real(R4P),    intent(IN)::           textCoo(1:NC_NN,1:dimm) !< Texture.
  integer(I4P)::                       E_IO              !< Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done.
  character(len=maxlen)::              s_buffer                !< Buffer string.
  integer(I4P)::                       rf                      !< Real file index.
  integer(I8P)::                       n1,n2                   !< Counters.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  rf = f
  if (present(cf)) then
    rf = cf ; f = cf
  endif
  select case(vtk(rf)%f)
  case(ascii)
    write(unit=vtk(rf)%u,fmt='(A,1X,'//FI4P//',1X,A)',iostat=E_IO)'TEXTURE_COORDINATES '//trim(varname),dimm,' float'
    write(s_buffer,fmt='(I1)',iostat=E_IO)dimm
    s_buffer='('//trim(s_buffer)//FR4P//')'
    write(unit=vtk(rf)%u,fmt=trim(s_buffer),iostat=E_IO)((textCoo(n1,n2),n2=1,dimm),n1=1,NC_NN)
  case(binary)
    write(s_buffer,fmt='(A,1X,'//FI4P//',1X,A)',iostat=E_IO)'TEXTURE_COORDINATES '//trim(varname),dimm,' float'
    write(unit=vtk(rf)%u,iostat=E_IO)trim(s_buffer)//end_rec
    write(unit=vtk(rf)%u,iostat=E_IO)((textCoo(n1,n2),n2=1,dimm),n1=1,NC_NN)
    write(unit=vtk(rf)%u,iostat=E_IO)end_rec
  endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_VAR_TEXT_R4
  !> @}

  !>Function for finalizing open file,  it has not inputs, @libvtk manages the file unit without the
  !>user's action.
  !> @note An example of usage is: \n
  !> @code ...
  !> E_IO=VTK_END()
  !> ... @endcode
  !> @return E_IO: integer(I4P) error flag
  !> @ingroup Lib_VTK_IOPublicProcedure
  function VTK_END(cf) result(E_IO)
  !---------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(I4P), intent(INOUT), optional:: cf   !< Current file index (for concurrent files IO).
  integer(I4P)::                          E_IO !< Input/Output inquiring flag: $0$ if IO is done, $> 0$ if IO is not done.
  integer(I4P)::                          rf   !< Real file index.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  rf = 1
  if (present(cf)) rf = cf
  close(unit=vtk(rf)%u,iostat=E_IO)
  call vtk_update(act='remove',cf=rf,Nvtk=Nvtk,vtk=vtk)
  f = rf
  if (present(cf)) cf = rf
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endfunction VTK_END
endmodule Lib_VTK_IO
