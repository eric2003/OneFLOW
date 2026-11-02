      program write_grid_unst
      use cgns
      implicit none
#include "cgnstypes_f03.h"
!
!   Creates simple 3-D unstructured grid and writes it to a
!   CGNS file.
!
!   Example compilation for this program is (change paths if needed!):
!   Note: when using the cgns module file, you must use the SAME fortran compiler
!   used to compile CGNS (see make.defs file)
!   ...or change, for example, via environment "setenv FC ifort"
!
!   ifort -I ../.. -c write_grid_unst.F90
!   ifort -o write_grid_unst write_grid_unst.o -L ../../lib -lcgns
!
!   (../../lib is the location where the compiled
!   library libcgns.a is located)
!
!   The following is no longer supported; now superseded by "use cgns":
!     include 'cgnslib_f.h'
!   Note Windows machines need to include cgnswin_f.h
!
      integer, parameter :: maxelemi=20*16*8
      integer, parameter :: maxelemj=1216
      real*8 x(21*17*9),y(21*17*9),z(21*17*9)
      integer i,j,k,ifirstnode,nbdyelem,ni,nj,nk
      integer index_file,index_section,ielem_no
      integer index_base,index_zone,index_coord
      integer ier,iset,iphysdim,icelldim
      integer(cgsize_t) isize(1,3),ielem(8,maxelemi),jelem(4,maxelemj)
      integer(cgsize_t) nelem_start,nelem_end
      character basename*32,zonename*32
!
      write(6,'('' Program write_grid_unst'')')
      if (CG_BUILD_64BIT) then
        write(6,'('' ...using 64-bit mode for particular integers'')')
      end if
!
!   create gridpoints for simple example:
      ni=21
      nj=17
      nk=9
      iset=0
      do k=1,nk
        do j=1,nj
          do i=1,ni
            iset=iset+1
            x(iset)=float(i-1)
            y(iset)=float(j-1)
            z(iset)=float(k-1)
          enddo
        enddo
      enddo
      write(6,'('' created simple 3-D grid points'')')
!
!   WRITE X, Y, Z GRID POINTS TO CGNS FILE
!   open CGNS file for write
      call cg_open_f('grid.cgns',CG_MODE_WRITE,index_file,ier)
      if (ier .ne. CG_OK) call cg_error_exit_f
!   create base (user can give any name)
      basename='Base'
      icelldim=3
      iphysdim=3
      call cg_base_write_f(index_file,basename,icelldim,iphysdim,index_base,ier)
!   define zone name (user can give any name)
      zonename = 'Zone  1'
!   vertex size
      isize(1,1)=ni*nj*nk
!   cell size
      isize(1,2)=(ni-1)*(nj-1)*(nk-1)
!   boundary vertex size (zero if elements not sorted)
      isize(1,3)=0
!   create zone
      call cg_zone_write_f(index_file,index_base,zonename,isize,               &
           CGNS_ENUMV(Unstructured),index_zone,ier)
!   write grid coordinates (user must use SIDS-standard names here)
      call cg_coord_write_f(index_file,index_base,index_zone,CGNS_ENUMV(RealDouble),       &
           'CoordinateX',x,index_coord,ier)
      call cg_coord_write_f(index_file,index_base,index_zone,CGNS_ENUMV(RealDouble),       &
           'CoordinateY',y,index_coord,ier)
      call cg_coord_write_f(index_file,index_base,index_zone,CGNS_ENUMV(RealDouble),       &
           'CoordinateZ',z,index_coord,ier)
!  set element connectivity:
!  ----------------------------------------------------------
!  do all the HEXA_8 elements (this part is mandatory):
!  maintain SIDS-standard ordering
      ielem_no=0
!  index no of first element
      nelem_start=1
      do k=1,nk-1
        do j=1,nj-1
          do i=1,ni-1
            ielem_no=ielem_no+1
!  in this example, due to the order in the node numbering, the
!  hexahedral elements can be reconstructed using the following
!  relationships:
            ifirstnode=i+(j-1)*ni+(k-1)*ni*nj
            ielem(1,ielem_no)=ifirstnode
            ielem(2,ielem_no)=ifirstnode+1
            ielem(3,ielem_no)=ifirstnode+1+ni
            ielem(4,ielem_no)=ifirstnode+ni
            ielem(5,ielem_no)=ifirstnode+ni*nj
            ielem(6,ielem_no)=ifirstnode+ni*nj+1
            ielem(7,ielem_no)=ifirstnode+ni*nj+1+ni
            ielem(8,ielem_no)=ifirstnode+ni*nj+ni
          enddo
        enddo
      enddo
!  index no of last element (=2560)
      nelem_end=ielem_no
      if (nelem_end .gt. maxelemi) then
        write(6,'('' Error, must increase maxelemi to at least '',             &
         i7)') nelem_end
        stop
      end if
!  unsorted boundary elements
      nbdyelem=0
!  write HEXA_8 element connectivity (user can give any name)
      call cg_section_write_f(index_file,index_base,index_zone,                &
           'Elem',CGNS_ENUMV(HEXA_8),nelem_start,nelem_end,nbdyelem,ielem,                 &
           index_section,ier)
!  ----------------------------------------------------------
!  do boundary (QUAD) elements (this part is optional,
!  but you must do it if you eventually want to define BCs
!  at element faces rather than at nodes):
!  maintain SIDS-standard ordering
!  INFLOW:
      ielem_no=0
!  index no of first element
      nelem_start=nelem_end+1
      i=1
      do k=1,nk-1
        do j=1,nj-1
          ielem_no=ielem_no+1
          ifirstnode=i+(j-1)*ni+(k-1)*ni*nj
          jelem(1,ielem_no)=ifirstnode
          jelem(2,ielem_no)=ifirstnode+ni*nj
          jelem(3,ielem_no)=ifirstnode+ni*nj+ni
          jelem(4,ielem_no)=ifirstnode+ni
        enddo
      enddo
!  index no of last element
      nelem_end=nelem_start+ielem_no-1
      if (ielem_no .gt. maxelemj) then
        write(6,'('' Error, must increase maxelemj to at least '',             &
         i7)') ielem_no
        stop
      end if
!  write QUAD element connectivity for inflow face (user can give any name)
      call cg_section_write_f(index_file,index_base,index_zone,                &
           'InflowElem',CGNS_ENUMV(QUAD_4),nelem_start,nelem_end,nbdyelem,                 &
           jelem,index_section,ier)
!  OUTFLOW:
      ielem_no=0
!  index no of first element
      nelem_start=nelem_end+1
      i=ni-1
      do k=1,nk-1
        do j=1,nj-1
          ielem_no=ielem_no+1
          ifirstnode=i+(j-1)*ni+(k-1)*ni*nj
          jelem(1,ielem_no)=ifirstnode+1
          jelem(2,ielem_no)=ifirstnode+1+ni
          jelem(3,ielem_no)=ifirstnode+ni*nj+1+ni
          jelem(4,ielem_no)=ifirstnode+ni*nj+1
        enddo
      enddo
!  index no of last element
      nelem_end=nelem_start+ielem_no-1
      if (ielem_no .gt. maxelemj) then
        write(6,'('' Error, must increase maxelemj to at least '',             &
         i7)') ielem_no
        stop
      end if
!  write QUAD element connectivity for outflow face (user can give any name)
      call cg_section_write_f(index_file,index_base,index_zone,                &
           'OutflowElem',CGNS_ENUMV(QUAD_4),nelem_start,nelem_end,nbdyelem,                &
           jelem,index_section,ier)
!  SIDEWALLS:
      ielem_no=0
!  index no of first element
      nelem_start=nelem_end+1
      j=1
      do k=1,nk-1
        do i=1,ni-1
          ielem_no=ielem_no+1
          ifirstnode=i+(j-1)*ni+(k-1)*ni*nj
          jelem(1,ielem_no)=ifirstnode
          jelem(2,ielem_no)=ifirstnode+ni*nj
          jelem(3,ielem_no)=ifirstnode+ni*nj+1
          jelem(4,ielem_no)=ifirstnode+1
        enddo
      enddo
      j=nj-1
      do k=1,nk-1
        do i=1,ni-1
          ielem_no=ielem_no+1
          ifirstnode=i+(j-1)*ni+(k-1)*ni*nj
          jelem(1,ielem_no)=ifirstnode+1+ni
          jelem(2,ielem_no)=ifirstnode+ni
          jelem(3,ielem_no)=ifirstnode+ni*nj+ni
          jelem(4,ielem_no)=ifirstnode+ni*nj+1+ni
        enddo
      enddo
      k=1
      do j=1,nj-1
        do i=1,ni-1
          ielem_no=ielem_no+1
          ifirstnode=i+(j-1)*ni+(k-1)*ni*nj
          jelem(1,ielem_no)=ifirstnode
          jelem(2,ielem_no)=ifirstnode+1
          jelem(3,ielem_no)=ifirstnode+1+ni
          jelem(4,ielem_no)=ifirstnode+ni
        enddo
      enddo
      k=nk-1
      do j=1,nj-1
        do i=1,ni-1
          ielem_no=ielem_no+1
          ifirstnode=i+(j-1)*ni+(k-1)*ni*nj
          jelem(1,ielem_no)=ifirstnode+ni*nj
          jelem(2,ielem_no)=ifirstnode+ni*nj+ni
          jelem(3,ielem_no)=ifirstnode+ni*nj+1+ni
          jelem(4,ielem_no)=ifirstnode+ni*nj+1
        enddo
      enddo
!  index no of last element
      nelem_end=nelem_start+ielem_no-1
      if (ielem_no .gt. maxelemj) then
        write(6,'('' Error, must increase maxelemj to at least '',             &
         i7)') ielem_no
        stop
      end if
!  write QUAD element connectivity for sidewall face (user can give any name)
      call cg_section_write_f(index_file,index_base,index_zone,                &
           'SidewallElem',CGNS_ENUMV(QUAD_4),nelem_start,nelem_end,nbdyelem,               &
           jelem,index_section,ier)
!  ----------------------------------------------------------
!   close CGNS file
      call cg_close_f(index_file,ier)
      write(6,'('' Successfully wrote unstructured grid to file'',             &
        '' grid.cgns'')')
      stop
      end
