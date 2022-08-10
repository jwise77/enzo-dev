#define NODEBUG_MRP
#include "fortran.def"
c=======================================================================
c//////////////////////  SUBROUTINE CIC_FLAG  \\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine cic_flag(irefflag, posx, posy, posz, ndim, npositions, 
     &                    ffield, leftedge, 
     &                    dim1, dim2, dim3, cellsize, buffersize)
#ifndef CONFIG_PFLOAT_16
c
c  PERFORMS 1/2/3D CLOUD-IN-CELL MARKING OF FLAGGING FIELD
c
c  written by: Greg Bryan
c  date:       September, 2003
c  modified1:
c
c  PURPOSE: This routine performs a 1,2 or 3 dimension setting of a
c           flagging field for must-refine particles.

c
c  INPUTS:
c     irefflag   - flag per particle marking it as must refine
c     ndim       - dimensionality
c     cellsize   - the cell size of field
c     dim1,2,3   - real dimensions of field
c     leftedge   - the left edge(s) of field
c     npositions - number of particles
c     posx,y,z   - particle positions
c     buffersize - how many cells to flag around the particle
c
c  OUTPUT ARGUMENTS: 
c     ffield      - field to be deposited to
c
c  EXTERNALS: 
c
c  LOCALS:
c
c-----------------------------------------------------------------------
c
      implicit NONE
#include "fortran_types.def"
c
c-----------------------------------------------------------------------
c
c  argument declarations
c
      INTG_PREC dim1, dim2, dim3, npositions, ndim
      INTG_PREC buffersize
      INTG_PREC irefflag(npositions)
      P_PREC posx(npositions), posy(npositions), posz(npositions), 
     &     leftedge(3)
      P_PREC cellsize
      INTG_PREC ffield(dim1, dim2, dim3)
c
c  locals
c
      INTG_PREC i1, j1, k1, n, ii, jj, kk
      INTG_PREC istart, iend, jstart, jend, kstart, kend
      R_PREC    xpos, ypos, zpos
      P_PREC edge1, edge2, edge3, half, fact
      parameter (half = 0.5001_PKIND)
      INTG_PREC nflag, dim
      R_PREC rmin(3), rmax(3), pmin(3), pmax(3)

      do dim = 1, 3
            rmin(dim) = 1
            rmax(dim) = 0
            pmin(dim) = 1
            pmax(dim) = 0
      enddo

c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////
c=======================================================================
c
      nflag = 0
      fact = 1._PKIND/cellsize
      edge1 = REAL(dim1,PKIND) - half
      edge2 = REAL(dim2,PKIND) - half
      edge3 = REAL(dim3,PKIND) - half
c
c     1D
c
      if (ndim .eq. 1) then
c
         do n=1, npositions
c
c           only do this for must-refine particles
c
            if (irefflag(n) .eq. 1) then
c
c             Compute the position of the central cell
c
              xpos = min(max((posx(n) - leftedge(1))*fact, half), edge1)
c
c             Convert this into an INTG_PREC index
c
              i1  = int(xpos + 0.5_RKIND,IKIND)
c
c             set flagging field
c
c             keep default option unrolled!
              if (buffersize .eq. 1) then
                 ffield(i1  ,1,1) = 1
                 ffield(i1+1,1,1) = 1
              else
                 istart = max(i1-buffersize+1, 1)
                 iend = min(i1+buffersize, dim1)
                 do ii = istart, iend
                    ffield(ii,1,1) = 1
                 enddo
              endif
c
	    endif
c
         enddo
c
      endif
c
c     2D
c
      if (ndim .eq. 2) then
c
         do n=1, npositions
c
c           only do this for must-refine particles
c
            if (irefflag(n) .eq. 1) then
c
c             Compute the position of the central cell
c
              xpos = min(max((posx(n) - leftedge(1))*fact, half), edge1)
              ypos = min(max((posy(n) - leftedge(2))*fact, half), edge2)
c
c             Convert this into an INTG_PREC index
c
              i1  = int(xpos + 0.5_RKIND,IKIND)
              j1  = int(ypos + 0.5_RKIND,IKIND)
c
c             Interpolate from field into sumfield
c
c             keep default option unrolled!
              if (buffersize .eq. 1) then
                 ffield(i1  ,j1  ,1) = 1
                 ffield(i1+1,j1  ,1) = 1
                 ffield(i1  ,j1+1,1) = 1
                 ffield(i1+1,j1+1,1) = 1
              else
                 istart = max(i1-buffersize+1, 1)
                 iend = min(i1+buffersize, dim1)
                 jstart = max(j1-buffersize+1, 1)
                 jend = min(j1+buffersize, dim2)
                 do jj = jstart, jend
                    do ii = istart, iend
                       ffield(ii,jj,1) = 1
                    enddo
                 enddo
              endif
c
            endif
c
         enddo
c
      endif
c
c     3D
c
      if (ndim .eq. 3) then
c
         do n=1, npositions
c
c           only do this for must-refine particles
c
            if (irefflag(n) .eq. 1) then
c
c             Compute the position of the central cell
c
              xpos = min(max((posx(n) - leftedge(1))*fact, half), edge1)
              ypos = min(max((posy(n) - leftedge(2))*fact, half), edge2)
              zpos = min(max((posz(n) - leftedge(3))*fact, half), edge3)

              pmin(1) = min(pmin(1), posx(n))
              pmin(2) = min(pmin(2), posy(n))
              pmin(3) = min(pmin(3), posz(n))
              pmax(1) = max(pmax(1), posx(n))
              pmax(2) = max(pmax(2), posy(n))
              pmax(3) = max(pmax(3), posz(n))

c  
c             Convert this into an INTG_PREC index
c  
              i1  = int(xpos + 0.5_RKIND,IKIND)
              j1  = int(ypos + 0.5_RKIND,IKIND)
              k1  = int(zpos + 0.5_RKIND,IKIND)

              rmin(1) = min(rmin(1), (i1-0.5)*cellsize+leftedge(1))
              rmin(2) = min(rmin(2), (j1-0.5)*cellsize+leftedge(2))
              rmin(3) = min(rmin(3), (k1-0.5)*cellsize+leftedge(3))
              rmax(1) = max(rmax(1), (i1+0.5)*cellsize+leftedge(1))
              rmax(2) = max(rmax(2), (j1+0.5)*cellsize+leftedge(2))
              rmax(3) = max(rmax(3), (k1+0.5)*cellsize+leftedge(3))

c
c             Set flagging field
c     
c             keep default option unrolled!
              if (buffersize .eq. 1) then
                 ffield(i1  ,j1  ,k1  ) = 1
                 ffield(i1+1,j1  ,k1  ) = 1
                 ffield(i1  ,j1+1,k1  ) = 1
                 ffield(i1+1,j1+1,k1  ) = 1
                 ffield(i1  ,j1  ,k1+1) = 1
                 ffield(i1+1,j1  ,k1+1) = 1
                 ffield(i1  ,j1+1,k1+1) = 1
                 ffield(i1+1,j1+1,k1+1) = 1
                 nflag = nflag + 1
            else
                 istart = max(i1-buffersize+1, 1)
                 iend = min(i1+buffersize, dim1)
                 jstart = max(j1-buffersize+1, 1)
                 jend = min(j1+buffersize, dim2)
                 kstart = max(k1-buffersize+1, 1)
                 kend = min(k1+buffersize, dim3)
c	         print*,'cic_flag: buffersize is:',buffersize
                 do kk = kstart, kend
                    do jj = jstart, jend
                       do ii = istart, iend
                          ffield(ii,jj,kk) = 1
                       enddo
                    enddo
                 enddo
              endif
c
	  endif
c
         enddo
c
      endif
c
#ifdef DEBUG_MRP
      if (npositions .eq. 17848184) then
            write(*,*) 'cic_flag:', rmin(1), rmin(2), rmin(3), '//', 
     &            rmax(1), 
     &            rmax(2), rmax(3)
            write(*,*) 'cic_flag:', pmin(1), pmin(2), pmin(3), '//', 
     &            pmax(1), 
     &            pmax(2), pmax(3)
      endif
      write(*,*) 'cic_flag: nflag = ', nflag, ', npart = ', npositions
#endif
      return
#endif
      end
