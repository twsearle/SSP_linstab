!A subroutine for transforming a complex array of numbers to real space
      SUBROUTINE RStransform(inarr, oarr, N, M, numXs, numYs, kx)
          implicit none
          integer :: N 
          integer :: M 
          integer :: numXs 
          integer :: numYs 
          integer :: ndx
          integer :: mdx
          integer :: xIndx
          integer :: yIndx
          real :: kx
          real :: xpos
          real :: ypos
          real :: pi
          real :: tmp3
          double complex :: tmp
          double complex :: tmp2
          double complex :: II = (0.0,1.0)
          double complex ,dimension((2*N+1)*M) :: inarr
          double complex ,dimension(numXs, numYs) :: oarr
          !f2py intent(in) inarr
          !f2py intent(out) oarr
          !f2py depend(numXs) oarr
          !f2py depend(numYs) oarr
          !f2py depend(N) inarr
          !f2py depend(M) inarr 



          pi = 3.141592653589793 

          do xIndx = 1,numXs
            do yIndx = 1,numYs
            xpos = (2*pi/kx) * ((real(xIndx-1)/real(numXs-1.0)) - 0.5)
            ypos = 2.0*real(yIndx-1)/(real(numYs)-1.0)-1.0
              do ndx = 0,(2*N)
                tmp2 = exp(II*cmplx(ndx-N)*cmplx(kx*xpos))
                do mdx = 1,M
                  tmp3 = cos((mdx-1)*acos(ypos))
                  tmp = inarr(ndx*M+mdx)*tmp2*cmplx(tmp3)
                  oarr(yIndx, xIndx) = oarr(yIndx, xIndx)+tmp
                end do
              end do
            end do
          end do

          END
