      REAL(KIND(0.0D0)) FUNCTION plgndr(l,m,x)
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN)      :: L, M
      REAL(DOUBLE), INTENT(IN) :: x
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!----------------------------------------------
      INTEGER      :: i, ll
      REAL(DOUBLE) :: fact,pll,pmm,pmmp1,somx2
!----------------------------------------------
!
!GG      if(m.lt.0.or.m.gt.l.or.abs(x).gt.1.d0)pause                      &
!GG                                              'bad arguments in plgndr'
      if(m.lt.0.or.m.gt.l.or.abs(x).gt.1.d0) then
         print*, "bad arguments in plgndr"
         stop
      end if
      pmm=1.d0
!                m
!     Compute   P
!                m
      if(m.gt.0) then
        somx2=sqrt((1.d0-x)*(1.d0+x))
        fact=1.d0
        do 11 i=1,m
          pmm=-pmm*fact*somx2
          fact=fact+2.d0
11      continue
      endif
      if(l.eq.m) then
        plgndr=pmm
      else
        pmmp1=x*(2*m+1)*pmm
!                m
!     Compute   P
!                m+1
        if(l.eq.m+1) then
          plgndr=pmmp1
        else
!                m
!     Compute   P  ,  l > m+1
!                l
          do 12 ll=m+2,l
            pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
            pmm=pmmp1
            pmmp1=pll
12        continue
          plgndr=pll
        endif
      endif
      return
      END FUNCTION plgndr
!  (C) Copr. 1986-92 Numerical Recipes Software oV)#S.
