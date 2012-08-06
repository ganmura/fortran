c
c     Numerical Analysis:
c     The Mathematics of Scientific Computing
c     Third Edition
c     D.R. Kincaid & E.W. Cheney
c     Brooks/Cole Publ., 2002
c     Copyright (c) 1996
c
c     Section 9.2
c
c     Example of solving a P.D.E using implicit method
c
c
c     file:  exs92.f
c
      parameter (n=10)
      dimension v(0:n+1),a(n),d(n),c(n),error(n)
      data M,hk/200,0.005125/  
c
      print *
      print *,' Implicit method example'
      print *,' Section 9.2, Kincaid-Cheney'
      print *
      print 7
c
      pi = 4.0*atan(1.0)
      h = 1.0/real(n+1)
      s = hk/h**2 
c
      do 2 i=1,n
         v(i) = sin(pi*real(i)*h) 
 2    continue
c
      t = 0.0
      print 8,0,t,v(0),v(1),v(10),v(11)
      print *
c
      do 3 i=1,n-1
         c(i) = -s
         a(i) = -s
 3    continue 
c
      do 6 j=1,M
         d(i) = 1.0 + 2.0*s
         call tri(n,a,d,c,v(1),v(1))
         t = real(j)*hk
         if (0 .eq. mod(j,10)) then 
            print 8,j,t,v(0),v(1),v(10),v(11)
            print *
         endif
c
         do 5 i=1,n
            x = real(i)*h
            error(i) = abs(sin(pi*x)/exp(pi*pi*t) - v(i))
 5       continue
c   
         if (0 .eq. mod(j,10)) then
            print *,' norm of error =', unorm(n,error)
            print *
         endif
 6    continue
c
 7    format (3x,'j',8x,'t',13x,'v(0)',11x,'v(1)',6x,'.....',
     +        5x,'v(10)',10x,'v(11)')
 8    format (1x,i3,2x,3(e13.6,2x),5x,2(e13.6,2x))
      stop
      end
c
      subroutine tri(n,a,d,c,b,x)     
c
c     solve a tridiagonal system
c
      dimension a(n),d(n),c(n),b(n),x(n) 
c
      do 2 i = 2,n
         xmult = a(i-1)/d(i-1) 
         d(i) = d(i) - xmult*c(i-1)    
         b(i) = b(i) - xmult*b(i-1)    
 2    continue    
c
      x(n) = b(n)/d(n)      
      do 3 i = n-1,1,-1     
         x(i) = (b(i) - c(i)*x(i+1))/d(i)
 3    continue    
c
      return      
      end 
c
      function unorm(n,u)
c
c     computes infinity norm of n component vector v
c
      dimension u(n)
      real max
c
      temp = 0.0
      do 2 i=1,n
         temp = max(temp,abs(u(i)))
 2    continue
      unorm = temp
c
      return
      end
