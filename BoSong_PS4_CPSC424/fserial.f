      parameter (MAXB=100000)
      integer MAXB


      double precision p(MAXB,3), v(MAXB,3), m(MAXB)
      double precision deltaf(3), f(MAXB,3)
      double precision dt, dt2, g, r, r2, vavg(3), pcom(3) 
      double precision ax, ay, az, mj, mx, my, mz, mtotal
      double precision time0, time1, cptime, wctime
      integer nb, nts, nt, i, j, thisbody, otherbody

      g = 1.0

      read (5,*) nb
      read (5,*) nts
      read (5,*) dt
      if (nb .eq. 0) then
         print *, 'Nothing to do'
         stop
      endif
      print *, 'N-Body Problem (Serial Run)  N = ',nb
      do i=1,nb
         read (5,*) m(i)
      enddo
c      print 992,(m(i),i=1,nb)
 992  format(' masses:'/(1x, 10g13.4))
      do i=1,nb
         read (5,*) p(i,1), p(i,2), p(i,3)
      enddo
c      print 993,(p(i,1),p(i,2),p(i,3),i=1,nb)
 993  format(' positions:'/(1x,3g20.8))
      do i=1,nb
         read (5,*) v(i,1), v(i,2), v(i,3)
      enddo
c      print 994,(v(i,1),v(i,2),v(i,3),i=1,nb)
 994  format(' velocities:'/(1x,3g20.8))

      mtotal = 0.
      do i=1,nb
         mtotal = mtotal + m(i)
      enddo
c      print *, 'mtotal = ', mtotal

      dt2 = dt / 2.
      
      call timing(time0,cptime)

c Timestep loop

      do nt=0,nts

c Output if needed

         if (mod(nt,128).eq.0 .or. nt.eq.nts) then
            mx = 0.
            my = 0.
            mz = 0.
            px = 0.
            vavg(1) = 0.
            vavg(2) = 0.
            vavg(3) = 0.
            do i=1,nb
               mx = mx + m(i)*p(i,1)
               my = my + m(i)*p(i,2)
               mz = mz + m(i)*p(i,3)
               vavg(1) = vavg(1) + v(i,1)
               vavg(2) = vavg(2) + v(i,2)
               vavg(3) = vavg(3) + v(i,3)
            enddo
            vavg(1) = vavg(1)/nb
            vavg(2) = vavg(2)/nb
            vavg(3) = vavg(3)/nb
c            print *,'directional mass sums:',mx, my, mz
            pcom(1) = mx/mtotal
            pcom(2) = my/mtotal
            pcom(3) = mz/mtotal
            
            time = nt * dt
            if (nt.eq.0) then
               print 990, time, pcom(1), pcom(2), pcom(3),
     1                    vavg(1), vavg(2), vavg(3)
            else 
               print 991, nt, time, pcom(1), pcom(2), pcom(3),
     1                    vavg(1), vavg(2), vavg(3)
            endif

 990        format(/' Initial Conditions (time = ', g16.8, '):'//
     1             '     Center of Mass:   (',
     2                      g16.8, ', ', g16.8, ', ', g16.8, ' )'/
     3             '     Average Velocity: (',
     4                      g16.8, ', ', g16.8, ', ', g16.8, ' )' )

 991        format(/' Conditions after timestep ', i5, 
     1             ' (time = ', g16.8, '):'//
     2             '     Center of Mass:   (',
     3                      g16.8, ', ', g16.8, ', ', g16.8, ' )'/
     4             '     Average Velocity: (',
     5                      g16.8, ', ', g16.8, ', ', g16.8, ' )' )
         endif

c Computation

         if (nt .lt. nts) then

c  Initialize forces
            do thisbody=1,nb
               f(thisbody,1) = 0.
               f(thisbody,2) = 0.
               f(thisbody,3) = 0.
            enddo

c  Compute all pairwise interbody forces
            do thisbody=1,nb
               do otherbody=thisbody+1,nb
                  call force(MAXB, thisbody, otherbody, deltaf, p, m, g)
                  f(thisbody,1) = f(thisbody,1) + deltaf(1)
                  f(thisbody,2) = f(thisbody,2) + deltaf(2)
                  f(thisbody,3) = f(thisbody,3) + deltaf(3)
                  f(otherbody,1) = f(otherbody,1) - deltaf(1)
                  f(otherbody,2) = f(otherbody,2) - deltaf(2)
                  f(otherbody,3) = f(otherbody,3) - deltaf(3)
               enddo
            enddo

c Move the bodies
            do thisbody=1,nb
               ax = f(thisbody,1) / m(thisbody)
               ay = f(thisbody,2) / m(thisbody)
               az = f(thisbody,3) / m(thisbody)
               vavgx = v(thisbody,1) + dt2*ax
               vavgy = v(thisbody,2) + dt2*ay
               vavgz = v(thisbody,3) + dt2*az
               p(thisbody,1) = p(thisbody,1) + dt*vavgx
               p(thisbody,2) = p(thisbody,2) + dt*vavgy
               p(thisbody,3) = p(thisbody,3) + dt*vavgz
               v(thisbody,1) = v(thisbody,1) + dt*ax
               v(thisbody,2) = v(thisbody,2) + dt*ay
               v(thisbody,3) = v(thisbody,3) + dt*az
            enddo
         endif
      enddo

      call timing(time1, cptime)
      wctime = time1 - time0
      print 999, nts, nb, wctime
 999  format(/' Elapsed time for ', i5, ' timesteps with ',
     1         i5, ' bodies: ', f9.4, ' seconds')

      stop
      end

      subroutine force(MAXB, body1, body2, deltaf, p, m, g)

      integer body1, body2
      real*8 deltaf(3), p(MAXB,3), m(MAXB)
      real*8 gmmr3, r, r2, dx, dy, dz
      real*8 g

      dx = p(body2,1) - p(body1,1)
      dy = p(body2,2) - p(body1,2)
      dz = p(body2,3) - p(body1,3)
      r2 = dx*dx + dy*dy + dz*dz
      r = sqrt(r2)

      if (r .le. 5.0) then
         gmmr3 = g * m(body1) * m(body2) / (r2 * r)
         deltaf(1) = gmmr3 * dx
         deltaf(2) = gmmr3 * dy
         deltaf(3) = gmmr3 * dz
      else
         deltaf(1) = 0.
         deltaf(2) = 0.
         deltaf(3) = 0.
      endif

      return
      end

      




