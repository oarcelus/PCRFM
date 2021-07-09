module pseudocoulomb

   ! This module contains functions to calculates the repulsive pseudo 
   ! electrostatic forces in each voxel belonging to the selected phase. 
   ! And it calculates trajectories of possitively charged particles inside
   ! the phase in order to identify particles.

   implicit none

   integer, save,   public                 :: npart
   integer, save,   public, allocatable    :: nodes(:)
   integer, save,   public, allocatable    :: Pmat(:,:), part_pospix(:,:)
   real,    save,   public, allocatable    :: pseudo_f(:,:)

   !Subroutines
   public  :: calc_f_field
   public  :: calc_particles

   private :: line_sight
   private :: is_phase

   contains

   function is_phase(x,y,X1,Y1,X2,Y2) result (res)

      use images

      integer, intent(in)    :: x, y, X1, Y1, X2, Y2
      logical                :: res, neigh_pix

      if ((x==(X1-1).and.y==Y1)     .or. &
          (x==(X1+1).and.y==Y1)     .or. &
          (x==X1.and.y==(Y1-1))     .or. &
          (x==X1.and.y==(Y1+1))     .or. &
          (x==(X1-1).and.y==(Y1-1)) .or. &
          (x==(X1-1).and.y==(Y1+1)) .or. &
          (x==(X1+1).and.y==(Y1-1)) .or. &
          (x==(X1+1).and.y==(Y1+1))) then
           
         neigh_pix = .true.

      end if


      if ((img(x,y) == phase).or.(x==X2).or.(y==Y2).or.neigh_pix) then
         res = .true.
      else
         res = .false.
      end if

   end function is_phase

   function line_sight(X1,Y1,X2,Y2) result (l)

      use images

      integer             :: X, Y, dX, dY, e, l, i
      integer, intent(in) :: X1, Y1, X2, Y2

      ! 1st quadrant
      if ((X2 >= X1) .and. (Y2 >= Y1)) then
         X = X1
         Y = Y1

         dX = X2 - X1
         dY = Y2 - Y1

         if (dX >= dY) then
            e = 3*dY - dX
            do i = 1, dX
               if (is_phase(X,Y,X1,Y1,X2,Y2)) then
                  l = 1
               else
                  l = 0
                  exit
               end if

               X = X + 1

               if (e > 0) then
                  Y = Y + 1

                  if (is_phase(X,Y,X1,Y1,X2,Y2)) then
                     l = 1
                  else
                     l = 0
                     exit
                  end if

                  e = e + 2*(dY - dX)

               else
                  e = e + 2*dY
               end if
            end do

         else if (dY > dX) then
            e = 3*dX - dY
            do i = 1, dY
               if (is_phase(X,Y,X1,Y1,X2,Y2)) then
                  l = 1
               else
                  l = 0
                  exit
               end if

               Y = Y + 1

               if (e > 0) then
                  X = X + 1

                  if (is_phase(X,Y,X1,Y1,X2,Y2)) then
                     l = 1
                  else
                     l = 0
                     exit
                  end if

                  e = e + 2*(dX - dY)

               else
                  e = e + 2*dX
               end if
            end do
         end if

      ! 2nd quadrant
      else if ((X2 < X1) .and. (Y2 >= Y1)) then
         X = X1
         Y = Y1

         dX = X1 - X2
         dY = Y2 - Y1

         if (dX >= dY) then
            e = 3*dY - dX
            do i = 1, dX
               if (is_phase(X,Y,X1,Y1,X2,Y2)) then
                  l = 1
               else
                  l = 0
                  exit
               end if

               X = X - 1

               if (e > 0) then
                  Y = Y + 1

                  if (is_phase(X,Y,X1,Y1,X2,Y2)) then
                     l = 1
                  else
                     l = 0
                     exit
                  end if

                  e = e + 2*(dY - dX)

               else
                  e = e + 2*dY
               end if
            end do

         else if (dY > dX) then
            e = 3*dX - dY
            do i = 1, dY
               if (is_phase(X,Y,X1,Y1,X2,Y2)) then
                  l = 1
               else
                  l = 0
                  exit
               end if

               Y = Y + 1

               if (e > 0) then
                  X = X - 1

                  if (is_phase(X,Y,X1,Y1,X2,Y2)) then
                     l = 1
                  else
                     l = 0
                     exit
                  end if

                  e = e + 2*(dX - dY)

               else
                  e = e + 2*dX
               end if
            end do
         end if

      ! 3rd quadrant
      else if ((X2 < X1) .and. (Y2 < Y1)) then
         X = X1
         Y = Y1

         dX = X1 - X2
         dY = Y1 - Y2

         if (dX >= dY) then
            e = 3*dY - dX
            do i = 1, dX
               if (is_phase(X,Y,X1,Y1,X2,Y2)) then
                  l = 1
               else
                  l = 0
                  exit
               end if

               X = X - 1

               if (e > 0) then
                  Y = Y - 1

                  if (is_phase(X,Y,X1,Y1,X2,Y2)) then
                     l = 1
                  else
                     l = 0
                     exit
                  end if

                  e = e + 2*(dY - dX)

               else
                  e = e + 2*dY
               end if
            end do

         else if (dY > dX) then
            e = 3*dX - dY
            do i = 1, dY
               if (is_phase(X,Y,X1,Y1,X2,Y2)) then
                  l = 1
               else
                  l = 0
                  exit
               end if

               Y = Y - 1

               if (e > 0) then
                  X = X - 1

                  if (is_phase(X,Y,X1,Y1,X2,Y2)) then
                     l = 1
                  else
                     l = 0
                     exit
                  end if

                  e = e + 2*(dX - dY)

               else
                  e = e + 2*dX
               end if
            end do
         end if

      ! 4th quadrant
      else if ((X2 >= X1) .and. (Y2 < Y1)) then
         X = X1
         Y = Y1

         dX = X2 - X1
         dY = Y1 - Y2

         if (dX >= dY) then
            e = 3*dY - dX
            do i = 1, dX
               if (is_phase(X,Y,X1,Y1,X2,Y2)) then
                  l = 1
               else
                  l = 0
                  exit
               end if

               X = X + 1

               if (e > 0) then
                  Y = Y - 1

                  if (is_phase(X,Y,X1,Y1,X2,Y2)) then
                     l = 1
                  else
                     l = 0
                     exit
                  end if

                  e = e + 2*(dY - dX)

               else
                  e = e + 2*dY
               end if
            end do

         else if (dY > dX) then
            e = 3*dX - dY
            do i = 1, dY
               if (is_phase(X,Y,X1,Y1,X2,Y2)) then
                  l = 1
               else
                  l = 0
                  exit
               end if

               Y = Y - 1

               if (e > 0) then
                  X = X + 1

                  if (is_phase(X,Y,X1,Y1,X2,Y2)) then
                     l = 1
                  else
                     l = 0
                     exit
                  end if

                  e = e + 2*(dX - dY)

               else
                  e = e + 2*dX
               end if
            end do
         end if
      end if        

   end function line_sight

   subroutine calc_f_field(param)

      use images

      integer             :: b, p, i, j, ip, jp, num, X1, Y1, X2, Y2, e
      integer             :: ppos(2), bpos(2)
      integer, intent(in) :: param
      real                :: dist_pb
      real                :: e_pb(2)

      ! Calculate pseudo Coulomb field for each voxel of the phase
      allocate(pseudo_f(phas_numpix,2))

      ! Loop over phase points
      pseudo_f = 0.0
      do p = 1, phas_numpix
         ! Loop over boundary voxels
         do b = 1, nbound
           
            ! Calculate distance between phase and boundary voxels
            ppos(:) = map(phas_pospix(p),:)
            bpos(:) = map(phas_bound(b),:)
 
            dist_pb = norm2(real(ppos) - real(bpos))

            ! Calculate unitary vector between phase and boundary voxels
            e_pb = (real(ppos) - real(bpos))/dist_pb

            X1 = ppos(1)
            Y1 = ppos(2)
                
            X2 = bpos(1)
            Y2 = bpos(2)
        
            pseudo_f(p,:) = pseudo_f(p,:) + line_sight(X1,Y1,X2,Y2)*e_pb(:)/(dist_pb)**param

         end do
      end do

      ! Correct divergent boundary force points (those which point out of the
      ! phase boundary)
      do p = 1, phas_numpix
         i = map(phas_pospix(p),1)
         j = map(phas_pospix(p),2)

         if (pseudo_f(p,1) < -1e-3) then
            ip = i - 1
         else if (pseudo_f(p,1) > 1e-3) then
            ip = i + 1
         else
            ip = i
         end if

         if (pseudo_f(p,2) < -1e-3) then
            jp = j - 1
         else if (pseudo_f(p,2) > 1e-3) then
            jp = j + 1
         else
            jp = j
         end if

         ! Check if ip and jp are part of a different phase
         if (img(ip,jp) /= phase) then
            pseudo_f(p,:) = 0.0
         end if

      end do
         
   end subroutine calc_f_field

   subroutine calc_particles()

      use images

      integer              :: x, xp, num, i, j, ip, jp, t, npath, &
                              it, jt, nid, p, ind, pos
      integer, allocatable :: tmp(:), traj(:), tmp1(:), non_assigned(:)                    
      real                 :: r

      allocate(Pmat(nx,ny))
      allocate(nodes(1))
      allocate(non_assigned(phas_numpix))
      Pmat = 0
      non_assigned = 0
      npart = 0
      call random_seed()
      do ! Trajectory start loop
         ! Check which pixels of the phase are still unassigned with a particle
         ! number
         ind = 1
         nid = 0
         do p = 1, phas_numpix
            i = map(phas_pospix(p),1)
            j = map(phas_pospix(p),2)

            if (Pmat(i,j) == 0) then
               non_assigned(ind) = phas_pospix(p)
               nid = nid + 1
               ind = ind + 1
            end if
         end do

         if (nid == 0) exit

         call random_number(r)
         num = int(r*nid) + 1 
         x = non_assigned(num)

         i = map(x,1)
         j = map(x,2)

         Pmat(i,j) = npart + 1
         ! Initialize trajectory
         allocate(traj(1))
         traj(1) = x
        
         npath = 1
         do ! Trajectory evolution loop
            ! Current position x
            i = map(x,1)
            j = map(x,2)
        
            pos = maxval(findloc(phas_pospix,x))

            ! Next position x'
            if (pseudo_f(pos,1) < -1e-3) then
               ip = i - 1
            else if (pseudo_f(pos,1) > 1e-3) then
               ip = i + 1
            else
               ip = i
            end if

            if (pseudo_f(pos,2) < -1e-3) then
               jp = j - 1
            else if (pseudo_f(pos,2) > 1e-3) then
               jp = j + 1
            else
               jp = j
            end if

            xp = pixel_no(ip,jp)

            ! Is this a new point?
            if ((img(ip,jp) == phase) .and. &
                (Pmat(ip,jp) == 0)    .and. &
                (x /= xp))             then
               print*, 'Evolve' 
               Pmat(ip,jp) = npart + 1
               x = xp      
               ! Extend trajetory allocation here
               allocate(tmp1(npath + 1))
               tmp1(1:npath) = traj(:)
               call move_alloc(tmp1,traj)
               npath = npath + 1
               traj(npath) = xp
  
            ! Particle detected?
            else if ((Pmat(ip,jp) == npart + 1) .or. &
                     (xp == x))          then
               print*, 'Detection'
               npart = npart + 1
               nodes(npart) = x
               !Allocate space for node
               allocate(tmp(npart + 1))
               tmp(1:npart) = nodes(:)
               call move_alloc(tmp,nodes)              
 
               deallocate(traj)
               ! End LOOP
               exit

            ! Does this trajectory intersect with a previous one?
            else if ((img(ip,jp) == phase)      .and. &
                     (Pmat(ip,jp) /= 0)         .and. &
                     (Pmat(ip,jp) /= npart + 1)) then
               ! Merge particle number for the full trajectory
               do t = 1, npath
                  it = map(traj(t),1)
                  jt = map(traj(t),2)

                  Pmat(it,jt) = Pmat(ip,jp)
               end do

               deallocate(traj)
               ! End LOOP
               exit
               print*, 'Merging'
            else
               print*, 'Error', ip, jp
               stop               
            end if
         end do

      end do

   end subroutine calc_particles

end module pseudocoulomb
