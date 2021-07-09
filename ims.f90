module images

!This module reads segmented images for the calculation
!The format of the text must include a header with image dimensions,
!number of phases included and the pixels such that the different phases
!are marked as, for example, 0: Pore phase, 1: AM, ect.
!Border pixels of the image are set to number 255, so that no phase is the same
!This will help define particles that are cut, at the cost of losing the
!information about the image in the borders.

   implicit none
 
   ! Parameters 
   integer    , save, public               :: nphas, phase, nbound, nx, ny, phas_numpix 
   integer    , save, public, allocatable  :: img(:,:), bound(:,:)
   integer    , save, public, allocatable  :: map(:,:)
   integer    , save, public, allocatable  :: phas_pospix(:), phas_bound(:)
   real       , save, public               :: lvox

   ! Subroutines
   public     :: get_image
   public     :: get_phases
   public     :: get_boundary_voxels

   ! Functions
   public     :: pixel_no

   contains

   function pixel_no(indx,indy) result (pos)

      integer, intent(in)   :: indx, indy
      integer               :: pos

      pos = (indx-1)*ny + indy

   end function pixel_no

   subroutine get_input()

      logical               :: io_find

      inquire(file='input.in', exist=io_find)
      if (io_find) then
         open(unit=1, file='input.in', status='old')
            read(1,*) phase
            read(1,*) lvox
         close(1)
      else
         print*, 'No input.in file found'
         stop
      end if

   end subroutine get_input

   subroutine get_image()
      logical               :: io_find
      integer               :: i, j, ind1

      inquire(file='img.dat', exist=io_find)
      if (io_find) then
         open(unit=1, file='img.dat', status='old')
            read(1,*) nphas
            read(1,*) nx, ny
         
            ! Allocate space for array
            allocate(img(nx,ny))

            do j = 1, ny
               read(1,*) img(:,j)
            end do

         close(1)

         ! Transform pixels at borders       
         img(:,1)  = 10000
         img(:,ny) = 10000
         img(1,:)  = 10000
         img(nx,:) = 10000

      else
         print*, 'No img.dat file found'
         stop
      end if

      ! Construct map for i,j,k positions
      allocate(map(nx*ny,2))       

      ind1 = 1
      do i = 1, nx
         do j = 1, ny
            map(ind1,1) = i
            map(ind1,2) = j
             
            ind1 = ind1 + 1
         end do
      end do
   
   end subroutine get_image

   subroutine get_phases()

      integer               :: i, j, ind1, ind2               

      phas_numpix = 0
      do i = 1, nx
         do j = 1, ny
                  
            if (img(i,j) == phase) then
               phas_numpix = phas_numpix + 1
            end if

         end do
      end do

      allocate(phas_pospix(phas_numpix))

      ind1 = 1
      ind2 = 1
      do i = 1, nx
         do j = 1, ny

            if (img(i,j) == phase) then
               phas_pospix(ind2) = ind1
               ind2 = ind2 + 1
            end if

            ind1 = ind1 + 1

         end do
      end do
               
   end subroutine get_phases

   subroutine get_boundary_voxels()

      ! construct a list of the voxel indices that are adjacent
      ! but outside of the selected phase

      integer               :: ip, jp, p, b, ind1
      integer, allocatable  :: dummy(:)

      ! Count number of boundary voxels of selected phase, and record position
      allocate(dummy(nx*ny))
      dummy = 0      

      ind1 = 1
      nbound = 0
      do p = 1, phas_numpix
         ip = map(phas_pospix(p),1)
         jp = map(phas_pospix(p),2)

         if (img(ip + 1,jp) /= phase) then
            dummy(ind1) = pixel_no(ip + 1,jp)
            ind1 = ind1 + 1
            nbound = nbound + 1
         end if
         if (img(ip - 1,jp) /= phase) then
            dummy(ind1) = pixel_no(ip - 1,jp)
            ind1 = ind1 + 1
            nbound = nbound + 1
         end if
         if (img(ip,jp + 1) /= phase) then
            dummy(ind1) = pixel_no(ip,jp + 1)
            ind1 = ind1 + 1 
            nbound = nbound + 1
         end if
         if (img(ip,jp - 1) /= phase) then
            dummy(ind1) = pixel_no(ip,jp - 1)
            ind1 = ind1 + 1
            nbound = nbound + 1
         end if
         if (img(ip + 1,jp + 1) /= phase) then
            dummy(ind1) = pixel_no(ip + 1,jp + 1)
            ind1 = ind1 + 1
            nbound = nbound + 1
         end if
         if (img(ip - 1,jp + 1) /= phase) then
            dummy(ind1) = pixel_no(ip - 1,jp + 1)
            ind1 = ind1 + 1
            nbound = nbound + 1
         end if
         if (img(ip - 1,jp - 1) /= phase) then
            dummy(ind1) = pixel_no(ip - 1,jp - 1)
            ind1 = ind1 + 1
            nbound = nbound + 1
         end if
         if (img(ip + 1,jp - 1) /= phase) then
            dummy(ind1) = pixel_no(ip + 1,jp - 1)
            ind1 = ind1 + 1
            nbound = nbound + 1
         end if
      end do

      ! Now that we know how much voxels we have we can reduce the dummy array
      allocate(phas_bound(nbound))
      phas_bound(:) = dummy(1:nbound)
      deallocate(dummy)      

      allocate(bound(nx,ny))
      bound = 0
      do b = 1, nbound
         ip = map(phas_bound(b),1)
         jp = map(phas_bound(b),2)

         bound(ip,jp) = 1
      end do

   end subroutine get_boundary_voxels

end module images
