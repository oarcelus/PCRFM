module postprocess

   ! This module applies some post-processing steps to the obtained particle
   ! trajectories. Among them we apply the following operations
   ! 1- Contiguous centers are merged
   ! 2- Removal of chessboard pattern
   ! 3- Removal of single voxel particles
   ! 4- Eliminate particles that verify d-PSD < c-PSD (another module)

   implicit none   

   integer, save, public, allocatable    :: Pmat_pp(:,:), part_pospix(:,:)
   integer, save, public, allocatable    :: chess_pospix(:), part_numpix(:), part_centroids(:)
   integer, save, public                 :: npart_pp, nchess

   public :: merge_contiguous
   public :: detect_chessboard
   public :: get_particle_pospix
   public :: calc_centroids
   public :: relabel_chessboard

   contains

   subroutine merge_contiguous()

      use images
      use pseudocoulomb

      integer                :: i, j, inod1, jnod1, inod2, jnod2, p, ip, jp
      integer, allocatable   :: tmp(:)

      ! Merge voxels that belong to contiguous centers of attraction
      ! Everythig in the 2 voxel vicinity

      allocate(Pmat_pp(nx,ny))
      npart_pp = npart

      do i = 1, npart_pp - 1
         do j = i + 1, npart_pp

            inod1 = map(nodes(i),1)
            jnod1 = map(nodes(i),2)

            inod2 = map(nodes(j),1)
            jnod2 = map(nodes(j),2)

            if (abs(inod1 - inod2) <= 2 .and. abs(jnod1 - jnod2) <= 2) then

               ! Merge both particles
               do p = 1, phas_numpix
                  ip = map(phas_pospix(p),1)
                  jp = map(phas_pospix(p),2)

                  if (Pmat(ip,jp) == j) then
                     Pmat(ip,jp) = i
                  end if
               end do

               do p = 1, phas_numpix
                  ip = map(phas_pospix(p),1)
                  jp = map(phas_pospix(p),2)

                  if (Pmat(ip,jp) > j) then
                     Pmat(ip,jp) = Pmat(ip,jp) - 1
                  end if
               end do

               npart_pp = npart_pp - 1

               ! Reallocate nodes array
               allocate(tmp(npart_pp))
               tmp(1:j-1) = nodes(1:j-1)
               tmp(j:npart_pp) = nodes(j+1:npart_pp+1)
               call move_alloc(tmp,nodes)

            end if
         end do
      end do

      Pmat_pp = Pmat

      deallocate(tmp)
      deallocate(Pmat)

   end subroutine merge_contiguous

   subroutine get_particle_pospix()

      use images
      use pseudocoulomb

      integer              :: part, p, ip, jp, countp, valmax
      integer, allocatable :: tmp(:,:)

      allocate(part_pospix(npart_pp,phas_numpix))
      allocate(part_numpix(npart_pp))

      part_pospix = 0
      part_numpix = 0
      countp = 0
      do part = 1, npart_pp
         do p = 1, phas_numpix
            ip = map(phas_pospix(p),1)
            jp = map(phas_pospix(p),2)
            
            if (Pmat_pp(ip,jp) == part) then
               countp = countp + 1
               part_pospix(part,countp) = phas_pospix(p)
            end if 
            
         end do

         part_numpix(part) = countp
         countp = 0
      end do 

      allocate(tmp(npart_pp,maxval(part_numpix)))
      tmp(:,1:maxval(part_numpix)) = part_pospix(:,1:maxval(part_numpix))
      call move_alloc(tmp,part_pospix)
      deallocate(tmp)

   end subroutine get_particle_pospix

   subroutine detect_chessboard()

      use images
      use pseudocoulomb

      integer              :: p, ip, jp, part
      integer, allocatable :: tmp(:)
      logical              :: chess

      allocate(chess_pospix(phas_numpix))     
      chess_pospix = 0
      nchess = 0

      ! Detect chessboard pattern
      do p = 1, phas_numpix
         ip = map(phas_pospix(p),1)
         jp = map(phas_pospix(p),2)

         !Create truth condition
         part = Pmat_pp(ip,jp)
         chess = (Pmat_pp(ip+1,jp) /= part) .and. (Pmat_pp(ip-1,jp) /= part)     .and. &
                 (Pmat_pp(ip,jp+1) /= part) .and. (Pmat_pp(ip,jp-1) /= part)     !.and. &
                 !(Pmat_pp(ip+1,jp+1) == part) .and. (Pmat_pp(ip+1,jp-1) == part) .and. &
                 !(Pmat_pp(ip-1,jp+1) == part) .and. (Pmat_pp(ip-1,jp-1) == part)

         if (chess) then
            nchess = nchess + 1
            chess_pospix(nchess) = phas_pospix(p)     
         end if

      end do
     
      ! Reallocate nodes array
      allocate(tmp(nchess))
      tmp(1:nchess) = chess_pospix(1:nchess)
      call move_alloc(tmp,chess_pospix)
      deallocate(tmp)

   end subroutine detect_chessboard

   subroutine calc_centroids()

      use images
      use pseudocoulomb

      integer                 :: part, p, ip, jp, ic, jc, m00, m01, m10, xc

      allocate(part_centroids(npart_pp))

      !Calculate centroids only if not part of chessboard pattern
      !calculate moments m10, m01, m00
      do part = 1, npart_pp
         m00 = 0
         m01 = 0
         m10 = 0
         do p = 1, part_numpix(part)
            ip = map(part_pospix(part,p),1)
            jp = map(part_pospix(part,p),2)

            if (.not.any(chess_pospix == part_pospix(part,p))) then 
               m00 = m00 + Pmat_pp(ip,jp)
               m10 = m10 + ip*Pmat_pp(ip,jp)
               m01 = m01 + jp*Pmat_pp(ip,jp)
            end if

         end do

         ic = int(m10/m00)
         jc = int(m01/m00)

         xc = pixel_no(ic,jc)

         part_centroids(part) = xc

      end do

   end subroutine calc_centroids

   subroutine relabel_chessboard()

      use images
      use pseudocoulomb

      integer                :: chess, ich, jch, ic, jc, part, minpart      
      real*8, allocatable    :: dist(:)
      real*8                 :: dist2

      allocate(dist(npart_pp))

      !Relabel pixels belonging to chessboard pattern to closests centroid label
      do chess = 1, nchess
         ich = map(chess_pospix(chess),1)
         jch = map(chess_pospix(chess),2)
         do part = 1, npart_pp
            ic = map(part_centroids(part),1)
            jc = map(part_centroids(part),2)
            
            dist2 = real((ic-ich)**2 + (jc-jch)**2)

            dist(part) = sqrt(dist2)                        
         end do

         minpart = minloc(dist, DIM=1)

         Pmat_pp(ich,jch) = minpart

      end do

   end subroutine relabel_chessboard

end module postprocess
