program main 

   use images
   use pseudocoulomb
   use postprocess

   implicit none

   integer                           :: i, j, p
   integer, allocatable              :: phase_ims(:,:), bound_ims(:,:)

   open(unit=33, file='std.out', status='replace')
      write(33,'(a)') '!=======================================================!'
      write(33,'(a)') '!          PPPP    CCCC RRRR  FFFFF MM     MM           !'
      write(33,'(a)') '!          P   P  C     R   R F     M M   M M           !'
      write(33,'(a)') '!          P   P C      R   R F     M  M M  M           !'
      write(33,'(a)') '!          PPPP  C      RRRR  FFFF  M   M   M           !'
      write(33,'(a)') '!          P     C      RR    F     M       M           !'
      write(33,'(a)') '!          P     C      R R   F     M       M           !'
      write(33,'(a)') '!          P      C     R  R  F     M       M           !'
      write(33,'(a)') '!          P       CCCC R   R F     M       M           !'
      write(33,'(a)') '!                                                       !'
      write(33,'(a)') '! This program executes the Pseudo-Coulomb Repulsive-   !'
      write(33,'(a)') '! Field Method, for pore or particle identification,    !'
      write(33,'(a)') '! and post-processing analysis.                         !'
      write(33,'(a)') '!                                                       !'
      write(33,'(a)') '! Author : Oier Arcelus                                 !'
      write(33,'(a)') '! Version: 0.0                                          !'
      write(33,'(a)') '! Date   : July 2020                                    !'
      write(33,'(a)') '!=======================================================!'
      write(33,'(a)') ' '


      write(33,'(a)') 'Reading input file...'
      write(33,'(a)') ' '
      !Read input file
      call get_input()    
      write(33,'(a,i2)')    'Phase to study : ', phase
      write(33,'(a,f12.6)') 'Voxel lenght   : ', lvox
      write(33,'(a)') ' '

      write(33,'(a)') 'Reading image...'
      !Read image
      call get_image() 
      write(33,'(a,i2)')    'Total No of phases : ', nphas
      write(33,'(a,3i)')    'Image dimension    : ', nx, ny
      write(33,'(a)') ' '

      allocate(phase_ims(nx,ny))
      allocate(bound_ims(nx,ny))

      write(33,'(a)') 'Getting information about phases and boundary voxels...'
      write(33,'(a)') ' '
      ! Read phases and boundary voxels
      call get_phases()
      call get_boundary_voxels()
      write(33,'(a)') 'Writting phase information to file...'
      
      ! Phase information
      phase_ims = 0
      do p = 1, phas_numpix
         i = map(phas_pospix(p),1)
         j = map(phas_pospix(p),2)
         phase_ims(i,j) = 255
      end do
 
      open(unit = 1, file = 'phase.dat', status = 'replace')
         write(1,'(2i)') nx, ny
         do i = 1, nx
            do j = 1, ny
               write(1,'(i)') phase_ims(i,j)
            end do
         end do
      close(1)

      write(33,'(a)') 'Writting boundary information to file...'

      ! Boundary information
      bound_ims = 0
      do p = 1, nbound
         i = map(phas_bound(p),1)
         j = map(phas_bound(p),2)
         bound_ims(i,j) = 255
      end do

      open(unit = 1, file = 'bound.dat', status = 'replace')
         write(1,'(2i)') nx, ny
         do i = 1, nx
            do j = 1, ny
               write(1,'(i)') bound_ims(i,j)
            end do
         end do
      close(1)


      write(33,'(a)') 'Computing pseudo coulomb field...'
      write(33,'(a)') ' '
      ! Compute pseudo-coulomb forces
      call calc_f_field(2)

      open(unit = 1, file = 'force.dat', status = 'replace')
         do p = 1, phas_numpix
            i = map(phas_pospix(p),1)
            j = map(phas_pospix(p),2)
            write(1,'(2i,2F12.6)') i, j, pseudo_f(p,1), pseudo_f(p,2)
         end do
      close(1)

      write(33,'(a)') 'Detecting particles...'
      write(33,'(a)') ' '
      call calc_particles()

      write(33,'(a)') 'Starting postprocessing...'
      write(33,'(a)') 'Mergin voxels of contiguous particles...'
      call merge_contiguous()
      write(33,'(a,i2)')    'Total particles detected : ', npart_pp
      write(33,'(a)') ' '
      call get_particle_pospix()
      write(33,'(a)') 'Detecting chessboard pattern pixels...'
      call detect_chessboard() 
      write(33,'(a)') 'Chessboard pattern detected'
      write(33,'(a)') ' '
      write(33,'(a)') 'Calculating particle centroids...'
      call calc_centroids()
      write(33,'(a)') 'Particle centroids calculated'
      write(33,'(a)') ' '
      write(33,'(a)') 'Relabeling chessboard pixels...'
      call relabel_chessboard()
      write(33,'(a)') 'Relabeling done'

      open(unit = 1, file = 'particle.dat', status = 'replace')
         write(1,'(2i)') nx, ny
         write(1,'(i)') npart_pp
         do i = 1, nx
            do j = 1, ny
               write(1,'(i)') Pmat_pp(i,j)
            end do
         end do
      close(1)
   close(33)

end program main
