MODULE loc_aper
    !-----------------------------------------------------
    ! Author:  	Lily Barta
    ! Contact: 	lily.barta@etu.u-bordeaux.fr
    ! 
    ! Supervisior: PD Dr. Denis Usvyat
    ! Contact:     denis.usvyat@hu-berlin.de
    !-----------------------------------------------------
    
    USE mp2mem
    
    PUBLIC
      
      character(len=50) :: ContextName
      integer :: ierr
      
      integer :: n_at_dom
      integer, pointer :: min_domain_list(:), min_domain_start(:) 
      integer, pointer :: domain_list(:), domain_start(:), domain_end(:)
      integer, pointer :: pair_list(:,:)
      integer :: n_pairs, n_strong,  n_weak, n_dist
      
      integer :: n_virt ! number of PAOs after removing the redundant ones
      real*8, pointer :: Fock_local_occ(:,:), Fock_local_virt(:,:), S_paos(:,:)
      real*8, pointer :: Transformation_Matrix_W(:,:), Transformation_Matrix_Q(:,:), Transformation_Matrix_V(:,:)
      real*8, pointer :: Int_p_ai(:,:,:), Int_p_ji(:,:,:), Int_p_ab(:,:,:)
      real*8, pointer :: Coulomb_Metrik_J_PQ(:,:)
    
    CONTAINS   
    SUBROUTINE my_loc(n_bas_cas, n_occ_fcidump, C_loc_mu_i)
    ! Orbital localization algorithm of Pipek and Mezey
    ! (J. Chem. Phys. 90, 4916(1989))
    ! with rotation angles determined as described by 
    ! Edmiston and Ruedenberg (Rev. Mod. Phys. 35, 457 (1963))
    
      USE fci_interface, ONLY : n_at_cas, S_proj, cas_mul_map   ! number of atoms, AO overlap matrix, map of atoms in the fragment
      
      IMPLICIT NONE
     
      integer, intent(in) :: n_bas_cas	! number of atomic orbitals
      integer, intent(in) :: n_occ_fcidump	! number of molecular orbitals
      real*8, intent(inout), dimension(n_bas_cas,n_occ_fcidump) :: C_loc_mu_i   ! occupied LMO coefficients 
      
      integer :: i_at, i_ao 
      integer :: i_mo, j_mo, k_mo
      integer :: i_opt, j_opt      ! index of MOs to be rotated
      integer :: n_mo_pairs	! number of MO pairs
      integer :: n_opt_pair	! index of the MO pair to be rotated
      integer :: counter
      integer :: at_in_cell
      integer :: n_iter_PM, n_max_iter_PM  
      double precision :: conv_thresh_PM
      real*8, pointer :: P_A_intermediate(:,:), P_A(:,:), delta_L(:), A_ij(:), B_ij(:), P_ii(:)
      double precision :: cos4a, sin4a, xx, x, y, test, cosa, sina, thresh, ci, cj, delta_L_opt, Q_ii_sq, Q_ii
       
      ContextName = 'my_loc'
       
      write(*,*) ' '
      write(*,*) '*******************************************************************************'
      write(*,*) '*                                                                             *'
      write(*,*) '*                          PIPEK MEZEY LOCALIZATION                           *'
      write(*,*) '*                                                                             *'
      write(*,*) '*******************************************************************************'
      write(*,*) ' ' 
      write(*,'(A30,I4)') 'Number of molecular orbitals: ', n_occ_fcidump
      write(*,'(A27,I4)') 'Number of atomic orbitals: ', n_bas_cas
      write(*,'(A17,I4)') 'Number of atoms: ', n_at_cas
       
      !call matprt(C_loc_mu_i,n_bas_cas,n_occ_fcidump,'Canonical orbitals')
       
      !-----------------------------------------------------------------------------
      !         FIND MO PAIR (i,j) WHICH YIELDS THE
      !         GREATEST INCREASE IN LOCALIZATION (LARGEST delta_L) 
      !-----------------------------------------------------------------------------
      ! 
      n_mo_pairs = n_occ_fcidump*(n_occ_fcidump-1)/2
      call mp2alloc(P_A_intermediate, 2, n_bas_cas, ContextName, ierr, 'P_A_intermediate')
      call mp2alloc(P_A, 2, 2, ContextName, ierr, 'P_A')
      call mp2alloc(delta_L, n_mo_pairs, ContextName, ierr, 'delta_L')
      call mp2alloc(A_ij, n_mo_pairs, ContextName, ierr, 'A_ij')
      call mp2alloc(B_ij, n_mo_pairs, ContextName, ierr, 'B_ij')
      call mp2alloc(P_ii, n_bas_cas, ContextName, ierr, 'P_ii')
      ! 
      !------- FOR EACH MO PAIR (i,j) ----------------------------------------------
      ! 
      counter = 1 
      do i_mo = 1, n_occ_fcidump-1
        do j_mo = i_mo+1, n_occ_fcidump
          ! 
          !----- COMPUTE A_ij AND B_ij ---------------------------------------------
          ! 
          A_ij(counter) = 0.0d0
          B_ij(counter) = 0.0d0
          call mxma(C_loc_mu_i(1,i_mo), (j_mo-i_mo)*n_bas_cas, 1, &
                     S_proj, 1, n_bas_cas, &
                      P_A_intermediate, 1, 2, &
                      2, n_bas_cas, n_bas_cas)
          do i_at = 1, n_at_cas
            call mxma(P_A_intermediate(1,cas_mul_map(i_at)), 1, 2, &
                      C_loc_mu_i(cas_mul_map(i_at),i_mo), 1, (j_mo-i_mo)*n_bas_cas, &
                      P_A(1,1), 1, 2, &
                      2, cas_mul_map(i_at+1)-cas_mul_map(i_at), 2)
            A_ij(counter) = A_ij(counter) + 0.25d0*((P_A(1,2)+P_A(2,1))**2.0d0 - (P_A(1,1)-P_A(2,2))**2.0d0)
            B_ij(counter) = B_ij(counter) + 0.5d0*(P_A(1,2)+P_A(2,1))*(P_A(1,1)-P_A(2,2))
          end do
          ! 
          !----- COMPUTE delta_L -----------------------------------------------------
          ! 
          delta_L(counter) = A_ij(counter) + dsqrt(A_ij(counter)**2+B_ij(counter)**2)
          counter = counter + 1 
        end do
      end do
      ! 
      !----- FIND LARGEST delta_L AND CORRSPONDING MO PAIR (i,j) ---------------------
      !      (only the number of the optimal mo pair)
      ! 
      n_opt_pair = 1
      do i_mo = 2, n_mo_pairs
        if (delta_L(i_mo).gt.delta_L(n_opt_pair)) then
          n_opt_pair = i_mo
        end if
      end do
      delta_L_opt = delta_L(n_opt_pair)
      ! 
      !*********** Begin convergence loop *********************************************
      ! 
      n_iter_PM = 0
      n_max_iter_PM = n_occ_fcidump**3
      conv_thresh_PM = 1.0d-10
      do while (conv_thresh_PM.lt.delta_L_opt)
      ! 
      ! 
      n_iter_PM = n_iter_PM +1
      if (n_iter_PM.gt.n_max_iter_PM) then
        write(*,*)
        write(*,'(A21,I4,A11)')'No convergence after ',n_iter_PM,' iterations'
        write(*,'(A,E10.3)')'Localization convergence threshold: ', conv_thresh_PM
        write(*,'(A,E10.3)')'Max change in localization: ', delta_L_opt
        stop
      end if
      ! 
      !--------------------------------------------------------------------------------
      !         FIND THE CORRESPONDING TRANSFORMATION
      !         AND THUS THE TWO NEW ORBITALS 
      !--------------------------------------------------------------------------------
      ! 
      !------ COMPUTE cos(4alpha) AND sin(4alpha) ------------------------------------- 
      ! 
      cos4a =-A_ij(n_opt_pair)/dsqrt(A_ij(n_opt_pair)**2+B_ij(n_opt_pair)**2)
      sin4a = B_ij(n_opt_pair)/dsqrt(A_ij(n_opt_pair)**2+B_ij(n_opt_pair)**2)
      ! 
      !------ COMPUTE THE TWO PAIRS (x1,y1) AND (x2,y2) -------------------------------
      !       KEEP THE PAIR WHICH SATISFIES: 
      !       4xy(x**2-y**2) = sin(4alpha)
      ! 
      xx = 0.5d0*(1.0d0+dsqrt(1.0d0-0.5d0*(1.0d0-cos4a)))
      x = dsqrt(xx)
      y = dsqrt(1-xx)
      test = 4.0d0*x*y*(2.0d0*xx-1.0d0)
      thresh = 1.0d-8
      if (dabs(sin4a-test).lt.thresh) then
        cosa = x
        sina = y
      else
        xx = 0.5d0*(1.0d0-dsqrt(1.0d0-0.5d0*(1.0d0-cos4a)))
        test = 4.0d0*x*y*(2.0d0*xx-1.0d0)
        if (dabs(sin4a-test).lt.thresh) then
          cosa = y
          sina = x
        else
          write(*,*) 'Problem: no valid rotation angle found'
          stop
        end if
      end if
      ! 
      !------ COMPUTE (i_opt,j_opt) FROM n_opt_pair -----------------------------------
      ! 
      counter = 0
      do k_mo = 1, n_occ_fcidump-1
        counter = counter + n_occ_fcidump - k_mo
        if (n_opt_pair.le.counter) then
          i_opt = k_mo
          exit
        end if
      end do
      j_opt = n_opt_pair - (i_opt-1)*n_occ_fcidump + i_opt*(i_opt+1)/2
      ! 
      !----- ROTATE THE CHOSEN MO PAIR (i_opt,j_opt) ----------------------------------
      ! 
      do i_ao = 1, n_bas_cas
        ci = C_loc_mu_i(i_ao,i_opt)
        cj = C_loc_mu_i(i_ao,j_opt)
        C_loc_mu_i(i_ao,i_opt) = cosa*ci + sina*cj
        C_loc_mu_i(i_ao,j_opt) =-sina*ci + cosa*cj
      end do
      ! 
      !--------------------------------------------------------------------------------
      !         COMPUTE NEW A_ij, B_ij AND delta_L ONLY FOR 
      !         THE PAIRS AFFECTED BY THE PREVIOUS ROTATION 
      !--------------------------------------------------------------------------------
      ! 
       do k_mo = 1, n_occ_fcidump
         if (k_mo.lt.i_opt) then 
           i_mo = k_mo 
           j_mo = i_opt
           counter = (i_mo-1)*n_occ_fcidump - i_mo*(i_mo+1)/2 + j_mo
           A_ij(counter) = 0.0d0
           B_ij(counter) = 0.0d0
           call mxma(C_loc_mu_i(1,i_mo), (j_mo-i_mo)*n_bas_cas, 1, &
                 S_proj, 1, n_bas_cas, &
                 P_A_intermediate, 1, 2, &
                 2, n_bas_cas, n_bas_cas)
           do i_at = 1, n_at_cas
             call mxma(P_A_intermediate(1,cas_mul_map(i_at)), 1, 2, &
                       C_loc_mu_i(cas_mul_map(i_at),i_mo), 1, (j_mo-i_mo)*n_bas_cas, &
                       P_A(1,1), 1, 2, &
                       2, cas_mul_map(i_at+1)-cas_mul_map(i_at), 2)
             A_ij(counter) = A_ij(counter) + 0.25d0*((P_A(1,2)+P_A(2,1))**2 - (P_A(1,1)-P_A(2,2))**2)
             B_ij(counter) = B_ij(counter) + 0.5d0*(P_A(1,2)+P_A(2,1))*(P_A(1,1)-P_A(2,2))
           end do
           delta_L(counter) = A_ij(counter) + dsqrt(A_ij(counter)**2+B_ij(counter)**2)
         else if (k_mo.gt.i_opt) then
           i_mo = i_opt
           j_mo = k_mo
           counter = (i_mo-1)*n_occ_fcidump - i_mo*(i_mo+1)/2 + j_mo
           A_ij(counter) = 0.0d0
           B_ij(counter) = 0.0d0
           call mxma(C_loc_mu_i(1,i_mo), (j_mo-i_mo)*n_bas_cas, 1, &
                 S_proj, 1, n_bas_cas, &
                 P_A_intermediate, 1, 2, &
                 2, n_bas_cas, n_bas_cas)
           do i_at = 1, n_at_cas
             call mxma(P_A_intermediate(1,cas_mul_map(i_at)), 1, 2, &
                       C_loc_mu_i(cas_mul_map(i_at),i_mo), 1, (j_mo-i_mo)*n_bas_cas, &
                       P_A(1,1), 1, 2, &
                       2, cas_mul_map(i_at+1)-cas_mul_map(i_at), 2)
             A_ij(counter) = A_ij(counter) + 0.25d0*((P_A(1,2)+P_A(2,1))**2 - (P_A(1,1)-P_A(2,2))**2)
             B_ij(counter) = B_ij(counter) + 0.5d0*(P_A(1,2)+P_A(2,1))*(P_A(1,1)-P_A(2,2))
           end do
           delta_L(counter) = A_ij(counter) + dsqrt(A_ij(counter)**2+B_ij(counter)**2)
         end if
         if (k_mo.lt.j_opt) then 
           i_mo = k_mo 
           j_mo = j_opt
           counter = (i_mo-1)*n_occ_fcidump - i_mo*(i_mo+1)/2 + j_mo
           A_ij(counter) = 0.0d0
           B_ij(counter) = 0.0d0
           call mxma(C_loc_mu_i(1,i_mo), (j_mo-i_mo)*n_bas_cas, 1, &
                 S_proj, 1, n_bas_cas, &
                 P_A_intermediate, 1, 2, &
                 2, n_bas_cas, n_bas_cas)
           do i_at = 1, n_at_cas
             call mxma(P_A_intermediate(1,cas_mul_map(i_at)), 1, 2, &
                       C_loc_mu_i(cas_mul_map(i_at),i_mo), 1, (j_mo-i_mo)*n_bas_cas, &
                       P_A(1,1), 1, 2, &
                       2, cas_mul_map(i_at+1)-cas_mul_map(i_at), 2)
             A_ij(counter) = A_ij(counter) + 0.25d0*((P_A(1,2)+P_A(2,1))**2 - (P_A(1,1)-P_A(2,2))**2)
             B_ij(counter) = B_ij(counter) + 0.5d0*(P_A(1,2)+P_A(2,1))*(P_A(1,1)-P_A(2,2))
           end do
           delta_L(counter) = A_ij(counter) + dsqrt(A_ij(counter)**2+B_ij(counter)**2)
         else if (k_mo.gt.j_opt) then
           i_mo = j_opt
           j_mo = k_mo
           counter = (i_mo-1)*n_occ_fcidump - i_mo*(i_mo+1)/2 + j_mo
           A_ij(counter) = 0.0d0
           B_ij(counter) = 0.0d0
           call mxma(C_loc_mu_i(1,i_mo), (j_mo-i_mo)*n_bas_cas, 1, &
                 S_proj, 1, n_bas_cas, &
                 P_A_intermediate, 1, 2, &
                 2, n_bas_cas, n_bas_cas)
           do i_at = 1, n_at_cas
             call mxma(P_A_intermediate(1,cas_mul_map(i_at)), 1, 2, &
                       C_loc_mu_i(cas_mul_map(i_at),i_mo), 1, (j_mo-i_mo)*n_bas_cas, &
                       P_A(1,1), 1, 2, &
                       2, cas_mul_map(i_at+1)-cas_mul_map(i_at), 2)
             A_ij(counter) = A_ij(counter) + 0.25d0*((P_A(1,2)+P_A(2,1))**2 - (P_A(1,1)-P_A(2,2))**2)
             B_ij(counter) = B_ij(counter) + 0.5d0*(P_A(1,2)+P_A(2,1))*(P_A(1,1)-P_A(2,2))
           end do
           delta_L(counter) = A_ij(counter) + dsqrt(A_ij(counter)**2+B_ij(counter)**2)
         end if
       end do
       ! 
       !------- FIND LARGEST delta_L AND CORRSPONDING MO PAIR (i,j) ---------------------
       !        (only the number of the optimal mo pair)
       ! 
       n_opt_pair = 1
       do i_mo = 2, n_mo_pairs
         if (delta_L(i_mo).gt.delta_L(n_opt_pair)) then
           n_opt_pair = i_mo
         end if
       end do
       delta_L_opt = delta_L(n_opt_pair)
       ! 
       !*********** End convergence loop ************************************************
       end do 
       ! 
       !call matprt(C_loc_mu_i,n_bas_cas,n_occ_fcidump,'Localized orbitals')
       write(*,*)
       write(*,'(A,I4,A12)')'Convergence reached after ',n_iter_PM,' iterations'
       write(*,'(A,E10.3)')'Localization convergence threshold: ', conv_thresh_PM
       write(*,'(A,E10.3)')'Max change in localization: ', delta_L_opt
       write(*,*)
       ! 
       call mp2dealloc(P_A_intermediate, 'P_A_intermediate')
       call mp2dealloc(P_A, 'P_A')
       call mp2dealloc(delta_L, 'delta_L')
       call mp2dealloc(A_ij, 'A_ij')
       call mp2dealloc(B_ij, 'B_ij')
       call mp2dealloc(P_ii, 'P_ii')    
              
   END SUBROUTINE my_loc
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   SUBROUTINE Frag_PAOs(n_bas_cas,n_occ_fcidump,density_mat,C_loc_mu_i,C_paos)
   
       USE fci_interface, ONLY : S_proj  ! AO overlap matrix
       ! 
       IMPLICIT NONE
       integer, intent(in) :: n_bas_cas		! number of AOs
       integer, intent(in) :: n_occ_fcidump	! number of MOs
       real*8, intent(in), dimension(n_bas_cas,n_bas_cas) :: density_mat		
       real*8, intent(in), dimension(n_bas_cas,n_occ_fcidump) :: C_loc_mu_i	! occupied LMO coefficients
       real*8, intent(inout), dimension(n_bas_cas,n_bas_cas) :: C_paos(:,:)	! virtual LMO coefficients (including redondant MOs)
       ! 
       integer :: i_ao, i_pao, j_pao
       real*8, pointer :: temp_mat(:,:)
       ! 
       !------ BUILD IDENTITY MATRIX --------------------------------------
       ! 
       C_paos = 0.0d0
       do i_ao = 1, n_bas_cas
         C_paos(i_ao,i_ao) = 1.0d0
       end do
       ! 
       !------ COMPUTE PAO COEFFICIENTS (I-DS) ----------------------------
       ! 
       call mxmbn(density_mat, 1, n_bas_cas, &
             S_proj, 1, n_bas_cas, &
             C_paos, 1, n_bas_cas, &
             n_bas_cas, n_bas_cas, n_bas_cas)
   
       !call matprt(C_paos,n_bas_cas,n_bas_cas,'PAO coefficients')
       ! 
       !------- BUILD PAO OVERLAP MATRIX ---------------------------------------------
       ! 
       call mp2alloc(temp_mat, n_bas_cas, n_bas_cas, ContextName, ierr, 'temp_mat')
       call mp2alloc(S_paos, n_bas_cas, n_bas_cas, ContextName, ierr, 'S_paos')
       
       call mxma(S_proj, 1, n_bas_cas, &
                 C_paos, 1, n_bas_cas, &
                 temp_mat, 1, n_bas_cas, &
                 n_bas_cas, n_bas_cas, n_bas_cas)
       
       call mxma(C_paos, n_bas_cas, 1, &
                 temp_mat, 1, n_bas_cas, &
                 S_paos, 1, n_bas_cas, &
                 n_bas_cas, n_bas_cas, n_bas_cas)
       
       call mp2dealloc(temp_mat,'temp_mat')
       ! 
       !------ NORMALIZE PAOs ---------------------------------------------
       ! 
       do i_pao = 1, n_bas_cas
         do j_pao = 1, n_bas_cas
           C_paos(j_pao,i_pao) = C_paos(j_pao,i_pao)/dsqrt(S_paos(i_pao,i_pao))
         end do
       end do
       !call matprt(C_paos,n_bas_cas,n_bas_cas,'normalized PAO coefficients')
       call mp2dealloc(S_paos,'S_paos')
   
   END SUBROUTINE Frag_PAOs
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   SUBROUTINE Frag_domains(n_bas_cas,n_occ_fcidump,C_loc_mu_i)
         
       USE fci_interface,  ONLY : n_at_cas, cas_mul_map, ats_cas, ats_cas_cells
       USE basato_module,  ONLY : xa, nat 
       USE gvect_module,   ONLY : xg
       USE mp2_fixindex,   ONLY : rext, fac_rcov, iext, max_dom_ats
       ! 
       IMPLICIT NONE
       integer, intent(in) :: n_bas_cas, n_occ_fcidump
       real*8, intent(in), dimension(n_bas_cas,n_occ_fcidump) :: C_loc_mu_i
       ! 
       integer :: i, j, k
       integer :: i_at, i_at_dom, i_at_frag, i_cell_frag
       integer :: i_mo, i_ao
       integer :: at, at_dom
       integer :: min_dom_list_size, max_n_at_dom, temp_list_size
       double precision :: dist, dist_cov, rcov_at, rcov_at_dom
       double precision :: thresh_min_domain 
       double precision :: x_at(3), x_at_dom(3)
       integer, pointer :: domain_map(:,:)
       integer, pointer :: temp_start(:), temp_list(:)
       integer, pointer :: new_temp_start(:), new_temp_list(:)
       REAL(FLOAT) :: RADCOV(0:92)
       REAL(FLOAT) :: RAYDAT(0:92)
       DATA RAYDAT/0._FLOAT,.68_FLOAT,1.47_FLOAT,1.65_FLOAT,1.18_FLOAT,0.93_FLOAT,&
       0.81_FLOAT,0.78_FLOAT,0.78_FLOAT,0.76_FLOAT,1.68_FLOAT,2.01_FLOAT, &
       1.57_FLOAT,1.50_FLOAT,1.23_FLOAT,1.15_FLOAT,1.09_FLOAT,1.05_FLOAT, &
       1.97_FLOAT,2.31_FLOAT,2.07_FLOAT,1.68_FLOAT,1.47_FLOAT,1.41_FLOAT, &
       1.47_FLOAT,1.47_FLOAT,1.47_FLOAT,1.41_FLOAT,1.41_FLOAT,1.41_FLOAT, &
       1.41_FLOAT,1.36_FLOAT,1.31_FLOAT,1.21_FLOAT,1.21_FLOAT,1.21_FLOAT, &
       2.10_FLOAT,2.31_FLOAT,2.10_FLOAT,1.94_FLOAT,1.63_FLOAT,1.52_FLOAT, &
       1.52_FLOAT,1.42_FLOAT,1.36_FLOAT,1.42_FLOAT,1.47_FLOAT,1.68_FLOAT, &
       1.62_FLOAT,1.62_FLOAT,1.52_FLOAT,1.52_FLOAT,1.47_FLOAT,1.47_FLOAT, &
       2.66_FLOAT,2.73_FLOAT,2.10_FLOAT,1.84_FLOAT,1.63_FLOAT,1.63_FLOAT, &
       1.63_FLOAT,1.63_FLOAT,1.63_FLOAT,1.63_FLOAT,1.63_FLOAT,1.63_FLOAT, &
       1.63_FLOAT,1.63_FLOAT,1.63_FLOAT,1.63_FLOAT,1.63_FLOAT,1.63_FLOAT, &
       1.63_FLOAT,1.52_FLOAT,1.42_FLOAT,1.42_FLOAT,1.36_FLOAT,1.42_FLOAT, &
       1.42_FLOAT,1.42_FLOAT,1.57_FLOAT,1.99_FLOAT,1.89_FLOAT,1.68_FLOAT, &
       1.42_FLOAT,1.42_FLOAT,1.62_FLOAT,2.94_FLOAT,1.68_FLOAT,2.05_FLOAT, &
       1.63_FLOAT,1.63_FLOAT,1.63_FLOAT/
       ! 
       ContextName = 'Frag_domains'
       ! 
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       !               Building the minimal domains
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! 
       !---- BUILD MINIMAL DOMAIN MAP ---------------------------------------------
       ! 
       call mp2alloc(domain_map, n_at_cas, n_occ_fcidump, ContextName, ierr, 'domain_map')
       thresh_min_domain = 0.1d0
       domain_map = 0
       min_dom_list_size = 0
       ! 
       do i_mo = 1, n_occ_fcidump
         do i_at = 1, n_at_cas
           do i_ao = cas_mul_map(i_at), cas_mul_map(i_at+1)-1
             if (dabs(C_loc_mu_i(i_ao,i_mo)).gt.thresh_min_domain) then
               domain_map(i_at,i_mo) = 1
               min_dom_list_size = 1 + min_dom_list_size
               exit
             end if
           end do
         end do
       end do
       ! 
       !---- PRINT MINIMAL DOMAIN MAP ---------------------------------------------
       ! 
       !write(*,*)
       !write(*,'(A)') 'domain_map'
       !write(*,'(5X, *(I4,1X))')(j, j=1, n_occ_fcidump)
       !write(*,'(5X, *(A5))')('-----', j=1, n_occ_fcidump)
       !do i = 1, n_at_cas
       !  write(*,'(I3," |", *(I4,1X))') i, (domain_map(i,j), j=1, n_occ_fcidump)
       !end do
       ! 
       !---- CONVERT DOMAIN MAP TO min_dom_list -----------------------------------
       ! 
       call mp2alloc(min_domain_list, min_dom_list_size, ContextName, ierr, 'min_domain_list')
       call mp2alloc(min_domain_start, n_occ_fcidump+1, ContextName, ierr, 'min_domain_start')
       call izero(min_domain_list, min_dom_list_size)
       call izero(min_domain_start, n_occ_fcidump+1)
       ! 
       n_at_dom = 1
       min_domain_start(1) = 1
       ! 
       do i_mo = 1, n_occ_fcidump
         do i_at = 1, n_at_cas
           if (domain_map(i_at,i_mo).eq.1) then
             min_domain_list(n_at_dom) = i_at
             n_at_dom = n_at_dom + 1
           end if
         end do
         min_domain_start(i_mo+1) = n_at_dom
       end do
       !
       !----- PRINT min_domain_list AND min_domain_start -------------------------
       !
       !write(*,*)
       !write(*,'(A)') 'min_domain_list'
       !do i = 1, min_dom_list_size
       !  write(*,'(I3," |", I4)') i, min_domain_list(i)
       !end do
       !write(*,*)
       !write(*,'(A)') 'min_domain_start'
       !do i = 1, n_occ_fcidump+1
       !  write(*,'(I3," |", I4)') i, min_domain_start(i)
       !end do
       !
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       !                       Extending the domains
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! 
       write(*,'(A,F10.5)') 'rext =', rext
       write(*,'(A,F10.5)') 'fac_rcov =', fac_rcov
       write(*,'(A,I4)') 'iext =', iext
       write(*,'(A,I4)') 'max_dom_ats =', max_dom_ats
       max_n_at_dom = n_at_cas*n_occ_fcidump
       write(*,'(A,I4)') 'max_n_at_dom =', max_n_at_dom
       ! 
       !----- ADD ATOMS WHICH ARE iext BONDS AWAY 
       !      FROM AN ATOM IN THE MINIMAL DOMAIN --------------------------------
       ! 
       temp_list_size = 10*n_occ_fcidump
       call mp2alloc(temp_list, temp_list_size, ContextName, ierr, 'temp_list')
       call mp2alloc(temp_start, n_occ_fcidump+1, ContextName, ierr, 'temp_start')
       call mp2alloc(new_temp_list, temp_list_size, ContextName, ierr, 'new_temp_list')
       call mp2alloc(new_temp_start, n_occ_fcidump+1, ContextName, ierr, 'new_temp_start')
       ! 
       !iext = 1
       !
       if (iext.gt.0) then
         RADCOV(0:92)=RAYDAT(0:92)/0.52918*fac_rcov
         n_at_dom = 1
         temp_start(1) = 1
         new_temp_start(1) = 1
         do i_mo = 1, n_occ_fcidump
           do i_at = 1, n_at_cas
             if (n_at_dom.gt.temp_list_size) then
               write(*,'(A)') 'ERROR: temp_list_size too small'
               write(*,'(A,I4)') 'MO ', i_mo
               write(*,'(A,I4)') 'Atom ', i_at
               write(*,'(A,I4)') 'temp_list_size = ', temp_list_size
               write(*,'(A,F10.5)') 'iext =', iext
               write(*,*) 'Increase temp_list_size or decrease iext'
               stop
             end if
             if (domain_map(i_at,i_mo).eq.0) then
               i_at_frag = ats_cas(i_at)
               i_cell_frag = ats_cas_cells(i_at)
               x_at(:) = XA(:,i_at_frag)+XG(:,i_cell_frag)
               at = MOD(NAT(i_at_frag),100)
               rcov_at = RADCOV(at)
               do j = min_domain_start(i_mo), min_domain_start(i_mo+1)-1
                 i_at_dom = min_domain_list(j)
                 i_at_frag = ats_cas(i_at_dom)
                 i_cell_frag = ats_cas_cells(i_at_dom)
                 x_at_dom(:) = XA(:,i_at_frag)+XG(:,i_cell_frag)
                 dist = dsqrt((x_at_dom(1)-x_at(1))**2 + &
                               (x_at_dom(2)-x_at(2))**2 + &
                               (x_at_dom(3)-x_at(3))**2)
                 at_dom = MOD(NAT(i_at_frag),100)
                 rcov_at_dom = RADCOV(at_dom)
                 dist_cov = rcov_at + rcov_at_dom
                 if (dist.lt.dist_cov) then
                   domain_map(i_at,i_mo) = 1
                   temp_list(n_at_dom) = i_at
                   n_at_dom = n_at_dom + 1
                   exit
                 end if
               end do
             end if 
           end do
           temp_start(i_mo+1) = n_at_dom
         end do
       end if
       if (iext.gt.1) then
         do i = 2, iext
           n_at_dom = 1
           do i_mo = 1, n_occ_fcidump
             do i_at = 1, n_at_cas
               if (n_at_dom.gt.temp_list_size) then
                 write(*,'(A)') 'ERROR: temp_list_size too small'
                 write(*,'(A,I4)') 'MO ', i_mo
                 write(*,'(A,I4)') 'Atom ', i_at
                 write(*,'(A,I4)') 'temp_list_size = ', temp_list_size
                 write(*,'(A,F10.5)') 'iext =', iext
                 write(*,*) 'Increase temp_list_size or decrease iext'
                 stop
               end if
               if (domain_map(i_at,i_mo).eq.0) then
                 i_at_frag = ats_cas(i_at)
                 i_cell_frag = ats_cas_cells(i_at)
                 x_at(:) = XA(:,i_at_frag)+XG(:,i_cell_frag)
                 at = MOD(NAT(i_at_frag),100)
                 rcov_at = RADCOV(at)
                 do j = temp_start(i_mo), temp_start(i_mo+1)-1
                   i_at_dom = temp_list(j)
                   i_at_frag = ats_cas(i_at_dom)
                   i_cell_frag = ats_cas_cells(i_at_dom)
                   x_at_dom(:) = XA(:,i_at_frag)+XG(:,i_cell_frag)
                   dist = dsqrt((x_at_dom(1)-x_at(1))**2 + &
                                 (x_at_dom(2)-x_at(2))**2 + &
                                 (x_at_dom(3)-x_at(3))**2)
                   at_dom = MOD(NAT(i_at_frag),100)
                   rcov_at_dom = RADCOV(at_dom)
                   dist_cov = rcov_at + rcov_at_dom
                   if (dist.lt.dist_cov) then
                     domain_map(i_at,i_mo) = 1
                     new_temp_list(n_at_dom) = i_at
                     n_at_dom = n_at_dom + 1
                     exit
                   end if
                 end do
               end if
             end do
             new_temp_start(i_mo+1) = n_at_dom
           end do
           temp_start(:) = new_temp_start(:)
           temp_list(:) = new_temp_list(:)
         end do
       end if
       ! 
       call mp2dealloc(temp_list,'temp_list')
       call mp2dealloc(temp_start,'temp_start')
       call mp2dealloc(new_temp_list,'new_temp_list')
       call mp2dealloc(new_temp_start,'new_temp_start')
       ! 
       !---- ADD ATOMS WHICH ARE AT A DISTANCE < rext 
       !     TO AN ATOM IN THE MINIMAL DOMAIN -------------------------------------
       ! 
       !rext = 0.0d0
       ! 
       if (rext.gt.0.0d0) then
         do i_mo = 1, n_occ_fcidump
           do i_at = 1, n_at_cas
             if (domain_map(i_at,i_mo).eq.0) then
               i_at_frag = ats_cas(i_at)
               i_cell_frag = ats_cas_cells(i_at)
               x_at(:) = XA(:,i_at_frag)+XG(:,i_cell_frag)
               do j = min_domain_start(i_mo), min_domain_start(i_mo+1)-1
                 i_at_dom = min_domain_list(j)
                 i_at_frag = ats_cas(i_at_dom)
                 i_cell_frag = ats_cas_cells(i_at_dom)
                 x_at_dom(:) = XA(:,i_at_frag)+XG(:,i_cell_frag)
                 dist = dsqrt((x_at_dom(1)-x_at(1))**2 + &
                               (x_at_dom(2)-x_at(2))**2 + &
                               (x_at_dom(3)-x_at(3))**2)
                 if (dist.le.rext) then
                   domain_map(i_at,i_mo) = 1
                   exit
                 end if
               end do
             end if
           end do
         end do
       end if
       ! 
       !---- PRINT EXTENDED DOMAIN MAP ---------------------------------------------
       ! 
       !write(*,*)
       !write(*,'(A)') 'extended domain_map'
       !write(*,'(5X, *(I4,1X))')(j, j=1, n_occ_fcidump)
       !write(*,'(5X, *(A5))')('-----', j=1, n_occ_fcidump)
       !do i = 1, n_at_cas
       !  write(*,'(I3," |", *(I4,1X))') i, (domain_map(i,j), j=1, n_occ_fcidump)
       !end do
       ! 
       !---- CONVERT DOMAIN MAP TO domain_list ------------------------------------
       ! 
       call mp2alloc(domain_list, max_n_at_dom, ContextName, ierr, 'domain_list')
       call mp2alloc(domain_start, n_occ_fcidump+1, ContextName, ierr, 'domain_start')
       ! 
       n_at_dom = 1
       domain_start(1) = 1
       ! 
       do i_mo = 1, n_occ_fcidump
         do i_at = 1, n_at_cas
           if (n_at_dom.gt.max_n_at_dom) then
             write(*,'(A)') 'ERROR: max_n_at_dom too small'
             write(*,'(A,I4)') 'MO ', i_mo
             write(*,'(A,I4)') 'Atom ', i_at
             write(*,'(A,I4)') 'max_n_at_dom = ', max_n_at_dom
             write(*,'(A,F10.5)') 'rext =', rext
             write(*,*) 'Increase max_n_at_dom or decrease rext'
             stop
           end if
           if (domain_map(i_at,i_mo).eq.1) then
             domain_list(n_at_dom) = i_at
             n_at_dom = n_at_dom + 1
           end if
         end do
         domain_start(i_mo+1) = n_at_dom
       end do 
       n_at_dom = n_at_dom-1
       ! 
       call mp2dealloc(domain_map,'domain_map')
       ! 
       !---- PRINT domain_list AND domain_start ---------------------------------
       !
       !write(*,*)
       !write(*,'(A)') 'domain_list'
       !do i = 1, n_at_dom
       !  write(*,'(I3," |", I4)') i, domain_list(i)
       !end do
       !write(*,*)
       !write(*,'(A)') 'domain_start'
       !do i = 1, n_occ_fcidump+1
       !  write(*,'(I3," |", I4)') i, domain_start(i)
       !end do
       !
       !---- REWRITTING THE LIST IN EC's FORMAT ---------------------------------
       ! 
       !call mp2alloc(domain_end, n_occ_fcidump, ContextName, ierr, 'domain_end')
       !call izero(domain_end, n_occ_fcidump)
       ! 
       !do i_mo = 1, n_occ_fcidump
       !  domain_end(i_mo) = domain_start(i_mo+1)-1
       !  domain_start(i_mo) = domain_start(i_mo)-1 
       !  do i_at = 1, n_at_dom
       !      domain_list(i_at) = domain_list(i_at)-1
       !  end do
       !end do
   
   END SUBROUTINE Frag_domains
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   SUBROUTINE Frag_pair_list(n_occ_fcidump)
   
     USE fci_interface,  ONLY : ats_cas, ats_cas_cells
     USE basato_module,  ONLY : xa
     USE gvect_module,   ONLY : xg
   
     IMPLICIT NONE
     integer, intent(in) :: n_occ_fcidump
   
     integer :: imo, jmo, i_at_dom, i_at, j_at_dom, j_at
     integer :: i_at_frag, j_at_frag, i_cell_frag, j_cell_frag
     double precision :: R_ij, R_ij_min, R_cutoff
     double precision :: x_i_at(3), x_j_at(3)
   
     ContextName = 'Frag_pair_list'
   
     call mp2alloc(pair_list, n_occ_fcidump, n_occ_fcidump, ContextName, ierr, 'pair_list')
     call izero(pair_list, n_occ_fcidump*n_occ_fcidump)
   
     n_pairs = 0
     n_strong = 0
     n_weak = 0
     n_dist = 0
     R_cutoff = 3.0d0
   
     do imo = 1, n_occ_fcidump
       do jmo = imo, n_occ_fcidump
         n_pairs = n_pairs + 1
         pair_list(jmo,imo) = n_pairs
         R_ij_min = 100.0d0
         do i_at_dom = min_domain_start(imo), min_domain_start(imo+1)-1
           i_at = min_domain_list(i_at_dom)
           i_at_frag = ats_cas(i_at)
           i_cell_frag = ats_cas_cells(i_at)
           x_i_at(:) = XA(:,i_at_frag)+XG(:,i_cell_frag)
           do j_at_dom = min_domain_start(jmo), min_domain_start(jmo+1)-1
             j_at = min_domain_list(j_at_dom)
             j_at_frag = ats_cas(j_at)
             j_cell_frag = ats_cas_cells(j_at)
             x_j_at(:) = XA(:,j_at_frag)+XG(:,j_cell_frag)
             R_ij = dsqrt((x_i_at(1)-x_j_at(1))**2 + &
                          (x_i_at(2)-x_j_at(2))**2 + &
                          (x_i_at(3)-x_j_at(3))**2)
             if (R_ij.lt.R_ij_min) R_ij_min = R_ij
           end do
         end do
         if (R_ij_min.le.R_cutoff) then
           n_strong = n_strong + 1
         else
           n_weak = n_weak + 1
         end if
       end do
     end do   
     
     call mp2dealloc(min_domain_list,'min_domain_list')
     call mp2dealloc(min_domain_start,'min_domain_start')
   
   END SUBROUTINE Frag_pair_list
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   SUBROUTINE Frag_Fock(Fock,n_occ_fcidump,n_bas_cas,C_loc_mu_i,C_paos)
   
     USE fragment_methods, ONLY : get_eigenvalues
     USE fci_interface,    ONLY : S_proj
   
     IMPLICIT NONE
     integer, intent(in) :: n_occ_fcidump, n_bas_cas
     real*8, intent(in)  :: Fock(n_bas_cas,n_bas_cas)
     real*8, intent(in)  :: C_loc_mu_i(n_bas_cas,n_occ_fcidump)
     real*8, intent(in)  :: C_paos(n_bas_cas,n_bas_cas)
   
     integer :: imo, jmo, amo, bmo, indx 
     double precision :: s 
     real*8, pointer :: temp_mat(:,:), Focc_eigenvector(:,:), Focc_eigenvalue(:)
     real*8, pointer :: SPAO_eigenvector(:,:), SPAO_eigenvalue(:)
     real*8, pointer :: Unitary_matrix_X(:,:), SPAO_eigenvector_Y(:,:), SPAO_eigenvalue_s(:)
     real*8, pointer :: Fock_prime_virt(:,:), Fvirt_eigenvector(:,:), Fvirt_eigenvalue(:)
   
     ContextName = 'Frag_Fock'
      
     !------ Transform Fock_occ matrix in AO basis into Fock_occ in LMO basis --------------------
      
     call mp2alloc(temp_mat, n_bas_cas, n_occ_fcidump, ContextName, ierr, 'temp_mat')
     call mp2alloc(Fock_local_occ, n_occ_fcidump, n_occ_fcidump, ContextName, ierr, 'Fock_local_occ')
      
     call mxma(Fock, 1, n_bas_cas, &
               C_loc_mu_i, 1, n_bas_cas, &
               temp_mat, 1, n_bas_cas, &
               n_bas_cas, n_bas_cas, n_occ_fcidump)
      
     call mxma(C_loc_mu_i, n_bas_cas, 1, &
               temp_mat, 1, n_bas_cas, &
               Fock_local_occ, 1, n_occ_fcidump, &
               n_occ_fcidump, n_bas_cas, n_occ_fcidump)
      
     call mp2dealloc(temp_mat,'temp_mat')
     !call matprt(Fock_local_occ,n_occ_fcidump,n_occ_fcidump,'Fock_local_occ')
      
     !------- Compute Transformation matrix W ---------------------------------------------------
      
     call mp2alloc(Focc_eigenvector, n_occ_fcidump, n_occ_fcidump, ContextName, ierr, 'Focc_eigenvector')
     call mp2alloc(Focc_eigenvalue, n_occ_fcidump, ContextName, ierr, 'Focc_eigenvalue')
     call mp2alloc(Transformation_Matrix_W, n_occ_fcidump, n_occ_fcidump, ContextName, ierr, 'Transformation_Matrix_W')
   
     call get_eigenvalues(Fock_local_occ,n_occ_fcidump,Focc_eigenvector,Focc_eigenvalue)
      
     do imo = 1, n_occ_fcidump
       do jmo = 1, n_occ_fcidump
         Transformation_Matrix_W(jmo,n_occ_fcidump-imo+1)=Focc_eigenvector(jmo,imo)
       end do
     end do
   
     call mp2dealloc(Focc_eigenvector,'Focc_eigenvector')
     call mp2dealloc(Focc_eigenvalue,'Focc_eigenvalue')
     !call matprt(Transformation_Matrix_W,n_occ_fcidump,n_occ_fcidump,'Transformation_Matrix_W')
   
     !------- Transform Fock matrix into Fock_virt in PAO basis --------------------------------
   
     call mp2alloc(temp_mat, n_bas_cas, n_bas_cas, ContextName, ierr, 'temp_mat')
     call mp2alloc(Fock_local_virt, n_bas_cas, n_bas_cas, ContextName, ierr, 'Fock_local_virt')
     
     !call matprt(Fock,n_bas_cas,n_bas_cas,'Fock')
     !call matprt(C_paos,n_bas_cas,n_bas_cas,'C_paos')
     call mxma(Fock, 1, n_bas_cas, &
               C_paos, 1, n_bas_cas, &
               temp_mat, 1, n_bas_cas, &
               n_bas_cas, n_bas_cas, n_bas_cas)
      
     call mxma(C_paos, n_bas_cas, 1, &
               temp_mat, 1, n_bas_cas, &
               Fock_local_virt, 1, n_bas_cas, &
               n_bas_cas, n_bas_cas, n_bas_cas)
     !call matprt(Fock_local_virt,n_bas_cas,n_bas_cas,'Fock_local_virt')
     call mp2dealloc(temp_mat,'temp_mat')
     !call matprt(Fock_local_virt,n_bas_cas,n_bas_cas,'Fock_local_virt')
   
     !------- Build PAO overlap matrix ---------------------------------------------------------
   
     call mp2alloc(temp_mat, n_bas_cas, n_bas_cas, ContextName, ierr, 'temp_mat')
     call mp2alloc(S_paos, n_bas_cas, n_bas_cas, ContextName, ierr, 'S_paos')
      
     call mxma(S_proj, 1, n_bas_cas, &
               C_paos, 1, n_bas_cas, &
               temp_mat, 1, n_bas_cas, &
               n_bas_cas, n_bas_cas, n_bas_cas)
      
     call mxma(C_paos, n_bas_cas, 1, &
               temp_mat, 1, n_bas_cas, &
               S_paos, 1, n_bas_cas, &
               n_bas_cas, n_bas_cas, n_bas_cas)
      
     call mp2dealloc(temp_mat,'temp_mat')
     !call matprt(S_paos,n_bas_cas,n_bas_cas,'S_paos')
   
     !-------- Diagonalize PAO overlap matrix and remove small eigenvalues ---------------------
   
     call mp2alloc(SPAO_eigenvector, n_bas_cas, n_bas_cas, ContextName, ierr, 'SPAO_eigenvector')
     call mp2alloc(SPAO_eigenvalue, n_bas_cas, ContextName, ierr, 'SPAO_eigenvalue')
   
     call get_eigenvalues(S_paos,n_bas_cas,SPAO_eigenvector,SPAO_eigenvalue)
     !call matprt(SPAO_eigenvector,n_bas_cas,n_bas_cas,'SPAO_eigenvector')
     !call matprt(SPAO_eigenvalue,1,n_bas_cas,'SPAO_eigenvalue')
   
     call mp2alloc(SPAO_eigenvector_Y, n_bas_cas, n_bas_cas, ContextName, ierr, 'SPAO_eigenvector_Y')
     call mp2alloc(SPAO_eigenvalue_s, n_bas_cas, ContextName, ierr, 'SPAO_eigenvalue_s')
   
     n_virt = 0
     do amo = 1, n_bas_cas
       if (SPAO_eigenvalue(amo).gt.1.0d-6) then
         n_virt = n_virt + 1
         SPAO_eigenvalue_s(n_virt) = SPAO_eigenvalue(amo)
         SPAO_eigenvector_Y(:,n_virt) = SPAO_eigenvector(:,amo)
       end if
     end do
   
     !call matprt(SPAO_eigenvector_Y,n_bas_cas,n_virt,'SPAO_eigenvector_Y')
     !call matprt(SPAO_eigenvalue_s,1,n_virt,'SPAO_eigenvalue_s')
   
     call mp2alloc(Unitary_matrix_X, n_bas_cas, n_virt, ContextName, ierr, 'Unitary_matrix_X')
   
     do amo = 1, n_virt
       s = sqrt(SPAO_eigenvalue_s(amo))
       Unitary_matrix_X(:,amo) = SPAO_eigenvector_Y(:,amo)/s
     end do
     !call matprt(Unitary_matrix_X,n_bas_cas,n_virt,'Unitary_matrix_X')
     
     call mp2dealloc(SPAO_eigenvector,'SPAO_eigenvector')
     call mp2dealloc(SPAO_eigenvalue,'SPAO_eigenvalue')
     call mp2dealloc(SPAO_eigenvector_Y,'SPAO_eigenvector_Y')
     call mp2dealloc(SPAO_eigenvalue_s,'SPAO_eigenvalue_s')
   
     !------- Build new Fock_local_virt ---------------------------------------------------------
   
     call mp2alloc(Fock_prime_virt, n_virt, n_virt, ContextName, ierr, 'Fock_prime_virt')
     call mp2alloc(temp_mat, n_virt, n_bas_cas, ContextName, ierr, 'temp_mat')
   
     call mxma(Unitary_matrix_X, n_bas_cas, 1, &
               Fock_local_virt, 1, n_bas_cas, &
               temp_mat, 1, n_virt, &
               n_virt, n_bas_cas, n_bas_cas)
   
     call mxma(temp_mat, 1, n_virt, &
               Unitary_matrix_X, 1, n_bas_cas, &
               Fock_prime_virt, 1, n_virt, &
               n_virt, n_bas_cas, n_virt)
   
     call mp2dealloc(temp_mat,'temp_mat')
     !call matprt(Fock_prime_virt,n_virt,n_virt,'Fock_prime_virt')
   
     !------- Diagonalize Fock_prime_virt ------------------------------------------------------
   
     call mp2alloc(Fvirt_eigenvector, n_virt, n_virt, ContextName, ierr, 'Fvirt_eigenvector')
     call mp2alloc(Fvirt_eigenvalue, n_virt, ContextName, ierr, 'Fvirt_eigenvalue')
   
     call get_eigenvalues(Fock_prime_virt,n_virt,Fvirt_eigenvector,Fvirt_eigenvalue)
     
     !------- Compute Transformation matrix Q --------------------------------------------------
     
     call mp2alloc(temp_mat, n_virt, n_bas_cas, ContextName, ierr, 'temp_mat')
     call mp2alloc(Transformation_Matrix_Q, n_bas_cas, n_virt, ContextName, ierr, 'Transformation_Matrix_Q')
   
     do amo = 1, n_virt
       do bmo = 1, n_virt
         indx = n_virt-amo+1
         temp_mat(bmo,amo) = Fvirt_eigenvector(bmo,indx)
       end do
     end do
   
     call mxma(Unitary_matrix_X, 1, n_bas_cas, &
               temp_mat, 1, n_virt, &
               Transformation_Matrix_Q, 1, n_bas_cas, &
               n_bas_cas, n_virt, n_virt)
     
     call mp2dealloc(temp_mat,'temp_mat')
     call mp2dealloc(Fock_prime_virt,'Fock_prime_virt')
     call mp2dealloc(Fvirt_eigenvector,'Fvirt_eigenvector')
     call mp2dealloc(Fvirt_eigenvalue,'Fvirt_eigenvalue')
     !call matprt(Transformation_Matrix_Q,n_bas_cas,n_virt,'Transformation_Matrix_Q')
   
     !------- Compute Transformation matrix V --------------------------------------------------
   
     call mp2alloc(Transformation_Matrix_V, n_bas_cas, n_bas_cas, ContextName, ierr, 'Transformation_Matrix_V')
   
     call mxma(Unitary_matrix_X, 1, n_bas_cas, &
               Unitary_matrix_X, n_bas_cas, 1, &
               Transformation_Matrix_V, 1, n_bas_cas, &
               n_bas_cas, n_virt, n_bas_cas)
               
     call mp2dealloc(Unitary_matrix_X,'Unitary_matrix_X')
     !call matprt(Transformation_Matrix_V,n_bas_cas,n_bas_cas,'Transformation_Matrix_V')
   
   END SUBROUTINE Frag_Fock  
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   SUBROUTINE Frag_integrals(C_loc_mu_i,C_paos,n_bas_cas,n_occ_fcidump, &
               n_all_aux,n_all_aux_full,n_pi_aux_at,n_pi_aux_at_red,needed_P,naf,max_pi_aux_atoms,Iunit_3ints_for_litf, &
               n_orb_cas,occ_shift,C_p_mn,C_p_mn_ao,Iunit_2ints_for_litf)
   
     USE INPOUT,      ONLY: leggi, apri,riavvolgi
   
     IMPLICIT NONE
     integer, intent(in) :: n_occ_fcidump, n_bas_cas, n_orb_cas, occ_shift
     integer, intent(in) :: n_all_aux, n_all_aux_full, naf, max_pi_aux_atoms, Iunit_3ints_for_litf, Iunit_2ints_for_litf
     integer, intent(in), dimension(naf) :: n_pi_aux_at,n_pi_aux_at_red
     real*8, intent(in), dimension(n_bas_cas,n_occ_fcidump)  :: C_loc_mu_i
     real*8, intent(in), dimension(n_bas_cas,n_bas_cas) :: C_paos
     double precision, intent(in), dimension(*)::C_p_mn
     double precision, intent(in), dimension(*)::C_p_mn_ao
     integer, intent(in), dimension(max_pi_aux_atoms,naf) :: needed_P
   
     integer :: P, Q, ind_P, ind_Q, iatom_P, iatom_Q, Q_on_at, P_on_at
     real*8, pointer :: Iunit_2buffer(:)
     
     ContextName = 'Frag_integrals'
   
     call mp2alloc(Int_p_ai,n_all_aux,n_bas_cas,n_occ_fcidump,ContextName,ierr,'Int_p_ai')
     call mp2alloc(Int_p_ji,n_all_aux,n_occ_fcidump,n_occ_fcidump,ContextName,ierr,'Int_p_ji')
   
     call transform_occ_litf(C_loc_mu_i,C_paos,  &
                   n_pi_aux_at,max_pi_aux_atoms,n_pi_aux_at_red,needed_P,Iunit_3ints_for_litf,  &
                   n_occ_fcidump,n_bas_cas,naf,                   &
                   n_orb_cas,n_bas_cas,occ_shift,C_p_mn,C_p_mn_ao,     &
                   Int_p_ai,Int_p_ji,n_all_aux)
   
     call mp2alloc(Int_p_ab,n_all_aux,n_bas_cas,n_bas_cas,ContextName,ierr,'Int_p_ab')
   
     call transform_virt_litf(C_paos,  &
                   n_pi_aux_at,max_pi_aux_atoms,n_pi_aux_at_red,needed_P,Iunit_3ints_for_litf,           &
                   n_bas_cas,naf,                                         &
                   n_orb_cas,n_bas_cas,occ_shift,C_p_mn,C_p_mn_ao,             &
                   Int_p_ab,n_all_aux)
   
     call mp2alloc(Coulomb_Metrik_J_PQ,n_all_aux,n_all_aux,ContextName,ierr,'Coulomb_Metrik_J_PQ')
     call mp2alloc(Iunit_2buffer,n_all_aux_full**2,ContextName,ierr,'Iunit_2buffer')
     
     call riavvolgi(Iunit_2ints_for_litf)
     call leggi(Iunit_2ints_for_litf,Iunit_2buffer,n_all_aux_full**2)
   
     ind_Q=0
     Q=0
     do iatom_Q=1,naf
       do Q_on_at=1,n_pi_aux_at(iatom_Q)
          Q=Q+1
          if (Needed_P(Q_on_at,iatom_Q).eq.0) cycle
          ind_Q=ind_Q+1
          ind_P=0
          P=0
          do iatom_P=1,naf
             do P_on_at=1,n_pi_aux_at(iatom_P)
                 P=P+1
                 if (Needed_P(P_on_at,iatom_P).eq.0) cycle
                 ind_P=ind_P+1
                 Coulomb_Metrik_J_PQ(ind_P,ind_Q) = Iunit_2buffer(P+(Q-1)*n_all_aux_full)
             end do
          end do
       end do
     end do
     
     call mp2dealloc(Iunit_2buffer,'Iunit_2buffer')
      
     do P=1,n_all_aux
       do Q=1,n_all_aux
        if((Coulomb_Metrik_J_PQ(P,Q)-Coulomb_Metrik_J_PQ(Q,P)).gt.1.d-10) then
         write(*,*)'no sym',P,Q,Coulomb_Metrik_J_PQ(P,Q),Coulomb_Metrik_J_PQ(Q,P)
         stop
        end if
       end do
      end do
      
   END SUBROUTINE Frag_integrals
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   
   SUBROUTINE write_frag_molpro_input(n_bas_cas,n_occ_fcidump,n_all_aux)  
   
     USE fci_interface,          ONLY: n_at_cas, cas_mul_map
   
     IMPLICIT NONE
     
     integer, intent(in) :: n_bas_cas		! number of atomic orbitals
     integer, intent(in) :: n_occ_fcidump		! number of molecular orbitals
     integer, intent(in) :: n_all_aux		! number of auxiliary basis functions
       
     INTEGER :: infostatus,imo,jmo,amo,bmo,P,Q,stream_input,stream_molpro,iAtom,iA,iP,atom_i,iA_Atom
     INTEGER :: atom_a,atom_b,iB_Atom,iDomain
     INTEGER :: nPao_on_atom
     integer :: n_all_aux_at
   
     ContextName = 'write_frag_molpro_input'
   
     stream_input=74
     stream_molpro=75
   
     open(unit=stream_input,file="fragment_itf_input.dat",status='new',action='write',iostat=infostatus)
   
     if(infostatus.ne.0) then
      write(*,*)'"fragment_itf_input.dat"-File can not be created!'
      stop
     end if
   
     write(stream_input,'(A)')'Input for CIS calculation in Molpro ITF.'
     write(stream_input,7399) n_bas_cas		    	
     write(stream_input,7400) n_occ_fcidump	        
     write(stream_input,7401) 0 !nAct_in_fragment   !Number of active orbitals in fragment
     write(stream_input,7402) 0 !nInt_in_fragment   !Number of internal orbitals in fragment
     write(stream_input,7403) n_virt		      	
     write(stream_input,7404) n_at_cas		    
     write(stream_input,7405) n_occ_fcidump*2	
     write(stream_input,7406) n_all_aux		     
     write(stream_input,7407) n_pairs
     write(stream_input,7408) n_strong
     write(stream_input,7409) n_weak
     write(stream_input,7410) n_dist
     write(stream_input,7411) n_at_dom
   
   7399     format('naos    =',i)
   7400     format('noccs   =',i)
   7401     format('nactv   =',i)
   7402     format('nintl   =',i)
   7403     format('nvirt   =',i)
   7404     format('natoms  =',i)
   7405     format('nelec   =',i)
   7406     format('nPfun   =',i)
   7407     format('npairs  =',i)
   7408     format('nstrong =',i)
   7409     format('nweak   =',i)
   7410     format('ndist   =',i)
   7411     format('ndoma   =',i)
   
     write(stream_input,'(A)')'PAOs'             
     do iAtom=1,n_at_cas         
       nPao_on_atom = cas_mul_map(iAtom+1)-cas_mul_map(iAtom)
       write(stream_input,'(i)')nPao_on_atom
     end do
   
     write(stream_input,'(A)')'Pfun'                  
     do iAtom=1,n_at_cas-1
       n_all_aux_at = n_all_aux/n_at_cas   
       write(stream_input,'(i)')n_all_aux_at
     end do
     n_all_aux_at = n_all_aux/n_at_cas + MOD(n_all_aux, n_at_cas)
     write(stream_input,'(i)')n_all_aux_at
   
     write(stream_input,'(A)')'listp_ij'        
     do imo=1,n_occ_fcidump
      do jmo=1,n_occ_fcidump
       write(stream_input,'(i)')pair_list(imo,jmo)
      end do
     end do
   
     write(stream_input,'(A)')'domain_start'    
     do imo=1,n_occ_fcidump				
      write(stream_input,'(i)')domain_start(imo)-1		
     end do
   
     write(stream_input,'(A)')'domain_end'     
     do imo=1,n_occ_fcidump				
      write(stream_input,'(i)')domain_start(imo+1)-1
     end do
   
     write(stream_input,'(A)')'domain_list'     
     do iDomain=1,n_at_dom				
      write(stream_input,'(i)')domain_list(iDomain)-1	
     end do
   
     call mp2dealloc(pair_list,'pair_list')        
     call mp2dealloc(domain_start,'domain_start')	
     !call mp2dealloc(domain_end,'domain_end')
     call mp2dealloc(domain_list,'domain_list')	
   
     write(stream_input,'(A)')'Fock_local_occ'   
     do imo = 1, n_occ_fcidump
      do jmo = 1, n_occ_fcidump
       write(stream_input,'(f)')Fock_local_occ(imo,jmo)
      end do
     end do
   
     call mp2dealloc(Fock_local_occ,'Fock_local_occ')
   
     write(stream_input,'(A)')'Fock_local_virt' 
     do amo = 1, n_bas_cas
      do bmo = 1, n_bas_cas
       write(stream_input,'(f)')Fock_local_virt(amo,bmo)
      end do
     end do
   
     call mp2dealloc(Fock_local_virt,'Fock_local_virt')
   
     write(stream_input,'(A)')'Overlap_local_virt'         
     do amo = 1, n_bas_cas
      do bmo = 1, n_bas_cas
       write(stream_input,'(f)')S_paos(amo,bmo)
      end do
     end do
   
     call mp2dealloc(S_paos,'S_paos')
   
     write(stream_input,'(A)')'Transformation_Matrix_Q'  
     do amo = 1, n_bas_cas
      do bmo = 1, n_virt
       write(stream_input,'(f)')Transformation_Matrix_Q(amo,bmo)
      end do
     end do
   
     call mp2dealloc(Transformation_Matrix_Q,'Transformation_Matrix_Q') 
   
     write(stream_input,'(A)')'Transformation_Matrix_W'
     do imo = 1, n_occ_fcidump
      do jmo = 1, n_occ_fcidump
       write(stream_input,'(f)')Transformation_Matrix_W(jmo,imo)
      end do
     end do
   
     call mp2dealloc(Transformation_Matrix_W,'Transformation_Matrix_W')
   
     write(stream_input,'(A)')'Transformation_Matrix_V'
     do amo = 1, n_bas_cas
      do bmo = 1, n_bas_cas
       write(stream_input,'(f)')Transformation_Matrix_V(bmo,amo)
      end do
     end do
   
     call mp2dealloc(Transformation_Matrix_V,'Transformation_Matrix_V')
   
     write(stream_input,'(A)')'iaP integral'                
     do imo = 1, n_occ_fcidump
       do amo=1, n_bas_cas
         do iP = 1, n_all_aux
          write(stream_input,'(f)')Int_p_ai(iP,amo,imo)
         end do
       end do
     end do
   
     call mp2dealloc(Int_p_ai,'Int_p_ai')
   
     write(stream_input,'(A)')'Coulomb_Metrik_J_PQ'
     do P = 1, n_all_aux
      do Q = 1, n_all_aux
       write(stream_input,'(f)')Coulomb_Metrik_J_PQ(P,Q)
      end do
     end do
   
     call mp2dealloc(Coulomb_Metrik_J_PQ,'Coulomb_Metrik_J_PQ')
   
     write(stream_input,'(A)')'ijP integral'
     do imo = 1, n_occ_fcidump
      do jmo = 1, n_occ_fcidump
        do iP = 1, n_all_aux
         write(stream_input,'(f)')Int_p_ji(iP,jmo,imo)
        end do
      end do
     end do
   
     call mp2dealloc(Int_p_ji,'Int_p_ji')
   
     write(stream_input,'(A)')'abP integral'
     do amo = 1, n_bas_cas
       do bmo = 1, n_bas_cas
         do iP = 1, n_all_aux
           write(stream_input,'(f)')Int_p_ab(iP,bmo,amo)
         end do
       end do
     end do
   
     call mp2dealloc(Int_p_ab,'Int_p_ab')
   
     close(unit=stream_input)
   
    END SUBROUTINE write_frag_molpro_input
   
   
   
   END MODULE Loc_aper
   