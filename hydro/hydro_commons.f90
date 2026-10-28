module hydro_commons
  use amr_parameters
  use hydro_parameters

  type hydro_kernel_t
     integer::iu1,iu2,ju1,ju2,ku1,ku2
     integer::if1,if2,jf1,jf2,kf1,kf2
     integer::io1,io2,jo1,jo2,ko1,ko2
     real(kind=8),dimension(:,:,:,:),allocatable::uloc
     real(kind=8),dimension(:,:,:,:),allocatable::uvirtual
     ! Left and right virtual-cell face states in each active direction
     real(kind=8),dimension(:,:,:,:,:,:),allocatable::uface
     ! Pointwise fluxes at virtual-cell interfaces
     real(kind=8),dimension(:,:,:,:,:),allocatable::fvirtual
     real(kind=8),dimension(:,:,:,:),allocatable::gloc
     real(kind=8),dimension(:,:,:,:),allocatable::qloc
     real(kind=8),dimension(:,:,:),allocatable::cloc
     real(kind=8),dimension(:,:,:,:,:),allocatable::flux
     real(kind=8),dimension(:,:,:,:,:),allocatable::tmp
     real(kind=8),dimension(:,:,:,:,:),allocatable::dq
     real(kind=8),dimension(:,:,:,:,:),allocatable::qm
     real(kind=8),dimension(:,:,:,:,:),allocatable::qp
     real(kind=8),dimension(:,:,:,:),allocatable::fx
     real(kind=8),dimension(:,:,:,:),allocatable::tx
     real(kind=8),dimension(:,:,:),allocatable::divu
     logical,dimension(:,:,:),allocatable::okloc
     integer,dimension(:,:,:),allocatable::cellloc
     integer,dimension(:,:,:),allocatable::childloc
     integer,dimension(:,:,:),allocatable::gridloc
     integer,dimension(:,:,:,:),allocatable::nborloc
#ifdef MHD
     real(kind=8),dimension(:,:,:,:),allocatable::bloc
     real(kind=8),dimension(:,:,:),allocatable::emfx
     real(kind=8),dimension(:,:,:),allocatable::emfy
     real(kind=8),dimension(:,:,:),allocatable::emfz
     real(kind=8),dimension(:,:,:),allocatable::Ex
     real(kind=8),dimension(:,:,:),allocatable::Ey
     real(kind=8),dimension(:,:,:),allocatable::Ez
     real(kind=8),dimension(:,:,:,:),allocatable::bf
     real(kind=8),dimension(:,:,:,:,:),allocatable::dbf
     real(kind=8),dimension(:,:,:,:,:),allocatable::qRT
     real(kind=8),dimension(:,:,:,:,:),allocatable::qRB
     real(kind=8),dimension(:,:,:,:,:),allocatable::qLT
     real(kind=8),dimension(:,:,:,:,:),allocatable::qLB
     integer,dimension(:,:,:,:),allocatable::nborsonloc
#endif
   contains
     procedure :: init => init_hydro_kernel
     procedure :: size => size_hydro_kernel
     procedure :: unpack_stencil => unpack_hydro_stencil
     procedure :: interpolate_faces => interpolate_hydro_faces
  end type hydro_kernel_t

  type hydro_workspace_t
     type(hydro_kernel_t)::kernel_1,kernel_2,kernel_4,kernel_8,kernel_16,kernel_32
  end type hydro_workspace_t

contains

  subroutine init_hydro_kernel(h,nn)
    use amr_parameters, only: ndim
    use hydro_parameters, only: nvar,nvar_per_minicell,dofs_per_cell
    use rt_parameters, only: nrtvar

    integer::nn,ndofsy,ndofsz
    class(hydro_kernel_t)::h

    h%iu1=-1; h%iu2=nn+2
    h%ju1=(1-ndim/2)-1*(ndim/2); h%ju2=(1-ndim/2)+(nn+2)*(ndim/2)
    h%ku1=(1-ndim/3)-1*(ndim/3); h%ku2=(1-ndim/3)+(nn+2)*(ndim/3)
    h%if1=1; h%if2=nn+1
    h%jf1=1; h%jf2=(1-ndim/2)+(nn+1)*(ndim/2)
    h%kf1=1; h%kf2=(1-ndim/3)+(nn+1)*(ndim/3)
    h%io1=0; h%io2=nn/2+1
    h%jo1=(1-ndim/2); h%jo2=(1-ndim/2)+(nn/2+1)*(ndim/2)
    h%ko1=(1-ndim/3); h%ko2=(1-ndim/3)+(nn/2+1)*(ndim/3)

    allocate(h%uloc (h%iu1:h%iu2,h%ju1:h%ju2,h%ku1:h%ku2,1:nvar))
    if(dofs_per_cell>1)then
       ndofsy=merge(dofs_per_cell,1,ndim>1)
       ndofsz=merge(dofs_per_cell,1,ndim>2)
       allocate(h%uvirtual(dofs_per_cell*h%iu1:dofs_per_cell*(h%iu2+1)-1, &
            & ndofsy*h%ju1:ndofsy*(h%ju2+1)-1, &
            & ndofsz*h%ku1:ndofsz*(h%ku2+1)-1,1:nvar_per_minicell))
       allocate(h%uface(lbound(h%uvirtual,1):ubound(h%uvirtual,1), &
            & lbound(h%uvirtual,2):ubound(h%uvirtual,2), &
            & lbound(h%uvirtual,3):ubound(h%uvirtual,3), &
            & 1:nvar_per_minicell,1:2,1:ndim))
       allocate(h%fvirtual(lbound(h%uvirtual,1):ubound(h%uvirtual,1), &
            & lbound(h%uvirtual,2):ubound(h%uvirtual,2), &
            & lbound(h%uvirtual,3):ubound(h%uvirtual,3), &
            & 1:nvar_per_minicell,1:ndim))
    endif
    allocate(h%gloc (h%iu1:h%iu2,h%ju1:h%ju2,h%ku1:h%ku2,1:ndim))
    allocate(h%qloc (h%iu1:h%iu2,h%ju1:h%ju2,h%ku1:h%ku2,1:nprim))
    allocate(h%cloc (h%iu1:h%iu2,h%ju1:h%ju2,h%ku1:h%ku2))
    allocate(h%okloc(h%iu1:h%iu2,h%ju1:h%ju2,h%ku1:h%ku2))
    allocate(h%dq   (h%iu1:h%iu2,h%ju1:h%ju2,h%ku1:h%ku2,1:nprim,1:ndim))
    allocate(h%qm   (h%iu1:h%iu2,h%ju1:h%ju2,h%ku1:h%ku2,1:nprim,1:ndim))
    allocate(h%qp   (h%iu1:h%iu2,h%ju1:h%ju2,h%ku1:h%ku2,1:nprim,1:ndim))
    allocate(h%fx   (h%iu1:h%iu2,h%ju1:h%ju2,h%ku1:h%ku2,1:nprim))
    allocate(h%tx   (h%iu1:h%iu2,h%ju1:h%ju2,h%ku1:h%ku2,1:2    ))

    allocate(h%flux(h%if1:h%if2,h%jf1:h%jf2,h%kf1:h%kf2,1:nprim,1:ndim))
    allocate(h%tmp (h%if1:h%if2,h%jf1:h%jf2,h%kf1:h%kf2,1:2    ,1:ndim))
    allocate(h%divu(h%if1:h%if2,h%jf1:h%jf2,h%kf1:h%kf2))

    allocate(h%childloc(h%io1:h%io2,h%jo1:h%jo2,h%ko1:h%ko2))
    allocate(h%gridloc (h%io1:h%io2,h%jo1:h%jo2,h%ko1:h%ko2))
    allocate(h%cellloc (h%io1:h%io2,h%jo1:h%jo2,h%ko1:h%ko2))
    allocate(h%nborloc (h%io1:h%io2,h%jo1:h%jo2,h%ko1:h%ko2,1:twondim))
#ifdef MHD
    allocate(h%bloc(h%iu1:h%iu2,h%ju1:h%ju2,h%ku1:h%ku2,1:6))
    allocate(h%emfx(h%if1:h%if2,h%jf1:h%jf2,h%kf1:h%kf2))
    allocate(h%emfy(h%if1:h%if2,h%jf1:h%jf2,h%kf1:h%kf2))
    allocate(h%emfz(h%if1:h%if2,h%jf1:h%jf2,h%kf1:h%kf2))
    allocate(h%Ex(h%iu1:h%iu2,h%ju1:h%ju2,h%ku1:h%ku2))
    allocate(h%Ey(h%iu1:h%iu2,h%ju1:h%ju2,h%ku1:h%ku2))
    allocate(h%Ez(h%iu1:h%iu2,h%ju1:h%ju2,h%ku1:h%ku2))
    allocate(h%bf(h%iu1:h%iu2+1,h%ju1:h%ju2+1,h%ku1:h%ku2+1,1:3))
    allocate(h%dbf(h%iu1:h%iu2+1,h%ju1:h%ju2+1,h%ku1:h%ku2+1,1:3,1:2))
    allocate(h%qRT(h%iu1:h%iu2,h%ju1:h%ju2,h%ku1:h%ku2,1:nprim,1:3))
    allocate(h%qRB(h%iu1:h%iu2,h%ju1:h%ju2,h%ku1:h%ku2,1:nprim,1:3))
    allocate(h%qLT(h%iu1:h%iu2,h%ju1:h%ju2,h%ku1:h%ku2,1:nprim,1:3))
    allocate(h%qLB(h%iu1:h%iu2,h%ju1:h%ju2,h%ku1:h%ku2,1:nprim,1:3))
    allocate(h%nborsonloc(h%io1:h%io2,h%jo1:h%jo2,h%ko1:h%ko2,1:twondim))
#endif
  end subroutine init_hydro_kernel

  subroutine unpack_hydro_stencil(h)
    use amr_parameters, only: ndim,dp
    use hydro_parameters, only: nvar,nvar_per_minicell,dofs_per_cell,unpack_parent_cell

    class(hydro_kernel_t),intent(inout)::h
    integer::i,j,k,ii,jj,kk,iminicell,ndofsy,ndofsz
    real(dp)::parent(1:nvar)
    real(dp)::minicells(1:nvar_per_minicell,1:dofs_per_cell**ndim)

    ndofsy=merge(dofs_per_cell,1,ndim>1)
    ndofsz=merge(dofs_per_cell,1,ndim>2)
    do k=h%ku1,h%ku2
       do j=h%ju1,h%ju2
          do i=h%iu1,h%iu2
             parent=h%uloc(i,j,k,:)
             call unpack_parent_cell(parent,minicells)
             do kk=0,ndofsz-1
                do jj=0,ndofsy-1
                   do ii=0,dofs_per_cell-1
                      iminicell=1+ii+dofs_per_cell*jj+dofs_per_cell**2*kk
                      h%uvirtual(dofs_per_cell*i+ii,ndofsy*j+jj,ndofsz*k+kk,:)=minicells(:,iminicell)
                   end do
                end do
             end do
          end do
       end do
    end do
  end subroutine unpack_hydro_stencil

  subroutine interpolate_hydro_faces(h,idim)
    use amr_parameters, only: ndim,dp
    use hydro_parameters, only: p3_center,p3_face_left,p3_face_right

    class(hydro_kernel_t),intent(inout)::h
    integer,intent(in)::idim
    integer::i,j,k,ii,jj,kk,jdim
    integer::reach(1:3)
    real(dp)::wl(-2:2,1:3),wr(-2:2,1:3),wleft,wright

    wl=0.0_dp
    wr=0.0_dp
    reach=0
    do jdim=1,ndim
       if(jdim==idim)then
          wl(:,jdim)=p3_face_left
          wr(:,jdim)=p3_face_right
          reach(jdim)=2
       else
          wl(-1:1,jdim)=p3_center
          wr(-1:1,jdim)=p3_center
          reach(jdim)=1
       endif
    end do
    do jdim=ndim+1,3
       wl(0,jdim)=1.0_dp
       wr(0,jdim)=1.0_dp
    end do

    do k=lbound(h%uvirtual,3)+reach(3),ubound(h%uvirtual,3)-reach(3)
       do j=lbound(h%uvirtual,2)+reach(2),ubound(h%uvirtual,2)-reach(2)
          do i=lbound(h%uvirtual,1)+reach(1),ubound(h%uvirtual,1)-reach(1)
             h%uface(i,j,k,:,1,idim)=0.0_dp
             h%uface(i,j,k,:,2,idim)=0.0_dp
             do kk=-reach(3),reach(3)
                do jj=-reach(2),reach(2)
                   do ii=-reach(1),reach(1)
                      wleft=wl(ii,1)*wl(jj,2)*wl(kk,3)
                      wright=wr(ii,1)*wr(jj,2)*wr(kk,3)
                      h%uface(i,j,k,:,1,idim)=h%uface(i,j,k,:,1,idim) &
                           & +wleft*h%uvirtual(i+ii,j+jj,k+kk,:)
                      h%uface(i,j,k,:,2,idim)=h%uface(i,j,k,:,2,idim) &
                           & +wright*h%uvirtual(i+ii,j+jj,k+kk,:)
                   end do
                end do
             end do
          end do
       end do
    end do
  end subroutine interpolate_hydro_faces

  function size_hydro_kernel(h)
    use amr_parameters, only: ndim
    use hydro_parameters, only: nvar
    integer::size_hydro_kernel
    class(hydro_kernel_t)::h

    integer::nint

    nint=0
    nint=nint+size(transfer(h%uloc,(/1/)))
    if(allocated(h%uvirtual))nint=nint+size(transfer(h%uvirtual,(/1/)))
    if(allocated(h%uface))nint=nint+size(transfer(h%uface,(/1/)))
    if(allocated(h%fvirtual))nint=nint+size(transfer(h%fvirtual,(/1/)))
    nint=nint+size(transfer(h%gloc,(/1/)))
    nint=nint+size(transfer(h%qloc,(/1/)))
    nint=nint+size(transfer(h%cloc,(/1/)))
    nint=nint+size(transfer(h%okloc,(/1/)))
    nint=nint+size(transfer(h%dq   ,(/1/)))
    nint=nint+size(transfer(h%qm   ,(/1/)))
    nint=nint+size(transfer(h%qp   ,(/1/)))
    nint=nint+size(transfer(h%fx   ,(/1/)))
    nint=nint+size(transfer(h%tx   ,(/1/)))

    nint=nint+size(transfer(h%flux ,(/1/)))
    nint=nint+size(transfer(h%tmp  ,(/1/)))
    nint=nint+size(transfer(h%divu ,(/1/)))

    nint=nint+size(transfer(h%childloc,(/1/)))
    nint=nint+size(transfer(h%gridloc ,(/1/)))
    nint=nint+size(transfer(h%cellloc ,(/1/)))
    nint=nint+size(transfer(h%nborloc ,(/1/)))
#ifdef MHD
    nint=nint+size(transfer(h%bloc,(/1/)))
    nint=nint+size(transfer(h%emfx,(/1/)))
    nint=nint+size(transfer(h%emfy,(/1/)))
    nint=nint+size(transfer(h%emfz,(/1/)))
    nint=nint+size(transfer(h%Ex  ,(/1/)))
    nint=nint+size(transfer(h%Ey  ,(/1/)))
    nint=nint+size(transfer(h%Ez  ,(/1/)))
    nint=nint+size(transfer(h%bf  ,(/1/)))
    nint=nint+size(transfer(h%dbf ,(/1/)))
    nint=nint+size(transfer(h%qRT ,(/1/)))
    nint=nint+size(transfer(h%qRB ,(/1/)))
    nint=nint+size(transfer(h%qLT ,(/1/)))
    nint=nint+size(transfer(h%qLB ,(/1/)))
    nint=nint+size(transfer(h%nborsonloc,(/1/)))
#endif
    size_hydro_kernel = nint

  end function size_hydro_kernel

end module hydro_commons
