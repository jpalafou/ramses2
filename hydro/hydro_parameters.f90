module hydro_parameters
  use amr_parameters, only: ndim
#if DOFS_PER_CELL>1
  use amr_parameters, only: dp
#endif

  ! Number of independant variables
#ifndef NENER
  integer,parameter::nener=0
#else
  integer,parameter::nener=NENER
#endif
  ! Number of hydro degrees of freedom per coordinate direction
#ifndef DOFS_PER_CELL
  integer,parameter::dofs_per_cell=1
#else
  integer,parameter::dofs_per_cell=DOFS_PER_CELL
#endif
  ! Number of packed variables and physical variables in one state
#ifndef NVAR
  integer,parameter::nvar=(5+nener)*dofs_per_cell**ndim
#else
  integer,parameter::nvar=NVAR
#endif
#if DOFS_PER_CELL>1
  integer,parameter::nminicell=dofs_per_cell**ndim
  integer,parameter::nvar_per_minicell=nvar/nminicell
#endif

#ifdef MHD
  integer,parameter::nprim=NVAR+3
  integer,parameter::ie=8
#else
  integer,parameter::nprim=NVAR
  integer,parameter::ie=5
#endif

#ifdef NION
  integer,parameter::nion=NION  ! # of ionization fractions species
#else
  integer,parameter::nion=1
#endif

  integer,parameter::solver_llf=1
  integer,parameter::solver_hll=2
  integer,parameter::solver_hllc=3
  integer,parameter::solver_hlld=4
  integer,parameter::solver_roe=5
  integer,parameter::solver_upwind=6
#ifdef _METAL
  integer,parameter::solver_uct_hlld=7
#endif

  integer,parameter::solver2d_llf=1
  integer,parameter::solver2d_hllf=2
  integer,parameter::solver2d_hlla=3
  integer,parameter::solver2d_hlld=4
  integer,parameter::solver2d_roe=5
  integer,parameter::solver2d_upwind=6

#if DOFS_PER_CELL>1
  ! Unlimited fourth-order conservative interpolation stencils
  real(dp),parameter::p3_center(-1:1)=(/-1,26,-1/)/24.0_dp
  real(dp),parameter::p3_face_left(-2:2)=(/-1,10,20,-6,1/)/24.0_dp
  real(dp),parameter::p3_face_right(-2:2)=(/1,-6,20,10,-1/)/24.0_dp
  real(dp),parameter::p3_integrate(-1:1)=(/1,22,1/)/24.0_dp

contains

  pure subroutine unpack_parent_cell(parent,minicells)
    real(dp),intent(in)::parent(1:nvar)
    real(dp),intent(out)::minicells(1:nvar_per_minicell,1:nminicell)
    integer::iminicell,ilo,ihi

    minicells(:,1)=parent(1:nvar_per_minicell)
    do iminicell=2,nminicell
       ilo=(iminicell-1)*nvar_per_minicell+1
       ihi=iminicell*nvar_per_minicell
       minicells(:,iminicell)=parent(1:nvar_per_minicell)+parent(ilo:ihi)
       minicells(:,1)=minicells(:,1)-parent(ilo:ihi)
    end do
  end subroutine unpack_parent_cell

  pure subroutine pack_parent_cell(minicells,parent)
    real(dp),intent(in)::minicells(1:nvar_per_minicell,1:nminicell)
    real(dp),intent(out)::parent(1:nvar)
    integer::iminicell,ilo,ihi

    parent(1:nvar_per_minicell)=sum(minicells,dim=2)/real(nminicell,kind=dp)
    do iminicell=2,nminicell
       ilo=(iminicell-1)*nvar_per_minicell+1
       ihi=iminicell*nvar_per_minicell
       parent(ilo:ihi)=minicells(:,iminicell)-parent(1:nvar_per_minicell)
    end do
  end subroutine pack_parent_cell
#endif

end module hydro_parameters

module const

  ! Some useful constant
  real(kind=8),parameter ::bigreal = 1.0d+30
  real(kind=8),parameter ::zero = 0.0
  real(kind=8),parameter ::one = 1.0
  real(kind=8),parameter ::two = 2.0
  real(kind=8),parameter ::three = 3.0
  real(kind=8),parameter ::four = 4.0
  real(kind=8),parameter ::two3rd = 0.6666666666666667
  real(kind=8),parameter ::half = 0.5
  real(kind=8),parameter ::third = 0.33333333333333333
  real(kind=8),parameter ::forth = 0.25
  real(kind=8),parameter ::sixth = 0.16666666666666667

end module const
