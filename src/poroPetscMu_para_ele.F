!     ------------------------------------------------------------------
!
!     porosplit.f
!
!
!     Created on Wed May 18 15:45:04 2011
!     Copyright (c) 2011 MyCompany. All rights reserved.
!
!
!     ------------------------------------------------------------------

      program phase
      implicit none

#include <petsc/finclude/petsc.h>

      integer maxnd,maxele,maxfc,maxwk,maxsol,maxedg,maxdbc
      integer maxcon,ndim,nnds,nints,ndpn,maxbdnd,maxprocs,csteps
      parameter (maxnd=100000,maxele=15000,maxfc=10000,
     &           maxcon=4,maxdbc=600000,maxwk=100,maxsol=40,maxedg=20,
     &           ndim=2,nnds=9,nints=9,ndpn=6,maxbdnd=100,
     &           maxprocs=1000,csteps=1000)
      real*8 node(maxnd,ndim),B(maxele,nnds,nints,ndim)
      real*8 det(maxele,nints),sumsquT,crackloc
      real*8 u(maxnd,ndim),mu(maxnd),NN(nnds,nints)
      real*8 gradv(maxele,nints,ndim,ndim,2),tau(maxele,nints,ndim,ndim)
      real*8 strn(maxele,nints,ndim,ndim,2),vint(maxele,nints,ndim)
      real*8 gradmu(maxele,nints,ndim),strs(maxele,nints,ndim,ndim)
      real*8 v(maxnd,ndim),psimu(maxele,nints),p(maxnd)
      real*8 gradp(maxele,nints,ndim),sumsqf,sumsqu
      real*8 ksi(maxele,nints,ndim),muint(maxele,nints,2)
      real*8 conAu(maxdbc,maxcon),conBu(maxdbc,maxcon)
      real*8 A(maxwk),dbcval(maxdbc),fbcval(maxfc)
      real*8 drhs(maxnd*ndpn),f(maxnd*ndpn)
      real*8 mu0,t,dt,beta,x,y,xold,fint(maxnd*ndpn)
      real*8 unitdbc(maxdbc),unitfbc(maxfc),e0
      real*8 mumin,mumax,xmin,xmax,xcrack,papp,uapp,vapp,muapp
      real*8 conAmu(maxdbc,maxcon),conBmu(maxdbc,maxcon)
      real*8 conAvp(maxdbc,maxcon),conBvp(maxdbc,maxcon)
      real*8 ngradv(maxnd,ndim,ndim),nstrn(maxnd,ndim,ndim)
      real*8 nstrs(maxnd,ndim,ndim),nv(maxnd,ndim),ngradp(maxnd,ndim)
      real*8 bsf(maxele,nints,ndim),bphiJ(maxele,nints)
      real*8 bp(maxele,nints,ndim),pint(maxele,nints,2)
      real*8 Jacob(maxele,nints,2),Finv(maxele,nints,ndim,ndim,2)
      real*8 dJacob(maxele,nints),dphi(maxele,nints)
      real*8 duint(maxele,nints,ndim),bv(maxele,nints)
      real*8 Cstiff(4),c11,c33,c13,c44,w0,kperm
      real*8 dstrn(maxele,nints,ndim,ndim),Hpsi(maxele,nints)
      real*8 ngradmu(maxnd,ndim),COD,uint(maxele,nints,ndim)
      real*8 conAvpu(maxdbc,maxcon),conBvpu(maxdbc,maxcon)
      real*8 fT(maxnd*(ndpn)),fintT(maxnd*(ndpn))
      real*8 drhsT(maxnd*ndpn),pgauge(2)
      real*8 conA(maxdbc,maxcon),conB(maxdbc,maxcon)
      real*8 lamda,sdot1,sdot2,s,stot,sumsqftemp
      real*8 betafactor,a0,LsqFint(maxnd*ndpn*ndim)
      real*8 FiJ(maxele,nints,ndim,ndim),phi(maxele,nints,2)
      real*8 visc_eff(maxele,nints),dHds(maxele,nints,ndim,ndim)
      real*8 Hpsi0(maxele,nints),sigmaHmax,sigmahmin,p0
      real*8 biot(maxele),biot_temp
      real*8 YM(maxele),poisson(maxele)
      real*8 YM_temp,poisson_temp,freestrain_temp,phi0_temp
      real*8 Imperm_temp(ndim,ndim),fluidsource_temp
      real*8 freestrain(maxele),phi0(maxele)
      real*8 fluidsource(maxele),muavg(maxele),open_thresh
      real*8 KDtN(maxbdnd*ndim,maxbdnd*ndim),norm_u,norm_mu,sv
      real*8 Imperm(maxele,ndim,ndim),maxJ(maxele)
      real*8 drhsTproc(maxnd*ndpn),drhs_normalize(maxnd*ndpn)
      real*8 maxU_u,maxU_v,maxU_p,maxU_mu,sumsqu_n
      real*8 maxF_u,maxF_v,maxF_p,maxF_mu,sumsqf_n
      real*8 l_maxF_u,l_maxF_v,l_maxF_p,l_maxF_mu
      real*8 Wp(maxele,nints),Wp_total,Wmu(maxele,nints),Wmu_total
      real*8 Ws(maxele,nints),Ws_total,dmuint(maxele,nints)
      real*8 kappa,viscosity,timerule,Qrate,delta_t
      real*8 u_temp(maxnd,ndim),v_temp(maxnd,ndim),p_temp(maxnd)
      real*8 mu_temp(maxnd),tt,vflux,gmuint(maxele,nints,ndim)
      real*8 stretch(maxele,nints),Cf,Cphi,diverged,xx,yy,area_sum
      real*8 BtB(nnds,ndim,ndim),Hphi(maxele,nints),Hphi0(maxele,nints)
      real*8 dHphids(maxele,nints,ndim,ndim),dHphidp(maxele,nints)
      real*8 val,val2,mass_bal,dHdp(maxele,nints)
      integer nodes,neles,ele(maxele*maxprocs,nnds),nudofs
      integer umpcnd(maxdbc,2),nmudbcs,nwkmu,nconu
      integer nvpdbcs,nwkvp,indices(maxnd*ndpn*ndim)
      integer i,j,fbcnd(maxfc),fbcdof(maxfc),nudbcs,nfbcs,ii,m
      integer gconu(maxnd,ndpn),dof,umpc(maxdbc,maxcon,2)
      integer maxa(maxnd*ndim+1),nwk,isol,k,cntr,control,l
      integer dbcnd(maxdbc),dbcdof(maxdbc),icont,ibcs
      integer nmudofs,gconmu(maxnd,ndpn),maxamu(maxnd+1)
      integer umpcnum(maxdbc),nconmu,jj,mumpcnum(maxdbc)
      integer mumpc(maxdbc,maxcon,2),mumpcnd(maxdbc,2)
      integer nvpdofs,gconvp(maxnd,ndpn),maxavp(maxnd*(ndim+1)+1)
      integer nconvp,vpmpcnum(maxdbc),solsteps,stepone
      integer vpmpc(maxdbc,maxcon,2),vpmpcnd(maxdbc,2),rank,sizeof
      integer vpumpcnum(maxdbc),vpumpc(maxdbc,maxcon,2)
      integer vpumpcnd(maxdbc,2),ndofs,ndbcs
      integer injele,DtNnodes
      integer mpcnum(maxdbc),mpc(maxdbc,maxcon,2)
      integer mpcnd(maxdbc,2),gcon(maxnd,ndpn),betaOn
      integer InjIdx(maxnd*(ndpn)),nlamda,constraineles,constraindofs
      integer Istart,Iend,isu,isvp,ismu,zeroRowcount,zeroRow(maxnd*ndpn)
      integer ele_per_node(maxnd),ele_at_node(maxnd,10)
      integer DtNnode(maxbdnd*ndim),DtNdof(maxbdnd*ndim)
      integer DtNidx(maxbdnd*ndim),DtNdofs
      integer tag,sendcount,recvcount,Ivpustart,Ivpuend
      integer Imustart,Imuend,IGstart,IGend,indicesV(maxnd*ndpn*ndim)
      integer t_resets,imod
      integer l_ele(maxele*maxprocs),local_eles,g_ele(maxele)
      character(5) fnum
      character(1) ou1,ou2,ou3
      logical linesearch,conv1,conv2,convAll,makeau,moduli_on,DtN_on
      logical tc_split_on,scale_visc,eleproc(maxele*maxprocs),artstiff
      logical true_tstep,div_reset

      integer stat(MPI_STATUS_SIZE)

      PetscScalar get(1)
      PetscOffset i_get

      PetscInt low
      PetscInt high
      PetscInt otherlow
      PetscInt otherhigh
      PetscInt proc_count
      PetscInt lows(maxprocs)
      PetscInt proc_counts(maxprocs)

      PetscErrorCode ierr
      Vec PetR
      Vec PetRs
      Vec PetRnorm
      Vec PetdU
      Vec PetdUnorm
      Vec LsqF
      Vec LsqU

      Mat PetA
      Mat PetAvpu
      Mat PetAmu
!      Mat PetG

      KSP ksp
!      KSP kspG

      call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
      CHKERRQ(ierr)
      call MPI_COMM_RANK(PETSC_COMM_WORLD,rank,ierr)
      call MPI_COMM_SIZE(PETSC_COMM_WORLD,sizeof,ierr)

      print*, 'Process ',rank+1,' of ',sizeof,' is active'
      tag=9999
      diverged=1.d20

      write(ou1,910)'1'
      write(ou2,910)'2'
      write(ou3,910)'3'

      moduli_on=.true.
      DtN_on=.true.
      tc_split_on=.true.
      scale_visc=.true.

      if (DtN_on) then
        open(unit=45,file='input/DtN_dofs')
        read(45,*)DtNnodes,DtNdofs
        close(unit=45)
        if (DtNdofs.eq.0) DtN_on=.false.
      endif

      if (scale_visc) then
        sv=1.d0
      else
        sv=0.d0
      endif

      xold=-1.d9
      
***  e0 is the characteristic strain, e0 = sqrt(Gc / E / l0)
      e0=0.001d0

***  a0=-2.d0 for quadratic, -1.d-5 for cubic, >0.d0 for bell curve
      a0=-1.d-5

*** Normalized elastic properties for v = 0.25
      c11=1.2d0
      c13=0.4d0
      c33=1.2d0
      c44=0.4d0

!      c11=1.0d0
!      c13=0.0d0
!      c33=1.0d0
!      c44=0.5d0

      Cstiff(1)=c11
      Cstiff(2)=c33      
      Cstiff(3)=c13
      Cstiff(4)=c44
      
!     Normalized by sigma_0
      viscosity=1.d-10*1.d0
      
!     Normalized by (l_0)^2
      kappa=1.d-13

!     Fluid Compressibilty sigma_0/Kf
      Cf=0.01364d0*0.d0

      Qrate=1.d-1*1.d0
      delta_t=1.0d0
      
      if (scale_visc) then
!       Permeability constant (1 over small number from paper)
        kperm=1.d3
!       open_thresh = crit_stretch - 1
        open_thresh=1.d0/((1.d0+1.d0/(kappa*kperm))**(1.d0/3.d0)-1.d0)
      else
        kperm=1.d6
        open_thresh=1.d0/((1.d0+1.d0/(kappa*kperm))**(1.d0/3.d0)-1.d0)
      endif

      if (rank.eq.0) then
        open(unit=51,file='results1/params')
        write(51,*),'tc_split:',tc_split_on
        write(51,*),'scale_visc:',scale_visc
        write(51,*),'e0:',e0
        write(51,*),'viscosity:',viscosity
        write(51,*),'kappa:',kappa
        write(51,*),'Q:',Qrate
        write(51,*),'dt:',delta_t
        write(51,*),'a0:',a0
        close(unit=51)
      endif

!     Normalizer for v
      w0=(kperm*kappa)/viscosity

      if (rank.eq.0) print*,open_thresh,kperm,w0
      
      mu0=1.d0

      sigmahmin=-1.665d0*0.d0-0.d0
      sigmaHmax=-1.895d0*0.d0-0.d0
      p0=1.562d0*0.d0+0.0d0
      
      icont=1

      pgauge(1)=-20.25d0
      pgauge(2)=0.0d0

      call readinp(ndim,nnds,ndpn,maxnd,maxele,maxfc,maxdbc,
     &             nodes,neles,nudbcs,nfbcs,nudofs,gconu,ele,fbcnd,
     &             fbcdof,dbcnd,dbcdof,node,fbcval,dbcval,umpc,conAu,
     &             conBu,umpcnum,nconu,umpcnd,nmudofs,gconmu,
     &             nmudbcs,maxcon,nconmu,conAmu,conBmu,mumpcnum,
     &             mumpc,mumpcnd,nvpdofs,gconvp,nvpdbcs,nconvp,
     &             conAvp,conBvp,vpmpcnum,vpmpc,vpmpcnd,conAvpu,
     &             conBvpu,vpumpcnum,vpumpc,vpumpcnd,
     &             conA,conB,mpcnum,mpc,mpcnd,gcon,freestrain,
     &             phi0,biot,YM,poisson,fluidsource,Imperm,
     &             moduli_on,rank,maxprocs)

      do 201 i=1,nodes
        ele_per_node(i)=0        
        do 202 j=1,neles
          do 205 m=1,4
            if (i.eq.ele(j,m)) then
              ele_per_node(i)=ele_per_node(i)+1
              if (ele_per_node(i).gt.10) then
                if (rank.eq.0) print*,'Increase ele_at_node size!'
                stop
              endif
              ele_at_node(i,ele_per_node(i))=j
            endif
 205      continue
 202    continue      
 201  continue

      ndofs=nudofs+nvpdofs+nmudofs
      ndbcs=nudbcs+nvpdbcs+nmudbcs

      if (rank.eq.0) then
      print*,'********************************'
      print*,'****** INPUT READ SUCCESS ******'
      print*,'********************************'
      print*
      print*,'======================================='
      print*,'|  nudofs=',nudofs
      print*,'|  nvpdofs=',nvpdofs
      print*,'|  nmudofs=',nmudofs
      print*,'|  total=',ndofs
      print*,'======================================='
      endif

********************************************************
************** Declare PETsc components ****************
********************************************************

      call MatCreate(PETSC_COMM_WORLD,PetA,ierr)
      call MatSetSizes(PetA,PETSC_DECIDE,
     &                      PETSC_DECIDE,
     &                      ndofs,
     &                      ndofs,ierr)
      if (sizeof.eq.1) then
         call MatSetType(PetA,MATSEQAIJ,ierr)
      else
         call MatSetType(PetA,MATMPIAIJ,ierr)
      endif
      call MatMPIAIJSetPreallocation(PetA,
     &        (ndpn)*nnds*10,PETSC_NULL_INTEGER,
     &        (ndpn)*nnds*5,PETSC_NULL_INTEGER,ierr)
      call MatSeqAIJSetPreallocation(PetA,
     &       (ndpn)*nnds*10,PETSC_NULL_INTEGER,ierr)
!      call MatCreateSeqAIJ(PETSC_COMM_SELF,
!     &       nvpdofs+nudofs+nmudofs+1,
!     &       nvpdofs+nudofs+nmudofs+1,
!     &       (ndpn)*nnds*7,PETSC_NULL_INTEGER,
!     &       PetA,ierr)
!      call MatCreate(PETSC_COMM_WORLD,PetA,ierr)
!      call MatSetSizes(PetA,PETSC_DECIDE,
!     &                      PETSC_DECIDE,
!     &           nudofs+nvpdofs+nmudofs,
!     &           nudofs+nvpdofs+nmudofs,ierr)
      call MatSetFromOptions(PetA,ierr)
      call MatSetUp(PetA,ierr)
      call MatGetOwnershipRange(PetA,Istart,Iend,ierr)


      call MatCreate(PETSC_COMM_WORLD,PetAvpu,ierr)
      call MatSetSizes(PetAvpu,PETSC_DECIDE,
     &                      PETSC_DECIDE,
     &                      ndofs,
     &                      ndofs,ierr)
      if (sizeof.eq.1) then
         call MatSetType(PetAvpu,MATSEQAIJ,ierr)
      else
         call MatSetType(PetAvpu,MATMPIAIJ,ierr)
      endif
      call MatMPIAIJSetPreallocation(PetAvpu,
     &        (ndpn)*nnds*10,PETSC_NULL_INTEGER,
     &        (ndpn)*nnds*5,PETSC_NULL_INTEGER,ierr)
      call MatSeqAIJSetPreallocation(PetAvpu,
     &       (ndpn)*nnds*10,PETSC_NULL_INTEGER,ierr)
      call MatSetFromOptions(PetAvpu,ierr)
      call MatSetUp(PetAvpu,ierr)
      call MatGetOwnershipRange(PetAvpu,Ivpustart,Ivpuend,ierr)


      call MatCreate(PETSC_COMM_WORLD,PetAmu,ierr)
      call MatSetSizes(PetAmu,PETSC_DECIDE,
     &                      PETSC_DECIDE,
     &                      ndofs,
     &                      ndofs,ierr)
      if (sizeof.eq.1) then
         call MatSetType(PetAmu,MATSEQAIJ,ierr)
      else
         call MatSetType(PetAmu,MATMPIAIJ,ierr)
      endif
      call MatMPIAIJSetPreallocation(PetAmu,
     &        (ndpn)*nnds*10,PETSC_NULL_INTEGER,
     &        (ndpn)*nnds*5,PETSC_NULL_INTEGER,ierr)
      call MatSeqAIJSetPreallocation(PetAmu,
     &       (ndpn)*nnds*10,PETSC_NULL_INTEGER,ierr)
      call MatSetFromOptions(PetAmu,ierr)
      call MatSetUp(PetAmu,ierr)
      call MatGetOwnershipRange(PetAmu,Imustart,Imuend,ierr)


!  Create vectors.  Note that we form 1 vector from scratch and
!  then duplicate as needed.
      call MatCreateVecs(PetA,PETSC_NULL_OBJECT,PetR,ierr)

!      call VecCreate(PETSC_COMM_WORLD,PetR,ierr)
!      call VecSetSizes(PetR,PETSC_DECIDE,
!     &             nudofs+nvpdofs+nmudofs+1,ierr)
      call VecSetFromOptions(PetR,ierr)

      call VecDuplicate(PetR,PetdU,ierr)
      call VecDuplicate(PetR,PetRs,ierr)
      call VecDuplicate(PetR,PetRnorm,ierr)
      call VecDuplicate(PetR,PetdUnorm,ierr)

      call KSPCreate(PETSC_COMM_WORLD,ksp,ierr)


!      call MatCreate(PETSC_COMM_WORLD,PetG,ierr)
!      call MatSetSizes(PetG,PETSC_DECIDE,
!     &                      PETSC_DECIDE,
!     &      ndim*(ndofs+ndbcs),
!     &      ndim*(ndofs+ndbcs),ierr)
!      if (sizeof.eq.1) then
!         call MatSetType(PetG,MATSEQAIJ,ierr)
!      else
!         call MatSetType(PetG,MATMPIAIJ,ierr)
!      endif
!      call MatMPIAIJSetPreallocation(PetG,
!     &        (ndpn)*nnds*ndim*7,PETSC_NULL_INTEGER,
!     &        (ndpn)*nnds*ndim*3,PETSC_NULL_INTEGER,ierr)
!      call MatSeqAIJSetPreallocation(PetG,
!     &       (ndpn)*nnds*ndim*7,PETSC_NULL_INTEGER,ierr)
!
!      call MatSetFromOptions(PetG,ierr)
!      call MatSetUp(PetG,ierr)
!      call MatGetOwnershipRange(PetG,IGstart,IGend,ierr)
!
!!  Create vectors.  Note that we form 1 vector from scratch and
!!  then duplicate as needed.
!      call MatCreateVecs(PetG,PETSC_NULL_OBJECT,LsqF, ierr )
!
!      call VecSetFromOptions(LsqF,ierr)
!
!      call VecDuplicate(LsqF,LsqU,ierr)
!
!      call KSPCreate(PETSC_COMM_WORLD,kspG,ierr)

********************************************************
************** Finished PETsc components ***************
********************************************************

      local_eles=0
      do 210 i=1,neles
        eleproc(i)=.true.
        do 211 j=1,nnds
          do 212 k=1,ndpn
            dof=gcon(ele(i,j),k)-1
            if (dof.le.ndofs-1) then
              if (dof.ge.Iend) eleproc(i)=.false.
              if (dof.lt.Istart) eleproc(i)=.false.
            endif
 212      continue
 211    continue
        if (.not.eleproc(i)) then
          do 213 j=1,nnds
            do 214 k=1,ndpn
              dof=gcon(ele(i,j),k)-1
              if (dof.lt.Iend) eleproc(i)=.true.
 214        continue
 213      continue
          do 215 j=1,nnds
            do 216 k=1,ndpn
              dof=gcon(ele(i,j),k)-1
              if (dof.lt.Istart) eleproc(i)=.false.
 216        continue
 215      continue
        endif
        if (eleproc(i)) then
          local_eles=local_eles+1
          g_ele(local_eles)=i
          l_ele(i)=local_eles
          if (local_eles.gt.maxele) then
            print*,'*** TOO MANY ELEMENTS ON PROC: ',rank
            stop
          endif
        endif
 210  continue

      print*,'|  Process: ',rank,' Elements: ',local_eles

      imod=20
      if (moduli_on) then
       open(unit=imod,file='input/moduli')
       do 501 i=1,neles
         if (ndim.eq.2) then
           read(imod,*)k,YM_temp,biot_temp,phi0_temp,
     &       poisson_temp,freestrain_temp,fluidsource_temp,
     &       Imperm_temp(1,1),Imperm_temp(1,2),
     &       Imperm_temp(2,1),Imperm_temp(2,2)
           if (eleproc(i)) then
             YM(l_ele(i))=YM_temp
             biot(l_ele(i))=biot_temp
             phi0(l_ele(i))=phi0_temp
             poisson(l_ele(i))=poisson_temp
             freestrain(l_ele(i))=freestrain_temp
             fluidsource(l_ele(i))=fluidsource_temp
             do 505 j=1,ndim
               do 506 k=1,ndim
                 Imperm(l_ele(i),j,k)=Imperm_temp(j,k)
 506           continue
 505         continue
           endif
         else
           print*,
     &      ' ERROR: Readinp does not support number of dimensions!'
           stop
         endif
 501   continue
       close(unit=imod)
      else
       do 502 i=1,local_eles
         YM(i)=1.d0
         biot(i)=0.d0
         phi0(i)=0.d0
         poisson(i)=0.25d0*1.d0
         freestrain(i)=0.d0
         fluidsource(i)=0.d0
         do 503 j=1,ndim
           do 504 k=1,ndim
             Imperm(i,j,k)=0.d0
             if (j.eq.k) Imperm(i,j,k)=1.d0
 504       continue
 503     continue
 502   continue
      endif

      call calcB(ndim,nnds,nints,maxnd,maxele,local_eles,ele,
     &           node,det,B,NN,g_ele,maxprocs)

!!    CONSTRUCT LEAST_SQUARES MATRIX FOR CALCULATING NODAL GRADIENTS
!
!      call stfLsq(nnds,nints,ndim,ndpn,B,NN,det,node,maxele,
!     &            maxnd,maxcon,maxdbc,nvpdofs,neles,ele,gcon,mpc,
!     &            mpcnum,conA,nvpdbcs,nmudofs,nmudbcs,nudofs,
!     &            nudbcs,PetG,rank,eleproc)
!
!      call MatAssemblyBegin(PetG,MAT_FINAL_ASSEMBLY,ierr)
!      call MatAssemblyEnd(PetG,MAT_FINAL_ASSEMBLY,ierr)
!!      call MatView(PetG,PETSC_VIEWER_STDOUT_WORLD,ierr) 
!
!      call KSPSetOperators(kspG,PetG,PetG,ierr)
!!      call KSPSetFromOptions(kspG,ierr)
!      call KSPSetType(kspG,KSPGMRES,ierr)

!    DONE CONSTRUCTING

      DtNdofs=0
      do 703 i=1,maxbdnd*ndim
        DtNidx(i)=-1
        do 704 j=1,maxbdnd*ndim
          KDtN(i,j)=0.d0
 704    continue
 703  continue

      if (DtN_on)
     &  call DtNstf(ndim,maxnd,maxbdnd,ndpn,node,gcon,KDtN,
     &                  DtNidx,DtNdofs,DtNnode,DtNdof)

*     initialize the nodal variables

      open(unit=25,file='input/initial')
      do 20 i=1,nodes
        u(i,1)=0.d0
        u(i,2)=0.d0
        v(i,1)=0.0d0
        v(i,2)=0.d0
        mu(i)=mu0
        p(i)=-0.0d0*(-1.d0+node(i,1))
        if (1.eq.0) then
          read(25,*) x,y,u(i,1),u(i,2),mu(i),p(i),v(i,1),v(i,2)
          u(i,1)=u(i,1)/e0
          u(i,2)=u(i,2)/e0
          v(i,1)=v(i,1)/w0
          v(i,2)=v(i,2)/w0
        endif
        u_temp(i,1)=u(i,1)
        u_temp(i,2)=u(i,2)
        v_temp(i,1)=v(i,1)
        v_temp(i,2)=v(i,2)
        mu_temp(i)=mu(i)
        p_temp(i)=p(i)        
 20   continue
      close(unit=25)

      dt=1.d0

      call gradvp(nnds,nints,ndim,maxele,maxnd,local_eles,
     &              ele,B,NN,gradv,vint,gradp,pint,v,p,g_ele,
     &              maxprocs)

      call gradm(nnds,nints,ndim,maxele,maxnd,local_eles,ele,B,
     &             NN,muint,gradmu,mu,dt,muavg,dmuint,g_ele,
     &             maxprocs)

      call gradu(ndim,nnds,nints,maxele,maxnd,local_eles,ele,
     &               B,NN,strn,u,Jacob,Finv,dJacob,dstrn,e0,dt,uint,
     &               duint,FiJ,muint,phi,open_thresh,g_ele,maxprocs,
     &               maxJ)

      call nodalstress(ndim,maxnd,maxele,nints,nnds,nodes,neles,
     &                     ele,gradp,gradmu,strn,strs,ngradp,ngradmu,
     &                     nstrn,nstrs,ele_per_node,ele_at_node,
     &                     eleproc,maxprocs,l_ele,det)

      call normmu(nnds,nints,ndim,maxele,maxnd,local_eles,
     &                 ele,B,NN,ngradmu,gmuint,g_ele,maxprocs)

!      do 22 i=1,neles
!        do 23 j=1,nints
!          Hpsi(i,j)=-1.d0
!          phi(i,j,1)=phi0(i)
! 23     continue
! 22   continue

      do 732 i=1,local_eles
        do 733 j=1,nints
          Hpsi0(i,j)=Hpsi(i,j)
          Jacob(i,j,2)=Jacob(i,j,1)
          muint(i,j,2)=muint(i,j,1)
          phi(i,j,2)=phi(i,j,1)
          pint(i,j,2)=pint(i,j,1)
          do 734 k=1,ndim
            do 735 l=1,ndim
              strn(i,j,k,l,2)=strn(i,j,k,l,1)
              gradv(i,j,k,l,2)=gradv(i,j,k,l,1)
              Finv(i,j,k,l,2)=Finv(i,j,k,l,1)
 735        continue
 734      continue
 733    continue
 732  continue

      call stresses(nints,ndim,maxele,local_eles,strn,gradv,
     &                    gradmu,pint,gradp,vint,strs,tau,bsf,bp,bv,
     &                    psimu,ksi,phi0,Jacob,Finv,Cstiff,w0,
     &                    Imperm,dJacob,dphi,dstrn,Hpsi,e0,duint,
     &                    isol,a0,phi,visc_eff,muint,g_ele,bphiJ,
     &                    dHds,Hpsi0,biot,dt,freestrain,YM,poisson,
     &                    tc_split_on,muavg,kperm,Wp,uint,Wmu,Ws,
     &                    dmuint,fluidsource,rank,open_thresh,
     &                    gmuint,scale_visc,stretch,maxprocs,Cf,
     &                    p0,maxJ,dHdp)

      do 30 i=1,nodes*ndpn
        indices(i)=i-1
        if (i-1.lt.Istart.or.i-1.ge.Iend) then
          indicesV(i)=-1
        else
          indicesV(i)=indices(i)
        endif
        drhs(i)=0.d0
 30   continue
      do 31 i=1,nodes*ndpn
        f(i)=0.d0
 31   continue

      do 556 i=1,nudbcs+nmudbcs+nvpdbcs
        unitdbc(i)=dbcval(i)
 556  continue
      do 557 i=1,nfbcs
        unitfbc(i)=fbcval(i)
 557  continue

      t=0.0d0
      stepone=1

!      if (rank.eq.0) then
!        open(unit=50,file='results1/crack_vs_time')
!        write(50,*),'t crack_length pressure opening'
!        close(unit=50)
!      endif

!      open(unit=50,file='results1/consolidation')
!      write(50,*),'t displacement pressure(right_side)'

!      open(unit=50,file='results1/pressure')
!      write(50,*),'time pressure'

      open(unit=50,file='results1/flux')

      open(unit=24,file='input/crackloc')
      read(24,*) crackloc
      close(unit=24)
     
      do 777 ibcs=0,200

        if (rank.eq.0) then
        print*
        print*,'======================================='
        print*,'|  LOAD STEP INFORMATION'
        print*,'|  ibcs=',ibcs
        endif

        if (ibcs.le.0) then
          muapp=1.0d0
          papp=0.d0+p0+0.d0*(0.9d0+1.d-2*dble(ibcs))
     &         /dsqrt(crackloc*dacos(-1.d0)*(1.d0-1.d0*0.25d0*0.25d0))
          uapp=0.d0
          vapp=0.d0
          makeau=.true.
          true_tstep=.true.
          timerule=1.0d0
        elseif (ibcs.le.0) then
          muapp=1.0d0
          papp=0.d0+p0
          uapp=1.d-3*dble(ibcs)
          vapp=Qrate*0.d-4
          makeau=.true.
          true_tstep=.true.
          timerule=1.d0
        elseif (ibcs.le.1000) then 
          muapp=1.0d0
          papp=0.d0+p0
          uapp=0.d-5*2.d0**dble(ibcs-0)
          vapp=Qrate*1.d0
          makeau=.true.
          true_tstep=.true.
          timerule=1.d0
        else 
          muapp=1.0d0
          papp=0.d0+p0
          uapp=1.d2
          vapp=Qrate*0.d0
          makeau=.true.
          true_tstep=.true.
          timerule=1.d0
        endif

        if (rank.eq.0) 
     &  print*,'|  ~ Opening Threshold: ',open_thresh

        if (ibcs.le.0) then
          dt=1.d30
        elseif (ibcs.le.0) then
          dt=delta_t*2.d4
        elseif (ibcs.le.1000) then
          dt=2.d0**(dble(ibcs)*0.15d0+1.d0)
     &       -2.d0**(dble(ibcs-1)*0.15d0+1.d0)
          dt=delta_t
        else
          dt=1.d20
        endif
 
        if (ibcs.le.0) then
          t=0.d0
        else
          t=t+dt
        endif

        if (ibcs.eq.-1) t=t+2.d0*delta_t-delta_t*2.d4

        t_resets=0

        if (rank.eq.0) then
        print*,'|  muapp =',muapp
        print*,'|  papp   =',papp
        print*,'|  uapp   =',uapp
        print*,'|  vapp   =',vapp
        print*,'|  t      =',t
        print*,'|  Time Integration Rule:',timerule
        print*,'======================================='
        print*
        print*,'======================================='
        print*,'|  SOLVER STATUS...'
        print*,'|'
        endif
 
        do 778 i=1,nfbcs
         if (fbcdof(i).eq.ndim) then
           fbcval(i)=unitfbc(i)*sigmaHmax
         elseif (fbcdof(i).eq.1) then
           fbcval(i)=unitfbc(i)*sigmahmin
         endif
 778    continue
 
        do 779 i=1,nudbcs+nmudbcs+nvpdbcs
         if (dbcdof(i).eq.ndpn-1) then
           dbcval(i)=unitdbc(i)*papp
         else
           if (dbcdof(i).le.2) then
             dbcval(i)=unitdbc(i)*uapp/e0
           else if (dbcdof(i).eq.6) then
             dbcval(i)=unitdbc(i)*muapp
           else
             dbcval(i)=unitdbc(i)*vapp/w0
           endif
         endif
 779    continue
        
****    global control loop
        betaOn=0
        beta=0.d0
        tt=t-dt
        div_reset=.false.
        do 900 control=0,csteps

          if (control.eq.csteps) dt=t-tt
          tt=tt+dt

          if (true_tstep) then
          if (.not.div_reset) then
            do 722 i=1,local_eles
              do 723 j=1,nints
                Hpsi0(i,j)=Hpsi(i,j)
                Jacob(i,j,2)=Jacob(i,j,1)
                muint(i,j,2)=muint(i,j,1)
                phi(i,j,2)=phi(i,j,1)
                pint(i,j,2)=pint(i,j,1)
                do 724 k=1,ndim
                  do 725 l=1,ndim
                    strn(i,j,k,l,2)=strn(i,j,k,l,1)
                    gradv(i,j,k,l,2)=gradv(i,j,k,l,1)
                    Finv(i,j,k,l,2)=Finv(i,j,k,l,1)
 725              continue
 724            continue
 723          continue
 722        continue
          endif
          endif

          if (div_reset) then
            do 422 i=1,local_eles
              do 423 j=1,nints
                Hpsi(i,j)=Hpsi0(i,j)
 423          continue
 422        continue
            do 120 i=1,nodes
              u(i,1)=u_temp(i,1)
              u(i,2)=u_temp(i,2)
              v(i,1)=v_temp(i,1)
              v(i,2)=v_temp(i,2)
              mu(i)=mu_temp(i)
              p(i)=p_temp(i)
 120        continue
          endif
  
          if (rank.eq.0) then
          print*,'|'
          print*,'|  CONTORL STEP: ',control
          print*,'|  Time Resets: ',t_resets
          print*,'|  t=',t,'tt=',tt
          print*,'|'
          endif

          do 600 i=1,ndpn*nodes
            drhs(i)=0.d0
 600      continue

          call calcvp(ndpn,ndim,maxnd,maxdbc,maxcon,nodes,nvpdofs,
     &              gconvp,dbcval,drhs,v,p,nudbcs,dbcnd,dbcdof,
     &              vpmpcnum,nvpdbcs,nconvp,vpmpc,vpmpcnd,
     &              conAvp,conBvp,nmudbcs)
     
          call gradvp(nnds,nints,ndim,maxele,maxnd,local_eles,
     &              ele,B,NN,gradv,vint,gradp,pint,v,p,g_ele,
     &              maxprocs)

          call calcmu(ndpn,ndim,maxnd,maxdbc,maxcon,nodes,nmudofs,
     &              gconmu,dbcval,drhs,mu,nudbcs,dbcnd,dbcdof,
     &              mumpcnum,nmudbcs,nconmu,mumpc,mumpcnd,
     &              conAmu,conBmu,nvpdbcs)
     
          call gradm(nnds,nints,ndim,maxele,maxnd,local_eles,ele,B,
     &              NN,muint,gradmu,mu,dt,muavg,dmuint,g_ele,
     &              maxprocs)

          call calcu(ndim,ndpn,maxnd,maxdbc,nodes,nudofs,nconu,
     &             gconu,dbcval,drhs,u,umpc,nudbcs,dbcnd,dbcdof,
     &             conAu,conBu,umpcnd,nmudbcs,nvpdbcs,umpcnum,maxcon)

          call gradu(ndim,nnds,nints,maxele,maxnd,local_eles,ele,
     &               B,NN,strn,u,Jacob,Finv,dJacob,dstrn,e0,dt,uint,
     &               duint,FiJ,muint,phi,open_thresh,g_ele,maxprocs,
     &               maxJ)

          call nodalstress(ndim,maxnd,maxele,nints,nnds,nodes,neles,
     &                     ele,gradp,gradmu,strn,strs,ngradp,ngradmu,
     &                     nstrn,nstrs,ele_per_node,ele_at_node,
     &                     eleproc,maxprocs,l_ele,det)

          call normmu(nnds,nints,ndim,maxele,maxnd,local_eles,
     &                 ele,B,NN,ngradmu,gmuint,g_ele,maxprocs)

          if (rank.eq.0) then
          do 602 i=1,nodes
            if (dabs(node(i,1)-pgauge(1)).lt.1.d-4) then
            if (dabs(node(i,2)-pgauge(2)).lt.1.d-4) then
              print*,'|        PRESSURE GAUGE: ',p(i)
            endif
            endif
 602      continue
          endif

          if (ibcs.eq.0.and.rank.eq.0) then             
          call paraview(0,ou2,maxnd,maxele,ndim,nodes,neles,nnds,
     &                  node,ele,u,mu,v,p,ngradp,nstrn,nstrs,e0,
     &                  ngradmu,nints,NN,phi,w0,maxprocs)
          endif
              
          if (control.eq.1) then
            solsteps=maxsol/2
          else
            solsteps=maxsol/2
          endif
          solsteps=200

          if (ibcs.ge.0) then
            conv1=.true.
            conv2=.true.
          else
            conv1=.false.
            conv2=.false.
          endif
          convAll=.false.
          do 100 isol=1,solsteps

            if (conv1.and.conv2) then
              isu=1
              isvp=1
              ismu=1
              if (ibcs.eq.-1) ismu=0
              if (rank.eq.0) then
              print*,'|  |-------------------------'
              print*,'|  |  isol=',(isol+1)/2
              print*,'|  |  DOFs : mu,u,vp'
              print*,'|  |-------------------------'
              endif
              convAll=.true.
            elseif (mod(isol,2).eq.0) then
              isu=0
              isvp=0
              ismu=1
              if (rank.eq.0) then
              print*,'|  |-------------------------'
              print*,'|  |  isol=',(isol+1)/2
              print*,'|  |  DOFs : mu'
              print*,'|  |-------------------------'
              endif
            elseif (mod(isol,2).eq.1) then
              isu=0
              if (makeau) isu=1
              isvp=1
              ismu=0
              if (rank.eq.0) then
              print*,'|  |-------------------------'
              print*,'|  |  isol=',(isol+1)/2
              print*,'|  |  DOFs : u,vp'
              print*,'|  |-------------------------'
              endif
            else
              if (rank.eq.0) print*,'You done messed up!!!'
              stop
            endif

            l_maxF_u=1.d0
            l_maxF_v=1.d0
            l_maxF_p=1.d0
            l_maxF_mu=1.d0

            maxU_u=0.d0
            maxU_v=0.d0
            maxU_p=0.d0
            maxU_mu=0.d0

            call applybcu(ndpn,ndim,maxfc,fbcnd,fbcdof,maxnd,nnds,
     &                    nfbcs,gconu,fbcval,f,nudofs,nodes,ele,
     &                    local_eles,g_ele,maxprocs,rank,maxele)

            do 532 i=1,nodes
              do 533 j=1,ndim
                fT(gcon(i,j))=f(gconu(i,j))*dble(isu)
                if (dabs(fT(gcon(i,j))).gt.l_maxF_u) 
     &            l_maxF_u=dabs(fT(gcon(i,j)))
 533          continue
 532        continue

            call applybcmu(nnds,nints,ndpn,ndim,maxfc,fbcnd,fbcdof,
     &                   maxnd,nodes,node,nfbcs,gconmu,fbcval,f,dt,
     &                   beta,maxele,neles,ele,det,NN,muint,nmudofs,
     &                   local_eles,g_ele,maxprocs,rank)

            do 534 i=1,nodes
              fT(gcon(i,6))=f(gconmu(i,1))*dble(ismu)
              if (dabs(fT(gcon(i,6))).gt.l_maxF_mu) 
     &          l_maxF_mu=dabs(fT(gcon(i,6)))
 534        continue
            
            do 101 ii=1,maxsol

            call VecSet(PetR,0.d0,ierr)
            call VecSet(PetRnorm,0.d0,ierr)
            call MPI_Barrier(MPI_COMM_WORLD,ierr)

            call applybcvp(nnds,nints,ndpn,ndim,maxfc,fbcnd,fbcdof,
     &                   maxnd,nodes,node,nfbcs,gconvp,fbcval,f,dt,
     &                   maxele,neles,ele,det,NN,mu,nmudofs,u,papp,e0,
     &                   unitfbc,vapp,fluidsource,w0,gradv,Finv,Jacob,
     &                   timerule,local_eles,g_ele,maxprocs,rank)

            do 530 i=1,nodes
              do 531 j=1,ndim+1
                fT(gcon(i,j+ndim))=f(gconvp(i,j))*dble(isvp)
                if (j.le.ndim) then
                  if (dabs(fT(gcon(i,j+ndim))).gt.l_maxF_v) 
     &              l_maxF_v=dabs(fT(gcon(i,j+ndim)))
                else
                  if (dabs(fT(gcon(i,j+ndim))).gt.l_maxF_p) 
     &              l_maxF_p=dabs(fT(gcon(i,j+ndim)))
                endif
 531          continue
 530        continue

            call stresses(nints,ndim,maxele,local_eles,strn,gradv,
     &                    gradmu,pint,gradp,vint,strs,tau,bsf,bp,bv,
     &                    psimu,ksi,phi0,Jacob,Finv,Cstiff,w0,
     &                    Imperm,dJacob,dphi,dstrn,Hpsi,e0,duint,
     &                    isol,a0,phi,visc_eff,muint,g_ele,bphiJ,
     &                    dHds,Hpsi0,biot,dt,freestrain,YM,poisson,
     &                    tc_split_on,muavg,kperm,Wp,uint,Wmu,Ws,
     &                    dmuint,fluidsource,rank,open_thresh,
     &                    gmuint,scale_visc,stretch,maxprocs,Cf,
     &                    p0,maxJ,dHdp)

            call intfvp(nnds,nints,ndim,ndpn,maxnd,nodes,maxele,
     &                  ele,gconvp,det,B,NN,bv,bp,tau,fint,dt,muint,
     &                  mu,bsf,Finv,Jacob,gradv,phi,w0,bphiJ,
     &                  timerule,local_eles,g_ele,maxprocs)

            do 540 i=1,nodes
              do 541 j=1,ndim+1
                fintT(gcon(i,j+ndim))=fint(gconvp(i,j))*dble(isvp)
                if (j.le.ndim) then
                  if (dabs(fintT(gcon(i,j+ndim))).gt.l_maxF_v) 
     &              l_maxF_v=dabs(fintT(gcon(i,j+ndim)))
                else
                  if (dabs(fintT(gcon(i,j+ndim))).gt.l_maxF_p) 
     &              l_maxF_p=dabs(fintT(gcon(i,j+ndim)))
                endif
 541          continue
 540        continue

            call intfu(nints,ndim,ndpn,nnds,maxnd,nodes,maxele,
     &                 ele,gconu,det,B,NN,strs,muint,pint,fint,e0,
     &                 gradmu,local_eles,g_ele,maxprocs)

!           DtN Forces 
            if (rank.eq.0) then
            do 70 i=1,DtNdofs
              do 71 j=1,DtNdofs
                fint(gconu(DtNnode(i),DtNdof(i)))=
     &            fint(gconu(DtNnode(i),DtNdof(i)))
     &            +KDtN(i,j)*u(DtNnode(j),DtNdof(j))                
 71           continue
 70         continue
            endif
!           End DtN Forces

            do 542 i=1,nodes
              do 543 j=1,ndim
                fintT(gcon(i,j))=fint(gconu(i,j))*dble(isu)
                if (dabs(fintT(gcon(i,j))).gt.l_maxF_u) 
     &            l_maxF_u=dabs(fintT(gcon(i,j)))
 543          continue
 542        continue

            call intfmu(nnds,nints,ndim,ndpn,maxnd,nodes,maxele,
     &          ele,gconmu,det,B,NN,psimu,ksi,fint,dt,beta,muint,
     &          local_eles,g_ele,maxprocs)

            do 544 i=1,nodes
              fintT(gcon(i,6))=fint(gconmu(i,1))*dble(ismu)
              if (dabs(fintT(gcon(i,6))).gt.l_maxF_mu) 
     &          l_maxF_mu=dabs(fintT(gcon(i,6)))
 544        continue
     
            do 590 i=1,nodes*ndpn
              drhs(i)=fT(i)-fintT(i)
 590        continue

            call MPI_Allreduce(l_maxF_u,maxF_u,1,MPI_REAL8,
     &                       MPI_MAX,MPI_COMM_WORLD,ierr)
            call MPI_Allreduce(l_maxF_v,maxF_v,1,MPI_REAL8,
     &                       MPI_MAX,MPI_COMM_WORLD,ierr)
            call MPI_Allreduce(l_maxF_p,maxF_p,1,MPI_REAL8,
     &                       MPI_MAX,MPI_COMM_WORLD,ierr)
            call MPI_Allreduce(l_maxF_mu,maxF_mu,1,MPI_REAL8,
     &                       MPI_MAX,MPI_COMM_WORLD,ierr)

            if (maxF_u.lt.1.d-14) maxF_u=1.d0
            if (maxF_v.lt.1.d-14) maxF_v=1.d0
            if (maxF_p.lt.1.d-14) maxF_p=1.d0
            maxF_mu=1.d0
!            if (maxF_mu.lt.1.d-14) maxF_mu=1.d0

            do 591 i=1,nodes
              do 592 j=1,ndpn
                dof=gcon(i,j)
                if (j.le.ndim) then
                  drhs_normalize(dof)=(fT(dof)-fintT(dof))/maxF_u
                elseif (j.le.ndim+ndim) then
                  drhs_normalize(dof)=(fT(dof)-fintT(dof))/maxF_v
                elseif (j.eq.ndim+ndim+1) then
                  drhs_normalize(dof)=(fT(dof)-fintT(dof))/maxF_p
                else
                  drhs_normalize(dof)=(fT(dof)-fintT(dof))/maxF_mu
                endif
 592          continue
 591        continue
 
            do 596 i=1,nconvp+nconu+nconmu
             j=nvpdofs+nvpdbcs+nudofs+nudbcs+nmudofs+nmudbcs+i
             do 597 jj=1,mpcnum(i)
               k=gcon(mpc(i,jj,1),mpc(i,jj,2))
               drhs(k)=drhs(k)+drhs(j)*conA(i,jj)
               drhs_normalize(k)=drhs_normalize(k)
     &                          +drhs_normalize(j)*conA(i,jj)
 597         continue
 596        continue


            call VecSetValues(PetR,nvpdofs+nudofs+nmudofs,
     &                     indices,drhs,ADD_VALUES,ierr)

            call VecSetValues(PetRnorm,nvpdofs+nudofs+nmudofs,
     &           indices,drhs_normalize,ADD_VALUES,ierr)

            call VecAssemblyBegin(PetR,ierr)
            call VecAssemblyEnd(PetR,ierr)

            call VecAssemblyBegin(PetRnorm,ierr)
            call VecAssemblyEnd(PetRnorm,ierr)
            
            if (ibcs.eq.-1.and.ii.eq.1)
     &      call VecView(PetR,PETSC_VIEWER_STDOUT_WORLD,ierr)
            call VecNorm(PetR,NORM_2,sumsqf,ierr)
            call VecNorm(PetRnorm,NORM_2,sumsqf_n,ierr)

            if (rank.eq.0) then
            print*,'|         ii=',ii
            print*,'|           sumsqf=',sumsqf,
     &             'relative=',sumsqf_n
            endif

            if (sumsqf.lt.1.d-10.and.ii.gt.1) then
              call VecSet(PetdU,0.d0,ierr)
              goto 911
            elseif (sumsqf.gt.diverged) then
              call VecSet(PetdU,0.d0,ierr)
              goto 912
            endif

            artstiff=.false.
            if (ii.le.10) artstiff=.false.

            if (ismu.eq.1) then
              if (isu.eq.1) then
 
                call stf(nnds,nints,ndim,B,NN,det,node,maxele,
     &          maxnd,maxcon,maxdbc,nvpdofs,ele,dt,Cstiff,
     &          gcon,mpc,mpcnum,conA,nvpdbcs,nmudofs,nmudbcs,
     &          nudofs,nudbcs,muint,Jacob,Finv,Imperm,w0,
     &          gradv,pint,tau,PetA,e0,dstrn,bsf,vint,strn,
     &          Hpsi,dJacob,dphi,phi,beta,a0,isu,isvp,ismu,
     &          PetR,PetdU,ksp,nodes,indices,gradp,dHdp,
     &          visc_eff,dHds,local_eles,p0,Cf,fluidsource,
     &          biot,phi0,freestrain,YM,poisson,KDtN,maxbdnd,
     &          DtNdofs,DtNidx,tc_split_on,muavg,open_thresh,
     &          kperm,rank,eleproc,ndpn,artstiff,dmuint,timerule,
     &          gmuint,stretch,g_ele,maxprocs)

              else

                call stf(nnds,nints,ndim,B,NN,det,node,maxele,
     &          maxnd,maxcon,maxdbc,nvpdofs,ele,dt,Cstiff,
     &          gcon,mpc,mpcnum,conA,nvpdbcs,nmudofs,nmudbcs,
     &          nudofs,nudbcs,muint,Jacob,Finv,Imperm,w0,
     &          gradv,pint,tau,PetAmu,e0,dstrn,bsf,vint,strn,
     &          Hpsi,dJacob,dphi,phi,beta,a0,isu,isvp,ismu,
     &          PetR,PetdU,ksp,nodes,indices,gradp,dHdp,
     &          visc_eff,dHds,local_eles,p0,Cf,fluidsource,
     &          biot,phi0,freestrain,YM,poisson,KDtN,maxbdnd,
     &          DtNdofs,DtNidx,tc_split_on,muavg,open_thresh,
     &          kperm,rank,eleproc,ndpn,.false.,dmuint,timerule,
     &          gmuint,stretch,g_ele,maxprocs)

              endif
            else

                call stf(nnds,nints,ndim,B,NN,det,node,maxele,
     &          maxnd,maxcon,maxdbc,nvpdofs,ele,dt,Cstiff,
     &          gcon,mpc,mpcnum,conA,nvpdbcs,nmudofs,nmudbcs,
     &          nudofs,nudbcs,muint,Jacob,Finv,Imperm,w0,
     &          gradv,pint,tau,PetAvpu,e0,dstrn,bsf,vint,strn,
     &          Hpsi,dJacob,dphi,phi,beta,a0,isu,isvp,ismu,
     &          PetR,PetdU,ksp,nodes,indices,gradp,dHdp,
     &          visc_eff,dHds,local_eles,p0,Cf,fluidsource,
     &          biot,phi0,freestrain,YM,poisson,KDtN,maxbdnd,
     &          DtNdofs,DtNidx,tc_split_on,muavg,open_thresh,
     &          kperm,rank,eleproc,ndpn,.false.,dmuint,timerule,
     &          gmuint,stretch,g_ele,maxprocs)

            endif

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!          Create the linear solver and set various options
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!  Set operators. Here the matrix that defines the linear system
!  also serves as the preconditioning matrix.

!            call KSPSetOperators(ksp,PetA,PetA,ierr)

!  Set runtime options, e.g.,
!      -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
!  These options will override those specified above as long as
!  KSPSetFromOptions() is called _after_ any other customization
!  routines.

            call KSPSetFromOptions(ksp,ierr)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!                      Solve the linear system
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

            call KSPSolve(ksp,PetR,PetdU,ierr)

!            call KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD,ierr)

 911        continue

!            call VecView(PetdU,PETSC_VIEWER_STDOUT_WORLD,ierr)

            call VecNorm(PetdU,NORM_2,sumsqu,ierr)

!            call VecGetValues(PetdU,nvpdofs+nudofs+nmudofs,
!     &                        indices,drhsT,ierr)

            call VecGetOwnershipRange(PetdU,low,high,ierr)
            call VecGetLocalSize(PetdU,proc_count,ierr)
            call VecGetArray(PetdU,get,i_get,ierr)

            do 306 i=1,proc_count
              drhsTproc(i)=get(i_get+i)
 306        continue

            call VecRestoreArray(PetdU,get,i_get,ierr)

!********** POPULATE drhsT ON EACH PROCESS

            call MPI_Allgather(proc_count,1,MPI_INTEGER,
     &                         proc_counts,1,MPI_INTEGER,
     &                         MPI_COMM_WORLD,ierr)

            call MPI_Allgather(low,1,MPI_INTEGER,
     &                         lows,1,MPI_INTEGER,
     &                         MPI_COMM_WORLD,ierr)

            call MPI_Gatherv(drhsTproc,proc_count,MPI_REAL8,
     &                       drhsT,proc_counts,lows,MPI_REAL8,
     &                       0,MPI_COMM_WORLD,ierr)

            call MPI_Bcast(drhsT,ndofs,MPI_REAL8,
     &                     0,MPI_COMM_WORLD,ierr)

            ! synchronization
!            call MPI_Barrier(MPI_COMM_WORLD,ierr)

!********** DONE

            do 550 i=1,nodes
              do 551 j=1,ndim+1
                if (gconvp(i,j).le.nvpdofs)
     &            drhs(gconvp(i,j))=drhsT(gcon(i,j+ndim))
 551          continue
 550        continue
            
            call calcvp(ndpn,ndim,maxnd,maxdbc,maxcon,nodes,nvpdofs,
     &                  gconvp,dbcval,drhs,v,p,nudbcs,dbcnd,dbcdof,
     &                  vpmpcnum,nvpdbcs,nconvp,vpmpc,vpmpcnd,
     &                  conAvp,conBvp,nmudbcs)
     
            call gradvp(nnds,nints,ndim,maxele,maxnd,local_eles,
     &                  ele,B,NN,gradv,vint,gradp,pint,v,p,g_ele,
     &                  maxprocs)

            do 554 i=1,nodes
              if (gconmu(i,1).le.nmudofs)
     &          drhs(gconmu(i,1))=drhsT(gcon(i,6))
 554        continue

            call calcmu(ndpn,ndim,maxnd,maxdbc,maxcon,nodes,nmudofs,
     &                  gconmu,dbcval,drhs,mu,nudbcs,dbcnd,dbcdof,
     &                  mumpcnum,nmudbcs,nconmu,mumpc,mumpcnd,
     &                  conAmu,conBmu,nvpdbcs)
                
            call gradm(nnds,nints,ndim,maxele,maxnd,local_eles,ele,B,
     &                 NN,muint,gradmu,mu,dt,muavg,dmuint,g_ele,
     &                 maxprocs)

            do 552 i=1,nodes
              do 553 j=1,ndim
                if (gconu(i,j).le.nudofs)
     &            drhs(gconu(i,j))=drhsT(gcon(i,j))
 553          continue
 552        continue

            call calcu(ndim,ndpn,maxnd,maxdbc,nodes,nudofs,nconu,
     &                 gconu,dbcval,drhs,u,umpc,nudbcs,dbcnd,dbcdof,
     &                 conAu,conBu,umpcnd,nmudbcs,nvpdbcs,umpcnum,
     &                 maxcon)

            call gradu(ndim,nnds,nints,maxele,maxnd,local_eles,ele,
     &               B,NN,strn,u,Jacob,Finv,dJacob,dstrn,e0,dt,uint,
     &               duint,FiJ,muint,phi,open_thresh,g_ele,maxprocs,
     &               maxJ)

            call nodalstress(ndim,maxnd,maxele,nints,nnds,nodes,neles,
     &                     ele,gradp,gradmu,strn,strs,ngradp,ngradmu,
     &                     nstrn,nstrs,ele_per_node,ele_at_node,
     &                     eleproc,maxprocs,l_ele,det)

            call normmu(nnds,nints,ndim,maxele,maxnd,local_eles,
     &                 ele,B,NN,ngradmu,gmuint,g_ele,maxprocs)

            do 611 i=1,nodes
              do 612 j=1,ndim
                if (dabs(u(i,j)).gt.maxU_u) maxU_u=dabs(u(i,j))
                if (dabs(v(i,j)).gt.maxU_v) maxU_v=dabs(v(i,j))
 612          continue
              if (dabs(p(i)).gt.maxU_v) maxU_v=dabs(p(i))
              if (dabs(mu(i)).gt.maxU_mu) maxU_mu=dabs(mu(i))
 611        continue

            if (maxU_u.lt.1.d-14) maxU_u=1.d0
            if (maxU_v.lt.1.d-14) maxU_v=1.d0
            if (maxU_p.lt.1.d-14) maxU_p=1.d0
            if (maxU_mu.lt.1.d-14) maxU_mu=1.d0

            do 613 i=1,nodes
              do 614 j=1,ndpn
                dof=gcon(i,j)
                if (j.le.ndim) then
                  drhs_normalize(dof)=(drhsT(dof))/maxU_u
                elseif (j.le.ndim+ndim+1) then
                  drhs_normalize(dof)=(drhsT(dof))/maxU_v
                elseif (j.eq.ndim+ndim+1) then
                  drhs_normalize(dof)=(drhsT(dof))/maxU_p
                else
                  drhs_normalize(dof)=(drhsT(dof))/maxU_mu
                endif
 614          continue
 613        continue

            sumsqu_n=0.d0
            do 621 i=1,ndofs
              sumsqu_n=sumsqu_n
     &                 +drhs_normalize(i)*drhs_normalize(i)
 621        continue
            sumsqu_n=dsqrt(sumsqu_n)

            if (rank.eq.0)
     &      call VecSetValues(PetdUnorm,nvpdofs+nudofs+nmudofs,
     &           indices,drhs_normalize,INSERT_VALUES,ierr)

            call VecAssemblyBegin(PetdUnorm,ierr)
            call VecAssemblyEnd(PetdUnorm,ierr)

!            call VecView(PetR,PETSC_VIEWER_STDOUT_WORLD,ierr)
!            call VecNorm(PetdUnorm,NORM_2,sumsqu_n,ierr)

            if (rank.eq.0) then
            print*,'|           sumsqu=',sumsqu,
     &             'relative=',sumsqu_n
            endif

            if (rank.eq.0) then
            do 601 i=1,nodes
              if (dabs(node(i,1)-pgauge(1)).lt.1.d-4) then
              if (dabs(node(i,2)-pgauge(2)).lt.1.d-4) then
                print*,'|        PRESSURE GAUGE: ',p(i)
              endif
              endif
 601        continue
            endif
           
!            if (betaOn.eq.1) then
!              if (ii.eq.1.and.sumsqu.lt.2.d0) then
!                betafactor=1.d-1
!              endif
!            endif

            if (isol.eq.solsteps-1) then
              conv1=.true.
              conv2=.true.
            elseif (ii.eq.1.and.isol.gt.2) then
              if (sumsqf_n.lt.1.d-2.or.sumsqf.lt.1.d-10) then
              if (sumsqu_n.lt.1.d-2.or.sumsqu.lt.1.d-10) then
                if (isu.eq.1) conv1=.true.
                if (ismu.eq.1) conv2=.true.
              endif
              endif
            endif


            if (sumsqf_n.lt.1.d-6.and.
     &          sumsqu_n.lt.1.d-6) then
              if (convall) then
                if (tt.ge.t) then
                  goto 901
                else
                  div_reset=.false.
                  goto 902
                endif
              endif
              goto 111
            elseif (sumsqf.lt.1.d-6.and.sumsqu.lt.1.d-6) then
              if (convall) then
                if (tt.ge.t) then
                  goto 901
                else
                  div_reset=.false.
                  goto 902
                endif
              endif
              goto 111
            elseif (ii.eq.maxsol) then
              if (rank.eq.0) print*,'|  MAXSOL HAS BEEN REACHED'
              
!              open(unit=51,file='results1/crack_vs_pressure'
!     &            ,access='append')
!              write(51,*),crackloc,papp,ibcs-10
!              close(unit=51)
              
              if (convall) then
                if (tt.ge.t) then
                  tt=tt-dt
                  dt=dt*0.5d0
                  t_resets=t_resets+1
                  div_reset=.true.
                  goto 903
                else
                  div_reset=.false.
                  goto 902
                endif
              endif
              goto 111
            endif

 912        continue
    
            if (sumsqf+sumsqu.gt.diverged) then
              if (rank.eq.0) print*,'|  SCHEME IS DIVERGENT'
              if (rank.eq.0) print*,'|'
              tt=tt-dt
              dt=dt*0.5d0
              t_resets=t_resets+1
              if (control.lt.csteps) then
              if (convall.and.t_resets.le.12) then
                div_reset=.true.
                goto 903
              endif
              endif
              goto 999
            endif

            if (1.eq.0.and.rank.eq.0) then             
            call paraview(ii,ou3,maxnd,maxele,ndim,nodes,neles,nnds,
     &                  node,ele,u,mu,v,p,ngradp,nstrn,nstrs,e0,
     &                  ngradmu,nints,NN,phi,w0,maxprocs)
            endif

 101        continue
 111        continue

!            call findcrack(ndim,maxnd,node,mu,xcrack,nodes)
!            if (rank.eq.0) print*,'|  ~Crack: ',xcrack
!            if (xcrack.gt.crackloc+10.d0) then
!              print*,'|       p_critical=',papp
!              print*,'|       error=',ibcs-10,'%'
!              goto 999
!            endif
            

            if (ibcs.eq.-1.and.rank.eq.0) then             
            call paraview(isol,ou2,maxnd,maxele,ndim,nodes,neles,nnds,
     &                  node,ele,u,mu,v,p,ngradp,nstrn,nstrs,e0,
     &                  ngradmu,nints,NN,phi,w0,maxprocs)
            endif

 100      continue
 902      continue

          if (control.gt.0.and.rank.eq.0) then
          call paraview(control,ou2,maxnd,maxele,ndim,nodes,neles,
     &                  nnds,node,ele,u,mu,v,p,ngradp,nstrn,nstrs,e0,
     &                  ngradmu,nints,NN,phi,w0,maxprocs)
          endif

          do 921 i=1,nodes
            u_temp(i,1)=u(i,1)
            u_temp(i,2)=u(i,2)
            v_temp(i,1)=v(i,1)
            v_temp(i,2)=v(i,2)
            mu_temp(i)=mu(i)
            p_temp(i)=p(i)        
 921      continue

 903    continue
      
 900    continue
 901    continue

        if (rank.eq.0)
     &  print*,'======================================='

        call stresses(nints,ndim,maxele,local_eles,strn,gradv,
     &                    gradmu,pint,gradp,vint,strs,tau,bsf,bp,bv,
     &                    psimu,ksi,phi0,Jacob,Finv,Cstiff,w0,
     &                    Imperm,dJacob,dphi,dstrn,Hpsi,e0,duint,
     &                    1,a0,phi,visc_eff,muint,g_ele,bphiJ,
     &                    dHds,Hpsi0,biot,dt,freestrain,YM,poisson,
     &                    tc_split_on,muavg,kperm,Wp,uint,Wmu,Ws,
     &                    dmuint,fluidsource,rank,open_thresh,
     &                    gmuint,scale_visc,stretch,maxprocs,Cf,
     &                    p0,maxJ,dHdp)

        Wp_total=0.d0
        Wmu_total=0.d0
        Ws_total=0.d0
        vflux=0.d0
        val=0.d0
        val2=0.d0
!        print*,' xx yy stretch(i,j) gmuint(i,j,1) gmuint(i,j,2)', 
!     &         ' gradmu(i,j,1) gradmu(i,j,2) Jacob(i,j,1)', 
!     &         ' FiJ(i,j,1,1) FiJ(i,j,1,2) FiJ(i,j,2,1) FiJ(i,j,2,2)'
        mass_bal=0.d0
        do 140 i=1,local_eles
          area_sum=0.d0
          do 142 j=1,nnds
            BtB(j,1,1)=0.d0
            BtB(j,1,2)=0.d0
            BtB(j,2,1)=0.d0
            BtB(j,2,2)=0.d0
 142      continue
!          val=0.d0
!          val2=0.d0
          do 141 j=1,nints
            xx=0.d0
            yy=0.d0
            do 773 k=1,nnds
              xx=xx+NN(k,j)*node(ele(i,k),1)
              yy=yy+NN(k,j)*node(ele(i,k),2)
              BtB(k,1,1)=BtB(k,1,1)+B(i,k,j,1)*B(i,k,j,1)*det(i,j)
              BtB(k,1,2)=BtB(k,1,2)+B(i,k,j,1)*B(i,k,j,2)*det(i,j)
              BtB(k,2,1)=BtB(k,2,1)+B(i,k,j,2)*B(i,k,j,1)*det(i,j)
              BtB(k,2,2)=BtB(k,2,2)+B(i,k,j,2)*B(i,k,j,2)*det(i,j)
 773        continue
            if (i.eq.13203) print*,i,xx,yy,B(i,3,j,1)
            if (i.eq.13204) print*,i,xx,yy,B(i,1,j,1)
            area_sum=area_sum+det(i,j)
            if (dabs(yy).lt.0.375d0) then
!              print*,i,xx,yy,Jacob(i,j,1),vint(i,j,1)*w0
!              print*,'Finv:',Finv(i,j,1,1,1),Finv(i,j,1,2,1)
!     &               ,Finv(i,j,2,1,1),Finv(i,j,2,2,1)
              if (j.eq.nints) then
!                if (i.eq.13203.or.i.eq.13204)
!     &            print*,area_sum,(BtB(k,1,1),k=1,nnds)
              endif
            endif
            Wp_total=Wp_total+Wp(i,j)*det(i,j)
            Wmu_total=Wmu_total+Wmu(i,j)*det(i,j)
            Ws_total=Ws_total+Ws(i,j)*det(i,j)
            if (yy.gt.0.d0.or.1.eq.0)
     &        vflux=vflux+vint(i,j,1)*det(i,j)*Jacob(i,j,1)
            if (yy.lt.0.d0)
     &        val=val+vint(i,j,1)*det(i,j)*Jacob(i,j,1)
            mass_bal=mass_bal+det(i,j)*(bv(i,j)+bphiJ(i,j))
            if (dabs(yy).lt.0.25d0) val=val+det(i,j)*bv(i,j)
            if (dabs(yy).lt.0.25d0) val2=val2+det(i,j)*bphiJ(i,j)
 141      continue
!          if (dabs(val).gt.1.d-15) print*,node(ele(i,4),1),val,val2
!          print*,Jacob(i,5,1),maxJ(i)
 140    continue
!        print*,'MASS: ',mass_bal

        if (rank.eq.0) print*,'|  Wp = ',Wp_total
        if (rank.eq.0) print*,'|  Wmu = ',Wmu_total
        if (rank.eq.0) print*,'|  2 Ws = ',2.d0*Ws_total
        if (rank.eq.0) print*,'|  Q = ',vflux/(2.d0)*w0

!        call LsqGrads(nints,ndim,ndpn,nnds,maxnd,maxdbc,maxcon,
!     &              nodes,maxele,nvpdofs,nudofs,nmudofs,nvpdbcs,
!     &              nudbcs,nmudbcs,ndofs,ndbcs,indices,neles,
!     &              ele,gcon,gconvp,gconu,gconmu,dbcval,dbcnd,
!     &              dbcdof,mpcnum,mpc,mpcnd,conA,conB,vpmpcnum,
!     &              nconvp,vpmpc,vpmpcnd,conAvp,conBvp,umpcnum,
!     &              nconu,umpc,umpcnd,conAu,conBu,mumpcnum,nconmu,
!     &              mumpc,mumpcnd,conAmu,conBmu,det,B,NN,
!     &              strn,gradv,gradp,gradmu,nstrn,ngradv,ngradp,
!     &              ngradmu,kspG,LsqF,LsqU)
!
!        call LsqGrads(nints,ndim,ndpn,nnds,maxnd,maxdbc,maxcon,
!     &              nodes,maxele,nvpdofs,nudofs,nmudofs,nvpdbcs,
!     &              nudbcs,nmudbcs,ndofs,ndbcs,indices,neles,
!     &              ele,gcon,gconvp,gconu,gconmu,dbcval,dbcnd,
!     &              dbcdof,mpcnum,mpc,mpcnd,conA,conB,vpmpcnum,
!     &              nconvp,vpmpc,vpmpcnd,conAvp,conBvp,umpcnum,
!     &              nconu,umpc,umpcnd,conAu,conBu,mumpcnum,nconmu,
!     &              mumpc,mumpcnd,conAmu,conBmu,det,B,NN,
!     &              strs,gradv,gradp,gradmu,nstrs,ngradv,ngradp,
!     &              ngradmu,kspG,LsqF,LsqU)

        call nodalstress(ndim,maxnd,maxele,nints,nnds,nodes,neles,
     &                     ele,gradp,gradmu,strn,strs,ngradp,ngradmu,
     &                     nstrn,nstrs,ele_per_node,ele_at_node,
     &                     eleproc,maxprocs,l_ele,det)

        call normmu(nnds,nints,ndim,maxele,maxnd,local_eles,
     &                 ele,B,NN,ngradmu,gmuint,g_ele,maxprocs)

        if (rank.eq.0.and.mod(ibcs,1).eq.0) then

        call paraview(ibcs,ou1,maxnd,maxele,ndim,nodes,neles,nnds,
     &                  node,ele,u,mu,v,p,ngradp,nstrn,nstrs,
     &                  e0,ngradmu,nints,NN,phi,w0,maxprocs)
        endif

        do 920 i=1,nodes
          u_temp(i,1)=u(i,1)
          u_temp(i,2)=u(i,2)
          v_temp(i,1)=v(i,1)
          v_temp(i,2)=v(i,2)
          mu_temp(i)=mu(i)
          p_temp(i)=p(i)        
 920    continue

        call MPI_Barrier(MPI_COMM_WORLD,ierr)
        
       
!        if (rank.eq.0) write(50,*),t,e0*u(5002,1),p(5003)
!        if (rank.eq.0) write(50,*),t,p(39284)

!        call findcrack(ndim,maxnd,node,mu,xcrack,nodes)
!        if (rank.eq.0) then
!          print*,'|  ~Crack: ',xcrack
!          open(unit=50,file='results1/crack_vs_time',access='append')
!          write(50,*),t,xcrack,p(1),(u(3,1)-u(2,1))*e0
!          close(unit=50)
!        endif

        if (rank.eq.0) write(50,*),uapp,vflux/(10.d0)*w0+0.d0*kappa
     &                            ,val/10.d0*w0

 777  continue

 999  continue

!      close(unit=50)

 910  format(A)

      if (rank.eq.0) then
      print*
      print*,'********************************'
      print*,'*********  COMPLETED  **********'
      print*,'********************************'
      print*
      endif

      call MPI_Barrier(MPI_COMM_WORLD,ierr)

      call VecDestroy(PetR,ierr)
      call VecDestroy(PetRs,ierr)
      call VecDestroy(PetdU,ierr)
      call VecDestroy(PetRnorm,ierr)
      call VecDestroy(PetdUnorm,ierr)
!      call VecDestroy(LsqF,ierr)
!      call VecDestroy(LsqU,ierr)
      call MatDestroy(PetA,ierr)
      call MatDestroy(PetAvpu,ierr)
      call MatDestroy(PetAmu,ierr)
!      call MatDestroy(PetG,ierr)
      call KSPDestroy(ksp,ierr)
!      call KSPDestroy(kspG,ierr)


      call PetscFinalize(ierr)
      CHKERRQ(ierr)

      end

***************   Starting line of subroutines   *********************

      subroutine readinp(ndim,nnds,ndpn,maxnd,maxele,maxfc,maxdbc,
     &             nodes,neles,nudbcs,nfbcs,nudofs,gconu,ele,fbcnd,
     &             fbcdof,dbcnd,dbcdof,node,fbcval,dbcval,umpc,conAu,
     &             conBu,umpcnum,nconu,umpcnd,nmudofs,gconmu,
     &             nmudbcs,maxcon,nconmu,conAmu,conBmu,mumpcnum,
     &             mumpc,mumpcnd,nvpdofs,gconvp,nvpdbcs,nconvp,
     &             conAvp,conBvp,vpmpcnum,vpmpc,vpmpcnd,conAvpu,
     &             conBvpu,vpumpcnum,vpumpc,vpumpcnd,
     &             conA,conB,mpcnum,mpc,mpcnd,gcon,freestrain,
     &             phi0,biot,YM,poisson,fluidsource,Imperm,
     &             moduli_on,rank,maxprocs)
      implicit none
      real*8 node(maxnd,ndim),dv,dbcval(maxdbc),fbcval(maxfc)
      real*8 conAu(maxdbc,maxcon),conBu(maxdbc,maxcon)
      real*8 conAmu(maxdbc,maxcon),conBmu(maxdbc,maxcon)
      real*8 conAvp(maxdbc,maxcon),conBvp(maxdbc,maxcon)
      real*8 conAvpu(maxdbc,maxcon),conBvpu(maxdbc,maxcon)
      real*8 conA(maxdbc,maxcon),conB(maxdbc,maxcon)
      real*8 YM(maxele),biot(maxele)
      real*8 phi0(maxele),poisson(maxele)
      real*8 freestrain(maxele),fluidsource(maxele)
      real*8 Imperm(maxele,ndim,ndim)
      integer fbcnd(maxfc),fbcdof(maxfc),nudofs,umpc(maxdbc,maxcon,2)
      integer dbcnd(maxdbc),dbcdof(maxdbc)
      integer gconmu(maxnd,ndpn),nmudofs,nmudbcs,umpcnum(maxdbc)
      integer gconvp(maxnd,ndpn),nvpdofs,nvpdbcs
      integer i,j,k,nodes,neles,nudbcs,nfbcs
      integer gconu(maxnd,ndpn),numvpu,numtot
      integer maxnd,maxele,maxdbc,maxfc,maxcon
      integer n,d,inds,iele,imat,idbc,ifbc,ndim,nnds,ndpn,imod
      integer icontr,nconu,num,umpcnd(maxdbc,2),nconmu
      integer mumpc(maxdbc,maxcon,2),mumpcnd(maxdbc,2)
      integer nconvp,vpmpcnum(maxdbc),vpmpc(maxdbc,maxcon,2)
      integer vpmpcnd(maxdbc,2),mumpcnum(maxdbc)
      integer vpumpcnum(maxdbc),vpumpc(maxdbc,maxcon,2)
      integer vpumpcnd(maxdbc,2),ele(maxele*maxprocs,nnds)
      integer mpcnum(maxdbc),mpc(maxdbc,maxcon,2)
      integer mpcnd(maxdbc,2),gcon(maxnd,ndpn),rank,maxprocs
      integer igcon,igconu,igconmu,igconvp
      logical moduli_on


      inds=15
      iele=16
      imat=17
      idbc=18
      ifbc=19
      imod=20
      icontr=31
      igcon=21
      igconu=22
      igconmu=23
      igconvp=24

***************************************************************
      if (rank.eq.0) print*,'...Reading nodes'
      open(unit=inds,file='input/nodes')
      read(inds,*)nodes
      if (nodes.gt.maxnd) then 
         print*,'too many nodes'
         print*,'nodes=',nodes,'maxnd=',maxnd
         stop
      endif
      do 250 i=1,nodes
        read(inds,*)k,(node(i,j),j=1,ndim)
        node(i,1)=node(i,1)*1.d0
        node(i,2)=node(i,2)*1.d0
 250  continue
      close(unit=inds)
      nudofs=nodes*ndim
      nmudofs=nodes
      nvpdofs=nodes*(ndim+1)

      if (rank.eq.0) print*,'...Reading elements'
      open(unit=iele,file='input/elements')
      read(iele,*)neles
      if (neles.gt.maxele*maxprocs) then 
         print*,'too many elements'
         print*,'neles',neles,'maxele',maxele*maxprocs
         stop
      endif
      do 300 i=1,neles
        read(iele,*)k,(ele(i,j),j=1,nnds)
 300  continue
      close(unit=iele)

      if (rank.eq.0) print*,'...Constructing Connectivity'
      do 390 i=1,nodes
        gconmu(i,1)=i
        do 391 j=1,ndim
          gconu(i,j)=(i-1)*ndim+j
          gconvp(i,j)=(i-1)*(ndim+1)+j
 391    continue
        gconvp(i,ndim+1)=i*(ndim+1)
        do 393 j=1,ndpn
          gcon(i,j)=(i-1)*(ndpn)+j
 393    continue
 390  continue

!      do 390 i=1,nodes
!        do 391 j=1,ndim
!          gconu(i,j)=(i-1)*ndim+j
!          gconvp(i,j)=(i-1)*ndim+j
!          gcon(i,j)=(i-1)*ndim+j
!          gcon(i,j+ndim)=ndim*nodes+(i-1)*ndim+j
! 391    continue
!        gconmu(i,1)=i
!        gconvp(i,ndim+1)=ndim*(nodes)+i
!        gcon(i,2*ndim+1)=2*ndim*nodes+i
!        gcon(i,ndpn)=(ndpn-1)*nodes+i
! 390  continue

      nmudbcs=0
      nvpdbcs=0
      open(unit=idbc,file='input/dbcs')
      read(idbc,*)nudbcs,nconu,nconvp,nconmu
      if (nudbcs.gt.maxdbc) then
        print*,'too many dbcs'
        stop
      endif
      if (nconu+nconmu+nconvp.gt.maxdbc) then
        print*,'too many multi-point constraints'
        stop
      endif

      if (rank.eq.0) print*,'...Reading dbcs'   
      do 400 i=1,nudbcs
        read(idbc,*)n,d,dv
        dbcnd(i)=n
        dbcdof(i)=d
        dbcval(i)=dv
        if (d.le.ndim) then
!          num=gconu(n,d)
!          numtot=gcon(n,d)
!          do 410 j=1,nodes
!            do 420 k=1,ndim
!              if(gconu(j,k).gt.num) gconu(j,k)=gconu(j,k)-1
! 420        continue
!            do 460 k=1,ndpn
!              if(gcon(j,k).gt.numtot) gcon(j,k)=gcon(j,k)-1
! 460        continue
! 410      continue
!          gconu(n,d)=nodes*ndim
!          gcon(n,d)=nodes*(ndpn)
          nudofs=nudofs-1
        else
          if (d.eq.ndpn) then
            nmudbcs=nmudbcs+1
!            num=gconmu(n,1)
!            numtot=gcon(n,d)
!            do 411 j=1,nodes
!              if(gconmu(j,1).gt.num) gconmu(j,1)=gconmu(j,1)-1
!              do 461 k=1,ndpn
!                if(gcon(j,k).gt.numtot) gcon(j,k)=gcon(j,k)-1
! 461          continue
! 411        continue
!            gconmu(n,1)=nodes
!            gcon(n,d)=nodes*(ndpn)
            nmudofs=nmudofs-1
          else
            nvpdbcs=nvpdbcs+1
!            num=gconvp(n,d-ndim)
!            numtot=gcon(n,d)
!            do 412 j=1,nodes
!              do 413 k=1,ndim+1
!                if(gconvp(j,k).gt.num) gconvp(j,k)=gconvp(j,k)-1
! 413          continue
!              do 462 k=1,ndpn
!                if(gcon(j,k).gt.numtot) gcon(j,k)=gcon(j,k)-1
! 462          continue
! 412        continue
!            gconvp(n,d-ndim)=nodes*(ndim+1)
!            gcon(n,d)=nodes*(ndpn)
            nvpdofs=nvpdofs-1
          endif
        endif
 400  continue
      nudbcs=nudbcs-nmudbcs-nvpdbcs
      
      if (rank.eq.0) print*,'...Reading umpcs'
      do 425 i=1,nconu
        read(idbc,*)umpcnum(i),umpcnd(i,1),umpcnd(i,2)
        if (umpcnum(i).gt.maxcon) then
          print*,'too much connectivity'
          stop
        endif
        do 426 j=1,umpcnum(i)
          read(idbc,*)umpc(i,j,1),umpc(i,j,2),conAu(i,j),conBu(i,j)
 426    continue
!        num=gconu(umpcnd(i,1),umpcnd(i,2))
!        numtot=gcon(umpcnd(i,1),umpcnd(i,2))
!        do 427 j=1,nodes
!          do 428 k=1,ndim
!            if(gconu(j,k).gt.num) gconu(j,k)=gconu(j,k)-1
! 428      continue
!          do 463 k=1,ndpn
!            if(gcon(j,k).gt.numtot) gcon(j,k)=gcon(j,k)-1
! 463      continue
! 427    continue
!        gconu(umpcnd(i,1),umpcnd(i,2))=nodes*ndim
!        gcon(umpcnd(i,1),umpcnd(i,2))=nodes*(ndpn)
        nudofs=nudofs-1
 425  continue

      if (rank.eq.0) print*,'...Reading vpmpcs'
      do 445 i=1,nconvp
        read(idbc,*)vpmpcnum(i),vpmpcnd(i,1),vpmpcnd(i,2)
        if (vpmpcnum(i).gt.maxcon) then
          print*,'too much connectivity'
          stop
        endif
        do 446 j=1,vpmpcnum(i)
          read(idbc,*)vpmpc(i,j,1),vpmpc(i,j,2),conAvp(i,j),conBvp(i,j)
 446    continue
!        num=gconvp(vpmpcnd(i,1),vpmpcnd(i,2)-ndim)
!        numtot=gcon(vpmpcnd(i,1),vpmpcnd(i,2))
!        do 449 j=1,nodes
!          do 450 k=1,ndim+1
!            if(gconvp(j,k).gt.num) gconvp(j,k)=gconvp(j,k)-1
! 450      continue
!          do 464 k=1,ndpn
!            if(gcon(j,k).gt.numtot) gcon(j,k)=gcon(j,k)-1
! 464      continue
! 449    continue
!        gconvp(vpmpcnd(i,1),vpmpcnd(i,2)-ndim)=nodes*(ndim+1)
!        gcon(vpmpcnd(i,1),vpmpcnd(i,2))=nodes*(ndpn)
        nvpdofs=nvpdofs-1
 445  continue

      if (rank.eq.0) print*,'...Reading mumpcs'
      do 435 i=1,nconmu
        read(idbc,*)mumpcnum(i),mumpcnd(i,1),mumpcnd(i,2)
        if (mumpcnum(i).gt.maxcon) then
          print*,'too much connectivity'
          stop
        endif
        do 436 j=1,mumpcnum(i)
         read(idbc,*)mumpc(i,j,1),mumpc(i,j,2),
     &               conAmu(i,j),conBmu(i,j)
 436    continue
!        num=gconmu(mumpcnd(i,1),1)
!        numtot=gcon(mumpcnd(i,1),mumpcnd(i,2))
!        do 439 j=1,nodes
!          if((gconmu(j,1)).gt.num) gconmu(j,1)=gconmu(j,1)-1
!          do 465 k=1,ndpn
!            if(gcon(j,k).gt.numtot) gcon(j,k)=gcon(j,k)-1
! 465      continue
! 439    continue
!        gconmu(mumpcnd(i,1),1)=nodes
!        gcon(mumpcnd(i,1),mumpcnd(i,2))=nodes*(ndpn)
        nmudofs=nmudofs-1
 435  continue

      open(unit=igcon,file='input/gcon')
      do 520 i=1,nodes
        read(igcon,*)(gcon(i,j),j=1,ndpn)
 520  continue
      close(unit=igcon)

      open(unit=igconu,file='input/gconu')
      do 521 i=1,nodes
        read(igconu,*)(gconu(i,j),j=1,ndim)
 521  continue
      close(unit=igconu)

      open(unit=igconmu,file='input/gconmu')
      do 522 i=1,nodes
        read(igconmu,*)gconmu(i,1)
 522  continue
      close(unit=igconmu)

      open(unit=igconvp,file='input/gconvp')
      do 523 i=1,nodes
        read(igconvp,*)(gconvp(i,j),j=1,ndim+1)
 523  continue
      close(unit=igconvp)

      do 510 i=1,nconu
        vpumpcnum(i)=umpcnum(i)
        vpumpcnd(i,1)=umpcnd(i,1)
        vpumpcnd(i,2)=umpcnd(i,2)
        mpcnum(i)=umpcnum(i)
        mpcnd(i,1)=umpcnd(i,1)
        mpcnd(i,2)=umpcnd(i,2)
        do 511 j=1,umpcnum(i)
          vpumpc(i,j,1)=umpc(i,j,1)
          vpumpc(i,j,2)=umpc(i,j,2)
          conAvpu(i,j)=conAu(i,j)
          conBvpu(i,j)=conBu(i,j)
          mpc(i,j,1)=umpc(i,j,1)
          mpc(i,j,2)=umpc(i,j,2)
          conA(i,j)=conAu(i,j)
          conB(i,j)=conBu(i,j)
 511    continue
 510  continue
      do 512 i=1,nconvp
        vpumpcnum(i+nconu)=vpmpcnum(i)
        vpumpcnd(i+nconu,1)=vpmpcnd(i,1)
        vpumpcnd(i+nconu,2)=vpmpcnd(i,2)
        mpcnum(i+nconu)=vpmpcnum(i)
        mpcnd(i+nconu,1)=vpmpcnd(i,1)
        mpcnd(i+nconu,2)=vpmpcnd(i,2)
        do 513 j=1,vpmpcnum(i)
          vpumpc(i+nconu,j,1)=vpmpc(i,j,1)
          vpumpc(i+nconu,j,2)=vpmpc(i,j,2)
          conAvpu(i+nconu,j)=conAvp(i,j)
          conBvpu(i+nconu,j)=conBvp(i,j)
          mpc(i+nconu,j,1)=vpmpc(i,j,1)
          mpc(i+nconu,j,2)=vpmpc(i,j,2)
          conA(i+nconu,j)=conAvp(i,j)
          conB(i+nconu,j)=conBvp(i,j)
 513    continue
 512  continue
      do 514 i=1,nconmu
        mpcnum(i+nconu+nconvp)=mumpcnum(i)
        mpcnd(i+nconu+nconvp,1)=mumpcnd(i,1)
        mpcnd(i+nconu+nconvp,2)=mumpcnd(i,2)
        do 515 j=1,mumpcnum(i)
          mpc(i+nconu+nconvp,j,1)=mumpc(i,j,1)
          mpc(i+nconu+nconvp,j,2)=mumpc(i,j,2)
          conA(i+nconu+nconvp,j)=conAmu(i,j)
          conB(i+nconu+nconvp,j)=conBmu(i,j)
 515    continue
 514  continue

      if (rank.eq.0) print*,'...Reading fbcs'
      open(unit=ifbc,file='input/fbcs')
      read(ifbc,*)nfbcs
      if (nfbcs.gt.maxfc) then
        print*,'|  *ERROR: too many force bcs'
        stop
      endif
      do 500 i=1,nfbcs
        read(ifbc,*)n,d,dv
        if (dv.eq.ndpn) then
          print*,'Force BCs not allowed on mu'
        endif
        fbcnd(i)=n
        fbcdof(i)=d
        fbcval(i)=dv
 500  continue
      close(unit=ifbc)


!      if (moduli_on) then
!       if (rank.eq.0) print*,'...Reading moduli'
!       open(unit=imod,file='input/moduli')
!       do 501 i=1,neles
!         if (ndim.eq.2) then
!           read(imod,*)k,YM(i),biot(i),phi0(i),poisson(i),
!     &               freestrain(i),fluidsource(i),
!     &               Imperm(i,1,1),Imperm(i,1,2),
!     &               Imperm(i,2,1),Imperm(i,2,2)
!         else
!           print*,
!     &      ' ERROR: Readinp does not support number of dimensions!'
!           stop
!         endif
! 501   continue
!       close(unit=imod)
!      else
!       do 502 i=1,neles
!         YM(i)=1.d0
!         biot(i)=0.d0
!         phi0(i)=0.d0
!         poisson(i)=0.25d0*1.d0
!         freestrain(i)=0.d0
!         fluidsource(i)=0.d0
!         do 503 j=1,ndim
!           do 504 k=1,ndim
!             Imperm(i,j,k)=0.d0
!             if (j.eq.k) Imperm(i,j,k)=1.d0
! 504       continue
! 503     continue
! 502   continue
!      endif



      return
      end

**********************************************************************
*************************************************************
*     This subroutine will construct the bele(8,4,2) for    *
*     each element and store them in B(maxele,8,4,2).       *
*************************************************************

      subroutine calcB(ndim,nnds,nints,maxnd,maxele,local_eles,
     &                 ele,node,det,B,NN,g_ele,maxprocs)
      implicit none
      real*8 node(maxnd,ndim),B(maxele,nnds,nints,ndim)
      real*8 det(maxele,nints),detp,ted,eleA
      real*8 rs(4),r,s,p(nnds,ndim),aj(ndim,ndim),aji(ndim,ndim)
      real*8 NN(nnds,nints),N9,ww(nints),w1,w2,rT(4),sT(4)
      integer ndim,nnds,nints,nc
      integer i,j,k,m,ele(maxele*maxprocs,nnds),maxnd,maxele
      integer k1,k2,k123,local_eles,maxprocs,g_ele(maxele)
      logical iserr

      
      iserr=.false.
      nc=idint(dsqrt(dble(nints)))
      if (nnds.eq.6.or.nnds.eq.3) then
        if (nints.eq.4) then
          rT(1)=1.d0/3.d0
          rT(2)=1.d0/5.d0
          rT(3)=1.d0/5.d0
          rT(4)=3.d0/5.d0
          sT(1)=1.d0/3.d0
          sT(2)=1.d0/5.d0
          sT(3)=3.d0/5.d0
          sT(4)=1.d0/5.d0
        else
          print*,'Number of integration stations is not supported'
          stop
        endif
      else
        if (nints.eq.9) then
          rs(1)=-dsqrt(0.6d0)
          rs(2)=0.d0
          rs(3)=-rs(1)
        elseif (nints.eq.16) then
          rs(1)=-0.861136311594053d0
          rs(2)=-0.339981043584856d0
          rs(3)=-rs(2)
          rs(4)=-rs(1)
        elseif (nints.eq.4) then
          rs(1)=-1.d0/dsqrt(3.d0)
          rs(2)=-rs(1)
        else
          print*,'Number of integration stations is not supported'
          stop
        endif
      endif

      if (nnds.eq.6.or.nnds.eq.3) then
        if (nints.eq.4) then
          ww(1)=(-27.d0/48.d0)/2.d0
          ww(2)=(25.d0/48.d0)/2.d0
          ww(3)=ww(2)
          ww(4)=ww(2)
        endif
      else
        if (nints.eq.9) then
          ww(1)=25.d0/81.d0
          ww(2)=40.d0/81.d0
          ww(3)=ww(1)
          ww(4)=ww(2)
          ww(5)=64.d0/81.d0
          ww(6)=ww(2)
          ww(7)=ww(1)
          ww(8)=ww(2)
          ww(9)=ww(1)
        elseif (nints.eq.16) then
          w1=0.347854845137454d0
          w2=0.652145154862546d0
          ww(1)=w1*w1
          ww(2)=w1*w2
          ww(3)=w1*w2
          ww(4)=w1*w1
          ww(5)=w1*w2
          ww(6)=w2*w2
          ww(7)=w2*w2
          ww(8)=w2*w1
          ww(9)=w2*w1
          ww(10)=w2*w2
          ww(11)=w2*w2
          ww(12)=w1*w2
          ww(13)=w1*w1
          ww(14)=w1*w2
          ww(15)=w1*w2
          ww(16)=w1*w1
       elseif (nints.eq.4) then
          ww(1)=1.d0
          ww(2)=1.d0
          ww(3)=1.d0
          ww(4)=1.d0
        endif
      endif
      
      if ((nnds.ne.8).and.(nnds.ne.9).and.(nnds.ne.4)
     &    .and.(nnds.ne.6)) then
        print*,'Number of nodes is not supported'
        stop
      endif
     
      if (nnds.eq.6.or.nnds.eq.3) then
        do 11 k1=1,nints
          r=rT(k1)
          s=sT(k1)
          k123=k1
          if (nnds.eq.6) then
            NN(1,k123)=r*(2.d0*r-1.d0)
            NN(2,k123)=s*(2.d0*s-1.d0)
            NN(3,k123)=(r+s-1.d0)*(2.d0*r+2.d0*s-1.d0)
            NN(4,k123)=4.d0*r*s
            NN(5,k123)=4.d0*r*(1.d0-r-s)
            NN(6,k123)=4.d0*s*(1.d0-r-s)
          endif
 11     continue
      else
        do 10 k2=1,nc
         s=rs(k2)
         do 20 k1=1,nc
          r=rs(k1)
             
          k123=nc*k2-nc+k1
              
          if (nnds.eq.8) then
            NN(1,k123)=-(1.d0-r)*(1.d0-s)*(1.d0+r+s)/4.d0
            NN(2,k123)=-(1.d0+r)*(1.d0-s)*(1.d0-r+s)/4.d0
            NN(3,k123)=-(1.d0-r)*(1.d0+s)*(1.d0+r-s)/4.d0
            NN(4,k123)=-(1.d0+r)*(1.d0+s)*(1.d0-r-s)/4.d0
            NN(5,k123)=(1.d0-r*r)*(1.d0-s)/2.d0
            NN(6,k123)=(1.d0-r)*(1.d0-s*s)/2.d0
            NN(7,k123)=(1.d0+r)*(1.d0-s*s)/2.d0
            NN(8,k123)=(1.d0-r*r)*(1.d0+s)/2.d0
          elseif (nnds.eq.4) then
            NN(1,k123)=(1.d0-r)*(1.d0-s)/4.d0
            NN(2,k123)=(1.d0+r)*(1.d0-s)/4.d0
            NN(3,k123)=(1.d0-r)*(1.d0+s)/4.d0
            NN(4,k123)=(1.d0+r)*(1.d0+s)/4.d0
          elseif (nnds.eq.9) then
            N9=(1.d0-s*s)*(1.d0-r*r)
            NN(1,k123)=r*s*(r-1.d0)*(s-1.d0)/4.d0
            NN(2,k123)=r*s*(r+1.d0)*(s-1.d0)/4.d0
            NN(3,k123)=r*s*(r-1.d0)*(s+1.d0)/4.d0
            NN(4,k123)=r*s*(r+1.d0)*(s+1.d0)/4.d0
            NN(5,k123)=s*(1.d0-r*r)*(s-1.d0)/2.d0
            NN(6,k123)=r*(1.d0-s*s)*(r-1.d0)/2.d0
            NN(7,k123)=r*(1.d0-s*s)*(r+1.d0)/2.d0
            NN(8,k123)=s*(1.d0-r*r)*(s+1.d0)/2.d0
            NN(9,k123)=N9
          elseif (nnds.eq.6) then
            NN(1,k123)=-(1.d0-r)*(1.d0-s)*(1.d0+r+s)/4.d0
            NN(2,k123)=-(1.d0+r)*(1.d0-s)*(1.d0-r+s)/4.d0
            NN(3,k123)=-(1.d0-r)*(1.d0+s)*(1.d0+r-s)/4.d0
     &                 -(1.d0+r)*(1.d0+s)*(1.d0-r-s)/4.d0
     &                 +(1.d0-r*r)*(1.d0+s)/2.d0
            NN(4,k123)=(1.d0-r*r)*(1.d0-s)/2.d0
            NN(5,k123)=(1.d0-r)*(1.d0-s*s)/2.d0
            NN(6,k123)=(1.d0+r)*(1.d0-s*s)/2.d0
          endif
          
 20      continue
 10     continue
      endif
      
      do 100 m=1,local_eles

        if (nnds.eq.6.or.nnds.eq.3) then

          eleA=0.d0
          do 311 k1=1,nints
            r=rT(k1)
            s=sT(k1)
            k123=k1

            if (nnds.eq.6) then
              p(1,1)=4.d0*r-1.d0
              p(2,1)=0.d0
              p(3,1)=4.d0*(r+s)-3.d0
              p(4,1)=4.d0*s
              p(5,1)=4.d0*(1.d0-2.d0*r-s)
              p(6,1)=-4.d0*s
              p(1,2)=0.d0
              p(2,2)=4.d0*s-1.d0
              p(3,2)=4.d0*(r+s)-3.d0
              p(4,2)=4.d0*r
              p(5,2)=-4.d0*r
              p(6,2)=4.d0*(1.d0-r-2.d0*s)
            endif
 
            do 403 i=1,ndim
              do 411 j=1,ndim
                aj(i,j)=0.d0
                do 421 k=1,nnds
                  aj(i,j)=aj(i,j)+p(k,i)*node(ele(g_ele(m),k),j)
 421            continue
 411          continue
 403        continue
 
            if (ndim.ne.2) then
              print*,'Number of dimensions is not supported'
            endif
            detp=aj(1,1)*aj(2,2)-aj(1,2)*aj(2,1)
            
            ted=1.d0/detp
            
            aji(1,1)=ted*aj(2,2)
            aji(1,2)=-ted*aj(1,2)
            aji(2,1)=-ted*aj(2,1)
            aji(2,2)=ted*aj(1,1)
            
            do 501 i=1,ndim
              do 511 j=1,nnds
                B(m,j,k123,i)=0.d0
                do 521 k=1,ndim
                  B(m,j,k123,i)=B(m,j,k123,i)+aji(i,k)*p(j,k)
 521            continue
 511          continue
 501        continue
            
            det(m,k123)=detp*ww(k123)

            eleA=eleA+det(m,k123)            

 311      continue

          if (eleA.lt.0.d0) then
            print*,'ERROR: Negative Area',
     &             ' Element: ',g_ele(m),eleA
            iserr=.true.
          endif


        else

          do 310 k2=1,nc
            s=rs(k2)
            do 320 k1=1,nc
              r=rs(k1)             
              k123=nc*k2-nc+k1
 
              if (nnds.eq.8) then
                p(1,1)=-((-1.d0+s)*(2.d0*r+s))/4.d0
                p(2,1)=-((2.d0*r-s)*(-1.d0+s))/4.d0
                p(3,1)=((2.d0*r-s)*(1.d0+s))/4.d0
                p(4,1)=((1.d0+s)*(2.d0*r+s))/4.d0
                p(5,1)=r*(-1.d0+s)
                p(6,1)=(-1.d0+s**2.d0)/2.d0
                p(7,1)=(1.d0-s**2.d0)/2.d0
                p(8,1)=-(r*(1.d0+s))
                p(1,2)=-((-1.d0+r)*(r+2.d0*s))/4.d0
                p(2,2)=-((1.d0+r)*(r-2.d0*s))/4.d0
                p(3,2)=((-1.d0+r)*(r-2.d0*s))/4.d0
                p(4,2)=((1.d0+r)*(r+2.d0*s))/4.d0
                p(5,2)=(-1.d0+r**2.d0)/2.d0
                p(6,2)=(-1.d0+r)*s
                p(7,2)=-((1.d0+r)*s)
                p(8,2)=(1.d0-r**2.d0)/2.d0
              elseif (nnds.eq.4) then
                p(1,1)=-(1.d0-s)/4.d0
                p(2,1)=(1.d0-s)/4.d0
                p(3,1)=-(1.d0+s)/4.d0
                p(4,1)=(1.d0+s)/4.d0
                p(1,2)=-(1.d0-r)/4.d0
                p(2,2)=-(1.d0+r)/4.d0
                p(3,2)=(1.d0-r)/4.d0
                p(4,2)=(1.d0+r)/4.d0
              elseif (nnds.eq.9) then
                p(1,1)=((-1.d0+2.d0*r)*(-1.d0+s)*s)/4.d0
                p(2,1)=((1.d0+2.d0*r)*(-1.d0+s)*s)/4.d0
                p(3,1)=((-1.d0+2.d0*r)*s*(1.d0+s))/4.d0
                p(4,1)=((1.d0+2.d0*r)*s*(1.d0+s))/4.d0
                p(5,1)=-(r*(-1.d0+s)*s)
                p(6,1)=-((-1.d0+2.d0*r)*(-1.d0+s**2))/2.d0
                p(7,1)=-((1.d0+2.d0*r)*(-1.d0+s**2))/2.d0
                p(8,1)=-(r*s*(1.d0+s))
                p(9,1)=2.d0*r*(-1.d0+s**2)
                p(1,2)=((-1.d0+r)*r*(-1.d0+2.d0*s))/4.d0
                p(2,2)=(r*(1.d0+r)*(-1.d0+2.d0*s))/4.d0
                p(3,2)=((-1.d0+r)*r*(1.d0+2.d0*s))/4.d0
                p(4,2)=(r*(1.d0+r)*(1.d0+2.d0*s))/4.d0
                p(5,2)=-((-1.d0+r**2)*(-1.d0+2.d0*s))/2.d0
                p(6,2)=-((-1.d0+r)*r*s)
                p(7,2)=-(r*(1.d0+r)*s)
                p(8,2)=-((-1.d0+r**2)*(1.d0+2.d0*s))/2.d0
                p(9,2)=2.d0*(-1.d0+r**2)*s
              elseif (nnds.eq.6) then
                p(1,1)=-((-1.d0+s)*(2.d0*r+s))/4.d0
                p(2,1)=-((2.d0*r-s)*(-1.d0+s))/4.d0
                p(3,1)=((2.d0*r-s)*(1.d0+s))/4.d0
     &                +((1.d0+s)*(2.d0*r+s))/4.d0
     &                -(r*(1.d0+s))
                p(4,1)=r*(-1.d0+s)
                p(5,1)=(-1.d0+s**2.d0)/2.d0
                p(6,1)=(1.d0-s**2.d0)/2.d0
                p(1,2)=-((-1.d0+r)*(r+2.d0*s))/4.d0
                p(2,2)=-((1.d0+r)*(r-2.d0*s))/4.d0
                p(3,2)=((-1.d0+r)*(r-2.d0*s))/4.d0
     &                +((1.d0+r)*(r+2.d0*s))/4.d0
     &                +(1.d0-r**2.d0)/2.d0
                p(4,2)=(-1.d0+r**2.d0)/2.d0
                p(5,2)=(-1.d0+r)*s
                p(6,2)=-((1.d0+r)*s)
              endif
 
              do 402 i=1,ndim
                do 410 j=1,ndim
                  aj(i,j)=0.d0
                  do 420 k=1,nnds
                    aj(i,j)=aj(i,j)+p(k,i)*node(ele(g_ele(m),k),j)
 420              continue
 410            continue
 402          continue
 
              if (ndim.ne.2) then
                print*,'Number of dimensions is not supported'
              endif
              detp=aj(1,1)*aj(2,2)-aj(1,2)*aj(2,1)
              
              ted=1.d0/detp
              
              aji(1,1)=ted*aj(2,2)
              aji(1,2)=-ted*aj(1,2)
              aji(2,1)=-ted*aj(2,1)
              aji(2,2)=ted*aj(1,1)
              
              do 500 i=1,ndim
                do 510 j=1,nnds
                  B(m,j,k123,i)=0.d0
                  do 520 k=1,ndim
                    B(m,j,k123,i)=B(m,j,k123,i)+aji(i,k)*p(j,k)
 520              continue
 510            continue
 500          continue
             
              det(m,k123)=detp*ww(k123)

              if (det(m,k123).lt.0.d0) then
                print*,'ERROR: Negative determinant in',
     &                 ' Element: ',g_ele(m),k123,det(m,k123)
                iserr=.true.
              endif
              
 320        continue
 310      continue

        endif
 
 100  continue

      if (iserr) stop

      
      return
      end

***********************************************************************

      subroutine applybcu(ndpn,ndim,maxfc,fbcnd,fbcdof,maxnd,nnds,
     &                    nfbcs,gconu,fbcval,f,nudofs,nodes,ele,
     &                    local_eles,g_ele,maxprocs,rank,maxele)
      implicit none
      real*8 fbcval(maxfc),f(maxnd*ndpn)
      integer maxfc,fbcnd(maxfc),fbcdof(maxfc),nnds
      integer nfbcs,ndof,gconu(maxnd,ndpn),maxnd
      integer i,nudofs,ndpn,nodes,ndim,rank,maxprocs,maxele
      integer local_eles,g_ele(maxele),ele(maxele*maxprocs,nnds)


      do 50 i=1,nodes*ndim
        f(i)=0.d0
 50   continue

      if (rank.eq.0) then
      do 100 i=1,nfbcs
        if (fbcdof(i).le.ndim) then
          ndof=gconu(fbcnd(i),fbcdof(i))
          f(ndof)=f(ndof)+fbcval(i)
        endif
 100  continue
      endif

      return
      end

***********************************************************************
      
      subroutine applybcmu(nnds,nints,ndpn,ndim,maxfc,fbcnd,fbcdof,
     &                   maxnd,nodes,node,nfbcs,gconmu,fbcval,f,dt,
     &                   beta,maxele,neles,ele,det,NN,muint,nmudofs,
     &                   local_eles,g_ele,maxprocs,rank)
      implicit none
      real*8 fbcval(maxfc),f(maxnd*ndpn),node(maxnd,ndim)
      real*8 NN(nnds,nints),muint(maxele,nints,2)
      real*8 fele(nnds),det(maxele,nints),dt,beta
      integer maxfc,fbcnd(maxfc)
      integer nfbcs,ndof,nodes,gconmu(maxnd,ndpn),maxnd
      integer i,j,k,l,nmudofs,fbcdof(maxfc),maxprocs
      integer maxele,neles,nnds,nints,ndim,ndpn,rank
      integer local_eles,g_ele(maxele),ele(maxele*maxprocs,nnds)

      
      do 50 i=1,nodes
        f(i)=0.d0
 50   continue

      if (rank.eq.0) then 
      do 100 i=1,nfbcs
        if (fbcdof(i).eq.ndpn) then
          ndof=gconmu(fbcnd(i),1)
          f(ndof)=f(ndof)+fbcval(i)
        endif
 100  continue
      endif

      do 200 i=1,local_eles

        do 210 j=1,nnds
 210      fele(j)=0.d0

        do 220 k=1,nints
         
          do 230 l=1,nnds
            fele(l)=fele(l)+det(i,k)*NN(l,k)*0.d0*beta/dt*muint(i,k,1)
 230      continue
 
 220    continue

        do 300 j=1,nnds
          ndof=gconmu(ele(g_ele(i),j),1)
          f(ndof)=f(ndof)+fele(j)
 300    continue
 
 200  continue


      return
      end

***********************************************************************
      
      subroutine applybcvp(nnds,nints,ndpn,ndim,maxfc,fbcnd,fbcdof,
     &                   maxnd,nodes,node,nfbcs,gconvp,fbcval,f,dt,
     &                   maxele,neles,ele,det,NN,mu,nmudofs,u,papp,e0,
     &                   unitfbc,vapp,fluidsource,w0,gradv,Finv,Jacob,
     &                   timerule,local_eles,g_ele,maxprocs,rank)
      implicit none
      real*8 fbcval(maxfc),f(maxnd*ndpn),node(maxnd,ndim)
      real*8 papp,e0,factor,unitfbc(maxfc),vapp,phiterm
      real*8 NN(nnds,nints),u(maxnd,ndim),length1,length2
      real*8 fele(nnds*(ndim+1)),det(maxele,nints),mu(maxnd),dt
      real*8 fluidsource(maxele),w0,gradv(maxele,nints,ndim,ndim,2)
      real*8 Finv(maxele,nints,ndim,ndim,2),Jacob(maxele,nints,2)
      real*8 timerule
      integer maxfc,fbcnd(maxfc),nbound(maxele*maxprocs,2),ncount
      integer nfbcs,ndof,nodes,gconvp(maxnd,ndpn),maxnd
      integer i,j,k,l,nmudofs,fbcdof(maxfc),a,b,maxprocs,m,n
      integer maxele,neles,nnds,nints,ndim,ndpn,rank,is9
      integer local_eles,g_ele(maxele),ele(maxele*maxprocs,nnds)

      is9=0
      if (nnds.eq.9) is9=1

      do 50 i=1,nodes*(ndim+1)
        f(i)=0.d0
 50   continue

      do 380 i=1,neles
        ncount=0
        do 381 j=1,(nnds-is9)/2
          do 382 k=1,nfbcs
            if (fbcdof(k).eq.3) then
            if (ele(i,j).eq.fbcnd(k)) then
              ncount=ncount+1
              nbound(i,ncount)=ele(i,j)
            endif
            endif
            if (ncount.gt.2) then
              print*,'Too many fbcs in ele: ',i
              stop
            endif
 382      continue
 381    continue
        if (ncount.eq.1) then
!         This will result in zero force...
          nbound(i,2)=nbound(i,1) 
        endif
 380  continue

      do 390 i=1,nfbcs
        if (fbcdof(i).eq.3) then
        factor=0.d0
        if (unitfbc(i).lt.1.d-10) then
          factor=-1.d0*papp*dabs(node(fbcnd(i),1))/node(fbcnd(i),1)
        else if (unitfbc(i).lt.1.d0+1.d-10) then
          factor=-0.9d0
        endif
        fbcval(i)=0.d0
        do 391 j=1,neles
          do 392 k=1,nnds-is9
            if (ele(j,k).eq.fbcnd(i)) then
              a=nbound(j,1)
              b=nbound(j,2)
              length1=dsqrt((node(a,2)+e0*u(a,2)
     &                      -node(b,2)-e0*u(b,2))**2.d0
     &                     +(node(a,1)+e0*u(a,1)
     &                      -node(b,1)-e0*u(b,1))**2.d0)
              if (k.le.(nnds-is9)/2) then
                fbcval(i)=fbcval(i)+length1/6.d0
              else
                fbcval(i)=fbcval(i)+2.d0*length1/3.d0
              endif
            endif
 392      continue
 391    continue
        fbcval(i)=fbcval(i)*factor
       endif
 390  continue

      do 480 i=1,neles
        ncount=0
        do 481 j=1,(nnds-is9)/2
          do 482 k=1,nfbcs
            if (fbcdof(k).eq.4) then
            if (ele(i,j).eq.fbcnd(k)) then
              ncount=ncount+1
              nbound(i,ncount)=ele(i,j)
            endif
            endif
            if (ncount.gt.2) then
              print*,'Too many fbcs in ele: ',i
              stop
            endif
 482      continue
 481    continue
        if (ncount.eq.1) then
!         This will result in zero force...
          nbound(i,2)=nbound(i,1) 
        endif
 480  continue

      do 490 i=1,nfbcs
        if (fbcdof(i).eq.4) then
        factor=0.d0
        if (unitfbc(i).lt.1.d-10) then
          factor=-1.d0*papp*dabs(node(fbcnd(i),2))/node(fbcnd(i),2)
        else if (unitfbc(i).lt.1.d0+1.d-10) then
          factor=-0.9d0
        endif
        fbcval(i)=0.d0
        do 491 j=1,neles
          do 492 k=1,nnds
            if (ele(j,k).eq.fbcnd(i)) then
              a=nbound(j,1)
              b=nbound(j,2)
              length1=dsqrt((node(a,2)+e0*u(a,2)
     &                      -node(b,2)-e0*u(b,2))**2.d0
     &                     +(node(a,1)+e0*u(a,1)
     &                      -node(b,1)-e0*u(b,1))**2.d0)
              if (k.le.(nnds-is9)/2) then
                fbcval(i)=fbcval(i)+length1/6.d0
              else
                fbcval(i)=fbcval(i)+2.d0*length1/3.d0
              endif
            endif
 492      continue
 491    continue
        fbcval(i)=fbcval(i)*factor
       endif
 490  continue 
 
      do 200 i=1,local_eles

        do 210 j=1,nnds*(ndim+1)
 210      fele(j)=0.d0

        do 220 k=1,nints 
          do 230 l=1,nnds
              fele(3*l)=fele(3*l)
     &           +det(i,k)*NN(l,k)*(-vapp)*fluidsource(i)/w0

            do 231 m=1,ndim
              do 232 n=1,ndim
                fele(3*l)=fele(3*l)+det(i,k)*NN(l,k)*
     &                    Jacob(i,k,2)*gradv(i,k,m,N,2)*Finv(i,k,N,m,2)
     &                    *(1.d0-timerule)
 232          continue
 231        continue

 230      continue
 220    continue

        do 300 j=1,nnds
          do 310 k=1,ndim+1
            ndof=gconvp(ele(g_ele(i),j),k)
            f(ndof)=f(ndof)+fele(3*(j-1)+k)
 310      continue
 300    continue

 200  continue

      if (rank.eq.0) then
      do 100 i=1,nfbcs
        if (fbcdof(i).gt.ndim.and.fbcdof(i).lt.ndpn) then
          ndof=gconvp(fbcnd(i),fbcdof(i)-ndim)
          f(ndof)=f(ndof)+fbcval(i)
        endif
 100  continue
      endif


      return
      end

***********************************************************************

      subroutine calcu(ndim,ndpn,maxnd,maxdbc,nodes,nudofs,nconu,
     &                 gconu,dbcval,drhs,u,umpc,nudbcs,dbcnd,dbcdof,
     &                 conAu,conBu,umpcnd,nmudbcs,nvpdbcs,umpcnum,
     &                 maxcon)
      implicit none
      real*8 drhs(maxnd*ndpn),u(maxnd,ndim)
      real*8 dbcval(maxdbc),conAu(maxdbc,maxcon),conBu(maxdbc,maxcon)
      integer maxnd,nodes,umpc(maxdbc,maxcon,2),dbcnd(maxdbc)
      integer dbcdof(maxdbc),maxdbc,maxcon,umpcnum(maxdbc)
      integer i,j,gconu(maxnd,ndpn),ndof,nconu,nudofs,nudbcs
      integer umpcnd(maxdbc,2),nmudbcs,ndim,ndpn,nvpdbcs


      do 100 i=1,nodes
        do 200 j=1,ndim
          ndof=gconu(i,j)
          if (ndof.le.nudofs) then
            u(i,j)=u(i,j)+drhs(ndof)
          endif
 200    continue
 100  continue

      do 250 i=1,nudbcs+nmudbcs+nvpdbcs
        ndof=dbcdof(i)
        if (ndof.le.ndim) then
          u(dbcnd(i),ndof)=dbcval(i)
        endif
 250  continue
 
      do 300 i=1,nconu
      
        if (umpcnd(i,2).eq.1) then
          u(umpcnd(i,1),1)=0.d0
          do 310 j=1,umpcnum(i)
            if(umpc(i,j,2).eq.1) u(umpcnd(i,1),1)
     &        =u(umpcnd(i,1),1)+conAu(i,j)*u(umpc(i,j,1),1)+conBu(i,j)
            if(umpc(i,j,2).eq.2) u(umpcnd(i,1),1)
     &        =u(umpcnd(i,1),1)+conAu(i,j)*u(umpc(i,j,1),2)+conBu(i,j)
 310      continue
        endif
        
        if (umpcnd(i,2).eq.2) then
          u(umpcnd(i,1),2)=0.d0
          do 320 j=1,umpcnum(i)
            if(umpc(i,j,2).eq.1) u(umpcnd(i,1),2)
     &        =u(umpcnd(i,1),2)+conAu(i,j)*u(umpc(i,j,1),1)+conBu(i,j)
            if(umpc(i,j,2).eq.2) u(umpcnd(i,1),2)
     &        =u(umpcnd(i,1),2)+conAu(i,j)*u(umpc(i,j,1),2)+conBu(i,j)
 320      continue
        endif
        
 300  continue

        
      return
      end

***********************************************************************

      subroutine calcmu(ndpn,ndim,maxnd,maxdbc,maxcon,nodes,nmudofs,
     &                 gconmu,dbcval,drhs,mu,nudbcs,dbcnd,dbcdof,
     &                 mumpcnum,nmudbcs,nconmu,mumpc,mumpcnd,
     &                 conAmu,conBmu,nvpdbcs)
      implicit none
      real*8 drhs(maxnd*ndpn),mu(maxnd)
      real*8 dbcval(maxdbc),conAmu(maxdbc,maxcon)
      real*8 conBmu(maxdbc,maxcon)
      integer maxnd,maxdbc,maxcon,nodes,dbcnd(maxdbc),dbcdof(maxdbc)
      integer i,gconmu(maxnd,ndpn),ndof,nmudofs,nudbcs
      integer mumpcnd(maxdbc,2),mumpc(maxdbc,maxcon,2)
      integer mumpcnum(maxdbc),ndpn,j,nvpdbcs,ndim,nmudbcs,nconmu


      do 100 i=1,nodes
        ndof=gconmu(i,1)
        if (ndof.le.nmudofs) then
          mu(i)=mu(i)+drhs(ndof)
        endif
 100  continue

      do 250 i=1,nudbcs+nmudbcs+nvpdbcs
        ndof=dbcdof(i)
        if (ndof.eq.ndpn) then
          mu(dbcnd(i))=dbcval(i)
        endif
 250  continue
 
      do 300 i=1,nconmu
      
       mu(mumpcnd(i,1))=0.d0
       do 340 j=1,mumpcnum(i)
        if(mumpc(i,j,2).eq.ndpn) 
     &     mu(mumpcnd(i,1))=mu(mumpcnd(i,1))
     &     +conAmu(i,j)*mu(mumpc(i,j,1))+conBmu(i,j)
 340   continue
        
 300  continue

        
      return
      end

**********************************************************************
      
      subroutine calcvp(ndpn,ndim,maxnd,maxdbc,maxcon,nodes,nvpdofs,
     1                  gconvp,dbcval,drhs,v,p,nudbcs,dbcnd,dbcdof,
     2                  vpmpcnum,nvpdbcs,nconvp,vpmpc,vpmpcnd,
     3                  conAvp,conBvp,nmudbcs)
      implicit none
      real*8 drhs(maxnd*ndpn),v(maxnd,ndim),p(maxnd)
      real*8 dbcval(maxdbc),conAvp(maxdbc,maxcon)
      real*8 conBvp(maxdbc,maxcon)
      integer maxnd,maxdbc,maxcon,nodes,dbcnd(maxdbc),dbcdof(maxdbc)
      integer i,gconvp(maxnd,ndpn),ndof,nvpdofs,nudbcs,nvpdbcs,nconvp
      integer vpmpcnd(maxdbc,2),vpmpc(maxdbc,maxcon,2)
      integer vpmpcnum(maxdbc),ndpn,j,nmudbcs,ndim


      do 100 i=1,nodes
        do 110 j=1,ndim
          ndof=gconvp(i,j)
          if (ndof.le.nvpdofs) then
            v(i,j)=v(i,j)+drhs(ndof)
          endif
 110    continue
        ndof=gconvp(i,3)
        if (ndof.le.nvpdofs) then
          p(i)=p(i)+drhs(ndof)
        endif
 100  continue

      do 250 i=1,nudbcs+nmudbcs+nvpdbcs
        ndof=dbcdof(i)-ndim
        if (ndof.gt.0.and.ndof.le.ndim) then
          v(dbcnd(i),ndof)=dbcval(i)
        endif
        if (ndof.eq.3) then
          p(dbcnd(i))=dbcval(i)
        endif
 250  continue
 
      do 300 i=1,nconvp
      
        if (vpmpcnd(i,2).eq.3) then
          v(vpmpcnd(i,1),1)=0.d0
          do 310 j=1,vpmpcnum(i)
            if(vpmpc(i,j,2).eq.3) v(vpmpcnd(i,1),1)
     &        =v(vpmpcnd(i,1),1)+conAvp(i,j)*v(vpmpc(i,j,1),1)+
     &         conBvp(i,j)
            if(vpmpc(i,j,2).eq.4) v(vpmpcnd(i,1),1)
     &        =v(vpmpcnd(i,1),1)+conAvp(i,j)*v(vpmpc(i,j,1),2)+
     &         conBvp(i,j)
 310      continue
        endif
        
        if (vpmpcnd(i,2).eq.4) then
          v(vpmpcnd(i,1),2)=0.d0
          do 320 j=1,vpmpcnum(i)
            if(vpmpc(i,j,2).eq.3) v(vpmpcnd(i,1),2)
     &        =v(vpmpcnd(i,1),2)+conAvp(i,j)*v(vpmpc(i,j,1),1)+
     &         conBvp(i,j)
            if(vpmpc(i,j,2).eq.4) v(vpmpcnd(i,1),2)
     &        =v(vpmpcnd(i,1),2)+conAvp(i,j)*v(vpmpc(i,j,1),2)+
     &         conBvp(i,j)
 320      continue
        endif
        
        if (vpmpcnd(i,2).eq.5) then
          p(vpmpcnd(i,1))=0.d0
          do 330 j=1,vpmpcnum(i)
            p(vpmpcnd(i,1))
     &        =p(vpmpcnd(i,1))+conAvp(i,j)*p(vpmpc(i,j,1))+
     &         conBvp(i,j)
 330      continue
        endif
        
 300  continue

        
      return
      end
      
**********************************************************************

      subroutine gradu(ndim,nnds,nints,maxele,maxnd,local_eles,ele,
     &                 B,NN,strn,u,Jacob,Finv,dJacob,dstrn,e0,dt,uint,
     &                 duint,FiJ,muint,phi,open_thresh,g_ele,maxprocs,
     &                 maxJ)
      implicit none
      real*8 strn(maxele,nints,ndim,ndim,2),dt,NN(nnds,nints)
      real*8 dstrn(maxele,nints,ndim,ndim),maxJ(maxele)
      real*8 u(maxnd,ndim),Jacob(maxele,nints,2)
      real*8 Finv(maxele,nints,ndim,ndim,2),e0
      real*8 B(maxele,nnds,nints,ndim),dJacob(maxele,nints)
      real*8 uint(maxele,nints,ndim),duint(maxele,nints,ndim)
      real*8 FiJ(maxele,nints,ndim,ndim),muint(maxele,nints,2)
      real*8 delta(ndim,ndim),isalt,phi(maxele,nints,2)
      real*8 h0,open_thresh
      integer maxele,maxnd,local_eles,ele(maxele*maxprocs,nnds)
      integer i,j,k,m,ndim,nnds,nints,maxprocs,g_ele(maxele),ii,jj
      integer Sint


      h0=open_thresh

      do 10 i=1,ndim
        do 11 j=1,ndim
          delta(i,j)=0.d0
          if (i.eq.j) delta(i,j)=1.d0
 11     continue
 10   continue

      do 100 m=1,local_eles
        
        maxJ(m)=0.d0
        Sint=5

        do 110 j=1,nints

          do 150 i=1,ndim
            uint(m,j,i)=0.d0
            do 151 k=1,ndim
              strn(m,j,i,k,1)=0.d0
 151        continue   
 150      continue

          do 200 i=1,nnds
            k=ele(g_ele(m),i)
            do 201 ii=1,ndim
              uint(m,j,ii)=uint(m,j,ii)+NN(i,j)*u(k,ii)
              do 202 jj=1,ndim
                strn(m,j,ii,jj,1)=strn(m,j,ii,jj,1)
     &                            +B(m,i,j,jj)*u(k,ii)
 202          continue
 201        continue
 200      continue

          if (strn(m,j,1,1,1)+strn(m,j,2,2,1).gt.maxJ(m)) then
            maxJ(m)=strn(m,j,1,1,1)+strn(m,j,2,2,1)
            Sint=j
          endif

 110    continue

        maxJ(m)=0.d0
        do 210 j=1,nints

          do 250 i=1,ndim
            do 251 k=1,ndim
              dstrn(m,j,i,k)=(strn(m,j,i,k,1)-strn(m,j,i,k,2))/dt
!              dstrn(m,j,i,k)=0.d0
              FiJ(m,j,i,K)=delta(i,K)+e0*strn(m,j,i,K,1)
 251        continue   
 250      continue

          Jacob(m,j,1)=FiJ(m,j,1,1)*FiJ(m,j,2,2)
     &              -FiJ(m,j,1,2)*FiJ(m,j,2,1)

          if (Jacob(m,j,1).gt.maxJ(m)) maxJ(m)=Jacob(m,j,1)
!          maxJ(m)=maxJ(m)+Jacob(m,j,1)
          if (Jacob(m,j,1).lt.1.d-15) then
!            print*,'Negative Jacobian: ',m,j
!            stop
          endif
          Finv(m,j,1,1,1)=FiJ(m,j,2,2)/Jacob(m,j,1)
          Finv(m,j,1,2,1)=-FiJ(m,j,1,2)/Jacob(m,j,1)
          Finv(m,j,2,1,1)=-FiJ(m,j,2,1)/Jacob(m,j,1)
          Finv(m,j,2,2,1)=FiJ(m,j,1,1)/Jacob(m,j,1)

          dJacob(m,j)=(Jacob(m,j,1)-Jacob(m,j,2))/dt

 210    continue

        maxJ(m)=maxJ(m)*0.3d0+Jacob(m,5,1)*0.7d0

 100  continue


      return
      end

********************************************************************** 
      
      subroutine gradm(nnds,nints,ndim,maxele,maxnd,local_eles,
     &                 ele,B,NN,muint,gradmu,mu,dt,muavg,dmuint,
     &                 g_ele,maxprocs)
      implicit none
      real*8 muint(maxele,nints,2),mu(maxnd)
      real*8 gradmu(maxele,nints,ndim)
      real*8 B(maxele,nnds,nints,ndim),NN(nnds,nints)
      real*8 dt,muavg(maxele),dmuint(maxele,nints)
      integer maxele,maxnd,local_eles,ele(maxele*maxprocs,nnds)
      integer i,j,k,m,nnds,nints,ndim,maxprocs,g_ele(maxele),ii


      do 100 m=1,local_eles
        
        muavg(m)=0.d0
        do 110 j=1,nints

          do 165 i=1,ndim
            gradmu(m,j,i)=0.d0
 165      continue

          muint(m,j,1)=0.d0

          do 200 i=1,nnds
            k=ele(g_ele(m),i)
            do 201 ii=1,ndim
              gradmu(m,j,ii)=gradmu(m,j,ii)+B(m,i,j,ii)*mu(k)
 201        continue
            muint(m,j,1)=muint(m,j,1)+NN(i,j)*mu(k)
 200      continue
          muavg(m)=muavg(m)+muint(m,j,1)

          dmuint(m,j)=(muint(m,j,1)-muint(m,j,2))/dt

 110    continue

        muavg(m)=muavg(m)/nints
 100  continue


      return
      end

********************************************************************** 
      
      subroutine normmu(nnds,nints,ndim,maxele,maxnd,local_eles,
     &                 ele,B,NN,ngradmu,gmuint,g_ele,maxprocs)
      implicit none
      real*8 B(maxele,nnds,nints,ndim),NN(nnds,nints)
      real*8 Finv(maxele,nints,ndim,ndim,2),ngradmu(maxnd,ndim)
      real*8 gmuint(maxele,nints,ndim),gmuintc(maxele,nints,ndim)
      real*8 magmubig,magmusmall
      integer maxele,maxnd,local_eles,ele(maxele*maxprocs,nnds)
      integer i,j,k,m,nnds,nints,ndim,p,q,maxprocs,g_ele(maxele)


      do 100 m=1,local_eles
      
        do 110 j=1,nints

          do 165 i=1,ndim
            gmuint(m,j,i)=0.d0
 165      continue

          do 200 i=1,nnds
            k=ele(g_ele(m),i)
            do 201 q=1,ndim
              gmuint(m,j,q)=gmuint(m,j,q)+NN(i,j)*ngradmu(k,q)
 201        continue
 200      continue

 110    continue
 100  continue


      return
      end

**********************************************************************     
      subroutine gradvp(nnds,nints,ndim,maxele,maxnd,local_eles,
     &                  ele,B,NN,gradv,vint,gradp,pint,v,p,g_ele,
     &                  maxprocs)
      implicit none
      real*8 vint(maxele,nints,ndim),v(maxnd,ndim)
      real*8 gradv(maxele,nints,ndim,ndim,2),p(maxnd)
      real*8 gradp(maxele,nints,ndim),pint(maxele,nints,2)
      real*8 B(maxele,nnds,nints,ndim),NN(nnds,nints)
      integer maxele,maxnd,local_eles,ele(maxele*maxprocs,nnds)
      integer i,j,k,l,m,nnds,nints,ndim,maxprocs,g_ele(maxele),II


      do 100 m=1,local_eles
      
        do 110 j=1,nints

          pint(m,j,1)=0.d0
          
          do 165 i=1,ndim
            vint(m,j,i)=0.d0
            gradp(m,j,i)=0.d0
            do 166 k=1,ndim
              gradv(m,j,i,k,1)=0.d0
 166        continue
 165      continue

          do 200 i=1,nnds
            k=ele(g_ele(m),i)
            pint(m,j,1)=pint(m,j,1)+NN(i,j)*p(k)
            do 210 l=1,ndim
              do 211 II=1,ndim
               gradv(m,j,l,II,1)=gradv(m,j,l,II,1)
     &                           +B(m,i,j,II)*v(k,l)
 211          continue
              vint(m,j,l)=vint(m,j,l)+NN(i,j)*v(k,l)
              gradp(m,j,L)=gradp(m,j,L)+B(m,i,j,L)*p(k)
 210        continue
 200      continue

 110    continue
 100  continue


      return
      end

**********************************************************************

      subroutine stresses(nints,ndim,maxele,local_eles,strn,gradv,
     &                    gradmu,pint,gradp,vint,strs,tau,bsf,bp,bv,
     &                    psimu,ksi,phi0,Jacob,Finv,Cstiff,w0,
     &                    Imperm,dJacob,dphi,dstrn,Hpsi,e0,duint,
     &                    isol,a,phi,visc_eff,muint,g_ele,bphiJ,
     &                    dHds,Hpsi0,biot,dt,freestrain,YM,poisson,
     &                    tc_split_on,muavg,kperm,Wp,uint,Wmu,Ws,
     &                    dmuint,fluidsource,rank,open_thresh,
     &                    gmuint,scale_visc,stretch,maxprocs,Cf,
     &                    p0,maxJ,dHdp)
      implicit none
      real*8 strn(maxele,nints,ndim,ndim,2),e0
      real*8 gradv(maxele,nints,ndim,ndim,2),bsf(maxele,nints,ndim)
      real*8 muint(maxele,nints,2),gradmu(maxele,nints,ndim)
      real*8 bp(maxele,nints,ndim),bv(maxele,nints)
      real*8 strs(maxele,nints,ndim,ndim),vint(maxele,nints,ndim)
      real*8 pint(maxele,nints,2),gradp(maxele,nints,ndim)
      real*8 psimu(maxele,nints),tau(maxele,nints,ndim,ndim)
      real*8 ksi(maxele,nints,ndim),Finv(maxele,nints,ndim,ndim,2)
      real*8 Jacob(maxele,nints,2),mu,Cstiff(4)
      real*8 Imperm(maxele,ndim,ndim),bphiJ(maxele,nints)
      real*8 c11,c13,c33,c44,w0,vbar(ndim,ndim),vkk
      real*8 dJacob(maxele,nints),dphi(maxele,nints),test
      real*8 dstrn(maxele,nints,ndim,ndim),Hpsi(maxele,nints)
      real*8 duint(maxele,nints,ndim),delta(ndim,ndim),factor,a
      real*8 dfactor,phi(maxele,nints,2)
      real*8 Cphi,p,theta,pi,nxn(ndim,ndim),ekk
      real*8 eplus(ndim,ndim),eminus(ndim,ndim),ekkplus,ekkminus
      real*8 strnplus(ndim,ndim),strnminus(ndim,ndim),Ba,h0
      real*8 visc_eff(maxele,nints),dHds(maxele,nints,ndim,ndim)
      real*8 Hpsi0(maxele,nints),Bc,phicrack,phibiot,dphidstrain
      real*8 biot(maxele),phi0(maxele),dphibiot
      real*8 freestrain(maxele),YM(maxele)
      real*8 poisson(maxele),dt,dkkplus
      real*8 aexp,bexp,phiF(2),darcyF,muavg(maxele),kperm,phieff(2)
      real*8 Wp(maxele,nints),uint(maxele,nints,ndim),DphiF(2)
      real*8 Wmu(maxele,nints),Ws(maxele,nints),dmuint(maxele,nints)
      real*8 fluidsource(maxele),crit_strain
      real*8 gmuint(maxele,nints,ndim),p0,dHdp(maxele,nints)
      real*8 HeavisideJ1,HeavisideJ2,istr,maxJ(maxele)
      real*8 stretch(maxele,nints),Cf,fdensity(2),gfdensity(ndim)
      real*8 s_numerator,s_denominator,open_thresh,dHJ1
      real*8 crit_stretch,Fstretch,Fstretch_norm
      integer maxele,local_eles,nints,ndim,isol,H_total
      integer m,l,i,j,k,n,q,rank,maxprocs,g_ele(maxele)
      logical tc_split_on,scale_visc


      c11=Cstiff(1)
      c33=Cstiff(2)
      c13=Cstiff(3)
      c44=Cstiff(4)

      pi=4.d0*datan(1.d0)

      do 10 i=1,ndim
        do 11 j=1,ndim
          delta(i,j)=0.d0
          eplus(i,j)=0.d0
          eminus(i,j)=0.d0
          if (i.eq.j) delta(i,j)=1.d0
 11     continue
 10   continue

      H_total=0
      
      do 100 m=1,local_eles
        do 200 l=1,nints

          c11=YM(m)*(1.d0-poisson(m))
     &        /((1.d0+poisson(m))*(1.d0-2.d0*poisson(m)))
          c33=c11
          c13=YM(m)*poisson(m)
     &        /((1.d0+poisson(m))*(1.d0-2.d0*poisson(m)))
          c44=YM(m)
     &        /(2.d0*(1.d0+poisson(m)))
          
          if (dabs(biot(m)-1.d0).lt.1.d-15) then
            Cphi=1.d20
          elseif (dabs(biot(m)-phi0(m)).lt.1.d-15) then
            Cphi=1.d20
          else 
            Cphi=(YM(m)/(e0*3.d0*(1.d0-2.d0*poisson(m))))
     &           /((biot(m)-phi0(m))*(1.d0-biot(m)))
          endif

          s_numerator=0.d0
          s_denominator=0.d0
          do 311 i=1,ndim
            s_numerator=s_numerator+gmuint(m,l,i)*gmuint(m,l,i)
            do 321 j=1,ndim
              do 331 k=1,ndim
                s_denominator=s_denominator
     &    +Finv(m,l,J,i,1)*Finv(m,l,K,i,1)*gmuint(m,l,J)*gmuint(m,l,K)
 331          continue             
 321        continue
 311      continue
          if (s_denominator.lt.1.d-15) then
            stretch(m,l)=1.d0
          else
            stretch(m,l)=dsqrt(s_numerator/s_denominator)
          endif
          stretch(m,l)=maxJ(m)*1.d0+Jacob(m,l,1)*0.d0
          crit_stretch=1.d0+open_thresh
          
          if (.not.scale_visc) then
            Fstretch=1.d0
          elseif (stretch(m,l).le.1.d0) then
            Fstretch=(crit_stretch/(crit_stretch-1.d0))**(3.d0)
          elseif (stretch(m,l).lt.2.d0*crit_stretch-1.d0) then
            Fstretch=((stretch(m,l)*(stretch(m,l)-2.d0)
     &                 +(1.d0-2.d0*crit_stretch)**2.d0)
     &                /(5.d0+stretch(m,l)*(stretch(m,l)-2.d0)
     &                  +4.d0*crit_stretch*(crit_stretch-2.d0)))**(3.d0)
          else
            Fstretch=(stretch(m,l)/(stretch(m,l)-1.d0))**(3.d0)
          endif
          Fstretch_norm=(crit_stretch/(crit_stretch-1.d0))**(3.d0)
!         Fstretch_norm turns out to be 1 + 1/(kperm*kappa) such that
!          the flow equation is normalized w0
 
          mu=muint(m,l,1)
          p=pint(m,l,1)

          if (a.le.0.d0) then
!          Cubic
           factor=-mu*mu*(2.d0*mu-3.d0+a*(mu-1.d0))
           dfactor=(6.d0*(1.d0-mu)+a*(2.d0-3.d0*mu))*mu
          else
!          Bell-Curve
           factor=(1.d0+a)*(1.d0-(a/(1.d0+a))**(mu*mu))
           dfactor=2.d0*(factor-(1.d0+a))*(mu)*dlog(a/(1.d0+a))
          endif

!         Fluid-mass indicator function for current and past time step
          phiF(1)=1.d0-3.d0*mu**2.d0+2.d0*mu**3.d0
          phiF(2)=1.d0-3.d0*muint(m,l,2)**2.d0+2.d0*muint(m,l,2)**3.d0
          DphiF(1)=-6.d0*mu*(1.d0-mu)
          DphiF(2)=-6.d0*muint(m,l,2)*(1.d0-muint(m,l,2))

!         Darcy degradation function
          aexp=4.d0
          bexp=4.d0
          darcyF=(1.d0-muavg(m)**aexp)**bexp+0.d-6
!          darcyF=(1.d0-mu**aexp)**bexp+0.d-6

!         Heaviside function for J
          if (Jacob(m,l,1).ge.1.d0) then
            HeavisideJ1=1.d0
            H_total=H_total+1
          else
            HeavisideJ1=1.d0
            if (tc_split_on) HeavisideJ1=0.d0
            if (fluidsource(m).gt.1.d-6) HeavisideJ1=1.d0
          endif
          dHJ1=0.d0

          if (Jacob(m,l,2).ge.1.d0) then
            HeavisideJ2=1.d0
          else
            HeavisideJ2=1.d0
            if (tc_split_on) HeavisideJ2=0.d0
            if (fluidsource(m).gt.1.d-6) HeavisideJ2=1.d0
          endif

          visc_eff(m,l)=Fstretch/(Fstretch_norm-1.d0)

          vkk=0.d0
          ekk=0.d0
          do 212 i=1,ndim
            do 214 j=1,ndim
            vkk=vkk+(gradv(m,l,j,I,1)+e0*dstrn(m,l,j,I)/w0)
     &              *Finv(m,l,I,j,1)
 214        continue
            ekk=ekk+strn(m,l,i,i,1)-freestrain(m)
 212      continue

          ekkplus=0.5d0*(ekk+dabs(ekk))
          ekkminus=0.5d0*(ekk-dabs(ekk))
          if (ekk.ge.0.d0) then
            dkkplus=1.d0
          else
            dkkplus=0.d0
          endif

          fdensity(1)=1.d0+Cf*p
          fdensity(2)=1.d0+Cf*pint(m,l,2)
          do 240 i=1,ndim
            gfdensity(I)=Cf*gradp(m,l,I)
 240      continue

          phicrack=1.0d0
          phibiot=phi0(m)/Jacob(m,l,1)
     &    +(biot(m)*(Jacob(m,l,1)-1.d0)+(p-p0)/Cphi)
     &     *(1.d0-phiF(1)*HeavisideJ1)
     &     /(Jacob(m,l,1)-phiF(1)*HeavisideJ1*(Jacob(m,l,1)-1.d0))

          phi(m,l,1)=phibiot

          phieff(1)=phi(m,l,1)+phiF(1)*(1.d0-phi(m,l,1))
     &                *((Jacob(m,l,1)-1.d0)*HeavisideJ1)/Jacob(m,l,1)

          phieff(2)=phi(m,l,2)+phiF(2)*(1.d0-phi(m,l,2))
     &                *((Jacob(m,l,2)-1.d0)*HeavisideJ2)/Jacob(m,l,2)
           
          dphi(m,l)=(phi(m,l,1)-phi(m,l,2))/dt

          do 213 i=1,ndim
            do 223 j=1,ndim
                if (i.eq.j) then
                  vbar(i,j)=-vkk/3.d0
                else
                  vbar(i,j)=0.d0
                endif
              do 233 k=1,ndim
                vbar(i,j)=vbar(i,j)
     &              +(gradv(m,l,i,K,1)+e0*dstrn(m,l,i,K)/w0)
     &               *Finv(m,l,K,j,1)
 233          continue
 223        continue
 213      continue

          if (dabs(strn(m,l,1,2,1)+strn(m,l,2,1,1)).lt.1.d-13) then
            if (strn(m,l,1,1,1).ge.strn(m,l,2,2,1)) then
              theta=0.d0
            else
              theta=pi/2.d0
            endif
          elseif (dabs(strn(m,l,1,1,1)-strn(m,l,2,2,1)).lt.1.d-13) then
            if (strn(m,l,1,2,1)+strn(m,l,2,1,1).ge.0.d0) then
              theta=pi/4.d0
            else
              theta=3.d0*pi/4.d0
            endif
          else
            theta=0.5d0*
     &                  datan2(strn(m,l,1,2,1)+strn(m,l,2,1,1),
     &                         strn(m,l,1,1,1)-strn(m,l,2,2,1))
          endif

          nxn(1,1)=dcos(theta)
          nxn(1,2)=-dsin(theta)
          nxn(2,1)=dsin(theta)
          nxn(2,2)=nxn(1,1)

          eplus(1,1)=0.5d0*(
     &                 0.5d0*(strn(m,l,1,1,1)+strn(m,l,2,2,1)
     &                                     -2.d0*freestrain(m)
     &                 +dsqrt((strn(m,l,1,1,1)-strn(m,l,2,2,1))**2.d0
     &                      +(strn(m,l,1,2,1)+strn(m,l,2,1,1))**2.d0))
     &           +dabs(0.5d0*(strn(m,l,1,1,1)+strn(m,l,2,2,1)
     &                                     -2.d0*freestrain(m)
     &                 +dsqrt((strn(m,l,1,1,1)-strn(m,l,2,2,1))**2.d0
     &                      +(strn(m,l,1,2,1)+strn(m,l,2,1,1))**2.d0)))
     &                     )
          eplus(2,2)=0.5d0*(
     &                 0.5d0*(strn(m,l,1,1,1)+strn(m,l,2,2,1)
     &                                     -2.d0*freestrain(m)
     &                 -dsqrt((strn(m,l,1,1,1)-strn(m,l,2,2,1))**2.d0
     &                      +(strn(m,l,1,2,1)+strn(m,l,2,1,1))**2.d0))
     &           +dabs(0.5d0*(strn(m,l,1,1,1)+strn(m,l,2,2,1)
     &                                     -2.d0*freestrain(m)
     &                 -dsqrt((strn(m,l,1,1,1)-strn(m,l,2,2,1))**2.d0
     &                      +(strn(m,l,1,2,1)+strn(m,l,2,1,1))**2.d0)))
     &                     )
          eminus(1,1)=0.5d0*(
     &                 0.5d0*(strn(m,l,1,1,1)+strn(m,l,2,2,1)
     &                                     -2.d0*freestrain(m)
     &                 +dsqrt((strn(m,l,1,1,1)-strn(m,l,2,2,1))**2.d0
     &                      +(strn(m,l,1,2,1)+strn(m,l,2,1,1))**2.d0))
     &           -dabs(0.5d0*(strn(m,l,1,1,1)+strn(m,l,2,2,1)
     &                                     -2.d0*freestrain(m)
     &                 +dsqrt((strn(m,l,1,1,1)-strn(m,l,2,2,1))**2.d0
     &                      +(strn(m,l,1,2,1)+strn(m,l,2,1,1))**2.d0)))
     &                      )
          eminus(2,2)=0.5d0*(
     &                 0.5d0*(strn(m,l,1,1,1)+strn(m,l,2,2,1)
     &                                     -2.d0*freestrain(m)
     &                 -dsqrt((strn(m,l,1,1,1)-strn(m,l,2,2,1))**2.d0
     &                      +(strn(m,l,1,2,1)+strn(m,l,2,1,1))**2.d0))
     &           -dabs(0.5d0*(strn(m,l,1,1,1)+strn(m,l,2,2,1)
     &                                     -2.d0*freestrain(m)
     &                 -dsqrt((strn(m,l,1,1,1)-strn(m,l,2,2,1))**2.d0
     &                      +(strn(m,l,1,2,1)+strn(m,l,2,1,1))**2.d0)))
     &                      )
          
          do 313 i=1,ndim
            do 323 j=1,ndim
              strnplus(i,j)=0.d0
              strnminus(i,j)=0.d0
              do 333 k=1,ndim
                do 343 q=1,ndim
                  strnplus(i,j)=strnplus(i,j)+
     &                          nxn(i,k)*eplus(k,q)*nxn(j,q)
                  strnminus(i,j)=strnminus(i,j)+
     &                           nxn(i,k)*eminus(k,q)*nxn(j,q)
 343            continue
 333          continue             
 323        continue
 313      continue

          if (.not.tc_split_on) then

          strs(m,l,1,1)=(c13*ekk+2.d0*c44*strn(m,l,1,1,1))*factor
     &      -p*((1.d0-phiF(1)*HeavisideJ1)*biot(m)
     &          +phiF(1)*HeavisideJ1)

          strs(m,l,2,2)=(c13*ekk+2.d0*c44*strn(m,l,2,2,1))*factor
     &      -p*((1.d0-phiF(1)*HeavisideJ1)*biot(m)
     &          +phiF(1)*HeavisideJ1) 

          strs(m,l,1,2)=(c44*(strn(m,l,1,2,1)+strn(m,l,2,1,1)))*factor
          strs(m,l,2,1)=strs(m,l,1,2)

          else

          strs(m,l,1,1)=(c13*ekkplus+2.d0*c44*strnplus(1,1))*factor
     &      +(c13*ekkminus+2.d0*c44*strnminus(1,1))
     &      -p*((1.d0-phiF(1)*HeavisideJ1)*biot(m)
     &          +phiF(1)*HeavisideJ1)

          strs(m,l,2,2)=(c13*ekkplus+2.d0*c44*strnplus(2,2))*factor
     &      +(c13*ekkminus+2.d0*c44*strnminus(2,2))
     &      -p*((1.d0-phiF(1)*HeavisideJ1)*biot(m)
     &          +phiF(1)*HeavisideJ1)

          strs(m,l,1,2)=(c44*(strnplus(1,2)+strnplus(2,1)))*factor
     &                 +(c44*(strnminus(1,2)+strnminus(2,1)))
          strs(m,l,2,1)=strs(m,l,1,2)

          endif
          
          do 210 i=1,ndim
            do 220 j=1,ndim
              tau(m,l,i,j)=(1.d0-mu)**2.d0*visc_eff(m,l)
     &                     *(vbar(i,j)+vbar(j,i))-p*delta(i,j)*1.d0
 220        continue
 210      continue

!          bphiJ(m,l)=-(phi(m,l,1)*Jacob(m,l,1)-phi(m,l,2)*Jacob(m,l,2))/dt
!     &            -0.5d0*(phiF(1)+phiF(2))*
!     &             ((1.d0-phi(m,l,1))-(1.d0-phi(m,l,2)))/dt
!     &            -(phiF(1)*(1.d0-phi(m,l,1))
!     &             *(dmax1(Jacob(m,l,1)-1.d0,0.d0))
!     &             -phiF(2)*(1.d0-phi(m,l,2))
!     &             *(dmax1(Jacob(m,l,2)-1.d0,0.d0)))/dt

          bphiJ(m,l)=-(phieff(1)*Jacob(m,l,1)*fdensity(1)
     &                -phieff(2)*Jacob(m,l,2)*fdensity(2))/dt

          bphiJ(m,l)=bphiJ(m,l)/(w0*fdensity(1))

          bv(m,l)=0.d0
          do 241 i=1,ndim
            bsf(m,l,i)=-(vint(m,l,i))*(1.d0-darcyF)
     &                 *Imperm(m,i,i)*kperm
            do 242 j=1,ndim
              bv(m,l)=bv(m,l)
     &        -(gfdensity(I)*vint(m,l,j)/fdensity(1)+gradv(m,l,j,I,1))
     &         *Jacob(m,l,1)*Finv(m,l,I,j,1)

              bsf(m,l,i)=bsf(m,l,i)
     &                   -Finv(m,l,J,i,1)*gradp(m,l,J)*0.d0
 242        continue
 241      continue

          if (.not.tc_split_on) then

          test=c44*(0.5d0*(strn(m,l,1,2,1)+strn(m,l,2,1,1))**2.d0+
     &              (strn(m,l,1,1,1)*strn(m,l,1,1,1))+
     &              (strn(m,l,2,2,1)*strn(m,l,2,2,1)))
     &         +c13*ekk*ekk/2.d0
     &         +HeavisideJ1*(p*ekk*(1.d0-biot(m))-0.5d0*p*p/(Cphi*e0))

          dHds(m,l,1,1)=0.d0
          dHds(m,l,2,2)=0.d0
          dHds(m,l,1,2)=0.d0
          dHds(m,l,2,1)=0.d0
          dHdp(m,l)=0.d0

          if (test.gt.Hpsi0(m,l).or.1.eq.1) then
            Hpsi(m,l)=test
            dHds(m,l,1,1)=(c13*ekk+2.d0*c44*strn(m,l,1,1,1))*dfactor
     &                    +(HeavisideJ1*p*(1.d0-biot(m)))*dfactor
            dHds(m,l,2,2)=(c13*ekk+2.d0*c44*strn(m,l,2,2,1))*dfactor
     &                    +(HeavisideJ1*p*(1.d0-biot(m)))*dfactor
            dHds(m,l,1,2)=(c44*strn(m,l,1,2,1)
     &                    +c44*strn(m,l,2,1,1))*dfactor
            dHds(m,l,2,1)=dHds(m,l,1,2)
            dHdp(m,l)=HeavisideJ1*((1.d0-biot(m))*ekk-p/(Cphi*e0))
          endif

          else

          test=c44*(0.5d0*(strnplus(1,2)+strnplus(2,1))**2.d0+
     &              (strnplus(1,1)*strnplus(1,1))+
     &              (strnplus(2,2)*strnplus(2,2)))
     &         +c13*ekkplus*ekkplus/2.d0
     &         +HeavisideJ1*(p*ekk*(1.d0-biot(m))-0.5d0*p*p/(Cphi*e0))

          dHds(m,l,1,1)=0.d0
          dHds(m,l,2,2)=0.d0
          dHds(m,l,1,2)=0.d0
          dHds(m,l,2,1)=0.d0
          dHdp(m,l)=0.d0

          if (test.gt.Hpsi0(m,l).or.1.eq.1) then
            Hpsi(m,l)=test
            dHds(m,l,1,1)=(c13*ekkplus+2.d0*c44*strnplus(1,1))*dfactor
     &                    +(HeavisideJ1*p*(1.d0-biot(m)))*dfactor
            dHds(m,l,2,2)=(c13*ekkplus+2.d0*c44*strnplus(2,2))*dfactor
     &                    +(HeavisideJ1*p*(1.d0-biot(m)))*dfactor
            dHds(m,l,1,2)=(c44*strnplus(1,2)+c44*strnplus(2,1))*dfactor
            dHds(m,l,2,1)=dHds(m,l,1,2)
            dHdp(m,l)=HeavisideJ1*((1.d0-biot(m))*ekk-p/(Cphi*e0))
          endif

          endif

          psimu(m,l)=-0.5d0*(1.d0-mu)+dfactor*(Hpsi(m,l))
!     &               -(p*ekk*HeavisideJ1*(1.d0-biot(m))
!     &                 -0.5d0*HeavisideJ1*p*p/(Cphi*e0))*DphiF(1)

          do 250 i=1,ndim
 250        ksi(m,l,i)=2.d0*gradmu(m,l,i)

          Wp(m,l)=0.d0
          Wmu(m,l)=0.25d0*(1.d0-mu)*(1.d0-mu)
          Ws(m,l)=(c44*(0.5d0*(strn(m,l,1,2,1)+strn(m,l,2,1,1))**2.d0+
     &              (strn(m,l,1,1,1)*strn(m,l,1,1,1))+
     &              (strn(m,l,2,2,1)*strn(m,l,2,2,1)))
     &         +c13*ekk*ekk/2.d0)*factor
          do 230 i=1,ndim
            Wp(m,l)=Wp(m,l)+p*gradmu(m,l,i)*uint(m,l,i)
            Wmu(m,l)=Wmu(m,l)+gradmu(m,l,i)*gradmu(m,l,i)
 230      continue        

 200    continue
 100  continue

      print*,H_total

      return
      end

**********************************************************************

      subroutine intfu(nints,ndim,ndpn,nnds,maxnd,nodes,maxele,
     &                 ele,gconu,det,B,NN,strs,muint,pint,fint,e0,
     &                 gradmu,local_eles,g_ele,maxprocs)
      implicit none
      real*8 det(maxele,nints),B(maxele,nnds,nints,ndim)
      real*8 NN(nnds,nints),fint(maxnd*ndpn),e0
      real*8 fele(nnds*ndim),strs(maxele,nints,ndim,ndim)
      real*8 muint(maxele,nints,2),pint(maxele,nints,2)
      real*8 gradmu(maxele,nints,ndim)
      integer maxnd,maxele,nodes,local_eles,g_ele(maxele)
      integer gconu(maxnd,ndpn),ele(maxele*maxprocs,nnds)
      integer i,j,k,l,m,ndof,nints,ndim,ndpn,nnds,maxprocs,ii,jj

      
      do 100 i=1,nodes*ndim
 100    fint(i)=0.d0

      do 200 i=1,local_eles

        do 210 j=1,nnds*ndim
 210      fele(j)=0.d0

        do 220 k=1,nints
          do 230 l=1,nnds
            m=ndim*l
            do 231 ii=1,ndim
              do 232 jj=1,ndim
                fele(m+ii-ndim)=fele(m+ii-ndim)
     &                +det(i,k)*B(i,l,k,jj)*strs(i,k,jj,ii)
 232          continue
 231        continue
 230      continue
 220    continue

        do 300 j=1,nnds
          do 310 k=1,ndim
            ndof=gconu(ele(g_ele(i),j),k)
            fint(ndof)=fint(ndof)+fele(ndim*(j-1)+k)
 310      continue
 300    continue

 200  continue


      return
      end

**********************************************************************
      
      subroutine intfmu(nnds,nints,ndim,ndpn,maxnd,nodes,maxele,
     &          ele,gconmu,det,B,NN,psimu,ksi,fint,dt,beta,muint,
     &          local_eles,g_ele,maxprocs)
      implicit none
      real*8 det(maxele,nints),B(maxele,nnds,nints,ndim)
      real*8 NN(nnds,nints),psimu(maxele,nints)
      real*8 fele(nnds),dt,beta,muint(maxele,nints,2)
      real*8 ksi(maxele,nints,ndim),fint(maxnd*ndpn)
      integer maxnd,maxele,nodes,local_eles,g_ele(maxele)
      integer gconmu(maxnd,ndpn),ele(maxele*maxprocs,nnds)
      integer i,j,k,l,ndof,nnds,nints,ndim,ndpn,maxprocs,ii

      
      do 100 i=1,nodes
 100    fint(i)=0.d0

      do 200 i=1,local_eles

        do 210 j=1,nnds
 210      fele(j)=0.d0

        do 220 k=1,nints
          do 230 l=1,nnds
            fele(l)=fele(l)+det(i,k)*NN(l,k)*psimu(i,k)+
     &               (0.d0*beta/dt)*det(i,k)*NN(l,k)*muint(i,k,1)
            do 231 ii=1,ndim
              fele(l)=fele(l)+det(i,k)*B(i,l,k,ii)*ksi(i,k,ii)
 231        continue
 230      continue
 220    continue

        do 300 j=1,nnds
          ndof=gconmu(ele(g_ele(i),j),1)
          fint(ndof)=fint(ndof)+fele(j)
 300    continue

 200  continue


      return
      end
      
**********************************************************************
      
      subroutine intfvp(nnds,nints,ndim,ndpn,maxnd,nodes,maxele,
     &                  ele,gconvp,det,B,NN,bv,bp,tau,fint,dt,muint,
     &                  mu,bsf,Finv,Jacob,gradv,phi,w0,bphiJ,
     &                  timerule,local_eles,g_ele,maxprocs)
      implicit none
      real*8 det(maxele,nints),B(maxele,nnds,nints,ndim)
      real*8 NN(nnds,nints),bv(maxele,nints),bphiJ(maxele,nints)
      real*8 bp(maxele,nints,ndim),tau(maxele,nints,ndim,ndim)
      real*8 fele(nnds*(ndim+1)),dt,muint(maxele,nints,2)
      real*8 fint(maxnd*ndpn),bsf(maxele,nints,ndim)
      real*8 mu(maxnd),Jacob(maxele,nints,2)
      real*8 Finv(maxele,nints,ndim,ndim,2),w0,timerule
      real*8 gradv(maxele,nints,ndim,ndim,2),phi(maxele,nints,2)
      integer maxnd,maxele,nodes,local_eles,g_ele(maxele)
      integer gconvp(maxnd,ndpn),ele(maxele*maxprocs,nnds)
      integer i,j,k,l,ndof,nnds,nints,ndim,ndpn,maxprocs,m,ii,jj,kk

      
      do 100 i=1,nodes*(ndim+1)
 100    fint(i)=0.d0

      do 200 i=1,local_eles

        do 210 j=1,nnds*(ndim+1)
 210      fele(j)=0.d0

        do 220 k=1,nints
 
          do 230 l=1,nnds
            m=(ndim+1)*l

            do 231 kk=1,ndim
              do 232 ii=1,ndim
                do 233 jj=1,ndim
                  fele(m-(ndim+1)+kk)=fele(m-(ndim+1)+kk)
     &             +det(i,k)*B(i,l,k,II)*Jacob(i,k,1)
     &              *tau(i,k,jj,kk)*Finv(i,k,II,jj,1)*timerule
 233            continue
 232          continue
              fele(m-(ndim+1)+kk)=fele(m-(ndim+1)+kk)
     &                            +det(i,k)*NN(l,k)*Jacob(i,k,1)
     &                            *(-bsf(i,k,kk))*timerule
 231        continue

            fele(m)=fele(m)+det(i,k)*NN(l,k)
     &                      *(bv(i,k)*timerule+bphiJ(i,k))
 230      continue
 220    continue

        do 300 j=1,nnds
          do 310 k=1,ndim+1
            ndof=gconvp(ele(g_ele(i),j),k)
            fint(ndof)=fint(ndof)+fele(3*(j-1)+k)
 310      continue
 300    continue

 200  continue


      return
      end

**********************************************************************
**********************************************************************

      subroutine stf(nnds,nints,ndim,B,NN,det,node,maxele,maxnd,
     &            maxcon,maxdbc,nvpdofs,ele,dt,Cstiff,
     &            gcon,mpc,mpcnum,conA,nvpdbcs,nmudofs,nmudbcs,
     &            nudofs,nudbcs,muint,Jacob,Finv,Imperm,w0,
     &            gradv,pint,tau,PetA,e0,dstrn,bsf,vint,strn,
     &            Hpsi,dJacob,dphi,phi,beta,a0,isu,isvp,ismu,
     &            PetR,PetdU,ksp,nodes,indices,gradp,dHdp,
     &            visc_eff,dHds,local_eles,p0,Cf,fluidsource,
     &            biot,phi0,freestrain,YM,poisson,KDtN,maxbdnd,
     &            DtNdofs,DtNidx,tc_split_on,muavg,open_thresh,
     &            kperm,rank,eleproc,ndpn,artstiff,dmuint,timerule,
     &            gmuint,stretch,g_ele,maxprocs)
      implicit none
#include <petsc/finclude/petsc.h>

      real*8 B(maxele,nnds,nints,ndim),det(maxele,nints)
      real*8 bele(nnds,nints,ndim),detele(nints)
      real*8 kele(nnds*(2*ndim+2),nnds*(2*ndim+2))
      real*8 Imperm(maxele,ndim,ndim)
      real*8 w0,dt,Cstiff(4),biot(maxele)
      real*8 gradv(maxele,nints,ndim,ndim,2),e0,beta
      real*8 phi0(maxele)
      real*8 node(maxnd,ndim),conA(maxdbc,maxcon)
      real*8 NN(nnds,nints),kcon,muint(maxele,nints,2)
      real*8 Jacob(maxele,nints,2),Finv(maxele,nints,ndim,ndim,2)
      real*8 pint(maxele,nints,2),tau(maxele,nints,ndim,ndim)
      real*8 dstrn(maxele,nints,ndim,ndim)
      real*8 bsf(maxele,nints,ndim)
      real*8 vint(maxele,nints,ndim),gradp(maxele,nints,ndim)
      real*8 strn(maxele,nints,ndim,ndim,2)
      real*8 Hpsi(maxele,nints),dJacob(maxele,nints)
      real*8 dphi(maxele,nints),a0,phi(maxele,nints,2)
      real*8 Petkele(nnds*ndpn*(maxcon+1),
     &               nnds*ndpn*(maxcon+1))
      real*8 KDtN(maxbdnd*ndim,maxbdnd*ndim)
      real*8 visc_eff(maxele,nints),dHds(maxele,nints,ndim,ndim)
      real*8 freestrain(maxele),YM(maxele)
      real*8 poisson(maxele),p0,Cf,fluidsource(maxele)
      real*8 muavg(maxele),open_thresh,kperm
      real*8 dmuint(maxele,nints),timerule,dHdp(maxele,nints)
      real*8 gmuint(maxele,nints,ndim),stretch(maxele,nints)
      integer maxele,maxnd,maxcon,maxdbc,local_eles,maxbdnd,ndpn
      integer ele(maxele*maxprocs,nnds),gcon(maxnd,ndpn)
      integer i,j,k,l,n,nvpdofs,nvpdbcs,nudofs,nudbcs
      integer dofi,dofj,kglb,nloci,nlocj,nmudofs,nmudbcs
      integer mpc(maxdbc,maxcon,2),mpcnum(maxdbc)
      integer ii,jj,coni,conj,nci,ncj,nnds,ndim,nints
      integer idx(nnds*ndpn*(maxcon+1)),Peti,Petj
      integer isdof,isu,isvp,ismu,ndofs,zeroRowcount
      integer zeroRow(maxnd*ndpn),nodes,maxprocs
      integer indices(maxnd*ndpn*ndim),g_ele(maxele)
      integer DtNidx(maxbdnd*ndim),DtNdofs,rank
      logical tc_split_on,eleproc(maxele),artstiff

      Mat PetA
      Vec PetR
      Vec PetRs
      Vec PetdU
      KSP ksp

      PetscErrorCode ierr
      PetscViewer Aview


      call MatZeroEntries(PetA,ierr)
      ndofs=nmudofs+nudofs+nvpdofs

      do 100 l=1,local_eles
 
        do 110 i=1,nnds
          do 111 k=1,nints
            do 112 j=1,ndim
              bele(i,k,j)=B(l,i,k,j)
 112        continue
 111      continue
 110    continue

        do 120 i=1,nints
          detele(i)=det(l,i)
 120    continue

        do 121 i=1,nnds*(maxcon+1)*(2*ndim+2)
          idx(i)=-1
          do 122 j=1,nnds*(maxcon+1)*(2*ndim+2)
            Petkele(i,j)=0.d0
 122      continue
 121    continue

        call stfe(nints,nnds,ndim,dt,maxele,l,kele,
     &                detele,bele,NN,pint,Cstiff,Imperm,Finv,
     &                Jacob,gradv,muint,w0,tau,e0,dstrn,bsf,
     &                vint,dJacob,strn,Hpsi,dphi,phi,beta,a0,
     &                visc_eff,dHds,g_ele,maxprocs,gradp,
     &                biot,phi0,freestrain,YM,p0,Cf,dHdp,
     &                poisson,tc_split_on,muavg,open_thresh,
     &                kperm,artstiff,dmuint,timerule,gmuint,
     &                stretch,fluidsource)

        do 130 i=1,nnds*(2*ndim+2)
          do 140 j=1,nnds*(2*ndim+2)

            nloci=(i-1)/(2*ndim+2)+1
            dofi=mod(i-1,2*ndim+2)+1
            nlocj=(j-1)/(2*ndim+2)+1
            dofj=mod(j-1,2*ndim+2)+1

            if (dofi.le.2) then
              if (dofj.le.2) then
                isdof=isu*isu
              elseif (dofj.le.5) then
                isdof=isu*isvp
              else
                isdof=isu*ismu
              endif
            elseif (dofi.le.5) then
              if (dofj.le.2) then
                isdof=isvp*isu
              elseif (dofj.le.5) then
                isdof=isvp*isvp
              else
                isdof=isvp*ismu
              endif
            else
              if (dofj.le.2) then
                isdof=ismu*isu
              elseif (dofj.le.5) then
                isdof=ismu*isvp
              else
                isdof=ismu*ismu
              endif
            endif

            dofi=gcon(ele(g_ele(l),nloci),dofi)
            dofj=gcon(ele(g_ele(l),nlocj),dofj)
            coni=dofi-nvpdofs-nvpdbcs-nudofs-nudbcs-nmudofs-nmudbcs
            conj=dofj-nvpdofs-nvpdbcs-nudofs-nudbcs-nmudofs-nmudbcs
            
            Peti=(i-1)*(maxcon+1)+1
            Petj=(j-1)*(maxcon+1)+1

            if (coni.gt.0) then
              nci=mpcnum(coni)
              do 150 ii=1,nci
                dofi=gcon(mpc(coni,ii,1),mpc(coni,ii,2))
                if (dofi.gt.ndofs) goto 150
                if(conj.gt.0) then
                  ncj=mpcnum(conj)              
                  do 155 jj=1,ncj
                    dofj=gcon(mpc(conj,jj,1),
     &                          mpc(conj,jj,2))
                    if (dofj.gt.ndofs) goto 155

                    Petkele(Petj+jj,Peti+ii)=kele(i,j)
     &                              *conA(coni,ii)*conA(conj,jj)
     &                *dble(isdof)

                    if (Peti+ii.eq.Petj+jj) then
                      idx(Peti+ii)=dofi-1
                    endif

 155              continue
                else
                  if (dofj.gt.ndofs) goto 150

                  Petkele(Petj,Peti+ii)=kele(i,j)*conA(coni,ii)
     &              *dble(isdof)

                endif
 150          continue
            elseif(conj.gt.0) then
              if (dofi.gt.ndofs) goto 140
              ncj=mpcnum(conj)
              do 160 jj=1,ncj
                dofj=gcon(mpc(conj,jj,1),mpc(conj,jj,2))
                if (dofj.gt.ndofs) goto 160

                Petkele(Petj+jj,Peti)=kele(i,j)*conA(conj,jj)
     &            *dble(isdof)

 160          continue
            else
              if (dofi.gt.ndofs) goto 140
              if (dofj.gt.ndofs) goto 140


              Petkele(Petj,Peti)=kele(i,j)
     &          *dble(isdof)

              if (Peti.eq.Petj) then
                idx(Peti)=dofi-1
                Petkele(Petj,Peti)=Petkele(Petj,Peti)
     &                            +dble(1-isdof)
              endif

            endif
 140      continue
 130    continue

        call MatSetValues(PetA,nnds*(2*ndim+2)*(maxcon+1),idx,
     &                         nnds*(2*ndim+2)*(maxcon+1),idx,
     &                    Petkele,ADD_VALUES,ierr)
      CHKERRQ(ierr)

 100  continue

      if (isu.eq.1) then
        if (rank.eq.0)
     &  call MatSetValues(PetA,maxbdnd*ndim,DtNidx,
     &                         maxbdnd*ndim,DtNidx,
     &                    KDtN,ADD_VALUES,ierr)
      endif

      call MatAssemblyBegin(PetA,MAT_FINAL_ASSEMBLY,ierr)
      call MatAssemblyEnd(PetA,MAT_FINAL_ASSEMBLY,ierr)
      call PetscViewerASCIIOpen(PETSC_COMM_WORLD,"mat.output",
     &                  Aview,ierr)
!      call PetscViewerSetFormat(Aview,PETSC_VIEWER_ASCII_DENSE
!     &                          ,ierr)
!      call MatView(PetA,Aview,ierr)

      zeroRowcount=0
      do 570 i=1,nodes
        do 571 j=1,ndim
          if (gcon(i,j).le.ndofs.and.isu.eq.0) then
            zeroRowcount=zeroRowcount+1
            zeroRow(zeroRowcount)=indices(gcon(i,j))
          endif
 571    continue
        do 572 j=1,ndim+1
          if (gcon(i,j+2).le.ndofs.and.isvp.eq.0) then
            zeroRowcount=zeroRowcount+1
            zeroRow(zeroRowcount)=indices(gcon(i,j+2))
          endif
 572    continue
        if (gcon(i,6).le.ndofs.and.ismu.eq.0) then
          zeroRowcount=zeroRowcount+1
          zeroRow(zeroRowcount)=indices(gcon(i,6))
        endif
 570  continue

!            do 1500 i=1,zeroRowcount
! 1500         print*,i,zeroRow(i)
    
!      call MatSetOption(PetA,MAT_KEEP_NONZERO_PATTERN,
!     &                  PETSC_TRUE,ierr)

      call MatZeroRows(PetA,zeroRowcount,zeroRow,
     &                 1.d0,PetR,PetdU,ierr)

      call MatSetOption(PetA,MAT_NEW_NONZERO_LOCATIONS,
     &                  PETSC_FALSE,ierr)

      call KSPSetOperators(ksp,PetA,PetA,ierr)

      call PetscViewerDestroy(Aview,ierr)
     
      return
      end

**********************************************************************

      subroutine stfe(nints,nnds,ndim,dt,maxele,iele,kele,
     &                detele,bele,nele,pint,Cstiff,imperm,Finv,
     &                Jacob,gradv,muint,w0,tau,e0,dstrn,bsf,
     &                vint,dJacob,strn,Hpsi,dphi,phi,beta,a0,
     &                visc_eff,dHds,g_ele,maxprocs,gradp,
     &                biot,phi0,freestrain,YM,p0,Cf,dHdp,
     &                poisson,tc_split_on,muavg,open_thresh,
     &                kperm,artstiff,dmuint,timerule,gmuint,
     &                stretch,fluidsource)
      implicit none
      real*8 kele(nnds*(2*ndim+2),nnds*(2*ndim+2))
      real*8 detele(nints),bele(nnds,nints,ndim)
      real*8 nele(nnds,nints),muint(maxele,nints,2),emu
      real*8 gradv(maxele,nints,ndim,ndim,2),egradv(ndim,ndim)
      real*8 dt,Cstiff(4),imperm(maxele,ndim,ndim),w0,e0,beta
      real*8 biot(maxele),ebiot,dHdp(maxele,nints)
      real*8 freestrain(maxele),estrn0
      real*8 Jacob(maxele,nints,2),eJac
      real*8 ebsf(ndim),phi0(maxele)
      real*8 pint(maxele,nints,2),ep,bsf(maxele,nints,ndim),ephi0
      real*8 Finv(maxele,nints,ndim,ndim,2),eFinv(ndim,ndim)
      real*8 visc_eff(maxele,nints),evisc_eff,eopen,eopena
      real*8 etau(ndim,ndim),bf(ndim),tau(maxele,nints,ndim,ndim)
      real*8 dstrn(maxele,nints,ndim,ndim),edstrn(ndim,ndim)
      real*8 vint(maxele,nints,ndim),gradp(maxele,nints,ndim)
      real*8 strn(maxele,nints,ndim,ndim,2)
      real*8 Hpsi(maxele,nints),dJacob(maxele,nints)
      real*8 dHds(maxele,nints,ndim,ndim),edHds(ndim,ndim)
      real*8 dphi(maxele,nints),edmu,phi(maxele,nints,2),edphi
      real*8 v(ndim),estrn(ndim,ndim),eHpsi,edJac,ephi,ephib
      real*8 YM(maxele),poisson(maxele),egradp(ndim)
      real*8 eYM,epoisson,kperm,ep_old,p0,Cf,fluidsource(maxele)
      real*8 dmuint(maxele,nints),gmuint(maxele,nints,ndim)
      real*8 egmu(ndim),stretch(maxele,nints),estretch
      real*8 uBBv(ndim,ndim,ndim,ndim),vBBv(ndim,ndim,ndim,ndim)
      real*8 pNBv(ndim,ndim),uBNv(ndim,ndim,ndim),pBNv(ndim,ndim)
      real*8 vNNv(ndim,ndim),uBBu(ndim,ndim,ndim,ndim)
      real*8 pNBu(ndim,ndim),uBNp(ndim,ndim),vBNp(ndim,ndim)
      real*8 oNNv(ndim),oNBv(ndim,ndim),oNBu(ndim,ndim),pNNp
      real*8 oBBo(ndim,ndim),uBNo(ndim,ndim),oNNo,pNNo,oNNp,a0
      real*8 pBNp(ndim),vNNp(ndim),edHdp
      real*8 muavg(maxele),emuavg,open_thresh,eimperm(ndim,ndim)
      real*8 ephi_old,eJac_old,emu_old,timerule,einject
      integer maxele,nints,nnds,ndim,g_ele(maxele),maxprocs
      integer iele,i,ia,k,l,i1,j1,j,Peti,Petj
      logical tc_split_on,artstiff


      bf(1)=0.d0
      bf(2)=0.d0

      do 100 i=1,nnds*(2*ndim+2)
        do 101 j=1,nnds*(2*ndim+2)
          kele(i,j)=0.d0
 101    continue
 100  continue
      
      do 200 ia=1,nints

        emu=muint(iele,ia,1)
        emu_old=muint(iele,ia,2)
        edphi=dphi(iele,ia)
        eJac=Jacob(iele,ia,1)
        eJac_old=Jacob(iele,ia,2)
        edJac=dJacob(iele,ia)
        ep=pint(iele,ia,1)
        ep_old=pint(iele,ia,2)
        eHpsi=Hpsi(iele,ia)
        ephi=phi(iele,ia,1)
        ephi_old=phi(iele,ia,2)
        edHdp=dHdp(iele,ia)
        evisc_eff=visc_eff(iele,ia)
        ephi0=phi0(iele)
        ebiot=biot(iele)
        estrn0=freestrain(iele)
        eYM=YM(iele)
        epoisson=poisson(iele)
        emuavg=muavg(iele)
        edmu=dmuint(iele,ia)
        estretch=stretch(iele,ia)
        einject=fluidsource(iele)

        do 300 i=1,ndim
          ebsf(i)=bsf(iele,ia,i)
          v(i)=vint(iele,ia,i)
          egmu(i)=gmuint(iele,ia,i)
          egradp(i)=gradp(iele,ia,i)
          do 301 j=1,ndim
            eFinv(i,j)=Finv(iele,ia,i,j,1)
            egradv(i,j)=gradv(iele,ia,i,j,1)
            etau(i,j)=tau(iele,ia,i,j)
            edstrn(i,j)=dstrn(iele,ia,i,j)
            estrn(i,j)=strn(iele,ia,i,j,1)
            edHds(i,j)=dHds(iele,ia,i,j)
            if (ia.eq.1) eimperm(i,j)=imperm(iele,i,j)
 301      continue
 300    continue            

        call matT(ndim,e0,emu,eFinv,eJac,evisc_eff,eimperm,Cstiff,
     &           etau,ep,ebsf,dt,v,egradv,bf,edstrn,estrn,edJac,eHpsi,
     &           edphi,beta,ephi,edHds,ep_old,ebiot,p0,Cf,einject,
     &           ephi0,uBBv,vBBv,pNBv,uBNv,vNNv,uBBu,pNBu,uBNp,
     &           vBNp,pBNv,oNNv,oNBv,oNBu,oBBo,uBNo,oNNo,pNNo,oNNp,
     &           pNNp,pBNp,vNNp,a0,estrn0,eYM,epoisson,tc_split_on,
     &           emuavg,open_thresh,kperm,artstiff,edmu,ephi_old,
     &           eJac_old,emu_old,w0,timerule,egmu,estretch,egradp,
     &           edHdp)
     
        do 220 i1=1,nnds
          do 230 j1=1,nnds
            
            do 240 i=1,ndim
              Peti=(i1-1)*(2*ndim+2)+i
              do 250 j=1,ndim
                Petj=(j1-1)*(2*ndim+2)+j
                do 251 k=1,ndim
                do 252 l=1,ndim
                  kele(Peti,Petj)=kele(Peti,Petj)
     &              +bele(j1,ia,k)*uBBu(j,k,l,i)
     &              *bele(i1,ia,l)*detele(ia)
 252            continue
 251            continue
 250          continue
              do 260 j=1,ndim
                Petj=(j1-1)*(2*ndim+2)+ndim+j
 260          continue
              Petj=(j1-1)*(2*ndim+2)+2*ndim+1
              do 263 k=1,ndim
                kele(Peti,Petj)=kele(Peti,Petj)
     &            +nele(j1,ia)*pNBu(k,i)
     &            *bele(i1,ia,k)*detele(ia)
 263          continue
              Petj=(j1-1)*(2*ndim+2)+2*ndim+2
              do 264 k=1,ndim
                kele(Peti,Petj)=kele(Peti,Petj)
     &            +nele(j1,ia)*oNBu(k,i)
     &            *bele(i1,ia,k)*detele(ia)
 264          continue
 240        continue

            do 270 i=1,ndim
              Peti=(i1-1)*(2*ndim+2)+ndim+i
              do 271 j=1,ndim
                Petj=(j1-1)*(2*ndim+2)+j
                do 272 k=1,ndim
                do 273 l=1,ndim
                  kele(Peti,Petj)=kele(Peti,Petj)
     &              +bele(j1,ia,k)*uBBv(j,k,l,i)
     &              *bele(i1,ia,l)*detele(ia)
 273            continue
                  kele(Peti,Petj)=kele(Peti,Petj)
     &              +bele(j1,ia,k)*uBNv(j,k,i)
     &              *nele(i1,ia)*detele(ia)
 272            continue
 271          continue
              do 280 j=1,ndim
                Petj=(j1-1)*(2*ndim+2)+ndim+j
                do 281 k=1,ndim
                do 282 l=1,ndim
                  kele(Peti,Petj)=kele(Peti,Petj)
     &              +bele(j1,ia,k)*vBBv(j,k,l,i)
     &              *bele(i1,ia,l)*detele(ia)
 282            continue
 281            continue
                kele(Peti,Petj)=kele(Peti,Petj)
     &            +nele(j1,ia)*vNNv(j,i)
     &            *nele(i1,ia)*detele(ia)
 280          continue
              Petj=(j1-1)*(2*ndim+2)+2*ndim+1
              do 283 k=1,ndim
                kele(Peti,Petj)=kele(Peti,Petj)
     &            +nele(j1,ia)*pNBv(k,i)
     &            *bele(i1,ia,k)*detele(ia)
                kele(Peti,Petj)=kele(Peti,Petj)
     &            +bele(j1,ia,k)*pBNv(k,i)
     &            *nele(i1,ia)*detele(ia)
 283          continue
              Petj=(j1-1)*(2*ndim+2)+2*ndim+2
              kele(Peti,Petj)=kele(Peti,Petj)
     &          +nele(j1,ia)*oNNv(i)
     &          *nele(i1,ia)*detele(ia)
              do 284 k=1,ndim
                kele(Peti,Petj)=kele(Peti,Petj)
     &            +nele(j1,ia)*oNBv(k,i)
     &            *bele(i1,ia,k)*detele(ia)
 284          continue
 270        continue

            Peti=(i1-1)*(2*ndim+2)+2*ndim+1
            do 290 j=1,ndim
              Petj=(j1-1)*(2*ndim+2)+j
              do 291 k=1,ndim
                kele(Peti,Petj)=kele(Peti,Petj)
     &            +bele(j1,ia,k)*uBNp(j,k)
     &            *nele(i1,ia)*detele(ia)
 291          continue
 290        continue
            do 292 j=1,ndim
              Petj=(j1-1)*(2*ndim+2)+ndim+j
              do 293 k=1,ndim
                kele(Peti,Petj)=kele(Peti,Petj)
     &            +bele(j1,ia,k)*vBNp(j,k)
     &            *nele(i1,ia)*detele(ia)
 293          continue
              kele(Peti,Petj)=kele(Peti,Petj)
     &          +nele(j1,ia)*vNNp(j)
     &          *nele(i1,ia)*detele(ia)
 292        continue
            Petj=(j1-1)*(2*ndim+2)+2*ndim+1
            kele(Peti,Petj)=kele(Peti,Petj)
     &        +nele(j1,ia)*pNNp
     &        *nele(i1,ia)*detele(ia)
            do 285 k=1,ndim
              kele(Peti,Petj)=kele(Peti,Petj)
     &          +bele(j1,ia,k)*pBNp(k)
     &          *nele(i1,ia)*detele(ia)
 285        continue
            Petj=(j1-1)*(2*ndim+2)+2*ndim+2
            kele(Peti,Petj)=kele(Peti,Petj)
     &        +nele(j1,ia)*oNNp
     &        *nele(i1,ia)*detele(ia)

            Peti=(i1-1)*(2*ndim+2)+2*ndim+2
            do 294 j=1,ndim
              Petj=(j1-1)*(2*ndim+2)+j
              do 295 k=1,ndim
                kele(Peti,Petj)=kele(Peti,Petj)
     &            +bele(j1,ia,k)*uBNo(j,k)
     &            *nele(i1,ia)*detele(ia)
 295          continue
 294        continue
            do 296 j=1,ndim
              Petj=(j1-1)*(2*ndim+2)+ndim+j
 296        continue
            Petj=(j1-1)*(2*ndim+2)+2*ndim+1
            kele(Peti,Petj)=kele(Peti,Petj)
     &        +nele(j1,ia)*pNNo
     &        *nele(i1,ia)*detele(ia)
            Petj=(j1-1)*(2*ndim+2)+2*ndim+2
            kele(Peti,Petj)=kele(Peti,Petj)
     &        +nele(j1,ia)*oNNo
     &        *nele(i1,ia)*detele(ia)
            do 297 k=1,ndim
            do 298 l=1,ndim
              kele(Peti,Petj)=kele(Peti,Petj)
     &          +bele(j1,ia,k)*oBBo(k,l)
     &          *bele(i1,ia,l)*detele(ia)
 298        continue
 297        continue

 230      continue
 220    continue
 200  continue


      return
      end

**********************************************************************


      subroutine matT(ndim,e0,mu,Finv,Jac,visc_eff,imperm,Cstiff,tau,
     &               p,bsf,dt,v,gradv,bf,dstrn,strn,dJac,Hpsi,dphi,
     &               beta,phi,dHds,p_old,biot,p0,Cf,einject,
     &               phi0,uBBv,vBBv,pNBv,uBNv,vNNv,uBBu,pNBu,uBNp,vBNp,
     &               pBNv,oNNv,oNBv,oNBu,oBBo,uBNo,oNNo,pNNo,oNNp,pNNp,
     &               pBNp,vNNp,a,strn0,YM,poisson,tc_split_on,muavg,
     &               open_thresh,kperm,artstiff,dmu,phi_old,Jac_old,
     &               mu_old,w0,timerule,gmu,stretch,gradp,dHdp)
      implicit none
      real*8 Cstiff(4),w0,imperm(ndim,ndim),mu,p,tau(ndim,ndim)
      real*8 Jac,dJac,c11,c33,c13,c44,bf(ndim),gradv(ndim,ndim)
      real*8 Finv(ndim,ndim),bsf(ndim),dt,e0,v(ndim),strn(ndim,ndim)
      real*8 vbar(ndim,ndim),vkk,Hpsi,dphi,beta,phi,YM,poisson
      real*8 dJFdu(ndim,ndim,ndim,ndim),dtaudv(ndim,ndim,ndim,ndim)
      real*8 dtaudu(ndim,ndim,ndim,ndim),dstrn(ndim,ndim)
      real*8 dbsfdv(ndim,ndim),dstrsdstrn(2,2,2,2)
      real*8 dstrsdp(ndim,ndim),delta(ndim,ndim),factor,dfactor,a
      real*8 dtaudmu(ndim,ndim),dbsfdmu(ndim),ddfactor
      real*8 biot,Bc,phicrack,Ba,dFdu(ndim,ndim,ndim,ndim)
      real*8 dstrsdmu(ndim,ndim),dHds(ndim,ndim),gradp(ndim)
      real*8 uBBv(ndim,ndim,ndim,ndim),vBBv(ndim,ndim,ndim,ndim)
      real*8 pNBv(ndim,ndim),uBNv(ndim,ndim,ndim),pBNv(ndim,ndim)
      real*8 vNNv(ndim,ndim),uBBu(ndim,ndim,ndim,ndim)
      real*8 pNBu(ndim,ndim),uBNp(ndim,ndim),vBNp(ndim,ndim)
      real*8 oNNv(ndim),oNBv(ndim,ndim),oNBu(ndim,ndim)
      real*8 oBBo(ndim,ndim),uBNo(ndim,ndim),oNNo,pNNo,oNNp,pNNp
      real*8 pBNp(ndim),vNNp(ndim),dstretchdJ
      real*8 Cphi,phi0,phiterm,ekk,dkkplus,dkkminus,dHdp
      real*8 ekkplus,ekkminus,z,etotal(ndim),dphidstrain
      real*8 eplus(ndim,ndim),eminus(ndim,ndim),theta,pi
      real*8 strnplus(ndim,ndim),strnminus(ndim,ndim),nxn(ndim,ndim)
      real*8 h0,strn0,open_thresh,kperm,p_old,p0,Cf,einject
      real*8 aexp,bexp,phiF(2),darcyF,DphiF(2),DdarcyF,muavg,dmu
      real*8 dJdu(ndim,ndim),DDphiF(2),dphibdmu,dphibdp,beta2
      real*8 dphibdJ,dphieffdmu,dphieffdp,dphieffdJ
      real*8 phieff(2),phi_old,Jac_old,mu_old,visc_eff,timerule
      real*8 crit_strain,HeavisideJ1,HeavisideJ2,dHJ1,istr
      real*8 stretch,vscale_factor,dstretchdF(ndim,ndim),gmu(ndim)
      real*8 s_numerator,s_denominator,crit_stretch,fdensity(2)
      integer i,j,k,l,m,n,ndim,q
      logical tc_split_on,no_tc_split,artstiff

      pi=4.d0*datan(1.d0)

      bf(1)=0.d0
      bf(2)=0.d0

      delta(1,1)=1.d0
      delta(1,2)=0.d0
      delta(2,1)=0.d0
      delta(2,2)=1.d0

      c11=YM*(1.d0-poisson)
     &    /((1.d0+poisson)*(1.d0-2.d0*poisson))
      c33=c11
      c13=YM*poisson
     &    /((1.d0+poisson)*(1.d0-2.d0*poisson))
      c44=YM
     &    /(2.d0*(1.d0+poisson))

      if (dabs(biot-1.d0).lt.1.d-15) then
        Cphi=1.d20
      elseif (dabs(biot-phi0).lt.1.d-15) then
        Cphi=1.d20
      else 
        Cphi=(YM/(e0*3.d0*(1.d0-2.d0*poisson)))
     &       /((biot-phi0)*(1.d0-biot))
      endif

      h0=open_thresh

      if (a.le.0.d0) then
        factor=-mu*mu*(2.d0*mu-3.d0+a*(mu-1.d0))
        dfactor=(6.d0*(1.d0-mu)+a*(2.d0-3.d0*mu))*mu
        ddfactor=2.d0*(3.d0+a-6.d0*mu-3.d0*a*mu)
      else
        factor=(1.d0+a)*(1.d0-(a/(1.d0+a))**(mu*mu))
        dfactor=2.d0*(factor-(1.d0+a))*mu*dlog(a/(1.d0+a))
        ddfactor=2.d0*(factor-(1.d0+a))*dlog(a/(1.d0+a))*
     &        (1.d0+2.d0*mu*mu*dlog(a/(1.d0+a)))
      endif

!      if (factor.lt.1.d-4.and.artstiff) factor=1.d-2

!      aexp=1.d0
!      bexp=2.d0
!      phiF(1)=(1.d0-mu**aexp)**bexp
!      phiF(2)=(1.d0-mu_old**aexp)**bexp
!
!      DphiF(1)=-aexp*bexp*(mu**(aexp-1.d0))
!     &     *(1.d0-mu**aexp)**(bexp-1.d0)
!      DphiF(2)=-aexp*bexp*(mu_old**(aexp-1.d0))
!     &     *(1.d0-mu_old**aexp)**(bexp-1.d0)
!
!      if (aexp.lt.1.1d0.and.bexp.lt.1.1d0) then
!        DDphiF(1)=0.d0
!        DDphiF(2)=0.d0
!      elseif (aexp.lt.1.1d0) then
!        DDphiF(1)=(bexp-1.d0)*bexp*(1.d0-mu)**(bexp-2.d0)
!        DDphiF(2)=(bexp-1.d0)*bexp*(1.d0-mu_old)**(bexp-2.d0)
!      elseif (bexp.lt.1.1d0) then
!        DDphiF(1)=-(aexp-1.d0)*aexp*mu**(aexp-2.d0)
!        DDphiF(2)=-(aexp-1.d0)*aexp*mu_old**(aexp-2.d0)
!      else
!        DDphiF(1)=aexp*bexp*mu**(aexp-2.d0)
!     &      *(1.d0-mu**aexp)**(bexp-2.d0)
!     &      *(1.d0-mu**aexp-aexp*(1.d0-bexp*mu**aexp))
!        DDphiF(2)=aexp*bexp*mu_old**(aexp-2.d0)
!     &      *(1.d0-mu_old**aexp)**(bexp-2.d0)
!     &      *(1.d0-mu_old**aexp-aexp*(1.d0-bexp*mu_old**aexp))
!      endif

      phiF(1)=1.d0-3.d0*mu**2.d0+2.d0*mu**3.d0
      phiF(2)=1.d0-3.d0*mu_old**2.d0+2.d0*mu_old**3.d0
      DphiF(1)=-6.d0*mu*(1.d0-mu)
      DphiF(2)=-6.d0*mu_old*(1.d0-mu_old)
      DDphiF(1)=-6.d0+12.d0*mu
      DDphiF(2)=-6.d0+12.d0*mu_old

      aexp=4.d0
      bexp=4.d0
      darcyF=(1.d0-muavg**aexp)**bexp+0.d-6
      DdarcyF=-aexp*bexp*(muavg**(aexp-1.d0))
     &       *(1.d0-muavg**aexp)**(bexp-1.d0)
!      darcyF=(1.d0-mu**aexp)**bexp+0.d-6
!      DdarcyF=-aexp*bexp*(mu**(aexp-1.d0))
!     &       *(1.d0-mu**aexp)**(bexp-1.d0)

      vkk=0.d0
      ekk=0.d0
      do 212 i=1,ndim
        vkk=vkk+(gradv(1,i)+e0*dstrn(1,i)/w0)*Finv(i,1)
     &         +(gradv(2,i)+e0*dstrn(2,i)/w0)*Finv(i,2)
        ekk=ekk+strn(i,i)-strn0
 212  continue
!      vkk=vkk+dJac/Jac
      
      ekkplus=0.5d0*(ekk+dabs(ekk))
      ekkminus=0.5d0*(ekk-dabs(ekk))
      if (ekk.ge.0.d0) then
        dkkplus=1.d0
        dkkminus=0.d0
      else
        dkkplus=0.d0
        dkkminus=1.d0
      endif

      crit_strain=e0*1.d0
     &     *dsqrt((1.d0+poisson)*(1.d0-2.d0*poisson)/3.d0)
      istr=-crit_strain*0.5d0

      if (Jac-(1.d0+istr).lt.0.d0) then
        HeavisideJ1=0.d0
        dHJ1=0.d0
      elseif (Jac-(1.d0+istr).ge.crit_strain) then
        HeavisideJ1=1.d0
        dHJ1=0.d0
      else
        HeavisideJ1=
     &  (3.d0/((crit_strain)**2.d0))
     &         *(Jac-(1.d0+istr))**2.d0
     & -(2.d0/((crit_strain)**3.d0))
     &         *(Jac-(1.d0+istr))**3.d0
        dHJ1=-6.d0*(Jac-(1.d0+istr))*(Jac-(1.d0+istr)-crit_strain)
     &        /((crit_strain)**3.d0)
        dHJ1=0.d0
      endif

      if (Jac_old-(1.d0+istr).lt.0.d0) then
        HeavisideJ2=0.d0
      elseif (Jac_old-(1.d0+istr).ge.crit_strain) then
        HeavisideJ2=1.d0
      else
        HeavisideJ2=
     &  (3.d0/((crit_strain)**2.d0))
     &         *(Jac_old-(1.d0+istr))**2.d0
     & -(2.d0/((crit_strain)**3.d0))
     &         *(Jac_old-(1.d0+istr))**3.d0
      endif

      dHJ1=0.d0
      if (Jac.ge.1.d0) then
        HeavisideJ1=1.d0
      else
        HeavisideJ1=1.d0
        if (tc_split_on) HeavisideJ1=0.d0
        if (einject.gt.1.d-6) HeavisideJ1=1.d0
      endif

      if (Jac_old.ge.1.d0) then
        HeavisideJ2=1.d0
      else
        HeavisideJ2=1.d0
        if (tc_split_on) HeavisideJ2=0.d0
        if (einject.gt.1.d-6) HeavisideJ2=1.d0
      endif

      do 213 i=1,ndim
        do 223 j=1,ndim
          if (i.eq.j) then
            vbar(i,j)=-vkk/3.d0
          else
            vbar(i,j)=0.d0
          endif
          do 233 k=1,ndim
            vbar(i,j)=vbar(i,j)+(gradv(i,k)+e0*dstrn(i,k)/w0)*Finv(k,j)
 233      continue
 223    continue
 213  continue

      crit_stretch=1.d0+h0
      if (stretch.gt.2.d0*crit_stretch-1.d0) then
        vscale_factor=-3.d0/(stretch*(stretch-1.d0))
      elseif (stretch.gt.1.d0) then
        vscale_factor=6.d0*(stretch-1.d0)*(
     &   1/((stretch-2.d0)*stretch+(1.d0-2.d0*crit_stretch)**2.d0)
     &  -1/(5.d0+(stretch-2.d0)*stretch
     &      +4.d0*(crit_stretch-2.d0)*crit_stretch))
      else
        vscale_factor=0.d0
      endif

      s_numerator=0.d0
      s_denominator=0.d0
      do 401 i=1,ndim
        s_numerator=s_numerator+gmu(i)*gmu(i)
        do 402 j=1,ndim
          do 403 k=1,ndim
            s_denominator=s_denominator
     &                    +Finv(j,i)*Finv(k,i)*gmu(j)*gmu(k)
 403      continue             
 402    continue
 401  continue
      
      do 404 i=1,ndim
        do 405 j=1,ndim
          dstretchdF(i,j)=0.d0
          dJdu(i,j)=e0*Jac*Finv(j,i)
          do 406 k=1,ndim
            if (s_denominator.lt.1.d-15) then
              dstretchdF(i,j)=0.d0
            else
              dstretchdF(i,j)=2.d0*gmu(i)*Finv(k,j)*gmu(k)*
     &      dsqrt(s_numerator)/(s_denominator*dsqrt(s_denominator))
            endif
            do 407 l=1,ndim
              dFdu(i,j,k,l)=-e0*Finv(i,k)*Finv(l,j)
 407        continue
 406      continue             
 405    continue
 404  continue

      dstretchdJ=1.d0

      phicrack=1.d0
      phieff(1)=phi+phiF(1)*(1.d0-phi)
     &                *((Jac-1.d0)*HeavisideJ1)/Jac
      phieff(2)=phi_old+phiF(2)*(1.d0-phi_old)
     &                *((Jac_old-1.d0)*HeavisideJ2)/Jac_old

      fdensity(1)=1.d0+Cf*p
      fdensity(2)=1.d0+Cf*p_old

      dphibdmu=-DphiF(1)*(biot*(Jac-1.d0)+(p-p0)/Cphi)*HeavisideJ1
     &          /((Jac-phiF(1)*(Jac-1.d0)*HeavisideJ1)**2.d0)
      dphibdp=(1.d0-phiF(1)*HeavisideJ1)
     &        /(Cphi*(Jac-(Jac-1.d0)*phiF(1)*HeavisideJ1))
      dphibdJ=-phi0/(Jac**2.d0)
     &        +(biot-(phi-phi0/Jac))*(1.d0-phiF(1)*HeavisideJ1)
     &                    /(Jac-(Jac-1.d0)*phiF(1)*HeavisideJ1)

      if (dabs(1.d0-phiF(1)*HeavisideJ1).gt.1.d-15)
     &  dphibdJ=dphibdJ+phi*dHJ1/(1.d0-phiF(1)*HeavisideJ1)

      dphieffdmu=dphibdmu+(Jac-1.d0)*HeavisideJ1
     &               *((1.d0-phi)*DphiF(1)-phiF(1)*dphibdmu)/Jac
      dphieffdJ=dphibdJ+phiF(1)*
     &             ((HeavisideJ1+(Jac-1.d0)*dHJ1)*(1.d0-phi)/Jac
     &                      -(Jac-1.d0)*HeavisideJ1*dphibdJ)/Jac
      dphieffdp=(1.d0-phiF(1)*(Jac-1.d0)*HeavisideJ1/Jac)*dphibdp

      do 100 i=1,ndim
       dbsfdmu(i)=0.d0
       do 101 j=1,ndim
        dbsfdmu(i)=dbsfdmu(i)+DdarcyF*imperm(i,j)*v(j)*kperm
        dbsfdv(i,j)=0.d0
        dstrsdp(i,j)=-((1.d0-phiF(1)*HeavisideJ1)*biot
     &                      +phiF(1)*HeavisideJ1)*delta(i,j)
        dtaudmu(i,j)=-2.d0*(1.d0-mu)*visc_eff*(vbar(i,j)+vbar(j,i))
        do 102 k=1,ndim
         dbsfdv(i,j)=dbsfdv(i,j)
     &               -(1.d0-darcyF)*imperm(i,k)*delta(k,j)*kperm
         do 103 l=1,ndim
          dJFdu(i,j,k,l)=e0*Jac*Finv(l,k)*Finv(i,j)
     &                  -e0*Jac*Finv(i,k)*Finv(l,j)
          dtaudv(i,j,k,l)=(1.d0-mu)**2.d0*visc_eff
     &                            *(delta(i,k)*delta(j,l)+
     &                              delta(j,k)*delta(i,l)-
     &                    2.d0/3.d0*delta(i,j)*delta(k,l))
          dtaudu(i,j,k,l)=(1.d0-mu)**2.d0*visc_eff*e0/dt
     &                            *(delta(i,k)*Finv(l,j)+
     &                              delta(j,k)*Finv(l,i)-
     &                    2.d0/3.d0*delta(i,j)*Finv(l,k))/w0
          do 104 m=1,ndim
           dtaudu(i,j,k,l)=dtaudu(i,j,k,l)-(1.d0-mu)**2.d0*e0*visc_eff*
     &         ((gradv(i,m)+e0*dstrn(i,m)/w0)*Finv(m,k)*Finv(l,j)
     &         +(gradv(j,m)+e0*dstrn(j,m)/w0)*Finv(m,k)*Finv(l,i))
           do 105 n=1,ndim
            dtaudu(i,j,k,l)=dtaudu(i,j,k,l)+(1.d0-mu)**2.d0*e0*visc_eff*
     &        2.d0/3.d0*delta(j,i)*(gradv(m,n)+e0*dstrn(m,n)/w0)
     &        *Finv(n,k)*Finv(l,m)

!            dtaudu(i,j,k,l)=dtaudu(i,j,k,l)+(tau(i,j)+p*delta(i,j))*
!     &                      vscale_factor*dstretchdF(m,n)*dFdu(m,n,k,l)

            dtaudu(i,j,k,l)=dtaudu(i,j,k,l)+(tau(i,j)+p*delta(i,j))*
     &                      vscale_factor*dstretchdJ*dJdu(k,l)
 105       continue
 104      continue
 103     continue
 102    continue
 101   continue
 100  continue

      phiterm=-(1.d0-biot)*phiF(1)*p

      if (dabs(strn(1,2)+strn(2,1)).lt.1.d-14) then
        if (strn(1,1).ge.strn(2,2)) then
          theta=0.d0
        else
          theta=pi/2.d0
        endif
      elseif (dabs(strn(1,1)-strn(2,2)).lt.1.d-14) then
        if (strn(1,2)+strn(2,1).ge.0.d0) then
          theta=pi/4.d0
        else
          theta=3.d0*pi/4.d0
        endif
      else
        theta=0.5d0*datan2(strn(1,2)+strn(2,1),strn(1,1)-strn(2,2))
      endif

      nxn(1,1)=dcos(theta)
      nxn(1,2)=-dsin(theta)
      nxn(2,1)=dsin(theta)
      nxn(2,2)=nxn(1,1)

      do 120 i=1,ndim
        do 121 j=1,ndim
          eplus(i,j)=0.d0
          eminus(i,j)=0.d0
 121    continue
 120  continue

      etotal(1)=0.5d0*(strn(1,1)+strn(2,2)-2.d0*strn0
     &             +dsqrt((strn(1,1)-strn(2,2))**2.d0
     &                  +(strn(1,2)+strn(2,1))**2.d0))
      etotal(2)=0.5d0*(strn(1,1)+strn(2,2)-2.d0*strn0
     &             -dsqrt((strn(1,1)-strn(2,2))**2.d0
     &                  +(strn(1,2)+strn(2,1))**2.d0))

      eplus(1,1)=0.5d0*(etotal(1)+dabs(etotal(1)))
      eplus(2,2)=0.5d0*(etotal(2)+dabs(etotal(2)))
      eminus(1,1)=0.5d0*(etotal(1)-dabs(etotal(1)))
      eminus(2,2)=0.5d0*(etotal(2)-dabs(etotal(2)))

      do 313 i=1,ndim
        do 323 j=1,ndim
          strnplus(i,j)=0.d0
          strnminus(i,j)=0.d0
          do 333 k=1,ndim
            do 343 q=1,ndim
              strnplus(i,j)=strnplus(i,j)+
     &                      nxn(i,k)*eplus(k,q)*nxn(j,q)
              strnminus(i,j)=strnminus(i,j)+
     &                       nxn(i,k)*eminus(k,q)*nxn(j,q)
 343        continue
 333      continue
 323    continue
 313  continue

      no_tc_split=.not.tc_split_on
      if (etotal(2).ge.0.d0.or.no_tc_split) then
        dstrsdstrn(1,1,1,1)=c11*factor+phiterm*dHJ1*dJdu(1,1)
        dstrsdstrn(1,1,2,2)=c13*factor+phiterm*dHJ1*dJdu(2,2)
        dstrsdstrn(1,1,1,2)=0.d0
        dstrsdstrn(1,1,2,1)=0.d0
        dstrsdstrn(2,2,1,1)=dstrsdstrn(1,1,2,2)
        dstrsdstrn(2,2,2,2)=c33*factor+phiterm*dHJ1*dJdu(2,2)
        dstrsdstrn(2,2,1,2)=0.d0
        dstrsdstrn(2,2,2,1)=0.d0
        dstrsdstrn(1,2,1,1)=0.d0
        dstrsdstrn(1,2,2,2)=0.d0
        dstrsdstrn(1,2,1,2)=c44*factor
        dstrsdstrn(1,2,2,1)=c44*factor
        dstrsdstrn(2,1,1,1)=0.d0
        dstrsdstrn(2,1,2,2)=0.d0
        dstrsdstrn(2,1,1,2)=dstrsdstrn(1,2,1,2)
        dstrsdstrn(2,1,2,1)=dstrsdstrn(1,2,2,1)
      elseif (etotal(1).lt.-1.d-13) then
        dstrsdstrn(1,1,1,1)=c11+phiterm*dHJ1*dJdu(1,1)
        dstrsdstrn(1,1,2,2)=c13+phiterm*dHJ1*dJdu(2,2)
        dstrsdstrn(1,1,1,2)=0.d0
        dstrsdstrn(1,1,2,1)=0.d0
        dstrsdstrn(2,2,1,1)=dstrsdstrn(1,1,2,2)
        dstrsdstrn(2,2,2,2)=c33+phiterm*dHJ1*dJdu(2,2)
        dstrsdstrn(2,2,1,2)=0.d0
        dstrsdstrn(2,2,2,1)=0.d0
        dstrsdstrn(1,2,1,1)=0.d0
        dstrsdstrn(1,2,2,2)=0.d0
        dstrsdstrn(1,2,1,2)=c44
        dstrsdstrn(1,2,2,1)=c44
        dstrsdstrn(2,1,1,1)=0.d0
        dstrsdstrn(2,1,2,2)=0.d0
        dstrsdstrn(2,1,1,2)=dstrsdstrn(1,2,1,2)
        dstrsdstrn(2,1,2,1)=dstrsdstrn(1,2,2,1)
      else
        z=dsqrt((strn(1,1)-strn(2,2))**2.d0
     &         +(strn(1,2)+strn(2,1))**2.d0)

        dstrsdstrn(1,1,1,1)=factor*
     &    (c13*dkkplus+2.d0*c44*
     &     ((-(strn(1,1)-strn0)**3.d0+3.d0*(strn(1,1)-strn0)*(z**2.d0)
     &       +2.d0*(z**3.d0)
     &       +((strn(1,1)-strn0)**2.d0)*(strn(2,2)-strn0)
     &       -(strn(2,2)-strn0)*(z**2.d0)
     &       +(strn(1,1)-strn0)*((strn(2,2)-strn0)**2.d0)
     &       -(strn(2,2)-strn0)**3.d0)/(4.d0*z**3.d0)))+
     &    (c13*dkkminus+2.d0*c44*
     &     (((strn(1,1)-strn0)**3.d0-3.d0*(strn(1,1)-strn0)*(z**2.d0)
     &       +2.d0*(z**3.d0)
     &       -((strn(1,1)-strn0)**2.d0)*(strn(2,2)-strn0)
     &       +(strn(2,2)-strn0)*(z**2.d0)
     &       -(strn(1,1)-strn0)*((strn(2,2)-strn0)**2.d0)
     &       +(strn(2,2)-strn0)**3.d0)/(4.d0*z**3.d0)))
     &    +phiterm*dHJ1*dJdu(1,1)

        dstrsdstrn(1,2,1,1)=factor*
     &    (2.d0*c44*
     &     ((strn(1,2)+strn(2,1))*((strn(1,2)+strn(2,1))**2.d0
     &      +2.d0*(strn(2,2)-strn0)*(strn(2,2)-strn(1,1))))
     &     /(4.d0*(z**3.d0)))+
     &    (2.d0*c44*
     &     ((strn(1,2)+strn(2,1))*(-(strn(1,2)+strn(2,1))**2.d0
     &      +2.d0*(strn(2,2)-strn0)*(strn(1,1)-strn(2,2))))
     &     /(4.d0*(z**3.d0)))

        dstrsdstrn(2,1,1,1)=dstrsdstrn(1,2,1,1)

        dstrsdstrn(2,2,1,1)=factor*
     &    (c13*dkkplus+2.d0*c44*
     &     ((-(strn(1,2)+strn(2,1))**2.d0
     &      *(strn(1,1)+strn(2,2)-2.d0*strn0))/(4.d0*z**3.d0)))+
     &    (c13*dkkminus+2.d0*c44*
     &     (((strn(1,2)+strn(2,1))**2.d0
     &      *(strn(1,1)+strn(2,2)-2.d0*strn0))/(4.d0*z**3.d0)))
     &    +phiterm*dHJ1*dJdu(1,1)

        dstrsdstrn(1,1,1,2)=dstrsdstrn(1,2,1,1)

        dstrsdstrn(1,2,1,2)=factor*
     &    (2.d0*c44*
     &     (((strn(1,1)-strn0)**3.d0+z**3.d0
     &       -(strn(1,1)-strn0)**2.d0*(strn(2,2)-strn0)
     &       -(strn(1,1)-strn0)*(strn(2,2)-strn0)**2.d0
     &       +(strn(2,2)-strn0)**3.d0)/(4.d0*z**3.d0)))+
     &    (2.d0*c44*
     &     ((-(strn(1,1)-strn0)**3.d0+z**3.d0
     &       +(strn(1,1)-strn0)**2.d0*(strn(2,2)-strn0)
     &       +(strn(1,1)-strn0)*(strn(2,2)-strn0)**2.d0
     &       -(strn(2,2)-strn0)**3.d0)/(4.d0*z**3.d0)))

        dstrsdstrn(2,1,1,2)=dstrsdstrn(1,2,1,2)

        dstrsdstrn(2,2,1,2)=factor*
     &    (2.d0*c44*
     &     (((strn(1,2)+strn(2,1))*((strn(1,2)+strn(2,1))**2.d0
     &       +2.d0*(strn(1,1)-strn0)*(strn(1,1)-strn(2,2))))/
     &      (4.d0*z**3.d0)))+
     &    (2.d0*c44*
     &     ((-(strn(1,2)+strn(2,1))*((strn(1,2)+strn(2,1))**2.d0
     &       +2.d0*(strn(1,1)-strn0)*(strn(1,1)-strn(2,2))))/
     &      (4.d0*z**3.d0)))

        dstrsdstrn(1,1,2,1)=dstrsdstrn(1,2,1,1)
        dstrsdstrn(1,2,2,1)=dstrsdstrn(1,2,1,2)
        dstrsdstrn(2,1,2,1)=dstrsdstrn(1,2,1,2)
        dstrsdstrn(2,2,2,1)=dstrsdstrn(2,2,1,2)
        dstrsdstrn(1,1,2,2)=dstrsdstrn(2,2,1,1)
        dstrsdstrn(1,2,2,2)=dstrsdstrn(2,2,1,2)
        dstrsdstrn(2,1,2,2)=dstrsdstrn(2,2,1,2)

        dstrsdstrn(2,2,2,2)=factor*
     &    (c13*dkkplus+2.d0*c44*
     &     ((-(strn(1,1)-strn0)**3.d0-(strn(1,1)-strn0)*(z**2.d0)
     &       +2.d0*(z**3.d0)
     &       +((strn(1,1)-strn0)**2.d0)*(strn(2,2)-strn0)
     &       +3.d0*(strn(2,2)-strn0)*(z**2.d0)
     &       +(strn(1,1)-strn0)*((strn(2,2)-strn0)**2.d0)
     &       -(strn(2,2)-strn0)**3.d0)/(4.d0*z**3.d0)))+
     &    (c13*dkkminus+2.d0*c44*
     &     (((strn(1,1)-strn0)**3.d0+(strn(1,1)-strn0)*(z**2.d0)
     &       +2.d0*(z**3.d0)
     &       -((strn(1,1)-strn0)**2.d0)*(strn(2,2)-strn0)
     &       -3.d0*(strn(2,2)-strn0)*(z**2.d0)
     &       -(strn(1,1)-strn0)*((strn(2,2)-strn0)**2.d0)
     &       +(strn(2,2)-strn0)**3.d0)/(4.d0*z**3.d0)))
     &    +phiterm*dHJ1*dJdu(2,2)

      endif


      if (.not.tc_split_on) then

      dstrsdmu(1,1)=(c13*ekk+2.d0*c44*strn(1,1))*dfactor
     &              -p*(-DphiF(1)*biot+DphiF(1))*HeavisideJ1
      dstrsdmu(2,2)=(c13*ekk+2.d0*c44*strn(2,2))*dfactor
     &              -p*(-DphiF(1)*biot+DphiF(1))*HeavisideJ1
      dstrsdmu(1,2)=(c44*strn(1,2)+c44*strn(2,1))*dfactor
      dstrsdmu(2,1)=dstrsdmu(1,2)

      else

      dstrsdmu(1,1)=(c13*ekkplus+2.d0*c44*strnplus(1,1))*dfactor
     &              -p*(-DphiF(1)*biot+DphiF(1))*HeavisideJ1
      dstrsdmu(2,2)=(c13*ekkplus+2.d0*c44*strnplus(2,2))*dfactor
     &              -p*(-DphiF(1)*biot+DphiF(1))*HeavisideJ1
      dstrsdmu(1,2)=(c44*strnplus(1,2)+c44*strnplus(2,1))*dfactor
      dstrsdmu(2,1)=dstrsdmu(1,2)

      endif

      beta2=0.d0
      if (artstiff) beta2=0.d0
          
      oNNo=0.5d0+ddfactor*Hpsi+(beta2+beta)/dt
!     &       -(ekk*HeavisideJ1*p*(1.d0-biot)
!     &       -0.5d0*HeavisideJ1*p*p/(Cphi*e0))*DDphiF(1)
      pNNo=dHdp
!              -(ekk*HeavisideJ1*(1.d0-biot)
!     &        -HeavisideJ1*p/(Cphi*e0))*DphiF(1)
      oNNp=-dphieffdmu*Jac/(dt*w0)
      pNNp=-dphieffdp*Jac/(dt*w0)
     &     -phieff(2)*Jac_old*fdensity(2)*Cf
     &      /((fdensity(1)**2.d0)*w0*dt)
      do 110 i=1,ndim
       oNNv(i)=-Jac*dbsfdmu(i)*timerule
       pBNp(i)=0.d0
       vNNp(i)=0.d0
       do 111 j=1,ndim
        pNBv(i,j)=-Jac*Finv(i,j)*1.d0*timerule
        pBNv(i,j)=Jac*Finv(i,j)*0.d0*timerule
        vNNv(i,j)=-Jac*dbsfdv(j,i)*timerule
        pNBu(i,j)=dstrsdp(j,i)
        vBNp(i,j)=-Jac*Finv(j,i)*timerule
        pBNp(i)=pBNp(i)-Cf*v(j)*Finv(i,j)*Jac/fdensity(1)
        vNNp(i)=vNNp(i)-Cf*gradp(j)*Finv(j,i)*Jac/fdensity(1)
        pNNp=pNNp+Cf*Cf*gradp(J)*v(i)*Finv(J,i)*Jac
     &            /(fdensity(1)**2.d0)
        uBNp(i,j)=-dJdu(j,i)*(dphieffdJ*Jac+phieff(1))/(dt*w0)
        oNBu(i,j)=dstrsdmu(j,i)
        uBNo(i,j)=dHds(j,i)
!                           -p*DphiF(1)*(delta(i,j)*HeavisideJ1
!     &                      +ekk*dHJ1*dJdu(i,j))*(1.d0-biot)
        oBBo(i,j)=2.d0*delta(j,i)
        oNBv(i,j)=0.d0
        do 112 k=1,ndim
         oNBv(i,j)=oNBv(i,j)+dtaudmu(k,j)*Jac*Finv(i,k)*timerule
         uBNv(i,j,k)=-e0*Jac*Finv(j,i)*(bsf(k)+bf(k))*timerule
         do 113 l=1,ndim
          uBNp(i,j)=uBNp(i,j)-dJFdu(k,l,i,j)*gradv(l,k)*timerule
          uBBu(i,j,k,l)=dstrsdstrn(l,k,i,j)
          uBBv(i,j,k,l)=0.d0
          vBBv(i,j,k,l)=0.d0
          do 114 m=1,ndim
           uBBv(i,j,k,l)=uBBv(i,j,k,l)
     &               +dJFdu(k,m,i,j)*tau(m,l)*timerule
     &               +Jac*Finv(k,m)*dtaudu(m,l,i,j)*timerule
           do 115 n=1,ndim
            vBBv(i,j,k,l)=vBBv(i,j,k,l)+Jac*Finv(k,m)*Finv(j,n)
     &                *dtaudv(m,l,i,n)*timerule
 115       continue
 114      continue
 113     continue
 112    continue
 111   continue
 110  continue


      return
      end

**********************************************************************
**********************************************************************

      subroutine nodalstress(ndim,maxnd,maxele,nints,nnds,nodes,neles,
     &                    ele,gradp,gradmu,strn,strs,ngradp,ngradmu,
     &                    nstrn,nstrs,ele_per_node,ele_at_node,
     &                    eleproc,maxprocs,l_ele,det)
      implicit none
#include <petsc/finclude/petsc.h>

      integer maxnd,maxele,ndim,nodes,neles,nints,nnds,maxprocs
      integer ele(maxele*maxprocs,nnds),i,j,k,l,ele_per_node(maxnd),m
      integer ele_at_node(maxnd,10),l_ele(maxele*maxprocs)
      real*8 gradp(maxele,nints,ndim),gradmu(maxele,nints,ndim)
      real*8 strn(maxele,nints,ndim,ndim,2),strs(maxele,nints,ndim,ndim)
      real*8 ngradp(maxnd,ndim),ngradmu(maxnd,ndim),det(maxele,nints)
      real*8 nstrn(maxnd,ndim,ndim),nstrs(maxnd,ndim,ndim)
      real*8 ngradp_local(maxnd,ndim),ngradmu_local(maxnd,ndim)
      real*8 nstrn_local(maxnd,ndim,ndim),nstrs_local(maxnd,ndim,ndim)
      real*8 gradmu_avg(ndim),area
      logical eleproc(maxele*maxprocs)

      PetscErrorCode ierr
      
      do 1 i=1,nodes
      
        do 3 k=1,ndim
!          ngradp(i,k)=0.d0
          ngradmu_local(i,k)=0.d0
          do 4 l=1,ndim
!            nstrn(i,k,l)=0.d0
!            nstrs(i,k,l)=0.d0
 4        continue
 3      continue
        
        do 2 j=1,ele_per_node(i)
         if (eleproc(ele_at_node(i,j))) then
          do 211 k=1,ndim
!           ngradp(i,k)=ngradp(i,k)+gradp(j,1,k)
           gradmu_avg(k)=0.d0
           do 221 l=1,ndim
!           nstrn(i,k,l)=nstrn(i,k,l)+strn(j,1,k,l,1)
!           nstrs(i,k,l)=nstrs(i,k,l)+strs(j,1,k,l)
 221       continue
 211      continue
          area=0.d0
          do 103 m=1,nints
           area=area+det(l_ele(ele_at_node(i,j)),m)
           do 111 k=1,ndim
!            ngradp(i,k)=ngradp(i,k)+gradp(j,1,k)
            gradmu_avg(k)=gradmu_avg(k)
     &                   +gradmu(l_ele(ele_at_node(i,j)),m,k)
     &                    *det(l_ele(ele_at_node(i,j)),m)
            do 121 l=1,ndim
!            nstrn(i,k,l)=nstrn(i,k,l)+strn(j,1,k,l,1)
!            nstrs(i,k,l)=nstrs(i,k,l)+strs(j,1,k,l)
 121        continue
 111       continue
 103      continue
          do 311 k=1,ndim
!            ngradp(i,k)=ngradp(i,k)+gradp(j,1,k)
           gradmu_avg(k)=gradmu_avg(k)/area
           ngradmu_local(i,k)=ngradmu_local(i,k)+gradmu_avg(k)
           do 321 l=1,ndim
!           nstrn(i,k,l)=nstrn(i,k,l)+strn(j,1,k,l,1)
!           nstrs(i,k,l)=nstrs(i,k,l)+strs(j,1,k,l)
 321       continue
 311      continue
         endif
 2      continue
 
        if (ele_per_node(i).gt.0) then
        do 19 k=1,ndim
!          ngradp(i,k)=ngradp(i,k)/ele_per_node(i)
!          if (abs(ngradp(i,k)).lt.1.d-16) ngradp(i,k)=0.d0
          call MPI_Allreduce(ngradmu_local(i,k),ngradmu(i,k),1,
     &                       MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,
     &                       ierr)
          ngradmu(i,k)=ngradmu(i,k)/ele_per_node(i)
          if (abs(ngradmu(i,k)).lt.1.d-16) ngradmu(i,k)=0.d0

!         FIX
!          if (k.eq.2) then
!            ngradmu(i,k)=1.d0
!          else
!            ngradmu(i,k)=0.d0
!          endif
!          do 29 l=1,ndim
!            nstrn(i,k,l)=nstrn(i,k,l)/ele_per_node(i)
!            if (abs(nstrn(i,k,l)).lt.1.d-16) nstrn(i,k,l)=0.d0
!            nstrs(i,k,l)=nstrs(i,k,l)/ele_per_node(i)
!            if (abs(nstrs(i,k,l)).lt.1.d-16) nstrs(i,k,l)=0.d0
! 29       continue
 19     continue
        endif
 
 1    continue

      if (nnds.eq.8.or.nnds.eq.9) then
       do 102 j=1,neles
        do 119 k=1,ndim
          ngradmu(ele(j,5),k)=
     &         0.5d0*(ngradmu(ele(j,1),k)+ngradmu(ele(j,2),k))
          ngradmu(ele(j,6),k)=
     &         0.5d0*(ngradmu(ele(j,1),k)+ngradmu(ele(j,3),k))
          ngradmu(ele(j,7),k)=
     &         0.5d0*(ngradmu(ele(j,2),k)+ngradmu(ele(j,4),k))
          ngradmu(ele(j,8),k)=
     &         0.5d0*(ngradmu(ele(j,3),k)+ngradmu(ele(j,4),k))
          if (nnds.eq.9) ngradmu(ele(j,9),k)=
     &        0.25d0*(ngradmu(ele(j,1),k)+ngradmu(ele(j,2),k))
     &       +0.25d0*(ngradmu(ele(j,3),k)+ngradmu(ele(j,4),k))
 119    continue
 102   continue
      endif
          
      return
      end
      
*********************************************************************

      subroutine paraview(num,ou,maxnd,maxele,ndim,nodes,neles,nnds,
     &                  node,ele,u,mu,v,p,ngradp,nstrn,nstrs,e0,
     &                  ngradmu,nints,NN,phi,w0,maxprocs)
      implicit none
      integer maxnd,maxele,ndim,nnds
      real*8 node(maxnd,ndim),u(maxnd,ndim),mu(maxnd)
      real*8 v(maxnd,ndim),p(maxnd),ngradp(maxnd,ndim),e0
      real*8 nstrn(maxnd,ndim,ndim),nstrs(maxnd,ndim,ndim)
      real*8 ngradmu(maxnd,ndim),xx,yy,NN(nnds,nints)
      real*8 phi(maxele,nints,2),w0
      integer ele(maxele*maxprocs,nnds),maxprocs
      integer nodes,neles,num,i,nints,j,k,is9
      character(5) fnum
      character(1) ou


****************   ParaView Post-Processing   ********************

          write(fnum,915)num
          
          open(unit=11,file='results'// ou //'/sol/sol'
     &                      //trim(adjustl(fnum)))
     
          write(11,*)'  X_node','  Y_node','  u_x','  u_y',
     &      '  mu','  p','  v_x','  v_y'

          do 850 i=1,nodes
            write(11,*) ,node(i,1),node(i,2),u(i,1)*e0,u(i,2)*e0,
     &                   mu(i),p(i),v(i,1)*w0,v(i,2)*w0
 850      continue

          close(unit=11)

          print*,'|  ...WRITING: '// ou //'/sol'
     &                      //trim(adjustl(fnum))//'.vtk'
          print*,'|'

          is9=0
          if (nnds.eq.9) is9=1
     
          open(unit=20,file='paraview'// ou //'/sol'
     &                      //trim(adjustl(fnum))//'.vtk')
          write(20,910) '# vtk DataFile Version 3.0'
          write(20,910) 'Solution'
          write(20,910) 'ASCII'
          write(20,910) 'DATASET UNSTRUCTURED_GRID'
          write(20,911) 'POINTS ',nodes,' float'
          do 201 i=1,nodes
           write(20,*)node(i,1),node(i,2),0.0d0
 201      continue
          write(20,913) 'CELLS ',neles,neles*(nnds+1-is9)
          do 202 i=1,neles
           if (nnds.eq.8.or.nnds.eq.9) then
            write(20,*)8,ele(i,1)-1,ele(i,2)-1,ele(i,4)-1,
     &       ele(i,3)-1,ele(i,5)-1,ele(i,7)-1,ele(i,8)-1,ele(i,6)-1
           elseif (nnds.eq.4) then
            write(20,*)nnds,ele(i,1)-1,ele(i,2)-1,ele(i,4)-1,
     &       ele(i,3)-1
           elseif (nnds.eq.6) then
            write(20,*)6,ele(i,1)-1,ele(i,2)-1,ele(i,3)-1,
     &       ele(i,4)-1,ele(i,6)-1,ele(i,5)-1
           endif
 202      continue
          write(20,914) 'CELL_TYPES ',neles
          do 203 i=1,neles
           if (nnds.eq.8.or.nnds.eq.9) then
            write(20,*)23
           elseif (nnds.eq.4) then
            write(20,*)9
           elseif (nnds.eq.6) then
            write(20,*)22
           endif
 203      continue
          write(20,914) 'POINT_DATA ',nodes
          write(20,910) 'SCALARS Phase-field float 1'
          write(20,910) 'LOOKUP_TABLE default'
          do 204 i=1,nodes
           write(20,916)mu(i)
 204      continue
          write(20,910) 'SCALARS Pressure float 1'
          write(20,910) 'LOOKUP_TABLE default'
          do 205 i=1,nodes
           write(20,916)p(i)
 205      continue
          write(20,910) 'VECTORS Displacement float'
          do 206 i=1,nodes
           write(20,912)u(i,1)*e0,u(i,2)*e0,0.d0
 206      continue
          write(20,910) 'VECTORS Fluid_Flow float'
          do 208 i=1,nodes
           write(20,912)v(i,1)*w0,v(i,2)*w0,0.0d0
 208      continue
          write(20,910) 'VECTORS grad_p float'
          do 207 i=1,nodes
           write(20,*)ngradp(i,1),ngradp(i,2),0.0d0
 207      continue
          write(20,910) 'VECTORS grad_mu float'
          do 211 i=1,nodes
           write(20,*)ngradmu(i,1),ngradmu(i,2),0.0d0
 211      continue
          write(20,910) 'TENSORS Strain float'
          do 209 i=1,nodes
           write(20,*)nstrn(i,1,1),nstrn(i,1,2),0.0d0
           write(20,*)nstrn(i,2,1),nstrn(i,2,2),0.0d0
           write(20,*)0.0d0,0.0d0,0.0d0
           write(20,*)
 209      continue
          write(20,910) 'TENSORS Stress float'
          do 210 i=1,nodes
           write(20,*)nstrs(i,1,1),nstrs(i,1,2),0.0d0
           write(20,*)nstrs(i,2,1),nstrs(i,2,2),0.0d0
           write(20,*)0.0d0,0.0d0,0.0d0
           write(20,*)
 210      continue

          close(unit=20)

          if (1.eq.0) then

          open(unit=21,file='paraview'// ou //'/norms'
     &                      //trim(adjustl(fnum))//'.vtk')
          write(21,910) '# vtk DataFile Version 3.0'
          write(21,910) 'Solution'
          write(21,910) 'ASCII'
          write(21,910) 'DATASET UNSTRUCTURED_GRID'
          write(21,911) 'POINTS ',neles*nints,' float'
          do 771 i=1,neles
           do 772 j=1,nints
            xx=0.d0
            yy=0.d0
            do 773 k=1,nnds
             xx=xx+NN(k,j)*node(ele(i,k),1)
             yy=yy+NN(k,j)*node(ele(i,k),2)
 773        continue
            write(21,*)xx,yy,0.0d0
 772       continue
 771      continue
          write(21,913) 'CELLS ',neles,neles*(nints+1)
          do 302 i=1,neles
            write(21,*)nints,((i-1)*nints+j, j=0,nints-1)
 302      continue
          write(21,914) 'CELL_TYPES ',neles
          do 303 i=1,neles
            write(21,*)2
 303      continue
          write(21,914) 'POINT_DATA ',neles*nints
          write(21,910) 'SCALARS Phi float 1'
          write(21,910) 'LOOKUP_TABLE default'
          do 304 i=1,neles
           do 305 j=1,nints
            write(21,916)phi(i,j,1)
 305       continue
 304      continue

          close(unit=21)

          endif

 
 910      format(A)
 911      format(A,I8,A)
 912      format(F30.15, F30.15, F30.15)
 913      format(A,I8, ' ',I8)
 914      format(A,I8)
 915      format(I5)
 916      format(F30.15)
  
!******************************************************************


      return
      end

***********************************************************************

      subroutine sols(ndpn,ndofs,maxwk,maxnd,maxa,a,drhs)
      real*8 a(maxwk),drhs(maxnd*(ndpn-1)),t,r,s
      integer maxa(maxnd*(ndpn-1)+1),ndofs,maxwk,ndpn
      integer je,m1,jb,j,jp,jp1,i1,ie,ij,i,ib,kp,k


      je=maxa(1)
      m1=2
      t=a(je)
      a(je)=1.d0/t

      do 10 j=m1,ndofs
       jb=je+1
       je=maxa(j)
       jp=je-j
       jp1=jp-1
       i1=jb-jp
       ie=maxa(i1)
       ij=jb

       do 20 i=i1+1,j-1
         ib=ie+1
         ie=maxa(i)
         kp=ie-i-jp
         ij=ij+1
         r=a(ij)

         do 30 k=max(jb,ib-kp),jp1+i
           r=r-a(k)*a(k+kp)
 30      continue

         a(ij)=r

 20    continue

       t=a(je)
       ij=jb

       do 40 i=i1,j-1
         r=a(ij)
         s=r*a(maxa(i))
         a(ij)=s
         ij=ij+1
         t=t-r*s
 40    continue

       a(je)=1.d0/t

 10   continue

      kp=maxa(1)

      do 50 j=2,ndofs
       r=drhs(j)
       kp=kp+1

       do 60 i=j+kp-maxa(j),j-1
         r=r-a(kp)*drhs(i)
         kp=kp+1
 60    continue

       drhs(j)=r

 50   continue

      do 70 i=1,ndofs
       drhs(i)=a(maxa(i))*drhs(i)
 70   continue

      kp=maxa(ndofs)

      do 80 j=ndofs,2,-1
       r=drhs(j)
       kp=kp-1

       do 90 i=j-1,j-kp+maxa(j-1),-1
         drhs(i)=drhs(i)-a(kp)*r
         kp=kp-1
 90    continue

 80   continue


      return
      end

*********************************************************************

      subroutine zerorates(maxele,nints,ndim,neles,dJacob,dstrn,dphi,
     &                     duint,dmuint)
      integer maxele,ndim,nints
      real*8 dJacob(maxele,nints),dphi(maxele,nints)
      real*8 dstrn(maxele,nints,ndim,ndim)
      real*8 duint(maxele,nints,ndim),dmuint(maxele,nints)
      integer i,j,k,l,neles

      
        do 68 i=1,neles
          do 69 j=1,nints
            dJacob(i,j)=0.d0
            dphi(i,j)=0.d0
            dmuint(i,j)=0.d0
            do 70 k=1,ndim
              duint(i,j,k)=0.d0
              do 71 l=1,ndim
                dstrn(i,j,k,l)=0.d0
 71          continue
 70        continue
 69      continue
 68    continue

 
      return
      end
      
***********************************************************************

      subroutine stfLsq(nnds,nints,ndim,ndpn,B,NN,det,node,maxele,
     &            maxnd,maxcon,maxdbc,nvpdofs,neles,ele,gcon,mpc,
     &            mpcnum,conA,nvpdbcs,nmudofs,nmudbcs,nudofs,
     &            nudbcs,PetG,rank,eleproc)
      implicit none
#include <petsc/finclude/petsc.h>

      real*8 B(maxele,nnds,nints,ndim),det(maxele,nints)
      real*8 bele(nnds,nints,ndim),detele(nints)
      real*8 kele(nnds,nnds)
      real*8 node(maxnd,ndim),conA(maxdbc,maxcon)
      real*8 NN(nnds,nints),kcon
      real*8 Petkele(nnds*(ndpn)*(maxcon+1)*ndim,
     &               nnds*(ndpn)*(maxcon+1)*ndim)
      integer maxele,maxnd,maxcon,maxdbc,neles,ndpn
      integer ele(maxele,nnds),gcon(maxnd,ndpn)
      integer i,j,k,l,n,nvpdofs,nvpdbcs,nudofs,nudbcs
      integer dofi,dofj,kglb,nloci,nlocj,nmudofs,nmudbcs
      integer mpc(maxdbc,maxcon,2),mpcnum(maxdbc)
      integer ii,jj,coni,conj,nci,ncj,nnds,ndim,nints
      integer idx(nnds*(ndpn)*(maxcon+1)*ndim),Peti,Petj
      integer ndofs,ndbcs,rank
      logical eleproc(maxele)

      Mat PetG
      PetscErrorCode ierr


      call MatZeroEntries(PetG,ierr)

      ndofs=nvpdofs+nudofs+nmudofs
      ndbcs=nvpdbcs+nudbcs+nmudbcs

      do 100 l=1,neles
      if (eleproc(l)) then

        do 110 i=1,nnds
          do 111 k=1,nints
            do 112 j=1,ndim
              bele(i,k,j)=B(l,i,k,j)
 112        continue
 111      continue
 110    continue

        do 120 i=1,nints
          detele(i)=det(l,i)
 120    continue

        do 121 i=1,nnds*(maxcon+1)*(ndpn)*ndim
          idx(i)=-1
          do 122 j=1,nnds*(maxcon+1)*(ndpn)*ndim
            Petkele(i,j)=0.d0
 122      continue
 121    continue

        call stfeLsq(nints,nnds,ndim,ndpn,maxele,l,kele,detele,bele,NN)

        do 130 i=1,nnds*(ndpn)
          do 140 j=1,nnds*(ndpn)

            nloci=(i-1)/(ndpn)+1
            dofi=mod(i-1,ndpn)+1
            nlocj=(j-1)/(ndpn)+1
            dofj=mod(j-1,ndpn)+1

            if (dofi.ne.dofj) goto 140            

            dofi=gcon(ele(l,nloci),dofi)
            dofj=gcon(ele(l,nlocj),dofj)
            coni=dofi-ndofs-ndbcs
            conj=dofj-ndofs-ndbcs
            
            Peti=(i-1)*(maxcon+1)+1
            Petj=(j-1)*(maxcon+1)+1

            if(coni.gt.0) then
              nci=mpcnum(coni)
              do 150 ii=1,nci
                dofi=gcon(mpc(coni,ii,1),mpc(coni,ii,2))
                if (dofi.gt.ndofs+ndbcs) goto 150
                if(conj.gt.0) then
                  ncj=mpcnum(conj)              
                  do 155 jj=1,ncj
                    dofj=gcon(mpc(conj,jj,1),
     &                          mpc(conj,jj,2))
                    if (dofj.gt.ndofs+ndbcs) goto 155

                    if (Peti+ii.eq.Petj+jj) then
                      idx(ndim*(Peti+ii-1)+1)=ndim*(dofi-1)
                      idx(ndim*(Peti+ii-1)+2)=ndim*(dofi-1)+1
                    endif

                    Petkele(ndim*(Petj+jj-1)+1,ndim*(Peti+ii-1)+1)=
     &                         kele(nloci,nlocj)
     &                         *conA(coni,ii)*conA(conj,jj)
                    Petkele(ndim*(Petj+jj-1)+2,ndim*(Peti+ii-1)+2)=
     &                         kele(nloci,nlocj)
     &                         *conA(coni,ii)*conA(conj,jj)


 155              continue
                else
                  if (dofj.gt.ndofs+ndbcs) goto 150

                  Petkele(ndim*(Petj-1)+1,ndim*(Peti+ii-1)+1)
     &             =kele(nloci,nlocj)*conA(coni,ii)
                  Petkele(ndim*(Petj-1)+2,ndim*(Peti+ii-1)+2)
     &             =kele(nloci,nlocj)*conA(coni,ii)

                endif
 150          continue
            elseif(conj.gt.0) then
              if (dofi.gt.ndofs+ndbcs) goto 140
              ncj=mpcnum(conj)
              do 160 jj=1,ncj
                dofj=gcon(mpc(conj,jj,1),mpc(conj,jj,2))
                if (dofj.gt.ndofs+ndbcs) goto 160

                Petkele(ndim*(Petj+jj-1)+1,ndim*(Peti-1)+1)
     &           =kele(nloci,nlocj)*conA(conj,jj)
                Petkele(ndim*(Petj+jj-1)+2,ndim*(Peti-1)+2)
     &           =kele(nloci,nlocj)*conA(conj,jj)

 160          continue
            else
              if (dofi.gt.ndofs+ndbcs) goto 140
              if (dofj.gt.ndofs+ndbcs) goto 140

              if (Peti.eq.Petj) then
                idx(ndim*(Peti-1)+1)=ndim*(dofi-1)
                idx(ndim*(Peti-1)+2)=ndim*(dofi-1)+1
              endif
                
              Petkele(ndim*(Petj-1)+1,ndim*(Peti-1)+1)
     &         =kele(nloci,nlocj)
              Petkele(ndim*(Petj-1)+2,ndim*(Peti-1)+2)
     &         =kele(nloci,nlocj)

            endif
 140      continue
 130    continue

        call MatSetValues(PetG,nnds*(ndpn)*(maxcon+1)*ndim,idx,
     &                         nnds*(ndpn)*(maxcon+1)*ndim,idx,
     &                    Petkele,ADD_VALUES,ierr)
      CHKERRQ(ierr)

      endif
 100  continue


      return
      end

**********************************************************************

      subroutine stfeLsq(nints,nnds,ndim,ndpn,maxele,iele,kele,
     &                detele,bele,nele)
      implicit none
      real*8 kele(nnds,nnds)
      real*8 detele(nints),bele(nnds,nints,ndim)
      real*8 nele(nnds,nints)
      integer maxele,nints,nnds,ndim,ndpn
      integer iele,i,ia,k,l,i1,j1,j,Peti,Petj


      do 100 i=1,nnds
        do 101 j=1,nnds
          kele(i,j)=0.d0
 101    continue
 100  continue
      
      do 200 ia=1,nints
     
        do 220 i1=1,nnds
          do 230 j1=1,nnds
            
            kele(i1,j1)=kele(i1,j1)
     &            +nele(j1,ia)*nele(i1,ia)*detele(ia)

 230      continue
 220    continue
 200  continue


      return
      end

**********************************************************************

      subroutine intfLsq(nints,ndim,ndpn,nnds,maxnd,nodes,maxele,neles,
     &                 ele,gcon,det,B,NN,ugrad,vgrad,pgrad,mugrad,
     &                 fint)
      implicit none
      real*8 det(maxele,nints),B(maxele,nnds,nints,ndim)
      real*8 NN(nnds,nints),fint(maxnd*ndpn*ndim),e0
      real*8 fele(nnds*ndpn*ndim),ugrad(maxele,nints,ndim,ndim)
      real*8 vgrad(maxele,nints,ndim,ndim)
      real*8 mugrad(maxele,nints,ndim)
      real*8 pgrad(maxele,nints,ndim)
      integer maxnd,maxele,nodes,neles,ele(maxele,nnds)
      integer gcon(maxnd,ndpn)
      integer i,j,k,l,m,ndof,nints,ndim,ndpn,nnds

      
      do 100 i=1,nodes*ndpn*ndim
 100    fint(i)=0.d0

      do 200 i=1,neles

        do 210 j=1,nnds*ndpn*ndim
 210      fele(j)=0.d0

        do 220 k=1,nints
          do 230 l=1,nnds
            m=ndim*ndpn*l
            fele(m-11)=fele(m-11)+det(i,k)*
     &           (NN(l,k)*ugrad(i,k,1,1))
            fele(m-10)=fele(m-10)+det(i,k)*
     &           (NN(l,k)*ugrad(i,k,1,2))
            fele(m-9)=fele(m-9)+det(i,k)*
     &           (NN(l,k)*ugrad(i,k,2,1))
            fele(m-8)=fele(m-8)+det(i,k)*
     &           (NN(l,k)*ugrad(i,k,2,2))
            fele(m-7)=fele(m-7)+det(i,k)*
     &           (NN(l,k)*vgrad(i,k,1,1))
            fele(m-6)=fele(m-6)+det(i,k)*
     &           (NN(l,k)*vgrad(i,k,1,2))
            fele(m-5)=fele(m-5)+det(i,k)*
     &           (NN(l,k)*vgrad(i,k,2,1))
            fele(m-4)=fele(m-4)+det(i,k)*
     &           (NN(l,k)*vgrad(i,k,2,2))
            fele(m-3)=fele(m-3)+det(i,k)*
     &           (NN(l,k)*pgrad(i,k,1))
            fele(m-2)=fele(m-2)+det(i,k)*
     &           (NN(l,k)*pgrad(i,k,2))
            fele(m-1)=fele(m-1)+det(i,k)*
     &           (NN(l,k)*mugrad(i,k,1))
            fele(m)=fele(m)+det(i,k)*
     &           (NN(l,k)*mugrad(i,k,2))
 230      continue
 220    continue

        do 300 j=1,nnds
          do 310 k=1,ndpn
            ndof=gcon(ele(i,j),k)
            fint(ndim*(ndof-1)+1)=fint(ndim*(ndof-1)+1)
     &             +fele(ndim*ndpn*(j-1)+ndim*(k-1)+1)
            fint(ndim*(ndof-1)+2)=fint(ndim*(ndof-1)+2)
     &             +fele(ndim*ndpn*(j-1)+ndim*(k-1)+2)
 310      continue
 300    continue

 200  continue


      return
      end

**********************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Set Nodal values of gradients. This is so long because gradient 
!!! shape functions mimic regular shape functions and therefor use
!!! mpc information.

      subroutine LsqGrads(nints,ndim,ndpn,nnds,maxnd,maxdbc,maxcon,
     &              nodes,maxele,nvpdofs,nudofs,nmudofs,nvpdbcs,
     &              nudbcs,nmudbcs,ndofs,ndbcs,indices,neles,
     &              ele,gcon,gconvp,gconu,gconmu,dbcval,dbcnd,
     &              dbcdof,mpcnum,mpc,mpcnd,conA,conB,vpmpcnum,
     &              nconvp,vpmpc,vpmpcnd,conAvp,conBvp,umpcnum,
     &              nconu,umpc,umpcnd,conAu,conBu,mumpcnum,nconmu,
     &              mumpc,mumpcnd,conAmu,conBmu,det,B,NN,
     &              ugrad,vgrad,pgrad,mugrad,nugrad,nvgrad,npgrad,
     &              nmugrad,kspG,LsqF,LsqU)
      implicit none
#include <petsc/finclude/petsc.h>

      real*8 det(maxele,nints),B(maxele,nnds,nints,ndim)
      real*8 NN(nnds,nints),LsqFint(maxnd*ndpn*ndim)
      real*8 drhs(maxnd*ndpn*ndim),drhsT(maxnd*ndpn*ndim)
      real*8 fele(nnds*ndpn*ndim),ugrad(maxele,nints,ndim,ndim)
      real*8 vgrad(maxele,nints,ndim,ndim),dbcval(maxdbc)
      real*8 mugrad(maxele,nints,ndim)
      real*8 pgrad(maxele,nints,ndim)
      real*8 nugrad(maxnd,ndim,ndim),nmugrad(maxnd,ndim)
      real*8 nvgrad(maxnd,ndim,ndim),npgrad(maxnd,ndim)
      real*8 udx(maxnd,ndim),udy(maxnd,ndim),pdx(maxnd),pdy(maxnd)
      real*8 vdx(maxnd,ndim),vdy(maxnd,ndim)
      real*8 mudx(2,maxnd),mudy(2,maxnd)
      real*8 conA(maxdbc,maxcon),conB(maxdbc,maxcon)
      real*8 conAu(maxdbc,maxcon),conBu(maxdbc,maxcon)
      real*8 conAvp(maxdbc,maxcon),conBvp(maxdbc,maxcon)
      real*8 conAmu(maxdbc,maxcon),conBmu(maxdbc,maxcon)
      integer maxnd,maxele,maxdbc,maxcon,nodes,neles
      integer i,j,k,l,m,jj,nints,ndim,ndpn,nnds,gcon(maxnd,ndpn)
      integer nudofs,umpc(maxdbc,maxcon,2),ele(maxele,nnds)
      integer umpcnd(maxdbc,2),nmudbcs,nconu,gconu(maxnd,ndpn)
      integer nvpdbcs,indices(maxnd*ndpn*ndim),nudbcs
      integer dbcnd(maxdbc),dbcdof(maxdbc)
      integer nmudofs,gconmu(maxnd,ndpn)
      integer umpcnum(maxdbc),nconmu,mumpcnum(maxdbc)
      integer mumpc(maxdbc,maxcon,2),mumpcnd(maxdbc,2) 
      integer nconvp,vpmpcnum(maxdbc),gconvp(maxnd,ndpn)
      integer vpmpc(maxdbc,maxcon,2),vpmpcnd(maxdbc,2),nvpdofs
      integer vpumpcnum(maxdbc),vpumpc(maxdbc,maxcon,2)
      integer vpumpcnd(maxdbc,2),ndofs,ndbcs
      integer mpcnum(maxdbc)
      integer mpcnd(maxdbc,2),mpc(maxdbc,maxcon,2)

      PetscErrorCode ierr
      Vec LsqF
      Vec LsqU
      KSP kspG


        call intfLsq(nints,ndim,ndpn,nnds,maxnd,nodes,maxele,neles,
     &                 ele,gcon,det,B,NN,ugrad,vgrad,pgrad,mugrad,
     &                 LsqFint)

        do 790 i=1,nodes*ndpn*ndim
          drhs(i)=LsqFint(i)
 790    continue
 
        do 796 i=1,nconvp+nconu+nconmu
          j=ndim*(ndofs+ndbcs+i-1)+1
          do 797 jj=1,mpcnum(i)
            k=ndim*(gcon(mpc(i,jj,1),mpc(i,jj,2))-1)+1
            drhs(k)=drhs(k)+drhs(j)*conA(i,jj)
            drhs(k+1)=drhs(k+1)+drhs(j+1)*conA(i,jj)
 797      continue
 796    continue

        call VecSetValues(LsqF,(ndofs+ndbcs)*ndim,
     &                     indices,drhs,INSERT_VALUES,ierr)
        call VecAssemblyBegin(LsqF,ierr)
        call VecAssemblyEnd(LsqF,ierr)
!        call VecView(LsqF,PETSC_VIEWER_STDOUT_WORLD,ierr)

        call KSPSolve(kspG,LsqF,LsqU,ierr)

!        call VecView(LsqU,PETSC_VIEWER_STDOUT_WORLD,ierr)
        call VecGetValues(LsqU,(ndofs+ndbcs)*ndim,
     &                     indices,drhsT,ierr)

        do 720 i=1,nodes
          udx(i,1)=0.d0
          udx(i,2)=0.d0
          vdx(i,1)=0.d0
          vdx(i,2)=0.d0
          mudx(1,i)=0.d0
          mudx(2,i)=0.d0  
          pdx(i)=0.d0
          udy(i,1)=0.d0
          udy(i,2)=0.d0
          vdy(i,1)=0.d0
          vdy(i,2)=0.d0
          mudy(1,i)=0.d0
          mudy(2,i)=0.d0  
          pdy(i)=0.d0
 720    continue


        do 760 i=1,nodes
          do 761 j=1,ndim+1
            if (gconvp(i,j).le.nvpdofs+nvpdbcs)
     &        drhs(gconvp(i,j))=drhsT(ndim*(gcon(i,j+ndim)-1)+1)
 761      continue
 760    continue

        call calcvp(ndpn,ndim,maxnd,maxdbc,maxcon,nodes,
     &              nvpdofs+nvpdbcs,gconvp,dbcval,drhs,vdx,pdx,0,
     &              dbcnd,dbcdof,vpmpcnum,0,nconvp,vpmpc,
     &              vpmpcnd,conAvp,conBvp,0)

        do 750 i=1,nodes
          do 751 j=1,ndim+1
            if (gconvp(i,j).le.nvpdofs+nvpdbcs)
     &        drhs(gconvp(i,j))=drhsT(ndim*(gcon(i,j+ndim)-1)+2)
 751      continue
 750    continue

        call calcvp(ndpn,ndim,maxnd,maxdbc,maxcon,nodes,
     &              nvpdofs+nvpdbcs,gconvp,dbcval,drhs,vdy,pdy,0,
     &              dbcnd,dbcdof,vpmpcnum,0,nconvp,vpmpc,
     &              vpmpcnd,conAvp,conBvp,0)

        do 752 i=1,nodes
          do 753 j=1,ndim
            if (gconu(i,j).le.nudofs+nudbcs)
     &        drhs(gconu(i,j))=drhsT(ndim*(gcon(i,j)-1)+1)
 753      continue
 752    continue

        call calcu(ndim,ndpn,maxnd,maxdbc,nodes,nudofs+nudbcs,nconu,
     &             gconu,dbcval,drhs,udx,umpc,0,dbcnd,dbcdof,
     &             conAu,conBu,umpcnd,0,0,umpcnum,maxcon)

        do 762 i=1,nodes
          do 763 j=1,ndim
            if (gconu(i,j).le.nudofs+nudbcs)
     &        drhs(gconu(i,j))=drhsT(ndim*(gcon(i,j)-1)+2)
 763      continue
 762    continue

        call calcu(ndim,ndpn,maxnd,maxdbc,nodes,nudofs+nudbcs,nconu,
     &             gconu,dbcval,drhs,udy,umpc,0,dbcnd,dbcdof,
     &             conAu,conBu,umpcnd,0,0,umpcnum,maxcon)

        do 754 i=1,nodes
          if (gconmu(i,1).le.nmudofs+nmudbcs)
     &      drhs(gconmu(i,1))=drhsT(ndim*(gcon(i,6)-1)+1)
 754    continue

        call calcmu(ndpn,ndim,maxnd,maxdbc,maxcon,nodes,
     &              nmudofs+nmudbcs,gconmu,dbcval,drhs,mudx,0,
     &              dbcnd,dbcdof,mumpcnum,0,nconmu,mumpc,
     &              mumpcnd,conAmu,conBmu,0)

        do 764 i=1,nodes
          if (gconmu(i,1).le.nmudofs+nmudbcs)
     &      drhs(gconmu(i,1))=drhsT(ndim*(gcon(i,6)-1)+2)
 764    continue

        call calcmu(ndpn,ndim,maxnd,maxdbc,maxcon,nodes,
     &              nmudofs+nmudbcs,gconmu,dbcval,drhs,mudy,0,
     &              dbcnd,dbcdof,mumpcnum,0,nconmu,mumpc,
     &              mumpcnd,conAmu,conBmu,0)

        do 721 i=1,nodes
          if (abs(udx(i,1)).lt.1.d-15) then
            nugrad(i,1,1)=0.d0
          else
            nugrad(i,1,1)=udx(i,1)
          endif
          if (abs(udx(i,2)).lt.1.d-15) then
            nugrad(i,2,1)=0.d0
          else
            nugrad(i,2,1)=udx(i,2)
          endif
          if (abs(vdx(i,1)).lt.1.d-15) then
            nvgrad(i,1,1)=0.d0
          else
            nvgrad(i,1,1)=vdx(i,1)
          endif
          if (abs(vdx(i,2)).lt.1.d-15) then
            nvgrad(i,2,1)=0.d0
          else
            nvgrad(i,2,1)=vdx(i,2)
          endif
          if (abs(mudx(2,i)).lt.1.d-15) then
            nmugrad(i,1)=0.d0
          else
            nmugrad(i,1)=mudx(2,i)
          endif
          if (abs(pdx(i)).lt.1.d-15) then
            npgrad(i,1)=0.d0
          else  
            npgrad(i,1)=pdx(i)
          endif
          if (abs(udy(i,1)).lt.1.d-15) then
            nugrad(i,1,2)=0.d0
          else
            nugrad(i,1,2)=udy(i,1)
          endif
          if (abs(udy(i,2)).lt.1.d-15) then
            nugrad(i,2,2)=0.d0
          else
            nugrad(i,2,2)=udy(i,2)
          endif
          if (abs(vdy(i,1)).lt.1.d-15) then
            nvgrad(i,1,2)=0.d0
          else
            nvgrad(i,1,2)=vdy(i,1)
          endif
          if (abs(vdy(i,2)).lt.1.d-15) then
            nvgrad(i,2,2)=0.d0
          else
            nvgrad(i,2,2)=vdy(i,2)
          endif
          if (abs(mudy(2,i)).lt.1.d-15) then
            nmugrad(i,2)=0.d0
          else
            nmugrad(i,2)=mudy(2,i)
          endif
          if (abs(pdy(i)).lt.1.d-15) then
            npgrad(i,2)=0.d0
          else  
            npgrad(i,2)=pdy(i)
          endif
 721    continue


      return
      end

*************************************************************************


      subroutine findcrack(ndim,maxnd,node,mu,xcrack,nodes)
      implicit none
      integer maxnd,i,ndim,nodes
      real*8 node(maxnd,ndim),mu(maxnd),xcrack,mumin
      real*8 mumax,xmin,xmax
      
         xmin=-100.d0
         xmax=1.d6
         mumin=0.d0
         mumax=1.d0
         do 851 i=1,nodes
          if (abs(node(i,1)-0.0d0).lt.1.d-3) then
           if(node(i,2).gt.xmin) then
            if(mu(i).ge.mumin.and.mu(i).le.0.1d0) then
             mumin=mu(i)
             xmin=node(i,2)
            endif
           endif
           if(node(i,2).lt.xmax) then
            if(mu(i).le.mumax.and.mu(i).ge.0.1d0) then
             mumax=mu(i)
             xmax=node(i,2)
            endif
           endif
          endif
 851     continue
         if(mumax.ne.mumin) then
          xcrack=((xmax-xmin)/(mumax-mumin))*(0.1d0-mumin)+xmin
         else
         xcrack=xmin
         endif
         
      return
      end
      
********************************************************************

      subroutine DtNstf(ndim,maxnd,maxbdnd,ndpn,node,gcon,KDtN,
     &                  DtNidx,DtNdofs,DtNnode,DtNdof)
      implicit none
      integer maxnd,ndim,ndpn,i,j,k,l,n,a,b,nFunc,maxbdnd
      integer DtNnodes,DtNdofs,kfold
      integer DtNnode(maxbdnd*ndim),DtNdof(maxbdnd*ndim)
      integer DtNidx(maxbdnd*ndim),gcon(maxnd,ndpn)
      real*8 YM,poisson,Pi,Ebar,vbar,kappa,dTheta,R
      real*8 alpha(2,2,4,4),DtNI,ThetaA,ThetaB,aConst,bConst
      real*8 KDtN(maxbdnd*ndim,maxbdnd*ndim),node(maxnd,ndim)
      real*8 asymfactor,bsymfactor


      kfold=2
      
      YM=1.d0
      poisson=0.25d0
      Pi=dacos(-1.d0)

      Ebar=YM/(1.d0-poisson**2.d0)
      vbar=poisson/(1.d0-poisson)

      open(unit=45,file='input/DtN_dofs')
      read(45,*)DtNnodes,DtNdofs
      do 2 i=1,DtNdofs
        read(45,*)DtNnode(i),DtNdof(i)
 2    continue
      close(unit=45)

      if (kfold.eq.1) then
        dTheta=2.d0*Pi/dble(DtNnodes)
      else
        dTheta=2.d0*Pi/dble((DtNnodes-1)*kfold)
      endif
      R=dsqrt(node(DtNnode(1),1)**2.d0+node(DtNnode(1),2)**2.d0)

      kappa=YM/(Pi*(1.d0+poisson/(1.d0-poisson))
     &            *(1.d0-poisson**2.d0)*R**2.d0)

      do 3 i=1,DtNdofs
        DtNidx(i)=gcon(DtNnode(i),DtNdof(i))-1
 3    continue

      do 11 nFunc=0,200

        n=nFunc

        do 21 i=1,2
         do 22 j=1,2
          do 23 k=1,4
           do 24 l=1,4
            alpha(i,j,k,l)=0.d0
 24        continue
 23       continue
 22      continue
 21     continue

        if (kfold.eq.1) then

        if (n.eq.0) then
          alpha(1,1,1,1)=kappa*0.5d0
          alpha(1,1,2,2)=kappa*0.5d0
          alpha(1,2,1,2)=kappa*0.5d0
          alpha(1,2,2,1)=-kappa*0.5d0
          do 31 k=1,4
           do 32 l=1,4
            alpha(2,1,k,l)=-alpha(1,2,k,l)
            alpha(2,2,k,l)=alpha(1,1,k,l)
 32        continue
 31       continue       
        elseif (n.eq.1) then
          alpha(1,1,1,1)=2.d0*kappa
          alpha(1,1,2,2)=4.d0*kappa
          alpha(1,1,3,3)=4.d0*kappa
          alpha(1,1,4,4)=2.d0*kappa
          alpha(1,1,2,3)=4.d0*kappa
          alpha(1,1,3,2)=4.d0*kappa

          alpha(1,2,1,2)=kappa*2.0d0
          alpha(1,2,1,3)=kappa*2.0d0
          alpha(1,2,2,1)=-kappa*2.0d0*0.d0
          alpha(1,2,3,1)=-kappa*2.0d0*0.d0
          alpha(1,2,2,4)=kappa*2.0d0*0.d0
          alpha(1,2,3,4)=kappa*2.0d0*0.d0
          alpha(1,2,4,2)=-kappa*2.0d0
          alpha(1,2,4,3)=-kappa*2.0d0
          do 33 k=1,4
           do 34 l=1,4
            alpha(2,1,k,l)=-alpha(1,2,k,l)
            alpha(2,2,k,l)=alpha(1,1,k,l)
 34        continue
 33       continue
        else
          aConst=kappa*(1.d0+2.d0*dble(n)-vbar)/(3.d0-vbar)
          bConst=kappa*(2.d0+(1.d0-vbar)*dble(n))/(3.d0-vbar)
          alpha(1,1,1,1)=aConst
          alpha(1,1,2,2)=aConst
          alpha(1,1,3,3)=aConst
          alpha(1,1,4,4)=aConst
          alpha(1,1,1,4)=-bConst
          alpha(1,1,2,3)=bConst
          alpha(1,1,3,2)=bConst
          alpha(1,1,4,1)=-bConst

          alpha(1,2,1,2)=aConst
          alpha(1,2,1,3)=bConst
          alpha(1,2,2,1)=-aConst
          alpha(1,2,3,1)=-bConst
          alpha(1,2,2,4)=bConst
          alpha(1,2,3,4)=aConst
          alpha(1,2,4,2)=-bConst
          alpha(1,2,4,3)=-aConst
          do 35 k=1,4
           do 36 l=1,4
            alpha(2,1,k,l)=-alpha(1,2,k,l)
            alpha(2,2,k,l)=alpha(1,1,k,l)
 36        continue
 35       continue
        endif

        elseif (kfold.eq.2) then

        if (n.eq.0) then
          alpha(1,1,1,1)=kappa*0.5d0
          alpha(1,2,1,2)=kappa*0.5d0
          alpha(2,1,2,1)=kappa*0.5d0
          alpha(2,2,2,2)=kappa*0.5d0      
        elseif (n.eq.1) then
          alpha(1,1,1,1)=2.d0*kappa
          alpha(1,1,4,4)=2.d0*kappa

          alpha(1,2,1,2)=kappa*2.0d0
          alpha(1,2,1,3)=kappa*2.0d0
          alpha(1,2,4,2)=-kappa*2.0d0
          alpha(1,2,4,3)=-kappa*2.0d0

          alpha(2,1,2,1)=kappa*2.0d0
          alpha(2,1,3,1)=kappa*2.0d0
          alpha(2,1,2,4)=-kappa*2.0d0
          alpha(2,1,3,4)=-kappa*2.0d0

          alpha(2,2,2,2)=4.d0*kappa
          alpha(2,2,3,3)=4.d0*kappa
          alpha(2,2,2,3)=4.d0*kappa
          alpha(2,2,3,2)=4.d0*kappa
        else
          aConst=kappa*(1.d0+2.d0*dble(n)-vbar)/(3.d0-vbar)
          bConst=kappa*(2.d0+(1.d0-vbar)*dble(n))/(3.d0-vbar)
          alpha(1,1,1,1)=aConst
          alpha(1,1,4,4)=aConst
          alpha(1,1,1,4)=-bConst
          alpha(1,1,4,1)=-bConst

          alpha(1,2,1,2)=aConst
          alpha(1,2,1,3)=bConst
          alpha(1,2,4,2)=-bConst
          alpha(1,2,4,3)=-aConst

          alpha(2,1,2,1)=aConst
          alpha(2,1,3,1)=bConst
          alpha(2,1,2,4)=-bConst
          alpha(2,1,3,4)=-aConst

          alpha(2,2,2,2)=aConst
          alpha(2,2,3,3)=aConst
          alpha(2,2,2,3)=bConst
          alpha(2,2,3,2)=bConst
        endif

        else
          print*,'ERROR: kfold'
          stop
        endif

          do 41 a=1,DtNdofs
           do 42 b=1,DtNdofs
            asymfactor=dsqrt(dble(kfold))
            bsymfactor=dsqrt(dble(kfold))
 
            if (kfold.eq.2) then
              if (node(DtNnode(a),2).lt.1.d-4)
     &          asymfactor=0.5d0*asymfactor
              if (node(DtNnode(b),2).lt.1.d-4)
     &          bsymfactor=0.5d0*bsymfactor
            endif
             
            ThetaA=datan2(node(DtNnode(a),2),node(DtNnode(a),1))
            ThetaB=datan2(node(DtNnode(b),2),node(DtNnode(b),1))
            do 43 k=1,4
             do 44 l=1,4
              KDtN(a,b)=KDtN(a,b)+alpha(DtNdof(a),DtNdof(b),k,l)
     &                  *asymfactor*DtNI(R,ThetaA,dTheta,n,k)
     &                  *bsymfactor*DtNI(R,ThetaB,dTheta,n,l)
 44          continue
 43         continue
 42        continue
 41       continue
        
 11   continue


      return
      end
      
********************************************************************

      function DtNI(R,Theta,dTheta,nFunc,j)
      implicit none
      integer nFunc,j
      real*8 DtNI
      real*8 Theta,dTheta,R,ans,n

      n=dble(nFunc)

      if (nFunc.eq.1) then
        if (j.eq.1) then
          ans=-(R*(-4.d0*dTheta**2.d0+dCos(2.d0*(dTheta-Theta))
     &             -2.d0*dCos(2.d0*Theta)+dCos(2.d0*(dTheta+Theta))))
     &         /(8.d0*dTheta)
        elseif (j.eq.2) then
          ans=(R*dCos(Theta)*dSin(dTheta)**2.d0*dSin(Theta))/dTheta
        elseif (j.eq.3) then
          ans=(R*dCos(Theta)*dSin(dTheta)**2.d0*dSin(Theta))/dTheta
        elseif (j.eq.4) then
          ans=(R*(4.d0*dTheta**2.d0+dCos(2.d0*(dTheta-Theta))
     &             -2.d0*dCos(2.d0*Theta)+dCos(2.d0*(dTheta+Theta))))
     &         /(8.d0*dTheta)
        else
          print*,'  DtN error!: j = 1,2,3,or 4'
          stop
        endif
      else
        if (j.eq.1) then
          ans=(-2.d0*R*(dCos(Theta)*dCos(n*Theta)*(-1.d0-n**2.d0
     &        +(1.d0+n**2.d0)*dCos(dTheta)*dCos(dTheta*n)
     &        +2.d0*n*dSin(dTheta)*dSin(dTheta*n))
     &        +(-2.d0*n+2.d0*n*dCos(dTheta)*dCos(dTheta*n)
     &        +(1.d0+n**2.d0)*dSin(dTheta)*dSin(dTheta*n))
     &        *dSin(Theta)*dSin(n*Theta)))
     &        /(dTheta*(-1.d0+n**2.d0)**2.d0)
        elseif (j.eq.2) then
          ans=(2.d0*R*(dCos(n*Theta)*(1.d0+n**2.d0-(1.d0+n**2.d0)
     &        *dCos(dTheta)*dCos(dTheta*n)
     &        -2.d0*n*dSin(dTheta)*dSin(dTheta*n))
     &        *dSin(Theta)+dCos(Theta)*(-2.d0*n
     &        +2.d0*n*dCos(dTheta)*Cos(dTheta*n)+(1.d0+n**2.d0)
     &        *dSin(dTheta)*dSin(dTheta*n))*Sin(n*Theta)))
     &        /(dTheta*(-1.d0+n**2.d0)**2.d0)
        elseif (j.eq.3) then
          ans=(2.d0*R*(dCos(n*Theta)*(-2.d0*n
     &        +2.d0*n*dCos(dTheta)*dCos(dTheta*n)
     &        +(1.d0+n**2.d0)*dSin(dTheta)*Sin(dTheta*n))
     &        *dSin(Theta)+dCos(Theta)*(1.d0+n**2.d0
     &        -(1.d0+n**2.d0)*dCos(dTheta)*dCos(dTheta*n)
     &        -2.d0*n*dSin(dTheta)*dSin(dTheta*n))*dSin(n*Theta)))
     &        /(dTheta*(-1.d0+n**2.d0)**2.d0)
        elseif (j.eq.4) then
          ans=(-2.d0*R*(dCos(Theta)*dCos(n*Theta)*(-2.d0*n
     &        +2.d0*n*dCos(dTheta)*dCos(dTheta*n)
     &        +(1.d0+n**2.d0)*dSin(dTheta)*dSin(dTheta*n))
     &        +(-1.d0-n**2.d0+(1.d0+n**2.d0)*dCos(dTheta)
     &        *dCos(dTheta*n)+2.d0*n*dSin(dTheta)*dSin(dTheta*n))
     &        *dSin(Theta)*dSin(n*Theta)))
     &        /(dTheta*(-1.d0+n**2.d0)**2.d0)
        else
          print*,'  DtN error!: j = 1,2,3,or 4'
          stop
        endif
      endif
             
      DtNI=ans
      return
      end
      
********************************************************************
