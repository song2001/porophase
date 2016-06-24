!     ------------------------------------------------------------------
!
!     meshabaqus.f
!
!
!     Created on Fri Jun 22 15:10:47 2012
!     Copyright (c) 2012 MyCompany. All rights reserved.
!
!
!     ------------------------------------------------------------------

      program meshgen
      implicit none
      integer maxnode,maxcon
      parameter (maxnode=1000000,maxcon=4)
      real*8 node(maxnode,2),dbcval(maxnode),fbcval(maxnode)
      real*8 conAu(maxnode,maxcon),conBu(maxnode,maxcon)
      real*8 conAvp(maxnode,maxcon),conBvp(maxnode,maxcon)
      real*8 conAmu(maxnode,maxcon),conBmu(maxnode,maxcon)
      real*8 ndistsq(maxnode),lx,ly,mx,my,x,y,h,L,cY1,cY2
      real*8 tempnode(maxnode,2),r,s,data_point,empty,crackloc
      real*8 YM(maxnode),biot(maxnode),phi0(maxnode),poisson(maxnode)
      real*8 freestrain(maxnode),hele,Rmax,rr,fluidsource(maxnode)
      real*8 Imperm(maxnode),data_u,data_mu
      integer i,j,k,ii,nodes,neles,ele(maxnode,9),ndbcs
      integer dbcnd(maxnode),dbcdof(maxnode),numpcs,umpcnd(maxnode,3)
      integer umpc(maxnode,maxcon,2),nmumpcs,mumpc(maxnode,maxcon,2)
      integer problem,order(maxnode),rank(maxnode)
      integer tempele(maxnode,9),mumpcnd(maxnode,3),CrackOn
      integer umunode(maxnode),pnode(maxnode)
      integer nvpmpcs,vpmpc(maxnode,maxcon,2),vpmpcnd(maxnode,3)
      integer isp(maxnode),isumu(maxnode),pcontrol(maxnode,maxcon)
      integer bnode,tnode,lnode,rnode,lbnode,ltnode,rbnode,rtnode
      integer cfound,umunodes,iflow,nbounds,bounds(2,2)
      integer fbcnd(maxnode),fbcnum(maxnode),nfbcs,fbcdof(maxnode)
      integer DtNnode(maxnode*2),DtNdof(maxnode*2),DtNdofs,DtNnodes
      integer gconmu(maxnode,6),gconvp(maxnode,6),gconu(maxnode,6)
      integer n,d,gcon(maxnode,6),numvpu,numtot,num,ndim,ndpn
      integer nconvp,nconmu,nconu

      character(3) fnum
      
 
 22   format(10(i9,2x))
      open(unit=7,file='mesh/nodes_in')
      open(unit=8,file='mesh/elements_in')
      open(unit=17,file='mesh/nodes')
      open(unit=18,file='mesh/elements')
      open(unit=19,file='mesh/dbcs')
      open(unit=21,file='mesh/fbcs')
      open(unit=22,file='mesh/moduli')
      open(unit=23,file='mesh/DtN_dofs')
      open(unit=24,file='mesh/gcon')
      open(unit=25,file='mesh/gconu')
      open(unit=26,file='mesh/gconmu')
      open(unit=27,file='mesh/gconvp')
      
      print*,"-> Reading inputs and sorting nodes..."

      read(7,*) nodes
      do 70 i=1,nodes
        read(7,*)ii,tempnode(i,1),tempnode(i,2)
 70   continue

      read(8,*) neles
      do 80 i=1,neles
        read(8,*)ii,tempele(i,4),tempele(i,3),tempele(i,1),
     1       tempele(i,2),tempele(i,8),tempele(i,6),tempele(i,5),
     2       tempele(i,7)
        nodes=nodes+1
        tempele(i,9)=nodes
        do 81 j=1,2
          tempnode(nodes,j)=(tempnode(tempele(i,1),j)
     &                      +tempnode(tempele(i,2),j)
     &                      +tempnode(tempele(i,3),j)
     &                      +tempnode(tempele(i,4),j))/4.d0
 81     continue
 80   continue

      do 41 i=1,nodes 
!!!!!!! SETS SORTING ORDER !!!!!!!!!!        
        ndistsq(i)=1.d0*((tempnode(i,1)+0.d0)**2)+
     1             1.d6*((tempnode(i,2)+0.d0)**2)
        
        rank(i)=i
        do 71 j=1,i-1
          if (ndistsq(i).lt.ndistsq(j)) then
            if (rank(j).lt.rank(i)) then
              rank(i)=rank(j)
            endif
            rank(j)=rank(j)+1
          endif
 71     continue
 41   continue
            
      do 40 i=1,nodes
        order(rank(i))=i
 40   continue 


      print*,"-> Writing outputs (nodes and elements)..." 
      write(17,*) nodes
      do 170 i=1,nodes
        node(i,1)=tempnode(order(i),1)*1.0d0
        node(i,2)=tempnode(order(i),2)*1.0d0
        write(17,*)i,node(i,1),node(i,2)
 170  continue
      
      write(18,*) neles
      do 180 i=1,neles
        do 181 j=1,9
          ele(i,j)=rank(tempele(i,j))
 181    continue
        write(18,22)i,ele(i,1),ele(i,2),ele(i,3),
     1       ele(i,4),ele(i,5),ele(i,6),ele(i,7),ele(i,8),
     2       ele(i,9)   
 180  continue
 
!     PROBLEM DESCRIPTION
!     Center Crack:          problem = 0
!     Pressurized Crack 1d:  problem = 1
!     1D Crack:              problem = 2
!     Center Crack (no sub): problem = 3 

      problem=0

      print*,"-> Applying mpcs and dbcs..." 
*********** Center Crack *********************************
      
      if (problem.eq.0) then
            
        ndbcs=0
        nfbcs=0
        numpcs=0
        nmumpcs=0
        nvpmpcs=0
        CrackOn=0
        DtNnodes=0

        lx=node(1,1)+1.d-4
        ly=node(1,2)+1.d-4
        mx=node(nodes,1)-1.d-4
        my=node(nodes,2)-1.d-4

        lx=-20.25d0+1.d-4
        ly=-20.25d0+1.d-4
        mx=20.25d0-1.d-4
        my=20.25d0-1.d-4

        Rmax=5.d2-1.d-4

        h=0.5d0
        L=lx+mx

        cY1=-12.5d0
        cY2=12.5d0
        crackloc=5.01d0

        open(unit=30,file='mesh/crackloc')
        read(30,*) crackloc
        close(unit=30)

        crackloc=h

        umunodes=0

        do 680 i=1,nodes
          isp(i)=1
          isumu(i)=0
 680    continue

        do 690 i=1,neles
          if (isumu(ele(i,1)).eq.0) then
            isumu(ele(i,1))=1            
          endif
          if (isumu(ele(i,2)).eq.0) then
            isumu(ele(i,2))=1
          endif
          if (isumu(ele(i,3)).eq.0) then
            isumu(ele(i,3))=1           
          endif
          if (isumu(ele(i,4)).eq.0) then
            isumu(ele(i,4))=1         
          endif
          isp(ele(i,5))=0
          pcontrol(ele(i,5),1)=ele(i,1)
          pcontrol(ele(i,5),2)=ele(i,2)
          isp(ele(i,6))=0
          pcontrol(ele(i,6),1)=ele(i,1)
          pcontrol(ele(i,6),2)=ele(i,3)
          isp(ele(i,7))=0
          pcontrol(ele(i,7),1)=ele(i,2)
          pcontrol(ele(i,7),2)=ele(i,4)
          isp(ele(i,8))=0
          pcontrol(ele(i,8),1)=ele(i,3)
          pcontrol(ele(i,8),2)=ele(i,4)
          isp(ele(i,9))=2
          isumu(ele(i,9))=2
          pcontrol(ele(i,9),1)=ele(i,1)
          pcontrol(ele(i,9),2)=ele(i,2)
          pcontrol(ele(i,9),3)=ele(i,3)
          pcontrol(ele(i,9),4)=ele(i,4)
 690    continue

        do 691 i=1,nodes

         if (isumu(i).eq.2) then

           numpcs=numpcs+1
           umpcnd(numpcs,1)=4
           umpcnd(numpcs,2)=i
           umpcnd(numpcs,3)=1
           umpc(numpcs,1,1)=pcontrol(i,1)
           umpc(numpcs,2,1)=pcontrol(i,2)
           umpc(numpcs,3,1)=pcontrol(i,3)
           umpc(numpcs,4,1)=pcontrol(i,4)
           umpc(numpcs,1,2)=1
           umpc(numpcs,2,2)=1
           umpc(numpcs,3,2)=1
           umpc(numpcs,4,2)=1
           conAu(numpcs,1)=0.25d0
           conBu(numpcs,1)=0.0d0
           conAu(numpcs,2)=0.25d0
           conBu(numpcs,2)=0.0d0
           conAu(numpcs,3)=0.25d0
           conBu(numpcs,3)=0.0d0
           conAu(numpcs,4)=0.25d0
           conBu(numpcs,4)=0.0d0

           numpcs=numpcs+1
           umpcnd(numpcs,1)=4
           umpcnd(numpcs,2)=i
           umpcnd(numpcs,3)=2
           umpc(numpcs,1,1)=pcontrol(i,1)
           umpc(numpcs,2,1)=pcontrol(i,2)
           umpc(numpcs,3,1)=pcontrol(i,3)
           umpc(numpcs,4,1)=pcontrol(i,4)
           umpc(numpcs,1,2)=2
           umpc(numpcs,2,2)=2
           umpc(numpcs,3,2)=2
           umpc(numpcs,4,2)=2
           conAu(numpcs,1)=0.25d0
           conBu(numpcs,1)=0.0d0
           conAu(numpcs,2)=0.25d0
           conBu(numpcs,2)=0.0d0
           conAu(numpcs,3)=0.25d0
           conBu(numpcs,3)=0.0d0
           conAu(numpcs,4)=0.25d0
           conBu(numpcs,4)=0.0d0

           nmumpcs=nmumpcs+1
           mumpcnd(nmumpcs,1)=4
           mumpcnd(nmumpcs,2)=i
           mumpcnd(nmumpcs,3)=6
           mumpc(nmumpcs,1,1)=pcontrol(i,1)
           mumpc(nmumpcs,2,1)=pcontrol(i,2)
           mumpc(nmumpcs,3,1)=pcontrol(i,3)
           mumpc(nmumpcs,4,1)=pcontrol(i,4)
           mumpc(nmumpcs,1,2)=6
           mumpc(nmumpcs,2,2)=6
           mumpc(nmumpcs,3,2)=6
           mumpc(nmumpcs,4,2)=6
           conAmu(nmumpcs,1)=0.25d0
           conBmu(nmumpcs,1)=0.0d0
           conAmu(nmumpcs,2)=0.25d0
           conBmu(nmumpcs,2)=0.0d0
           conAmu(nmumpcs,3)=0.25d0
           conBmu(nmumpcs,3)=0.0d0
           conAmu(nmumpcs,4)=0.25d0
           conBmu(nmumpcs,4)=0.0d0

           isumu(i)=0

         elseif (isumu(i).eq.0) then

           numpcs=numpcs+1
           umpcnd(numpcs,1)=2
           umpcnd(numpcs,2)=i
           umpcnd(numpcs,3)=1
           umpc(numpcs,1,1)=pcontrol(i,1)
           umpc(numpcs,2,1)=pcontrol(i,2)
           umpc(numpcs,1,2)=1
           umpc(numpcs,2,2)=1
           conAu(numpcs,1)=0.5d0
           conBu(numpcs,1)=0.0d0
           conAu(numpcs,2)=0.5d0
           conBu(numpcs,2)=0.0d0

           numpcs=numpcs+1
           umpcnd(numpcs,1)=2
           umpcnd(numpcs,2)=i
           umpcnd(numpcs,3)=2
           umpc(numpcs,1,1)=pcontrol(i,1)
           umpc(numpcs,2,1)=pcontrol(i,2)
           umpc(numpcs,1,2)=2
           umpc(numpcs,2,2)=2
           conAu(numpcs,1)=0.5d0
           conBu(numpcs,1)=0.0d0
           conAu(numpcs,2)=0.5d0
           conBu(numpcs,2)=0.0d0

           nmumpcs=nmumpcs+1
           mumpcnd(nmumpcs,1)=2
           mumpcnd(nmumpcs,2)=i
           mumpcnd(nmumpcs,3)=6
           mumpc(nmumpcs,1,1)=pcontrol(i,1)
           mumpc(nmumpcs,2,1)=pcontrol(i,2)
           mumpc(nmumpcs,1,2)=6
           mumpc(nmumpcs,2,2)=6
           conAmu(nmumpcs,1)=0.5d0
           conBmu(nmumpcs,1)=0.0d0
           conAmu(nmumpcs,2)=0.5d0
           conBmu(nmumpcs,2)=0.0d0
          
         endif

         if (isp(i).eq.0) then !TAYLOR HOOD

           nvpmpcs=nvpmpcs+1
           vpmpcnd(nvpmpcs,1)=2
           vpmpcnd(nvpmpcs,2)=i
           vpmpcnd(nvpmpcs,3)=5
           vpmpc(nvpmpcs,1,1)=pcontrol(i,1)
           vpmpc(nvpmpcs,2,1)=pcontrol(i,2)
           vpmpc(nvpmpcs,1,2)=5
           vpmpc(nvpmpcs,2,2)=5
           conAvp(nvpmpcs,1)=0.5d0
           conBvp(nvpmpcs,1)=0.0d0
           conAvp(nvpmpcs,2)=0.5d0
           conBvp(nvpmpcs,2)=0.0d0

         elseif (isp(i).eq.2) then

           nvpmpcs=nvpmpcs+1
           vpmpcnd(nvpmpcs,1)=4
           vpmpcnd(nvpmpcs,2)=i
           vpmpcnd(nvpmpcs,3)=5
           vpmpc(nvpmpcs,1,1)=pcontrol(i,1)
           vpmpc(nvpmpcs,2,1)=pcontrol(i,2)
           vpmpc(nvpmpcs,3,1)=pcontrol(i,3)
           vpmpc(nvpmpcs,4,1)=pcontrol(i,4)
           vpmpc(nvpmpcs,1,2)=5
           vpmpc(nvpmpcs,2,2)=5
           vpmpc(nvpmpcs,3,2)=5
           vpmpc(nvpmpcs,4,2)=5
           conAvp(nvpmpcs,1)=0.25d0
           conBvp(nvpmpcs,1)=0.0d0
           conAvp(nvpmpcs,2)=0.25d0
           conBvp(nvpmpcs,2)=0.0d0
           conAvp(nvpmpcs,3)=0.25d0
           conBvp(nvpmpcs,3)=0.0d0
           conAvp(nvpmpcs,4)=0.25d0
           conBvp(nvpmpcs,4)=0.0d0
           
           isp(i)=0

         endif

 691    continue

        open(unit=52,file='mesh/sol_diamond22')

        do 692 i=1,nodes
          x=node(i,1)
          y=node(i,2)
          rr=dsqrt(x*x+y*y)
!          read(52,*)empty,empty,empty,data_u,data_mu,empty,empty,empty
          if (isumu(i).eq.1) then
           if (x.lt.lx.or.x.gt.mx) then
            ndbcs=ndbcs+1
            dbcnd(ndbcs)=i
            dbcdof(ndbcs)=1
            dbcval(ndbcs)=0.d0
           endif
!           ndbcs=ndbcs+1
!           dbcnd(ndbcs)=i
!           dbcdof(ndbcs)=2
!           dbcval(ndbcs)=1.d0*data_u/100.d0
           if (y.lt.ly) then
            ndbcs=ndbcs+1
            dbcnd(ndbcs)=i
            dbcdof(ndbcs)=2
            dbcval(ndbcs)=0.d0
           elseif (y.gt.my) then
            ndbcs=ndbcs+1
            dbcnd(ndbcs)=i
            dbcdof(ndbcs)=2
            dbcval(ndbcs)=0.d0
!           elseif (dabs(y).lt.1.d-2) then
!            ndbcs=ndbcs+1
!            dbcnd(ndbcs)=i
!            dbcdof(ndbcs)=2
!            dbcval(ndbcs)=0.0d0
           endif
          endif
          if (x.lt.lx.or.x.gt.mx) then
           ndbcs=ndbcs+1
           dbcnd(ndbcs)=i
           dbcdof(ndbcs)=3
           dbcval(ndbcs)=0.d0
          endif
          if (y.lt.ly+0.0d0.or.y.gt.my-0.0d0) then
           ndbcs=ndbcs+1
           dbcnd(ndbcs)=i
           dbcdof(ndbcs)=4
           dbcval(ndbcs)=0.d0
          endif
!          elseif (dabs(y).lt.1.d-3) then
!           ndbcs=ndbcs+1
!           dbcnd(ndbcs)=i
!           dbcdof(ndbcs)=3
!           dbcval(ndbcs)=0.d0
!          endif    
          if (isp(i).eq.1) then
!           if (x.lt.lx) then
!            ndbcs=ndbcs+1
!            dbcnd(ndbcs)=i
!            dbcdof(ndbcs)=5
!            dbcval(ndbcs)=0.1d0*(5.d0-x)
!           elseif (x.gt.mx) then
!            ndbcs=ndbcs+1
!            dbcnd(ndbcs)=i
!            dbcdof(ndbcs)=5
!            dbcval(ndbcs)=0.d0
!           endif
          endif
          if (isumu(i).eq.1) then
           if (x.lt.-10.24d0) then
            if (dabs(y).lt.0.5d0*h+1.d-2) then
             ndbcs=ndbcs+1
             dbcnd(ndbcs)=i
             dbcdof(ndbcs)=6
             dbcval(ndbcs)=0.d0
            endif
           endif
           if (dabs(x-10.d0).lt.0.5d0*h+1.d-2) then
            if (dabs(y).lt.10.5d0*h+1.d-2) then
             ndbcs=ndbcs+1
             dbcnd(ndbcs)=i
             dbcdof(ndbcs)=6
             dbcval(ndbcs)=0.d0
            endif
           endif
          endif
          if (x.gt.mx) then
!           if (isumu(i).eq.1) then
!            nfbcs=nfbcs+1
!            fbcnd(nfbcs)=i
!            fbcdof(nfbcs)=1
!            fbcval(nfbcs)=0.d0
!            fbcnum(i)=nfbcs
!           endif
!           nfbcs=nfbcs+1
!           fbcnd(nfbcs)=i
!           fbcdof(nfbcs)=3
!           fbcval(nfbcs)=0.d0
!           fbcnum(i)=nfbcs
          elseif (x.lt.lx) then
!           if (isumu(i).eq.1) then
!            nfbcs=nfbcs+1
!            fbcnd(nfbcs)=i
!            fbcdof(nfbcs)=1
!            fbcval(nfbcs)=0.d0
!            fbcnum(i)=nfbcs
!           endif
!           nfbcs=nfbcs+1
!           fbcnd(nfbcs)=i
!           fbcdof(nfbcs)=3
!           fbcval(nfbcs)=0.d0
!           fbcnum(i)=nfbcs
          endif
!          if (y.gt.my) then
!           if (isumu(i).eq.1) then
!            nfbcs=nfbcs+1
!            fbcnd(nfbcs)=i
!            fbcdof(nfbcs)=2
!            fbcval(nfbcs)=0.d0
!            fbcnum(i)=nfbcs
!           endif
!           nfbcs=nfbcs+1
!           fbcnd(nfbcs)=i
!           fbcdof(nfbcs)=4
!           fbcval(nfbcs)=0.d0
!           fbcnum(i)=nfbcs
!          endif

          if (rr.gt.Rmax.and.isumu(i).eq.1) then
            DtNnodes=DtNnodes+1
            DtNdofs=DtNdofs+1
            DtNnode(DtNdofs)=i
            DtNdof(dtNdofs)=1
            if (y.gt.1.d-4) then
              DtNdofs=DtNdofs+1
              DtNnode(DtNdofs)=i
              DtNdof(dtNdofs)=2
            endif
          endif

 692    continue

        close(unit=52)

        do 390 i=1,neles
          nbounds=0
          do 393 k=1,nfbcs
            if (fbcdof(k).eq.1) then
              do 392 j=1,4
                if (ele(i,j).eq.fbcnd(k)) then
                  nbounds=nbounds+1
                  bounds(nbounds,1)=ele(i,j)
                  bounds(nbounds,2)=k
                endif
                if (nbounds.gt.2) then
                  print*,'Too many fbcs in ele: ',i
                  stop
                endif
 392          continue
            endif
 393      continue
          if (nbounds.eq.2) then
            hele=dsqrt((node(bounds(1,1),2)-node(bounds(2,1),2))**2.d0
     &             +(node(bounds(1,1),1)-node(bounds(2,1),1))**2.d0)
            fbcval(bounds(1,2))=fbcval(bounds(1,2))
     &       +hele*0.5d0*dabs(node(bounds(1,1),1))/node(bounds(1,1),1)
            fbcval(bounds(2,2))=fbcval(bounds(2,2))
     &       +hele*0.5d0*dabs(node(bounds(2,1),1))/node(bounds(2,1),1)
          endif

          nbounds=0
          do 493 k=1,nfbcs
            if (fbcdof(k).eq.2) then
              do 492 j=1,4
                if (ele(i,j).eq.fbcnd(k)) then
                  nbounds=nbounds+1
                  bounds(nbounds,1)=ele(i,j)
                  bounds(nbounds,2)=k
                endif
                if (nbounds.gt.2) then
                  print*,'Too many fbcs in ele: ',i
                  stop
                endif
 492          continue
            endif
 493      continue
          if (nbounds.eq.2) then
            hele=dsqrt((node(bounds(1,1),2)-node(bounds(2,1),2))**2.d0
     &             +(node(bounds(1,1),1)-node(bounds(2,1),1))**2.d0)
            fbcval(bounds(1,2))=fbcval(bounds(1,2))
     &       +hele*0.5d0*dabs(node(bounds(1,1),2))/node(bounds(1,1),2)
            fbcval(bounds(2,2))=fbcval(bounds(2,2))
     &       +hele*0.5d0*dabs(node(bounds(2,1),2))/node(bounds(2,1),2)
          endif

          x=0.d0
          y=0.d0
          do 140 j=1,8
            x=x+node(ele(i,j),1)
            y=y+node(ele(i,j),2)
 140      continue
          x=x/8.d0
          y=y/8.d0
!          if (dabs(x-12.5d0).lt.(h*0.499d0).and.y.lt.89.99d0) then
!            YM(i)=1.1230769d-3
!            biot(i)=1.d0
!            phi0(i)=0.3d0
!            poisson(i)=0.4d0
!            freestrain(i)=73.d0*1.d0
!            fluidsource(i)=0.d0
!            Imperm(i)=1.d-3
!          else
            YM(i)=1.d0
            biot(i)=0.d0
            phi0(i)=0.d0
            poisson(i)=0.25d0
            freestrain(i)=0.d0
            fluidsource(i)=0.d0
            Imperm(i)=1.d0
            if (dabs(x+40.d0*h).lt.0.5d0*h+1.d-4) then
             if (dabs(y-0.0d0*h).lt.0.5d0*h+1.d-4) then
              print*,' Injection Element: ',i
              fluidsource(i)=1.d0
             endif
            endif
!            if (dabs(y-0.d0*h).lt.0.5d0*h+1.d-2) then
!             if (dabs(x-30.d0*h).lt.0.5d0*h+1.d-2) then
!              print*,' Injection Element: ',i
!              fluidsource(i)=1.d0
!             endif
!            endif

!            if (dabs(x+10.2d0).lt.0.5d0*h-1.d-4) then
!             if (dabs(y).lt.h*0.5d0-1.d-4) then
!            if (dabs(x).lt.crackloc*0.d0+0.51d0) then
!             if (dabs(y).lt.crackloc*0.d0+0.51d0) then
!              print*,' Injection Element: ',i
!              fluidsource(i)=1.d0
!             elseif (dabs(y).lt.h*0.5d0-1.d-4) then
!              print*,' Injection Element: ',i
!              fluidsource(i)=1.d0
!             elseif (abs(x+cY2).lt.h*0.5d0-1.d-4) then
!              fluidsource(i)=1.d0
!             endif
!            endif

!          endif

 390    continue

*******************************************************************
      
!      write(fnum,905)nx

      endif

      print*,"-> Checking for errors..." 

      do 801 i=1,ndbcs
        do 802 j=1,ndbcs
          if (i.ne.j) then
            if (dbcnd(i).eq.dbcnd(j)) then
              if (dbcdof(i).eq.dbcdof(j)) then
                print*,'ERROR: Repeated dbcs!',' node: ',dbcnd(i),
     1                 dbcdof(i),node(dbcnd(i),1),node(dbcnd(i),2)
              endif
            endif
          endif
 802    continue
        if (dbcdof(i).le.2) then
          do 803 j=1,numpcs
            if (dbcnd(i).eq.umpcnd(j,2)) then
              if (dbcdof(i).eq.umpcnd(j,3)) then
                print*,'ERROR: Cannot have dbc and mpc on same dof!',
     1                 'node: ',dbcnd(i),dbcdof(i),
     2                 node(dbcnd(i),1),node(dbcnd(i),2)
              endif
            endif
 803      continue
        else if (dbcdof(i).le.5) then
          do 804 j=1,nvpmpcs
            if (dbcnd(i).eq.vpmpcnd(j,2)) then
              if (dbcdof(i).eq.vpmpcnd(j,3)) then
                print*,'ERROR: Cannot have dbc and mpc on same dof!',
     1                 'node: ',dbcnd(i),dbcdof(i),
     2                 node(dbcnd(i),1),node(dbcnd(i),2)
              endif
            endif
 804      continue
        else if (dbcdof(i).eq.6) then
          do 805 j=1,nmumpcs
            if (dbcnd(i).eq.mumpcnd(j,2)) then
              if (dbcdof(i).eq.mumpcnd(j,3)) then
                print*,'ERROR: Cannot have dbc and mpc on same dof!',
     1                 'node: ',dbcnd(i),dbcdof(i),
     2                 node(dbcnd(i),1),node(dbcnd(i),2)
              endif
            endif
 805      continue
        else
          print*,'ERROR: Undefined dof #',dbcnd(i),dbcdof(i)
        endif
 801  continue

      do 810 i=1,numpcs
        do 811 j=1,numpcs
          do 812 k=1,umpcnd(i,1)
            if (umpc(i,k,1).eq.umpcnd(j,2)) then
              if (umpc(i,k,2).eq.umpcnd(j,3)) then
                print*,'ERROR: mpcnds cannot be mpcs!',' node: ',
     1                 umpc(i,k,1),umpc(i,k,2),
     2                 node(umpc(i,k,1),1),node(umpc(i,k,1),2)
              endif
            endif
 812      continue
 811    continue
 810  continue

      do 813 i=1,nvpmpcs
        do 814 j=1,nvpmpcs
          do 815 k=1,vpmpcnd(i,1)
            if (vpmpc(i,k,1).eq.vpmpcnd(j,2)) then
              if (vpmpc(i,k,2).eq.vpmpcnd(j,3)) then
                print*,'ERROR: mpcnds cannot be mpcs!',' node: ',
     1                 vpmpc(i,k,1),vpmpc(i,k,2),
     2                 node(vpmpc(i,k,1),1),node(vpmpc(i,k,1),2)

              endif
            endif
 815      continue
 814    continue
 813  continue

      do 816 i=1,nmumpcs
        do 817 j=1,nmumpcs
          do 818 k=1,mumpcnd(i,1)
            if (mumpc(i,k,1).eq.mumpcnd(j,2)) then
              if (mumpc(i,k,2).eq.mumpcnd(j,3)) then
                print*,'ERROR: mpcnds cannot be mpcs!',' node: ',
     1                 mumpc(i,k,1),mumpc(i,k,2),
     2                 node(mumpc(i,k,1),1),node(mumpc(i,k,1),2)

              endif
            endif
 818      continue
 817    continue
 816  continue

      print*,"-> Writing outputs (vtk, dbcs, fbcs)..." 

         
      write(19,*) ndbcs,numpcs,nvpmpcs,nmumpcs
      do 91 i=1,ndbcs
        write(19,*)dbcnd(i),dbcdof(i),dbcval(i)
 91   continue
      do 95 i=1,numpcs
        write(19,*)umpcnd(i,1),umpcnd(i,2),umpcnd(i,3)
        do 96 j=1,umpcnd(i,1)
          write(19,*)umpc(i,j,1),umpc(i,j,2),conAu(i,j),conBu(i,j)
 96     continue
 95   continue
      do 395 i=1,nvpmpcs
        write(19,*)vpmpcnd(i,1),vpmpcnd(i,2),vpmpcnd(i,3)
        do 396 j=1,vpmpcnd(i,1)
          write(19,*)vpmpc(i,j,1),vpmpc(i,j,2),conAvp(i,j),conBvp(i,j)
 396    continue
 395  continue
      do 295 i=1,nmumpcs
        write(19,*)mumpcnd(i,1),mumpcnd(i,2),mumpcnd(i,3)
        do 296 j=1,mumpcnd(i,1)
          write(19,*)mumpc(i,j,1),mumpc(i,j,2),conAmu(i,j),conBmu(i,j)
 296    continue
 295  continue

      write(21,*) nfbcs
      do 391 i=1,nfbcs
        write(21,*)fbcnd(i),fbcdof(i),fbcval(i)
 391  continue

      do 495 i=1,neles
        write(22,*)i,YM(i),biot(i),phi0(i),poisson(i),freestrain(i),
     &             fluidsource(i),Imperm(i),0.d0,0.d0,Imperm(i)
 495  continue

      if (problem.eq.3.or.problem.eq.4) then
       write(23,*) DtNnodes,DtNdofs
       do 496 i=1,DtNdofs
         write(23,*)DtNnode(i),DtNdof(i)
 496   continue
      endif
       
      open(unit=20,file='mesh/mesh.vtk')
      write(20,900) '# vtk DataFile Version 3.0'
      write(20,900) 'Mesh'
      write(20,900) 'ASCII'
      write(20,900) 'DATASET UNSTRUCTURED_GRID'
      write(20,901) 'POINTS ',nodes,' float'
      do 201 i=1,nodes
        write(20,902)node(i,1),node(i,2),0.0
 201  continue
      write(20,*)
      write(20,903) 'CELLS ',neles,neles*9
      do 202 i=1,neles
        write(20,*)8,ele(i,1)-1,ele(i,2)-1,ele(i,4)-1,ele(i,3)-1,
     2    ele(i,5)-1,ele(i,7)-1,ele(i,8)-1,ele(i,6)-1
 202  continue
      write(20,*)
      write(20,904) 'CELL_TYPES ',neles
      do 203 i=1,neles
        write(20,*)23
 203  continue

      print*,'-> Constructing Connectivity'
      ndim=2
      ndpn=6
      nconvp=nvpmpcs
      nconu=numpcs
      nconmu=nmumpcs

      do 60 i=1,nodes
        gconmu(i,1)=i
        do 61 j=1,ndim
          gconu(i,j)=(i-1)*ndim+j
          gconvp(i,j)=(i-1)*(ndim+1)+j
 61     continue
        gconvp(i,ndim+1)=i*(ndim+1)
        do 63 j=1,ndpn
          gcon(i,j)=(i-1)*(ndpn)+j
 63     continue
 60   continue

      print*,'...dbcs'   
      do 400 i=1,ndbcs
        n=dbcnd(i)
        d=dbcdof(i)
        if (d.le.ndim) then
          num=gconu(n,d)
          numtot=gcon(n,d)
          do 410 j=1,nodes
            do 420 k=1,ndim
              if(gconu(j,k).gt.num) gconu(j,k)=gconu(j,k)-1
 420        continue
            do 460 k=1,ndpn
              if(gcon(j,k).gt.numtot) gcon(j,k)=gcon(j,k)-1
 460        continue
 410      continue
          gconu(n,d)=nodes*ndim
          gcon(n,d)=nodes*(ndpn)
        else
          if (d.eq.ndpn) then
            num=gconmu(n,1)
            numtot=gcon(n,d)
            do 411 j=1,nodes
              if(gconmu(j,1).gt.num) gconmu(j,1)=gconmu(j,1)-1
              do 461 k=1,ndpn
                if(gcon(j,k).gt.numtot) gcon(j,k)=gcon(j,k)-1
 461          continue
 411        continue
            gconmu(n,1)=nodes
            gcon(n,d)=nodes*(ndpn)
          else
            num=gconvp(n,d-ndim)
            numtot=gcon(n,d)
            do 412 j=1,nodes
              do 413 k=1,ndim+1
                if(gconvp(j,k).gt.num) gconvp(j,k)=gconvp(j,k)-1
 413          continue
              do 462 k=1,ndpn
                if(gcon(j,k).gt.numtot) gcon(j,k)=gcon(j,k)-1
 462          continue
 412        continue
            gconvp(n,d-ndim)=nodes*(ndim+1)
            gcon(n,d)=nodes*(ndpn)
          endif
        endif
 400  continue
      
      print*,'...umpcs'
      do 425 i=1,nconu
        num=gconu(umpcnd(i,1+1),umpcnd(i,2+1))
        numtot=gcon(umpcnd(i,1+1),umpcnd(i,2+1))
        do 427 j=1,nodes
          do 428 k=1,ndim
            if(gconu(j,k).gt.num) gconu(j,k)=gconu(j,k)-1
 428      continue
          do 463 k=1,ndpn
            if(gcon(j,k).gt.numtot) gcon(j,k)=gcon(j,k)-1
 463      continue
 427    continue
        gconu(umpcnd(i,1+1),umpcnd(i,2+1))=nodes*ndim
        gcon(umpcnd(i,1+1),umpcnd(i,2+1))=nodes*(ndpn)
 425  continue

      print*,'...vpmpcs'
      do 445 i=1,nconvp
        num=gconvp(vpmpcnd(i,1+1),vpmpcnd(i,2+1)-ndim)
        numtot=gcon(vpmpcnd(i,1+1),vpmpcnd(i,2+1))
        do 449 j=1,nodes
          do 450 k=1,ndim+1
            if(gconvp(j,k).gt.num) gconvp(j,k)=gconvp(j,k)-1
 450      continue
          do 464 k=1,ndpn
            if(gcon(j,k).gt.numtot) gcon(j,k)=gcon(j,k)-1
 464      continue
 449    continue
        gconvp(vpmpcnd(i,1+1),vpmpcnd(i,2+1)-ndim)=nodes*(ndim+1)
        gcon(vpmpcnd(i,1+1),vpmpcnd(i,2+1))=nodes*(ndpn)
 445  continue

      print*,'...mumpcs'
      do 435 i=1,nconmu
        num=gconmu(mumpcnd(i,1+1),1)
        numtot=gcon(mumpcnd(i,1+1),mumpcnd(i,2+1))
        do 439 j=1,nodes
          if((gconmu(j,1)).gt.num) gconmu(j,1)=gconmu(j,1)-1
          do 465 k=1,ndpn
            if(gcon(j,k).gt.numtot) gcon(j,k)=gcon(j,k)-1
 465      continue
 439    continue
        gconmu(mumpcnd(i,1+1),1)=nodes
        gcon(mumpcnd(i,1+1),mumpcnd(i,2+1))=nodes*(ndpn)
 435  continue

      do 49 i=1,nodes
         write(24,*)gcon(i,1),gcon(i,2),gcon(i,3),gcon(i,4),
     &              gcon(i,5),gcon(i,6)
         write(25,*)gconu(i,1),gconu(i,2)
         write(26,*)gconmu(i,1)
         write(27,*)gconvp(i,1),gconvp(i,2),gconvp(i,3)
 49   continue
      

      close(unit=7)
      close(unit=8)
      close(unit=17)
      close(unit=18)
      close(unit=19)
      close(unit=20)
      close(unit=21)
      close(unit=22)
      close(unit=23)
      close(unit=24)
      close(unit=25)
      close(unit=26)
      close(unit=27)

      print*,"-> MESH GENERATED!" 
 
 900  format(A)
 901  format(A, I8, A)
 902  format(F15.5, F15.5, F15.5)
 903  format(A, I8, ' ', I8)
 904  format(A, I8)
 905  format(I3.3)
      stop
      end
