      program newpdf
      implicit real*8(a-h,o-z)
      parameter(nspec_max=10,num_ion_max=200)
      character*2 xa_name(nspec_max),xb_name(nspec_max)
      character*1     xx_name(nspec_max)
      dimension num_ion_a(nspec_max),num_ion_b(nspec_max),
     $     num_ion_x(nspec_max)
      dimension bscat_a(nspec_max),bscat_b(nspec_max),
     $     bscat_x(nspec_max)
      dimension therm_a(nspec_max),therm_b(nspec_max),
     $     therm_x(nspec_max)
      dimension val_r_a(nspec_max,nspec_max),
     $     val_r_b(nspec_max,nspec_max),
     $     val_N_a(nspec_max,nspec_max),val_N_b(nspec_max,nspec_max)
      dimension size(3)
      dimension posion_a(3,num_ion_max,nspec_max),
     $     posion_b(3,num_ion_max,nspec_max),
     $     posion_x(3,num_ion_max,nspec_max)
      

      call read_input(num_uc,num_spec_a,num_spec_b,num_spec_x,
     $     num_ion_a,num_ion_b,num_ion_x,xa_name,xb_name,xx_name,
     $     therm_a,therm_b,therm_x,bscat_a,bscat_b,bscat_x,
     $     posion_a,posion_b,posion_x,size,val_r_a,val_r_b,
     $     val_N_a,val_N_b)



      call get_pdf(num_uc,num_spec_a,num_spec_b,num_spec_x,
     $     num_ion_a,num_ion_b,num_ion_x,xa_name,xb_name,xx_name,
     $     bscat_a,bscat_b,bscat_x,therm_a,therm_b,therm_x,
     $     posion_a,posion_b,posion_x,size)

      call ao_analyze(num_spec_a,num_spec_x,
     $     num_ion_a,num_ion_x,xa_name,xx_name,
     $     val_r_a,val_N_a,
     $     posion_a,posion_b,posion_x,size)

      call bo_analyze(num_spec_b,num_spec_x,
     $     num_ion_b,num_ion_x,xb_name,xx_name,
     $     val_r_b,val_N_b,
     $     posion_a,posion_b,posion_x,size)

      call oxygen_analyze(num_spec_a,num_spec_b,num_spec_x,
     $     num_ion_a,num_ion_b,num_ion_x,xa_name,xb_name,xx_name,
     $     val_r_a,val_r_b,val_N_a,val_N_b,
     $     posion_a,posion_b,posion_x,size)



      end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine read_input(num_uc,num_spec_a,num_spec_b,num_spec_x,
     $     num_ion_a,num_ion_b,num_ion_x,xa_name,xb_name,xx_name,
     $     therm_a,therm_b,therm_x,bscat_a,bscat_b,bscat_x,
     $     posion_a,posion_b,posion_x,size,val_r_a,val_r_b,
     $     val_N_a,val_N_b)
      implicit real*8(a-h,o-z)
      parameter(nspec_max=10,num_ion_max=200)
      character*2 xa_name(nspec_max),xb_name(nspec_max)
      character*1     xx_name(nspec_max)
      dimension num_ion_a(nspec_max),num_ion_b(nspec_max),
     $     num_ion_x(nspec_max)
      dimension bscat_a(nspec_max),bscat_b(nspec_max),
     $     bscat_x(nspec_max)
      dimension therm_a(nspec_max),therm_b(nspec_max),
     $     therm_x(nspec_max)
      dimension val_r_a(nspec_max,nspec_max),
     $     val_r_b(nspec_max,nspec_max),
     $     val_N_a(nspec_max,nspec_max),val_N_b(nspec_max,nspec_max)
      dimension size(3)
      dimension posion_a(3,num_ion_max,nspec_max),
     $     posion_b(3,num_ion_max,nspec_max),
     $     posion_x(3,num_ion_max,nspec_max)

      open(99,file='Input.dat')
      open(100,file='Echo.dat')
      read(99,*) num_uc,num_spec_a,num_spec_b,num_spec_x
      write(100,*) num_uc,num_spec_a,num_spec_b,num_spec_x
      nsum_a=0
      do iaspec=1,num_spec_a
         read(99,*) num_ion_a(iaspec),xa_name(iaspec),
     $        bscat_a(iaspec),therm_a(iaspec)
         write(100,*) num_ion_a(iaspec),xa_name(iaspec),
     $        bscat_a(iaspec),therm_a(iaspec)
c         call flush(100,iostat)
         nsum_a=nsum_a+num_ion_a(iaspec)

      enddo
      if(nsum_a.gt.num_uc) then
         write(100,*) "Number of A-sites greater than 
     $        Number of unit cells",nsum_a,num_uc
         stop
      endif

      nsum_b=0
      do ibspec=1,num_spec_b
         read(99,*) num_ion_b(ibspec),xb_name(ibspec),
     $        bscat_b(ibspec),therm_b(ibspec)
         write(100,*) num_ion_b(ibspec),xb_name(ibspec),
     $        bscat_b(ibspec),therm_b(ibspec)
c         call flush(100,iostat)
         nsum_b=nsum_b+num_ion_b(ibspec)
      enddo
      if(nsum_b.gt.num_uc) then
         write(100,*) "Number of B-sites greater than 
     $        Number of unit cells"
         stop
      endif


      nsum_x=0
      do ixspec=1,num_spec_x
         read(99,*) num_ion_x(ixspec),xx_name(ixspec),
     $        bscat_x(ixspec),therm_x(ixspec)
         write(100,*) num_ion_x(ixspec),xx_name(ixspec),
     $        bscat_x(ixspec),therm_x(ixspec)
c         call flush(100,iostat)

         nsum_x=nsum_x+num_ion_x(ixspec)
      enddo
      if(nsum_x.gt.num_uc*3) then
         write(100,*) "Number of X-sites greater than 
     $        3*Number of unit cells"
         stop
      endif

      do iaspec=1,num_spec_a
         do ixspec=1,num_spec_x
            read(99,*) val_r_a(iaspec,ixspec),
     $           val_N_a(iaspec,ixspec)
         write(100,* ) val_r_a(iaspec,ixspec),
     $           val_N_a(iaspec,ixspec)
c         call flush(100,iostat)

         enddo
      enddo

      do ibspec=1,num_spec_b
         do ixspec=1,num_spec_x
            read(99,*) val_r_b(ibspec,ixspec),
     $           val_N_b(ibspec,ixspec)
            write(100,*) val_r_b(ibspec,ixspec),
     $           val_N_b(ibspec,ixspec)
c            call flush(100,iostat)
         enddo
      enddo

      read(99,*) size(1),size(2),size(3)
      do iaspec=1,num_spec_a
         do i_ion=1,num_ion_a(iaspec)
            read(99,*) posion_a(1,i_ion,iaspec),
     $           posion_a(2,i_ion,iaspec),posion_a(3,i_ion,iaspec)
            do idir=1,3
               if(posion_a(idir,i_ion,iaspec).gt.1.0) 
     $              posion_a(idir,i_ion,iaspec)=
     $              posion_a(idir,i_ion,iaspec)-1.0
               if(posion_a(idir,i_ion,iaspec).lt.-1.0) 
     $              posion_a(idir,i_ion,iaspec)=
     $              posion_a(idir,i_ion,iaspec)+1.0
               posion_a(idir,i_ion,iaspec)=
     $              posion_a(idir,i_ion,iaspec)*size(idir)
            enddo
            write(100,*) posion_a(1,i_ion,iaspec),
     $           posion_a(2,i_ion,iaspec),posion_a(3,i_ion,iaspec)
c            call flush(100,iostat)
           
         enddo
      enddo
      

      do ibspec=1,num_spec_b
         do i_ion=1,num_ion_b(ibspec)
            read(99,*) posion_b(1,i_ion,ibspec),
     $           posion_b(2,i_ion,ibspec),posion_b(3,i_ion,ibspec)
            do idir=1,3
               if(posion_b(idir,i_ion,ibspec).gt.1.0) 
     $              posion_b(idir,i_ion,ibspec)=
     $              posion_b(idir,i_ion,ibspec)-1.0
               if(posion_b(idir,i_ion,ibspec).lt.-1.0) 
     $              posion_b(idir,i_ion,ibspec)=
     $              posion_b(idir,i_ion,ibspec)+1.0
               posion_b(idir,i_ion,ibspec)=
     $              posion_b(idir,i_ion,ibspec)*size(idir)
            enddo
            write(100,*) posion_b(1,i_ion,ibspec),
     $           posion_b(2,i_ion,ibspec),posion_b(3,i_ion,ibspec)
c            call flush(100,iostat)

         enddo
      enddo


      do ixspec=1,num_spec_x
         do i_ion=1,num_ion_x(ixspec)
            read(99,*) posion_x(1,i_ion,ixspec),
     $           posion_x(2,i_ion,ixspec),posion_x(3,i_ion,ixspec)
            do idir=1,3
               if(posion_x(idir,i_ion,ixspec).gt.1.0) 
     $              posion_x(idir,i_ion,ixspec)=
     $              posion_x(idir,i_ion,ixspec)-1.0
               if(posion_x(idir,i_ion,ixspec).lt.-1.0) 
     $              posion_x(idir,i_ion,ixspec)=
     $              posion_x(idir,i_ion,ixspec)+1.0
               posion_x(idir,i_ion,ixspec)=
     $              posion_x(idir,i_ion,ixspec)*size(idir)
            enddo
            write(100,*) posion_x(1,i_ion,ixspec),
     $           posion_x(2,i_ion,ixspec),posion_x(3,i_ion,ixspec)
c            call flush(100,iostat)

         enddo
      enddo

C     Cput in a loop to put all position s betwee 0 and 1

      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine get_pdf(num_uc,num_spec_a,num_spec_b,num_spec_x,
     $     num_ion_a,num_ion_b,num_ion_x,xa_name,xb_name,xx_name,
     $     bscat_a,bscat_b,bscat_x,therm_a,therm_b,therm_x,
     $     posion_a,posion_b,posion_x,size)
      implicit real*8(a-h,o-z)
      parameter(nspec_max=10,num_ion_max=200)
      character*2 xa_name(nspec_max),xb_name(nspec_max),xa2
      character*1     xx_name(nspec_max)
      dimension num_ion_a(nspec_max),num_ion_b(nspec_max),
     $     num_ion_x(nspec_max)
      dimension bscat_a(nspec_max),bscat_b(nspec_max),
     $     bscat_x(nspec_max)
      dimension therm_a(nspec_max),therm_b(nspec_max),
     $     therm_x(nspec_max)
      dimension val_r_a(nspec_max,nspec_max),
     $     val_r_b(nspec_max,nspec_max),
     $     val_N_a(nspec_max,nspec_max),val_N_b(nspec_max,nspec_max)
      dimension size(3)
      character*4 xa4
      dimension posion_a(3,num_ion_max,nspec_max),
     $     posion_b(3,num_ion_max,nspec_max),
     $     posion_x(3,num_ion_max,nspec_max)
      parameter(numpdf=1000)
      dimension pdftot(numpdf),apdftot(numpdf),apdf(numpdf),
     $     total_pdf(numpdf)
      dimension partial_pdf_aa(numpdf,nspec_max,nspec_max),
     $     partial_pdf_ab(numpdf,nspec_max,nspec_max),
     $     partial_pdf_ax(numpdf,nspec_max,nspec_max),
     $     partial_pdf_bb(numpdf,nspec_max,nspec_max),
     $     partial_pdf_bx(numpdf,nspec_max,nspec_max),
     $     partial_pdf_xx(numpdf,nspec_max,nspec_max)
      dimension posion_i(3,num_ion_max),posion_j(3,num_ion_max)
      
      sizemesh=0.01

      do idist=1,1000
         pdftot(idist)=0.0
      enddo

c     calculate A-A partial pdfs
      do iaspec=1,num_spec_a
         do jaspec=1,num_spec_a
            do m=1,numpdf
               apdf(m)=0
            enddo
            
            itype=1
            jtype=1
            
            call prepare_arrays(itype,jtype,iaspec,jaspec,
     $           posion_i,posion_j,posion_a,posion_b,posion_x,
     $           num_ion_i,num_ion_j,num_ion_a,num_ion_b,num_ion_x,
     $           therm_i,therm_j,therm_a,therm_b,therm_x)

            call  evaluate_partial_pdf(posion_i,posion_j,
     $           num_ion_i,num_ion_j,therm_i,therm_j,size,apdf)

            open(111,
     $           file=xa_name(iaspec)//xa_name(jaspec)//'pdf')

c     partial pdf is apdf times the neutron scattering factor bscat
            do m=1,numpdf
               r=sizemesh*m
               partial_pdf_aa(m,jaspec,iaspec)=apdf(m)
     $              *sqrt(bscat_a(iaspec)*bscat_a(jaspec))
               write(111,*) r,
     $              partial_pdf_aa(m,jaspec,iaspec)
               total_pdf(m)=total_pdf(m)+partial_pdf_aa(m,jaspec,iaspec)
            enddo

         enddo
      enddo

c     calculate A-B partial pdfs
      do iaspec=1,num_spec_a
         do jbspec=1,num_spec_b
            do m=1,1000
               apdf(m)=0
            enddo
            
            itype=1
            jtype=2
            
            call prepare_arrays(itype,jtype,iaspec,jbspec,
     $           posion_i,posion_j,posion_a,posion_b,posion_x,
     $           num_ion_i,num_ion_j,num_ion_a,num_ion_b,num_ion_x,
     $           therm_i,therm_j,therm_a,therm_b,therm_x)

            call  evaluate_partial_pdf(posion_i,posion_j,
     $           num_ion_i,num_ion_j,therm_i,therm_j,size,apdf)
            open(111,
     $           file=xa_name(iaspec)//xb_name(jbspec)//'pdf')
            do m=1,numpdf
               r=sizemesh*m
               partial_pdf_ab(m,jbspec,iaspec)=apdf(m)
     $              *sqrt(bscat_a(iaspec)*bscat_b(jbspec))
               write(111,*) r, partial_pdf_ab(m,jbspec,iaspec)
               total_pdf(m)=total_pdf(m)+partial_pdf_ab(m,jbspec,iaspec)
            enddo
            close(111)
         enddo
      enddo
      
c     calculate A-X partial pdfs
      do iaspec=1,num_spec_a
         do jxspec=1,num_spec_x
            do m=1,1000
               apdf(m)=0
            enddo
            
            itype=1
            jtype=3
            
            call prepare_arrays(itype,jtype,iaspec,jxspec,
     $           posion_i,posion_j,posion_a,posion_b,posion_x,
     $           num_ion_i,num_ion_j,num_ion_a,num_ion_b,num_ion_x,
     $           therm_i,therm_j,therm_a,therm_b,therm_x)

            call  evaluate_partial_pdf(posion_i,posion_j,
     $           num_ion_i,num_ion_j,therm_i,therm_j,size,apdf)
            open(111,
     $           file=xa_name(iaspec)//xx_name(jxspec)//'pdf')
            do m=1,numpdf
               r=m*sizemesh
               partial_pdf_ax(m,jxspec,iaspec)=apdf(m)
     $              *sqrt(bscat_a(iaspec)*bscat_x(jxspec))
               write(111,*) r,partial_pdf_ax(m,jxspec,iaspec)
               total_pdf(m)=total_pdf(m)+partial_pdf_ax(m,jxspec,iaspec)
            enddo
            close(111)
         enddo
      enddo


c     calculate B-B partial pdfs
      do ibspec=1,num_spec_b
         do jbspec=ibspec,num_spec_b
            do m=1,1000
               apdf(m)=0
            enddo
            
            itype=2
            jtype=2
            
            call prepare_arrays(itype,jtype,ibspec,jbspec,
     $           posion_i,posion_j,posion_a,posion_b,posion_x,
     $           num_ion_i,num_ion_j,num_ion_a,num_ion_b,num_ion_x,
     $           therm_i,therm_j,therm_a,therm_b,therm_x)

            call  evaluate_partial_pdf(posion_i,posion_j,
     $           num_ion_i,num_ion_j,therm_i,therm_j,size,apdf)

            open(111,
     $           file=xb_name(ibspec)//xb_name(jbspec)//'pdf')
            do m=1,numpdf
               r=sizemesh*m 
               partial_pdf_aa(m,jbspec,ibspec)=apdf(m)
     $              *sqrt(bscat_b(ibspec)*bscat_b(jbspec))
               write(111,*) r,partial_pdf_bb(m,jbspec,ibspec)
               total_pdf(m)=total_pdf(m)+partial_pdf_bb(m,jbspec,ibspec)
            enddo
            close(111)
         enddo
      enddo

c     calculate B-X partial pdfs
      do ibspec=1,num_spec_b
         do jxspec=1,num_spec_x
            do m=1,1000
               apdf(m)=0
            enddo
            
            itype=2
            jtype=3
            
            call prepare_arrays(itype,jtype,ibspec,jxspec,
     $           posion_i,posion_j,posion_a,posion_b,posion_x,
     $           num_ion_i,num_ion_j,num_ion_a,num_ion_b,num_ion_x,
     $           therm_i,therm_j,therm_a,therm_b,therm_x)

            call  evaluate_partial_pdf(posion_i,posion_j,
     $           num_ion_i,num_ion_j,therm_i,therm_j,size,apdf)

            open(111,
     $           file=xb_name(ibspec)//xx_name(jxspec)//'pdf')
            do m=1,numpdf
               r=sizemesh*m
               partial_pdf_bx(m,jbspec,ibspec)=apdf(m)
     $              *sqrt(bscat_b(ibspec)*bscat_x(jxspec))
               write(111,*) r,partial_pdf_bx(m,jxspec,ibspec)
               total_pdf(m)=total_pdf(m)+partial_pdf_bx(m,jxspec,ibspec)
            enddo
            close(111)
         enddo
      enddo

c     calculate X-X partial pdfs
      do ixspec=1,num_spec_x
         do jxspec=1,num_spec_x
            do m=1,1000
               apdf(m)=0
            enddo
            
            itype=3
            jtype=3
            
            call prepare_arrays(itype,jtype,ixspec,jxspec,
     $           posion_i,posion_j,posion_a,posion_b,posion_x,
     $           num_ion_i,num_ion_j,num_ion_a,num_ion_b,num_ion_x,
     $           therm_i,therm_j,therm_a,therm_b,therm_x)

            call  evaluate_partial_pdf(posion_i,posion_j,
     $           num_ion_i,num_ion_j,therm_i,therm_j,size,apdf)

            open(111,
     $           file=xx_name(ixspec)//xx_name(jxspec)//'pdf')
            do m=1,numpdf
               r=sizemesh*m
               partial_pdf_xx(m,jbspec,ibspec)=apdf(m)
     $              *sqrt(bscat_x(ixspec)*bscat_x(jxspec))
               write(111,*) r,partial_pdf_xx(m,jxspec,ixspec)
               total_pdf(m)=total_pdf(m)+partial_pdf_xx(m,jxspec,ixspec)
            enddo
            close(111)
         enddo
      enddo

c     print out the total pdf
      open(100,file='Totalpdf')
      do m=1,1000,10
         r=m*sizemesh
         write(100,*) r,total_pdf(m),m
      enddo
    
      end


      subroutine prepare_arrays(itype,jtype,ispec,jspec,
     $     posion_i,posion_j,posion_a,posion_b,posion_x,
     $           num_ion_i,num_ion_j,num_ion_a,num_ion_b,num_ion_x,
     $           therm_i,therm_j,therm_a,therm_b,therm_x)
      implicit real*8(a-h,o-z)
      parameter(nspec_max=10,num_ion_max=200)
      dimension num_ion_a(nspec_max),num_ion_b(nspec_max),
     $     num_ion_x(nspec_max)
      dimension therm_a(nspec_max),therm_b(nspec_max),
     $     therm_x(nspec_max)
      dimension posion_a(3,num_ion_max,nspec_max),
     $     posion_b(3,num_ion_max,nspec_max),
     $     posion_x(3,num_ion_max,nspec_max)
      dimension posion_i(3,num_ion_max),posion_j(3,num_ion_max)

c     copy the appropriate arrays for the positions of the two elements used in 
c     partial pdf computation 
      if(itype.eq.1.and.jtype.eq.1) then
         therm_i=therm_a(ispec)
         therm_j=therm_a(jspec)
         num_ion_i=num_ion_a(ispec)
         num_ion_j=num_ion_a(jspec)
         do i_ion=1,num_ion_i
            do idir=1,3
               posion_i(idir,i_ion)=posion_a(idir,i_ion,ispec)
            enddo
         enddo
         do j_ion=1,num_ion_j
            do idir=1,3
               posion_j(idir,j_ion)=posion_a(idir,j_ion,jspec)
            enddo
         enddo
      elseif(itype.eq.1.and.jtype.eq.2) then
         therm_i=therm_a(ispec)
         therm_j=therm_b(jspec)
         num_ion_i=num_ion_a(ispec)
         num_ion_j=num_ion_b(jspec)
         do i_ion=1,num_ion_i
            do idir=1,3
               posion_i(idir,i_ion)=posion_a(idir,i_ion,ispec)
            enddo
         enddo
         do j_ion=1,num_ion_j
            do idir=1,3
               posion_j(idir,j_ion)=posion_b(idir,j_ion,jspec)
            enddo
         enddo
      elseif(itype.eq.1.and.jtype.eq.3) then
         therm_i=therm_a(ispec)
         therm_j=therm_x(jspec)
         num_ion_i=num_ion_a(ispec)
         num_ion_j=num_ion_x(jspec)
         do i_ion=1,num_ion_i
            do idir=1,3
               posion_i(idir,i_ion)=posion_a(idir,i_ion,ispec)
            enddo
         enddo
         do j_ion=1,num_ion_j
            do idir=1,3
               posion_j(idir,j_ion)=posion_x(idir,j_ion,jspec)
            enddo
         enddo
      elseif(itype.eq.2.and.jtype.eq.2) then
         therm_i=therm_b(ispec)
         therm_j=therm_b(jspec)
         num_ion_i=num_ion_b(ispec)
         num_ion_j=num_ion_b(jspec)
         do i_ion=1,num_ion_i
            do idir=1,3
               posion_i(idir,i_ion)=posion_b(idir,i_ion,ispec)
            enddo
         enddo
         do j_ion=1,num_ion_j
            do idir=1,3
               posion_j(idir,j_ion)=posion_b(idir,j_ion,jspec)
            enddo
         enddo     
      elseif(itype.eq.2.and.jtype.eq.3) then
         therm_i=therm_b(ispec)
         therm_j=therm_x(jspec)
         num_ion_i=num_ion_b(ispec)
         num_ion_j=num_ion_x(jspec)
         do i_ion=1,num_ion_i
            do idir=1,3
               posion_i(idir,i_ion)=posion_b(idir,i_ion,ispec)
            enddo
         enddo
         do j_ion=1,num_ion_j
            do idir=1,3
               posion_j(idir,j_ion)=posion_x(idir,j_ion,jspec)
            enddo
         enddo
      elseif(itype.eq.3.and.jtype.eq.3) then
         therm_i=therm_x(ispec)
         therm_j=therm_x(jspec)
         num_ion_i=num_ion_x(ispec)
         num_ion_j=num_ion_x(jspec)
         do i_ion=1,num_ion_i
            do idir=1,3
               posion_i(idir,i_ion)=posion_x(idir,i_ion,ispec)
            enddo
         enddo
         do j_ion=1,num_ion_j
            do idir=1,3
               posion_j(idir,j_ion)=posion_x(idir,j_ion,jspec)
            enddo
         enddo
      endif

         


      end


      subroutine evaluate_partial_pdf(posion_i,posion_j,
     $     num_ion_i,num_ion_j,therm_i,therm_j,size,apdf)

      implicit real*8(a-h,o-z)
      parameter(nspec_max=10,num_ion_max=200)
      parameter(numpdf=1000)
      dimension apdf(numpdf),size(3),posion_i(3,num_ion_max),
     $     posion_j(3,num_ion_max),distc(3),add(3)

c     loop over all ions 
      sizemesh=0.01
      rcut=10
      w=(therm_i+therm_j)/2
      do i_ion=1,num_ion_i
         do j_ion=1,num_ion_j
            
c  this takes care of the periodic images of the ions  in the 26 cells 
c  bordering the unit cell in the input file
            do ix=-1,1
               add(1)=ix*size(1)
               do iy=-1,1
                  add(2)=iy*size(2)
                  do iz=-1,1
                     add(3)=iz*size(3)
                     dist_tot=0                     
                     do idir=1,3
                        distc(idir)=posion_i(idir,i_ion)-
     $                       posion_j(idir,j_ion)+add(idir)
                        dist_tot=dist_tot+distc(idir)*distc(idir)
                        write(100,*) distc(idir),posion_i(idir,i_ion),
     $                       posion_j(idir,j_ion),i_ion,j_ion
                     enddo

                     dist_tot=sqrt(dist_tot)
c                     write(100,*) "DIST",i_ion,j_ion,dist_tot
c                     call flush(100,iostat)
                     
                     if(dist_tot.gt.0.and.dist_tot.lt.rcut) then
                        sqrtwpi=sqrt(3.14159*w)  
                        iigauss=nint(sqrt(w/sizemesh)*40)+1
                        
                        sqrtwpi=sqrt(3.14159*w)
                        mcenter=nint(dist_tot/sizemesh)
                        mstart=mcenter-iigauss
                        mend=mcenter+iigauss
                        if(mstart.lt.1) mstart=1
                        if(mend.gt.1000) mend=1000

c   make the pdf by gaussian smearing of the distance between the ions                         
                        do m=mstart,mend    
                           r=m*sizemesh-dist_tot
                           apdf(m)=apdf(m)+exp(-r*r/w)/sqrtwpi
                        enddo
                     endif
                  enddo
               enddo
            enddo
                  
         enddo
      enddo

c     normalize the pdf by the factor of the sphere surface area
      pi=acos(-1.0)
      do m=1,numpdf
         r=m*sizemesh
         apdf(m)=apdf(m)/(4*pi*r*r)
      enddo
      do m=1,numpdf,10

        write(1001,*) m*sizemesh,apdf(m)
      enddo
      end
      


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine ao_analyze(num_spec_a,num_spec_x,
     $     num_ion_a,num_ion_x,xa_name,xx_name,
     $     val_r_a,val_N_a,
     $     posion_a,posion_b,posion_x,size)

      implicit real*8(a-h,o-z)
      parameter(nspec_max=10,num_ion_max=200,num_max_dist=100)
      character*2 xa_name(nspec_max)
      dimension num_ion_a(nspec_max),num_ion_x(nspec_max),
     $     posion_a(3,num_ion_max,nspec_max),
     $     posion_x(3,num_ion_max,nspec_max)
      dimension size(3),dist_AO(num_max_dist),mxspec_AO(num_max_dist),
     $     mion_AO(num_max_dist),center_x(3),pos_x(3),disp(3)
      dimension  val_r_a(nspec_max,nspec_max),
     $     val_N_a(nspec_max,nspec_max)
      dimension add(3),vctr_sum(4),distc(3)

      cutoff_AO=5.0
      ncut_AO=12

      do iaspec=1,num_spec_a
         open(1002,file=xa_name(iaspec)//'-O_complex_data')
         do i_ion=1,num_ion_a(iaspec)

c  assemble all the A-X distances less than cutoff_AO (in Ang) for  a given A-site cation
            icnt=0
            do jxspec=1,num_spec_x
               do j_ion=1,num_ion_x(jxspec)

c  this takes care of the periodic images of the ions  in the 26 cells 
c  bordering the unit cell in the input file
                  do ix=-1,1
                     add(1)=ix*size(1)
                     do iy=-1,1
                        add(2)=iy*size(2)
                        do iz=-1,1
                           add(3)=iz*size(3)
                           dist_tot=0
                           
                           do idir=1,3
                              distc(idir)=posion_a(idir,i_ion,iaspec)-
     $                             posion_x(idir,j_ion,jxspec)+add(idir)
                              dist_tot=dist_tot+distc(idir)*distc(idir)
                           enddo
                           dist_tot=sqrt(dist_tot)
                           if(dist_tot.lt.cutoff_AO) then
                              icnt=icnt+1
                              dist_AO(icnt)=dist_tot
                              mxspec_AO(icnt)=jxspec
                              mion_AO(icnt)=j_ion
                           endif
                        enddo
                     enddo
                  enddo
               enddo
            enddo

c     sort the distances from shortest to longest
            call sort(icnt,dist_AO,mxspec_AO,mion_AO)           
c  print the shortest ncut_AO+4 A-X distances and calculate the bond valence of the A-cation
            val_tot=0
            write(1002,*) " "
            write(1002,*) xa_name(iaspec)//' ion ',i_ion
c            call flush(1002,iostat)
            do idir=1,3
               vctr_sum(idir)=0
            enddo
            do m=1,ncut_AO+4
               ixspec=mxspec_AO(m)
               val=(dist_AO(m)
     $              /val_r_a(iaspec,ixspec))**
     $              (-val_N_a(iaspec,ixspec))
               write(1002,*) dist_AO(m),val,mxspec_AO(m),m
               val_tot=val_tot+val
            enddo

            write(1002,*) "VAL",val_tot,i_ion

c      calculate the center of the O12 cage and then use it to calculate the off-center displacement of the A-site cation
            do idir=1,3
               center_x(idir)=0
            enddo
            do m=1,ncut_AO
               jxspec=mxspec_AO(m)
               j_ion=mion_AO(m)
               dist_tot=0
               do idir=1,3
                  pos_x(idir)=posion_x(idir,j_ion,jxspec)
                  distc(idir)=posion_a(idir,i_ion,iaspec)
     $                 -posion_x(idir,j_ion,jxspec)
                  if(distc(idir).gt.size(idir)/2) pos_x(idir)=
     $                 pos_x(idir)+size(idir)
                  if(distc(idir).lt.-size(idir)/2) pos_x(idir)=
     $                 pos_x(idir)-size(idir)
                  center_x(idir)=center_x(idir)+pos_x(idir)
                  distc(idir)=posion_a(idir,i_ion,iaspec)-pos_x(idir)
                  dist_tot=dist_tot+distc(idir)*distc(idir)
               enddo
               dist_tot=sqrt(dist_tot)
               val=((dist_tot)
     $              /val_r_a(iaspec,ixspec))**
     $              (-val_N_a(iaspec,ixspec))
               
c     calculate bond vector and add it to the bond vector sum
               do idir=1,3
                  vctr_sum(idir)=vctr_sum(idir)+
     $                 distc(idir)/dist_tot*val
               enddo
               
            enddo
            
            disp_tot=0
            do idir=1,3
               center_x(idir)=center_x(idir)/ncut_AO
               disp(idir)=posion_a(idir,i_ion,iaspec)-center_x(idir)
               disp_tot=disp_tot+disp(idir)*disp(idir)
            enddo
            disp_tot=sqrt(disp_tot)              

            write(1002,200) disp(1),disp(2),disp(3),disp_tot,i_ion
            vctr_sum(4)=vctr_sum(1)**2+vctr_sum(2)**2+vctr_sum(3)**2
            vctr_sum(4)=sqrt(vctr_sum(4))
            write(1002,202) vctr_sum(1),vctr_sum(2),vctr_sum(3),
     $           vctr_sum(4),i_ion
            write(1002,203) vctr_sum(1)/val_tot*3,vctr_sum(2)/val_tot*3,
     $           vctr_sum(3)/val_tot*3,vctr_sum(4)/val_tot,i_ion
         enddo
         close(1002)
      enddo
 200  format(1x,'Disp',4f15.6,i1)
 202  format(1x,'Vctr',4f15.6,i1)
 203  format(1x,'Asym',4f15.6,i1)
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine bo_analyze(num_spec_b,num_spec_x,
     $     num_ion_b,num_ion_x,xb_name,xx_name,
     $     val_r_b,val_N_b,
     $     posion_a,posion_b,posion_x,size)

      implicit real*8(a-h,o-z)
      parameter(nspec_max=10,num_ion_max=200,num_max_dist=100,max=6)
      character*2 xb_name(nspec_max)
      dimension num_ion_b(nspec_max),num_ion_x(nspec_max),
     $     posion_b(3,num_ion_max,nspec_max),
     $     posion_x(3,num_ion_max,nspec_max)
      dimension size(3),dist_BO(num_max_dist),mxspec_BO(num_max_dist),
     $     mion_BO(num_max_dist),center_x(3),pos_x(3,max),disp(3)
      dimension  val_r_b(nspec_max,nspec_max),
     $     val_N_b(nspec_max,nspec_max)
      dimension add(3),distc(3),ind_x(3),ind_y(3),ind_z(3)
      dimension center_x_cov(3),disp_cov(3),tilt(3),glazer(3),
     $     deformation(3),pos_b(3)
      dimension vctr_sum(4)

      cutoff_BO=5.0
      ncut_BO=6
      pi=acos(-1.0)

      do ibspec=1,num_spec_b
         open(1002,file=xb_name(ibspec)//'-O_complex_data')
         do i_ion=1,num_ion_b(ibspec)
c  assemble all the B-X distances less than cutoff_BO (in Ang) for  a given B-site cation
            icnt=0
            do jxspec=1,num_spec_x
               do j_ion=1,num_ion_x(jxspec)

c  this takes care of the periodic images of the ions  in the 26 cells 
c  bordering the unit cell in the input file
                  do ix=-1,1
                     add(1)=ix*size(1)
                     do iy=-1,1
                        add(2)=iy*size(2)
                        do iz=-1,1
                           add(3)=iz*size(3)
                           dist_tot=0
                           
                           do idir=1,3
                              distc(idir)=posion_b(idir,i_ion,ibspec)-
     $                             posion_x(idir,j_ion,jxspec)
     $                             +add(idir)
                              dist_tot=dist_tot+distc(idir)*distc(idir)
                           enddo
                           dist_tot=sqrt(dist_tot)
                           if(dist_tot.lt.cutoff_BO) then
                              icnt=icnt+1
                              dist_BO(icnt)=dist_tot
                              mxspec_BO(icnt)=jxspec
                              mion_BO(icnt)=j_ion
                           endif
                        enddo
                     enddo
                  enddo
               enddo
            enddo

c     sort the distances from shortest to longest
            call sort(icnt,dist_BO,mxspec_BO,mion_BO)           
c  print the shortest ncut_BO+4 B-X distances and calculate the bond valence of the B-cation
            val_tot=0
            write(1002,*) " "
            write(1002,*) xb_name(ibspec)//' ion ',i_ion
c            call flush(1002,iostat)
            do m=1,ncut_BO+4
               ixspec=mxspec_BO(m)
               val=(dist_BO(m)
     $              /val_r_b(ibspec,ixspec))**
     $              (-val_N_b(ibspec,ixspec))
               jj=mion_BO(m)
               jx=mxspec_BO(m)
               write(1002,*) dist_BO(m),val,mxspec_BO(m),m,mion_BO(m),
     $             posion_x(1,jj,jx),posion_x(2,jj,jx),posion_x(3,jj,jx)
               write(1002,*) 
               val_tot=val_tot+val
            enddo

            write(1002,*) "VAL",val_tot,i_ion
            
c      calculate the center of the O6 cage and then use it to calculate the off-center displacement of the B-site cation
c      also calculate the covalent displacement.  Here the center of the O6 cage along a cartesian direction is only computed using the two O atoms located along that cartesian direction, see PRB 69, 144118  (2004).

            do idir=1,3
               center_x(idir)=0
               center_x_cov(idir)=0
               vctr_sum(idir)=0
            enddo
            ixcount=0
            iycount=0
            izcount=0
            do m=1,ncut_BO
               jxspec=mxspec_BO(m)
               j_ion=mion_BO(m)
               dist_tot=0
               do idir=1,3
                  pos_x(idir,m)=posion_x(idir,j_ion,jxspec)
                  distc(idir)=posion_b(idir,i_ion,ibspec)
     $                 -posion_x(idir,j_ion,jxspec)
                  if(distc(idir).gt.size(idir)/2)  then
                     pos_x(idir,m)=pos_x(idir,m)+
     $                    size(idir)
                     distc(idir)=distc(idir)-size(idir)
                  endif
                  if(distc(idir).lt.-size(idir)/2) then
                     pos_x(idir,m)=pos_x(idir,m)
     $                 -size(idir)
                     distc(idir)=distc(idir)+size(idir)
                  endif
                  center_x(idir)=center_x(idir)+pos_x(idir,m)
                  dist_tot=dist_tot+distc(idir)*distc(idir)
              enddo
              dist_tot=sqrt(dist_tot)
              val=((dist_tot)
     $              /val_r_b(ibspec,ixspec))**
     $              (-val_N_b(ibspec,ixspec))

c calculate bond vector and add it to the bond vector sum
              do idir=1,3
                 vctr_sum(idir)=vctr_sum(idir)+
     $                distc(idir)/dist_tot*val
              enddo


c   mark the O atoms as belonging to either x, y or z axis of the octahedron
               if(acos(abs(distc(1))/dist_tot).lt.pi/6)  then
                  ixcount=ixcount+1
                  ind_x(ixcount)=m
                  center_x_cov(1)=center_x_cov(1)
     $                 +pos_x(1,m)/2
               endif
               if(acos(abs(distc(2))/dist_tot).lt.pi/6) then
                  iycount=iycount+1
                  ind_y(iycount)=m
                  center_x_cov(2)=center_x_cov(2)+
     $                 pos_x(2,m)/2
               endif
               if(acos(abs(distc(3))/dist_tot).lt.pi/6)   then
                  izcount=izcount+1
                  ind_z(izcount)=m
                  center_x_cov(3)=center_x_cov(3)+
     $                 pos_x(3,m)/2
               endif

            enddo


            disp_tot=0
            disp_cov_tot=0
            do idir=1,3
               center_x(idir)=center_x(idir)/ncut_BO
               disp(idir)=posion_b(idir,i_ion,ibspec)-center_x(idir)
               disp_tot=disp_tot+disp(idir)*disp(idir)

               disp_cov(idir)=posion_b(idir,i_ion,ibspec)
     $              -center_x_cov(idir)
c               write(1002,*) posion_b(idir,i_ion,ibspec),
c     $              center_x_cov(idir),idir,i_ion
               disp_cov_tot=disp_cov_tot+disp_cov(idir)*disp_cov(idir)
            enddo
            disp_tot=sqrt(disp_tot)
            disp_cov_tot=sqrt(disp_cov_tot)
 200        format(1x,'DispTot',4f15.6,i1)
 201        format(1x,'DispCov',4f15.6,i1)
 202        format(1x,'%COV',4f15.6,i1)
 206        format(1x,'Vctr',4f15.6,i1)
 207        format(1x,'Asym',4f15.6,i1)
            write(1002,200) disp(1),disp(2),disp(3),disp_tot,i_ion
            write(1002,201) disp_cov(1),disp_cov(2),disp_cov(3),
     $           disp_cov_tot,i_ion
            write(1002,202) disp_cov(1)/disp(1)*100,
     $           disp_cov(2)/disp(2)*100,
     $           disp_cov(3)/disp(3)*100,disp_cov_tot/disp_tot*100
     $           ,i_ion
            vctr_sum(4)=vctr_sum(1)**2+vctr_sum(2)**2+vctr_sum(3)**2
            vctr_sum(4)=sqrt(vctr_sum(4))
            write(1002,206) vctr_sum(1),vctr_sum(2),vctr_sum(3),
     $           vctr_sum(4),i_ion
            write(1002,207) vctr_sum(1)/val_tot*3,vctr_sum(2)/val_tot*3,
     $           vctr_sum(3)/val_tot*3,vctr_sum(4)/val_tot,i_ion

c     calculate three tilt angles, tilts away from cartesian axis, glazer angle tilts around a cartesian axis, and the deformation angle between the three octahedral axes
            
            do idir=1,3
               pos_b(idir)=posion_b(idir,i_ion,ibspec)
            enddo
            call compute_tilt(ind_x,ind_y,ind_z,pos_x,size,
     $           tilt,glazer,deformation,pos_b)
 203        format(1x,'Tilt',3f15.6,i1)
 204        format(1x,'GlazerAngle',3f15.6,i1)
 205        format(1x,'DeformAngle',3f15.6,i1)
            write(1002,203) tilt(1),tilt(2),tilt(3),i_ion
            glazerX=(tilt(2)+tilt(3))/2
            glazerY=(tilt(1)+tilt(3))/2
            glazerZ=(tilt(1)+tilt(2))/2
            write(1002,204) glazer(1),glazer(2),glazer(3),i_ion
            write(1002,205) deformation(1),deformation(2),
     $           deformation(3),i_ion
         enddo
      enddo

c      call flush(1002,iostat)
c 1002 format('Glazer',3f14.5) 
c 1003 format('glazer',3f14.5) 
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine oxygen_analyze(num_spec_a,num_spec_b,num_spec_x,
     $     num_ion_a,num_ion_b,num_ion_x,xa_name,xb_name,xx_name,
     $     val_r_a,val_r_b,val_N_a,val_N_b,
     $     posion_a,posion_b,posion_x,size)

      implicit real*8(a-h,o-z)
      parameter(nspec_max=10,num_ion_max=200,num_max_dist=100,max=6)
      character*1 xx_name(nspec_max)
      character*2 xb_name(nspec_max),xa_name(nspec_max)
      dimension num_ion_b(nspec_max),num_ion_x(nspec_max),
     $     num_ion_a(nspec_max),
     $     posion_a(3,num_ion_max,nspec_max),
     $     posion_b(3,num_ion_max,nspec_max),
     $     posion_x(3,num_ion_max,nspec_max)
      dimension size(3),dist_BO(num_max_dist),mxspec_BO(num_max_dist),
     $     mion_BO(num_max_dist),center_x(3),pos_x(3,max),disp(3)
      dimension val_r_b(nspec_max,nspec_max),
     $     val_r_a(nspec_max,nspec_max),
     $     val_N_a(nspec_max,nspec_max),val_N_b(nspec_max,nspec_max)
      dimension add(3),distc(3),ind_x(3),ind_y(3),ind_z(3)
      dimension center_x_cov(3),disp_cov(3),tilt(3),glazer(3),
     $     deformation(3),pos_b(3)
      dimension maspec_AO(num_max_dist),mbspec_BO(num_max_dist)
      dimension mion_AO(num_max_dist),dist_AO(num_max_dist)
      dimension vctr_sum(4),dvec(3,num_max_dist),ipoint(num_max_dist)

      cutoff_BO=5.0
      cutoff_AO=5.01
      ncut_BO=2
      ncut_AO=4
      pi=acos(-1.0)
      write(1001,*) "IN HERE",num_ion_a(1)


      do ixspec=1,num_spec_x
         open(1002,file=xx_name(ixspec)//'-AB_complex_data')
         write(1001,*) "IN HERE",num_spec_x,ixspec

         do i_ion=1,num_ion_x(ixspec)
           write(1002,*) " "
           write(1002,*) xx_name(ixspec)//' ion ',i_ion,ixspec
c            call flush(1002,iostat)
            val_tot=0 
            do idir=1,3
               vctr_sum(idir)=0
            enddo
c  assemble all the A-X distances less than cutoff_AO (in Ang) for  a given A-site cation
            icnt=0
            write(1002,*) num_spec_a,num_ion_a(iaspec)
            do jaspec=1,num_spec_a
c               write(1002,*) "GOT HERE OK",num_ion_a(iaspec),
c     $              iaspec
               do j_ion=1,num_ion_a(jaspec)
                  call find_dist(dist_AO,dvec,maspec_AO,mion_AO,
     $                 size,posion_x,posion_a,cutoff_AO,icnt,
     $                 ixspec,i_ion,jaspec,j_ion)
c                  write(1002,*) dist_AO(icnt),icnt
               enddo
            enddo

c     sort the distances from shortest to longest
            do ipt=1,num_max_dist
               ipoint(ipt)=ipt
            enddo
            call sort3(icnt,dist_AO,maspec_AO,mion_AO,ipoint)           
c  print the shortest ncut_AO+4 A-X distances and calculate the bond valence of the A-cation
c            call flush(1002,iostat)
            do m=1,ncut_AO
               iaspec=maspec_AO(m)
               val=(dist_AO(m)
     $              /val_r_a(iaspec,ixspec))**
     $              (-val_N_a(iaspec,ixspec))
               write(1002,*) dist_AO(m),val,maspec_AO(m),m
               val_tot=val_tot+val
               do idir=1,3
                  distc(idir)=dvec(idir,ipoint(m))
                  vctr_sum(idir)=vctr_sum(idir)+
     $                 distc(idir)/dist_AO(m)*val
               enddo

            enddo

c  assemble all the B-X distances less than cutoff_BO (in Ang) for  a given B-site cation
            icnt=0
            do jbspec=1,num_spec_b
               do j_ion=1,num_ion_b(jbspec)

                  call find_dist(dist_BO,dvec,mbspec_BO,mion_BO,
     $                 size,posion_x,posion_b,cutoff_BO,icnt,
     $                 ixspec,i_ion,jbspec,j_ion)
               enddo
            enddo

c     sort the distances from shortest to longest
            do ipt=1,num_max_dist
               ipoint(ipt)=ipt
            enddo
            call sort3(icnt,dist_BO,mbspec_BO,mion_BO,ipoint)           
c  print the shortest ncut_BO+4 B-X distances and calculate the bond valence of the B-cation
            do m=1,ncut_BO
               ibspec=mbspec_BO(m)
               val=(dist_BO(m)
     $              /val_r_b(ibspec,ixspec))**
     $              (-val_N_b(ibspec,ixspec))
                  write(1002,*) dist_BO(m),val,mbspec_BO(m),m
               val_tot=val_tot+val
               do idir=1,3
                  distc(idir)=dvec(idir,ipoint(m))
                  vctr_sum(idir)=vctr_sum(idir)+
     $                 distc(idir)/dist_AO(m)*val
               enddo
            enddo

            write(1002,*) "VAL",val_tot,i_ion
            vctr_sum(4)=vctr_sum(1)**2+vctr_sum(2)**2+vctr_sum(3)**2
            vctr_sum(4)=sqrt(vctr_sum(4))
c            write(1002,202) vctr_sum(1),vctr_sum(2),vctr_sum(3),
c     $           vctr_sum(4),i_ion
c            write(1002,203) vctr_sum(1)/val_tot*3,vctr_sum(2)/val_tot*3,
c     $           vctr_sum(3)/val_tot*3,vctr_sum(4)/val_tot,i_ion

         enddo
      enddo

c 202  format(1x,'Vctr',4f15.6,i)
c 203  format(1x,'Asym',4f15.6,i)
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine compute_tilt(ind_x,ind_y,ind_z,pos_x,size,
     $           tilt,glazer,deformation,pos_b)
      parameter(max=6)
      implicit real*8(a-h,o-z)
      dimension tilt(3),glazer(3), size(3),distc(3)
      dimension ind_x(2),ind_y(2),ind_z(2)
      dimension pos_x(3,max),vec(3,3)
      dimension deformation(3),pos_b(3)

      pi=acos(-1.0)
      do idir=1,3
         glazer(idir)=0
         tilt(idir)=0
      enddo
      m1=ind_x(1)
      m2=ind_x(2)
      dist_tot=0
      idir=1

      dist1=pos_x(idir,m1)-pos_b(idir)
      dist2=pos_x(idir,m2)-pos_b(idir)

      if(dist1.gt.0.and.dist2.gt.0) then
         if(dist1.gt.dist2) pos_x(idir,m1)=
     $        pos_x(idir,m1)-size(idir)
         if(dist2.gt.dist1) pos_x(idir,m2)=
     $        pos_x(idir,m2)-size(idir)
      elseif(dist1.lt.0.and.dist2.lt.0) then
         if(dist1.lt.dist2) pos_x(idir,m1)=
     $        pos_x(idir,m1)+size(idir)
         if(dist2.lt.dist1) pos_x(idir,m2)=
     $        pos_x(idir,m2)+size(idir)
      endif



c     calculate the vector from the left  O atom to the right O atom along the x-axis of the octahedron
      itop=m1
      ibot=m2
      if(pos_x(idir,m1).lt.pos_x(idir,m2)) then
         itop=m2
         ibot=m1
      endif
      do idir=1,3
         distc(idir)=pos_x(idir,itop)-pos_x(idir,ibot)
         dist_tot=dist_tot+distc(idir)*distc(idir)
         vec(idir,1)=distc(idir)
      enddo

c     calculate the tilt of the octahedral axis away from x-axis, and the contribution of this O6 axis to the glazer tilt around the y and z cartesian axes
      dist_tot=sqrt(dist_tot)
      tilt(1)=acos(distc(1)/dist_tot)*180/pi
      if(tilt(1).gt.90) tilt(1)=tilt(1)-180
      glazer(2)=atan(distc(3)/distc(1))/2
      glazer(3)=atan(distc(2)/distc(1))/2


      m1=ind_y(1)
      m2=ind_y(2)
      dist_tot=0
      idir=2
      dist1=pos_x(idir,m1)-pos_b(idir)
      dist2=pos_x(idir,m2)-pos_b(idir)

      if(dist1.gt.0.and.dist2.gt.0) then
         if(dist1.gt.dist2) pos_x(idir,m1)=
     $        pos_x(idir,m1)-size(idir)
         if(dist2.gt.dist1) pos_x(idir,m2)=
     $        pos_x(idir,m2)-size(idir)
      elseif(dist1.lt.0.and.dist2.lt.0) then
         if(dist1.lt.dist2) pos_x(idir,m1)=
     $        pos_x(idir,m1)+size(idir)
         if(dist2.lt.dist1) pos_x(idir,m2)=
     $        pos_x(idir,m2)+size(idir)
      endif


c     calculate the vector from the front  O atom to the back O atom along the y-axis of the octahedron
      itop=m1
      ibot=m2
      if(pos_x(idir,m1).lt.pos_x(idir,m2)) then
         itop=m2
         ibot=m1
      endif
      do idir=1,3
         distc(idir)=pos_x(idir,itop)-pos_x(idir,ibot)
         dist_tot=dist_tot+distc(idir)*distc(idir)
         vec(idir,2)=distc(idir)
      enddo

c     calculate the tilt of the octahedral axis away from y-axis, and the contribution of this O6 axis to the glazer tilt around the x and z cartesian axes
      dist_tot=sqrt(dist_tot)
      tilt(2)=acos(distc(2)/dist_tot)*180/pi
      if(tilt(2).gt.90) tilt(2)=tilt(2)-180
      glazer(1)=atan(distc(3)/distc(2))/2
      glazer(3)=glazer(3)-atan(distc(1)/distc(2))/2
    
c     calculate the vector from the top  O atom to the bottom O atom along the z-axis of the octahedron
      m1=ind_z(1)
      m2=ind_z(2)
      dist_tot=0
      idir=3
      dist1=pos_x(idir,m1)-pos_b(idir)
      dist2=pos_x(idir,m2)-pos_b(idir)

      if(dist1.gt.0.and.dist2.gt.0) then
         if(dist1.gt.dist2) pos_x(idir,m1)=
     $        pos_x(idir,m1)-size(idir)
         if(dist2.gt.dist1) pos_x(idir,m2)=
     $        pos_x(idir,m2)-size(idir)
      elseif(dist1.lt.0.and.dist2.lt.0) then
         if(dist1.lt.dist2) pos_x(idir,m1)=
     $        pos_x(idir,m1)+size(idir)
         if(dist2.lt.dist1) pos_x(idir,m2)=
     $        pos_x(idir,m2)+size(idir)
      endif


      itop=m1
      ibot=m2
      if(pos_x(idir,m1).lt.pos_x(idir,m2)) then
         itop=m2
         ibot=m1
      endif
      do idir=1,3
         distc(idir)=pos_x(idir,itop)-pos_x(idir,ibot)
         dist_tot=dist_tot+distc(idir)*distc(idir)
         vec(idir,3)=distc(idir)
      enddo


c     calculate the tilt of the octahedral axis away from z-axis, and the contribution of this O6 axis to the glazer tilt around the x and y cartesian axes
      dist_tot=sqrt(dist_tot)
      tilt(3)=acos(distc(3)/dist_tot)*180/pi
      if(tilt(3).gt.90) tilt(3)=tilt(3)-180
      glazer(1)=glazer(1)-atan(distc(2)/distc(3))/2
      glazer(2)=glazer(2)-atan(distc(1)/distc(3))/2

      do idir=1,3
         glazer(idir)=glazer(idir)*180/pi
      enddo

c calculate the deformation angles.  A perfect O6 octahedron has 90 degree angles between all three octahedral axes.
      icount=0
      do ivec=1,3
         do jvec=ivec+1,3
            dot=0
            do idir=1,3
               dot=dot+vec(idir,ivec)*vec(idir,jvec)
            enddo
            angle=dot/sqrt(vec(1,ivec)**2+vec(2,ivec)**2
     $           +vec(3,ivec)**2)/sqrt(vec(1,jvec)**2
     $           +vec(2,jvec)**2+vec(3,3)**jvec)
            angle=acos(angle)*180/pi
            if(angle.lt.0) angle=angle+90
            if(angle.gt.0) angle=angle-90
            icount=icount+1
           deformation(icount)=angle
         enddo
      enddo
      
      end
      

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine sort(n,arr,mxspec,mion)
      implicit real*8(a-h,o-z)  
      parameter(num_max_dist=100)
      dimension arr(num_max_dist),mion(num_max_dist),
     $     mxspec(num_max_dist)
      
c  this routine is based on the shell's straight insertion sorting method, chapter 8.1 in Numerical Recipes in Fortran 77: the Art of Scientific Computing

      do 12 j=2,n
         a=arr(j)
         ind_xion=mion(j)      
         ind_xspec=mxspec(j)
         do 11 i=j-1,1,-1
            if(arr(i).le.a) goto 10
            arr(i+1)=arr(i)
            mion(i+1)=mion(i)
            mxspec(i+1)=mxspec(i)
 11      continue
         i=0
 10      arr(i+1)=a
         mion(i+1)=ind_xion
         mxspec(i+1)=ind_xspec
 12   continue
      return
      end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine sort3(n,arr,mxspec,mion,ipoint)
      implicit real*8(a-h,o-z)  
      parameter(num_max_dist=100)
      dimension arr(num_max_dist),mion(num_max_dist),
     $     mxspec(num_max_dist),ipoint(num_max_dist)
      
c  this routine is based on the shell's straight insertion sorting method, chapter 8.1 in Numerical Recipes in Fortran 77: the Art of Scientific Computing

      do 12 j=2,n
         a=arr(j)
         ind_xion=mion(j)      
         ind_xspec=mxspec(j)
         ind_ipt=ipoint(j)
         do 11 i=j-1,1,-1
            if(arr(i).le.a) goto 10
            arr(i+1)=arr(i)
            mion(i+1)=mion(i)
            mxspec(i+1)=mxspec(i)
            ipoint(i+1)=ipoint(i)
 11      continue
         i=0
 10      arr(i+1)=a
         mion(i+1)=ind_xion
         mxspec(i+1)=ind_xspec
         ipoint(i+1)=ind_ipt
 12   continue
      return
      end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine sort2(n,arr,at,mxspec,mion)
      implicit real*8(a-h,o-z)  
      parameter(num_max_dist=100)
      dimension arr(num_max_dist),mion(num_max_dist),
     $     mxspec(num_max_dist)
      dimension at(num_max_dist,3)
      
c  this routine is based on the shell's straight insertion sorting method, chapter 8.1 in Numerical Recipes in Fortran 77: the Art of Scientific Computing

      do 12 j=2,n
         a=arr(j)
         ind_xion=mion(j)      
         ind_xspec=mxspec(j)
         at1=at(j,1)
c         at2=at(j,2)
c         at3=at(j,3)
         do 11 i=j-1,1,-1
            if(arr(i).le.a) goto 10
            arr(i+1)=arr(i)
            at(i+1,i)=at(i,1)
c            at(i+1,i)=at(i,2)
c            at(i+1,i)=at(i,3)
            mion(i+1)=mion(i)
            mxspec(i+1)=mxspec(i)
 11      continue
         i=0
 10      arr(i+1)=a
c         at(i+1,1)=at1
c         at(i+1,2)=at2
c         at(i+1,3)=at3
         mion(i+1)=ind_xion
         mxspec(i+1)=ind_xspec
 12   continue
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine find_dist(dist_BO,dvec,mbspec_BO,mion_BO,
     $     size,posion_x,posion_b,cutoff_BO,icnt,
     $     ixspec,i_ion,jbspec,j_ion)
      implicit real*8(a-h,o-z)
      parameter(nspec_max=10,num_ion_max=200,num_max_dist=100,max=6)
      dimension posion_b(3,num_ion_max,nspec_max),
     $     posion_x(3,num_ion_max,nspec_max)
      dimension size(3),dist_BO(num_max_dist),mbspec_BO(num_max_dist),
     $     mion_BO(num_max_dist),distc(3),add(3)
      dimension dvec(3,num_max_dist),ipoint(num_max_dist)

c      write(1002, *) "In find dist",cutoff_BO
      do ix=-1,1
         add(1)=ix*size(1)
         do iy=-1,1
            add(2)=iy*size(2)
            do iz=-1,1
               add(3)=iz*size(3)
               dist_tot=0
               
               do idir=1,3
                  distc(idir)=posion_x(idir,i_ion,ixspec)-
     $                 posion_b(idir,j_ion,jbspec)
     $                 +add(idir)
                  dist_tot=dist_tot+distc(idir)*distc(idir)
               enddo
               dist_tot=sqrt(dist_tot)
c               write(1002, *) dist_tot,icnt,cutoff_BO
               if(dist_tot.lt.cutoff_BO) then
                  icnt=icnt+1
                  do idir=1,3
                     dvec(idir,icnt)=distc(idir)
                  enddo
                  dist_BO(icnt)=dist_tot
                  mbspec_BO(icnt)=jbspec
                  mion_BO(icnt)=j_ion
               endif
            enddo
         enddo
      enddo
      end
      

