C **********************************************************
C     coherence_analysis_toolik.f        DATE: AUGUST 2017
C **********************************************************
C
c intended to analyze toolik 4-channel data and
c write out three output files suitable for plotting with
c gnuplot.
c 
c currently hard-coded to take eight consecutive 32768-pt FFTs 
c computes two cross-spectra,  center-south and center-west
c converts the phase diffs to azimuth and elevation using 
c cable length delays and antenna positions
c
c this version is run by a master file, called filename
c format is multiple pairs of lines like this:
c 20170827-0956-1000-TLK-INT.dat
c 20170827-0956-1000-TLK-INT.dat 30 75 1250 2500
c ...
c
c means, start writing data at 0956:30, end at 0957:15
c write data only for 1250-2500 kHz
c
c file filename and all data files must be in current directory
c
c output files are in gnuplot format
c time freq data
c
c gnuplot commands to plot:
c set palette rgbformulae 22,13,-31   (for rainbow)
c set cbrange [-150:100]              (for power)
c set cbrange [-180:180]              (for azimuth)
c set cbrange [45:90]              (for elevation)
c plot [30:75][1250:2500] "20170827-0956-1000-TLK-INT.dat-pow.gnuplot" using 1:2:3 with image
c (for power; similar for az or el)
c
c
        character*30 file1
        character*42 outfile1, outfile2, outfile3
	real XMAG(32768),XR(327682),XI(32768)
        real gxx(32768), gyy(32768), gzz(32768)
        real coh12(32768),coh13(32768),phase12(32768),phase13(32768)
        real phase1(32768), phase2(32768), power(32768)
        real data(8,32768)
        integer*1 year,month,day,hour,minute,second,idiff,izero
        integer mincoh(32768)
        integer itime(7),endhour,endminute
        integer*1 byte1, byte2, byte3, byte4, byte5, byte6
        complex cdata(3,32768),gxy(32768),gxz(32768)
c needed for two-byte binary read:
      character c2(2)
      integer*2 i2
      INTEGER fgetc,ibb,status
c needed for fseek
      integer :: stati
c needed for ftell
      integer*8 spencer, nbytes, isample
	real wfn(32768)

c needed for two-byte binary read:
      equivalence (i2,c2)


      izero=0
      zero=0.0
      PI = 3.141592654
      cc=300000000.0
      isample=10000000
      sample=float(isample)
c
c this special version only:
c      write(6,*)"enter selected bin (1-4096)"
c      read(5,*)iselect
c      fselect=float(iselect-1)*(5000.0/4096.0)
c      write(6,*)"you selected ",fselect," kHz"
c      write(6,*)"your output will be in file fort.28"

c  cable length corrections used to compute DOAs
c  "phase12" is the arctangent of GXY, correction1 is subtracted from it
c  "phase13" is the arctangent of GXZ, correction2 is subtracted from it
C         GXY(k)=GXY(k)+(cdata(1,k)*conjg(cdata(2,k)))/bin
C         GXZ(k)=GXZ(k)+(cdata(1,k)*conjg(cdata(3,k)))/bin
C
c  I *think* correction1, 2 should be the extra delay due to additional
c  cabling to X1, Y1 at sondrestrom. hence it is subtracted from the 
c  measured delay/phase between those signals and reference
c  need to do equivalent at toolik: i.e., correction1 should be the 
c  "additional" delay to center beyond that to south (might be a negative
c  number if cable to south is longer) 
c  correction2 should be similar for west vs south.
c  south appears to have the shortest delay of these three, meaning
c  no negative numbers. The delays are very small---effort was made
c  to match the cable lengths. 
c  (If near were ever to be used, both these statements would be false)
c      correction1=25.0e-9
c      correction2=27.5e-9
        correction1=0.0
        correction2=0.0

c  antenna positions x=positve east, y=positive north
c  d1 points from ref to X1
c  d2 points from ref to Y1
c  at toolik, south is ref (first channel read), 
c             center is X1 (second channel read),
c             west is Y1 (third channel read)
c for comparison, at sondrestrom they are read in order 
c  work1  = reference
c  work2  = X1 data
c  work3  = Y1 data
      d1y=17.90
      d1x=44.13
      d2y=66.49
      d2x=23.34

c default number of pts, ensembles per coherence calculation
      npts=32768
      rpts=float(npts)
      nbin=8
      bin=float(nbin)
c calculate time increment in ms
      tinc=(bin*32768.0)/(10000.0)
      idiff=int(tinc)

c open master file, start massive loop over all files
      open(13,file="filename",status="old")
      do mymassiveloop=1,999999

c get the date and start time from the filename on line one
      read(13,821,end=9999)(itime(i), i=1,7)
 821  format(i4,i2,i2,1x,i2,i2,1x,i2,i2)
      itime(1)=itime(1)-2000
      endhour=itime(6)
      endminute=itime(7)
      nbytes=((endhour-itime(4))*60)+(endminute-itime(5))
      nbytes=nbytes*60*4*2*isample
      write(6,*)(itime(i), i=1,7),endhour,endminute,nbytes
c get filename, seconds to start, seconds to end, flow and fhigh 
c from line two
      read(13,*,end=9999)file1,sec1,sec2,flow,fhigh
c  open input file
      open(11,file=file1,status="old")
c create output files and open them
      write(unit=outfile1,fmt="(a28,'-pow.gnuplot')")file1
      write(unit=outfile2,fmt="(a28,'-azi.gnuplot')")file1
      write(unit=outfile3,fmt="(a28,'-ele.gnuplot')")file1
      open(27,file=outfile1)
      open(28,file=outfile2)
      open(29,file=outfile3)
      write(6,*)"opening files :",file1
      write(6,*)"               ",outfile1
      write(6,*)"               ",outfile2
      write(6,*)"               ",outfile3

c  calculate a bunch of stuff
c      istart=int((sec1*sample)/(bin*rpts))
c     iend=int((sec2*sample)/(bin*rpts))
      istart=int(sec1)
      iend=int((sec2+4)/4)
      nlow=int((flow*rpts/sample)+0.5)
      nhigh=int((fhigh*rpts/sample)+0.5)
      write(6,*)"istart/iend/nlow/nhigh= ",istart,iend,nlow,nhigh

c  open output file 
c      OPEN(20,FILE="coheredata",FORM="UNFORMATTED",
c     *     ACCESS="DIRECT",recl=1,ERR=1000)

c  set weighting fct: (0) boxcar, (1) Hanning, (2) 10% Cosine Bell'
      iwgt=0
      call weight(npts,wfn,iwgt)

c loop for spectra starts here
      do jklmno=1,iend
    
c zero arrays for cross-spectra
      do K=1,NPTS
         GXX(k)=0.0
         GYY(K)=0.0
         GZZ(K)=0.0
         GXY(K)=cmplx(0.0,0.0)
         GXZ(K)=cmplx(0.0,0.0)
      END do
c
c read in nbin sets of data for each of three units,
      do kkk=1,nbin
        do I=1,NPTS
          status=fgetc(11,c2(1))
          status=fgetc(11,c2(2))
          ibb =i2
          data(1,i)=(float(ibb))
          status=fgetc(11,c2(1))
          status=fgetc(11,c2(2))
          ibb =i2
          data(2,i)=(float(ibb))
          status=fgetc(11,c2(1))
          status=fgetc(11,c2(2))
          ibb =i2
          dummy=(float(ibb))
          status=fgetc(11,c2(1))
          status=fgetc(11,c2(2))
          ibb =i2
          data(3,i)=(float(ibb))
        enddo
        do iunit=1,3
          do i=1,npts
            XR(I) = data(iunit,i)*wfn(i)
            XI(I) = zero*wfn(i)
          END do
          IN = 0
c take fft
          CALL FFT842(IN,NPTS,XR,XI)
c store real and imaginary parts to arrays
          do K=1, NPTS/2
            CDATA(IUNIT,K)=cmplx(XR(K),XI(K))
          enddo
        enddo
c formulate and accumulate cross-spectral parameters
        do k=1,npts/2
          GXX(k)=GXX(k)+( (abs(cdata(1,k))**2) /bin)
          GYY(k)=GYY(k)+((abs(cdata(2,k))**2)/bin)
          GZZ(k)=GZZ(k)+((abs(cdata(3,k))**2)/bin)
          GXY(k)=GXY(k)+(cdata(1,k)*conjg(cdata(2,k)))/bin
          GXZ(k)=GXZ(k)+(cdata(1,k)*conjg(cdata(3,k)))/bin
        enddo
      END do

c formulate coherences and phases
      do k=1,npts/2
        coh12(k)=(abs(GXY(k)))/((GXX(K)*GYY(K))**0.5)
        coh13(k)=(abs(GXZ(k)))/((GXX(K)*GZZ(K))**0.5)
        phase12(k)=atan2(aimag(GXY(K)),real(GXY(K)))
        phase12(k)=180.*phase12(k)/pi
        phase13(k)=atan2(aimag(GXZ(K)),real(GXZ(K)))
        phase13(k)=180.*phase13(k)/pi
      enddo

c if this is the very first cross-spectrum, encode exact time
c and write to output, and write zero to output 
      if(jklmno.eq.1)then
c        write(20, rec=1)itime(1)
c        write(20, rec=2)itime(2)
c        write(20, rec=3)itime(3)
c        write(20, rec=4)itime(4)
c        write(20, rec=5)itime(5)
c        write(20, rec=6)itime(6)
c        write(20,rec=7)izero
c        kount=8
      else
c in all other cases, write the time increment in milliseconds
c        write(20,rec=kount)idiff
c        kount=kount+1
      endif

c write the coherences and phases to new arrays, correctly ordered
c      DO Kkk=1,npts/2
c        k=npts/2+kkk-1
c        i12=int(coh12(k)*100.)
c        i13=int(coh13(k)*100.)
c        if(i13.lt.i12)then
c          mincoh(KKK)=i13
c        else
c          mincoh(KKK)=i12
c        endif
c        phase1(KKK)=phase12(k)
c        phase2(KKK)=phase13(k)
c        power(KKK)=gxx(k)
c      ENDDO
c      DO Kkk=1,npts/2
c        k=kkk
c        i12=int(coh12(k)*100.)
c        i13=int(coh13(k)*100.)
c        if(i13.lt.i12)then
c          mincoh(KKK+(npts/2))=i13
c        else
c          mincoh(KKK+(npts/2))=i12
c        endif
c        phase1(KKK+(npts/2))=phase12(k)
c        phase2(KKK+(npts/2))=phase13(k)
c        power(KKK+(npts/2))=gxx(k)
c      ENDDO
      DO k=1,npts/2
        i12=int(coh12(k)*100.)
        i13=int(coh13(k)*100.)
        if(i13.lt.i12)then
          mincoh(k)=i13
        else
          mincoh(k)=i12
        endif
        phase1(k)=phase12(k)
        phase2(k)=phase13(k)
        power(k)=gxx(k)
      ENDDO

c now calculate frequencies and angles, compress and write these

      do j=nlow,nhigh
c compute frequency:
        rfreq=float(j-1)*(10000000.0/rpts)
c only do any calculations if coherence is above 0.50.
c otherwise set angle and elevation to out-of-range defaults
        if(mincoh(j).lt.50)then
          angle=-1000.0
          elevation=-95.0
        else
c apply cable length corrections:
  	  phase1(j)=phase1(j)-360.0*rfreq*correction1
  	  phase2(j)=phase2(j)-360.0*rfreq*correction2
c force phases into range -180 to 180
c		  phase1(j)=mod(phase1(j)+180.0,360.0)-180.0
c		  phase2(j)=mod(phase2(j)+180.0,360.0)-180.0
          do k=1,2
            if(phase1(j).gt.180.0)phase1(j)=phase1(j)-360.0
            if(phase1(j).lt.-180.0)phase1(j)=360.0+phase1(j)
          enddo
          do k=1,2
            if(phase2(j).gt.180.0)phase2(j)=phase2(j)-360.0
            if(phase2(j).lt.-180.0)phase2(j)=360.0+phase2(j)
          enddo
c         phase2(j)=phase2(j)+360.0
c        phase1(j)=phase1(j)-360.0
c calculate the x and y components of the k-vector
	  rkx=(d2y*phase1(j) - d1y*phase2(j))/(d1x*d2y - d2x*d1y)
          rky=(d2x*phase1(j) - d1x*phase2(j))/(d2x*d1y - d2y*d1x)
c atan2(y,x) gives angle off of x-axis (which is east)
          angle=atan2(rky,rkx)
c adjust for around north rather than around east, also convert rad-->degr
          angle=90.0 - ((angle/3.14159)*180.)
c subtract 180 in order to get DOA rather than k-direction
          angle=angle-180.0
c force into range -180 to 180
c		  angle = mod(angle+180.0,360.0)-180.0
	  do k=1,2
            if(angle.gt.180.0)angle=angle-360.0
	    if(angle.lt.-180.0)angle=angle+360.0
	  enddo
c calculate elevation
          rkk=360.0*rfreq/cc
          if((rkk**2 - rkx**2 - rky**2).ge.0.0)then
            rkz=-(rkk**2 - rkx**2 - rky**2)**0.5
            theta=180.0*acos(rkz/rkk)/pi
            elevation=theta-90.0
c set to out-of-range default if not in range 0-90
            if(elevation.lt.0.0.or.elevation.gt.90.0)elevation=-97.0
          else
            elevation=999
            angle=999
          endif
        endif

c write the three bytes to file, but only for the middle 90% of the band;
c the bandwidth is 833.333 kHz, but only 750 kHz is good
c        if(j.gt.npts/20.and.j.le.(19*npts)/20)then
c for toolik: write all values
        
c          i1=int(mincoh(j))
c          byte1=int(i1)
c          write(20,rec=kount)byte1
c          kount=kount+1
          
c          i22=int((angle+180.0)/1.5)
c          byte2=int(i22)
c          write(20,rec=kount)byte2
c          kount=kount+1
          
c          i3=int(elevation)
c          byte3=int(i3)
c          write(20,rec=kount)byte3
c          kount=kount+1

          if (power(j).le.0.0) then
          	power(j)=-100.0
          else
		 	power(j)=100.0*ALOG10(power(j))
		  endif
          i4=int(power(j))
c          if(i4.gt.255)i4=255
c          if(i4.lt.0)i4=0
c          byte4=int(i4)
c          write(20,rec=kount)byte4
c          kount=kount+1

c          i5=int((phase1(j)+180.0)/1.5)
c          byte5=int(i5)
c          write(20,rec=kount)byte5
c          kount=kount+1
          
c          i6=int((phase2(j)+180.0)/1.5)
c          byte6=int(i6)
c          write(20,rec=kount)byte6
c          kount=kount+1
c        endif
c          if(j.eq.iselect) then
c            if(jklmno.gt.13732)then
            if(jklmno.gt.istart)then
              tttt=float(jklmno-1)*4.0
c              freq=float(j-1)*(5000.0/4096.0)
              write(27,833)file1,tttt,rfreq/1000.0,i4,
     *                     int(angle+0.5),int(elevation+0.5)
c              write(28,*)file1,tttt,rfreq/1000.0
              write(28,*)file1,tttt,rfreq/1000.0,phase1(j),phase2(j)
c              write(29,*)file1,tttt,rfreq/1000.0
 833          format(a30,1x,f7.1,1x,f6.1,1x,i4,1x,i4,1x,i4)
           endif
c          endif
c          if(jklmno.eq.1)then
c            freq=float(j-1)*(5000.0/4096.0)
c            write(29,*)freq,mincoh(j),power(j),angle,elevation
c          endif

      enddo

c     increment to the start of the next four second interval
c     that is, read the remaining words until 1M samples of each of four
c     channels have been read.
      do i=1,1000000-(nbin*npts)
         do j=1,4
            status=fgetc(11,c2(1))
            status=fgetc(11,c2(2))
         enddo
      enddo

c check if at end of file
      spencer=ftell(11)
      if(float(jklmno/100).eq.float(jklmno)/100.0)then
        write(6,*)tttt
      endif
c      write(6,*)spencer,nbytes
      if(spencer.gt.nbytes)goto 888
c return to do another spectrum
      enddo

c come to here when finished with one data file, close input, output files
 888  continue
      close(11)
      close(27)
      close(28)
      close(29)
      enddo

c  come to here when massive loop is done, close master file
 9999 continue
      CLOSE(13)

 1000 STOP
      END

C-------------------------------------------------------------------
      SUBROUTINE FFT842(IN,N,X,Y)
C----------------------------------------------------------------
C SUBROUTINE:   FFT842
C FAST FOURIER TRANSFORM FOR N=2**M
C COMPLEX INPUT
C
C
C THIS PROGRAM REPLACES THE VECTOR Z=X+iY BY ITS FINITE
C DISCRETE FOURIER TRANSFORM IF IN=0 AND THE INVERSE
C TRANSFORM IF IN=1. IT PERFORMS AS MANY BASE 8 ITERATIONS
C AS POSSIBLE AND THEN FINISHES WITH A BASE 4 OR A BASE 2
C ITERATION IF NEEDED.
C
      DIMENSION X(*),Y(*),L(15)
      COMMON/CON2/ PI2,P7
      EQUIVALENCE (L15,L(1)),(L14,L(2)),(L13,L(3)),(L12,L(4)),
     *     (L11,L(5)),(L10,L(6)),(L9,L(7)),(L8,L(8)),(L7,L(9)),
     *     (L6,L(10)),(L5,L(11)),(L4,L(12)),(L3,L(13)),(L2,L(14)),
     *     (L1,L(15))
C
C IW IS A MACHINE DEPENDENT WRITE DEVICE NUMBER
C
      IW = 6
C
      PI2 = 8.*ATAN(1.)
      P7 = 1./SQRT(2.)
      DO 10 I=1,15
         M = I
         NT = 2**I
         IF(N.EQ.NT)GOTO 20
   10 CONTINUE
      WRITE(IW,9999)
 9999 FORMAT(35H N IS NOT A POWER OF TWO FOR FFT842)
      STOP
   20 N2POW = M
      NTHPO = N
      FN = NTHPO
      IF(IN.EQ.1)GOTO 40
      DO 30 I=1,NTHPO
         Y(I) = -Y(I)
   30 CONTINUE
   40 N8POW = N2POW/3
      IF(N8POW.EQ.0)GOTO 60
C
C RADIX 8 PASSES IF ANY
C
      DO 50 IPASS=1,N8POW
         NXTLT = 2**(N2POW-3*IPASS)
         LENGT = 8*NXTLT
         CALL R8TX(NXTLT,NTHPO,LENGT,X(1),X(NXTLT+1),X(2*NXTLT+1),
     *        X(3*NXTLT+1),X(4*NXTLT+1),X(5*NXTLT+1),X(6*NXTLT+1),
     *        X(7*NXTLT+1),Y(1),Y(NXTLT+1),Y(2*NXTLT+1),Y(3*NXTLT+1),
     *        Y(4*NXTLT+1),Y(5*NXTLT+1),Y(6*NXTLT+1),Y(7*NXTLT+1))
   50 CONTINUE
C
C IS THERE A FOUR FACTOR LEFT
C
   60 IF(N2POW-3*N8POW-1) 90,70,80
C
C GO THROUGH THE BASE 2 ITERATION
C
   70 CALL R2TX(NTHPO,X(1),X(2),Y(1),Y(2))
      GOTO 90
C
C GO THROUGH THE BASE 4 ITERATION
C
   80 CALL R4TX(NTHPO,X(1),X(2),X(3),X(4),Y(1),Y(2),Y(3),Y(4))
C
   90 DO 110 J=1,15
         L(J) = 1
         IF(J-N2POW) 100,100,110
  100    L(J) = 2**(N2POW+1-J)
  110 CONTINUE
      IJ = 1
      DO 130 J1=1,L1
      DO 130 J2=J1,L2,L1
      DO 130 J3=J2,L3,L2
      DO 130 J4=J3,L4,L3
      DO 130 J5=J4,L5,L4
      DO 130 J6=J5,L6,L5
      DO 130 J7=J6,L7,L6
      DO 130 J8=J7,L8,L7
      DO 130 J9=J8,L9,L8
      DO 130 J10=J9,L10,L9
      DO 130 J11=J10,L11,L10
      DO 130 J12=J11,L12,L11
      DO 130 J13=J12,L13,L12
      DO 130 J14=J13,L14,L13
      DO 130 JI=J14,L15,L14
         IF(IJ-JI) 120,130,130
  120    R = X(IJ)
         X(IJ) = X(JI)
         X(JI) = R
         FI = Y(IJ)
         Y(IJ) = Y(JI)
         Y(JI) = FI
  130    IJ =IJ + 1
      IF(IN.EQ.1)GOTO 150
      DO 140 I=1,NTHPO
         Y(I) = -Y(I)
      X(I) = X(I)/FN
      Y(I) = Y(I)/FN
  140 CONTINUE
      GOTO 170
  150 DO 160 I=1,NTHPO
         X(I) = X(I)
         Y(I) = Y(I)
  160 CONTINUE
 170  CONTINUE
      RETURN
      END
C
C------------------------------------------------------------------
C SUBROUTINE:   R2TX
C RADIX 2 ITERATION SUBROUTINE
C------------------------------------------------------------------
C
      SUBROUTINE R2TX(NTHPO,CR0,CR1,CI0,CI1)
      DIMENSION CR0(*),CR1(*),CI0(*),CI1(*)
      DO 10 K=1,NTHPO,2
         R1 = CR0(K) + CR1(K)
         CR1(K) = CR0(K) - CR1(K)
         CR0(K) = R1
         FI1 = CI0(K) + CI1(K)
         CI1(K) = CI0(K) - CI1(K)
         CI0(K) = FI1
   10 CONTINUE
      RETURN
      END
C
C-----------------------------------------------------------------
C SUBROUTINE:   R4TX
C RADIX 4 ITERATION SUBROUTINE
C-----------------------------------------------------------------
C
      SUBROUTINE R4TX(NTHPO,CR0,CR1,CR2,CR3,CI0,CI1,CI2,CI3)
      DIMENSION CR0(*),CR1(*),CR2(*),CR3(*),CI0(*),CI1(*),
     *      CI2(*),CI3(*)
      DO 10 K=1,NTHPO,4
         R1 = CR0(K) + CR2(K)
         R2 = CR0(K) - CR2(K)
         R3 = CR1(K) + CR3(K)
         R4 = CR1(K) - CR3(K)
         FI1 = CI0(K) + CI2(K)
         FI2 = CI0(K) - CI2(K)
         FI3 = CI1(K) + CI3(K)
         FI4 = CI1(K) - CI3(K)
         CR0(K) = R1 + R3
         CI0(K) = FI1 + FI3
         CR1(K) = R1 - R3
         CI1(K) = FI1 - FI3
         CR2(K) = R2 - FI4
         CI2(K) = FI2 + R4
         CR3(K) = R2 + FI4
         CI3(K) = FI2 - R4
   10 CONTINUE
      RETURN
      END
C
C-----------------------------------------------------------------
C SUBROUTINE:   R8TX
C RADIX 8 ITERATION ROUTINE
C-----------------------------------------------------------------
C
      SUBROUTINE R8TX(NXTLT,NTHPO,LENGT,CR0,CR1,CR2,CR3,CR4,
     *    CR5,CR6,CR7,CI0,CI1,CI2,CI3,CI4,CI5,CI6,CI7)
      DIMENSION CR0(*),CR1(*),CR2(*),CR3(*),CR4(*),CR5(*),
     *    CR6(*),CR7(*),CI0(*),CI1(*),CI2(*),CI3(*),
     *    CI4(*),CI5(*),CI6(*),CI7(*)
      COMMON/CON2/ PI2,P7
C
      SCALE = PI2/FLOAT(LENGT)
      DO 30 J=1,NXTLT
         ARG = FLOAT(J-1)*SCALE
         C1 = COS(ARG)
         S1 = SIN(ARG)
         C2 = C1**2 - S1**2
         S2 = C1*S1 + C1*S1
         C3 = C1*C2 - S1*S2
         S3 = C2*S1 + S2*C1
         C4 = C2**2 - S2**2
         S4 = C2*S2 + C2*S2
         C5 = C2*C3 - S2*S3
         S5 = C3*S2 + S3*C2
         C6 = C3**2 - S3**2
         S6 = C3*S3 + C3*S3
         C7 = C3*C4 - S3*S4
         S7 = C4*S3 + S4*C3
         DO 20 K=J,NTHPO,LENGT
            AR0 = CR0(K) + CR4(K)
            AR1 = CR1(K) + CR5(K)
            AR2 = CR2(K) + CR6(K)
            AR3 = CR3(K) + CR7(K)
            AR4 = CR0(K) - CR4(K)
            AR5 = CR1(K) - CR5(K)
            AR6 = CR2(K) - CR6(K)
            AR7 = CR3(K) - CR7(K)
            AI0 = CI0(K) + CI4(K)
            AI1 = CI1(K) + CI5(K)
            AI2 = CI2(K) + CI6(K)
            AI3 = CI3(K) + CI7(K)
            AI4 = CI0(K) - CI4(K)
            AI5 = CI1(K) - CI5(K)
            AI6 = CI2(K) - CI6(K)
            AI7 = CI3(K) - CI7(K)
            BR0 = AR0 + AR2
            BR1 = AR1 + AR3
            BR2 = AR0 - AR2
            BR3 = AR1 - AR3
            BR4 = AR4 - AI6
            BR5 = AR5 - AI7
            BR6 = AR4 + AI6
            BR7 = AR5 + AI7
            BI0 = AI0 + AI2
            BI1 = AI1 + AI3
            BI2 = AI0 - AI2
            BI3 = AI1 - AI3
            BI4 = AI4 + AR6
            BI5 = AI5 + AR7
            BI6 = AI4 - AR6
            BI7 = AI5 - AR7
            CR0(K) = BR0 + BR1
            CI0(K) = BI0 + BI1
            IF(J.LE.1)GOTO 10
            CR1(K) = C4*(BR0-BR1) - S4*(BI0-BI1)
            CI1(K) = C4*(BI0-BI1) + S4*(BR0-BR1)
            CR2(K) = C2*(BR2-BI3) - S2*(BI2+BR3)
            CI2(K) = C2*(BI2+BR3) + S2*(BR2-BI3)
            CR3(K) = C6*(BR2+BI3) - S6*(BI2-BR3)
            CI3(K) = C6*(BI2-BR3) + S6*(BR2+BI3)
            TR = P7*(BR5-BI5)
            TI = P7*(BR5+BI5)
            CR4(K) = C1*(BR4+TR) - S1*(BI4+TI)
            CI4(K) = C1*(BI4+TI) + S1*(BR4+TR)
            CR5(K) = C5*(BR4-TR) - S5*(BI4-TI)
            CI5(K) = C5*(BI4-TI) + S5*(BR4-TR)
            TR = -P7*(BR7+BI7)
            TI = P7*(BR7-BI7)
            CR6(K) = C3*(BR6+TR) - S3*(BI6+TI)
            CI6(K) = C3*(BI6+TI) + S3*(BR6+TR)
            CR7(K) = C7*(BR6-TR) - S7*(BI6-TI)
            CI7(K) = C7*(BI6-TI) + S7*(BR6-TR)
            GOTO 20
   10       CR1(K) = BR0 - BR1
            CI1(K) = BI0 - BI1
            CR2(K) = BR2 - BI3
            CI2(K) = BI2 + BR3
            CR3(K) = BR2 + BI3
            CI3(K) = BI2 - BR3
            TR = P7*(BR5-BI5)
            TI = P7*(BR5+BI5)
            CR4(K) = BR4 + TR
            CI4(K) = BI4 + TI
            CR5(K) = BR4 - TR
            CI5(K) = BI4 - TI
            TR = -P7*(BR7+BI7)
            TI = P7*(BR7-BI7)
            CR6(K) = BR6 + TR
            CI6(K) = BI6 + TI
            CR7(K) = BR6 - TR
            CI7(K) = BI6 - TI
   20    CONTINUE
   30 CONTINUE
      RETURN
      END
c-------------------------------------------------------------
      subroutine sort(n,p1)
C-------------------------------------------------------------
      dimension p1(32768),w1(32768)
C
      npts = n
      nhalf = npts/2
C
C  D.C. IS REPLACED BY AVERAGE OF THE TWO ADAJACENT FREQUENCIES
c     P1(1)=(P1(2)+P1(NPTS))/2.
C
      do I=1,NHALF-1
      II=I+NHALF
      W1(I)=P1(II+1)
      W1(II-1)=P1(I)
      END do
      W1(NPTS-1)=P1(NHALF)
      W1(NPTS)=P1(NHALF+1)
      do I=1,NPTS
      P1(I)=W1(I)
      END do
      RETURN
      END

C ------------------------------------------------------------
C
      subroutine weight(npts,wfn,iwgt)
C
C ------------------------------------------------------------
C
      real wfn(32768)
      PI=3.14159265358979
      FN=FLOAT(npts)
      PIFN=PI/FN
c      write(6,*)
c     *' Enter: (0) boxcar, (1) Hanning, (2) 10% Cosine Bell'
c      read(5,*)iwgt
C
      if(iwgt.eq.0)then
        do i=1,npts
          WFN(I)=1.0
        enddo
        CONST1=1.0
        BW=1.0
      endif
C
      IF(IWGT.EQ.1)then
        do i=1,npts
          wfn(i)=0.5+0.5*COS((2.0*FLOAT(I)+FN)*PIFN)
        enddo
        CONST1=2.6666666667
        BW=1.5
      endif
C
      IF(IWGT.EQ.2)then
        LWGT=npts/10
        NLWGT=npts-LWGT
        PIWT=PI/FLOAT(LWGT)
        do i=1,npts
          wfn(i)=1.0
          IF(I.LT.LWGT)wfn(i)=0.5*(1.0-COS(FLOAT(I)*PIWT))
          IF(I.GT.NLWGT)wfn(i)=0.5*(1.0-COS(FLOAT(npts-I)*PIWT))
        enddo
        BW=1.1
        CONST1=1.0/0.875
      endif
C
C   Recombine using the following???
C
C      TRAIN=FN/16.
C      F(I)=F(I)/FN
C      PWR1=2.0*F(I)*CONJG(F(I))
C      S(I)=PWR1*CONST1*TRAIN/BW
C
      return
      end
