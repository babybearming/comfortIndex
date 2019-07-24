!****************************************
!write by xiongmingming
!in order to contour the comfortableIdex
!****************************************

program main
integer,parameter :: mi=12,h=87,sy=1960,ey=2009,y=ey-sy+1
real,parameter :: aa=0.06,a=0.13,b=0.625
real,parameter :: Pi=3.1415926,PP=180/Pi
real D,M,L,S,F,W,DN !D为计算点经度的度值，M为计算点的分值，DN为积日订正值，S=12，F=0
real N0,x,Deta,TB,Fi,r,x1,w0,p2,I0,T0
integer id,year
real,dimension(366,y) :: TA,Qn
integer,dimension(mi,y) :: mon
real,dimension(mi,y) :: wind,rzss,rzbfl,t,rh !wind:风速;rzbfl:日照百分率;t:平均气温;f:平均相对湿度
real,dimension(mi,y) :: thi,wci,icl,ci !thi:温湿指数;wci:风寒指数;icl:着衣指数;ci:综合舒适指数
real,dimension(mi,y) :: q0,q,RW
real lonlat(2,5)
real lon,lat
real alpha
real ccc,cinian,cimonth
integer id2(5)
integer i,j,mm,k,flag
integer month_day1(13),month_day2(13)
data month_day1/0,31,60,91,121,152,182,213,244,274,305,335,366/
data month_day2/0,31,59,90,120,151,181,212,243,273,304,334,365/


open(1,file='F:\work\comfortIndex\Data\Hohhot.txt')
open(2,file='F:\work\comfortIndex\Data\latlon.txt')
   
do i=1,y
   do j=1,mi
   read(1,*)id,year,mon(j,i),wind(j,i),rzss(j,i),rzbfl(j,i),t(j,i),rh(j,i)
   end do
end do


do i=1,5
   read(2,*)id2(i),lonlat(1,i),lonlat(2,i)
   if(id2(i)==id) then 
   print *,id,lonlat(1,i),lonlat(2,i)
   lon=floor(lonlat(2,i)/100)+(lonlat(2,i)-floor(lonlat(2,i)/100)*100)/60
   lat=floor(lonlat(1,i)/100)+(lonlat(1,i)-floor(lonlat(1,i)/100)*100)/60
   end if
end do

print *,'conting the index'
print *,id,lon,lat
!单位换算
wind=wind*0.1  !m/s
rzss=rzss*0.1 !h
rzbfl=rzbfl*0.01 !%
t=t*0.1 !°C
rh=rh*0.01 !%

open(3,file='thi.txt')
open(4,file='thiidex.txt')
!计算温湿指数THI=(1.8t+32)-0.55(1-f)*(1.8t-26)
do i=1,y
   do j=1,mi
   thi(j,i)=(1.8*t(j,i)+32)-0.55*(1-rh(j,i))*(1.8*t(j,i)-26)
   write(3,*)i+1959,j,thi(j,i)
   if(thi(j,i)<40) then
   thi(j,i)=1
   else if(thi(j,i)<45) then
   thi(j,i)=3
   else if(thi(j,i)<55) then
   thi(j,i)=5
   else if(thi(j,i)<60) then
   thi(j,i)=7
   else if(thi(j,i)<65) then
   thi(j,i)=9
   else if(thi(j,i)<70) then
   thi(j,i)=7
   else if(thi(j,i)<75) then
   thi(j,i)=5
   else if(thi(j,i)<80) then
   thi(j,i)=3
   else
   thi(j,i)=1
   end if
   write(4,*)i+1959,j,thi(j,i)
   end do
end do


open(5,file='wci.txt')
open(6,file='wciidex.txt')
!计算风寒指数
do i=1,y
   do j=1,mi
   wci(j,i)=(33-t(j,i))*(9.0+10.9*sqrt(wind(j,i))-wind(j,i))
       write(5,*)i+1959,j,wci(j,i)  
   if(wci(j,i)>=1000) then
   wci(j,i)=1
   else if(wci(j,i)>=800) then
   wci(j,i)=3 
   else if(wci(j,i)>=600) then
   wci(j,i)=5
   else if(wci(j,i)>=300) then
   wci(j,i)=7
   else if(wci(j,i)>=200)  then
   wci(j,i)=9
   else if(wci(j,i)>=50) then
   wci(j,i)=7
   else if(wci(j,i)>=-80) then
   wci(j,i)=5
   else if(wci(j,i)>=-160) then
   wci(j,i)=3
   else
   wci(j,i)=1
   end if
   write(6,*)i+1959,j,wci(j,i)
   end do
end do 

open(7,file='icl.txt')
open(8,file='iclindex.txt')
!***********************************************************************************************
!计算着衣指数ICL
!***********************************************************************************************
!===============step1计算日天文辐射值========================================
D=floor(lon)
M=(lon-D)*60.0
print *,D,M
Fi=lat/PP
r=34.0/60.0
r=r/PP
S=12
F=0
W=S+F/60
L=-(D+M/60.0)/15.0
DN=(W-L)/24.0
I0=0.001367  !太阳常数
T0=24*60*60

mm=0
do j=SY,EY
   N0=79.6764+0.2422*(j-1985)-INT(0.25*(j-1985))
   mm=mm+1
   do i=0,365
   x=2*Pi*57.3*(i+DN-N0)/365.2422
   x=x/PP
   Deta=0.3723+23.2567*sin(x)+0.1149*sin(2.0*x)-0.1712*sin(3.0*x)-0.7580*cos(x)+0.3656*cos(2.0*x)+0.0201*cos(3.0*x)
   Deta=Deta/PP   !太阳赤纬
   x1=sin(45.0/pp+(Fi-Deta+r)/2)*sin(45.0/PP-(Fi-Deta-r)/2)/(cos(Fi)*cos(Deta))
   TB=2*asin(sqrt(x1))*PP
   !TA=2*TB
   w0=acos(-tan(Fi)*tan(Deta))  !日落时角
   p2=1.000423+0.032359*sin(x)+0.000086*sin(2.0*x)-0.008349*cos(x)+0.000115*cos(2.0*x)  !日地相对距离
   
   TA(i+1,mm)=TB/15.0*2  !可照时数
   Qn(i+1,mm)=T0*I0*(w0*sin(Fi)*sin(Deta)+cos(Fi)*cos(Deta)*sin(w0))/(Pi*p2)
   end do 
end do 

!======================step2计算月太阳辐射总量=====================================
mm=0
flag=0
do k=sy,ey
   mm=mm+1
   do j=1,12
   Q0(j,mm)=0.0
    Q(j,mm)=0.0
	if(mod(k,400)==0 .or. (mod(k,4)==0 .and. mod(k,100)/=0)) then
    do i=month_day1(j)+1,month_day1(j+1)
    Q0(j,mm)=Q0(j,mm)+Qn(i,mm)
	flag=1
	end do
  else
  do i=month_day2(j)+1,month_day2(j+1)
	 Q0(j,mm)=Q0(j,mm)+Qn(i,mm)
	 flag=0
  end do
    end if
	Q(j,mm)=Q0(j,mm)*(a+b*rzbfl(j,mm))
	if(flag==1) then
	RW(j,mm)=Q(j,mm)/((month_day1(j+1)-month_day1(j))*24*60*60)!将Q转换成icl中的R
	else
    RW(j,mm)=Q(j,mm)/((month_day2(j+1)-month_day2(j))*24*60*60)
	end if
   end do
end do

!=================step3计算ICL==========================================
ccc=23.0+26.0/60.0
do i=1,y
   do j=1,mi
   if(j>=6 .and. j<=8) then
   alpha=90-lat+ccc
   else 
   if(j==12 .or. j==1 .or. j==2) then
   alpha=90-lat-ccc
   else
   alpha=90-lat
   end if
   end if
   alpha=alpha/pp
   icl(j,i)=(33-t(j,i))/(0.155*h)-(h+aa*rw(j,i)*cos(alpha))/((0.62+19.0*sqrt(wind(j,i)))*h)
   write(7,*)i+1959,j,icl(j,i)
   if(icl(j,i)>=2.5) then
   icl(j,i)=1
   else if(icl(j,i)>=1.8) then
    icl(j,i)=3
   else if(icl(j,i)>=1.5) then
   icl(j,i)=5
   else if(icl(j,i)>=1.3) then
   icl(j,i)=7
   else if(icl(j,i)>=0.7) then
   icl(j,i)=9
   else if(icl(j,i)>=0.5) then
   icl(j,i)=7
   else if(icl(j,i)>=0.3) then
   icl(j,i)=5
   else if(icl(j,i)>=0.1) then
   icl(j,i)=3
   else
   icl(j,i)=1
   end if
   write(8,*)i+1959,j,icl(j,i)
   end do
end do

open(9,file='ci.txt')
open(10,file='cinian.txt')
open(11,file='cimonth.txt')
!=================计算综合舒适度指数ci==========================================
do i=1,y
   cinian=0.0
   do j=1,mi
   ci(j,i)=0.6*thi(j,i)+0.3*wci(j,i)+0.1*icl(j,i)
!   write(9,*)i+1959,j,ci(j,i)
   cinian=cinian+ci(j,i)
   end do
   write(10,'(i,f7.2)')i+1959,cinian
end do

do j=1,mi
   do i=1,y
   write(9,'(f6.2)',advance='no')ci(j,i)
   end do
   write(9,*)
end do

do j=1,mi
   cimonth=0.0
   do i=1,y
   cimonth=cimonth+ci(j,i)
   end do
   cimonth=cimonth/y
   write(11,'(i,f7.2)')j,cimonth
end do

end 
