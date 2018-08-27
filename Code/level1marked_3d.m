clc;
clear all;
close all;
load level2_fore;
load Iilevel2_fore;
load Iqlevel2_fore;
load alevel2_fore;
load g1_fore;
load g0_fore;
load Iy_fore;
load E_fore;

[m n p]=size(g0);
[mB nB pB]=size(g1);
p1=zeros(mB,nB, pB);
flag1=zeros(mB,nB, pB);
num=10;
for k=1:num
    for i=2:mB-2
    for j=2:nB-2
        x=2*i; y=2*j;
        for r=max(1,k-1):min(num,k+1)
            for p=x-2:x+2
            for q=y-2:y+2
                if (x>=1 && x<=m && y>=1 && y<=n)   
                p1(i,j,k)=p1(i,j,k)+a(p,q,r)*exp(-50*g1(i,j,k));
            end
        end
        end
    end
    end
    end
end
% p1(:,:,1)=0;
% p1(:,:,5)=0;
[maxB,ind] = max(p1(:));
[mm,nn,t] = ind2sub(size(p1),ind);
flag1(mm,nn,t)=1;

tic
while(maxB>0.1)
x=2*mm; y=2*nn;
for k=max(1,t-1):min(num,t+1)
for i=max(1,x-1):min(x+1,m)
    for j=max(1,y-1):min(y+1,n)
        if (a(i,j,k)~=1)           
phi=0;
Ci=0;
Cq=0;

r=[max(1,t-1):min(num,t+1)];
p=[max(1,x-2):min(m,x+2)];
q=[max(1,y-2):min(n,y+2)];

phi=sum(sum(sum(a(p,q,r).*exp(-50*(abs(Iy(p,q,r)-Iy(i,j,k))+0.5*E(p,q,r))))));
Ci=sum(sum(sum(a(p,q,r).*exp(-50*(abs(Iy(p,q,r)-Iy(i,j,k))+0.5*E(p,q,r))).*Ii(p,q,r))));
Cq=sum(sum(sum(a(p,q,r).*exp(-50*(abs(Iy(p,q,r)-Iy(i,j,k))+0.5*E(p,q,r))).*Iq(p,q,r))));

Ii(i,j,k)=Ci/phi;
Iq(i,j,k)=Cq/phi;
a(i,j,k)=1;
        end
    end
end
end
% Priority updation
for k=max(1,t-1):min(num,t+1)
    for i=max(1,mm-11):min(mB,mm+11)
    for j=max(1,nn-11):min(nB,nn+11)
        x=2*i; y=2*j;
        p1(i,j,k)=0;
        
        r=[max(1,k-1):min(num,k+1)];
        p=[max(1,x-2):min(m,x+2)];
        q=[max(1,y-2):min(n,y+2)];
        if(flag1(i,j,k)~=1)
        p1(i,j,k)=sum(sum(sum(a(p,q,r).*exp(-50*g1(i,j,k)))));
        end
    end
    end
end
% p1(:,:,1)=0;
% p1(:,:,5)=0;
[maxB,ind] = max(p1(:));
[mm,nn,t] = ind2sub(size(p1),ind);
flag1(mm,nn,t)=1;
end
toc;

close all
for k=1:num
    Icol1(:,:,1)=Iy(:,:,k);
    Icol1(:,:,2)=Ii(:,:,k);
    Icol1(:,:,3)=Iq(:,:,k);
    Icol1=ntsc2rgb(Icol1);
   colvol(:,:,:,k)=Icol1;
end

implay(colvol)
save level1_fore colvol
save Iilevel1_fore Ii
save Iqlevel1_fore Iq
save alevel1_fore a