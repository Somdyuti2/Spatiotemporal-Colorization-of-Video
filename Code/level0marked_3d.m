clc;
clear all;
close all;
load level1_fore;
load alevel1_fore;
load Iilevel1_fore;
load Iqlevel1_fore;
load g0_fore;
load Iy_fore;
load E_fore;
[m n p]=size(g0);
[mC nC pC]=size(g0);
p0=zeros(mC,nC,pC);

%% Level 0 priorities
num=10;
for k=1:num
for i=2:mC-2
    for j=2:nC-2
        x=i;
        y=j;
        for r=max(1,k-1):min(num,k+1)
        for p=max(1,x-1):min(m,x+1)
            for q=max(1,y-1):min(n,y+1)
                if (a(i,j,k)~=1)
                if (x>=1 && x<=m && y>=1 && y<=n)
                p0(i,j,k)=p0(i,j,k)+a(p,q,r)*exp(-50*(abs(Iy(p,q,r)-Iy(i,j,k))+0.2*E(p,q,r)));
%                   p0(i,j)=p0(i,j)+a(p,q)*exp(-50*g0(i,j));
                end
                end
        end
        end
        end
    end
end
end
% p0(:,:,1)=0;
% p0(:,:,5)=0;
k=find(a~=1);
s=size(k);

tic
%% Level 0 propagation
for t=1:s  
    [maxC,ind] = max(p0(:));
    [m1,n1,p1] = ind2sub(size(p0),ind);
    if (m1==1 || n1==1 || m1==m || n1==n)
    p0(m1,n1,p1)=0;
    a(m1,n1,p1)=1;
    continue;
    end
phi=0;
Ci=0;
Cq=0;

r=[max(1,p1-1):min(num,p1+1)];
p=[max(1,m1-1):min(m,m1+1)];
q=[max(1,n1-1):min(n,n1+1)];
phi=sum(sum(sum(a(p,q,r).*exp(-50*(abs(Iy(p,q,r)-Iy(m1,n1,p1))+0.2*E(p,q,r))))));
Ci=sum(sum(sum(a(p,q,r).*exp(-50*(abs(Iy(p,q,r)-Iy(m1,n1,p1))+0.2*E(p,q,r))).*Ii(p,q,r))));
Cq=sum(sum(sum(a(p,q,r).*exp(-50*(abs(Iy(p,q,r)-Iy(m1,n1,p1))+0.2*E(p,q,r))).*Iq(p,q,r))));
Ii(m1,n1,p1)=Ci/phi;
Iq(m1,n1,p1)=Cq/phi;
a(m1,n1,p1)=1;
p0(m1,n1,p1)=0;


% Priority Updation
for k=max(1,p1-1):min(num,p1+1)
 for i=m1-1:m1+1
    for j=n1-1:n1+1
        x=i; y=j;
        p0(i,j,k)=0;
        r=[max(1,p1-1):min(num,p1+1)];
        p=[max(1,x-1):min(m,x+1)];
        q=[max(1,y-1):min(n,y+1)];
                if (a(i,j,k)~=1)
                 p0(i,j,k)=sum(sum(sum(a(p,q,r).*exp(-50*(abs(Iy(p,q,r)-Iy(i,j,k))+0.2*E(p,q,r))))));
                end
                end
            end
end
% p0(:,:,1)=0;
% p0(:,:,5)=0;
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
save colorized_video colvol