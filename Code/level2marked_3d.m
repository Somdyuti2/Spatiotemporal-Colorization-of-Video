clc;
clear all;
close all;
path = '../Foreman_gray_frames';
I=imread(strcat(path,'/foreman1_marked1.bmp'));
% I=imresize(I,0.5);
I1=rgb2ntsc(I);
Ii(:,:,1)=I1(:,:,2);
Iq(:,:,1)=I1(:,:,3);
I=imread(strcat(path,'/foreman1.bmp'));
% I=imresize(I,0.5);
vol(:,:,1)=im2double(I);
% I=rgb2ntsc(I);
Iy(:,:,1)=im2double(I);
num=10;
for k=2:num-1
    I=imread(strcat(path,'/foreman',num2str(k),'.bmp'));
%     I=imresize(I,0.5);
    vol(:,:,k)=im2double(I);
%     I=rgb2ntsc(I);
    Iy(:,:,k)=im2double(I);
    Ii(:,:,k)=0;
    Iq(:,:,k)=0;
end

I=imread(strcat(path,'/foreman10_marked1.bmp'));
% I=imresize(I,0.5);
I1=rgb2ntsc(I);
Ii(:,:,num)=I1(:,:,2);
Iq(:,:,num)=I1(:,:,3);
I=imread('../Foreman_gray_frames/foreman10.bmp');
% I=imresize(I,0.5);
vol(:,:,num)=im2double(I);
% I=rgb2ntsc(I);
Iy(:,:,num)=im2double(I);
[m n p]=size(vol);
a=zeros(m,n,p);

for i=1:num
    a(:,:,i)=Ii(:,:,i)~=0;
end


load g2_fore;
load E_fore;

[mA nA pA]=size(g2);
p2=zeros(mA,nA,pA);
flag2=zeros(mA,nA,pA);

%% Level 2 priorities
for k=1:num
for i=2:mA-2
    for j=2:nA-2
        x=4*i; y=4*j;
        for r=max(1,k-1):min(num,k+1)
        for p=x-3:x+3
            for q=y-3:y+3
                if (x>=1 && x<=m && y>=1 && y<=n)
                p2(i,j,k)=p2(i,j,k)+a(p,q,r)*exp(-50*g2(i,j,k));
            end
        end
        end
    end
end
end
end
% p2(:,:,1)=0;
% p2(:,:,5)=0;

[maxA,ind] = max(p2(:));
[m1,n1,p1] = ind2sub(size(p2),ind);
flag2(m1,n1,p1)=1;
tic

while(maxA>0.1)
x=4*m1; y=4*n1;

for k=max(1,p1-1):min(num,p1+1)
for i=max(1,x-2):min(x+2,m)
    for j=max(1,y-2):min(y+2,n)
        if a(i,j,k)~=1            
phi=0;
Ci=0;
Cq=0;

r=[max(1,p1-1):min(num,p1+1)];
p=[max(1,x-3):min(m,x+3)];
q=[max(1,y-3):min(n,y+3)];

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

for k=max(1,p1-1):min(num,p1+1)
for i=max(1,m1-11):min(m1+11,mA)
    for j=max(1,n1-11):min(n1+11,nA)
        x=4*i; y=4*j;
        p2(i,j,k)=0;
        
        r=[max(1,k-1):min(num,k+1)];
        p=[max(1,x-3):min(m,x+3)];
        q=[max(1,y-3):min(n,y+3)];
        
        if (flag2(i,j,k)~=1)
                p2(i,j,k)=sum(sum(sum(a(p,q,r)*exp(-50*g2(i,j,k)))));
        end
    end
end
end
% p2(:,:,1)=0;
% p2(:,:,5)=0;
[maxA,ind] = max(p2(:));
[m1,n1,p1] = ind2sub(size(p2),ind);
flag2(m1,n1,p1)=1;
end
toc

close all;
for k=1:num
    Icol(:,:,1)=Iy(:,:,k);
    Icol(:,:,2)=Ii(:,:,k);
    Icol(:,:,3)=Iq(:,:,k);
    Icol=ntsc2rgb(Icol);
    colvol(:,:,:,k)=Icol;
end

implay(colvol)
save level2_fore colvol
save Iilevel2_fore Ii
save Iqlevel2_fore Iq
save alevel2_fore a
save Iy_fore Iy