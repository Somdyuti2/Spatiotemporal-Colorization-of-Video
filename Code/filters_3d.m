clc;
clear all;
close all;
% addpath('frames');
tic

% % 9 Tap Inphase
% g1=[0.0084 0.1025 0.3540 -0.0537 -0.8230 -0.0537 0.3540 0.1025 0.0084]; %even symmetry
% g2=[0.0008 0.0176 0.1660 0.6383 1.0000 0.6383 0.1660 0.0176 0.0008]; %even symmetry
% g3=[-0.0034 -0.0582 -0.3662 -0.7039 0 0.7039 0.3662 0.0582 0.0034]; %odd symmetry
% g4=[-0.0020 -0.0354 -0.2225 -0.4277 0 0.4277 0.2225 0.0354 0.0020]; %odd symmetry

%% 13 Tap Inphase
g1=[0.0017265 0.01827 0.10551 0.30359 0.30275 -0.32046 -0.8230  -0.32046 0.30275 0.30359 0.10551 0.01827 0.0017265];
g2=[0.0001234 0.00193 0.018315 0.10540 0.36788 0.77880 1.0000 0.77880 0.36788 0.10540 0.018315 0.00193 0.0001234];
g3=[-0.00060937 -0.007943 -0.06029 -0.26022 -0.60550 -0.64092 0 0.64092 0.60550 0.26022 0.06029 0.007943 0.00060937];
g4=[-0.0003702 -0.004826 -0.03663 -0.158099 -0.36788 -0.38940 0 0.38940 0.36788 0.158099 0.03663 0.004826 0.0003702]; 

g1=single(g1);
g2=single(g2);
g3=single(g3);
g4=single(g4);

% for i=1:length(g3)
% G_yt(i,:)=g2*g3(i);
% end
% for j=1:length(g3)
% G(:,:,j)=g1'*G_yt(j,:);
% end

%% Inphase
for i=1:length(g3)
G_yt(i,:)=g2*g2(i);
end
for j=1:length(g3)
G_2a(:,:,j)=g1'*G_yt(j,:);
end

for i=1:length(g3)
G_yt(i,:)=g4*g2(i);
end
for j=1:length(g3)
G_2b(:,:,j)=g3'*G_yt(j,:);
end

for i=1:length(g3)
G_yt(i,:)=g1*g2(i);
end
for j=1:length(g3)
G_2c(:,:,j)=g2'*G_yt(j,:);
end

for i=1:length(g3)
G_yt(i,:)=g2*g4(i);
end
for j=1:length(g3)
G_2d(:,:,j)=g3'*G_yt(j,:);
end

for i=1:length(g3)
G_yt(i,:)=g3*g4(i);
end
for j=1:length(g3)
G_2e(:,:,j)=g2'*G_yt(j,:);
end

for i=1:length(g3)
G_yt(i,:)=g2*g1(i);
end
for j=1:length(g3)
G_2f(:,:,j)=g2'*G_yt(j,:);
end

% %% 9 Tap Quadrature Phase
% h1=[-0.0088 -0.0554 0.0895 0.6776 0 -0.6776 -0.0895 0.0554 0.0088]; %odd
% h2=[0.0043 0.0508 0.1522 -0.1695 -0.6595 -0.1695 0.1522 0.0508 0.0043]; %even
% h3=[0.0008 0.0176 0.1660 0.6383 1.0000 0.6383 0.1660 0.0176 0.0008]; %even
% h4=[-0.0018 -0.0310 -0.1953 -0.3754 0 0.3754 0.1953 0.0310 0.0018];
% h5=[-0.0020 -0.0354 -0.2225 -0.4277 0 0.4277 0.2225 0.0354 0.0020];

%% 13 Tap Quadrature Phase

h1=[-0.002192 -0.016928 -0.05614 0.000551 0.40493 0.68498 0 -0.68498 -0.40493 -0.000551 0.05614 0.016928 0.002192];
h2=[0.0008935 0.0093175 0.05223 0.13865 0.080298 -0.342717 -0.6595 -0.342717 0.080298 0.13865 0.05223  0.0093175 0.0008935];
h3=[0.0001234 0.0019304 0.018315 0.105399 0.36788 0.77880 1 0.77880 0.36788 0.105399 0.018315 0.0019304 0.0001234];
h4=[-0.000325 -0.0042302 -0.032154 -0.138775 -0.322915 -0.341806 0 0.341806 0.322915 0.138775 0.032154 0.0042302 0.000325];
h5=[-0.0003702 -0.004826 -0.03663 -0.158099 -0.367878 -0.38940 0 0.38940 0.367878 0.158099 0.03663 0.004826 0.0003702];

h1=single(h1);
h2=single(h2);
h3=single(h3);
h4=single(h4);
h5=single(h5);

for i=1:length(g3)
G_yt(i,:)=h3*h3(i);
end
for j=1:length(g3)
H_2a(:,:,j)=h1'*G_yt(j,:);
end

for i=1:length(g3)
G_yt(i,:)=h5*h3(i);
end
for j=1:length(g3)
H_2b(:,:,j)=h2'*G_yt(j,:);
end

for i=1:length(g3)
G_yt(i,:)=h2*h3(i);
end
for j=1:length(g3)
H_2c(:,:,j)=h5'*G_yt(j,:);
end

for i=1:length(g3)
G_yt(i,:)=h1*h3(i);
end
for j=1:length(g3)
H_2d(:,:,j)=h3'*G_yt(j,:);
end


for i=1:length(g3)
G_yt(i,:)=h3*h5(i);
end
for j=1:length(g3)
H_2e(:,:,j)=h2'*G_yt(j,:);
end

for i=1:length(g3)
G_yt(i,:)=h5*h5(i);
end
for j=1:length(g3)
H_2f(:,:,j)=h4'*G_yt(j,:);
end

for i=1:length(g3)
G_yt(i,:)=h2*h5(i);
end
for j=1:length(g3)
H_2g(:,:,j)=h3'*G_yt(j,:);
end

for i=1:length(g3)
G_yt(i,:)=h3*h2(i);
end
for j=1:length(g3)
H_2h(:,:,j)=h5'*G_yt(j,:);
end


for i=1:length(g3)
G_yt(i,:)=h5*h2(i);
end
for j=1:length(g3)
H_2i(:,:,j)=h3'*G_yt(j,:);
end
for i=1:length(g3)
G_yt(i,:)=h3*h1(i);
end
for j=1:length(g3)
H_2j(:,:,j)=h3'*G_yt(j,:);
end

%% Lowpass
lo = single((1/16)*[1 4 6 4 1]);

for i=1:length(lo)
G_low(i,:)=lo*lo(i);
end
for j=1:length(lo)
LP(:,:,j)=lo'*G_low(j,:);
end

for i=5:14
    I=imread(strcat('../Foreman_gray_frames/foreman',num2str(i-4),'.bmp'));
    I=histeq(I);
    V(:,:,i)=mat2gray(I);
end

%% Adding pad for first and last frames
I=imread(strcat('../Foreman_gray_frames/foreman1.bmp'));;
I=histeq(I);
V(:,:,1)=mat2gray(I);
V(:,:,2)=mat2gray(I);
V(:,:,3)=mat2gray(I);
V(:,:,4)=mat2gray(I);
I=imread(strcat('../Foreman_gray_frames/foreman10.bmp'));
I=histeq(I);
V(:,:,15)=mat2gray(I);
V(:,:,16)=mat2gray(I);
V(:,:,17)=mat2gray(I);
V(:,:,18)=mat2gray(I);


%% Level 0
V_G2a=convn(V,G_2a,'same');
V_G2b=convn(V,G_2b,'same');
V_G2c=convn(V,G_2c,'same');
V_G2d=convn(V,G_2d,'same');
V_G2e=convn(V,G_2e,'same');
V_G2f=convn(V,G_2f,'same');

V_H2a=convn(V,H_2a,'same');
V_H2b=convn(V,H_2b,'same');
V_H2c=convn(V,H_2c,'same');
V_H2d=convn(V,H_2d,'same');
V_H2e=convn(V,H_2e,'same');
V_H2f=convn(V,H_2f,'same');
V_H2g=convn(V,H_2g,'same');
V_H2h=convn(V,H_2h,'same');
V_H2i=convn(V,H_2i,'same');
V_H2j=convn(V,H_2j,'same');


tic

k=1;
for theta=0:20:360
    theta_r=theta*pi/180;
    for phi=-90:20:90
        phi_r=phi*pi/180;
        alpha=cos(theta)*sin(phi);
        beta=sin(theta)*sin(phi);
        gamma=cos(phi);
        G=alpha^2*V_G2a+2*alpha*beta*V_G2b+beta^2*V_G2c+2*alpha*gamma*V_G2d+2*beta*gamma*V_G2e+gamma^2*V_G2f;
        H=alpha^3*V_H2a+ 3*alpha^2*beta*V_H2b+ 3*alpha*beta^2*V_H2c+ beta^3*V_H2d+ 3*alpha^2*gamma*V_H2e+ 6*alpha*beta*gamma*V_H2f+ 3*beta^2*gamma*V_H2g+ 3*alpha*gamma^2*V_H2h+ 3*beta*gamma^2*V_H2i+ gamma^3*V_H2j;
        V0{k}=complex(G,H);
        k=k+1;
    end
end

V=convn(V,LP,'same');

V=V(1:2:end,1:2:end,:);

%% Level1
V_G2a=convn(V,G_2a,'same');
V_G2b=convn(V,G_2b,'same');
V_G2c=convn(V,G_2c,'same');
V_G2d=convn(V,G_2d,'same');
V_G2e=convn(V,G_2e,'same');
V_G2f=convn(V,G_2f,'same');

V_H2a=convn(V,H_2a,'same');
V_H2b=convn(V,H_2b,'same');
V_H2c=convn(V,H_2c,'same');
V_H2d=convn(V,H_2d,'same');
V_H2e=convn(V,H_2e,'same');
V_H2f=convn(V,H_2f,'same');
V_H2g=convn(V,H_2g,'same');
V_H2h=convn(V,H_2h,'same');
V_H2i=convn(V,H_2i,'same');
V_H2j=convn(V,H_2j,'same');




k=1;
for theta=0:20:360
    theta_r=theta*pi/180;
    for phi=-90:20:90
        phi_r=phi*pi/180;
        alpha=cos(theta)*sin(phi);
        beta=sin(theta)*sin(phi);
        gamma=cos(phi);
        G=alpha^2*V_G2a+2*alpha*beta*V_G2b+beta^2*V_G2c+2*alpha*gamma*V_G2d+2*beta*gamma*V_G2e+gamma^2*V_G2f;
        H=alpha^3*V_H2a+ 3*alpha^2*beta*V_H2b+ 3*alpha*beta^2*V_H2c+ beta^3*V_H2d+ 3*alpha^2*gamma*V_H2e+ 6*alpha*beta*gamma*V_H2f+ 3*beta^2*gamma*V_H2g+ 3*alpha*gamma^2*V_H2h+ 3*beta*gamma^2*V_H2i+ gamma^3*V_H2j;
        V1{k}=complex(G,H);
        k=k+1;
    end
end

V=convn(V,LP,'same');

V=V(1:2:end,1:2:end,:);

%% Level 2
V_G2a=convn(V,G_2a,'same');
V_G2b=convn(V,G_2b,'same');
V_G2c=convn(V,G_2c,'same');
V_G2d=convn(V,G_2d,'same');
V_G2e=convn(V,G_2e,'same');
V_G2f=convn(V,G_2f,'same');

V_H2a=convn(V,H_2a,'same');
V_H2b=convn(V,H_2b,'same');
V_H2c=convn(V,H_2c,'same');
V_H2d=convn(V,H_2d,'same');
V_H2e=convn(V,H_2e,'same');
V_H2f=convn(V,H_2f,'same');
V_H2g=convn(V,H_2g,'same');
V_H2h=convn(V,H_2h,'same');
V_H2i=convn(V,H_2i,'same');
V_H2j=convn(V,H_2j,'same');

k=1;
for theta=0:20:360
    theta_r=theta*pi/180;
    for phi=-90:20:90
        phi_r=phi*pi/180;
        alpha=cos(theta)*sin(phi);
        beta=sin(theta)*sin(phi);
        gamma=cos(phi);
        G=alpha^2*V_G2a+2*alpha*beta*V_G2b+beta^2*V_G2c+2*alpha*gamma*V_G2d+2*beta*gamma*V_G2e+gamma^2*V_G2f;
        H=alpha^3*V_H2a+ 3*alpha^2*beta*V_H2b+ 3*alpha*beta^2*V_H2c+ beta^3*V_H2d+ 3*alpha^2*gamma*V_H2e+ 6*alpha*beta*gamma*V_H2f+ 3*beta^2*gamma*V_H2g+ 3*alpha*gamma^2*V_H2h+ 3*beta*gamma^2*V_H2i+ gamma^3*V_H2j;
        V2{k}=complex(G,H);
        k=k+1;
    end
end

toc;        

% save V0 V0
% save V1 V1
% save V2 V2

%% Dominant orientation

 k=1;
[m2 n2 p2]=size(V2{1});
e2=zeros(m2,n2,p2);
[m1 n1 p1]=size(V1{1});
e1=zeros(m1,n1,p1);
[m0 n0 p0]=size(V0{1});
e0=zeros(m0,n0,p0);

for theta=0:20:360
    for phi=-90:20:90
edge_level2=abs(V2{k});
edge_level1=abs(V1{k});
edge_level0=abs(V0{k});
e2(e2<edge_level2)=edge_level2(e2<edge_level2);
e1(e1<edge_level1)=edge_level1(e1<edge_level1);
e0(e0<edge_level0)=edge_level0(e0<edge_level0);
k=k+1;
    end
end

e2=e2./max(e2(:));
e1=e1./max(e1(:));
e0=e0./max(e0(:));

g2=e2(:,:,5:14);
g1=e1(:,:,5:14);
g0=e0(:,:,5:14);

% x=V0{65};
% implay(abs(x));
toc;

[m n p]=size(g0);
for k=1:p
    g0(:,:,k)=histeq(g0(:,:,k));
    g1(:,:,k)=histeq(g1(:,:,k));
    g2(:,:,k)=histeq(g2(:,:,k));
end

save g0_fore g0
save g1_fore g1
save g2_fore g2
