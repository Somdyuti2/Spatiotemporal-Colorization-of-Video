clc
clear all
close all

for k=1:10
    I=imread(strcat('../Foreman_gray_frames/foreman',num2str(k),'.bmp'));
%     I=imresize(I,0.5);
%     I=rgb2gray(I);
    I=histeq(I);
    E(:,:,k)=edge(I,'canny');
end

save E_fore E