clc
clear all
close all
%% Comparación suavizado de imágenes 
[ima,MAP]=imread("MRI_pseudo_colored.jpg");
double(ima);

mask3=(1/9)*[ones(3,3)];
mask5=(1/25)*[ones(5,5)];
mask7=(1/49)*[ones(7,7)];

ima_res3=imfilter(ima,mask3);
ima_res5=imfilter(ima,mask5);
ima_res7=imfilter(ima,mask7);
figure(1)
subplot(1,3,1);imagesc(uint8(ima_res3));
subplot(1,3,2);imagesc(uint8(ima_res5));
subplot(1,3,3);imagesc(uint8(ima_res7));


