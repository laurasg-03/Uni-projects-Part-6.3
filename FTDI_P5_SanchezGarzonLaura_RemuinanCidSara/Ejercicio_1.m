close all
clear all
clc

[ima, map] = imread('MRI_pseudo_colored.jpg'); 
ima=double(ima);
matriz=[zeros(1,7);zeros(1,7);zeros(1,7);7,6,4,1,4,6,7;zeros(1,7);zeros(1,7);zeros(1,7)]/35;
ima_res1=imfilter(ima,matriz); 
ima_res2=imfilter(ima,matriz'); 
ima_res_sym=imfilter(ima,matriz,'symmetric'); 
ima_res_rep=imfilter(ima,matriz,'replicate'); 
ima_res_circ=imfilter(ima,matriz,'circular'); 

e_ima=calcular_energia(ima);
e_ima_res1=calcular_energia(ima);
e_ima_res2=calcular_energia(ima);
e_ima_res_sym=calcular_energia(ima_res_sym);
e_ima_res_rep=calcular_energia(ima_res_rep);
e_ima_res_circ=calcular_energia(ima_res_circ);

figure('Name', 'Aplicación de operadores locales lineales');
subplot(3,3,2); imshow(uint8(ima)); title(sprintf('Ima original; E=%g', e_ima)); colorbar;
subplot(3,3,4); imshow(uint8(ima_res1)); title(sprintf('Ima suav, X; E=%g', e_ima_res1)); colorbar;
subplot(3,3,6); imshow(uint8(ima_res2)); title(sprintf('Ima suavT, X; E=%g', e_ima_res2)); colorbar;
subplot(3,3,7); imshow(uint8(ima_res_sym)); title(sprintf('Im suav, sim; E=%g', e_ima_res_sym)); colorbar;
subplot(3,3,8); imshow(uint8(ima_res_rep)); title(sprintf('Ima suav, repl; E=%g', e_ima_res_rep)); colorbar;
subplot(3,3,9); imshow(uint8(ima_res_circ)); title(sprintf('Ima suav, circ; E=%g', e_ima_res_circ)); colorbar;

%% Filtros habituales
%% Suavizado con filtro de media 
[ima, map] = imread('MRI_pseudo_colored.jpg'); 
ima=double(ima);
mask = [ones(1,3);ones(1,3);ones(1,3)]/9;
ima_res=imfilter(ima,mask); 
r1 = double(ima_res(:,:,1)); 
g1 = double(ima_res(:,:,2)); 
b1 = double(ima_res(:,:,3)); 

r2 = double(ima(:,:,1)); 
g2 = double(ima(:,:,2)); 
b2 = double(ima(:,:,3)); 

rdif = (r1-r2).^2; 
gdif = (g1-g2).^2; 
bdif = (b1-b2).^2; 

e_ima=calcular_energia(ima);
e_ima_res=calcular_energia(ima_res);
e_ima_rdif=calcular_energia(rdif);
e_ima_gdif=calcular_energia(gdif);
e_ima_bdif=calcular_energia(bdif);

figure('Name', 'Ejercicio 1: suavizado con filtro de media. Visualice la imagen original y la resultante tras aplicar el filtro');
subplot(2,1,1); imshow(uint8(ima)); title(sprintf('Ima original; E=%g', e_ima)); colorbar;
subplot(2,1,2); imshow(uint8(ima_res)); title(sprintf('Ima suav, X; E=%g', e_ima_res)); colorbar;

figure('Name', 'Ejercicio 1: suavizado con filtro de media. Imágenes diferencia (una por canal)');
subplot(1,3,1); imagesc(rdif); title(sprintf('Ima suavT, X; E=%g', e_ima_rdif)); colorbar;
subplot(1,3,2); imagesc(gdif); title(sprintf('Im suav, sim; E=%g', e_ima_gdif)); colorbar;
subplot(1,3,3); imagesc(bdif); title(sprintf('Ima suav, repl; E=%g', e_ima_bdif)); colorbar;

%% Comparación suavizado de imágenes
[ima, map] = imread('MRI_pseudo_colored.jpg'); 
ima=double(ima);
mask = [ones(1,3);ones(1,3);ones(1,3)]/9;
ima_res=imfilter(ima,mask); 




function energia = calcular_energia(imagen)

imagen=double(imagen); % para evitar desbordamientos en caso de unit, logical, ...
energia = sum(sum(imagen .* imagen));

end