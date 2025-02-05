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
clc
close all
clear all

[ima, map] = imread('MRI_pseudo_colored.jpg'); 
ima=double(ima);

mask3=(1/9)*[ones(3,3)];
mask5=(1/25)*[ones(5,5)];
mask7=(1/49)*[ones(7,7)];

ima_res3=imfilter(ima,mask3);
ima_res5=imfilter(ima,mask5);
ima_res7=imfilter(ima,mask7);

e_ima_res3=calcular_energia(ima_res3);
e_ima_res5=calcular_energia(ima_res5);
e_ima_res7=calcular_energia(ima_res7);

figure('Name', 'comparación suavizado de imágenes');
subplot(1,3,1);imagesc(uint8(ima_res3));title(sprintf('Ima suavT, X; E=%g', e_ima_res3)); colorbar;
subplot(1,3,2);imagesc(uint8(ima_res5));title(sprintf('Ima suavT, X; E=%g', e_ima_res5)); colorbar;
subplot(1,3,3);imagesc(uint8(ima_res7));title(sprintf('Ima suavT, X; E=%g', e_ima_res7)); colorbar;

r1 = double(ima(:,:,1)); 
g1 = double(ima(:,:,2)); 
b1 = double(ima(:,:,3)); 

r23 = double(ima_res3(:,:,1)); 
g23 = double(ima_res3(:,:,2)); 
b23 = double(ima_res3(:,:,3)); 

r25 = double(ima_res5(:,:,1)); 
g25 = double(ima_res5(:,:,2)); 
b25 = double(ima_res5(:,:,3)); 

r27 = double(ima_res7(:,:,1)); 
g27 = double(ima_res7(:,:,2)); 
b27 = double(ima_res7(:,:,3)); 


rdif3 = (r1-r23).^2; 
gdif3 = (g1-g23).^2; 
bdif3 = (b1-b23).^2; 

rdif5 = (r1-r25).^2; 
gdif5 = (g1-g25).^2; 
bdif5 = (b1-b25).^2; 

rdif7 = (r1-r27).^2; 
gdif7 = (g1-g27).^2; 
bdif7 = (b1-b27).^2; 

e_ima_rdif3=calcular_energia(rdif3);
e_ima_gdif3=calcular_energia(gdif3);
e_ima_bdif3=calcular_energia(bdif3);
e_ima_rdif5=calcular_energia(rdif5);
e_ima_gdif5=calcular_energia(gdif5);
e_ima_bdif5=calcular_energia(bdif5);
e_ima_rdif7=calcular_energia(rdif7);
e_ima_gdif7=calcular_energia(gdif7);
e_ima_bdif7=calcular_energia(bdif7);


figure('Name', 'Ejercicio 2: comparación suavizado de imágenes. Imágenes diferencia (una por canal)');
subplot(3,3,1); imagesc(rdif3); title(sprintf('Ima suavT, X; E=%g', e_ima_rdif3)); colorbar;
subplot(3,3,2); imagesc(gdif3); title(sprintf('Im suav, sim; E=%g', e_ima_gdif3)); colorbar;
subplot(3,3,3); imagesc(bdif3); title(sprintf('Ima suav, repl; E=%g', e_ima_bdif3)); colorbar;
subplot(3,3,4); imagesc(rdif5); title(sprintf('Ima suavT, X; E=%g', e_ima_rdif5)); colorbar;
subplot(3,3,5); imagesc(gdif5); title(sprintf('Im suav, sim; E=%g', e_ima_gdif5)); colorbar;
subplot(3,3,6); imagesc(bdif5); title(sprintf('Ima suav, repl; E=%g', e_ima_bdif5)); colorbar;
subplot(3,3,7); imagesc(rdif7); title(sprintf('Ima suavT, X; E=%g', e_ima_rdif7)); colorbar;
subplot(3,3,8); imagesc(gdif7); title(sprintf('Im suav, sim; E=%g', e_ima_gdif7)); colorbar;
subplot(3,3,9); imagesc(bdif7); title(sprintf('Ima suav, repl; E=%g', e_ima_bdif7)); colorbar;

%% Utilización de filtros binomiales
clc
close all
clear all

[ima, map] = imread('MRI_pseudo_colored.jpg'); 
ima=double(ima);

ima_res3bi=imfilter_binomial(ima,3);
ima_res5bi=imfilter_binomial(ima,5);
ima_res7bi=imfilter_binomial(ima,7);

e_ima_res3bi=calcular_energia(ima_res3bi);
e_ima_res5bi=calcular_energia(ima_res5bi);
e_ima_res7bi=calcular_energia(ima_res7bi);

figure('Name', 'comparación suavizado de imágenes. Filtro binomial.');
subplot(1,3,1);imagesc(uint8(ima_res3bi));title(sprintf('Ima suavT, X; E=%g', e_ima_res3bi)); colorbar;
subplot(1,3,2);imagesc(uint8(ima_res5bi));title(sprintf('Ima suavT, X; E=%g', e_ima_res5bi)); colorbar;
subplot(1,3,3);imagesc(uint8(ima_res7bi));title(sprintf('Ima suavT, X; E=%g', e_ima_res7bi)); colorbar;

r1 = double(ima(:,:,1)); 
g1 = double(ima(:,:,2)); 
b1 = double(ima(:,:,3)); 

r23bi = double(ima_res3bi(:,:,1)); 
g23bi = double(ima_res3bi(:,:,2)); 
b23bi = double(ima_res3bi(:,:,3)); 

r25bi = double(ima_res5bi(:,:,1)); 
g25bi = double(ima_res5bi(:,:,2)); 
b25bi = double(ima_res5bi(:,:,3)); 

r27bi = double(ima_res7bi(:,:,1)); 
g27bi = double(ima_res7bi(:,:,2)); 
b27bi = double(ima_res7bi(:,:,3)); 


rdif3bi = (r1-r23bi).^2; 
gdif3bi = (g1-g23bi).^2; 
bdif3bi = (b1-b23bi).^2; 

rdif5bi = (r1-r25bi).^2; 
gdif5bi = (g1-g25bi).^2; 
bdif5bi = (b1-b25bi).^2; 

rdif7bi = (r1-r27bi).^2; 
gdif7bi = (g1-g27bi).^2; 
bdif7bi = (b1-b27bi).^2; 

e_ima_rdif3bi=calcular_energia(rdif3bi);
e_ima_gdif3bi=calcular_energia(gdif3bi);
e_ima_bdif3bi=calcular_energia(bdif3bi);
e_ima_rdif5bi=calcular_energia(rdif5bi);
e_ima_gdif5bi=calcular_energia(gdif5bi);
e_ima_bdif5bi=calcular_energia(bdif5bi);
e_ima_rdif7bi=calcular_energia(rdif7bi);
e_ima_gdif7bi=calcular_energia(gdif7bi);
e_ima_bdif7bi=calcular_energia(bdif7bi);


figure('Name', 'Ejercicio 2: comparación suavizado de imágenes. Filtro binomial. Imágenes diferencia (una por canal)');
subplot(3,3,1); imagesc(rdif3bi); title(sprintf('Ima suavT, X; E=%g', e_ima_rdif3bi)); colorbar;
subplot(3,3,2); imagesc(gdif3bi); title(sprintf('Im suav, sim; E=%g', e_ima_gdif3bi)); colorbar;
subplot(3,3,3); imagesc(bdif3bi); title(sprintf('Ima suav, repl; E=%g', e_ima_bdif3bi)); colorbar;
subplot(3,3,4); imagesc(rdif5bi); title(sprintf('Ima suavT, X; E=%g', e_ima_rdif5bi)); colorbar;
subplot(3,3,5); imagesc(gdif5bi); title(sprintf('Im suav, sim; E=%g', e_ima_gdif5bi)); colorbar;
subplot(3,3,6); imagesc(bdif5bi); title(sprintf('Ima suav, repl; E=%g', e_ima_bdif5bi)); colorbar;
subplot(3,3,7); imagesc(rdif7bi); title(sprintf('Ima suavT, X; E=%g', e_ima_rdif7bi)); colorbar;
subplot(3,3,8); imagesc(gdif7bi); title(sprintf('Im suav, sim; E=%g', e_ima_gdif7bi)); colorbar;
subplot(3,3,9); imagesc(bdif7bi); title(sprintf('Ima suav, repl; E=%g', e_ima_bdif7bi)); colorbar;

%% Extracción de bordes
% Lea la imagen ‘Hernia_disco.jpg4. Obtenga los bordes de la imagen utilizando los filtros de Prewitt
% 3x3. Primero debe calcular el gradiente en la dirección x e y. Posteriormente obtenga la intensidad de bordes
% total combinando ambos gradientes. Visualice la imagen original, las imágenes de intensidad de los
% gradientes calculados y la imagen de intensidad de bordes total.
% Se recomienda representar las imágenes de gradiente y bordes con imagesc y un mapa de colores en
% escala de grises. Para ello, utilice la instrucción colormap(gray) antes o después de la instrucción imagesc.
% Opere siempre con las imágenes de tipo double. Calcule el gradiente G total como:
% G =sqrt(Gx.^2 + Gy.^2) % Gx Gradiente en x / Gy Gradiente en y
%% Prewitt
clc
close all
clear all

[I,map]=imread("Hernia_disco.jpg");
%I=double(I);
Prewsd=[1 1 1; 0 0 0; -1 -1 -1]/6;%Prewsd=double(Prewsd);
PreswdT=[-1 0 1; -1 0 1; -1 0 1]/6;%PreswdT=double(PreswdT);

% Sea I la imagen a filtrar y s y s2 los filtros de contornos, la intensidad de borde se calcula como el
% módulo de las dos componentes: 
FX=imfilter(I,Prewsd);FX=double(FX);
FY=imfilter(I,PreswdT);FY=double(FY);
F=sqrt(FX.^2 + FY.^2);

e_F=calcular_energia(F);
e_FX=calcular_energia(FX);
e_FY=calcular_energia(FY);

figure('Name', 'Ejercicio 3: Extracción de bordes. Prewitt');
subplot(1,3,1); colormap(gray);imagesc(uint8(F)); title(sprintf('Ima suavT, X; E=%g', e_F)); colorbar;
subplot(1,3,2); colormap(gray);imagesc(uint8(FX)); title(sprintf('Im suav, sim; E=%g', e_FX)); colorbar;
subplot(1,3,3); colormap(gray);imagesc(uint8(FY)); title(sprintf('Ima suav, repl; E=%g', e_FY)); colorbar;

%% Sobel 
clc
close all
clear all

[I,map]=imread("Hernia_disco.jpg");
I=double(I);
Sobeld=[1 2 1; 0 0 0; -1 -2 -1]/8;Sobeld=double(Sobeld);
SobeldT=[-1 0 1; -2 0 2; -1 0 1]/8;SobeldT=double(SobeldT);

% Sea I la imagen a filtrar y s y s2 los filtros de contornos, la intensidad de borde se calcula como el
% módulo de las dos componentes: 
FX=imfilter(I,Sobeld);FX=double(FX);
FY=imfilter(I,SobeldT);FY=double(FY);
F=sqrt(FX.^2 + FY.^2);

e_F=calcular_energia(F);
e_FX=calcular_energia(FX);
e_FY=calcular_energia(FY);

figure('Name', 'Ejercicio 4: Extracción de bordes. Sobel');
subplot(1,3,1); colormap(gray);imagesc(uint8(F)); title(sprintf('Ima suavT, X; E=%g', e_F)); colorbar;
subplot(1,3,2); colormap(gray);imagesc(uint8(FX)); title(sprintf('Im suav, sim; E=%g', e_FX)); colorbar;
subplot(1,3,3); colormap(gray);imagesc(uint8(FY)); title(sprintf('Ima suav, repl; E=%g', e_FY)); colorbar;

%% Roberts
clc
close all
clear all

[I,map]=imread("Hernia_disco.jpg");
I=double(I);
Robertsd=[1 0; 0 -1];Robertsd=double(Robertsd);
RobertsdT=[0 1; -1 0];RobertsdT=double(RobertsdT);

% Sea I la imagen a filtrar y s y s2 los filtros de contornos, la intensidad de borde se calcula como el
% módulo de las dos componentes: 
FX=imfilter(I,Robertsd);FX=double(FX);
FY=imfilter(I,RobertsdT);FY=double(FY);
F=sqrt(FX.^2 + FY.^2);

e_F=calcular_energia(F);
e_FX=calcular_energia(FX);
e_FY=calcular_energia(FY);

figure('Name', 'Ejercicio 4: Extracción de bordes. Roberts');
subplot(1,3,1); colormap(gray);imagesc(uint8(F)); title(sprintf('Ima suavT, X; E=%g', e_F)); colorbar;
subplot(1,3,2); colormap(gray);imagesc(uint8(FX)); title(sprintf('Im suav, sim; E=%g', e_FX)); colorbar;
subplot(1,3,3); colormap(gray);imagesc(uint8(FY)); title(sprintf('Ima suav, repl; E=%g', e_FY)); colorbar;




function energia = calcular_energia(imagen)

imagen=double(imagen); % para evitar desbordamientos en caso de unit, logical, ...
energia = sum(sum(imagen .* imagen));

% Si la imagen tiene múltiples canales (por ejemplo, RGB), sumar las energías de cada canal
if size(imagen, 3) > 1
    energia = sum(energia(:));
end

end


function ima_res=imfilter_binomial(ima,orden)

triang = cell(1,orden);

triang{1} = [1 1];

for i=2:orden
    v_ant = triang{i-1};
    v=ones(1, i+1);
    v(1)=1;
    v(end)=1;
    for j=2:length(v)-1
        v(j)=v_ant(j-1) + v_ant(j);
    end
    triang{i}=v;
end
    
v= triang{orden};

% Mascara del filtro
v=double(v);
W = v'*v;

% Normalizacion
C = sum(sum(W));
W = W * (1/C);

% Filtrado
ima_res = imfilter(ima,W);
end