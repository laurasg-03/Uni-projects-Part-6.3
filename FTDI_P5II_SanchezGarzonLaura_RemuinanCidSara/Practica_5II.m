%% Ejercicio 1: Dilatación y erosión por correlación con el elemento estructurante
clc
clear all
close all

[ima,map]=imread("bandas.bmp");
mask=[1 1 1 1 1;1 1 1 1 1; 1 1 1 0 0;1 1 1 0 0; 1 1 1 0 0];
mask=uint8(mask);

ima_res=imfilter_dilate(ima,mask);
e_ima = calcular_energia(ima);
e_ima_res_dil = calcular_energia(ima_res);
figure('Name', 'Aplicación de dilatador');
subplot(2,1,1); imshow(ima); title(sprintf('Ima original; E=%g', e_ima)); colorbar;
subplot(2,1,2); imshow(ima_res); title(sprintf('Ima suav, X; E=%g', e_ima_res_dil)); colorbar;

ima_res_erode=imfilter_erode(ima,mask);
e_ima_res_erode = calcular_energia(ima_res_erode);
figure('Name', 'Aplicación de erosionador');
subplot(2,1,1); imshow(ima); title(sprintf('Ima original; E=%g', e_ima)); colorbar;
subplot(2,1,2); imshow(ima_res_erode); title(sprintf('Ima suav, X; E=%g', e_ima_res_erode)); colorbar;

%% Ejercicio 2: Dilatación y erosión por desplazamiento de la señal

[ima,map]=imread("bandas.bmp");
mask=[1 1 1 1 1;1 1 1 1 1; 1 1 1 0 0;1 1 1 0 0; 1 1 1 0 0];
mask=uint8(mask);

ima_resD=imfilter_dilateD(ima,mask);
e_ima = calcular_energia(ima);
e_ima_res_dilD = calcular_energia(ima_resD);
figure('Name', 'Aplicación de dilatador por desplazamiento de señal');
subplot(2,1,1); imshow(ima); title(sprintf('Ima original; E=%g', e_ima)); colorbar;
subplot(2,1,2); imshow(ima_resD); title(sprintf('Ima suav, X; E=%g', e_ima_res_dilD)); colorbar;

ima_res_erodeD=imfilter_erodeD(ima,mask);
e_ima_res_erodeD = calcular_energia(ima_res_erodeD);
figure('Name', 'Aplicación de erosionador por desplazamiento');
subplot(2,1,1); imshow(ima); title(sprintf('Ima original; E=%g', e_ima)); colorbar;
subplot(2,1,2); imshow(ima_res_erodeD); title(sprintf('Ima suav, X; E=%g', e_ima_res_erodeD)); colorbar;

% Comparar energías
x=e_ima_res_dil-e_ima_res_dilD;
y=e_ima_res_erode-e_ima_res_erodeD;

% Comparar tiempos
tic
for i = 1:10
    ima_res=imfilter_dilate(ima,mask);
end
toc

tic
for i = 1:10
    ima_resD=imfilter_dilateD(ima,mask);
end
toc

tic
for i = 1:10
    ima_res=imfilter_erode(ima,mask);
end
toc

tic
for i = 1:10
    ima_resD=imfilter_erodeD(ima,mask);
end
toc

%% Gradiente morfológico
clc
clear all
close all

[ima,map]=imread("Tools.bmp");
se = strel('square',3); % creación de un elemento estructurante 3x3 cuadrado
ima_d=imdilate(ima,se); % dilatacion de la imagen ima
ima_e=imerode (ima,se); % erosion de la imagen ima

% Gradiente por dilatación
G_dil=ima_d-ima;
G_er=ima-ima_e;
G_morf=ima_d-ima_e;

e_ima = calcular_energia(ima);
e_G_dil = calcular_energia(G_dil);
e_G_er = calcular_energia(G_er);
e_G_morf = calcular_energia(G_morf);

figure('Name', 'Gradiente morfológico');
subplot(2,2,1); imshow(ima); title(sprintf('Ima original; E=%g', e_ima)); colorbar;
subplot(2,2,2); imshow(G_dil); title(sprintf('Gradiente por dilatación; E=%g', e_G_dil)); colorbar;
subplot(2,2,3); imshow(G_er); title(sprintf('Gradiente por erosión; E=%g', e_G_er)); colorbar;
subplot(2,2,4); imshow(G_morf); title(sprintf('Gradiente morfológico; E=%g', e_G_morf)); colorbar;

%% Restauración de un cuadro por apertura
% ??? Tiene sentido usar la apertura?
clc
clear all
close all

[ima,map]=imread("MRI_gray_garab.jpg");
se = strel('sphere', 2); 
ima_e = imerode(ima,se); G_er=ima-ima_e; apertura=imdilate(G_er,se);
imshow(apertura)

%%
clc
clear all
close all

[ima,map]=imread("verjanegra.jpg");
se = strel('line', 20, 180); 
ima_d = imdilate(ima,se); cierre=imerode(ima_d,se); 
imshow(cierre)



%%
function ima_res=imfilter_dilate(ima,mask)

ima=uint8(ima); % Para asegurarse de que el tipo es uint8
[image_h,image_w,~] = size(ima);

[M,N] = size(mask);
margen_y = (M-1) / 2;
margen_x = (N-1) / 2;

mask=uint8(255*fliplr(flipud(mask))); 	% En la función ‘imfilter_dilate’-->cambiar para imfilter_erode

ima_res = uint8(zeros(image_h,image_w));

for f=margen_y+1:image_h-margen_y,
    for c=margen_x+1:image_w-margen_x,
        
        simage = ima( f-margen_y:f+margen_y , c-margen_x:c+margen_x);
        
        and_image=bitand(simage,mask);		% En la función ‘imfilter_dilate’ -->cambiar para imfilter_erode
        
        ima_res(f,c) = max(max(and_image)); % En la función ‘imfilter_dilate’-->cambiar para imfilter_erode

    end;
end;
% sustituir los valores no calculados
ima_res(1:margen_y,:) = ima(1:margen_y,:);
ima_res(end-margen_y+1:end,:) = ima(end-margen_y+1:end,:);
ima_res(:,1:margen_x) = ima(:,1:margen_x);
ima_res(:,end-margen_x+1:end) = ima(:,end-margen_x+1:end);
end

function ima_res=imfilter_erode(ima,mask)

ima=uint8(ima); % Para asegurarse de que el tipo es uint8
[image_h,image_w,~] = size(ima);

[M,N] = size(mask);
margen_y = (M-1) / 2;
margen_x = (N-1) / 2;


mask=uint8(255*(~mask)); 	% En la función ‘imfilter_dilate’-->cambiar para imfilter_erode

ima_res = uint8(zeros(image_h,image_w));

for f=margen_y+1:image_h-margen_y,
    for c=margen_x+1:image_w-margen_x,
        
        simage = ima( f-margen_y:f+margen_y , c-margen_x:c+margen_x);
        
        or_image=bitor(simage,mask);		% En la función ‘imfilter_dilate’ -->cambiar para imfilter_erode
        
        ima_res(f,c) = min(min(or_image)); % En la función ‘imfilter_dilate’-->cambiar para imfilter_erode

    end;
end;
% sustituir los valores no calculados
ima_res(1:margen_y,:) = ima(1:margen_y,:);
ima_res(end-margen_y+1:end,:) = ima(end-margen_y+1:end,:);
ima_res(:,1:margen_x) = ima(:,1:margen_x);
ima_res(:,end-margen_x+1:end) = ima(:,end-margen_x+1:end);
end

function energia = calcular_energia(imagen)

imagen=double(imagen); % para evitar desbordamientos en caso de unit, logical, ...
energia = sum(sum(imagen .* imagen));

% Si la imagen tiene múltiples canales (por ejemplo, RGB), sumar las energías de cada canal
if size(imagen, 3) > 1
    energia = sum(energia(:));
end

end

function ima_res=imfilter_dilateD(ima,mask)
ima=uint8(ima); 
[image_h,image_w,~] = size(ima);
[M,N] = size(mask);
% El centro del kernel
mask_center_y =ceil((M) / 2);
mask_center_x=ceil((N) / 2);
margen_y=(M-1)/2;
margen_x=(N-1)/2;
mask=uint8(255*fliplr(flipud(mask)));
ima_anterior=0.*ima;
for i=1:M
    for j=1:N
        if mask(i,j)==255
            shift_vector=[mask_center_y-i, mask_center_x-j];
            ima_despl = circshift(ima, shift_vector);
            ima_res=max(ima_anterior,ima_despl);
            ima_anterior=ima_res;
        end
    end
end
% sustituir los valores no calculados
ima_res(1:margen_y,:) = ima(1:margen_y,:);
ima_res(end-margen_y+1:end,:) = ima(end-margen_y+1:end,:);
ima_res(:,1:margen_x) = ima(:,1:margen_x);
ima_res(:,end-margen_x+1:end) = ima(:,end-margen_x+1:end);
end

function ima_res=imfilter_erodeD(ima,mask)
ima=uint8(ima); 
[image_h,image_w,~] = size(ima);
[M,N] = size(mask);
% El centro del kernel
mask_center_y =ceil((M) / 2);
mask_center_x=ceil((N) / 2);
margen_y=(M-1)/2;
margen_x=(N-1)/2;
mask=uint8(255*mask);
ima_anterior=uint8(255*ones(image_h,image_w));
for i=1:M
    for j=1:N
        if mask(i,j)==255
            shift_vector=[mask_center_y-i, mask_center_x-j];
            ima_despl = circshift(ima, shift_vector);
            ima_res=min(ima_anterior,ima_despl);
            ima_anterior=ima_res;
        end
    end
end
% sustituir los valores no calculados
ima_res(1:margen_y,:) = ima(1:margen_y,:);
ima_res(end-margen_y+1:end,:) = ima(end-margen_y+1:end,:);
ima_res(:,1:margen_x) = ima(:,1:margen_x);
ima_res(:,end-margen_x+1:end) = ima(:,end-margen_x+1:end);
end