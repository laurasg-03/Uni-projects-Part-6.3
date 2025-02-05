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
