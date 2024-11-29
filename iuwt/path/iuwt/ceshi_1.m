
load('exampleimages.mat');
% I=Original;
I=Dirtymap;
index = 1;

for sigma = 0.6:0.1:6
    
    gausFilter = fspecial('gaussian', [10, 10], sigma);
    
    adaptiveImg{index} = imfilter(I, gausFilter, 'replicate');
    
    index = index+1;
    
end

for index = 1:55
    
    adaptiveSub{index} = (I - adaptiveImg{index}).^2;
    
    temp = max(max(adaptiveImg{index}));
    
end

adaptiveI = I;

for i = 1:size(I, 1)
    
    for j = 1:size(I, 2)
        
        min = 100.0;
        
        sigma = 0.6;
        
        for index = 1:55
            
            temp = 1/sigma + adaptiveSub{index}(i,j);
            
            if temp < min
                
                min = temp;
                
                sigmaMap(i,j) = sigma;
                
                indexMap(i,j) = index;
                
                adaptiveI(i,j) = adaptiveImg{index}(i,j);
                
            end
            
            sigma = sigma + 0.1;
            
        end
        
    end
    
end

figure;

subplot(1,3,1);

imshow(I);

title('original image');

subplot(1,3,2)

imshow(adaptiveI);

title('adaptive image');

subplot(1,3,3);

imshow(adaptiveImg{20});

title('blur image with sigma = 2.5');


