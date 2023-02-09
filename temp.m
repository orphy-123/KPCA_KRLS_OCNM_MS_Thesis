Images = abs(0.5 + 0.25.*randn(10,16)); %10 images in bmp form, each as a row vector
Novel = []; %Novel image indices
Compress = []; %Novel image compressions

for i=1:size(Images,1)
    for j=1:size(Images,1)
        rho_numerator = std(Images(i,:)) + std(Images(j,:)) - std([Images(i,:) Images(j,:)]);
        rho_denominator = std(Images(i,:)) + std(Images(j,:));
        rho(i,j) = rho_numerator/rho_denominator;
    end
end

