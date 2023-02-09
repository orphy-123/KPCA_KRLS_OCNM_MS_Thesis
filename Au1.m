Images = abs(0.5 + 0.25.*randn(10,16)); %10 images in bmp form, each as a row vector
Novel = []; %Novel image indices
matchCounts = []; %Novel image matchCounts
Compress = []; %Novel image compressions
tau = 0.5; %threshold for similarity metric rho
matchFound = 0; %Boolean

for i=1:size(Images,1)
    
    I0 = Images(i,:) %Assign arriving image to I_nought
    matchFound = 0;
    I0 = I0 %Convert JPG image to BMP
    C0 = std(I0); %Compress image I0
    
    if i==1 %If I0 is first image
        Novel = [Novel; I0]; %Initialize novel image list
        matchFound = 1; %Similarity found
        matchCounts = [matchCounts; 1]; %Initialize matchCounts with first novel image
        Compress = [Compress; C0]; %Save C0
    else
        for n=1:size(Novel,1)
            I1 = Novel(n,:); %Assign n-th novel image to I_one
            C1 = std(I1); %Compress image I1
            I01 = [I0 I1]; %Concatenate (I0,I1)
            C01 = std(I01); %Compress concatenated image I01 = (I0,I1)
            
%            rho = (size(C0)+size(C1)-size(C01)) /(size(C0)+size(C1));
            rho = (C0+C1-C01)/(C0+C1);
            
            if rho<tau %dissimilar
            %i.e. arriving image is deemed novel when compared to n-th novel image,
            %so continue checking with rest of novel image list
                ;
            else %rho>=tau, meaning I1 from similar
                matchFound = 1;
                %matchCounts(n)=matchCounts(n)+1; %Increment matchCount for n-th novel image
                %Sort Novel by matchCount;
                break; %for n=1:size(Novel,1)
            end %if rho<tau
        end %for n=1:size(Novel,1)
    end %if i==1
    
    if matchFound==0 %I0 not first image, and found dissimilar to all novel images
        Novel = [Novel; I0];
        matchCounts = [matchCounts; 1];
        Compress = [Compress; C0]; %Save C0
    end %if matchFound==0
    
end %for i=1:size(Images,1)

