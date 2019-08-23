function mat = makeBowl(size)
%This creates a binary mask to black out unknown parts of an image.

mat = zeros(size);

[xGrid, yGrid] = meshgrid(1:size,1:size);
mat(sqrt((xGrid-((size+1)/2)).^2 + (yGrid-((size+1)/2)).^2) <=(size/2)+1)=1;
       