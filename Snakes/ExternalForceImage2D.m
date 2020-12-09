function Eextern = ExternalForceImage2D(I,Wline, Wedge, Wterm,Sigma,verbose)
% Eextern = ExternalForceImage2D(I,Wline, Wedge, Wterm,Sigma)
% 
% inputs, 
%  I : The image
%  Sigma : Sigma used to calculated image derivatives 
%  Wline : Attraction to lines, if negative to black lines otherwise white
%          lines
%  Wedge : Attraction to edges
%  Wterm : Attraction to terminations of lines (end points) and corners
%
% outputs,
%  Eextern : The energy function described by the image
%
% Function is written by D.Kroon University of Twente (July 2010)

Ix=ImageDerivatives2D(I,Sigma,'x');
Iy=ImageDerivatives2D(I,Sigma,'y');
Ixx=ImageDerivatives2D(I,Sigma,'xx');
Ixy=ImageDerivatives2D(I,Sigma,'xy');
Iyy=ImageDerivatives2D(I,Sigma,'yy');


Eline = imgaussian(I,Sigma);
Eterm = (Iyy.*Ix.^2 -2*Ixy.*Ix.*Iy + Ixx.*Iy.^2)./((1+Ix.^2 + Iy.^2).^(3/2));
Eedge = -sqrt(Ix.^2 + Iy.^2); 

% Eline=mat2gray(Eline);
% Eterm=mat2gray(Eterm);
% Eedge=mat2gray(Eedge);

Eextern= (Wline*Eline + Wedge*Eedge -Wterm * Eterm); 

if verbose
%     figure
%     subplot(1,4,1)
%     imagesc(Wline*Eline), colormap(gray), title('E_{line}','Interpreter','latex'), axis off
%     subplot(1,4,2)
%     imagesc( -Wterm *Eterm), colormap(gray), title('E_{term}'), axis off
%     subplot(1,4,3)
%     imagesc(Wedge*Eedge), colormap(gray), title('E_{edge}'), axis off
%     subplot(1,4,4)
%     imagesc(Eextern), colormap(gray), title('E_{img}'), axis off


   figure
    subplot(2,2,1)
    imagesc(Wline*Eline), colormap(gray), title('E_{line}','Interpreter','latex'), axis off
    axis equal
    subplot(2,2,2)
    imagesc(-Wterm *Eterm), colormap(gray), title('E_{term}'), axis off
    axis equal
    subplot(2,2,3)
    imagesc(Wedge*Eedge), colormap(gray), title('E_{edge}'), axis off
    axis equal
    subplot(2,2,4)
    imagesc(Eextern), colormap(gray), title('E_{img}'), axis off
    axis equal
%     
%    figure
%    imagesc(Eextern), colormap(gray), 
%    axis off
%    axis equal
%     
%   
%    figure
%    imshow(Eextern)
end