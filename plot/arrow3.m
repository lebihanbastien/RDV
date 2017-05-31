function [ligne, pointe]=arrow3(Xi,Xf,col,lw,loF,laF)
% ARROW3 3D arrow for 3D MATLAB plots.
%
% Emilien Fabacher 2015
%
% Log:
% 03/2016 Minor modifications for Matlab 2013a compatibility

N = 20;
t = (0:30/N:360)'/180*pi;
X = [0; cos(t)];
Y = [0; sin(t)];
Z = [1; 0*X(1:end-1)];


% loF = norm(Xi-Xf);
X = X*laF*loF;
Y = Y*laF*loF;
Z = Z*loF-loF;

XYZ = rotation([X Y Z],vrrotvec2mat(vrrotvec([0 0 1],Xf-Xi)));
XYZ = XYZ+ones(size(XYZ,1),1)*(Xf);
XYZbase = XYZ(2:end,:);

vert = [XYZ];
fac = [ones(1,length(vert)-1); 2:length(vert); 3:length(vert) 2]';
fac = [fac NaN*ones(size(fac,1),length(vert)-3); 1:length(vert)];
pointe = patch('Faces',fac,'Vertices',vert,'FaceColor',col, 'FaceAlpha', 1, 'EdgeColor', col);
hold on

aux = [Xi; Xf-loF*(Xf-Xi)/norm(Xf-Xi)];
ligne = plot3(aux(:,1),aux(:,2),aux(:,3),'Color',col,'LineWidth',lw);

axis equal
