function XYZrot = rotation(XYZ,P)
% Rotation function on tensor XYZ, with a rotation matrix P
%
% Emilien Fabacher 2015
XYZrot = XYZ;
for i =1:size(XYZ,1)
    XYZrot(i,:) = (P*XYZrot(i,:)')';
end