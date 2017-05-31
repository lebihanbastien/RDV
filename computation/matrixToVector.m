function y = matrixToVector(y, M, nr, nc, shift)
% MATRIXTOVECTOR custom routine to set a matrix into a vector with a given
% shift.
%
% MATRIXTOVECTOR(Y, M, NR, NC, SHIFT) puts the NR x NC matrix M in the
% vector Y(SHIFT+1:end) 
%
% See also VECTORTOMATRIX
%
% BLB 2015

for i = 1 : nr
    for j = 1 : nc
        m = nc*(i-1) + j;
        y(shift+m) = M(i,j);
    end
end

end


