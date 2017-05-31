function M = vectorToMatrix(y, nr, nc, shift)
% VECTORTOMATRIX custom routine to set a vector into a matrix with a given
% shift.
%
% VECTORTOMATRIX(Y, M, NR, NC, SHIFT) puts the vector Y(SHIFT+1:end) in the
% NR x NC matrix M.
%
% See also MATRIXTOVECTOR
%
% BLB 2015

M = eye(nr, nc);
for i = 1 : nr
    for j = 1 : nc
        m = nc*(i-1) + j;
        M(i,j) = y(shift+m);
    end
end

end

