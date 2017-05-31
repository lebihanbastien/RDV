function f = LRT(t,X,A_v)

A_LRT = reshape(A_v,6,6);

f = A_LRT*X;


end