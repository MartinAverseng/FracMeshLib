% element-wise assembly of hypersingular operator on a generalized mesh. 
function [A] = bemAssembly(M)


[J,Nf] = localToGlobal(M,0);
Nelt = M.nelt;
A = zeros(Nf);
A = CppAssembly(Nf,Nelt,J-1,M.vtx,M.elt-1,A); % -1 for 0-indexing in C++.
A = (A + A')/2;


end

