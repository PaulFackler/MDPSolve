% cellsparse creates a class to work with cell arrays of sparse matrices
% Often matrices need to be built up blocks of columns at a time
% Each submatrix has the same number of rows and, in principle, the whole matrix
% could be obtained using [A{:}]. Due to memory considerations, however, it may 
% be better to just keep the matrix stored in block form and operate on it
% without having to have the whole matrix in contiguous memory.
% To create a cellsparse object use
%   B=cellsparse(A);
% Note that 
%   A=cellsparse(A);
% can also be used. The individual elements of A can still be obtained using
% the {} operator. For the most part A{..} and B{..} should operate identically.
%
% Currently the methods that can be used woth cellsparse objects are
%   * (mtimes)
%   size
%   disp
%   indexed extraction (e.g. B(1:2,:), B(1,[1 4 6]), etc.)
% disp and extraction may fail if B is large.
classdef cellsparse
   properties
     m
     n
     s
     data
   end
   methods
     function a=cellsparse(A)
       if ~iscell(A)
         error('cellsparse objects are created with a cell array of matrices.')
       end
       if isempty(A), a.m=0; 
       else           a.m=size(A{1},1); 
       end
       a.n=length(A);
       a.s=zeros(1,a.n+1);
       for j=1:a.n
         if size(A{j},1)~=a.m
           error('All submatrices in cellsparse class must have the same number of rows.')
         end
         a.s(j+1)=a.s(j)+size(A{j},2);
       end
       a.data=A;
     end
   end
end