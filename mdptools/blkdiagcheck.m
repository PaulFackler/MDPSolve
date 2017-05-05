% blkdiagcheck Determines block diagonal form of a matrix
% USAGE
%   [blocksr,blocksc,AA]=blkdiagcheck(A);
% INPUT
%   A : m x n matrix
% OUTPUTS
%   blocksr : k-element cell array composed of vectors of integers on {1,...,m}
%               indicating the rows with non-zero values in each block
%               (columns with all 0s are assigned to the last block and
%                  the rows assigned to this block are the all 0 rows)
%   blocksc : k-element cell array composed of vectors of integers on {1,...,n}
%               indicating the columns with non-zero values in each block
%               Note: each column appears in a single block so
%                 length([blocksc{:}]) is equal to n
%   AA      : m x n block diagonal form of A
%
% k is the number of blocks in the block diagonal form (could be 1)
% Within blocks the columns retain the original order.
% All zero blocks are inserted at the end.
%
% To obtain the ith block use A(blocksr{i},blocksc{i})
% To obtain AA use A([blocksr{:}],[blocksc{:}])
%
% If A is block diagonal then AA will equal A.

% MDPSOLVE: MATLAB tools for solving Markov Decision Problems
% Copyright (c) 2014, Paul L. Fackler (paul_fackler@ncsu.edu)
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without  
% modification, are permitted provided that the following conditions are met:
% 
%    * Redistributions of source code must retain the above copyright notice, 
%        this list of conditions and the following disclaimer.
%    * Redistributions in binary form must reproduce the above copyright notice, 
%        this list of conditions and the following disclaimer in the 
%        documentation and/or other materials provided with the distribution.
%    * Neither the name of the North Carolina State University nor of Paul L. 
%        Fackler may be used to endorse or promote products derived from this 
%        software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
% 
% For more information, see the Open Source Initiative OSI site:
%   http://www.opensource.org/licenses/bsd-license.php

function [blocksr,blocksc,AA]=blkdiagcheck(A)
n=size(A,2);
blocksr{1}=A(:,1)~=0;
blocksc{1}=1;
% check if column is all 0s
if any(blocksr{1})
  zblock=0;
else
  zblock=1;
end
for j=2:n
  jcol=A(:,j)~=0;
  if any(jcol)
    jblock=0;
    for i=length(blocksr):-1:1
      if any(jcol(blocksr{i})) % check if any common elements
        jcol = jcol | blocksr{i};
        if jblock==0 % add j to a previous block
          blocksc{i}=[blocksc{i} j];
        else  % merge j's block with a previous block
          blocksc{i}=[blocksc{i} blocksc{jblock}];
          blocksr(jblock)=[];
          blocksc(jblock)=[];
          if zblock>jblock, zblock=zblock-1; end
        end
        blocksr{i} = jcol;
        jblock=i;
      end
    end
    if jblock==0 % create a new block
      blocksr{end+1}=A(:,j)~=0;  %#ok<*AGROW>
      blocksc{end+1}=j; 
    end
  else % handle all 0 blocks
    if zblock>0, blocksc{zblock}=[blocksc{zblock} j];
    else         blocksc{end+1}=j; blocksr{end+1}=[]; zblock=length(blocksc);
    end
  end
end

% move empty block to the end
if zblock>0 
  if zblock<length(blocksc)
    blocksc{end+1}=blocksc{zblock}; blocksr{end+1}=blocksr{zblock};
    blocksc(zblock)=[]; blocksr(zblock)=[];
  end
  blocksr{end}=~blocksr{1};
  for i=2:length(blocksr)-1
    blocksr{end}=blocksr{end} & ~blocksr{i};
  end
end

for i=1:length(blocksc)
  blocksc{i}=sort(blocksc{i});
  blocksr{i}=find(blocksr{i})';
end

% produce block diagonal matrix
if nargout>2
  AA=A([blocksr{:}],[blocksc{:}]);
end