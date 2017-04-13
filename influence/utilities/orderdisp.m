% computes a string that summarizes the order of operations
% If no output is requested it displays this string.
% Input is the order information returned by getvarorder

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

function s=orderdisp(V,names,tex)
m=size(V,1)+1;
if nargin<3 || isempty(tex)
  tex=false;
end
if nargin<2 || isempty(names)
  names=cell(1,m);
  if tex
    for i=1:m, names{i}=['F^{' num2str(i) '}']; end
  else
    for i=1:m, names{i}=['F' num2str(i)]; end
  end
end
if tex, open='_{'; else open='{'; end
F=cell(m,1);
for i=1:m-1
  o1=V{i,1}; 
  o2=V{i,2};
  if max([V{:,3:4}])>9, form='%1i,'; else form='%1i'; end
  if isempty(F{o1})
    if length(V{i,3})>1
      s1=[names{o1} open sprintf(form,sort(V{i,3}(1:end-1))) num2str(V{i,3}(end)) '}'];
    else
      s1=[names{o1} open sprintf('%1i',V{i,3})  '}'];
    end
  else
    s1= F{o1};
  end
  if isempty(F{o2})
    if length(V{i,4})>1
      s2=[names{o2} open sprintf(form,sort(V{i,4}(1:end-1))) num2str(V{i,4}(end)) '}'];
    else
      s2=[names{o2} open sprintf('%1i',V{i,4})  '}'];
    end
  else
    s2= F{o2};
  end
  outvars=unique([V{i,3} V{i,4}]);
  if form(end)==',', sd=','; else sd=''; end
  if tex
    ss='';
    for j=1:length(outvars)
      if ismember(outvars(j),V{i,7}), ss=[ss '\underline{' num2str(outvars(j)) '}' sd ];
      else                            ss=[ss  num2str(outvars(j)) sd];
      end
    end
    if ss(end)==',', ss=ss(1:end-1); end
  end
  if isempty(find([V{:,1}]==o2,1)) || find([V{:,1}]==o1,1)<find([V{:,1}]==o2,1) || 1
    if tex
      if length(outvars)>1
        F{o1}=['\underbrace{\left('  s1  s2 '\right)}_{' ss '}'];
      else
        F{o1}=['\underbrace{\left('  s1  s2 '\right)}_{' num2str(outvars) '}'];
      end
      
    else
      F{o1}=['('  s1 ',' s2 ')'];
    end
  else
    if tex
      F{o1}=['\left('  s2  s1 '\right)'];
    else
      F{o1}=['('  s2 ',' s1 ')'];
    end
  end
end
s=F{o1};
if nargout<1 
  if tex
    figure; clf; 
    h=text(0.5,0.5,['$$' s '$$'],'Interpreter','LaTex','HorizontalAlignment','center');
    set(h,'FontSize',2*get(h,'FontSize'),'FontWeight','bold')
    set(gca,'Xcolor',[1 1 1],'YColor',[1 1 1])
  else
    disp(s); 
  end
else
  if tex
    s=['\scriptsize\[\hspace*{-.75in}' s '\]\normalsize'];
  end
end
return
