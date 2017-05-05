function varargout = subsref(a,S)
switch S.type
  case '()'
    if numel(S.subs)==2
      indr=S.subs{1};
      indc=S.subs{2};
      if strcmp(indc,':')
        B=[a.data{:}];
        B=B(indr,:);
      elseif isempty(indc)
        B=zeros(a.m,0);
        B=B(indr,:);
      else
        B=cell(1,a.n);
        for k=1:a.n
          if indc(end)<a.s(k), break; end
          B{k}=a.data{k}(indr,ismember(a.s(k)+1:a.s(k+1),indc));
        end
        B=[B{:}];
      end
    else
      error('Single indexing of cellsparse objects is not supported')
    end
    varargout{1}=B;
  case '{}'
    varargout=a.data(S.subs{:});
  case '.'
    error('??? Attempt to reference field of non-structure array.')
  otherwise
    error('unknown error')
end