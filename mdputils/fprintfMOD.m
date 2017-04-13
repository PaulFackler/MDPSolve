% fixes the problem of printing sparse inputs with fprintf
function varargout = fprintf(varargin)
% adjust input arguments here
for k=1:nargin
    if issparse(varargin{k})
      try
        varargin{k} = full(varargin{k});
      catch
        error('cannot convert input to full for display')
      end
    end
end
% call the built-in function
[varargout{1:nargout}] = builtin('fprintf',varargin{:});
end
