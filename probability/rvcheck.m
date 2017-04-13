function rvcheck(rv,parents)
if ~isstruct(rv)
  error('rv must be an rv structure generated by rvdef')
end
if ~isempty(rv)
  if ~isempty(rv.values) 
    if length(rv.values)~=rv.size
      error('rv size and number of values does not match')
    end
  end
  if ~isempty(parents) && ~isempty(rv.weights)
    if prod(D.sizes(pind))~=size(rv.weights,2)
        error('Number of columns in CPT weights does not match the number of parent values')
    end
  end
  if isa(rv.parameters,'function_handle')
    if ~isempty(parents)      
      if nargin(rv.parameters)~=length(parents)
        error('Number of inputs for rv parameter function does not match the number of parents')
      end
    end
  end
end
