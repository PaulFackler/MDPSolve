function printD(DD,print,dname,typestring)
d=length(DD.names);
if nargin<3, dname      = 'D';        end    % name given to diagram
if nargin<4, typestring = 'sadfurcp'; end    % valid letters to designate types
parents=getparents(DD);

  if print>1
    % print dummy variables for CPDs
    % so code will run even though cpds are not defined
    for i=1:d      
      if ~isempty(DD.cpdnames{i})
        fprintf([DD.cpdnames{i} '=[];\n'])
      end
    end
  end
  k=0;
  fprintf([dname '=[];\n']);
  for i=1:d      
    if ischar(DD.types{i}), typei = DD.types{i};
    else                    typei = typestring(DD.types{i});
    end
    fprintf([dname '=add2diagram(' dname ',''' DD.names{i} ...
        ''',''' typei ''','])  
    if DD.obs(i), fprintf('true,{')
    else fprintf('false,{')
    end
    for j=1:length(parents{i})
      if j>1, fprintf(','); end
      fprintf(['''' DD.names{parents{i}(j)} ''''])
    end
    if isempty(DD.cpdnames{i})
      fprintf('},[],')
    else
      fprintf(['},' DD.cpdnames{i} ','])
    end
    if isempty(DD.locs)
      fprintf('[]');
    else
      fprintf('[%1.4f,%1.4f]',DD.locs(i,1),DD.locs(i,2));
    end

    if isempty(DD.attachments)
      fprintf(');\n')
    else
      fprintf(',[')
      ii=find(DD.attachments(:,2)==i);
      for j=1:length(parents{i})
        k=find(DD.attachments(ii,1)==parents{i}(j)); 
        if j==1
          fprintf('%1i %1i',  DD.attachments(ii(k),3), DD.attachments(ii(k),4))
        else
          fprintf(';%1i %1i', DD.attachments(ii(k),3), DD.attachments(ii(k),4))
        end
      end
      fprintf(']);\n')
    end
  end
  
  % location and attachment info now in node line
  if 0
    disp([dname '.locs=[ ...'])
    for j=1:2
      for i=1:size(DD.locs,1)
        if i>1, fprintf(' '); end
        fprintf('%5.3f',DD.locs(i,j));
      end
      if j<2, fprintf(';\n')
      else    fprintf(']'';\n')
      end
    end
    if ~isempty(DD.attachments)
      disp([dname '.attachments=[ ...'])
      for j=1:4
        for i=1:size(DD.attachments,1)
          if i>1, fprintf(' '); end
          fprintf('%2i',DD.attachments(i,j));
        end
        if j<4, fprintf(';\n')
        else    fprintf(']'';\n')
        end
      end
    end
  end
