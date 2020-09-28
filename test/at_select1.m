%function sel=at_select(atnam,at_res,reslst,atlst,atdim,offs),
%---------------------------------------------------
%  df-aug-03    df-dec-97
%	select atoms according to atnam & reslst
%   standard lists: 'all', 'heavy', 'bb', 'nh'
%   atdim is required if atlist in not a standard one
%---------------------------------------------------
%if nargin<6, offs=2; end		%default
%if nargin<5, atdim=2; end		%default
atlst = ' HN'
offs = 1;
atdim = 2;
atsel_flag=1;
if strcmp(atlst,'heavy')|strcmp(atlst,'HEAVY')
    at='NCOS'; atlen=1;
elseif strcmp(atlst,'bb')|strcmp(atlst,'BB')
    at='N CAC O '; atlen=2;
elseif strcmp(atlst,'nh')|strcmp(atlst,'NH') 
   %not now 
   at='HN'; atlen=2;
elseif strcmp(atlst,'all')|strcmp(atlst,'ALL')|isempty(atlst), 
   %get all atoms
   %at='HNCOS';atlen=1; 
   atsel_flag=0; atlen=1;
else
    atlen=atdim;
    at=atlst;
end

atnampos =[7+offs:7+offs+atlen-1];
atletter = atnampos + 6;
nat=size(atnam,1);
indsel=NaN*ones(nat,1);
for ii=1:nat,
  if isempty(reslst),			%take all resid
    if atsel_flag==0,           	%all atoms
        indsel(ii)=ii;
    else                        	%select atoms
        if ~isempty(findstr(atnam(ii,atnampos),at)),
          if ~isempty(findstr(atnam(ii,atletter),'A'))
            indsel(ii)=ii; 
          end  
        end
    end
  else					        %res. selection
   if ~isempty(find(at_res(ii,2)==reslst)),
     if atsel_flag==0,			%all atoms
        indsel(ii)=ii; 
     else				      %select atoms
        if ~isempty(findstr(atnam(ii,atnampos),at)),
          indsel(ii)=ii;
        end
     end
   end
  end
end
sel=indsel(~isnan(indsel));
%return
%===================================================

