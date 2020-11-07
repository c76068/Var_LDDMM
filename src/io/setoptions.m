function v = setoptions(struct, fieldname, default_value, possible_value)
% setoptions(struct, fieldname, default_value, possible_value) - given a struct,
% check if the field 'fieldname' exists. If not it create the field with 
% the default value. If it exists and possible_value is 
% given, it checks that the existing value is in possible value.
%
% Input
%  struct : a structure
%  fieldname : a string with the fieldname
%  default_value : a default value 
% Output
%  struct : the checked structure.
%

v = struct;

if ~isfield(struct, fieldname) && (nargin ==2)
    
    error([fieldname,' is mandatory!'])
    
elseif ~isfield(struct, fieldname) && (nargin > 2) 
    
    v.(fieldname) = default_value;

elseif nargin == 4 && ischar(default_value)
    
    if sum(strcmp(possible_value(:),v.(fieldname))) == 0
        error(['Possible values for ',fieldname,' is ',cell2mat(possible_value)])
    end
    
end
        