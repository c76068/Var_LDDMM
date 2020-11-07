function []=export_var_vtk(pos,vec,name,encod_type)
% export_mom_vtk(pos,mom,fname) create a vtk file 'fname' containing the initial momentums.
% Vectors should be in matrix mom (nb_points x dim_ambient_space) and position in  pos (matrix nb_points x dim_ambient_space). 
%
% encod_type is deprecated.
%
% Note: If the file fname already exists,  it is OVERWRITTEN without any warning.
%
% Input :
%   pos: matrix nb_points x dim_ambient_space
%   mom:  matrix nb_points x dim_ambient_space
%   name : string typically ('/path/to/myfile.vtk')
%
% See also : export_fshape_vtk, export_fshape_ply, export_atlas_HT, export_atlas_free, jnfmean_tan_free
%


if nargin ==3
    encod_type = 'ascii';
end

if size(pos,2)<=2
    pos = [pos,zeros(size(pos,1),3-size(pos,2))];
end

for k=1:length(vec)
    if size(vec{k},2)<=2
        vec{k} = [vec{k},zeros(size(vec{k},1),3-size(vec{k},2))];
    end
end

nb_points = size(pos,1);

fid = fopen(name, 'w'); 


%-------------
%  header
%-------------

%ASCII file header
fprintf(fid, '# vtk DataFile Version 3.0\n');
fprintf(fid, 'VTK from fshapesTk\n');
if strcmp(encod_type,'ascii')
    fprintf(fid, 'ASCII\n\n');
else
    fprintf(fid, 'BINARY\n\n');
end

%-------------
%  Position
%-------------

%ASCII sub header
fprintf(fid, 'DATASET STRUCTURED_GRID\n');
fprintf(fid, ['DIMENSIONS ',num2str(nb_points),' 1 1\n']);
%Record position
fprintf(fid, ['POINTS ',num2str(nb_points),' float\n']); 
if strcmp(encod_type,'ascii')
    fprintf(fid,'%G %G %G\n',pos');
else
    fwrite(fid,pos','float','b');
end

%-------------
%  vectors
%-------------

%ASCII sub header
fprintf(fid, ['\nPOINT_DATA ',num2str(nb_points),'\n']);
%Record vectors
for k=1:length(vec)
    
    fprintf(fid,['FRAME VECTOR momentum',num2str(k),'float\n']);
    if strcmp(encod_type,'ascii')
        fprintf(fid,'%G %G %G\n',vec{k}');
    else
        fwrite(fid, vec{k}','float','b');
    end
    
end

fclose(fid);

fprintf('\nFile saved in %s',name)

end
