
function [vertex,mom] = import_mom_vtk(filename)

% read_vtk - read data from VTK file. (Based on a work Mario Richtsfeld)
%
%   [vertex,face] = read_vtk(filename, verbose);
%
%   'vertex' is a 'nb.vert x 3' array specifying the position of the vertices.
%   'face' is a 'nb.face x 2' (POLYGONS) or 'nb.face x 2' (LINES) array specifying the connectivity of the mesh.
%
% Authors : this file is part of the fshapesTk by B. Charlier, N. Charon, A. Trouve (2012-2014)  


fid = fopen(filename,'r');

%---------------
% read header 
%----------------

if( fid==-1 )
    error('Can''t open the file.');
end

str = fgets(fid);  % -1 if eof
if ~strcmp(str(3:5), 'vtk')
    error('The file is not a valid VTK one.');    
end

%jump 3 lines
[~] = fgets(fid);
[~] = fgets(fid);
[~] = fgets(fid);

%----------------
% read vertices
%----------------
   [~] = fgets(fid);[~] = fgets(fid);
str = fgets(fid);
info = sscanf(str,'%s %*s %*s', 6);

if strcmp(info,'POINTS')
    nvert = sscanf(str,'%*s %d %*s', 1);
    [A,cnt] = fscanf(fid,'%f ', 3*nvert);
    if cnt~=3*nvert
        warning('Problem in reading vertices.');
    end
    A = reshape(A, 3, cnt/3);
    vertex = A';
end

%----------------
% read polygons
%----------------

str = fgets(fid);
temp = fgets(fid);
info = sscanf(temp,'%s %*s %*s');

if strcmp(info,'VECTORS') 
     nvect = sscanf(str,'%*s %d %*s', 1);
    [A,cnt] = fscanf(fid,'%f ', 3*nvect);
    if cnt~=3*nvect
        warning('Problem in reading momentums.');
    end
    A = reshape(A, 3, cnt/3);
    mom = A';
    
else
    error('Problem in reading momentums.')
end


fclose(fid);

end
