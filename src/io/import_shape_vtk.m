function fshape = import_shape_vtk(filename,d)
% IMPORT_VTK(filename) import a vtk file as a fshape structure. (This function uses a function written by M. Richtsfeld).
%
%Input :
%   path : string path to file
%Output :
%   fshape : struct('x',...,'G',...,'f',...)
%
% See also : import_shape_ply, import_shape_obj
%

if nargin == 1
    d=3;
end

[fshape.x,fshape.G,fshape.f] = read_vtk(filename,d);


end

function [vertex,face,signal] = read_vtk(filename,d)

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

str = fgets(fid);
info = sscanf(str,'%s %*s %*s', 6);

if strcmp(info,'POINTS') 
    nvert = sscanf(str,'%*s %d %*s', 1);
    [A,cnt] = fscanf(fid,'%G ', d*nvert);
    if cnt~=d*nvert
        warning('Problem in reading vertices.');
    end
    A = reshape(A, d, cnt/d);
    vertex = A';
end

%----------------
% read polygons
%----------------

str = fgets(fid);
info = sscanf(str,'%s %*s %*s');

if strcmp(info,'POLYGONS')  || strcmp(info,'LINES')|| strcmp(info,'VERTICES')
    nface = sscanf(str,'%*s %d %d', 2);
    [A,cnt] = fscanf(fid,'%d ', nface(2));
    if cnt~=nface(2)
        warning('Problem in reading faces.');
    end
    A = reshape(A,nface(2)/nface(1),nface(1));
    face = (A(2:end,:)+1)';
else
    error('Problem in reading faces.')
end

%-------------------
% read signal
%-------------------

try
    str = fgets(fid);
    info = sscanf(str,'%s %*s', 10);
    
    if strcmp(info,'POINT_DATA')
        
        nvertex = sscanf(str,'%*s %d', 1);
        [~] = fgets(fid);
        [~] = fgets(fid);
        
        [signal,cnt] = fscanf(fid,'%G', nvertex);
        if cnt~=nvertex
            error('Problem in reading signal.');
        end
    elseif strcmp(info,'CELL_DATA')
        nface = sscanf(str,'%*s %d', 1);
        [~] = fgets(fid);
        [~] = fgets(fid);
        
        [signal,cnt] = fscanf(fid,'%G', nface);
        if cnt~=nface
            error('Problem in reading signal.');
        end
    else
        error('Problem in reading signal.');
    end
catch
    signal = zeros( [size(vertex,1) ,1] );
end

    
fclose(fid);

end
