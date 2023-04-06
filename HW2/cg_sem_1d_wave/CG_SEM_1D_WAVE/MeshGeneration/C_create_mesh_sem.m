function [Region] = C_create_mesh_sem(Dati)
%% [Region] = C_create_mesh(Dati)
%==========================================================================
% Creates triangular or quadrilateral mesh
%==========================================================================
%    called in C_main1D.m
%
%    INPUT:
%          Dati        : (struct)  see C_dati.m
%
%    OUTPUT:
%          Region      : (struct) having fields: dimension
%                                                domain
%                                                mesh size
%                                                number of vertices
%                                                number of elements
%                                                coordinates
%                                                boundary points
%                                                connectivity
%

x0 = Dati.domain(1);
xL = Dati.domain(2);


if (strcmp(Dati.fem,'P1'))
    npdx = 2;
elseif (strcmp(Dati.fem,'P2'))
    npdx = 3;
elseif (strcmp(Dati.fem,'P3'))
    npdx = 4;
elseif (strcmp(Dati.fem,'P4'))
    npdx = 5;
elseif (strcmp(Dati.fem,'P5'))
    npdx = 6;
else
    disp('case not implemented')
end
    

%================================================
% GEOMETRICAL INFO
nEl = 2^Dati.nRefinement;
MeshSize = (xL-x0)./nEl;

xy = zeros(nEl,npdx);
xx = zeros(2,nEl);

i = 0; 
for ie = 1 : nEl
    xb_ie = x0 + ie*MeshSize;
    xa_ie = xb_ie - MeshSize;
    
    [xp,wp] = xwlgl(npdx,xa_ie,xb_ie); 
    p(i+1:npdx+i) = xp;
    i = i + npdx ;
end

p = unique(p);
nVert = size(p,2); 

if (strcmp(Dati.fem,'P1'))
    t = [[1:npdx-1:nVert-npdx+1]' [2:npdx-1:nVert-npdx+2]']';
elseif (strcmp(Dati.fem,'P2'))
    t = [[1:npdx-1:nVert-npdx+1]' [2:npdx-1:nVert-npdx+2]' [3:npdx-1:nVert-npdx+3]']';
elseif (strcmp(Dati.fem,'P3'))
    t = [[1:npdx-1:nVert-npdx+1]' [2:npdx-1:nVert-npdx+2]' [3:npdx-1:nVert-npdx+3]' [4:npdx-1:nVert-npdx+4]']';
elseif (strcmp(Dati.fem,'P4'))
    t = [[1:npdx-1:nVert-npdx+1]' [2:npdx-1:nVert-npdx+2]' [3:npdx-1:nVert-npdx+3]' [4:npdx-1:nVert-npdx+4]' [5:npdx-1:nVert-npdx+5]']';
elseif (strcmp(Dati.fem,'P5'))
    t = [[1:npdx-1:nVert-npdx+1]' [2:npdx-1:nVert-npdx+2]' [3:npdx-1:nVert-npdx+3]' [4:npdx-1:nVert-npdx+4]' [5:npdx-1:nVert-npdx+5]' [6:npdx-1:nVert-npdx+6]']';
else
    disp('case not implemented')
end

%================================================

% struttura dati della mesh
Region = struct('dim',1,...
    'domain',Dati.domain,...
    'h',MeshSize,...
    'nvert',nVert,...
    'ne',nEl,...
    'coord',p',...
    'boundary_points',[x0,xL],...
    'connectivity',t);




