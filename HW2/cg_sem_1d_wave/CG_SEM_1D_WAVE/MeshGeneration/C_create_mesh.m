function [Region] = C_create_mesh(Dati)
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

%================================================
% GEOMETRICAL INFO
 nEl = 2^Dati.nRefinement; 
 nVert = nEl + 1;
 p = linspace(x0,xL,nVert);
 t = [[1:nVert-1]' [2:nVert]']';
 MeshSize = (xL-x0)./nEl;
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
           
           

           noe=nov(npdx,ne);

xy=zeros(noe,1); ww=zeros(noe,1); jacx=zeros(ne,1);
xx=zeros(2,ne);
H=(xb-xa)/ne;
for ie=1:ne
xb_ie=xa+ie*H;
xa_ie=xb_ie-H;
xx(1:2,ie)=[xa_ie;xb_ie];
jacx(ie)=.5*(xb_ie-xa_ie);
for i=1:npdx
xy(nov(i,ie))=x(i)*jacx(ie)+.5*(xb_ie+xa_ie);
ww(nov(i,ie))=ww(nov(i,ie))+wx(i)*jacx(ie);
end
end