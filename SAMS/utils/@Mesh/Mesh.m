classdef Mesh < handle
%MESH Summary of this function goes here
%   Detailed explanation goes here

properties
    F           %matrix of faces
    nF          %number of faces
    Nf          %normal vectors of faces
    V           %matrix of vertices
    nV          %number of vertices
    Nv          %vertex normals, derived from face normals
    F2V         %Map between faces to vertices
    V2E         %Map between vertices to edges
    A           %Adjacency matrix
    E           %Edge list, not used often
    nE          %number of edges
    BE          %Boundary edges, not often saved
    BV          %Boundary vertices, not often saved
    Aux         %Auxillary variables, mostly stored curvature information
    V2V         %Deprecated
    E2E         %Deprecated
end

methods
    function obj = Mesh(varargin)
        %% Constructor if you already put in a mesh, just repeat input
        if (length(varargin)==1) && isa(varargin{1},'Mesh')
            obj.V=varargin{1}.V;
            obj.nV=varargin{1}.nV;
            obj.Nv=varargin{1}.Nv;
            obj.F=varargin{1}.F;
            obj.nF=varargin{1}.nF;
            obj.Nf=varargin{1}.Nf;
            obj.F2V=varargin{1}.F2V;
            obj.V2E=varargin{1}.V2E;
            obj.A=varargin{1}.A;
            obj.E=varargin{1}.E;
            obj.nE=varargin{1}.nE;
            obj.BE=varargin{1}.BE;
            obj.BV=varargin{1}.BV;
            obj.Aux=varargin{1}.Aux;
            obj.V2V=varargin{1}.V2V;
            obj.E2E=varargin{1}.E2E;
        elseif length(varargin)>=2
            %% Else read from file
            switch(varargin{1})
                case 'ply'
                    [obj.V,obj.F] = obj.read_ply(varargin{2});
                    if size(obj.V,1) ~= 3
                        obj.V = obj.V';
                    end
                    if size(obj.F,1) ~= 3
                        obj.F = obj.F';
                    end
                case 'off'
                    [obj.V,obj.F] = obj.read_off(varargin{2});
                case 'obj'
                    [obj.V,obj.F] = obj.read_obj(varargin{2});
                case 'VF'
                    if size(varargin{2},1)==2
                        obj.V=[varargin{2};zeros(1,size(varargin{2},2))];
                    else
                        obj.V=varargin{2};
                    end
                    obj.F=varargin{3};
            end
            %% Compute some remaining information
            obj.F2V = obj.ComputeF2V;
            obj.V2E = obj.ComputeV2E;
            obj.nV = size(obj.V,2);
            obj.nF = size(obj.F,2);
            obj.nE = size(obj.V2E,2);
        else
            %% Empty mesh constructor
            obj.F=[];
            obj.V=[];
        end
    end
end
methods(Static)
    [V,F,Fs] = read_off(filename)
    [V,F,Fs] = read_obj(filename)
    [V,F,Fs] = read_ply(filename)
    v = getoptions(options, name, v, mendatory)
    [D,S,Q] = PerformFrontPropagation(vertex, faces, W,start_points,end_points, nb_iter_max, H, L, values, dmax)
end
end


