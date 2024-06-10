classdef matvecH_aug
    %
    % matvecH_aug class
    %
    % A matvecH_aug object is used to represent a matrix [Hall 0] in matrix-vector
    % multiplication where the matrix is stacked by H_i.mat horizontally
    % and zero
    %
    % the matvecH has input(s):
    %   n - number observation 
    %   path - path where H matrices are located
    %   nbeta - number of added zero columns
    %
    % and is based on a structure with the following fields:
    %
    % Calling Syntax:
    % P = matvecH(n,path,nbeta)
    %
    % T.Cho, 10/29/2020
    % Modified by M. Sabate Landman 6/7/2024
    
    properties
        n
        nbeta
        H
        transpose
    end % properties

    methods

        function P = matvecH_aug(varargin) % constructor
            switch nargin
                case 3
                    P.transpose = false;
                    P.n = varargin{1};
                    P.nbeta = varargin{2};
                    P.H = varargin{3};
                otherwise
                    error('Incorrect number of input arguments')
            end % switch
        end % constructor

        function P = ctranspose(P) % Overload transpose
            P.transpose = not(P.transpose); % switches boolean transpose flag
        end % transpose

        function y = mtimes(A,x)
            %load('H_final.mat');
            [mA,nA] = size(A.H);

            if A.transpose % transpose
                [mx,nx] = size(x);
                if mx ~= mA
                    error('Invalid size of x')
                end
                y = [A.H'*x; zeros(A.nbeta,1)];
            else % no transpose
                y = A.H*x(1:end-A.nbeta);
            end 
        end % mtimes
        
        function varargout = size(A,dim)
            load('H_final.mat');
            [mA,nA] = size(H);
            d(1) = mA;
            d(2) = nA*A.n + A.nbeta;
            if nargout == 1 || nargout == 0
                if nargin >1
                    varargout{1} = d(dim);
                else
                    varargout{1} = d;
                end
            else
                varargout{1} = d(1);
                varargout{2} = d(2);
            end
        end % size
        
        function l = length(A)
            load('H_final.mat');
            [mA,nA] = size(H);
            if nA*A.n + A.nbeta > mA
                l = nA*A.n + A.nbeta;
            else
                l = mA; 
            end
        end % length
        
    end % methods
    
end % end classdef