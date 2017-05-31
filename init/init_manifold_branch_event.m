function [ manifold_branch ] = init_manifold_branch_event(stability, way, varargin)
% INIT_MANIFOLD_BRANCH_EVENT is a DEPRECATED version of the routine
% INIT_MANIFOLD_BRANCH. It has been kept for backwards compatibility 
% but does exactly the same thing as INIT_MANIFOLD_BRANCH.
%
% See also INIT_MANIFOLD_BRANCH
%
% BLB 2016

    manifold_branch = init_manifold_branch(stability, way, varargin{:});
end




