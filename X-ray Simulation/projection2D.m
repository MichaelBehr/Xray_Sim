function [proj2d] = projection2D(data3d,param)
% projection2D: Computes the 2D projection of an 3D array
% Uses a perspective projection implementation, in order to simulate a
% conical x-ray from a point source. Utilizes 2D interpolation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AUTHOR        Michael G. Behr
% CONTACT       MichaelBehr13@gmail.com
% INSTITUTION   Toronto Rehabilitation Institute, University Health Network
%
%
% USAGE         [proj2d] = XraySim(data3d,param)
%
% INPUTS
%
%               data3d - 3 dimensional array storing the 3D object
%               param  - the imaging parameters necessary for perspective 
%               projection
%
% OUTPUT     
%
%               proj2d - 2 dimensional array storing the 2D projection in
%               the Z direction (+ to -)
%
% EXAMPLE
%               Create the 2D perspective projection of the 3D array,
%               OBJECT
%               >> [PROJ] = projection2D(OBJECT, PARAMETERS);
%
%
% IMPORTANT NOTES
%
%   - Implementation adapted from projection.m to project from z-direction
%   instead of y-direction
%   - The ReadMe.txt file contains details into the math behind the
%   projection
%   - The projections are computed using 2D linear interpolation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create a 2D array to store the projections
proj2d = (zeros(param.nu,param.nv,'single'));

% Creates mesh sampling grids based on the geometric sampling vectors
% of the detector plate and the 3D object
[uu,vv] = meshgrid(param.us,param.vs);
[~,yy] = meshgrid(param.xs,param.ys);
[xx,~] = meshgrid(param.xs,param.zs);
data3d(isnan(data3d))=0;

% Go slice by slice in the Z-direction, computing each projection utilizing
% 2D, linearly interpolated mesh re-sampling grids
for iz = 1:param.nz
    
    % Distance ratio to determine projection size change (geometry). It
    % changes based on which z-slice you are currently in.
    Ratio = (param.zs(iz)+param.DSO)/(param.DSD);
    
    % Scale the mesh re-sampling grids by the size change ratio
    pu = uu*Ratio;
    pv = vv*Ratio;    
    
    % Shift the points back towards the origin
    pu = (pu - xx(1,1))/(param.dx)+1; 
    pv = (pv - yy(1,1))/(param.dy)+1; 
    
    
    % 2D interpolation, with the detector grids. The scaled and shifted 
    % mesh grids pu and pv are interpolated with the current z-slice
    tmp = (interp2((single(data3d(:,:,iz))),(single(pv)),(single(pu)),param.interp));
    
    tmp(isnan(tmp))=0;
    proj2d = proj2d + tmp';
end

% End FUNCTION: projection2D.m
end





