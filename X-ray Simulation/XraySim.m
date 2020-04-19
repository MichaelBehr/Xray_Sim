function [Xray_Bitmap] = XraySim(STLFILE,varargin)
% XraySim: Computes the simulated X-ray of a 3D object saved as a STL file. 
% Uses VOXELISE.m and projection2D.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AUTHOR        Michael G. Behr
% CONTACT       MichaelBehr13@gmail.com
% INSTITUTION   Toronto Rehabilitation Institute, University Health Network
%
%
% USAGE        [Xray_Bitmap] = XraySim('STLFILE.stl')
%        or... [Xray_Bitmap] = XraySim('STLFILE.stl',gridx,gridy,gridz)
%        or... [Xray_Bitmap] = XraySim('STLFILE.stl',gridx,gridy,gridz,'FILENAME.bmp')
%        or... [Xray_Bitmap] = XraySim('STLFILE.stl',gridx,gridy,gridz,'FILENAME.bmp','plot',true)
%        or... [Xray_Bitmap] = XraySim('STLFILE.stl',gridx,gridy,gridz,'FILENAME.bmp','plot',true,'D_Manual',Distance)
%
% INPUTS
%
%     STLFILE - Mandatory  - string       - Filename of the STL file.
%
%     gridx   - optional - 1xP array      - List of the grid X coordinates. 
%                           OR an integer - Number of voxels in the grid in 
%                                           X direction. Default is 100.
%                                           Maximum value is 209.
%
%     gridy   - optional - 1xQ array      - List of the grid Y coordinates.
%                           OR an integer - Number of voxels in the grid in 
%                                           Y direction. Default is 100.
%                                           Maximum value is 209.
%
%     gridz   - optional - 1xR array      - List of the grid Z coordinates.
%                           OR an integer - Number of voxels in the grid in 
%                                           Z direction. Default is 100.
%                                           Maximum value is 209.
%     
%     'FILENAME.bmp' - optional - string  - Name of output bitmap simulated
%                                           x-ray file
%     NAME - PAIR INPUT ARGS
%
%     'plot'  - optional - string         - Flag to output the 3D mesh
%                                           plot. Default is false.
% 
%     true    - optional - boolean        - Pairs with plot. Set to either 
%                                           true or false depending on 
%                                           whether you want 3D mesh plot
%
%     'D_Manual' - optional - string      - Flag to set object distance 
%                                           manually
%
%     Distance   - optional - Single or   - Value of the object distance
%                             Double        assumed to be positive                               
%
% OUTPUT
%      
%     Xray_Bitmap - Mandatory - 2D array of type single - 2D X-ray 
%                                                         projection
% 
% EXAMPLES
%     
%     Only the STL file input:
%     >>  [Xray_Bitmap] = XraySim('STLFILE.stl')
% 
%     Defining the voxelization grid manually
%     >>  [Xray_Bitmap] = XraySim('STLFILE.stl',gridx,gridy,gridz)
% 
%     File output is called Xray.bmp and saved in current folder
%     >>  [Xray_Bitmap] = XraySim('STLFILE.stl',gridx,gridy,gridz,'Xray.bmp')
%
%     Display a 3D mesh plot of the object
%     >>  [Xray_Bitmap] = XraySim('STLFILE.stl',gridx,gridy,gridz,'Xray.bmp','plot',true)
%
%     Change the distance of the object manually
%     >>  [Xray_Bitmap] = XraySim('STLFILE.stl',gridx,gridy,gridz,'Xray.bmp','plot',true,'D_Manual',500)
% IMPORTANT NOTES
%
%   - The ReadMe.txt file contains a detailed breakdown of the thought
%   process behind the overall coding implementation
%   - STLFILE can be both binary or ascii format
%   - The mesh (STLFILE) must be properly closed (watertight)in order
%   for proper voxelization into a 3D image matrix by use of VOXELISE.m
%   - It is recommended to at least use the default amount of grid size
%   when doing voxelization as you do not want to leave out 3D detail
%   - Projection parameters are set automatically. There are more details 
%   in the ReadMe.txt file about the assumptions
%   - 2D x-ray projections are computed by use of project2D.m, adapted from 
%   'projection.m' in the attached appendix code folder
%   - the distance of the x-ray source to the 3D object is automatically 
%   set based on cone beam geometry, and 'real' object size (matrix). This
%   is in order to make sure the object fits inside the conical beam. For 
%   example if the 'real' size of the object is increased in all 3 
%   cartesian directions the projection will stay nearly identical, as the
%   object will be placed further from the source in the simulation (to 
%   ensure fit).
%   -Warnings are given if the object is too large to fit inside the cone
%   beam geometry; it is proportionally scaled down to compensate
%
% REFERENCES
%
%   - Mesh voxelisation toolbox by: Adam A.
%     Specifically the two functions: READ_stl.m and VOXELISE.m which read
%     in STL files and perform mesh voxelization using xyz raycasting
%     respectively.
%     https://www.mathworks.com/matlabcentral/fileexchange/27390-mesh-voxelisation
%
%   - 3D Cone beam CT (CBCT) projection backprojection FDK, iterative 
%     reconstruction toolbox by: Kyungsang Kim
%     projection2D.m is adapted using code from projection.m
%     XraySim.m has code adapted from ParamSetting.m
%     https://www.mathworks.com/matlabcentral/fileexchange/35548-3d-cone-beam-ct-cbct-projection-backprojection-fdk-iterative-reconstruction-matlab-examples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Read STL file into MATLAB for Voxelization                         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read STL object into matlab to display 3D mesh
[Vertices] = READ_stl( STLFILE );

x = squeeze( Vertices(:,1,:) )';
y = squeeze( Vertices(:,2,:) )';
z = squeeze( Vertices(:,3,:) )';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Check Input Parameters                                             %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Manual distance flag set to false
Flag = false;

% Check for grid parameters + violations of max grid size 209x209x209
if nargin>1
    if ((varargin{1,1}<210)||(varargin{1,2}<210)||(varargin{1,2}<210))
        gridx = varargin{1,1};
        gridy = varargin{1,2};
        gridz = varargin{1,3};
    else
        warning('Dimensions too large. Reducing size to 200x200x200')
        gridx = 200;
        gridy = 200;
        gridz = 200; 
    end
else
    gridx = 100;
    gridy = 100;
    gridz = 100;
end

% Check for Image filename input
if nargin > 4
    filename = varargin{1,4};
else
    filename = 'Simulated X-ray Image.bmp';
end

% Check for mesh plot parameter
if (nargin>5)
    if (strcmp(varargin{1,5},'plot')) && (varargin{1,6}==true)
        figure;patch(x,y,z,'b');
    else
        
    end
else
end

% Check for manual distance parameter
if (nargin>7)
    if (strcmp(varargin{1,7},'D_Manual')) && (isa((varargin{1,8}),'double')||isa((varargin{1,8}),'single'))
        param.DSO = varargin{1,8};
        Flag = true;
    else
        
    end
else
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Use VOXELISE.m to voxelize the 3D mesh into a 3D image array       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Use the voxelise function to obtain a 3D matrix representation of the
% object using xyz raycasting projections
Object = single(VOXELISE( gridx, gridy, gridz, STLFILE,'xyz'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. Set Parameters - adapted from ParamSetting.m                       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Organize parameters into a structure variable

% 3D grid parameters
param.nx = gridx; 
param.ny = gridy;
param.nz = gridz;

% Object real size.
param.sx = round(gridx); % mm 
param.sy = round(gridy); % mm
param.sz = round(gridz); % mm

% The detector panel pixel density (number of pixels).
param.nu = 256;     % pixels
param.nv = 256;

% Detector real size, close to that of an actual flat panel detector.
% Works out to almost 300 PPI^2 with these settings.
param.su = 375;     % mm (real size)
param.sv = 375;     % mm

% Distance of X-ray point source to detector
param.DSD = 1000;   % mm (real distance)

% if no manual distance input
if Flag == false
    % Distance of X-ray point source to object origin - calculated by trig 
    % geometry to fit the conical ray around the xy matrix square. This 
    % guarantees projection is not cutoff and encapsulates full 3d matrix

    % tan(theta) = opp/adj - determines angle of maximum x-ray that fits panel
    Ray_angle = (param.su/2) / param.DSD;

    % Find minimum radius of circle that encompasses XY of 3D array
    Circle_radius = sqrt(( param.sx/2 )^2 + (param.sy/2)^2 );

    % Use radius to compute distance to object (adjusted by half of the
    % z-dimension)
    param.DSO = round( (Circle_radius/Ray_angle) + param.sz/2 );	

end

% Single voxel size
param.dx = param.sx/param.nx; 
param.dy = param.sy/param.ny;
param.dz = param.sz/param.nz;
param.du = param.su/param.nu;
param.dv = param.sv/param.nv;

% Calculate geometric sampling vectors for the object
param.xs = (-(param.nx-1)/2:1:(param.nx-1)/2)*param.dx;
param.ys = (-(param.ny-1)/2:1:(param.ny-1)/2)*param.dy;
param.zs = (-(param.nz-1)/2:1:(param.nz-1)/2)*param.dz;

% Calculate geometric sampling vectors for the detector
param.us = (-(param.nu-1)/2:1:(param.nu-1)/2)*param.du;
param.vs = (-(param.nv-1)/2:1:(param.nv-1)/2)*param.dv;

% Type of interpolation to use when computing 2D projection
param.interp = 'linear';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5. Compute the 2D projection using Conical X-rays                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Xray_Bitmap] = mat2gray(projection2D(Object,param));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 6. Save X-ray projection as 8 bit Bitmap file                         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Write out image
imwrite(Xray_Bitmap, filename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END FUNCTION: XraySim.m                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

