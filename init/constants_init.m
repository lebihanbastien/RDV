function cst = constants_init()
% CONSTANT_INIT defines various constants (coordinates, environement
% constants...). 
% 
% Should be called priori to any computation.
%       
% BLB 2016

%-------------------------------------------------------------------------%
% Deprecated (true and false are already defined in native MATLAB)
%-------------------------------------------------------------------------%
cst.TRUE  = 1;
cst.FALSE = 0;

%-------------------------------------------------------------------------%
%Type of coordinates
%-------------------------------------------------------------------------%
cst.coord.NC   = 0; % Normalized-Centered coordinates coordinates (x, p) (canonical)
cst.coord.SYS  = 1; % Three-Body Problem coordinates (x, p) (canonical) with the US convention (the Earth has a negative abscissa)
cst.coord.VSYS = 2; % Three-Body Problem (EM or SEM) coordinates (x, v) (non canonical) again with the US convention (the Earth has a negative abscissa)
cst.coord.VSEM = 3; % SEM coordinates (x, v) (non canonical) again with the US convention (the Earth has a negative abscissa)

%-------------------------------------------------------------------------%
%Type of framework
%-------------------------------------------------------------------------%
cst.fwrk.EM  = 0; %Earth-Moon focus
cst.fwrk.SEM = 1; %Sun-(Earth+Moon) focus

%-------------------------------------------------------------------------%
%Type of model
%-------------------------------------------------------------------------%
cst.model.RTBP = 0; %Restricted Three-Body model
cst.model.QBCP = 1; %Quasi-Bicircular model

%-------------------------------------------------------------------------%
%Type of manifold
%-------------------------------------------------------------------------%
cst.mantype.MAN_CENTER = 0;
cst.mantype.MAN_CENTER_S = 1;
cst.mantype.MAN_CENTER_U = 2;
cst.mantype.MAN_CENTER_US = 3;

%-------------------------------------------------------------------------%
% Environment
%-------------------------------------------------------------------------%
cst.env.G = 6.67428e-11; %gravitational constant
cst.env.AU = 1.49597871e8;        %astronomical unit in km
cst.env.julian.y2015 = 2457023.5; %Julian date of 01/01/2015 at 00:00
cst.env.julian.y2000 = 2451544.5; %Julian date of 01/01/2015 at 00:00
cst.env.hours = 3600;             %in seconds
cst.env.days  = 86400;            %in seconds
cst.env.years = 31556926;         %in seconds (true value)

%-------------------------------------------------------------------------%
% Orbit
%-------------------------------------------------------------------------%
% Family
cst.orbit.family.NORTHERN = 'NORTHERN';
cst.orbit.family.SOUTHERN = 'SOUTHERN';
cst.orbit.family.PLANAR   = 'PLANAR';

% Computation type used to built the orbit ('EMPTY' initially)
cst.orbit.EMPTY        =  'EMPTY';
cst.orbit.APPROXIMATED =  'APPROXIMATED';
cst.orbit.REAL         =  'REAL';

% Initialisation of the STM as the identity matrix
cst.orbit.STM0 = eye(6);

%Type
cst.orbit.type.HALO  = 'HALO';
cst.orbit.type.VLYAP = 'VERTICAL_LYAPUNOV';
cst.orbit.type.PLYAP = 'PLANAR_LYAPUNOV';

%-------------------------------------------------------------------------%
% Type of differential corrector
%-------------------------------------------------------------------------%
cst.corr.X0_FIXED = 1;
cst.corr.Z0_FIXED = 2;
cst.corr.MIN_NORM = 3;

%-------------------------------------------------------------------------%
% Type of plot
%-------------------------------------------------------------------------%
cst.plot.ADIM = 0;
cst.plot.DIM  = 1;

%-------------------------------------------------------------------------%
%Manifold
%-------------------------------------------------------------------------%
% EXTERIOR == towards east, INTERIOR == towards west
cst.manifold.EXTERIOR = 1;
cst.manifold.INTERIOR = -1;

% Type
cst.manifold.STABLE   = 1;
cst.manifold.UNSTABLE = -1;

% Type of termination (by event or freely)
cst.manifold.event.type.FREE      = 'FREE';                       %no termination event
cst.manifold.event.type.X_SECTION = 'X_SECTION';                  %termination on a section x = cst
cst.manifold.event.type.Y_SECTION = 'Y_SECTION';                  %termination on a section y = cst
cst.manifold.event.type.Z_SECTION = 'Z_SECTION';                  %termination on a section z = cst
cst.manifold.event.type.ANGLE_SECTION = 'ANGLE_SECTION';          %termination on a section phi = cst, with phi a given angle
cst.manifold.event.type.FLIGHT_PATH_ANGLE = 'FLIGHT_PATH_ANGLE';  %termination on a section fpa = cst

% Is the event terminal?
cst.manifold.event.isterminal.YES = 1;
cst.manifold.event.isterminal.NO = 0;

% Which zeros trigger the event?
cst.manifold.event.direction.INCREASING =  1;
cst.manifold.event.direction.DECREASING = -1;
cst.manifold.event.direction.ALL        =  0;


%-------------------------------------------------------------------------%
%Computation
%-------------------------------------------------------------------------%
cst.computation.MATLAB = 0;    %use MATLAB routines only
cst.computation.MEX = 1;       %use MEX routines

end


