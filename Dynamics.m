function Dynamics(block)
%MSFUNTMPL_BASIC A Template for a Level-2 MATLAB S-Function
%   The MATLAB S-function is written as a MATLAB function with the
%   same name as the S-function. Replace 'msfuntmpl_basic' with the 
%   name of your S-function.
%
%   It should be noted that the MATLAB S-function is very similar
%   to Level-2 C-Mex S-functions. You should be able to get more
%   information for each of the block methods by referring to the
%   documentation for C-Mex S-functions.
%
%   Copyright 2003-2010 The MathWorks, Inc.

%%
%% The setup method is used to set up the basic attributes of the
%% S-function such as ports, parameters, etc. Do not add any other
%% calls to the main body of the function.
%%
setup(block);

%endfunction

%% Function: setup ===================================================
%% Abstract:
%%   Set up the basic characteristics of the S-function block such as:
%%   - Input ports
%%   - Output ports
%%   - Dialog parameters
%%   - Options
%%
%%   Required         : Yes
%%   C-Mex counterpart: mdlInitializeSizes
%%
function setup(block)

% Register number of ports
block.NumInputPorts  = 2;
block.NumOutputPorts = 1;

% Setup port properties to be inherited or dynamic (input and output ports
% inherit their compiled properties from the model)
block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;

% Override input port properties
block.InputPort(1).Dimensions        = 4;
block.InputPort(1).DatatypeID        = 0;  % double
block.InputPort(1).Complexity        = 'real';
block.InputPort(1).DirectFeedthrough = false;

block.InputPort(2).Dimensions        = 6;
block.InputPort(2).DatatypeID        = 0;  % double
block.InputPort(2).Complexity        = 'real';
block.InputPort(2).DirectFeedthrough = false;

% Override output port properties
block.OutputPort(1).Dimensions       = 12;
block.OutputPort(1).DatatypeID       = 0; % double
block.OutputPort(1).Complexity       = 'real';

% Register parameters
block.NumDialogPrms     = 9;

% Register sample times
%  [0 offset]            : Continuous sample time
%  [positive_num offset] : Discrete sample time
%
%  [-1, 0]               : Inherited sample time
%  [-2, 0]               : Variable sample time
block.SampleTimes = [0 0];
block.NumContStates = 12;

% Number 
% Specify the block simStateCompliance. The allowed values are:
%    'UnknownSimState', < The default setting; warn and assume DefaultSimState
%    'DefaultSimState', < Same sim state as a built-in block
%    'HasNoSimState',   < No sim state
%    'CustomSimState',  < Has GetSimState and SetSimState methods
%    'DisallowSimState' < Error out when saving or restoring the model sim state
block.SimStateCompliance = 'DefaultSimState';

  % -----------------------------------------------------------------
  % Options
  % -----------------------------------------------------------------
  % Specify if Accelerator should use TLC or call back to the 
  % MATLAB file
  block.SetAccelRunOnTLC(false);

%% -----------------------------------------------------------------
%% The MATLAB S-function uses an internal registry for all
%% block methods. You should register all relevant methods
%% (optional and required) as illustrated below. You may choose
%% any suitable name for the methods and implement these methods
%% as local functions within the same file. See comments
%% provided for each function for more information.
%% -----------------------------------------------------------------

block.RegBlockMethod('InitializeConditions', @InitializeConditions);
block.RegBlockMethod('Outputs',              @Outputs);     % Required
block.RegBlockMethod('Derivatives',          @Derivatives);
% block.RegBlockMethod('Update',               @Update);    

%end setup


%%
%% InitializeConditions:
%%   Functionality    : Called at the start of simulation and if it is 
%%                      present in an enabled subsystem configured to reset 
%%                      states, it will be called when the enabled subsystem
%%                      restarts execution to reset the states.
%%   Required         : No
%%   C-MEX counterpart: mdlInitializeConditions
%%
function InitializeConditions(block)

%%%%%%%%%% Set Initial Conditions
    block.ContStates.Data(1) = 1;                % x1
    block.ContStates.Data(2) = 0;                % x2
    block.ContStates.Data(3) = 0;                % x3
    block.ContStates.Data(4) = 0;                % x4 
    block.ContStates.Data(5) = 0;                % x5
    block.ContStates.Data(6) = 0;                % x6
    block.ContStates.Data(7) = 0;                % x7 
    block.ContStates.Data(8) = 0;                % x8 
    block.ContStates.Data(9) = 0;                % x9 
    block.ContStates.Data(10) = 0;               % x10 
    block.ContStates.Data(11) = 0;               % x11 
    block.ContStates.Data(12) = 0;               % x12 
   
    for i=1:12
        block.OutputPort(1).Data(i) = 0;
    end
%end InitializeConditions

%%
%% Outputs:
%%   Functionality    : Called to generate block outputs in
%%                      simulation step
%%   Required         : Yes
%%   C-MEX counterpart: mdlOutputs
%%
function Outputs(block)
% fileID = fopen('Debug.txt','w');
% fprintf(fileID,'%E\t%E\t%E\t%E\t',block.ContStates.Data(1), block.ContStates.Data(2), block.ContStates.Data(3), block.ContStates.Data(4));
% fprintf(fileID,'%E\t%E\t%E\t%E\t',block.ContStates.Data(5), block.ContStates.Data(6), block.ContStates.Data(7), block.ContStates.Data(8));
% fprintf(fileID,'%E\t%E\t%E\t%E\t',block.ContStates.Data(9), block.ContStates.Data(10), block.ContStates.Data(11), block.ContStates.Data(12));
block.OutputPort(1).Data = block.ContStates.Data;               % Vector of 12 outputs
% block.OutputPort(1).Data


%%
%% Derivatives:
%%   Functionality    : Called to update derivatives of
%%                      continuous states during simulation step
%%   Required         : No
%%   C-MEX counterpart: mdlDerivatives
%%
function Derivatives(block)

%%------Parameters--------------
m = block.DialogPrm(1).Data;
Ix = block.DialogPrm(2).Data;
Iy = block.DialogPrm(3).Data;
Iz = block.DialogPrm(4).Data;
l = block.DialogPrm(5).Data;
Jr = block.DialogPrm(6).Data;
g = block.DialogPrm(7).Data;
Kt = block.DialogPrm(8).Data;
d = block.DialogPrm(9).Data;

%%-----Inputs-----------------
input = block.InputPort(1).Data;
Omega1= input(1);
Omega2= input(2);
Omega3= input(3);
Omega4= input(4);

disturbance = block.InputPort(2).Data;
DistX = disturbance(1);
DistY = disturbance(2);
DistZ = disturbance(3);
DistPhi = disturbance(4);

%% ------ Calculate U's --------
U1 = Kt*(Omega1^2 + Omega2^2 + Omega3^2 + Omega4^2);        % I think omega should be in rev/sec!! Check
U2 = Kt*(Omega4^2 - Omega2^2);
U3 = Kt*(Omega3^2 - Omega1^2);
U4 = d*(Omega1^2 + Omega3^2 - Omega2^2 - Omega4^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Omega_r = Omega1 + Omega3 - Omega2 - Omega4;

%% ---------states------------
%     x1 = block.ContStates.Data(1);             % x              
%     x2 = block.ContStates.Data(2);             % xdot    
%     x3 = block.ContStates.Data(3);             % y
%     x4 = block.ContStates.Data(4);             % ydot  
%     x5 = block.ContStates.Data(5);             % z 
%     x6 = block.ContStates.Data(6);             % zdot 
%     x7 = block.ContStates.Data(7);             % phi 
%     x8 = block.ContStates.Data(8);             % phi_dot  
%     x9 = block.ContStates.Data(9);             % theta
%     x10 = block.ContStates.Data(10);           % theta_dot 
%     x11 = block.ContStates.Data(11);           % psi
%     x12 = block.ContStates.Data(12);           % psi_dot
    
    x        = block.ContStates.Data(1);             % x              
    xdot     = block.ContStates.Data(2);             % xdot    
    y        = block.ContStates.Data(3);             % y
    ydot     = block.ContStates.Data(4);             % ydot  
    z        = block.ContStates.Data(5);             % z 
    zdot     = block.ContStates.Data(6);             % zdot 
    phi      = block.ContStates.Data(7);             % phi 
    phidot   = block.ContStates.Data(8);             % phi_dot  
    theta    = block.ContStates.Data(9);             % theta
    thetadot = block.ContStates.Data(10);           % theta_dot 
    psi      = block.ContStates.Data(11);           % psi
    psidot   = block.ContStates.Data(12);           % psi_dot
    
    %% ------------------ State Equations --------------------
%     x1_dot = x2;                                                            % x dot
%     x2_dot = (cos(x11)*sin(x9)*cos(x7) + sin(x11)*sin(x7))*(U1/m);          % x ddot
%     x3_dot = x4;                                                            % y dot
%     x4_dot = (sin(x11)*sin(x9)*cos(x7) - cos(x11)*sin(x7))*(U1/m);          % y ddot
%     x5_dot = x6;                                                            % z dot
%     x6_dot = (cos(x9)*cos(x7))*(U1/m) - g;                                  % z ddot
%     x7_dot = x8;
%     x8_dot = (1/Ix)*((Iy - Iz)*x12*x10 - Jr*x10*Omega_r + U2*l);
%     x9_dot = x10;
%     x10_dot = (1/Iy)*((Iz - Ix)*x12*x8 + Jr*x8*Omega_r + U3*l);
%     x11_dot = x12;
%     x12_dot = (1/Iz)*((Ix - Iy)*x10*x8 + U4);
    
    x1_dot = xdot;                                                            % x dot
    x2_dot = (cos(psi)*sin(theta)*cos(phi) + sin(psi)*sin(phi))*(U1/m) + DistX;          % x ddot
    x3_dot = ydot;                                                            % y dot
    x4_dot = (sin(psi)*sin(theta)*cos(phi) - cos(psi)*sin(phi))*(U1/m) + DistY;          % y ddot
    x5_dot = zdot ;                                                            % z dot
    x6_dot = (cos(theta)*cos(phi))*(U1/m) - g + DistZ;                                  % z ddot
    x7_dot = phidot;
    x8_dot =  (1/Iy)*((Iy - Iz)*psidot*thetadot - Jr*thetadot*Omega_r + U2*l);         % pi/30 to convert from rev/min to rad/sec !  (1/Ix)*((Iy - Iz)*psidot*thetadot - Jr*thetadot*Omega_r +
    x9_dot = thetadot;
    x10_dot = (1/Iy)*((Iz - Ix)*psidot*phidot + Jr*phidot*Omega_r + U3*l);
    x11_dot = psidot;
    x12_dot = (1/Iz)*((Ix - Iy)*thetadot*phidot + U4);
    %%%%%%%%%%%%%%%%%%%%
    % Rough rule to impose a "ground" boundary...could easily be improved...
%     if ((z <= 0) && (x5_dot <= 0)) % better  version then before?
%         x5_dot = 0;
%         block.ContStates.Data(5) = 0;
%     end
%     if(phi > 1.57079)
%         block.ContStates.Data(7) = 1.57;
%         block.ContStates.Data(8) = 0;
%     elseif (phi < -1.57079)
%         block.ContStates.Data(7) = -1.57;
%         block.ContStates.Data(8) = 0;
%     end
%     if(theta > 1.57079)
%         block.ContStates.Data(9) = 1.57;
%         block.ContStates.Data(10) = 0;
%     elseif (theta < -1.57079)
%         block.ContStates.Data(9) = -1.57;
%         block.ContStates.Data(10) = 0;
%      end
%     %%%
    block.Derivatives.Data = [x1_dot; x2_dot; x3_dot; x4_dot; x5_dot; x6_dot; x7_dot; x8_dot; x9_dot; x10_dot; x11_dot; x12_dot];

%end Derivatives
%%
%%
%function Update(block)
