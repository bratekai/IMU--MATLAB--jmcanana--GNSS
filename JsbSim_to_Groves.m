function [out_profile] = JsbSim_to_Groves(jsb_profile)
%JsbSim_to_Groves - Converts motions profile from jsbSim to Groves profile
%
% Software for use with "An Open Source Flight Dynamics Model
% and IMUSignal Simulator."
%
% This function created 1/4/2018 by James McAnanama
%
% Inputs:
%   jsb_profile T,P,V,A from jsbSim UDP port.
%               See https://github.com/jmcanana/JGM_PLANS_2018 for jsbSim setup.
%
% Outputs:
%   out_profile   Navigation solution as a motion profile array
%
%% Format of jsbSim profile:
%  Column 1: time (sec)
%  Column 2: latitude (geod, deg)
%  Column 3: longitude (deg)
%  Column 4: height (ft)
%  Column 5: north velocity (fpd)
%  Column 6: east velocity (fps)
%  Column 7: down velocity (fps)
%  Column 8: roll angle of body w.r.t NED (rad)
%  Column 9: pitch angle of body w.r.t NED (rad)
%  Column 10: yaw angle of body w.r.t NED (rad)

%
% Format of Groves profile:
%  Column 1: time (sec)
%  Column 2: latitude (rad)
%  Column 3: longitude (rad)
%  Column 4: height (m)
%  Column 5: north velocity (m/s)
%  Column 6: east velocity (m/s)
%  Column 7: down velocity (m/s)
%  Column 8: roll angle of body w.r.t NED (rad)
%  Column 9: pitch angle of body w.r.t NED (rad)
%  Column 10: yaw angle of body w.r.t NED (rad)
%
deg_to_rad = 0.01745329252;
ft_to_m    = 0.3048;
[no_epochs,no_columns] = size(jsb_profile);
if no_columns ~= 10
    disp('Input profile has the wrong number of columns')
else
    out_profile = jsb_profile;
    out_profile(:, 2:3) = out_profile(:, 2:3) * deg_to_rad;
    out_profile(:, 4:7) = out_profile(:, 4:7) * ft_to_m;
end

