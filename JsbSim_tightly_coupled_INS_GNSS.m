function [out_profile,out_errors,out_IMU_bias_est,out_clock,out_KF_SD] =...
                                          JsbSim_tightly_coupled_INS_GNSS(lport)
%JsbSim_tightly_coupled_INS_GNSS - Simulates inertial navigation using ECEF
% navigation equations and kinematic model, GNSS and tightly coupled
% INS/GNSS integration, with jsbSim flight dynamics as the input profile.
%
% Software for use with "Principles of GNSS, Inertial, and Multisensor
% Integrated Navigation Systems," Second Edition, extended to use jsbSim
% flight dynamics data as described in "An Open Source Flight Dynamics Model
% and IMUSignal Simulator."
%
% Dependencies:
% tcp_udp_ip: A UDP to matlab utility, copyright (c) 2015, Peter Rydesäter.
% jsbSim: JSBSim: An open source flight dynamics model.
% See https://github.com/jmcanana/JGM_PLANS_2018 for sources.
%
%
% This function created 1/4/2018 by James McAnanama, based on
% the Tightly_coupled_INS_GNSSInertial_navigation function
% created 12/4/2012 by Paul Groves
%
% Inputs:
%   lport       UDP port of incoming jsbSim data.
%               See https://github.com/jmcanana/JGM_PLANS_2018 for jsbSim setup.
%
%
% Assumptions for the following are taken as the tactical-grade IMU model used
% in the INS_GNSS_Demo_10 from this repo.
%   initialization_errors
%     .delta_r_eb_n     position error resolved along NED (m)
%     .delta_v_eb_n     velocity error resolved along NED (m/s)
%     .delta_eul_nb_n   attitude error as NED Euler angles (rad)
%   IMU_errors
%     .delta_r_eb_n     position error resolved along NED (m)
%     .b_a              Accelerometer biases (m/s^2)
%     .b_g              Gyro biases (rad/s)
%     .M_a              Accelerometer scale factor and cross coupling errors
%     .M_g              Gyro scale factor and cross coupling errors
%     .G_g              Gyro g-dependent biases (rad-sec/m)
%     .accel_noise_root_PSD   Accelerometer noise root PSD (m s^-1.5)
%     .gyro_noise_root_PSD    Gyro noise root PSD (rad s^-0.5)
%     .accel_quant_level      Accelerometer quantization level (m/s^2)
%     .gyro_quant_level       Gyro quantization level (rad/s)
%   GNSS_config
%     .epoch_interval     Interval between GNSS epochs (s)
%     .init_est_r_ea_e    Initial estimated position (m; ECEF)
%     .no_sat             Number of satellites in constellation
%     .r_os               Orbital radius of satellites (m)
%     .inclination        Inclination angle of satellites (deg)
%     .const_delta_lambda Longitude offset of constellation (deg)
%     .const_delta_t      Timing offset of constellation (s)
%     .mask_angle         Mask angle (deg)
%     .SIS_err_SD         Signal in space error SD (m)
%     .zenith_iono_err_SD Zenith ionosphere error SD (m)
%     .zenith_trop_err_SD Zenith troposphere error SD (m)
%     .code_track_err_SD  Code tracking error SD (m)
%     .rate_track_err_SD  Range rate tracking error SD (m/s)
%     .rx_clock_offset    Receiver clock offset at time=0 (m)
%     .rx_clock_drift     Receiver clock drift at time=0 (m/s)
%   TC_KF_config
%     .init_att_unc           Initial attitude uncertainty per axis (rad)
%     .init_vel_unc           Initial velocity uncertainty per axis (m/s)
%     .init_pos_unc           Initial position uncertainty per axis (m)
%     .init_b_a_unc           Initial accel. bias uncertainty (m/s^2)
%     .init_b_g_unc           Initial gyro. bias uncertainty (rad/s)
%     .init_clock_offset_unc  Initial clock offset uncertainty per axis (m)
%     .init_clock_drift_unc   Initial clock drift uncertainty per axis (m/s)
%     .gyro_noise_PSD         Gyro noise PSD (rad^2/s)
%     .accel_noise_PSD        Accelerometer noise PSD (m^2 s^-3)
%     .accel_bias_PSD         Accelerometer bias random walk PSD (m^2 s^-5)
%     .gyro_bias_PSD          Gyro bias random walk PSD (rad^2 s^-3)
%     .clock_freq_PSD         Receiver clock frequency-drift PSD (m^2/s^3)
%     .clock_phase_PSD        Receiver clock phase-drift PSD (m^2/s)
%     .pseudo_range_SD        Pseudo-range measurement noise SD (m)
%     .range_rate_SD          Pseudo-range rate measurement noise SD (m/s)
%
% Outputs:
%   out_profile        Navigation solution as a motion profile array
%   out_errors         Navigation solution error array
%   out_IMU_bias_est   Kalman filter IMU bias estimate array
%   out_clock          GNSS Receiver clock estimate array
%   out_KF_SD          Output Kalman filter state uncertainties
%
% Format of motion profiles:
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
% Format of error array:
%  Column 1: time (sec)
%  Column 2: north position error (m)
%  Column 3: east position error (m)
%  Column 4: down position error (m)
%  Column 5: north velocity error (m/s)
%  Column 6: east velocity error (m/s)
%  Column 7: down velocity (error m/s)
%  Column 8: attitude error about north (rad)
%  Column 9: attitude error about east (rad)
%  Column 10: attitude error about down = heading error  (rad)
%
% Format of output IMU biases array:
%  Column 1: time (sec)
%  Column 2: estimated X accelerometer bias (m/s^2)
%  Column 3: estimated Y accelerometer bias (m/s^2)
%  Column 4: estimated Z accelerometer bias (m/s^2)
%  Column 5: estimated X gyro bias (rad/s)
%  Column 6: estimated Y gyro bias (rad/s)
%  Column 7: estimated Z gyro bias (rad/s)
%
% Format of receiver clock array:
%  Column 1: time (sec)
%  Column 2: estimated clock offset (m)
%  Column 3: estimated clock drift (m/s)
%
% Format of KF state uncertainties array:
%  Column 1: time (sec)
%  Column 2: X attitude error uncertainty (rad)
%  Column 3: Y attitude error uncertainty (rad)
%  Column 4: Z attitude error uncertainty (rad)
%  Column 5: X velocity error uncertainty (m/s)
%  Column 6: Y velocity error uncertainty (m/s)
%  Column 7: Z velocity error uncertainty (m/s)
%  Column 8: X position error uncertainty (m)
%  Column 9: Y position error uncertainty (m)
%  Column 10: Z position error uncertainty (m)
%  Column 11: X accelerometer bias uncertainty (m/s^2)
%  Column 12: Y accelerometer bias uncertainty (m/s^2)
%  Column 13: Z accelerometer bias uncertainty (m/s^2)
%  Column 14: X gyro bias uncertainty (rad/s)
%  Column 15: Y gyro bias uncertainty (rad/s)
%  Column 16: Z gyro bias uncertainty (rad/s)
%  Column 17: clock offset uncertainty (m)
%  Column 18: clock drift uncertainty (m/s)

% Copyright 2012, Paul Groves
% Copyright 2018, James McAnanama
% License: BSD; see license.txt for details


  % Constants
  deg_to_rad = 0.01745329252;
  rad_to_deg = 1/deg_to_rad;
  micro_g_to_meters_per_second_squared = 9.80665E-6;

  %jsbSim streaming at 50Hz, the following caps the data history in the plot:
  data_rate = 50;
  seconds_of_data = data_rate * 300;


  use_old_GNSS_data = false;
  prev_no_GNSS_meas = 0;


  % CONFIGURATION
  % Output motion profile and error filenames
  output_profile_name = 'JsbSim_INS_GNSS_Demo_TC_Profile.csv';
  output_errors_name = 'JsbSim_INS_GNSS_Demo_TC_Errors.csv';

  % Attitude initialization error (deg, converted to rad; @N,E,D)
  initialization_errors.delta_eul_nb_n = [-0.05;0.04;1]*deg_to_rad; % rad

  % Accelerometer biases (micro-g, converted to m/s^2; body axes)
  IMU_errors.b_a = [900;-1300;800] * micro_g_to_meters_per_second_squared;
  % Gyro biases (deg/hour, converted to rad/sec; body axes)
  IMU_errors.b_g = [-9;13;-8] * deg_to_rad / 3600;
  % Accelerometer scale factor and cross coupling errors (ppm, converted to
  % unitless; body axes)
  IMU_errors.M_a = [500, -300, 200;...
                   -150, -600, 250;...
                   -250,  100, 450] * 1E-6;
  % Gyro scale factor and cross coupling errors (ppm, converted to unitless;
  % body axes)
  IMU_errors.M_g = [400, -300,  250;...
                      0, -300, -150;...
                      0,    0, -350] * 1E-6;
  % Gyro g-dependent biases (deg/hour/g, converted to rad-sec/m; body axes)
  IMU_errors.G_g = [0.9, -1.1, -0.6;...
                   -0.5,  1.9, -1.6;...
                    0.3,  1.1, -1.3] * deg_to_rad / (3600 * 9.80665);
  % Accelerometer noise root PSD (micro-g per root Hz, converted to m s^-1.5)
  IMU_errors.accel_noise_root_PSD = 100 *...
      micro_g_to_meters_per_second_squared;
  % Gyro noise root PSD (deg per root hour, converted to rad s^-0.5)
  IMU_errors.gyro_noise_root_PSD = 0.01 * deg_to_rad / 60;
  % Accelerometer quantization level (m/s^2)
  IMU_errors.accel_quant_level = 1E-2;
  % Gyro quantization level (rad/s)
  IMU_errors.gyro_quant_level = 2E-4;

  % Interval between GNSS epochs (s)
  GNSS_config.epoch_interval = 0.5;

  % Initial estimated position (m; ECEF)
  GNSS_config.init_est_r_ea_e = [0;0;0];

  % Number of satellites in constellation
  GNSS_config.no_sat = 30;
  % Orbital radius of satellites (m)
  GNSS_config.r_os = 2.656175E7;
  % Inclination angle of satellites (deg)
  GNSS_config.inclination = 55;
  % Longitude offset of constellation (deg)
  GNSS_config.const_delta_lambda = 0;
  % Timing offset of constellation (s)
  GNSS_config.const_delta_t = 0;

  % Mask angle (deg)
  GNSS_config.mask_angle = 10;
  % Signal in space error SD (m) *Give residual where corrections are applied
  GNSS_config.SIS_err_SD = 1;
  % Zenith ionosphere error SD (m) *Give residual where corrections are applied
  GNSS_config.zenith_iono_err_SD = 2;
  % Zenith troposphere error SD (m) *Give residual where corrections are applied
  GNSS_config.zenith_trop_err_SD = 0.2;
  % Code tracking error SD (m) *Can extend to account for multipath
  GNSS_config.code_track_err_SD = 1;
  % Range rate tracking error SD (m/s) *Can extend to account for multipath
  GNSS_config.rate_track_err_SD = 0.02;
  % Receiver clock offset at time=0 (m);
  GNSS_config.rx_clock_offset = 10000;
  % Receiver clock drift at time=0 (m/s);
  GNSS_config.rx_clock_drift = 100;

  % Initial attitude uncertainty per axis (deg, converted to rad)
  TC_KF_config.init_att_unc = degtorad(1);
  % Initial velocity uncertainty per axis (m/s)
  TC_KF_config.init_vel_unc = 0.1;
  % Initial position uncertainty per axis (m)
  TC_KF_config.init_pos_unc = 10;
  % Initial accelerometer bias uncertainty per instrument (micro-g, converted
  % to m/s^2)
  TC_KF_config.init_b_a_unc = 1000 * micro_g_to_meters_per_second_squared;
  % Initial gyro bias uncertainty per instrument (deg/hour, converted to rad/sec)
  TC_KF_config.init_b_g_unc = 10 * deg_to_rad / 3600;
  % Initial clock offset uncertainty per axis (m)
  TC_KF_config.init_clock_offset_unc = 10;
  % Initial clock drift uncertainty per axis (m/s)
  TC_KF_config.init_clock_drift_unc = 0.1;

  % Gyro noise PSD (deg^2 per hour, converted to rad^2/s)
  TC_KF_config.gyro_noise_PSD = (0.02 * deg_to_rad / 60)^2;
  % Accelerometer noise PSD (micro-g^2 per Hz, converted to m^2 s^-3)
  TC_KF_config.accel_noise_PSD = (200 *...
      micro_g_to_meters_per_second_squared)^2;
  % Accelerometer bias random walk PSD (m^2 s^-5)
  TC_KF_config.accel_bias_PSD = 1.0E-7;
  % Gyro bias random walk PSD (rad^2 s^-3)
  TC_KF_config.gyro_bias_PSD = 2.0E-12;
  % Receiver clock frequency-drift PSD (m^2/s^3)
  TC_KF_config.clock_freq_PSD = 1;
  % Receiver clock phase-drift PSD (m^2/s)
  TC_KF_config.clock_phase_PSD = 1;

  % Pseudo-range measurement noise SD (m)
  TC_KF_config.pseudo_range_SD = 2.5;
  % Pseudo-range rate measurement noise SD (m/s)
  TC_KF_config.range_rate_SD = 0.1;

  % Seeding of the random number generator for reproducability. Change
  % this value for a different random number sequence (may not work in Octave).
  RandStream.setGlobalStream(RandStream('mt19937ar','seed',1));

  %If you want to look at data after the fact, use no_epochs and increment
  %the epoch through the main loop.  Default is to just plot the results.
  no_epochs = 1;
  epoch = 1;

  % Initialize profile records and errors record
  in_profile = zeros(no_epochs,10);
  out_profile = zeros(no_epochs,10);
  out_errors = zeros(no_epochs,10);

  % Begins

  %% Open udp port
  % Add default argument
  if nargin<1, lport=1138; end
  % Write help message
      disp ' ';
      disp 'Run jsbSim to send motion profile.';
      disp(sprintf('  ./$JSBSIM_ROOT/./linux-build/src/JSBSim --script=scripts/c310_plans_racetrack.xml --realtime --logdirectivefile=data_output/matlabSocket.xml '));
      disp ' ';
  % Open  udpsocket and bind udp port adress to it.
  udp=pnet('udpsocket',lport);

  init_nav = true;

  % Open figure window
  f  = figure(1);
  f2 = figure(2);
  f3 = figure(3);

  if 1%try
      while epoch <= no_epochs
          % Wait/Read udp packet to reed buffer
          len=pnet(udp,'readpacket');
          if len <= 0
              continue
          end
          % readline and convert from char to doubles
          data=pnet(udp,'readline');
          C = strsplit(data, ',');
          X =  str2double(C);
          if length(X) ~= 10
              disp('Input profile has the wrong number of columns')
              continue
          end
          in_profile(epoch, :) = JsbSim_to_Groves(X);
          if init_nav
              init_nav = false;
              % Initialize true navigation solution
              min_time = in_profile(epoch,1);
              old_time = in_profile(epoch,1) - min_time;
              true_L_b = in_profile(epoch,2);
              true_lambda_b = in_profile(epoch,3);
              true_h_b = in_profile(epoch,4);
              true_v_eb_n = in_profile(epoch,5:7)';
              true_eul_nb = in_profile(epoch,8:10)';
              true_C_b_n = Euler_to_CTM(true_eul_nb)';
              [old_true_r_eb_e,old_true_v_eb_e,old_true_C_b_e] =...
                  NED_to_ECEF(true_L_b,true_lambda_b,true_h_b,true_v_eb_n,true_C_b_n);

              % Determine satellite positions and velocities
              [sat_r_es_e,sat_v_es_e] = Satellite_positions_and_velocities(old_time,...
                  GNSS_config);

              % Initialize the GNSS biases. Note that these are assumed constant throughout
              % the simulation and are based on the initial elevation angles. Therefore,
              % this function is unsuited to simulations longer than about 30 min.
              GNSS_biases = Initialize_GNSS_biases(sat_r_es_e,old_true_r_eb_e,true_L_b,...
                  true_lambda_b,GNSS_config);

              % Generate GNSS measurements
              [GNSS_measurements,no_GNSS_meas] = Generate_GNSS_measurements(old_time,...
                  sat_r_es_e,sat_v_es_e,old_true_r_eb_e,true_L_b,true_lambda_b,...
                  old_true_v_eb_e,GNSS_biases,GNSS_config);

              % Determine Least-squares GNSS position solution
              [old_est_r_eb_e,old_est_v_eb_e,est_clock] = GNSS_LS_position_velocity(...
                  GNSS_measurements,no_GNSS_meas,GNSS_config.init_est_r_ea_e,[0;0;0]);
              [old_est_L_b,old_est_lambda_b,old_est_h_b,old_est_v_eb_n] =...
                  pv_ECEF_to_NED(old_est_r_eb_e,old_est_v_eb_e);
              est_L_b = old_est_L_b;

              % Initialize estimated attitude solution
              old_est_C_b_n = Initialize_NED_attitude(true_C_b_n,initialization_errors);
              [temp1,temp2,old_est_C_b_e] = NED_to_ECEF(old_est_L_b,...
                  old_est_lambda_b,old_est_h_b,old_est_v_eb_n,old_est_C_b_n);

              % Initialize output profile record and errors record
              out_profile = zeros(no_epochs,10);
              out_errors = zeros(no_epochs,10);

              % Generate output profile record
              out_profile(1,1) = old_time;
              out_profile(1,2) = old_est_L_b;
              out_profile(1,3) = old_est_lambda_b;
              out_profile(1,4) = old_est_h_b;
              out_profile(1,5:7) = old_est_v_eb_n';
              out_profile(1,8:10) = CTM_to_Euler(old_est_C_b_n')';

              % Determine errors and generate output record
              [delta_r_eb_n,delta_v_eb_n,delta_eul_nb_n] = Calculate_errors_NED(...
                  old_est_L_b,old_est_lambda_b,old_est_h_b,old_est_v_eb_n,old_est_C_b_n,...
                  true_L_b,true_lambda_b,true_h_b,true_v_eb_n,true_C_b_n);
              out_errors(1,1) = old_time;
              out_errors(1,2:4) = delta_r_eb_n';
              out_errors(1,5:7) = delta_v_eb_n';
              out_errors(1,8:10) = delta_eul_nb_n';

              % Initialize Kalman filter P matrix and IMU bias states
              P_matrix = Initialize_TC_P_matrix(TC_KF_config);
              est_IMU_bias = zeros(6,1);

              % Initialize IMU quantization residuals
              quant_residuals = [0;0;0;0;0;0];

              % Generate IMU bias and clock output records
              out_IMU_bias_est(1,1) = old_time;
              out_IMU_bias_est(1,2:7) = est_IMU_bias';
              out_clock(1,1) = old_time;
              out_clock(1,2:3) = est_clock;

              % Generate KF uncertainty record
              out_KF_SD(1,1) = old_time;
              for i =1:17
                  out_KF_SD(1,i+1) = sqrt(P_matrix(i,i));
              end % for i

              % Initialize GNSS model timing
              time_last_GNSS = old_time;
              GNSS_epoch = 1;

              figure(f);
              clf
              h_truth = animatedline(in_profile(1, 3) * rad_to_deg, ...
                                     in_profile(1, 2) * rad_to_deg, ...
                                    'Color','g','LineWidth',1, ...
                                    'LineStyle', '--', ...
                                    'Marker', 'o', ...
                                    'MarkerSize', 4, ...
                                    'MaximumNumPoints', seconds_of_data);
              h_nav  = animatedline(out_profile(1, 3) * rad_to_deg, ...
                                     out_profile(1, 2) * rad_to_deg, ...
                                    'Color','b','LineWidth',1, ...
                                    'LineStyle', '--', ...
                                    'Marker', 'x', ...
                                    'MarkerSize', 4, ...
                                    'MaximumNumPoints', seconds_of_data);
              legend('truth', 'INS');
              xlabel('Longitude');
              ylabel('Latitude');
              title('INS Output');

              figure(f2);
              clf
              subplot(2,1,1);
              h_err_pn = animatedline(0, 0, ...
                                    'Color','k','LineWidth',1, ...
                                    'LineStyle', '-', ...
                                    'Marker', 'x', ...
                                    'MarkerSize', 4, ...
                                    'MaximumNumPoints', seconds_of_data);
              h_err_pe = animatedline(0, 0, ...
                                    'Color','k','LineWidth',1, ...
                                    'LineStyle', '-', ...
                                    'Marker', 'o', ...
                                    'MarkerSize', 4, ...
                                    'MaximumNumPoints', seconds_of_data);
              h_err_pd = animatedline(0, 0, ...
                                    'Color','k','LineWidth',1, ...
                                    'LineStyle', '-', ...
                                    'Marker', '^', ...
                                    'MarkerSize', 4, ...
                                    'MaximumNumPoints', seconds_of_data);

              title('INS Errors');
              legend('Error_N', 'Error_E', 'Error_D');
              ylabel('Position Error [m]');

              subplot(2,1,2);
              h_err_vn = animatedline(0, 0, ...
                                    'Color','k','LineWidth',1, ...
                                    'LineStyle', '-', ...
                                    'Marker', 'x', ...
                                    'MarkerSize', 4, ...
                                    'MaximumNumPoints', seconds_of_data);

              h_err_ve = animatedline(0, 0, ...
                                    'Color','k','LineWidth',1, ...
                                    'LineStyle', '-', ...
                                    'Marker', 'o', ...
                                    'MarkerSize', 4, ...
                                    'MaximumNumPoints', seconds_of_data);

              h_err_vd = animatedline(0, 0, ...
                                    'Color','k','LineWidth',1, ...
                                    'LineStyle', '-', ...
                                    'Marker', '^', ...
                                    'MarkerSize', 4, ...
                                    'MaximumNumPoints', seconds_of_data);


              legend('Error_N', 'Error_E', 'Error_D');
              ylabel('Velocity Error [m/s]');
              xlabel('Time (s)');


              figure(f3);
              clf
              h_err_psi = animatedline(0, 0, ...
                                    'Color','k','LineWidth',1, ...
                                    'LineStyle', '-', ...
                                    'Marker', '^', ...
                                    'MarkerSize', 4, ...
                                    'MaximumNumPoints', seconds_of_data);

              title('INS Heading Errors');
              legend('Error_{Heading}');
              ylabel('Heading Error [mrad]');
              xlabel('Time (s)');


          else
              % Input data from motion profile
              time = in_profile(epoch, 1) - min_time;
              true_L_b = in_profile(epoch, 2);
              true_lambda_b = in_profile(epoch, 3);
              true_h_b = in_profile(epoch, 4);
              true_v_eb_n = in_profile(epoch, 5:7)';
              true_eul_nb = in_profile(epoch, 8:10)';
              true_C_b_n = Euler_to_CTM(true_eul_nb)';
              [true_r_eb_e,true_v_eb_e,true_C_b_e] =...
                  NED_to_ECEF(true_L_b,true_lambda_b,true_h_b,true_v_eb_n,true_C_b_n);

              % Time interval
              tor_i = time - old_time;

              % Calculate specific force and angular rate
              [true_f_ib_b,true_omega_ib_b] = Kinematics_ECEF(tor_i,true_C_b_e,...
                  old_true_C_b_e,true_v_eb_e,old_true_v_eb_e,old_true_r_eb_e);

              % Simulate IMU errors
              [meas_f_ib_b,meas_omega_ib_b,quant_residuals] = IMU_model(tor_i,...
                  true_f_ib_b,true_omega_ib_b,IMU_errors,quant_residuals);

              % Correct IMU errors
              meas_f_ib_b = meas_f_ib_b - est_IMU_bias(1:3);
              meas_omega_ib_b = meas_omega_ib_b - est_IMU_bias(4:6);

              % Update estimated navigation solution
              [est_r_eb_e,est_v_eb_e,est_C_b_e] = Nav_equations_ECEF(tor_i,...
                  old_est_r_eb_e,old_est_v_eb_e,old_est_C_b_e,meas_f_ib_b,...
                  meas_omega_ib_b);

              % Determine whether to update GNSS simulation and run Kalman filter
              update_error_plot = false;
              if (time - time_last_GNSS) >= GNSS_config.epoch_interval
                  GNSS_epoch = GNSS_epoch + 1;
                  update_error_plot = true;
                  tor_s = time - time_last_GNSS;  % KF time interval
                  time_last_GNSS = time;

                  % Determine satellite positions and velocities
                  [sat_r_es_e,sat_v_es_e] = Satellite_positions_and_velocities(time,...
                      GNSS_config);

                  % Generate GNSS measurements
                  [GNSS_measurements,no_GNSS_meas] = Generate_GNSS_measurements(...
                      time,sat_r_es_e,sat_v_es_e,true_r_eb_e,true_L_b,true_lambda_b,...
                      true_v_eb_e,GNSS_biases,GNSS_config);

                  if (use_old_GNSS_data && prev_no_GNSS_meas ~= 0)
                      this_GNSS_measurements = prev_GNSS_measurements;
                      this_no_GNSS_meas      = prev_no_GNSS_meas;
                  else
                      this_GNSS_measurements = GNSS_measurements;
                      this_no_GNSS_meas      = no_GNSS_meas;
                  end
                  % Run Integration Kalman filter
                  [est_C_b_e,est_v_eb_e,est_r_eb_e,est_IMU_bias,est_clock,P_matrix] =...
                      TC_KF_Epoch(this_GNSS_measurements,this_no_GNSS_meas,tor_s,...
                      est_C_b_e, est_v_eb_e,est_r_eb_e,est_IMU_bias,est_clock,P_matrix,...
                      meas_f_ib_b,est_L_b,TC_KF_config);

                  prev_GNSS_measurements = GNSS_measurements;
                  prev_no_GNSS_meas = no_GNSS_meas;

                  % Generate IMU bias and clock output records
                  out_IMU_bias_est(GNSS_epoch,1) = time;
                  out_IMU_bias_est(GNSS_epoch,2:7) = est_IMU_bias';
                  out_clock(GNSS_epoch,1) = time;
                  out_clock(GNSS_epoch,2:3) = est_clock;

                  % Generate KF uncertainty output record
                  out_KF_SD(GNSS_epoch,1) = time;
                  for i =1:17
                      out_KF_SD(GNSS_epoch,i+1) = sqrt(P_matrix(i,i));
                  end % for i

              end % if time

              % Convert navigation solution to NED
              [est_L_b,est_lambda_b,est_h_b,est_v_eb_n,est_C_b_n] =...
                  ECEF_to_NED(est_r_eb_e,est_v_eb_e,est_C_b_e);

              % Generate output profile record
              out_profile(epoch,1) = time;
              out_profile(epoch,2) = est_L_b;
              out_profile(epoch,3) = est_lambda_b;
              out_profile(epoch,4) = est_h_b;
              out_profile(epoch,5:7) = est_v_eb_n';
              out_profile(epoch,8:10) = CTM_to_Euler(est_C_b_n')';

              % Determine errors and generate output record
              [delta_r_eb_n,delta_v_eb_n,delta_eul_nb_n] = Calculate_errors_NED(...
                  est_L_b,est_lambda_b,est_h_b,est_v_eb_n,est_C_b_n,true_L_b,...
                  true_lambda_b,true_h_b,true_v_eb_n,true_C_b_n);
              out_errors(epoch,1) = time;
              out_errors(epoch,2:4) = delta_r_eb_n';
              out_errors(epoch,5:7) = delta_v_eb_n';
              out_errors(epoch,8:10) = delta_eul_nb_n';

              % Reset old values
              old_time = time;
              old_true_r_eb_e = true_r_eb_e;
              old_true_v_eb_e = true_v_eb_e;
              old_true_C_b_e = true_C_b_e;
              old_est_r_eb_e = est_r_eb_e;
              old_est_v_eb_e = est_v_eb_e;
              old_est_C_b_e = est_C_b_e;


              if  update_error_plot
                  addpoints(h_truth, in_profile(epoch, 3) * rad_to_deg, ...
                      in_profile(epoch, 2) * rad_to_deg);
                  addpoints(h_nav,  out_profile(epoch, 3) * rad_to_deg, ...
                      out_profile(epoch, 2) * rad_to_deg);


                  addpoints(h_err_pn, out_errors(epoch, 1), ...
                                                          out_errors(epoch, 2));
                  addpoints(h_err_pe, out_errors(epoch, 1), ...
                                                          out_errors(epoch, 3));
                  addpoints(h_err_pd, out_errors(epoch, 1), ...
                                                          out_errors(epoch, 4));
                  addpoints(h_err_vn, out_errors(epoch, 1), ...
                                                          out_errors(epoch, 5));
                  addpoints(h_err_ve, out_errors(epoch, 1), ...
                                                          out_errors(epoch, 6));
                  addpoints(h_err_vd, out_errors(epoch, 1), ...
                                                          out_errors(epoch, 7));

                  addpoints(h_err_psi, out_errors(epoch, 1), ...
                                                          out_errors(epoch, 10) * 10^3);


              end

              drawnow;
              if no_epochs > 1
                  epoch =  epoch + 1;
              end
          end %if init_nav
      end %while loop
  end %try
  % On break or error close udpconnection and figure window.
  pnet(udp,'close');
  delete(f);
  delete(f2);
  delete(f3);
  return;



% Ends
