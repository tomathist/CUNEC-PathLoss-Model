%% Parameters for alpha-beta model
% LOS
alpha_LOS               =  1.58;
Delta_LOS               =  6;
sigma_S_LOS             =  1.2;
corr_dist_LOS           =  9;

% NLOS
alpha_NLOS              =  6.3;
Delta_NLOS              = -56.2;
sigma_S_NLOS            =  11.5;
corr_dist_NLOS          =  1;

%% Parameters for CUNEC model
% Free space path loss
FSPL_1m_3pt5GHz         =  43.32;

% Environment parameters (constant)
corner_dist_1_threshold = 70;

% Environment parameters (mean)
building_mean           = 29;
height_mean             = 23;
width_mean              = 19;

% Environment parameters (standard deviation)
building_stdev          = 13;
height_stdev            = 10;
width_stdev             = 5;

% Zeroth (mean)
a_0_mean                =  1.56;
Delta_0_mean            =  6.27;
sigma_S_0_AP_mean       =  1.05;
sigma_S_0_UE_mean       =  0.91;
corr_dist_0_AP_mean     =  17.8;
corr_dist_0_UE_mean     =  12.4;

mu_0 = [
    building_mean;
    width_mean;
    height_mean;
    a_0_mean;
    Delta_0_mean;
    sigma_S_0_AP_mean;
    sigma_S_0_UE_mean;
    corr_dist_0_AP_mean;
    corr_dist_0_UE_mean
];

% Zeroth (standard deviation)
a_0_stdev               =  0.05;
Delta_0_stdev           =  0.64;
sigma_S_0_AP_stdev      =  0.11;
sigma_S_0_UE_stdev      =  0.11;
corr_dist_0_AP_stdev    =  6.1;
corr_dist_0_UE_stdev    =  2.8;

sigma_0 = [
    building_stdev;
    width_stdev;
    height_stdev;
    a_0_stdev;
    Delta_0_stdev;
    sigma_S_0_AP_stdev;
    sigma_S_0_UE_stdev;
    corr_dist_0_AP_stdev;
    corr_dist_0_UE_stdev
];

% First (mean)
a_1_mean                =  1.4;
Delta_1_mean            =  29.7;
sigma_S_1_AP_mean       =  4.78;
sigma_S_1_UE_mean       =  4.66;
corr_dist_1_AP_mean     =  9.9;
corr_dist_1_UE_mean     =  14.4;
kappa_1_mean            =  0.037;
offset_1_mean           =  9.179;

mu_1 = [
    building_mean;
    height_mean;
    a_1_mean;
    Delta_1_mean;
    kappa_1_mean;
    offset_1_mean;
    sigma_S_1_AP_mean;
    sigma_S_1_UE_mean;
    corr_dist_1_AP_mean;
    corr_dist_1_UE_mean
];

% First (standard deviation)
a_1_stdev               =  0.65;
Delta_1_stdev           =  10.99;
sigma_S_1_AP_stdev      =  2.6;
sigma_S_1_UE_stdev      =  1.63;
corr_dist_1_AP_stdev    =  6.7;
corr_dist_1_UE_stdev    =  9.5;
kappa_1_stdev           =  0.02;
offset_1_stdev          =  4.63;

sigma_1 = [
    building_stdev;
    height_stdev;
    a_1_stdev;
    Delta_1_stdev;
    kappa_1_stdev;
    offset_1_stdev;
    sigma_S_1_AP_stdev;
    sigma_S_1_UE_stdev;
    corr_dist_1_AP_stdev;
    corr_dist_1_UE_stdev
];

% Second (mean)
a_2_mean                =  1.3;
sigma_S_2_AP_mean       =  7.02;
sigma_S_2_UE_mean       =  7.01;
corr_dist_2_AP_mean     =  10.4;
corr_dist_2_UE_mean     =  16.0;
% 
mu_2 = [
    building_mean;
    a_2_mean;
    sigma_S_2_AP_mean;
    sigma_S_2_UE_mean;
    corr_dist_2_AP_mean;
    corr_dist_2_UE_mean
];
% 
% Second (standard deviation)
a_2_stdev               =  0.58;
sigma_S_2_AP_stdev      =  2.40;
sigma_S_2_UE_stdev      =  1.91;
corr_dist_2_AP_stdev    =  10.4;
corr_dist_2_UE_stdev    =  9.2;
% 
sigma_2 = [
    building_stdev;
    a_2_stdev;
    sigma_S_2_AP_stdev;
    sigma_S_2_UE_stdev;
    corr_dist_2_AP_stdev;
    corr_dist_2_UE_stdev
];
