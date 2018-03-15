% Gesamttest für das Matlab-Toolbox-Repo
% 
% Führt alle verfügbaren Modultests aus um die Funktionalität
% sicherzustellen

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-03
% (C) Institut für mechatronische Systeme, Universität Hannover

clc; clear; close all;

this_repo_path = fullfile(fileparts(which('trajectory_toolbox_path_init.m')));
addpath(fullfile(this_repo_path, 'examples_tests'));

cone_rotation_trajectory_example
trapezoid_trajectory_example

clc
close all
fprintf('Alle Testfunktionen dieses Repos ausgeführt\n');