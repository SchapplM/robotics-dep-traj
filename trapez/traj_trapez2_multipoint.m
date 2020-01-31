% Berechne eine Trajektorie aus mehreren Trapezprofilen zwischen gegebenen
% Punkten
% Die Trajektorienzeitpunkte sind nicht äquidistant, sind aber als
% Simulink-Eingabe geeignet.
% 
% Input:
% QL [PXN]
%   P different Joint Positions
% vmax [1x1] oder [1xN]
%   Max Velocity: Same for every axis (1x1) or separate for each axis                
% T2 [1x1] oder [1xN]
%   Anstiegszeit der Geschwindigkeit (bestimmt maximale Beschleunigung)
% T3 [1x1] oder [1xN]
%   Anstiegszeit der Beschleunigung (bestimmt maximalen Ruck)
% T_Abt [1x1]
%   Sample Time
% T_pause [1x1] oder [Px1]
%   Pause nach jeder Teil-Trajektorie (außer der letzten)
% 
% Output:
% q_Arm_Full [MxN]
%   Position Trajectory
% qD_Arm_Full [MxN]
%   Velocity Trajectory 
% qDD_Arm_Full [MxN]
%   Acceleration Trajectory
% t
%   Zeit-Stützstellen der Trajektorie
% IL [NP-1x1]
%   Indizes der eingegebenen Eckwerte aus QL
% 
% TODO: Unterschiedliche Geschwindigkeiten etc. zwischen allen Teilstücken
% (vektorielle Eingabe der Parameter)

% Moritz Schappler, schappler@irt.uni-hannover.de, 2016-10
% (C) Institut für Regelungstechnik, Universität Hannover

function [Q,QD,QDD,t,IL] = traj_trapez2_multipoint(QL, ...
  vmax, T2, T3, T_Abt, T_pause)

Np = size(QL,1);
if all(size(T_pause) == [1 1])
  T_pause = T_pause * ones(Np,1);
end

Q = [];
QD = [];
QDD = [];
IL = NaN(Np-1,1);
t = [];
for i = 2:Np
  % Gelenkwinkeltrajektorie generieren
  [Traj_i_q, Traj_i_qD, Traj_i_qDD, Traj_i_t] = traj_trapez2(QL(i-1,:), QL(i,:), ...
    vmax, T2, T3, T_Abt);
  if i > 2 % Wartezeit für die nächste Trajektorie anhängen
    % Mindestens Abtastzeit, sonst sind zwei Zeitschritte gleich
    Traj_i_t = Traj_i_t + max(T_pause(i-1),T_Abt)+t(end);
  end
  IL(i-1) = size(Q,1)+1;
  % Trajektorie anhängen
  Q = [Q; Traj_i_q]; %#ok<AGROW>
  QD = [QD; Traj_i_qD]; %#ok<AGROW>
  QDD = [QDD; Traj_i_qDD]; %#ok<AGROW>
  t = [t; Traj_i_t]; %#ok<AGROW>
end