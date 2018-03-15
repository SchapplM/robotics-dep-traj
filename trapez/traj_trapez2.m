% Berechne Beschleunigungs-Trapezprofil für die Bewegung mehrerer Achsen
% Erhält Randbedingungen für N Achsen und berechnet daraus M Zeitwerte mit
% einem den RB entsprechenden Trapezprofil zweiter Ordnung.
% 
% 
% Input:
% q_Start [1xN]
%   Start Joint Positions [rad]
% q_End [1xN]
%   End Joint Positions [rad]
% vmax [1x1] oder [1xN]
%   Max Velocity: Same for every axis (1x1) or separate for each axis                
% T2 [1x1] oder [1xN]
%   Anstiegszeit der Geschwindigkeit (bestimmt maximale Beschleunigung)
% T3 [1x1] oder [1xN]
%   Anstiegszeit der Beschleunigung (bestimmt maximalen Ruck)
% T_Abt [1x1]
%   Sample Time
% 
% Output:
% q_Arm_Full [MxN]
%   Position Trajectory
% qD_Arm_Full [MxN]
%   Velocity Trajectory 
% qDD_Arm_Full [MxN]
%   Acceleration Trajectory

% Moritz Schappler, schappler@irt.uni-hannover.de, 2016-02
% (c) Institut für Regelungstechnik, Universität Hannover

function [q_Arm_Full,qD_Arm_Full,qDD_Arm_Full,w_t]= traj_trapez2(q_Start, q_End, vmax, T2, T3, T_Abt)
%% Init
% Coder Information
assert(isa(q_Start,'double') && isreal(q_Start) && size(q_Start,1) == 1, ...
  'traj_trapez2: q_Start has to be 1xN double.');
assert(isa(q_End,'double') && isreal(q_End) && size(q_End,1) == 1, ...
  'traj_trapez2: q_End has to be 1xN double.');
assert(isa(T_Abt,'double') && isreal(T_Abt) && all(size(T_Abt) == [1,1]), ...
  'traj_trapez2: T_Abt has to be 1x1 double.');


N = length(q_Start);


if length(vmax) == 1
  vmax = repmat(vmax, 1,N);
elseif length(vmax) ~= N
  error('traj_trapez2: Maximalgeschwindigkeit vmax muss [1xN] sein');
end
if length(T2) == 1
  T2 = repmat(T2, 1,N);
elseif length(T2) ~= N
  error('traj_trapez2: Verschliffzeit T2 muss [1xN] sein');
end
if length(T3) == 1
  T3 = repmat(T3, 1,N);
elseif length(T3) ~= N
  error('traj_trapez2: Verschliffzeit T3 muss [1xN] sein');
end

% calculate trapezoid trajectory
[~, ~, w_t, w_Q] = Trapez_nAbl_Multi([q_Start; zeros(2,N)], ...
[q_End; zeros(2,N)], [q_End; vmax; vmax./T2; vmax./(T2.*T3)],T_Abt, false);

% Ausgabegrößen zuweisen
q_Arm_Full = NaN(length(w_t), N);
qD_Arm_Full = NaN(length(w_t), N);
qDD_Arm_Full = NaN(length(w_t), N);  

for i = 1:N
  q_Arm_Full(:,i) = w_Q(:,1,i);
  qD_Arm_Full(:,i) = w_Q(:,2,i);
  qDD_Arm_Full(:,i) = w_Q(:,3,i);
end  
