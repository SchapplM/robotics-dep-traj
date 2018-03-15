% Rotationstrajektorie in Form eines Kegels mit Achse-Winkel-Notation
% Der Kegel wird ausgehend vom Koordinatenursprung durch die
% Einheitsvektoren des Koordinatensystems beschrieben
% 
% Eingabe:
% traj_settings
%   Struktur mit Parametern der zu berechnenden Trajektorie
%   (Standardeinstellungen siehe unten)
%   .delta_theta [1x1] [rad]
%     Winkel des Kegels gegenüber der Drehachse. Die erste Drehung
%     im "start_axis" dreht um diesen Winkel
%   .arclength [1x1] [rad]
%     Bogenlänge des Kegels (2*pi ist ein vollständiger Kegel)
%   .cone_axis [3x1]
%     Achse, um die der Kegel erzeugt wird
%     Bezogen auf KS 0
%   .start_axis [3x1]
%     Achse, um die die erste Verkippung zum Starten des Kegels erzeugt wird
%     Sollte Senkrecht auf der Drehachse des Kegels stehen
%     Bezogen auf KS 0
%   .omega_max [1x3] [rad/s]
%     Maximale Winkelgeschwindigkeiten in den Phasen Anfahrt, Kegel, Rückfahrt
%   .T2 [1x3] [s]
%     Verschliffzeit der Geschwindigkeit (in den drei Phasen). Siehe traj_trapez2.m
%   .T3 [1x3] [s]
%     Verschliffzeit der Beschleunigung (in den drei Phasen). Siehe traj_trapez2.m
%   .Ts [s]
%     Abtastrate der Trajektorie
%   .T_pause_hin [s]
%     Pause nach der Hinfahrt
%   .T_pause_rueck [s]
%     Pause vor der Rückfahrt
% 
% Ausgabe:
% traj_ges
%   Struktur mit der Orientierungs-Trajektorie. Felder
%     .t      Zeit [s]
%     .R      Rotationsmatrix i
%     .omega  Winkelgeschwindigkeit [rad/s]
%     .omegaD noch nicht belegt [rad/s^2]
% traj_settings
%   Struktur mit den Einstellungen. Gibt die Standard-Einstellungen zurück,
%   falls die Funktion mit traj_rot_cone_angaxis([]) aufgerufen wird.
% 
% Koordinatensysteme
% 0: Ursprüngliches KS vor Trajektorie
% 1: Nach Anfahrt
% 2: Nach Kegel-Trajektorie
% 3: Nach Rückfahrt

% Moritz Schappler, schappler@irt.uni-hannover.de, 2017-03
% (c) Institut für Regelungstechnik, Universität Hannover

function [traj_ges, traj_settings] = traj_rot_cone_angaxis(traj_settings)

%% Eingabe verarbeiten
if isempty(traj_settings)
  traj_settings = struct( ...
    'delta_theta', 10*pi/180, ... % Winkel des Kegels gegenüber der Drehachse
    'arclength', 2*pi, ... % Bogenlänge des Kegels (2*pi ist ein vollständiger Kegel)
    'cone_axis', [0;0;1], ... % Achse, um die der Kegel erzeugt wird
    'start_axis', [1;0;0], ... % Achse, um die die erste Verkippung zum Starten des Kegels erzeugt wird
    'omega_max', ones(1,3)*1.0, ...
    'T2', ones(1,3)*50e-3, ...
    'T3', ones(1,3)*5e-3, ...
    'Ts', 1e-3, ...
    'T_pause_hin', 0.5, ...
    'T_pause_rueck', 0.5);
  return
end

Ts = traj_settings.Ts;

%% Anfahrt zum Startpunkt
traj_hin = struct('t', [], 'R', [], 'omega', [], 'omegaD', []);
vec1 = traj_settings.start_axis;
if ~any(isnan(vec1))
  % Daten für Trajektorie
  dv = traj_settings.delta_theta;
  vmax = traj_settings.omega_max(1);
  T2 = traj_settings.T2(1);
  T3 = traj_settings.T3(1);

  % Trajektorie: Winkelbeschleunigungs-Trapez
  [theta1,theta1D,~,w_t1]= traj_trapez2(0, dv, vmax, T2, T3, Ts);

  % Anfahrt-Trajektorie als Struktur
  traj_hin.t = w_t1;
  traj_hin.R = NaN(3,3,length(w_t1));
  traj_hin.omega = NaN(length(w_t1),3);
  traj_hin.omegaD = NaN(length(w_t1),3);
  for i = 1:length(w_t1)
    % Interpoliere Rotationsmatrix aus Achse-Winkel-Darstellung
    R0i = angvec2r(theta1(i), vec1);
    traj_hin.R(:,:,i) = R0i;

    % Berechne Winkelgeschwindigkeit
    traj_hin.omega(i,:) = angvecD2omega_sym(vec1, theta1(i), theta1D(i));
  end
  R01 = R0i;
else
  R01 = eye(3);
end


%% Warten bei Übergang nach Hinfahrt
traj_warte_hin = struct('t', [], 'R', [], 'omega', [], 'omegaD', []);
T_pause_hin = traj_settings.T_pause_hin;
if ~isnan(T_pause_hin) && T_pause_hin > 0
  traj_warte_hin.t = (0:Ts:T_pause_hin)';
  traj_warte_hin.R = NaN(3,3,length(traj_warte_hin.t));
  for i = 1:length(traj_warte_hin.t)
    traj_warte_hin.R(:,:,i) = R01;
  end
  traj_warte_hin.omega = zeros(length(traj_warte_hin.t),3);
  traj_warte_hin.omegaD = zeros(length(traj_warte_hin.t),3);
end

%% Kegel-Bewegung
% Trajektorie für Bewegung

dv = traj_settings.arclength;
vmax = traj_settings.omega_max(2);
T2 = traj_settings.T2(2);
T3 = traj_settings.T3(2);

[theta2,theta2D,~,w_t2]= traj_trapez2(0, dv, vmax, T2, T3, Ts);
% Drehachse in KS1 umrechnen
k_1 = R01'*traj_settings.cone_axis;

% Als Struktur
traj_kegel = struct('t', [], 'R', [], 'omega', [], 'omegaD', []);
% Anfahrt-Trajektorie als Struktur
traj_kegel.t = w_t2;
traj_kegel.R = NaN(3,3,length(w_t2));
traj_kegel.omega = NaN(length(w_t2),3);
traj_kegel.omegaD = NaN(length(w_t2),3);
for i = 1:length(w_t2)
  R1i = angvec2r(theta2(i), k_1);
  R0i = R01*R1i;
  traj_kegel.R(:,:,i) = R0i;
  % Berechne Winkelgeschwindigkeit
  omega_1 = angvecD2omega_sym(k_1, theta2(i), theta2D(i));
  % Rotation in korrektes Koordinatensystem
  traj_kegel.omega(i,:) = R01*omega_1;
end
R02 = R0i; % Ist identisch mit R01, wenn einmal im Kreis gefahren wird

%% Warten bei Übergang nach Hinfahrt
traj_warte_rueck = struct('t', [], 'R', [], 'omega', [], 'omegaD', []);
T_pause_rueck = traj_settings.T_pause_rueck;
if ~isnan(T_pause_rueck) && T_pause_rueck > 0
  traj_warte_rueck.t = (0:Ts:T_pause_hin)';
  traj_warte_rueck.R = NaN(3,3,length(traj_warte_rueck.t));
  for i = 1:length(traj_warte_rueck.t)
    traj_warte_rueck.R(:,:,i) = R01;
  end
  traj_warte_rueck.omega = zeros(length(traj_warte_rueck.t),3);
  traj_warte_rueck.omegaD = zeros(length(traj_warte_rueck.t),3);
end

%% Rückfahrt zum Startpunkt
traj_rueck = struct('t', [], 'R', [], 'omega', [], 'omegaD', []);
if ~isnan(T_pause_rueck)
  % Daten für Trajektorie
  dv = traj_settings.delta_theta;
  vmax = traj_settings.omega_max(3);
  T2 = traj_settings.T2(3);
  T3 = traj_settings.T3(3);
  
  % TODO: Das hier funktioniert nur, wenn die Drehung 2*pi war!
  vec3 = traj_settings.start_axis;

  % Trajektorie
  [theta3,theta3D,~,w_t3]= traj_trapez2(dv, 0, vmax, T2, T3, Ts);

  % Anfahrt-Trajektorie als Struktur
  traj_rueck.t = w_t3;
  traj_rueck.R = NaN(3,3,length(w_t3));
  traj_rueck.omega = NaN(length(w_t3),3);
  traj_rueck.omegaD = NaN(length(w_t3),3);
  for i = 1:length(w_t3)
    R0i = angvec2r(theta3(i), vec3);
    traj_rueck.R(:,:,i) = R0i;

    % Berechne Winkelgeschwindigkeit
    traj_rueck.omega(i,:) = angvecD2omega_sym(vec3, theta3(i), theta3D(i));
  end 
end

%% Gesamttrajektorie
traj_ges = timestruct_append(traj_hin, traj_warte_hin, traj_hin.t(end)+Ts);
traj_ges = timestruct_append(traj_ges, traj_kegel, traj_ges.t(end)+Ts);
traj_ges = timestruct_append(traj_ges, traj_warte_rueck, traj_ges.t(end)+Ts);
traj_ges = timestruct_append(traj_ges, traj_rueck, traj_ges.t(end)+Ts);