% Beispielskript zum Erstellen einer Orientierungstrajektorie in Form eines Kegels
% Ruft die Funktion traj_rot_cone_angaxis() auf.

% Moritz Schappler, schappler@irt.uni-hannover.de, 2017-03
% (c) Institut für Regelungstechnik, Universität Hannover

clear
clc
%% Funktionsaufruf
traj_settings = struct( ...
  'delta_theta', 30*pi/180, ...
  'arclength', 2*pi, ...
  'cone_axis', [0;0;1], ...
  'start_axis', [1;0;0], ...
  'omega_max', [0.5 1 0.3],...
  'T2', ones(1,3)*50e-3, ...
  'T3', ones(1,3)*5e-3, ...
  'Ts', 1e-3, ...
  'T_pause_hin', 0.2, ...
  'T_pause_rueck', 0.3);
[traj_ges] = traj_rot_cone_angaxis(traj_settings);

%% Test: Integration der Winkelgeschwindigkeit
% Ergebnis: Die Trajektorie muss wieder am Ursprungsort ankommen
% omega_ges = traj_ges.omega;
R_t0 = traj_ges.R(:,:,1);

[R_int_test, t_out] = angvel_int_sl(traj_ges.t, traj_ges.omega, 1e-3, R_t0);
sl = struct('t', t_out, 'R_int', R_int_test);
clear t_out R_int_test

% Vergleiche Aufintegrierte Winkelgeschwindigkeit mit der vorgegebenen
% Rotationsmatrix
for i = 1:length(traj_ges.t)
  [~,I_out] = min(abs(sl.t-traj_ges.t(i)));
  % Vergleiche korrekte Integration
  R_int_i = sl.R_int(:,:,I_out);
  test3=(R_int_i\traj_ges.R(:,:,i)) - eye(3);
  if any(abs(test3(:))>1e-2)
    error('t=%fs. Winkelgeschwindigkeit stimmt nicht mit Rotationsmatrix überein.', traj_ges.t(i));
  end
end

%% Zeichne Vergleich
figure(1);clf;
for ix = 1:3
  for iy = 1:3
    subplot(3,3,sprc2no(3,3,ix,iy));hold on;
    plot(sl.t, squeeze(sl.R_int(ix,iy,:)));
    plot(traj_ges.t, squeeze(traj_ges.R(ix,iy,:)), '--');
    if ix==3 && iy==3
      legend({'int,SL','Soll'});
    end
    ylabel(sprintf('R_{%d,%d}',ix,iy));
    grid on;
  end
end
linkxaxes
% linkxaxes

figure(2);clf;
for i = 1:3
  subplot(3,1,i);hold on;
  plot(traj_ges.t, traj_ges.omega(:,i));
  ylabel(sprintf('\\omega_{%s}',char(119+i)));grid on;
  if i == 3
    legend({'soll'});
  end
end
linkxaxes

% Verlauf der Einheitsvektoren
scal = 0.05;
Ex = NaN(length(traj_ges.t),3);Ey=Ex;Ez=Ex;
for i = 1:length(traj_ges.t)
  Ex(i,:) = traj_ges.R(:,1,i)*scal;
  Ey(i,:) = traj_ges.R(:,2,i)*scal;
  Ez(i,:) = traj_ges.R(:,3,i)*scal;
end
figure(3);clf;hold on
% Verlauf
plot3(Ex(:,1),Ex(:,2),Ex(:,3),'r-');
plot3(Ey(:,1),Ey(:,2),Ey(:,3),'g-');
plot3(Ez(:,1),Ez(:,2),Ez(:,3),'b-');
% Start
plot3(Ex(1,1),Ex(1,2),Ex(1,3),'r^');
plot3(Ey(1,1),Ey(1,2),Ey(1,3),'g^');
plot3(Ez(1,1),Ez(1,2),Ez(1,3),'b^');
% Ende
plot3(Ex(end,1),Ex(end,2),Ex(end,3),'rv');
plot3(Ey(end,1),Ey(end,2),Ey(end,3),'gv');
plot3(Ez(end,1),Ez(end,2),Ez(end,3),'bv');
grid on;
xlabel('x');ylabel('y');zlabel('z');
view(3);