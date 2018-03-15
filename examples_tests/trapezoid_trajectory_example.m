% Beispielskript für die Erstellung von Trapez-Trajektorien

% Moritz Schappler, schappler@irt.uni-hannover.de, 2015-12
% (c) Institut für Regelungstechnik, Leibniz Universität Hannover

%% Beispiel für einachsige Bewegung
% Start- und Endposition
s_start = 0;
s_end = 1;

% Maximalwerte für Zeitableitungen
v_max = 2;
a_max = v_max / 200e-3;
j_max = a_max / 50e-3;
jD_max = j_max/20e-3;
jDD_max = jD_max/5e-3;

% Anzahl der Randbedingungen
% Das Wegprofil ist n-1 mal stetig differenzierbar.
n = 5;

% Eingabedaten füllen
z0 = [s_start; zeros(n-1,1)];
zT = [s_end; zeros(n-1,1)];
zmax = [s_end; v_max; a_max; j_max; jD_max; jDD_max];
zmax = zmax(1:n+1);

% Abtastrate der Trajektorie
T_Abt = 1e-3;

% Funktionsaufruf
[ew_t, ew_z, w_z, w_t] = Trapez_nAbl(z0, zT, 0, zmax, T_Abt, false);

% Ergebnisse zeichnen
figure(1);clf;
for i = 1:n+1
  subplot(n+1,1,i);
  if i < n+1
    plot(w_t, w_z(:, i));
  else
    stairs(w_t, w_z(:, i));
  end
  grid on;
  hold all;
  plot(ew_t, ew_z(:, i), 'o');
  
  plot([ew_t(1); ew_t(end)], zmax(i)*[1;1], 'r--');
  plot([ew_t(1); ew_t(end)], -zmax(i)*[1;1], 'r--');
  
  ylabel(sprintf('D%d s(t)', i-1));
end
linkxaxes