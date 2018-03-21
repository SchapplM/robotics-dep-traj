% Beispielskript für die Erstellung von Polynom-Trajektorien

% Moritz Schappler, schappler@irt.uni-hannover.de, 2018-03
% (C) Institut für mechatronische Systeme, Leibniz Universität Hannover

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
jDDD_max = jDD_max/2e-3;

% Anzahl der Randbedingungen
% Das Wegprofil ist n-1 mal stetig differenzierbar.
for n = 2:6


  % Eingabedaten füllen
  z0 = [s_start; zeros(n-1,1)];
  zT = [s_end; zeros(n-1,1)];
  zmax = [s_end; v_max; a_max; j_max; jD_max; jDD_max; jDDD_max];
  zmax = zmax(1:n+1);

  % Abtastrate der Trajektorie
  T_Abt = 1e-3;

  % Funktionsaufruf
  [w_z, w_t] = traj_polynomN_single(z0, zT, 0, zmax, T_Abt);

  % Ergebnisse zeichnen
  figure(n);clf;axhdl = NaN(n+1,1);
  for i = 1:n+1
    axhdl(i) = subplot(n+1,1,i);
    plot(w_t, w_z(:, i));

    grid on;
    hold all;
%     plot(w_t([1,end]), -zmax(i)*[1;1], 'r--');
%     plot(w_t([1,end]),  zmax(i)*[1;1], 'r--');
    ylabel(sprintf('D%d s(t)', i-1));
    
    if i == 1
      title(sprintf('Polynomprofil %d mal stetig Diff.bar', n-1));
      set(n, 'Name', sprintf('Stetigkeit%d', n-1), 'NumberTitle', 'off');
    end
  end
  linkxaxes
  remove_inner_labels(axhdl, 1);
end

dockall