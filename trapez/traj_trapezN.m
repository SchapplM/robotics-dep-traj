% Trapez-Trajektorien für mehrere Achsen gleichzeitig
% Die Synchronisierung erfolgt auf die langsamste Achse
% 
% * n Achsen
% * m vorgegegebene Randbedingungen: Positionen und Zeitableitungen
% * m+1 vorzugebene Grenzwerte
% * Die Stetigkeit des Positionsverlaufs der Trajektorie ist damit m
% 
% Eingabe:
% q0 [mxn]
%   Anfangswerte: m Zeitableitungen, n Gelenkwinkel
% qT [mxn]
%   Endwerte: m Zeitableitungen, n Gelenkwinkel
% qmax [(m+1)xn]
%   Endwerte: m Zeitableitungen, n Gelenkwinkel
% T_Abt     
%   Abtastzeit des Trapez-Profils. Die Eckzeiten sind ganzzahlige vielfache
%   davon
% debug     Zusätzliche Informationen ausgeben
% 
% Ausgabe:
% ew_t [px1]
%   Eckwerte der Zeit. p=2^m Zeitschritte
% ew_q [px(m+1)xn]
%   Eckwerte für Position und alle Zeitableitungen für alle n Achsen
% 
% TODO:
% * Modus mit alle Achsen schnellstmöglich hinzufügen
% * assert-Befehle für Dimensionen hinzufügen.

% Moritz Schappler, schappler@irt.uni-hannover.de, 2015-02
% (c) Institut für Regelungstechnik, Universität Hannover

function [ew_t, ew_Q, w_t, w_Q] = traj_trapezN(q0, qT, qmax, T_Abt, debug)

%% Init
n = size(q0,2); % Anzahl der Achsen
m = size(q0,1); % Stetigkeitsordnung

net = 2^m; % Anzahl der Eckzeiten
tmp_t = NaN(net,n); % Zwischenspeichern aller Eckzeiten für die Einzelachstrajektorien

ew_Q = NaN(net,m+1,n);
%% Berechne einzelne Trapeztrajektorie für jedes Gelenk
for i = 1:n
  % Berechne Trapezprofil für Achse i
  [ew_t_i, ew_z_i] = traj_trapezN_single(q0(:,i), qT(:,i), 0, qmax(:,i), T_Abt, debug);
  
  % Speichern
  tmp_t(:,i) = ew_t_i;
  ew_Q(:,:,i) = ew_z_i;
end

%% Skaliere die Zeit auf die langsamste Trajektorie
% langsamste Trajektorie finden
[t_langsamste, I_langsamste] = max(tmp_t(end,:));

% Alle Trajektorien skalieren
for i = 1:n
  % Skalierungsfaktor zur verringerung der Geschwindigkeiten
  Skalierung = tmp_t(end,i)/t_langsamste;
  
  % alle Zeitableitungen verlangsamen
  ew_Q(:,2:end,i) = ew_Q(:,2:end,i) * Skalierung;
end

% benutze die langsamste Zeit für die gesamte Trajektorie
ew_t = tmp_t(:,I_langsamste);

% Ausgabe der feinaufgelösten Trajektorie
w_t = (0:T_Abt:t_langsamste)';

% Neuberechnung der mit T_Abt aufgelösten Trajektorie
w_Q = NaN(length(w_t), m+1, n);

% Trajektorie für alle Achsen neu berechnen
for i = 1:n
  for iter = 1:3
    % Eckwerte mit neu Berechnen, da bei der Skalierung anscheinend nicht ein
    % konsistenter Positionsverlauf entsteht
    % zur Neuberechnung Eckwerte NaN setzen
    % Ersten Zeitpunkt behalten
    ew_Q(2:end,1:end-1,i) = NaN;

    % Berechnung der Eckwerte und feinaufgelösten Trajektorien
    [w_Q(:,:,i), ew_Q(:,:,i)] = traj_trapezN_values(ew_t, ew_Q(:,:,i), w_t);

    % Prüfe die Endposition der Trajektorie, kann sich durch Integration
    % der Numerik-Fehler verändern
    if abs(ew_Q(end,1,i)-qT(1,i))>1e-6
      if debug
        fprintf(['Achse %d: Der zurückgelegte Weg ist falsch. Endposition Soll: %f, Ist: %f. ', ...
          'Verändere letzten Grenzwert.\n'], i, qT(1,i), ew_Q(end,1,i));
      end
      % Die Endposition weicht ab. Neu Berechnen. Letzte Zeitableitung
      % entsprechend verändern.

      % Soll-Strecke geteilt durch Ist-Strecke
      Skalierung = abs( (qT(1,i)-ew_Q(1,1,i)) / (ew_Q(end,1,i)-ew_Q(1,1,i)) );
      ew_Q(:,end,i) = Skalierung * ew_Q(:,end,i);
      if debug
        fprintf('Achse %d: Grenzwert Abl. %d neu gesetzt: Skalierungsfaktor %f\n', ...
          i, size(ew_Q(end,:,i),1), Skalierung);
      end
    end
  end
end