% Bestimme ein Geschwindigkeitsprofil für eine abgebremste Nullraumbewegung
% Grundlage ist eine Rast-zu-Rast-Bewegung für eine Aufgabenbewegung. Die 
% Methode ist hauptsächlich bei Trapez-Profilen mit Abbremsphase am Ende
% sinnvoll. Vor der Verzögerungsphase der Aufgabenbewegung wird die Null-
% raumbewegung auf Null gebremst. Wird dies nicht gemacht, ist bei Rast der
% Aufgabenbewegung immer noch eine Nullraumgeschwindigkeit vorhanden.
% Diese Funktion ist daher nur bei Trajektorien für redundante Mehrachs-
% systeme sinnvoll (z.B. 6FG-Roboter für 5FG-Aufgaben).
% 
% Eingabe:
% t [NT x 1]
%   Zeit-Stützstellen der Trajektorie
% I_rest [N+1 x 1]
%   Indizes für Rast-Punkte in `t`. Die Trajektorie besteht aus `N` PTP-
%   Bewegungen, zwischen denen eine Rast stattfindet.
% Tv
%   Abbauzeit der Geschwindigkeit in Trapezprofil (Dauer der Beschleuni-
%   gungsphase)
% Tv_ns
%   Abbauzeit der Nullraumgeschwindigkeit
% Ts
%   Abtastzeit der Trajektorie (Annahme: äquidistante Schritte)
% 
% Ausgabe:
% nullspace_maxvel_interp [2xM]
%   Interpolations-Stützstellen (1. Zeile) und Skalierungswerte (2. Zeile)
%   für die Nullraumgeschwindigkeit. Die Skalierung ist normiert auf 0-1.
%   Zwischen den Stützstellen kann bspw. mit interp1 linear interpoliert
%   werden.
% 
% Siehe auch: traj_trapez2_multipoint.m, invkin_traj-Funktionen (Robotik)

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2022-02
% (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

function nullspace_maxvel_interp = nullspace_maxvel_from_tasktraj(t, ...
  I_rest, Tv, Tv_ns, Ts)
% Zeit-Stützstellen für die Rastpunkte
t_rest = t(I_rest);

% Berechne die Zeitpunkte der einzelnen Phasen der Bewegung aus den Rast-
% punkten. Zusätzliche Prüfung bei sehr kurzen Bewegungsphasen, damit die
% Beschleunigungs- und Bremsphasen nicht unlogisch überlappen.
% Der letzte Rastpunkt in t_rest entspricht dem Ende. Indizierung entspr.
% Zeitpunkte, an denen die Beschleunigung der Aufgaben-PTP-Bewegung endet.
t_acc_done = t_rest(1:end-1) + Tv; % Verschliffzeit von traj_trapez2_multipoint
% Zeitpunkte, an denen die Verzögerung (Aufgaben-PTP) beginnt.
t_dec_start = t_rest(2:end) - Tv;
% Start der Nullraum-Beschleunigung: Sofort, gemeinsam mit Aufgabe
t_ns_acc_start = t_rest(1:end-1); %Alternative für Abwarten der Aufgaben-Beschleunigung: t_acc_done;
% Ende der Rampe für Nullraum-Beschleunigung: Sofort die max. Nullraum-
% geschwindigkeit erlauben.
t_ns_acc_finish = t_rest(1:end-1)+0.5e-3;
% Beginn der Nullraum-Verzögerungsphase.
t_ns_dec_start = max([t_dec_start-Tv_ns, t_ns_acc_finish+Ts/2], [], 2);
% Ende der Nullraum-Verzögerungsphase
t_ns_dec_finish = min([t_ns_dec_start + Tv_ns, [t_ns_acc_start(2:end)-Ts/2;inf]], [], 2);
% Falls der letzte Wert nach der Trajektorie liegt, setze auf Ende
t_ns_dec_finish(t_ns_dec_finish>t(end)) = t(end);
% Zeitdauer der Verzögerungsphasen (Debug)
% T_dec = t_ns_dec_finish-t_ns_dec_start;
% Prüfe die Plausibilität (Reihenfolge der einzelnen Schritte). Wird durch
% min-/max-Befehle oben eigentlich sichergestellt.
assert(all(t_ns_dec_finish(1:end-1) < t_ns_acc_start(2:end)), ['Start ', ...
  'der Nullraumbewegung muss nach Ende der vorherigen sein']);
assert(all(t_ns_acc_start(1:end) < t_ns_acc_finish(1:end)), ['Logik-', ...
  'Fehler bei Nullraum-Begrenzung (t_ns_acc_start)']);
assert(all(t_ns_acc_finish(1:end) < t_ns_dec_start(1:end)), ['Logik-', ...
  'Fehler bei Nullraum-Begrenzung (t_ns_acc_finish)']);
assert(all(t_ns_dec_start(1:end) < t_ns_dec_finish(1:end)), ['Logik-', ...
  'Fehler bei Nullraum-Begrenzung (t_ns_dec_start)']);
nullspace_maxvel_interp = [ ... % zunächst nach Kategorie sortiert sammeln
  [t_ns_acc_start'; zeros(1,length(t_ns_acc_start))], ...
  [t_ns_acc_finish'; ones(1,length(t_ns_acc_finish))], ...
  [t_ns_dec_start'; ones(1,length(t_ns_dec_start))], ...
  [t_ns_dec_finish'; zeros(1,length(t_ns_dec_finish))], ...
  ];
% Letzten Zeitschritt hinzufügen, damit dieser auch definiert ist
if ~any(nullspace_maxvel_interp(1,:)==t(end))
  nullspace_maxvel_interp = [nullspace_maxvel_interp, [t(end); 0]];
end
% anschließend sortieren und nochmals prüfen
assert(all(nullspace_maxvel_interp(1,:)<=t(end)), ['Zeitschritte dürfen ', ...
  'nicht weiter als Trajektorie gehen'])
[~,I_sort] = sort(nullspace_maxvel_interp(1,:));
nullspace_maxvel_interp = nullspace_maxvel_interp(:,I_sort);
assert(all(diff(nullspace_maxvel_interp(1,:))>0), ['Zeit-Stützstellen ', ...
  'müssen aufsteigend sein']);
assert(size(nullspace_maxvel_interp,2)==length(unique( ...
  nullspace_maxvel_interp(1,:))), 'Zeit-Stützstellen nicht eindeutig');
