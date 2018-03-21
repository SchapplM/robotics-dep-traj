% Berechne Polynom-Trajektorie zwischen zwei Punkten mit Randbedingungen
%
% Die Berechnung erfolgt mit Bernstein-Polynomen. Diese sind numerisch
% robuster und in der Berechnung schneller als die Berechnung des Polynoms
% mit einem linearen Gleichungssystem.
%
% Eingaben:
% z0        Anfangswerte (beliebige Geschwindigkeiten etc.)
%           z0(1): Position
%           z0(2): Geschwindigkeit
%           z0(2): Beschleunigung (usw.)
% zT        Endwerte (genauso)
% zmax      Maximalwerte für Weg, Geschwindigkeit etc. (nur positiv)
% T_Abt     Abtastzeit, mit der die Zeitwerte w_t ausgegeben werden
%
% Ausgaben:
% ew_t      Zeit-Eckwerte für Minimal- und Maximalstellen der Polynome
%           aller Ableitungen
% ew_z      Werte aller Ableitungen zu den Eckzeiten ew_t
% A         Normierte Polynomkoeffizienten für den Weg und alle Ableitungen
%           normiert: Auf relative Zeit tau bezogen (0...1)
%           Höchster Koeffizient zuerst (zur Eingabe in polyval)
% T         Trajektoriendauer [s] für das Polynom (minimal)
% BegrAblNr
%           Nummer des Zustands z, der die Dauer begrenzt
% w_t       Zeitstützstellen des kompletten Polynoms (Abtastzeit T_Abt)
% w_z       Werte des Polynoms und seiner Zeitableitungen zu w_t
%
% Quelle:
% [1] Biagiotti, Melchiorri: Trajectory Planning for Automatic Machines
% and Robots (2006)

% MA Moritz Schappler, schapplm@stud.uni-hannover.de, 2013-12
% Institut für mechatronische Systeme, Universität Hannover
% Betreuer: Daniel Beckmann, Daniel.Beckmann@imes.uni-hannover.de

function [ew_t, ew_z, A, T, BegrAblNr, w_t, w_z] = ...
    traj_polynom_bernstein(z0, zT, zmax, T_Abt)
%% Initialisierung

nz = length(z0); % Anzahl der Randbedingungen (nciz = nz, ncf = nz)
n = 2*nz-1; % Grad des Polynoms (n+1 Koeffizienten)


% z0 und zT um eine Ordnung erweitern, damit Bernstein-Polynom gebildet
% werden kann
z0 = [z0; 0]; zT = [zT; 0];
h = zT(1)-z0(1);
%% Berechnung der Koeffizienten des Bernstein-Basis-Polynoms
if all([z0(2:end); zT(2:end)] == 0)
    % Null-Randbedingungen, einfachere Berechnung
    % [1], S. 36
    p = [zeros(nz,1); ones(nz,1)];
else
    % Gleichungssystem für Polynom-Koeffizienten erzeugen
    % Start [1], S.32, u.a. (2.12)
    M0 = zeros(nz, nz);
    M0(1,1)=1;
    for k = 1:nz-1
        for i = 0:k
            M0(k+1, i+1) = nchoosek(k, i) * (-1)^(k+i);
        end
    end
    b0 = zeros(nz, 1);
    for k = 1:nz-1
        b0(k+1) = z0(k+1) / factorial(n)/factorial(n-k);
    end
    p0 = M0\b0; % Koeffizienten berechnen.  [1], S.33 (2.14)
    % Ende [1], S.33, u.a. (2.13)
    M1 = zeros(nz, nz);
    M1(1,1)=1;
    for k = 1:nz-1
        for i = 0:k
            M1(k+1, i+1) = nchoosek(k, i) * (-1)^(i);
        end
    end
    b1 = ones(nz, 1);
    for k = 1:nz-1
        b1(k+1) = zT(k+1) / factorial(n)/factorial(n-k);
    end
    % Koeffizienten berechnen. Achtung: p1 enthält die Koeffizienten in
    % umgekehrter Reihenfolge. [1], S.33 (2.15)
    p1 = M1\b1;
    p = [p0; flipud(p1)]; % [1], S.34 oben
end

%% Berechnung der Koeffizienten des normierten Polynoms
% Vergleiche Polynomkoeffizienten zum Testen mit [1], S. 37, Tabelle 2.1
% [1], S.32, Formel (2.10)
a = zeros(1, n+1);
for j = 0:n
    K = factorial(n)/factorial(n-j);
    for i = 0:j
        a(j+1) = a(j+1) + K*(-1)^(i+j)/(factorial(i)*factorial(j-i))*p(i+1);
    end
end


%% Berechne Zahlenwerte des normierten Polynoms und der Ableitungen Eckzeiten
ew_tau = ones(1, nz); % Zeitpunkte von Maximalwerten aller Ableitungen
ew_z = zeros(nz, nz); % Werte der Ableitungen des nicht-normierten Polynoms
T_min = zeros(nz, 1); % Minimalmögliche Zeit aufgrund der zmax

% Berechne Koeffizienten der Polynome der Ableitungen in normierter Zeit tau
A = zeros(nz+1, n+1);
A(1, :) = fliplr(a);
% Bilde Ableitungen der normierten Polynome
% Kontrolle der ersten beiden Ableitungen mit [1], S. 459, Tabelle A.2, A.3
for i = 2:nz+1
    ablcoeff = polyder(A(i-1, :)); %polyder: Höchste Ableitung zuerst
    A(i, :) = [zeros(1, n+1-length(ablcoeff)), ablcoeff];
    if i>nz, continue;end
    % Maximalwert und Maximalstelle dieser Ableitung
   [pmax, pmin, maxstelle, minstelle] = polymax(A(i, :), 0, 1);
   [wmax, im] = max(abs([pmin, pmax])); % Maximalwert des Polynoms
   if im == 1
       ew_tau(i) = minstelle;
   else
       ew_tau(i) = maxstelle;
   end
   % Berechne Minimalzeit. Umformung von [1], S. 230, (5.6)
   T_min(i) = (abs(h)/(zmax(i)) * ...
        norm(wmax))^(1/(i-1)); % Minimalzeit aufgrund der Begrenzung

end

% Minimal mögliche Zeit für die Trajektorie berechnen
[T, BegrAblNr] = max([T_min; 1e-6]);

ew_t = ew_tau*T;
% Berechne Werte des Polynoms für alle gefundenen Stellen
for i = 1:nz
    ew_z(i, :) = h/T^(i-1)*polyval(A(i, :), ew_t/T);
end
ew_z(1, :) = ew_z(1, :)+z0(1);
% Rückgabe der Polynomkoeffizienten der höchsten Ableitung
% PolyKoeffnz_tau = A(end, :)*h;
% A=A;

% Berechne die Trajektorie
if nargout > 5
    w_t = (0:T_Abt:T)';
    w_z = NaN(length(w_t), nz+1);
    for i = 1:nz+1
        w_z(:, i) = h/T^(i-1)*polyval(A(i, :), w_t/T);
    end
    w_z(:, 1) = w_z(:, 1) + z0(1);
end

%% Debug: Polynom zeichnen
return; %#ok<*UNRCH>
figure(1);clf;
t_w = linspace(0, T, 1000);
z_w = zeros(nz, 1000);
clf;
for i = 1:nz
    z_w(i, :) = h/T^(i-1)*polyval(fliplr(A(i, :)), t_w/T);
    if i == 1
        z_w(1, :) = z_w(1, :) + z0(1);
    end
    subplot(nz, 1, i);
    hold on;
    plot(t_w', z_w(i, :)');
    plot(ew_t, ew_z(i, :), 'ro');
    plot([0; T], [zmax(i); zmax(i)], 'r--');
end

