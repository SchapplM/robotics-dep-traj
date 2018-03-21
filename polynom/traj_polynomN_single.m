% Berechnung einer Polynom-Trajektorie für die n. Ableitung einer Größe.
%
% Eingangsgrößen
% z0        Anfangswert für die Größe z(1) und alle Ableitungen (z(2), ...)
% zT        Endwert ...
% t0        Anfangszeit
% zmax      Maximalwert für die Größe und alle Ableitungen
% T_Abt     Abtastzeit der Trajektorie
% T_Dauer   Geforderte Dauer des Polynoms. Überschreibt Dauer aufgrund von
%           Maximalwerten zmax
%
% Ausgabe
% w_z       Zwischenwerte
% w_t       Zwischenzeiten
%

% MA Moritz Schappler, schapplm@stud.uni-hannover.de, 2014-01
% Institut für mechatronische Systeme, Universität Hannover
% Betreuer: Daniel Beckmann, Daniel.Beckmann@imes.uni-hannover.de


function [w_z, w_t] = traj_polynomN_single(z0, zT, t0, zmax, T_Abt, T_Dauer)
nz = size(z0, 1);

if nargin < 6
    T_Dauer = [];
end

% Polynomkoeffizienten (immer gleich) und Trajektoriendauer (teilw.
% unterschiedlich) berechnen.
[~, ~, PolyKoeff_norm, T] = traj_polynom_bernstein(z0(1:nz), zT(1:nz), zmax(1:nz));
% Polynomdauer an Abtastzeit anpassen (vergrößern)
T = ceil(T/T_Abt)*T_Abt;

% Feste Polynomdauer vorgeben (z.B. zum Vergleich mit anderen Polynomen)
if ~isempty(T_Dauer)
    if T_Dauer < T
        warning(['Polynomdauer von %1.4f s gefordert. Aufgrund von ', ...
            'Beschränkungen aber mindestens %1.4 s erforderlich. Nehme ', ...
            'geforderte Dauer von %1.4f s'], T_Dauer, T, T_Dauer);
    else
        warning('Nehme längere Polynomdauer als aufgrund von Beschränkungen benötigt.');
    end
    T = T_Dauer;
end

% Werte berechnen
w_t = (t0:T_Abt:t0+T)';

% Schrittweite
h = zT(1)-z0(1);

% Werte für relative Zeit berechnen
w_tau = (w_t-t0)/T;

% Polynom-Werte für jede Ableitung berechnen
w_z = NaN(length(w_t), nz+1);
for iz=1:nz+1
    w_z(:, iz) = h/T^(iz-1) * polyval(PolyKoeff_norm(iz, :), w_tau);
end
% Anfangs-Weg aufaddieren
w_z(:, 1) = w_z(:, 1) + z0(1);

% End-Zustand überschreiben, damit Rundungsfehler nicht Geschwindigkeit etc
% größer Null lässt
w_z(end, 1:nz) = zT';
