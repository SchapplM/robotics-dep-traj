% Berechne beliebige Zwischenwerte für die Trapez-Trajektorie der n. Abl.
% 
% Eingabe
% ew_t      Eckwerte der Zeiten
% ew_z      Eckwerte der Größe z(1) und ihren Ableitungen z(2:end)
% w_t       Zwischenwerte der Zeiten
% 
% Ausgabe
% 
% w_z       Zwischenwerte der Größe z(1) und ihren Ableitungen z(2:end) zu
%           den Zeiten w_t
% ew_z      Mit Integration neu berechnete Eckwerte zu den Zeiten ew_t,
%           ausgehend von den Zeiten und den Werten für die letzte Ableitung
% 
% 
% MA Moritz Schappler, schapplm@stud.uni-hannover.de, 2014-01
% Institut für mechatronische Systeme, Universität Hannover
% Betreuer: Daniel Beckmann, Daniel.Beckmann@imes.uni-hannover.de



function [w_z, ew_z] = Trapez_nAbl_Werte(ew_t, ew_z, w_t)


%% Initialisierung 

%#codegen
assert(isa(ew_z, 'double') && isreal(ew_z) && size(ew_z, 2) > 1 && size(ew_z, 1) > 1);
assert(isa(ew_t, 'double') && isreal(ew_t) && size(ew_t, 1) > 1);
assert(isa(w_t, 'double') && isreal(w_t) && size(w_t, 1) > 0);

if nargin == 2
    w_t = []; % Keine Berechnung von Zwischenwerten, nur Eckwerte mit Integration berechnen
end

coder.extrinsic('polyint')

w_z = NaN(length(w_t), size(ew_z, 2));
nz = size(ew_z, 2)-1;


%% Zwischenwerte iterativ durch Integration berechnen
for it = 1:length(ew_t)-1
    % Indizes der Zeitpunkte in diesem Segment
    % Werte der letzten <=-Bedingung werden in nächster Iteration überschrieben
    II = (ew_t(it) <= w_t) & (w_t <= ew_t(it+1)); 
    if ~any(II) && ~any(isnan(ew_z(it+1, :)))
        continue % keine Zeitwerte in diesem Bereich. Keine Integration notwendig
    end
    
    % Polynomkoeffizienten des Polynoms auf diesem Segment berechnen
    polykoeff = zeros(nz+1, nz+1); % Polynomkoeffizienten: höchster zuerst
    polykoeff(nz+1, nz+1) = ew_z(it, nz+1); % letzte Ableitung ist konstant
    % Integration von der höchsten Wegableitung z(nz+1) zum Weg z(1)
    for iz = nz:-1:1
        polykoeff(iz, :) = polyint(polykoeff(iz+1, 2:nz+1), ew_z(it, iz));
    end
    
    % Zwischenwerte berechnen
    for iz = 1:nz+1
        w_z(II, iz) = polyval(polykoeff(iz, :), w_t(II)-ew_t(it))';
    end
    
    % Eckwerte des folgenden Zeitschritts berechnen (falls nicht gesetzt)
    for iz = 1:nz+1
        if isnan(ew_z(it+1, iz))
            % nächsten Eckwert durch Integration im aktuellen Segment
            % berechnen
            ew_z(it+1, iz) = polyval(polykoeff(iz, :), ew_t(it+1)-ew_t(it));
        end
    end 
end

%% Zeiten außerhalb der Trajektorie abfangen
II = w_t < ew_t(1); % vor Beginn
w_z(II, :) = repmat(ew_z(1, :), sum(II), 1);
w_z(II, end) = 0;

II = w_t >= ew_t(end); % nach Ende
w_z(II, :) = repmat(ew_z(end, :), sum(II), 1);
w_z(II, end) = 0;