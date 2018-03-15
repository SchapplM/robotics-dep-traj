% Berechnung einer Trapez-Trajektorie f�r die n. Ableitung einer Gr��e.
% 
% Benutze unterschiedliche Algorithmen, je nach Randbedingungen.
% Bis zur dritten Ordnung kann analytisch gerechnet werden (schneller und
% genauer). Danach mit Faltungsalgorithmus
% 
% Eingabe:
% z0        Anfangswert f�r die Gr��e z(1) und alle Ableitungen (z(2), ...)
% zT        Endwert ...
% t0        Anfangszeit
% zmax      Maximalwert f�r die Gr��e und alle Ableitungen 
% T_Abt     Abtastzeit des Trapez-Profils. Die Eckzeiten sind ganzzahlige
% debug     Zus�tzliche Informationen ausgeben
% 
% Ausgabe:
% ew_t      Eckwerte der Zeiten [Nx1]
% ew_z      Eckwerte der Gr��e z(1) und ihren Ableitungen z(2:end) 
%           [N x (nz+1)]
% w_z       Alle Werte der Gr��e z(1) und ihren Ableitungen z(2:end) 
%           zu den Zwischenzeiten w_t. [n x (nz+1)]
% w_t       Alle Zeitschritte von t0 bis zum Ende der Trajektorie mit dem
%           Abstand T_Abt. [n x 1]


% MA Moritz Schappler, schapplm@stud.uni-hannover.de, 2014-01
% Institut f�r mechatronische Systeme, Universit�t Hannover
% Betreuer: Daniel Beckmann, Daniel.Beckmann@imes.uni-hannover.de

function [ew_t, ew_z, w_z, w_t] = Trapez_nAbl(z0, zT, t0, zmax, T_Abt, debug)
%% Vorbereiten
nz = length(z0);

if length(zmax) < nz + 1
    error('nz=%d Randbedingungen gegeben. %d Maximalwerte f�r z erforderlich. Nur %d gegeben.', ...
        nz, nz+1, length(zmax));
end

if length(z0) ~= length(zT)
    error('Links und Rechts unterschiedliche Randbedingungen vorgegeben. Nicht implementiert!');
end


Tv = zmax(1:end-1)./zmax(2:end);



%% Berechnen
if nz < 5 && nz ~= 2
    [ew_t, ew_z] = Trapez_nAbl_Analytisch(z0, zT, t0, zmax, T_Abt, debug);
else
%     try
%         [ew_t, ew_z] = Trapez_nAbl_Iterativ(z0, zT, t0, zmax, Tmin, T_Abt, 0);
%     catch e
%         warning('Trapez_nAbl:Iterativ_Fehler', 'Fehler beim iterativen L�sen');
%         getReport (e)
        [ew_t, ew_z] = Trapez_nAbl_Faltung(z0, zT, t0, zmax, T_Abt, debug);
%     end
end

%% Pr�fen

% Pr�fe, ob die Eckwerte f�r alle Ableitungen korrekt sind.

for it = 1:length(ew_t)-1
    % Polynomkoeffizienten des Polynoms auf diesem Segment berechnen
    polykoeff = zeros(nz+1, nz+1); % Polynomkoeffizienten: h�chster zuerst
    polykoeff(nz+1, nz+1) = ew_z(it, nz+1); % letzte Ableitung ist konstant
    % Integration von der h�chsten Wegableitung z(nz+1) zum Weg z(1)
    for iz = nz:-1:1
        polykoeff(iz, :) = polyint(polykoeff(iz+1, 2:nz+1), ew_z(it, iz));
    end
    
    % Pr�fe, ob berechneter n�chster Eckwert mit tats�chlichem
    % �bereinstimmt
    for iz = 1:nz
        ew_iz_berechnet = polyval(polykoeff(iz, :), ew_t(it+1)-ew_t(it)); % mit Polynom-Integration berechnet
        ew_iz_vorgabe = ew_z(it+1, iz); % Aus Faltung oder analytischer Berechnung vorgegeben.
        Toleranz = 10^(-7+iz); % Toleranz erh�hen mit iz
        if abs(ew_iz_berechnet-ew_iz_vorgabe) > Toleranz
            error('z%d(t=%1.4f) stimmt nicht. Berechnet: %1.5e, Vorgegeben: %1.5e. Diff: %1.5e. (it=%d)', ...
                iz, ew_t(it+1), ew_iz_berechnet, ew_iz_vorgabe, ew_iz_berechnet-ew_iz_vorgabe, it+1);
        end
    end
    
    % Vergleiche
end




% Pr�fe, ob die Eckpunkte richtig sind
FF = abs(ew_z(end, 1:nz)-zT(1:nz)');
II = FF > 1e-9;
if any(II)
    I = find(II,1);
    error(['Der Fehler der Trajektorienberechnung ist gr��er als die Toleranz 1e-9. ', ...
        'Ist: z%d(T)=%f; Soll: z%d(T)=%f'], I, ew_z(end, I) ,I, zT(I));
else
    ew_z(end, 1:nz) = zT(1:nz)';
end

FF = abs(ew_z(1, 1:nz)-z0(1:nz)');
II = FF > 1e-9;
if any(II)
    I = find(II,1);
    error(['Der Fehler der Trajektorienberechnung ist gr��er als die Toleranz 1e-9. ', ...
        'Ist: z%d(0)=%f; Soll: z%d(0)=%f'], I, ew_z(1,I), I, z0(I));
else
    ew_z(1, 1:nz) = z0(1:nz)';
end

% Berechne Zwischenwerte, falls gefordert
if nargout >= 3
    w_t = (t0:T_Abt:ew_t(end))';
    w_z = Trapez_nAbl_Werte(ew_t, ew_z, w_t);
end