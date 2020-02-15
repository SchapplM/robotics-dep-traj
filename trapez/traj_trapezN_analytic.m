% Berechnung einer Trapez-Trajektorie für die n. Ableitung einer Größe.
% Die Berechnung erfolgt analytisch bis zur 4. Ordnung
% 
% 
% Eingabe
% z0        Anfangswert für die Größe z(1) und alle Ableitungen (z(2), ...)
% zT        Endwert ...
% t0        Anfangszeit
% zmax      Maximalwert für die Größe und alle Ableitungen 
% T_Abt     Abtastzeit des Trapez-Profils. Die Eckzeiten sind ganzzahlige
% debug     Zusätzliche Informationen ausgeben
% 
% Ausgabe
% ew_t      Eckwerte der Zeiten [Nx1]
% ew_z      Eckwerte der Größe z(1) und ihren Ableitungen z(2:end) 
%           [N x (nz+1)]
% 
% 
% Quelle:
% [1] Lambrechts 2003 - Trajectory planning and feedforward design for
% electromechanical motion systems
% https://de.mathworks.com/matlabcentral/fileexchange/16352-advanced-setpoints-for-motion-systems

% MA Moritz Schappler, schapplm@stud.uni-hannover.de, 2014-01
% Institut für Mechatronische Systeme, Universität Hannover
% Betreuer: Daniel Beckmann, Daniel.Beckmann@imes.uni-hannover.de


function [ew_t, ew_z] = traj_trapezN_analytic(z0, zT, t0, zmax, T_Abt, debug)
%% Initialisierung
nz = length(z0);


% Auf negativen Weg prüfen
if zT(1) < z0(1) % Berechnung für negative Berechnung analog zu positiver. Geschwindigkeiten etc. einfach invertieren
    minus = true;
    [z0, zT] = deal(zT, z0); % tausche Anfang und Ende
elseif zT(1) == z0(1) % Keine Bewegung: sofortiges Beenden
    % Die Anzahl der Eckzeiten sollte der gewünschten Trajektorienordnung
    % entsprechen
    nt = 2^nz;
    ew_z = [repmat(z0', nt, 1), zeros(nt,1)];
    ew_t = repmat(t0, nt, 1);
    % die Trajektorie darf nicht Dauer 0 haben, ansonsten Fehler in folgenden Funktionen
    ew_t(end) = t0+T_Abt;
    return;
else
    minus = false;
end

%% Berechnung
if nz == 1 % Geschwindigkeits-Rechteck
    h = zT(1)-z0(1);
    T = h/zmax(2);
    T = ceil(T/T_Abt)*T_Abt;
    v = h/T;
    ew_t = [0; T];
    ew_z = [[0; v], [h; 0]]';
elseif nz == 2 % Geschwindigkeits-Trapez
    error('Keine Berechnung einer Trajektorie mit nur Geschwindigkeits-RB vorgesehen.');
%     h = zT(1)-z0(1);
%     T1 = zmax(2)/zmax(3);
%     T1 = ceil(T1/T_Abt)*T_Abt;
%     T2 = h/zmax(2)-zmax(3)*T1^2;
%     T2 = ceil(T2/T_Abt)*T_Abt;
%     ew_t = [0; T1; T1+T2; 2*T1+T2];
%     ew_z = [0;0;];
elseif nz == 3 % Beschleunigungs-Trapez (Trapez-2)
    [t,jd]=make3(zT(1)-z0(1), zmax(2), zmax(3), zmax(4), T_Abt);
    if debug, figure(1);clf;end
    [jj,~,~,~,~,~,tt]=profile3(t,jd(1),T_Abt,debug);
    ew_t = tt';
    ew_z = NaN(length(ew_t), nz+1); % Werte vorbelegen
    ew_z(:, end) = jj(end, 2:2:end);  % Ruck setzen
    
    ew_z([4 5], 3) = 0;
elseif nz == 4 % Ruck-Trapez (Trapez-3)
    % Segmentdauern und Ruckableitung berechnen
    [t,dd]=make4(zT(1)-z0(1), zmax(2), zmax(3), zmax(4), zmax(5), T_Abt);
    % Eckzeiten berechnen
    if debug, figure(1);clf;end
    [dj,~,~,~,~,~,~,tt]=profile4(t,dd(1),T_Abt,debug);
    ew_t = tt';
    
    % restliche Eckwerte berechnen
    ew_z = NaN(length(ew_t), nz+1); % Werte vorbelegen
    ew_z(:, end) = dj(end, 2:2:end);  % Ruckableitungen setzen
    
%     
    % Bekannte Stellen Null setzen (verringert Integrationsfehler)
    % dabei wird von 16 Eck-Zeiten ausgegangen (Null-Segmente bleiben in
    % make4(...) erhalten)
%     ew_z(4, [4 5 8 9 12 13]) = 0; % Ruck
%     ew_z(3, [8 9]) = 0; % Beschleunigung
    
    
    
elseif nz>4
    error('Berechnung von Trapez-Trajektorien höherer Ordnung als 3 (Gefordert: %d) nicht möglich.', nz-1);
end

if minus % Prüfe, ob Verlauf negativ ist
    ew_z(:, 2:end) = - ew_z(:, 2:end); % invertieren der Geschwindigkeiten
    [z0, zT] = deal(zT, z0); % tausche Anfang und Ende zurück
end

% Anfangswert setzen
ew_z(1, 1:nz) = z0(1:nz)'; 

% Endwert setzen
ew_z(end, 1:nz) = zT(1:nz)'; 

% restliche Werte durch Integration    
[~, ew_z] = traj_trapezN_values(ew_t, ew_z, 0); 

ew_t = ew_t + t0;




%% Zeichnen    
if debug
    figure(2);clf;
    w_t_fein = (t0:1e-5:ew_t(end))';
    [w_z_fein, ~] = traj_trapezN_values(ew_t, ew_z, w_t_fein);
%     for iz=nz+1:-1:1
%         ax(iz)=subplot(nz+1, 1, iz);hold all;
%         plot(tt, dj(end, 2:2:end), 'o')
    for iz=nz+1:-1:1
        subplot(nz+1, 1, iz);hold all;
        if iz ~= nz+1
            hdl=plot(w_t_fein, w_z_fein(:, iz));
            legend_append(hdl, 'w_z_fein');
%             hdl=plot(w_t_fein, w_z_neu_fein(iz, :), '--m');
%             legend_append(hdl, 'w_z_neu_fein');
        else
            stairs(ew_t, ew_z(:, iz));
        end
        plot(ew_t, ew_z(:, iz), 'ro');
%             plot(ew_t, ew_z_int(iz, :), 'go');
        plot(ew_t([1, end])', repmat(zmax(iz), 2, 1), 'r--');
        ylabel(sprintf('z%d = d^%d y_f / dt^%d', iz, iz-1, iz-1));

    grid on
    end
    linkxaxes;
end
