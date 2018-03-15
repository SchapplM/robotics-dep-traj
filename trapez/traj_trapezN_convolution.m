% Berechnung einer Trapez-Trajektorie für die n. Ableitung einer Größe.
% 
% Algorithmus basierend auf der mehrfachen Faltung der Zeitverläufe
% 
% Eingabe:
% z0        Anfangswert für die Größe z(1) und alle Ableitungen (z(2), ...)
% zT        Endwert ...
% t0        Anfangszeit
% zmax      Maximalwert für die Größe und alle Ableitungen 
% T_Abt     Abtastzeit des Trapez-Profils. Die Eckzeiten sind ganzzahlige
%           Vielfache der Abtastzeit
% debug     Ausgabe zusätzlicher Informationen und Bilder
% 
% Ausgabe:
% ew_t      Eckwerte der Zeiten
% ew_z      Eckwerte der Größe z(1) und ihren Ableitungen z(2:end)
% w_z       Alle Werte der Größe z(1) und ihren Ableitungen z(2:end) 
%           zu den Zwischenzeiten w_t. [n x (nz+1)]
% w_t       Alle Zeitschritte von t0 bis zum Ende der Trajektorie mit dem
%           Abstand T_Abt. [n x 1]
% 
% Quelle:
% [1] Lee 2013 - Convolution-Based Trajectory Generation Methods Using Physical
% System Limits 


% MA Moritz Schappler, schapplm@stud.uni-hannover.de, 2014-01
% Institut für mechatronische Systeme, Universität Hannover
% Betreuer: Daniel Beckmann, Daniel.Beckmann@imes.uni-hannover.de

function [ew_t, ew_z, w_z, w_t] = traj_trapezN_convolution(z0, zT, t0, zmax, T_Abt, debug)

% Beispiel aus [1], S.7
% z0=[0, 0, 0, 0, 0,0];
% zT=[8, 0, 0, 0, 0,0];
% zmax=[8, 4, 4, 16,64];
% load('test.mat');
% save('test.mat');
% debug=0;



if zT(1) < z0(1) % Berechnung für negative Berechnung analog zu positiver. Geschwindigkeiten etc. einfach invertieren
    minus = true;
    [z0, zT] = deal(zT,z0);% tausche Anfang und Ende
elseif zT(1) == z0(1) % Keine Bewegung: sofortiges Beenden
    ew_z = [[z0',0]; [zT',0]]; % die Trajektorie darf nicht Dauer 0 haben,
    ew_t = [t0; t0+T_Abt]; % ...ansonsten Fehler in folgenden Funktionen
    return;
else
    minus = false;
end


% Zusätzliche RB hinzufügen mit minimaler Verschliffzeit
% z0 = [z0; 0];
% zT = [zT; 0];
% zmax = [zmax; zmax(end)/T_Abt];

nz = length(z0);


% Maximalwerte der höchsten Ableitung anpassen
if length(zmax) == nz
    zmax = [zmax; zmax(nz)/T_Abt(nz)];
else
    % Die letzte Ableitung wurde als Zahlenwert angegeben. Vergleiche mit
    % gegebener Abtastrate. Mindestwert beachten
    Tv_gegeben = zmax(nz)/zmax(nz+1);
    if Tv_gegeben < max(T_Abt)
        zmax(nz+1) = zmax(nz)/max(T_Abt);
    end
    zmax = zmax(1:nz+1);
end


% zmax_orig=zmax; % ursprüngliche Max-Werte speichern
% zmin_orig = -zmax_orig;
% zmin_orig(1) = z0(1);

% Weg-Maximalwert anpassen
zmax(1) = zT(1)-z0(1);
zmin = -zmax;

%%
% [1], Seite 3; Es werden für Ableitungen der Geschwindigkeit bis zur
% Ordnung n Grenzen vorgegeben

% Randbedingungen bis zur (nz-2). Ableitungen der Geschwindigkeit.
%  -> Maximalwerte zu nz-3 Ableitungen der Geschwindigkeit
n = nz-3; 
% [1], Formel (5)
t = NaN(n+1,1);
for k = 0:n+2
    i = k+1;
    t(i) = zmax(k+1) / zmax(k+2); % Probe: k=0 -> Weg/Geschwindigkeit
    
    % Zeit an Abtastzeit anpassen und höhere Ableitung korrigieren
    t(i) = max(ceil(t(i)/T_Abt)*T_Abt, T_Abt);
    %if i>1
        zmax(k+2) = zmax(k+1)/t(i);
    %end
    if t(i) == T_Abt
        warning('Eine Verschliffzeit (t(%d)) entspricht der Abtastzeit. Damit keine Faltung wirksam! Verdoppele.', i);
        t(i) = 2*T_Abt;
    end
end

% [1], S. 4 oben:
% Die Zeitlängen müssen ganzzahlige Vielfache der Abtastzeit sein!
t = max(ceil(t/T_Abt)*T_Abt, T_Abt);
% t(end)=2*t(end);
% t(3)=10e-3;
% [1], Formel (6)
for l = n+2:-1:0%0:n+2
    i = l+1;
    if t(i) < sum(t(i+1:end))
        fprintf('Verschliffzeiten: [%s] ms\n', disp_array(1e3*t', '%1.3f'));
        error('Die Zeit t%d (%f) ist kleiner als die Summe t%d bis t%d (%f). Passe an!', ...
            l, t(i), l+1, n+2, sum(t(i+1:n+3)));
        % t(i) = sum(t(i+1:end)) + T_Abt;
    end
end

% Grenzwerte aktualisieren
for i = 2:nz
    zmvh = zmax(i+1);
    zmax(i+1) = zmax(i)/t(i);
    if debug && zmvh~=zmax(i+1)
        fprintf('zmax(%d): %1.4e -> %1.4e\n', i+1, zmvh, zmax(i+1));
    end
end


n_Zusatz = nz-1;
w_t = (0:T_Abt:sum(t))';
w_t = [w_t; w_t(end)+T_Abt*(1:nz-1+n_Zusatz)'];
nS = length(w_t);

%% Faltung für alle Ableitungen durchführen
% [1], Bild 1
y_vorh = zeros(nS, nz+1);
% Setze y0
y_vorh(:, 2) = repmat(zmax(2), nS, 1);%[zeros(nz-1+8, 1); repmat(zmax(2), nS-(nz-1+8), 1)];
% Funktion y0 nach t0 ausblenden
t_Anf = n_Zusatz*T_Abt;
t_End = t_Anf+t(1);
y_vorh(w_t<t_Anf, 2)=0;
y_vorh(w_t>=t_End, 2)=0;

y_nach = y_vorh;
% Schleife über durchführung der Faltung. 
for i = 1:n+3
    %% Faltung der Geschwindigkeit und aller verfügbaren Ableitungen durchführen
    if i>1
        % Schleife über Ableitungen der Geschwindigkeit
        % [1], Formel (9); Zähler i statt n
        for iz = 2:i
            % zuerst: Berechnung von yn mit n=1 (y0 ist konstant)

            mi = t(i)/T_Abt;
            for k = 1:nS
                Y1 = y_vorh(k, iz);
                if k-mi < 1
                    Y2 = 0;
                else
                    Y2 = y_vorh(floor(k-mi), iz);
                end
                if k<2
                    Y3 = 0;
                else
                    Y3 = y_nach(k-1, iz);
                end
                y_nach(k, iz) = (Y1-Y2)/mi + Y3;
            end
        end
    else
        mi = 0;
    end
    %% Weg integrieren
    % gefaltete Geschwindigkeit integrieren (wird als richtig
    % angesehen) (numerischer Fehler, da Geschw. Polynom höheren Grades)
    y_nach(:, 1) = cumtrapz(w_t, y_nach(:, 2));
    %% Gefaltete Ableitungen korrigieren
    % Alle Ableitungen der Geschwindigkeit, die im vorherigen Schritt
    % gefaltet wurden, müssen korrigiert werden.
    for iz = 3:i
        % gefaltete Ableitung bis zum Weg integrieren
        s_Falt = y_nach(:, iz);
        for ia = 1:iz-1
            s_Falt = cumtrapz(w_t, s_Falt);
        end
        % Korrekturfaktor aus dem Verhältnis Ist- und Soll-Weg bestimmen
        K_Korr = y_nach(end, 1) / s_Falt(end);
        % Korrigieren
        y_nach(:, iz) = K_Korr*y_nach(:, iz);
        continue;
        figure(1);clf;hold all;plot(s_Falt);plot(y_nach(:, 1))
    end
    %% Nächste Verfügbare Ableitung der Geschwindigkeit bilden
    % Ableitung einer linearen Funktion = Rechteck (Ableitungen mit diff bilden)
    if i>1
        y_nach(:, i+1) = [diff(y_nach(:, i)); 0]/T_Abt; % ; 0
        % y_nach(:, i+1) = gradient(y_nach(:, i))/T_Abt;
        
        % Ableitungen bereinigen (Rundungsbedingte Impulse entfernen)
        y_nach(abs(y_nach(:, i+1))<0.1*max(y_nach(:, i+1)), i+1) = 0;
        
        % Ableitungen bereinigen: Darf nur einen konstanten Wert haben.
        % Setze konstante Höhe des Rechtecks ab 90% der Rechteckhöhe
        I_Rechteck = abs(y_nach(:, i+1))>0.9*max(y_nach(:, i+1));
        y_nach(I_Rechteck, i+1) = sign(y_nach(I_Rechteck, i+1)) * max(y_nach(:, i+1));
        
        % Es dürften keine Werte mehr übrig sein, die keiner
        % Rechteckfunktion entsprechen
        if any( abs(y_nach(:, i+1)) >= 0.1*max(abs(y_nach(:, i+1))) & ...
                abs(y_nach(:, i+1)) <= 0.9*max(abs(y_nach(:, i+1))) )
            warning('Nach der Differentiation der gefalteten Funktion sind Werte übrig, die keiner Rechteckfunktion entsprechen! Fehler!');
        end
    end
    %% Zeichnen
    if debug
        figure(6+i);clf;
        for iz=nz+1:-1:1

            ax(iz)=subplot(nz+1, 1, iz);hold all;
            if iz == 1
                title(sprintf(['Trajektorie: Nach Faltungsschritt %d. ', ...
                    'Verschliffen mit %1.1fms bzw. %d Schritten'], i-1, t(i)*1e3, mi)); 
            end
            % neue Werte
            hdl=stairs(w_t, y_nach(:, iz));
            if iz <= i
                legend_append(hdl, sprintf('s^{(%d)}_{%d}=s^{(%d)}_{%d}*', iz-1,i, iz-1,i-1));
            elseif iz == i+1
                legend_append(hdl, sprintf('s^{(%d)}_{%d}=ds^{(%d)}_{%d}/dt', iz-1,i, iz-2,i));
            else
                legend_append(hdl, sprintf('s^{(%d)}_{%d}: nicht berechnet', iz-1,i));
            end
            % alte Werte
            hdl = stairs(w_t, y_vorh(:, iz), '--');
            legend_append(hdl, sprintf('s^{(%d)}_{%d}', iz-1,i-1));
            % Grenzwert
            plot([min(w_t), max(w_t)], zmax(iz)*[1,1], 'r--'); % _orig
            plot([min(w_t), max(w_t)], zmin(iz)*[1,1], 'r--'); % _orig
            grid on;
            ylabel(sprintf('z%d', iz));
        end
        subplot_expand(gcf, nz+1, 1);%, ax);
    end    
    % Variablen tauschen (für nächste Iteration)
    y_vorh = y_nach;

end

%% Eckwerte berechnen

% Eckzeiten extrahieren
I = logical([true;diff(y_nach(:, nz+1))] ~= 0);
ew_t = w_t(I);


ew_z = NaN(sum(I), nz+1); % Vorbelegen

% letzte Ableitung setzen
ew_z(:, nz+1) = [y_nach(I(2:end-1), nz+1); 0];

% andere Ableitungen auch übernehmen
% ew_z(:, 1:nz) = [y_nach(I(2:end-1), 1:nz); NaN(1,nz)];

% die letzte Ableitung immer auf den gleichen Wert setzen.
% ACHTUNG: Funktioniert nicht!
% Dann ist der zurückgelegte Weg proportional zu diesem Wert.
% ew_z(end, :) = sign(ew_z(end, :)) * abs(ew_z(end, 1));

% Anfangswerte setzen
ew_z(1, 1:nz) = zeros(nz,1);%z0(1:nz); % Anfangsposition und -ableitungen

% Fehlende Eckwerte ergänzen (Segmente mit Länge Null)
% Die letzte Ableitung muss abwechselnd 1, 0, -1, 0, ... besitzen
% problematisch, da zwei aufeinanderfolgende Rechtecke auch verschmelzen
% können. Dann sind die Eckpunkte nicht mehr zuordnugnsfähig.
% if length(ew_t) ~= 2^(nz)
%     if debug
%         fprintf('%d Segmente. Bei %d Randbedingungen: %d Segmente erwartet.\n', ...
%             length(ew_t)-1, nz, 2^(nz)-1);
%     end
%     ii = 1;
%     ew_t_neu = NaN(1, 2^(nz)-1);
%     ew_z_neu = NaN(nz+1, 2^(nz)-1);
%     for it = 1:2^(nz) % alle neuen Eck-Zeitpunkte durchgehen
%         if mod(it, 4) == 1 || mod(it, 4) == 3 % positive Flanke
%             if ew_z(end, ii) ~= 0
%                 iO = 1; % in Ordnung
%             else
%                 error('Flanke erwartet und nicht bekommen');
%             end
%         else
%             if ew_z(end, ii) == 0
%                 iO = 1;
%             else
%                 % Null-Segment einfügen, was vorher nicht da war
%                 ew_t_neu(it) = ew_t(ii); % Zeitdauer Null
%                 ew_z_neu(end,it) = 0;
%                 iO = 0; % Ursprüngliche Daten nicht weiterzählen
%             end
%         end
%         if iO
%             ew_t_neu(it) = ew_t(ii);
%             ew_z_neu(:,it) = ew_z(:, ii);
%             ii = ii+1; % Index der berechneten Werte ew_t
%         end
%     end
%     % neu berechnete Trajektorie übernehmen
%     ew_z = ew_z_neu;
%     ew_t = ew_t_neu;
% end

% Bekannte Null-Eckwerte setzen (funktioniert nicht richtig, da durch
% Ableitung und Falten die Trajektorie nicht symmetrisch sein kann.)
% ii = 1;
% for iz = nz:-1:nz % funktioniert nicht mit niedrigeren Ableitungen
%     % Anfangen mit Trapezprofil bei nz (alle 2^2-1=3 Eckwerte Null), dann
%     % nächsthöheres: alle 2^3-1=7 Eckwerte Null, usw...)
%     ii = ii+1;
%     I1 = (0:(2^(ii)):2^nz);
% 
%     ew_z(iz,  I1(2:end)) = 0;
%     ew_z(iz,  I1(1:end-1)+1) = 0;
% 
% end

% Verschiebung der gefalteten Funktionen ausgleichen
y_nach2 = y_nach;
y_nach2(:, 1:nz-1) = [y_nach2([2:end,end], 1:nz-1)];


if debug 
    % Zwischenwerte zum Zeichnen mit berechnen
    w_t_fein = (0:1e-3:w_t(end))';
    [w_z_fein, ~] = traj_trapezN_values(ew_t, ew_z, w_t_fein);
%     [w_z_neu_fein, ew_z_neu_int] = traj_trapezN_values(ew_t, ew_z_neu, w_t_fein);
else
    [~, ew_z] = traj_trapezN_values(ew_t, ew_z, 0);
end
% w_t_fein = w_t;
% w_z_fein = y_nach;

if debug
    figure(3);clf;
    for iz=1:nz+1
        % axes (ax(iz)); %#ok<LAXES>
        subplot(nz+1,1,iz);hold all
        if iz ~= nz+1
            hdl1=stairs(1e3*w_t_fein, w_z_fein(:, iz), '-c');
            %hdl2=stairs(1e3*w_t, y_nach(:, iz), '--k');
            hdl3=stairs(1e3*w_t, y_nach2(:, iz), '--b');
        else
            hdl1=stairs(1e3*w_t_fein, w_z_fein(:, iz), '-c');
            %hdl2=stairs(1e3*w_t, y_nach(:, iz), '--k');
            hdl3=stairs(1e3*w_t, y_nach2(:, iz), '--b');
        end
        legend_append(hdl1, 'w_z_fein');
        %legend_append(hdl2, 'y_nach');
        legend_append(hdl3, 'y_nach2');
        hdl=plot(1e3*ew_t, ew_z(:, iz), 'go');
        legend_append(hdl, 'ew_z');
        grid on;
        if iz == 1
            plot(1e3*w_t_fein(end), zT(1)-z0(1), 'rx');
        end
    end
    subplot_expand(gcf,nz+1,1);
    xlabel('t [ms]')
end

% ew_z = ew_z_int;

% Anfangswerte aufaddieren
ew_t = ew_t + t0;
w_t = w_t + t0;

%% Ergebnis prüfen und korrigieren
% Nullen am Anfang und Ende entfernen
% I = ew_z(end, :) == 0







% Fehler korrigieren, der evtl durch Rundung entsteht
for i = 1:3
    if abs(ew_z(end,1)-zmax(1))>1e-6
        debug=1;
        if debug
            fprintf('Der zurückgelegte Weg ist falsch. Soll: %f, Ist: %f. Verändere letzten Grenzwert.\n', ...
                zmax(1), ew_z(end,1));
        end
        zm_vorher = zmax(nz+1);
        zmax(nz+1) = abs((zmax(1)-ew_z(1,1))/(ew_z(end,1)-ew_z(1, 1))) * zmax(nz+1);
        % Alle Eckzeiten lassen, alle Eckpunkte ändern (manuell).
        ew_z(:,nz+1) = ew_z(:,nz+1) * zmax(nz+1)/zm_vorher;
        ew_z(2:end,1:nz) = NaN; % Eckwerte zurücksetzen
        [~, ew_z] = traj_trapezN_values(ew_t, ew_z, 0);
        if debug
            fprintf('Grenzwert Abl. %d neu gesetzt: %1.4e -> %1.4e\n', ...
                nz, zm_vorher, zmax(nz+1));
        end
    end
end

% Anfangs-Weg wieder draufaddieren
ew_z(:, 1) = ew_z(:, 1) + z0(1);

if minus % Prüfe, ob Verlauf negativ ist
    ew_z(:,2:end) = - ew_z(:,2:end); % invertieren der Geschwindigkeiten
    ew_z(:,1) = zT(1)+z0(1)-ew_z(:,1); % umkehren der Weg-Positionen. Ohne fliplr. Das führt zu Fehler
    % tausche Anfang und Ende zurück, damit im nächten Schritt geprüft
    % werden kann
    [z0, zT] = deal(zT,z0);% tausche Anfang und Ende zurück
    
%     plot(ew_t,ew_z(1, :))
%     plot(ew_t,fliplr(ew_z(1, :)))
%     hold all;plot(ew_t,zT(1)-ew_z(1, :))
end





% Prüfe, ob die Eckpunkte richtig sind

FF = abs(ew_z(end,1:nz)'-zT(1:nz));
Tol = 10.^(-6 +(1:nz))'; % bei höheren Ableitungen größere Toleranz
II = FF > Tol;
if any(II)
    I = find(II,1);
    error(['Der Fehler der Trajektorienberechnung bei Faltung ist größer als die Toleranz %e. ', ...
        'Ist: z%d(T)=%f; Soll: z%d(T)=%f. Diff=%e'],Tol(I), I, ew_z(end,I) ,I, zT(I),FF(I));
end
% exakten Wert einsetzen, um Rundungsfehler zu eliminieren
ew_z(end,1:nz) = zT(1:nz);

% Berechne Zwischenwerte, falls gefordert
if nargout >= 3
    w_t = t0:T_Abt:ew_t(end);
    w_z = traj_trapezN_values(ew_t, ew_z, w_t);
end

% Berechne die korrekten Eckwerte
[~, ew_z] = traj_trapezN_values(ew_t, ew_z, 0);
