% POLYMAX Suche Maximalwert- und Minimalwert eines Polynoms in einem Bereich
%
% Eingaben:
% polycoeff: Koeffizienten, beginnend beim höchsten
% xmin:      Beginn des betrachteten Bereichs
% xmax:      Ende des ...

% MA Moritz Schappler, schapplm@stud.uni-hannover.de, 2013-11
% Institut für mechatronische Systeme, Universität Hannover
% Betreuer: Daniel Beckmann, Daniel.Beckmann@imes.uni-hannover.de

function [maxval, minval, maxstelle, minstelle] = polymax(polycoeff, xmin, xmax)

% Ableiten
polyablcoeff = polyder(polycoeff);
% Null setzen
stellen = real(roots(polyablcoeff)); % rundungsbedingte Imaginärteile entfernen
% Randwerte hinzufügen
stellen = [stellen; xmin; xmax];

% Alle möglichen Extremstellen prüfen
stellen = stellen((stellen>=xmin) & (stellen<=xmax));

% Werte berechnen
polyextr = polyval(polycoeff, stellen);

[maxval, imax] = max(polyextr);
[minval, imin] = min(polyextr);

% Stellen der Maximal- und Minimalwerte
maxstelle = stellen(imax);
minstelle = stellen(imin);
