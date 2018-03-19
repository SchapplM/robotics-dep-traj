# Trajektorien-Toolbox

Dieses Repo enthält eine Sammlung von Matlab-Funktionen zum Thema Trajektorien.

Moritz Schappler, schappler@imes.uni-hannover.de, 2018-03  
(C) Institut für mechatronische Systeme, Universität Hannover

## Benutzung

Für die Benutzung des Repos muss die Pfadinitialisierung durchgeführt werden mit dem Skript `trajectory_toolbox_path_init.m` aus dem Hauptordner.

### Abhängigkeiten

Einige Funktionen benötigen das Repo `matlab-ext` mit externen Funktionen und Toolboxen.

### Modultests

Die einzelnen Funktionen des Repos können über die Modultests aus dem Ordner `examples_tests` gestartet werden. Die Funktionsfähigkeit des Gesamt-Repos kann mit dem Gesamt-Testskript `trajectory_toolbox_test_repo.m` aus dem Hauptordner geprüft werden.

## Inhalt

Hier erfolgt eine kurze Übersicht der im Repo enthaltenen Funktionen. Die Funktionsweise ist in den jeweiligen Funktionsdateien und Modultests erklärt.

* Trapez-Profile beliebiger stetiger Differenzierbarkeit: `traj_trapezN.m`, `traj_trapezN_single.m`
* Kegel-Trajektorie im Arbeitsraum: `traj_rot_cone_angaxis.m`

## Erweiterung

Weitere allgemeine Funktionen sollen unter Verwendung des Benennungsschemas und grundsätzlicher Programmierkonventionen ergänzt werden und damit allen Institutsangehörigen zur Verfügung gestellt werden.

