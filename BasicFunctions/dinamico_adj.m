function adj = dinamico_adj(screw) % optimized on 30.05.2022
adj = [0 -screw(3) screw(2) 0 0 0;screw(3) 0 -screw(1) 0 0 0;-screw(2) screw(1) 0 0 0 0;...
       0 -screw(6) screw(5) 0 -screw(3) screw(2);screw(6) 0 -screw(4) screw(3) 0 -screw(1);-screw(5) screw(4) 0 -screw(2) screw(1) 0];