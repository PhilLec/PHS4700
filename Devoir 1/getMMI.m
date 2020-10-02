 # Retourne la matrice du moment d'inertie d'un objet dans le référentiel global
   # mi : moment d'inertie de l'objet dans son propre référentiel
   # m : masse de l'objet
   # cm : centre de masse de l'objet dans le référentiel du patineur
 function [mmi] = getMMI(pos, mi, m, cm)
   mi = [mi([1]), 0, 0; 0, mi([2]), 0; 0, 0, mi([3])];
   d = pos - cm;
   
   mmi = mi + m * [d([2])^2 + d([3])^2, -d([1]) * d([2]), -d([1]) * d([3]);
          -d([2]) * d([1]), d([1])^2 + d([3])^2, -d([2]) * d([3]);
          -d([3]) * d([1]), d([3]) * d([2]), d([1])^2 + d([2])^2];
 end