# libOL version 1
Localisation géométrique rapide par arbre octal

# Objectif
La libOL permet de générer un arbre de localisation autour d'un nuage de points ou d'un ensemble de triangles.  
Elle retourne une simple étiquette associée à chaque maillage et permet par la suite d'effectuer des requêtes très rapidement afin de retrouver l'entité la plus proche d'une certaine coordonnée ou bien la liste des entités incluses dans une boîte rectangulaire.  
L'implémentation qui en est faite favorise une très grande compacité mémoire afin de gérer de très gros maillages sur un simple ordinateur portable.

# Compilation
Entrez les commandes dans l'ordre suivant :
- désarchivez le fichier ZIP
- mkdir build
- cd build
- cmake .
- make
- make install

# Utilisation
Il s'agit d'un simple fichier écrit en ANSI C et d'un header associé qu'il suffit d'inclure et de compiler avec son propre code.  
Elle peut être utilisée à partir des langages C, C++, F77 et F90.  
Elle a été testée sous les systèmes Linux et Mac OS X.
