# Objectifs :

Faire un programme qui récupère les informations d'un fichier SAM.

## Contraintes
- Input: fichier.sam
- Output: fichiers .txt (tabulé) et .png (graphique)
- readme (ce fichier doit contenir toutes les informations pour installer et utiliser votre script)
- docstring (help)
- rapport (10 pages)

## Comment sera structuré le programme ?

- Module ouverture fichier
	Parser le fichier :
	- récupérer header (logiciel mapping + version, nb de reads total)	> sortie ????
	- récupérer sur chaques lignes (N° read 1et2, flag 1et2, cigar 1et2 autre ???) > sortie ????
	
- Module Flag
	- récupère les lignes reads en entrée
	- travaille avec Nom read et flag
	- prérempli un tableau de sortie
	
- Module Cigar
	- récupère les lignes reads en entrée
	- travaille avec nom read et cigar
	- complète le tableau de sortie (cf flag)


- Main qui affiche les stats Flag et stats Cigar


