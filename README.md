# Projekt zaliczeniowy Genomika Porównawcza

### Instrukcja uruchomienia programu
Wymagane zainstalowane pakiety: Python: Biopython, numpy, matplotlib, requests, copy, sys, re, json, os; R: TreeTools, ape.
Wtmagane zainstalowane programy: muscle, mmseqs2, duptree.

#### Dane
Są dwie opcje podania danych o gatunkach do skryptu:

1. plik z nazwami gatunków przedzielony nowymi liniami, np.:

Aspergillus nidulans

Candida parapsilosis

Magnaporthe grisea\\
Debaryomyces hansenii\\
Candida tenuis
Torulaspora delbrueckii
Lachancea thermotolerans
Tetrapisispora blattae
Aspergillus fumigatus
Saccharomyces arboricola
Schizosaccharomyces japonicus

2. plik z nazwami gatunków oraz ID proteomów z bazy InterPro przedzieliony nowymi liniami, np.:
Aspergillus nidulans,UP000000560
Candida parapsilosis,UP000005221
Schizosaccharomyces japonicus,UP000001744

Pierwsza opcja służy do przypadków, gdy mamy bardzo dużo gatunków i nieefektywne byłoby wyszukiwanie wszystkich ID proteomów, natomiast może ona się mylić. Ten sposób zakłada, że po przeszukaniu bazy proteomów InterPro danym gatunkiem, pierwszy wynik jest poprawnym proteomem. Zwykle tak się zdarza, ale niestety nie zawsze, dlatego korzystając z tej opcji, trzeba przeprowadzić walidację, czy pobrane proteomy są poprawne, np. poprzez sprawdzenie, czy któryś proteom nie ma zaskakująco małej wielkości.

Druga opcja zapewnia, że zostaną konkretne wybrane proteomy, jednak wymaga podanie ich na wejściu.

Skorzystanie z pierwszej opcji:

Na początku skryptu należy ustawić zmienną SPECIES na True oraz podać nazwę pliku oraz zmienną IDS pozostawić ustawioną na False.

Skorzystanie z drugiej opcji:

Na początku skryptu ustawić zmienną IDS na True oraz podać nazwę pliku oraz zmienną SPECIES pozostawić ustawioną na False.

#### Pozostałe parametry

Na początku skryptu należy również podać szereg innych parametrów do programu:

working_dir = *str: ścieżka do folderu, w którym mamy plik wejściowy i w którym będą pojawiać się wyniki*

merged_fasta_output_filename = *str: nazwa pliku, który połączy wszystkie proteomy*

path_to_mmseq = *str: ścieżka do programu mmseqs2*

IDENT = *float: minimalny poziom identyczności do klastrowania*

COVER = *float: procent pokrycia do klastrowania*

COV_MODE = *int: typ stosowanego pokrycia w klastrowaniu, "coverage mode"*

PREF = *str: prefix plików wyjściowych z klastrowania*

path_to_muscle = *str: ścieżka do programu muscle*

TREE_CONSTRUCTION_METHOD = *str: metoda budowania drzew w Phylo.TreeConstruction, np. 'nj' lub 'mp'*

DISTANCE_METHOD = *str: metoda liczenia odległości między sekwencjami do budowania drzew, np. 'blosum62', 'pam70'*

PATH_RSCRIPT = *str: ścieżka do Rscript w celu wywołania skryptu R, w przypadku korzystania z systemu operacyjnego linux, prawdopodobnie wystarczy podać napis 'Rscript'*

PATH_DUPTREE = *ścieżka do programu duptree*

BOOTSTRAP_REP = *int: liczba drzew generowanych za pomocą metody bootstrap*

Aby uruchomić program należy mieć w katalogu, w którym pracujemy skrypty: projekt_gp.py, consensus.r, download_proteome.py oraz uruchomić skrypt projekt_gp.py.
