# **Alignement de Séquences ADN - Distance d'Édition** 🧬
    Ce projet implémente plusieurs algorithmes pour calculer la distance d'édition entre deux séquences ADN en utilisant la programmation dynamique. L'objectif est d'optimiser l'accès mémoire avec des approches itératives, cache-aware, et cache-oblivious.

## 🚀 **Fonctionnalités**
✅ Alignement de séquences ADN avec l'algorithme de Needleman-Wunsch
✅ Implémentation de trois versions :

    - Itératif : Complexité mémoire réduite à O(N + M)
    - Cache-Aware : Optimisé pour la hiérarchie mémoire
    - Cache-Oblivious : Approche récursive sans connaissance des caches
✅ Analyse des performances et défauts de cache avec Valgrind

## 🛠️ **Technologies utilisées**
    - Langage : C
    - Compilation : GCC & Makefile
    - Analyse mémoire : Valgrind (Cachegrind)
    - Données : Fichiers FASTA du NCBI
## 📁 Structure du projet
```bash
ADN-Distance-Edition/
├── src/         # Code source (.c et .h)
├── srcperg      # Code source mesure
├── tests/       # Scripts et fichiers de test
├── bin/         # Exécutables
├── README.md    # Documentation du projet
└── Makefile     # Script de compilation
```
## 🔧 Installation et exécution
### 1️⃣ Compilation
```bash
    make
```
### 2️⃣ Exécution
```bash
    ./bin/distanceEdition <seq1.fna> <pos1> <len1> <seq2.fna> <pos2> <len2>
```
Exemple :

```bash
    ./bin/distanceEdition data/wuhan_hu_1.fasta 116 30331 data/ba52_recent_omicron.fasta 153 30183
```
## 🏆 Performance et Analyse
Simulation des défauts de cache avec Valgrind
Comparaison des temps d'exécution et consommation mémoire
Estimation du coût énergétique
## ✨ Améliorations possibles
Implémentation parallèle pour accélérer l’exécution
Optimisation SIMD pour réduire les accès mémoire
Export des résultats en format standardisé (CSV, JSON)
## 👥 Équipe : 
Mouez Jaafoura, Zine Ouadii

