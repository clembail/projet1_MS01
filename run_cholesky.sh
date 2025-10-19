#!/bin/bash
#=============================================================================
# SCRIPT DE SOUMISSION SLURM POUR TESTS POISSON (JACOBI / GAUSS-SEIDEL)
#
# Lance des tests de scalabilité forte et faible.
#
# HYPOTHÈSE IMPORTANTE :
# Les exécutables acceptent 4 arguments : <Nx> <Ny> <MaxIter> <Tolerance>
# Et impriment leur sortie sous la forme :
#   Temps: X.XXX
#   Iterations: YYYY
#   Error: Z.ZZZe-ZZ
#
# Utilisation : sbatch ce_script.sh
#=============================================================================

# --- Configuration Slurm (A ADAPTER) ---
#SBATCH --job-name=poisson_scalability   # Nom du job
#SBATCH --partition=cpu_partition        # <-- A CHANGER: partition du cluster
#SBATCH --account=votre_projet           # <-- A CHANGER: nom de votre projet/compte
#SBATCH --nodes=8                        # <-- A CHANGER: 8 nœuds * 64 cœurs/nœud = 512
#SBATCH --ntasks-per-node=64             # <-- A CHANGER: Nb de cœurs par nœud
#SBATCH --ntasks=512                     # <-- A CHANGER: Nb TOTAL de processus (Nœuds * ntasks-per-node)
#SBATCH --cpus-per-task=1                # 1 CPU par tâche MPI
#SBATCH --time=04:00:00                  # <-- A CHANGER: temps max (HH:MM:SS)
#SBATCH --output=slurm_job_scalability_%j.out
#SBATCH --error=slurm_job_scalability_%j.err

# --- Configuration des Tests (A ADAPTER) ---

# Fichier de sortie
OUTPUT_FILE="scalability_results_$(date +%Y%m%d_%H%M).csv"

# Exécutables (adaptez les chemins si nécessaire)
J_SEQ="./jacobi_sequentiel"
J_PAR="./jacobi_parallele"
G_SEQ="./gauss_seidel_sequentiel"
G_PAR_RB="./gauss_seidel_parallele_rb" # Version Rouge/Noir

# Paramètres de simulation
MAX_ITER=50000
TOLERANCE=1e-6

# Liste des nombres de processus 'p' à tester
P_VALUES="2 4 8 16 32 64 128 256 512" # Doit être <= --ntasks

# --- Paramètres de Scalabilité ---

# 1. SCALABILITE FORTE (Problème Fixe, 'p' augmente)
NX_STRONG=2048
NY_STRONG=2048

# 2. SCALABILITE FAIBLE (Travail/processeur Fixe)
# On fixe Ny, et on fait augmenter Nx proportionnellement à 'p'
# (Nx = p * NX_PER_PROC_WEAK)
# Cela correspond à votre décomposition en bandes.
NY_WEAK=2048
NX_PER_PROC_WEAK=256 # Taille Nx de base par processeur

# --- Préparation de l'environnement ---

echo "Début des tests de scalabilité sur $(hostname)"
echo "Allocation Slurm : $SLURM_JOB_ID"
echo "Nb tâches max demandées : $SLURM_NTASKS"
echo ""

# Charger les modules nécessaires (ex: compilateur, bibliothèque MPI)
# module purge
# module load gcc/11.2.0        # <-- A DECOMMENTER / ADAPTER
# module load openmpi/4.1.1     # <-- A DECOMMENTER / ADAPTER

# Créer le fichier CSV et écrire l'en-tête
echo "Methode,Version,Test_Type,Nx,Ny,Processus_p,Temps_s,Iterations,Erreur_Finale" > $OUTPUT_FILE

# --- Fonction d'exécution ---
# Lance un test et enregistre les résultats
# Args: 1:EXE, 2:METHODE, 3:VERSION, 4:TEST_TYPE, 5:NX, 6:NY, 7:P
run_and_log() {
    local EXE=$1
    local METHODE=$2
    local VERSION=$3
    local TEST_TYPE=$4
    local NX=$5
    local NY=$6
    local P=$7

    echo "--- Lancement: $METHODE $VERSION ($TEST_TYPE) | N=${NX}x${NY} | p=$P ---"
    
    # Construction de la commande avec les 4 arguments
    local CMD="srun -n $P --exclusive $EXE $NX $NY $MAX_ITER $TOLERANCE"

    # Exécution : capture stdout et stderr (2>&1) dans la variable OUTPUT
    local OUTPUT
    OUTPUT=$($CMD 2>&1)

    # --- PARSING DES RÉSULTATS (basé sur les 'Nouveaux cout') ---
    
    local TIME_S
    TIME_S=$(echo "$OUTPUT" | grep -o 'Temps: [0-9.]*' | awk '{print $2}')
    [ -z "$TIME_S" ] && TIME_S="ERREUR_TEMPS"

    local ITERS
    ITERS=$(echo "$OUTPUT" | grep -o 'Iterations: [0-9]*' | awk '{print $2}')
    [ -z "$ITERS" ] && ITERS="NA"

    local ERROR
    ERROR=$(echo "$OUTPUT" | grep -o 'Error: [0-9.]*[eE]\{0,1\}[-+]\{0,1\}[0-9]*' | awk '{print $2}')
    [ -z "$ERROR" ] && ERROR="NA"

    # Écrire dans le CSV
    echo "$METHODE,$VERSION,$TEST_TYPE,$NX,$NY,$P,$TIME_S,$ITERS,$ERROR" >> $OUTPUT_FILE
    
    echo "Temps: $TIME_S s, Iters: $ITERS, Erreur: $ERROR"
    echo "----------------------------------------------------"
    echo ""
}

# --- Lancement des Tests ---

# 1. BASELINE SÉQUENTIEL (sur la taille "Forte")
echo "--- DEBUT TEST BASELINE (p=1, N=${NX_STRONG}x${NY_STRONG}) ---"
run_and_log $J_SEQ "Jacobi" "Sequentiel" "Baseline" $NX_STRONG $NY_STRONG 1
run_and_log $G_SEQ "Gauss-Seidel" "Sequentiel" "Baseline" $NX_STRONG $NY_STRONG 1


# 2. SCALABILITÉ FORTE (Problème Fixe)
echo "--- DEBUT SCALABILITE FORTE (N=${NX_STRONG}x${NY_STRONG} fixe) ---"
for p in $P_VALUES; do
    run_and_log $J_PAR "Jacobi" "Parallele" "Forte" $NX_STRONG $NY_STRONG $p
    run_and_log $G_PAR_RB "Gauss-Seidel" "Parallele (R/N)" "Forte" $NX_STRONG $NY_STRONG $p
done


# 3. SCALABILITÉ FAIBLE (Taille proportionnelle à 'p')
echo "--- DEBUT SCALABILITE FAIBLE (Nx = p * $NX_PER_PROC_WEAK, Ny = $NY_WEAK) ---"
for p in $P_VALUES; do
    # Calcule la taille Nx pour ce 'p'
    CURRENT_NX=$(( $NX_PER_PROC_WEAK * $p ))
    
    run_and_log $J_PAR "Jacobi" "Parallele" "Faible" $CURRENT_NX $NY_WEAK $p
    run_and_log $G_PAR_RB "Gauss-Seidel" "Parallele (R/N)" "Faible" $CURRENT_NX $NY_WEAK $p
done

echo "--- Tous les tests sont terminés. ---"
echo "Résultats enregistrés dans $OUTPUT_FILE"