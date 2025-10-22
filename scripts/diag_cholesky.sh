#!/bin/bash
#=============================================================================
# SCRIPT DE DIAGNOSTIC SLURM (v4 - "Safe Debug")
#
# Objectif: Voir ce que le script capture VRAIMENT.
#=============================================================================

# --- Configuration Slurm (VERSION DEBUG RAPIDE) ---
#SBATCH --job-name=poisson_DIAGNOSTIC     # Nom différent
#SBATCH --partition=debug_partition      # <-- IMPORTANT: Utilisez une partition 'debug' si elle existe
#SBATCH --account=votre_projet           # (Inchangé)
#SBATCH --nodes=1                        # 1 seul nœud
#SBATCH --ntasks=8                       # 8 cœurs au total (p max = 8)
#SBATCH --time=00:05:00                  # 5 MINUTES
#SBATCH --output=slurm_diag_%j.out
#SBATCH --error=slurm_diag_%j.err

# --- Configuration des Tests (VERSION DEBUG RAPIDE) ---
OUTPUT_FILE="diag_results.csv"
J_SEQ="./jacobi_sequentiel"
J_PAR="./jacobi_parallele"
G_SEQ="./gauss_seidel_sequentiel"
G_PAR_RB="./gauss_seidel_parallele_rb"

MAX_ITER=100
TOLERANCE=1e-2
P_VALUES="2 4 8"

NX_STRONG=100
NY_STRONG=100
NY_WEAK=50
NX_PER_PROC_WEAK=50

# --- Préparation de l'environnement ---

echo "--- DEBUT DU SCRIPT DE DIAGNOSTIC ---"
echo "Job ID: $SLURM_JOB_ID"
echo "Nœuds: $SLURM_JOB_NODELIST"
echo ""

# =========================================================================
# !! IMPORTANT !! DÉCOMMENTEZ ET ADAPTEZ CES LIGNES
# Si les modules ne sont pas chargés, mpirun ne fonctionnera pas !
#
# module purge
# module load gcc/11.2.0        
# module load openmpi/4.1.1     
#
echo "Modules chargés:"
module list
# =========================================================================

echo "Methode,Version,Test_Type,Nx,Ny,Processus_p,Temps_s,Iterations,Erreur_Finale" > $OUTPUT_FILE

# --- Fonction de Diagnostic ---
run_and_log() {
    local EXE=$1
    local METHODE=$2
    local VERSION=$3
    local TEST_TYPE=$4
    local NX=$5
    local NY=$6
    local P=$7

    echo "--- Lancement: $METHODE $VERSION (N=${NX}x${NY} p=$P) ---"
    
    # Commande mpirun "silencieuse"
    local CMD="mpirun -n $P $EXE $NX $NY $MAX_ITER $TOLERANCE"
    
    # Pour les tests séquentiels (p=1), on n'utilise pas mpirun
    if [ $P -eq 1 ]; then
        CMD="$EXE $NX $NY $MAX_ITER $TOLERANCE"
    fi

    # Exécution : capture stdout et stderr (2>&1)
    local OUTPUT
    OUTPUT=$($CMD 2>&1)

    # --- SECTION DE DIAGNOSTIC ---
    echo "1. SORTIE BRUTE DU PROGRAMME (ce que voit le script) :"
    echo "--- DEBUT SORTIE BRUTE ---"
    echo "$OUTPUT"
    echo "--- FIN SORTIE BRUTE ---"
    echo ""

    # 2. VÉRIFICATION DU GREP (ce que trouve le script) :
    # On n'utilise pas 'awk' pour l'instant, on veut la ligne entière.
    local TIME_LINE
    TIME_LINE=$(echo "$OUTPUT" | grep 'Temps:')
    local ITERS_LINE
    ITERS_LINE=$(echo "$OUTPUT" | grep 'Iterations:')
    local ERROR_LINE
    ERROR_LINE=$(echo "$OUTPUT" | grep 'Error:')

    echo "3. RÉSULTAT DU PARSING :"
    echo "   Ligne Temps trouvée  : '$TIME_LINE'"
    echo "   Ligne Iters trouvée  : '$ITERS_LINE'"
    echo "   Ligne Error trouvée  : '$ERROR_LINE'"
    echo "----------------------------------------------------"
    echo ""

    # Écrire des marqueurs de debug dans le CSV
    echo "$METHODE,$VERSION,$TEST_TYPE,$NX,$NY,$P,DEBUG,DEBUG,DEBUG" >> $OUTPUT_FILE
}

# --- Lancement des Tests ---

echo "Lancement Baseline (p=1)..."
run_and_log $J_SEQ "Jacobi" "Sequentiel" "Baseline" $NX_STRONG $NY_STRONG 1

echo "Lancement Scalabilité Forte (p=2)..."
run_and_log $J_PAR "Jacobi" "Parallele" "Forte" $NX_STRONG $NY_STRONG 2

echo "Lancement Scalabilité Faible (p=2)..."
local NX_WEAK_2=$(( $NX_PER_PROC_WEAK * 2 ))
run_and_log $J_PAR "Jacobi" "Parallele" "Faible" $NX_WEAK_2 $NY_WEAK 2

echo "--- Diagnostic terminé. ---"