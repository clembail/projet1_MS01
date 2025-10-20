#!/bin/bash
#=========================================================
# SCRIPT WORKER (ESCLAVE) - LANCE UN SEUL TEST
#
# Logique:
# 1. Reçoit 'P_PROCS' (ex: 8) en argument $1.
# 2. Demande P_PROCS cœurs à Slurm via #SBATCH --ntasks.
# 3. Lance mpirun SANS '-n'. mpirun lira SLURM_NTASKS.
#=========================================================

# --- 1. Récupération des arguments ---
P_PROCS=$1     # Ex: 8
EXE_PATH=$2    # Ex: ./buildClem/jacobi_parallele
TEST_TYPE=$3   # Ex: "Forte"
METHODE=$4     # Ex: "Jacobi"
NX=$5          # Ex: 2048
NY=$6          # Ex: 2048

# Fichier de sortie CSV (un par job)
OUTPUT_FILE="result_${METHODE}_${TEST_TYPE}_p${P_PROCS}.csv"

# --- 2. Configuration Slurm ---
# C'est ici que l'on fixe le nombre de processus
#SBATCH --job-name=CB_${METHODE}_${TEST_TYPE}_p${P_PROCS}
#SBATCH --partition=cpu_test
#SBATCH --account=ams301
#SBATCH --ntasks=${P_PROCS}           # Demande EXACTEMENT le bon nb. de cœurs
#SBATCH --time=04:00:00
#SBATCH --output=slurm_log_${METHODE}_p${P_PROCS}_%j.out
#SBATCH --error=slurm_log_${METHODE}_p${P_PROCS}_%j.err

# --- 3. Environnement ---
# !! ETAPE CRUCIALE !!
# (Pour que le bon 'mpirun' soit trouvé)
#
module purge
module cmake
module load gcc
module load openmpi


echo "Job $SLURM_JOB_ID: Lancement $METHODE $TEST_TYPE"
echo "  - Processus demandés (ntasks): ${P_PROCS}"
echo "  - Nœuds: $SLURM_JOB_NODELIST"
module list

# --- 4. Exécution ---
MAX_ITER=50000
TOLERANCE=1e-6

# =================================================================
# MODIFICATION CLÉ : 'mpirun' est appelé sans '-n' ou '-np'.
# Il va lire ${P_PROCS} depuis la variable $SLURM_NTASKS.
CMD="${EXE_PATH} ${NX} ${NY} ${MAX_ITER} ${TOLERANCE}"
# =================================================================

# Pour un code séquentiel (si P_PROCS=1)
if [ $P_PROCS -eq 1 ]; then
    echo "  - Mode séquentiel détecté (p=1)"
    OUTPUT=$($CMD 2>&1)
else
# Pour un code MPI
    echo "  - Mode MPI détecté (p > 1)"
    OUTPUT=$(mpirun $CMD 2>&1)
fi

# --- 5. Parsing et écriture CSV ---
echo "Methode,Version,Test_Type,Nx,Ny,Processus_p,Temps_s,Iterations,Erreur_Finale" > $OUTPUT_FILE

TIME_S=$(echo "$OUTPUT" | grep -o 'Temps: [0-9.]*' | awk '{print $2}')
[ -z "$TIME_S" ] && TIME_S="ERREUR_TEMPS"

ITERS=$(echo "$OUTPUT" | grep -o 'Iterations: [0-9]*' | awk '{print $2}')
[ -z "$ITERS" ] && ITERS="NA"

ERROR=$(echo "$OUTPUT" | grep -o 'Error: [0-9.]*[eE]\{0,1\}[-+]\{0,1\}[0-9]*' | awk '{print $2}')
[ -z "$ERROR" ] && ERROR="NA"

echo "$METHODE,Parallele,${TEST_TYPE},${NX},${NY},${P_PROCS},${TIME_S},${ITERS},${ERROR}" >> $OUTPUT_FILE

echo "Test $SLURM_JOB_ID terminé. Résultat dans $OUTPUT_FILE"