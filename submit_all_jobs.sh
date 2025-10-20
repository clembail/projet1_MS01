#!/bin/bash
#=========================================================
# SCRIPT MASTER (MAÎTRE) - LANCE TOUS LES JOBS
#
# A lancer depuis le nœud frontal:
# ./submit_all_jobs.sh
#=========================================================

echo "Lancement de la campagne de tests..."

# --- Chemins des exécutables ---
J_SEQ="./buildClem/jacobi"
J_PAR="./buildClem/jacobi_para"
G_SEQ="./buildClem/gauss-seidel"
G_PAR_RB="./buildClem/gauss-seidel_para"

# --- Paramètres de Scalabilité ---
P_VALUES="2 4 8 16 32 64 128" # Liste des 'p' à tester

# Scalabilité Forte
NX_STRONG=2048
NY_STRONG=2048

# Scalabilité Faible
NY_WEAK=2048
NX_PER_PROC_WEAK=256

# --- Soumission des jobs ---

# 0. Baseline séquentiel (p=1)
echo "Soumission des jobs 'Baseline' (p=1)..."
sbatch run_single_test.sbatch 1 $J_SEQ "Baseline" "Jacobi" $NX_STRONG $NY_STRONG
sbatch run_single_test.sbatch 1 $G_SEQ "Baseline" "Gauss-Seidel" $NX_STRONG $NY_STRONG

# 1. Scalabilité Forte
echo "Soumission des jobs 'Forte'..."
for p in $P_VALUES; do
    echo "  - Soumission Jacobi Forte pour p=$p"
    sbatch run_single_test.sbatch $p $J_PAR "Forte" "Jacobi" $NX_STRONG $NY_STRONG
    
    echo "  - Soumission Gauss-S Forte pour p=$p"
    sbatch run_single_test.sbatch $p $G_PAR_RB "Forte" "Gauss-Seidel-RN" $NX_STRONG $NY_STRONG
done

# 2. Scalabilité Faible
echo "Soumission des jobs 'Faible'..."
for p in $P_VALUES; do
    CURRENT_NX=$(( $NX_PER_PROC_WEAK * $p ))
    
    echo "  - Soumission Jacobi Faible pour p=$p (N=${CURRENT_NX}x${NY_WEAK})"
    sbatch run_single_test.sbatch $p $J_PAR "Faible" "Jacobi" $CURRENT_NX $NY_WEAK
    
    echo "  - Soumission Gauss-S Faible pour p=$p (N=${CURRENT_NX}x${NY_WEAK})"
    sbatch run_single_test.sbatch $p $G_PAR_RB "Faible" "Gauss-Seidel-RN" $CURRENT_NX $NY_WEAK
done

echo "Toutes les soumissions sont terminées. Utilisez 'squeue -u $USER' pour voir vos jobs."