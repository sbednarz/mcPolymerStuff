#!/bin/bash -l

#SBATCH --mem=79000

## Liczba węzłów
#SBATCH -N 1

## Ilość zadań na węzeł
#SBATCH --ntasks-per-node=1

## Maksymalny czas trwania zlecenia (format HH:MM:SS)
#SBATCH --time=01:00:00 

## Nazwa grantu do rozliczenia zużycia zasobów
#SBATCH -A sbedg03

## Specyfikacja partycji
#SBATCH -p plgrid-short

## Plik ze standardowym wyjściem
#SBATCH --output="output.out"

## Plik ze standardowym wyjściem błędów
#SBATCH --error="error.out"
 
 
## przejscie do katalogu z ktorego wywolany zostal sbatch
cd $SLURM_SUBMIT_DIR

srun echo "Job ID: $SLURM_JOB_ID `hostname` `date`"
srun --mem-per-cpu=79000 mcPolymer experiment.tcl
srun date
srun echo "Happy end."
