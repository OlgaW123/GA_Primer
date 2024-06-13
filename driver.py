import subprocess
import os

def run_experiment(pe_values, pm_values):
    for Pe, Pm in zip(pe_values, pm_values):
        try:
            # Construct the command to run the original script with the current Pe and Pm values
            command = [
                "python3", "code.py",
                f"--Pe={Pe}",
                f"--Pm={Pm}"
            ]

            # Run the command
            subprocess.run(command, check=True)

            print(f"Completed for Pm={Pm} and Pe={Pe}")

        except subprocess.CalledProcessError as e:
            print(f"Error for Pm={Pm} and Pe={Pe}: {e}")

# Define ranges for Pe and Pm
pe_values = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5]
pm_values = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]

# Call run_experiment function with defined ranges
run_experiment(pe_values, pm_values)

if os.path.exists("hg38_cut.fa"):
    os.remove("hg38_cut.fa")
if os.path.exists("initial_population.txt"):
    os.remove("initial_population.txt")
