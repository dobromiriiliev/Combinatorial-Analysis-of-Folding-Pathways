import os
import pandas as pd

class Protein:
    class Residue:
        def __init__(self):
            self.native = False  # Indicates if the residue is in the native state

    def __init__(self, chain_size, energy, temp, barrier):
        self.chain = [self.Residue() for _ in range(chain_size)]  # Initialize a chain of residues
        self.size = chain_size  # Size of the protein chain
        self.interaction_energy = energy  # Interaction energy parameter
        self.temperature = temp  # Temperature for kinetics modeling
        self.barrier_height = barrier  # Barrier height for kinetics modeling

    def set_sequence(self, sequence):
        if len(sequence) != self.size:
            print("Error: Sequence length does not match chain size.")
            return False

        for i in range(self.size):
            self.chain[i].native = sequence[i] in ['M', 'L', 'I', 'V']  # Assume specific residues are native
        return True

    def calculate_free_energy_profile(self):
        # Calculate energy based on whether residues are in a native state
        return [-self.interaction_energy if residue.native else 0.0 for residue in self.chain]

    def write_free_energy_profile_to_file(self, filename, free_energy_profile):
        try:
            with open(filename, 'w') as output_file:
                output_file.write("Position,Free Energy\n")
                for i, energy in enumerate(free_energy_profile):
                    output_file.write(f"{i},{energy}\n")
            print(f"Free energy profile has been written to {filename}.")
        except IOError as e:
            print(f"Error writing to file: {e}")

if __name__ == "__main__":
    interaction_energy = 0.5  # Predefined interaction energy
    temperature = 300.0  # Predefined temperature
    barrier_height = 10.0  # Predefined barrier height

    # Read input from CSV file
    csv_path = "/Users/dobromiriliev/Documents/GitHub/CombinatoricsPathways/final_unip_seqs_dict.csv"
    dataset = pd.read_csv(csv_path)

    for index, entry in dataset.iterrows():
        protein = Protein(len(entry['Sequence']), interaction_energy, temperature, barrier_height)
        if protein.set_sequence(entry['Sequence']):
            free_energy_profile = protein.calculate_free_energy_profile()

            # Construct the output file path
            output_filename = os.path.join("/Users/dobromiriliev/Documents/GitHub/CombinatoricsPathways/FreeEnergyProfiles", f"{entry['Entry']}_free_energy_profile.csv")
            print(f"Output filename: {output_filename}")

            protein.write_free_energy_profile_to_file(output_filename, free_energy_profile)
