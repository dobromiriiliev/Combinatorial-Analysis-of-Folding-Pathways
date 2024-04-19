class Residue:
    def __init__(self, type):
        self.type = type  # Type of amino acid
        self.native = False  # Indicates if the residue is in the native state
        # Additional properties
        self.hydrophobic = self.set_hydrophobicity()
        self.charge = self.set_charge()

    def set_hydrophobicity(self):
        # Example hydrophobicity values
        hydrophobic = {'I': True, 'V': True, 'L': True, 'M': True, 'A': False}
        return hydrophobic.get(self.type, False)

    def set_charge(self):
        # Example charge values
        charge = {'K': +1, 'R': +1, 'D': -1, 'E': -1}
        return charge.get(self.type, 0)

class Protein:
    def __init__(self, sequence, energy):
        self.chain = [Residue(x) for x in sequence]
        self.interaction_energy = energy

    def calculate_interaction_energy(self):
        total_energy = 0
        length = len(self.chain)
        for i in range(length):
            # Energy contribution from being in native state
            if self.chain[i].native:
                total_energy -= self.interaction_energy

            # Adjacent interactions
            if i < length - 1 and self.chain[i].hydrophobic and self.chain[i + 1].hydrophobic:
                total_energy -= 0.2  # Example energy value for hydrophobic interactions

            # Long-range interactions example (arbitrary distance set to 4 residues apart)
            if i < length - 4:
                if self.chain[i].charge != 0 and self.chain[i + 4].charge != 0:
                    # Electrostatic interactions considering charges
                    total_energy += self.chain[i].charge * self.chain[i + 4].charge * 0.1

        return total_energy

if __name__ == "__main__":
    protein_sequence = "MLKIADEV"
    protein = Protein(protein_sequence, 0.5)  # Interaction energy for native state
    interaction_energy = protein.calculate_interaction_energy()
    print(f"Total interaction energy: {interaction_energy}")
