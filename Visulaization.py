import pandas as pd
import os
import math
import networkx as nx
import matplotlib.pyplot as plt

class Protein:
    def __init__(self, chain_size, energy, temp, barrier):
        self.size = chain_size
        self.interaction_energy = energy
        self.temperature = temp
        self.barrier_height = barrier
        self.chain = [False] * chain_size  # Chain of residues, initialized to False

    def set_sequence(self, sequence):
        if len(sequence) != self.size:
            print("Error: Sequence length does not match chain size.")
            return
        for i in range(self.size):
            self.chain[i] = sequence[i] in ['M', 'L', 'I', 'V']  # Example: Specific residues contributing to interaction

    def calculate_free_energy_profile(self):
        # Negative energy for active residues, 0 for others
        return [-self.interaction_energy if residue else 0.0 for residue in self.chain]

    def model_kinetics(self):
        folding_time = math.exp(self.barrier_height / self.temperature) 
        print(f"Estimated folding time: {folding_time} seconds")

def read_dataset(filename):
    data = pd.read_csv(filename)
    return [(row['Entry'], row['Sequence']) for index, row in data.iterrows()]

def build_combinatorial_assembly_graph(dataset):
    G = nx.Graph()
    for entry in dataset:
        G.add_node(entry[0])  # Add protein entry as a node
    # Adding edges arbitrarily; in practice, this could be based on some criteria
    for i in range(len(dataset)):
        for j in range(i + 1, len(dataset)):
            G.add_edge(dataset[i][0], dataset[j][0])
    return G

def visualize_graph(G):
    nx.draw(G, with_labels=True, font_weight='bold')
    plt.show()

def main():
    # Paths and parameters
    csv_path = "/Users/dobromiriliev/Documents/GitHub/CombinatoricsPathways/final_unip_seqs_dict.csv"
    interaction_energy = 0.5
    temperature = 300
    barrier_height = 10

    dataset = read_dataset(csv_path)
    dataset = dataset[:8000]  # Limiting to the first 100 proteins for display
    for entry, sequence in dataset:
        protein = Protein(len(sequence), interaction_energy, temperature, barrier_height)
        protein.set_sequence(sequence)
        print(f"Free energy profile for {entry}: ", protein.calculate_free_energy_profile())
        protein.model_kinetics()

    G = build_combinatorial_assembly_graph(dataset[:10])  # Visualize only a small part for performance
    visualize_graph(G)

if __name__ == "__main__":
    main()
