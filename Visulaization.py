import pandas as pd
import math
import networkx as nx
import plotly.graph_objects as go
from plotly.offline import plot
import plotly.express as px

class Protein:
    def __init__(self, entry, sequence, interaction_energy, temperature, barrier_height):
        self.entry = entry
        self.sequence = sequence
        self.interaction_energy = interaction_energy
        self.temperature = temperature
        self.barrier_height = barrier_height
        self.chain = [residue in ['M', 'L', 'I', 'V'] for residue in sequence]

    def calculate_free_energy_profile(self):
        return [-self.interaction_energy if residue else 0.0 for residue in self.chain]

    def model_kinetics(self):
        folding_time = math.exp(self.barrier_height / self.temperature)
        return folding_time

def read_dataset(filename):
    return pd.read_csv(filename)

def build_combinatorial_assembly_graph(proteins):
    G = nx.Graph()
    for protein in proteins:
        G.add_node(protein.entry, sequence=protein.sequence, energy=protein.interaction_energy)
        for other in proteins:
            if protein != other:
                G.add_edge(protein.entry, other.entry)
    return G

def visualize_graph_with_plotly(G):
    pos = nx.spring_layout(G)
    edge_x = []
    edge_y = []
    for edge in G.edges():
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]
        edge_x.extend([x0, x1, None])
        edge_y.extend([y0, y1, None])

    edge_trace = go.Scatter(
        x=edge_x, y=edge_y,
        line=dict(width=0.5, color='#888'),
        hoverinfo='none',
        mode='lines')

    node_x = []
    node_y = []
    text = []
    colors = []  # List to store unique colors for each node
    color_palette = px.colors.qualitative.Plotly  # Using Plotly's qualitative color palette

    for i, node in enumerate(G.nodes()):
        x, y = pos[node]
        node_x.append(x)
        node_y.append(y)
        text.append(f"{node}: {G.nodes[node]['energy']}")
        colors.append(color_palette[i % len(color_palette)])  # Cycle through color palette

    node_trace = go.Scatter(
        x=node_x, y=node_y,
        mode='markers',
        hoverinfo='text',
        marker=dict(
            size=10,
            color=colors,  # Apply unique color for each node
            line=dict(width=2, color='DarkSlateGrey')
        ),
        text=text
    )

    fig = go.Figure(data=[edge_trace, node_trace],
                    layout=go.Layout(
                        showlegend=False,
                        hovermode='closest',
                        margin=dict(b=0, l=0, r=0, t=0),
                        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False))
                    )
    plot(fig, filename='network.html')

def main():
    csv_path = "final_unip_seqs_dict.csv"
    interaction_energy = 0.5
    temperature = 300
    barrier_height = 10

    data = read_dataset(csv_path)
    proteins = [Protein(row['Entry'], row['Sequence'], interaction_energy, temperature, barrier_height) for index, row in data.iterrows()]
    limited_proteins = proteins[:10]  # Limit to first 10 proteins for visualization

    for protein in limited_proteins:
        energy_profile = protein.calculate_free_energy_profile()
        folding_time = protein.model_kinetics()
        print(f"Free energy profile for {protein.entry}: {energy_profile}")
        print(f"Estimated folding time for {protein.entry}: {folding_time} seconds")

    G = build_combinatorial_assembly_graph(limited_proteins)
    visualize_graph_with_plotly(G)

if __name__ == "__main__":
    main()
