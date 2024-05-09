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
        G.add_node(protein.entry, sequence=protein.sequence, energy=protein.interaction_energy, size=len(protein.sequence))
        for other in proteins:
            if protein != other:
                G.add_edge(protein.entry, other.entry)
    return G

def visualize_graph_with_plotly(G):
    pos = nx.spring_layout(G, k=0.1, iterations=50)  # Improved layout
    edge_x = []
    edge_y = []
    for edge in G.edges():
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]
        edge_x.extend([x0, x1, None])
        edge_y.extend([y0, y1, None])

    edge_trace = go.Scatter(
        x=edge_x, y=edge_y,
        line=dict(width=0.5, color='Grey'),
        hoverinfo='none',
        mode='lines')

    node_x = []
    node_y = []
    texts = []
    colors = []  # List to store unique colors for each node
    node_size = []
    color_palette = px.colors.qualitative.D3  # Professional color palette

    for i, node in enumerate(G.nodes()):
        x, y = pos[node]
        node_x.append(x)
        node_y.append(y)
        node_info = f"{node}<br>Energy: {G.nodes[node]['energy']}<br>Size: {G.nodes[node]['size']}"
        texts.append(node_info)
        colors.append(color_palette[i % len(color_palette)])  # Cycle through color palette
        node_size.append(10 + G.nodes[node]['size'] / 5)  # Dynamic node size

    node_trace = go.Scatter(
        x=node_x, y=node_y,
        mode='markers+text',  # Add text to the nodes
        text=[node for node in G.nodes()],  # Labels are the protein entries
        textposition="middle center",  # Position text in the center of the nodes
        hoverinfo='text',
        marker=dict(
            size=node_size,
            color=colors,
            line=dict(width=2, color='DarkSlateGrey')
        ),
        textfont=dict(
            family="Arial",
            size=12,
            color='black'),
        hovertext=texts
    )

    fig = go.Figure(data=[edge_trace, node_trace],
                    layout=go.Layout(
                        title='Protein Interaction Graph',
                        titlefont_size=16,
                        showlegend=False,
                        hovermode='closest',
                        margin=dict(b=20, l=5, r=5, t=40),
                        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                        paper_bgcolor='white',
                        plot_bgcolor='white',
                    ))
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
