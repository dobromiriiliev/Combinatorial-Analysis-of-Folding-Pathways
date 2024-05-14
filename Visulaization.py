import pandas as pd
import math
import networkx as nx
import plotly.graph_objects as go
from plotly.offline import plot
import plotly.express as px
import nglview as nv
import ipywidgets

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
    pos = nx.spring_layout(G, k=0.3, iterations=100)  
    edge_x = []
    edge_y = []
    for edge in G.edges():
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]
        edge_x.extend([x0, x1, None])
        edge_y.extend([y0, y1, None])

    edge_trace = go.Scatter(
        x=edge_x, y=edge_y,
        line=dict(width=1, color='#888', dash='solid'),
        opacity=0.7,
        hoverinfo='none',
        mode='lines')

    node_x = []
    node_y = []
    texts = []
    colors = []
    node_size = []
    color_palette = px.colors.qualitative.Plotly

    for i, node in enumerate(G.nodes()):
        x, y = pos[node]
        node_x.append(x)
        node_y.append(y)
        node_info = f"{node}<br>Energy: {G.nodes[node]['energy']}<br>Size: {G.nodes[node]['size']}"
        texts.append(node_info)
        colors.append(color_palette[i % len(color_palette)])
        node_size.append(10 + G.nodes[node]['size'] / 5)  

    node_trace = go.Scatter(
        x=node_x, y=node_y,
        mode='markers+text',
        text=[node for node in G.nodes()],
        textposition="top center",
        hoverinfo='text',
        marker=dict(
            size=node_size,
            color=colors,
            line=dict(width=2, color='DarkSlateGrey')
        ),
        textfont=dict(
            family="Arial",
            size=10,
            color='black'),
        hovertext=texts
    )

    fig = go.Figure(data=[edge_trace, node_trace],
                    layout=go.Layout(
                        title='Protein-Protein Interaction Network',
                        titlefont_size=20,
                        showlegend=False,
                        hovermode='closest',
                        margin=dict(b=20, l=5, r=5, t=40),
                        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                        paper_bgcolor='white',
                        plot_bgcolor='white',
                    ))
    plot(fig, filename='network.html')

def visualize_proteins_with_nglview(proteins):
    views = []
    for protein in proteins:
        structure = nv.TextStructure(protein.sequence)
        view = nv.show_text(structure)
        view.add_cartoon(color='sstruc')
        views.append(view)

    tab = ipywidgets.Tab()
    children = [view._remote_repr_jupyter_() for view in views]
    tab.children = children
    for i, protein in enumerate(proteins):
        tab.set_title(i, protein.entry)
    display(tab)

def main():
    csv_path = "/Users/dobromiriliev/Documents/GitHub/CombinatoricsPathways/final_unip_seqs_dict.csv"
    interaction_energy = 0.5
    temperature = 300
    barrier_height = 10

    data = read_dataset(csv_path)
    proteins = [Protein(row['Entry'], row['Sequence'], interaction_energy, temperature, barrier_height) for index, row in data.iterrows()]
    limited_proteins = proteins[:10]

    for protein in limited_proteins:
        energy_profile = protein.calculate_free_energy_profile()
        folding_time = protein.model_kinetics()
        print(f"Free energy profile for {protein.entry}: {energy_profile}")
        print(f"Estimated folding time for {protein.entry}: {folding_time} seconds")

    G = build_combinatorial_assembly_graph(limited_proteins)
    visualize_graph_with_plotly(G)
    visualize_proteins_with_nglview(limited_proteins)

if __name__ == "__main__":
    main()
