from retrobiocat_web.retro.visualisation import node_information,\
    edge_information, colour_nodes, colour_edges, node_positions

import jinja2
from IPython.display import HTML
from pathlib import Path

from retrobiocat_web.retro.visualisation import rdkit_images
from retrobiocat_web.retro.generation import node_analysis

class Visualiser(object):

    def __init__(self, network, graph=None, fix_target=False):
        self.network = network

        if graph is not None:
            self.graph = graph
        else:
            self.graph = self.network.graph

        self.fix_target=fix_target

        self.nodes, self.edges = self._get_nodes_and_edges()

        self._set_node_positions()
        self._add_node_colour()
        self._add_node_information()
        self._add_edge_colour()
        self._add_edge_information()

    def html(self, height='750px', width='1000px', options=None):
        if options == None:
            options = {}

        path = str(Path(__file__).parents[0]) + '/visjs_notebook_template.html'
        with open(path) as template_html:
            content = template_html.read()

        template = jinja2.Template(content)
        rendered = template.render(nodes=self.nodes, edges=self.edges, height=height, width=width, options=options)
        html_object = HTML(rendered)

        return html_object

    def nodes_edges(self):
        return self.nodes, self.edges

    def _get_nodes_and_edges(self):
        graph = self.graph
        molSize = self.network.settings['molSize']

        nodes = []
        edges = []

        for node in node_analysis.get_substrate_nodes(graph):
            url = rdkit_images.smile_to_svg_url(node, size=molSize)

            nodes.append({'id': node,
                          'shape': 'circularImage',
                          'size': 50,
                          'borderWidth': 3,
                          'borderWidthSelected': 6,
                          'color': 'blue',
                          'image': url,
                          'label': '',
                          'title': node,
                          'type': 'substrate'})

        for node in node_analysis.get_reaction_nodes(graph):
            nodes.append({'id': node,
                          'size': 10,
                          'borderWidth': 2,
                          'borderWidthSelected': 4,
                          'color': 'darkgrey',
                          'label': graph.nodes[node]['attributes']['name'],
                          'title': node,
                          'shape': 'dot',
                          'type': 'reaction'
                          })

        for e in graph.edges:
            edges.append({'id': 'from ' + str(e[1]) + 'to ' + str(e[0]),
                          'from': e[1],
                          'to': e[0],
                          'arrows': "to",
                          'color': 'grey'})

        return nodes, edges

    def _add_node_information(self):
        self.nodes = node_information.add_substrate_info(self.nodes, self.graph)
        self.nodes = node_information.add_enzyme_info(self.nodes, self.graph)

    def _add_edge_information(self):
        if self.network.settings['display_cofactors'] == True:
            self.edges = edge_information.add_cofactors(self.graph, self.edges)

    def _add_edge_colour(self):
        if self.network.settings['colour_arrows'] == 'Complexity change':
            self.edges = colour_edges.colour_edges_by_change_in_complexity(self.graph, self.edges)

    def _add_node_colour(self):

        if self.network.settings['colour_substrates'] == 'Starting material':
            self.nodes = colour_nodes.colour_substrate_nodes_by_building_block(self.nodes, self.graph)
        elif self.network.settings['colour_substrates'] == 'Relative complexity':
            self.nodes = colour_nodes.colour_substrates_by_relative_complexity(self.nodes, self.graph,
                                                                               self.network.substrate_nodes,
                                                                               self.network.reaction_nodes)

        if self.network.settings['colour_reactions'] == 'Substrate specificity':
            self.nodes = colour_nodes.colour_reactions_by_activity(self.nodes, self.graph,
                                                                   show_negative=self.network.settings['show_negative_enzymes'])
        elif self.network.settings['colour_reactions'] == 'Complexity change':
            self.nodes = colour_nodes.colour_reactions_by_change_in_complexity(self.nodes, self.graph)

        self.nodes = colour_nodes.mark_node(self.nodes, self.network.target_smiles, 'orange')

    def _set_node_positions(self):
        if self.fix_target == True:
            node_positions.fix_node_position(self.nodes, self.network.target_smiles)


if __name__ == '__main__':
    from retrobiocat_web.retro.generation.network_generation.network import Network
    network = Network()
    network.generate('CCCCCO', 5)
    network.visualise()



