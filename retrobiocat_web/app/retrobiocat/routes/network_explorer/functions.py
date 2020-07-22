def add_new(original, new):
    ids = []
    unique = []
    for item in new:
        unique.append(item)
        ids.append(item['id'])

    for item in original:
        if item['id'] not in ids:
            ids.append(item['id'])
            unique.append(item)

    return unique


def get_ids_to_delete(nodes_to_delete):
    id_to_delete = []
    for node in nodes_to_delete:
        id_to_delete.append(node['id'])
    return id_to_delete

def delete_nodes_and_edges(nodes_to_delete, nodes, edges):
    new_nodes = [n for n in nodes if n['id'] not in nodes_to_delete]

    new_edges = []
    for edge in edges:
        if edge['to'] not in nodes_to_delete:
            if edge['from'] not in nodes_to_delete:
                new_edges.append(edge)

    return new_nodes, new_edges

