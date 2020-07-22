

def fix_node_position(list_nodes, target_node):
    for i, node in enumerate(list_nodes):
        if list_nodes[i]['id'] == target_node:
            list_nodes[i]['fixed'] = True

    return list_nodes

