
    function addNodes(data, newNodes, newEdges) {
        for (let i = 0; i < newNodes.length; i++) {
            data['nodes'].update(newNodes[i])
        }
        for (let i = 0; i < newEdges.length; i++) {
            data['edges'].update(newEdges[i])
        }
    }

    function remove_all_nodes(data) {
        data['nodes'].clear()
        data['edges'].clear()
    }

    function removeNodes(data, nodesToRemove) {

        data['nodes'].remove(nodesToRemove)

        var edges = data['edges'].getIds({
          filter: function(item) {
            return (nodesToRemove.indexOf(item.to)   !== -1) ||
                   (nodesToRemove.indexOf(item.from) !== -1);
            }
          });

        data['edges'].remove(edges)
    }
