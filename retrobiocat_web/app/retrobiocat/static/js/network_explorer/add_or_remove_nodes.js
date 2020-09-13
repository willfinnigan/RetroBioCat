window.pause_add_nodes = false;

function pause_off() {
    window.pause_add_nodes = false;
    document.getElementById("interaction").innerHTML = " ";
}

function pause_on() {
    window.pause_add_nodes = true;
    document.getElementById("interaction").innerHTML = "waiting for server";
}

document.getElementById("reset").addEventListener("click", function () {
    window.pause_add_nodes = false;
    document.getElementById("interaction").innerHTML = "";
});



function default_step(params, network_data) {
    if (window.pause_add_nodes === false) {
        pause_on()
        $.post(Flask.url_for("retrobiocat.step"), {
            smiles: params.nodes[0],
            x: params.pointer.canvas.x,
            y: params.pointer.canvas.y,
            max_reactions: document.getElementById("max_pathways").value,

        }).done(function(response_data) {
            if (response_data.result.mode === 'default') {
                let new_nodes = response_data.result.nodes
                let new_edges = response_data.result.edges
                let nodesToRemove = response_data.result.to_delete
                removeNodes(network_data, nodesToRemove)
                addNodes(network_data, new_nodes, new_edges)
                pause_off()
            }
        });
    }
}

function aizynth_step(params, network_data) {
    if (window.pause_add_nodes === false) {
        pause_on()
        $.post(Flask.url_for("retrobiocat.aizynth_step"), {
            smiles: params.nodes[0],
            x: params.pointer.canvas.x,
            y: params.pointer.canvas.y,
            max_reactions: document.getElementById("max_reactions").value,

        }).done(function(response_data) {
            if (response_data.result.mode === 'default') {
                let new_nodes = response_data.result.nodes
                let new_edges = response_data.result.edges
                let nodesToRemove = response_data.result.to_delete
                removeNodes(network_data, nodesToRemove)
                addNodes(network_data, new_nodes, new_edges)
                pause_off()
            }
        });
    }
}