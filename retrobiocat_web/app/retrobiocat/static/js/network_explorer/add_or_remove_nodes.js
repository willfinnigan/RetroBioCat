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



function step(params, network_data) {
    if (window.pause_add_nodes === false) {
        pause_on()
        $.post(Flask.url_for("retrobiocat.step"), {
            smiles: params.nodes[0],
            x: params.pointer.canvas.x,
            y: params.pointer.canvas.y,
            max_reactions: document.getElementById("max_pathways").value,
            mode: document.getElementById("Retrorules").checked

        }).done(function(response_data) {
            if (response_data.result.mode === 'default') {
                let new_nodes = response_data.result.nodes
                let new_edges = response_data.result.edges
                let nodesToRemove = response_data.result.to_delete
                removeNodes(network_data, nodesToRemove)
                addNodes(network_data, new_nodes, new_edges)
                pause_off()
            } else {
                get_retrorules_Status(response_data.result.response.data.task_id)
            }
        });
    }
}


function get_retrorules_Status(taskID) {
    var network_data = data
    $.get('_retrorules_step_status/' + taskID, {
    }).done(function(response) {
        const taskStatus = response.data.task_status;
        const taskProgress = response.data.task_progress;

        document.getElementById("rr_que").innerHTML = "queuing";

        if (taskStatus === 'finished') {
            document.getElementById("rr_que").innerHTML = "done";
            let new_nodes = response.data.task_result.nodes
            let new_edges = response.data.task_result.edges
            window.graph_dict = response.data.task_result.graph_dict
            window.attr_dict = response.data.task_result.attr_dict
            addNodes(network_data, new_nodes, new_edges)
        } else if (taskStatus === 'failed') {
            document.getElementById("rr_que").innerHTML = "error";
            window.pause_add_nodes = false;
            document.getElementById("interaction").innerHTML = "";
            return false;
        } else {
            if (taskStatus === 'started') {
                document.getElementById("rr_que").innerHTML = "started";
            }
            setTimeout(function() {
                get_retrorules_Status(response.data.task_id);
                }, 1000);
        }
    })
}