window.selected_node = ""

function set_reaction_info(data) {
    document.getElementById("enzyme_select").innerHTML = '';
        for (option of data.enzyme_options) {
            document.getElementById("enzyme_select").innerHTML += "<option value=" + option + ">" + option + "</option>";
            }
    document.getElementById("enzyme_select").selectedIndex = data.enzyme_select_index;
}

function get_node_info() {
    console.log('Get node info')
    document.getElementById("node_info_body").innerHTML = "<p> Please wait..  getting info.. </p>"
    $.post($SCRIPT_ROOT + '/_node_info', {
        node: window.selected_node,
        task_id: window.task_id

    }).done(function(data) {
        var node_type = data.result.type
        document.getElementById("node_info_body").innerHTML = data.result.html

        if (node_type === "Reaction Node") {
            set_reaction_info(data.result.data);
        } else if (node_type === "Substrate Node") {
        }
    })
}

function get_reaction_node_info() {
    console.log('Get reaction node info')
    $.post($SCRIPT_ROOT + '/_node_info_reaction', {
        node: window.selected_node,
        task_id: window.task_id

    }).done(function(data) {
        var node_type = data.result.type

        if (node_type === "Reaction Node") {
            set_reaction_info(data.result.data);
        } else if (node_type === "Substrate Node") {
        }
    })
}

function select_enzyme(network_data) {
    console.log('Select enzyme')
    $.post($SCRIPT_ROOT + '/_change_enzyme', {
        selected_node: window.selected_node,
        selected_enzyme:  document.getElementById("enzyme_select").value,
        task_id: window.task_id

    }).done(function(data) {
        let new_nodes = data.result.nodes
        let new_edges = data.result.edges
        addNodes(network_data, new_nodes, new_edges)
    })
}

if (document.getElementById("enzyme_select") !== null) {
    document.getElementById("enzyme_select").addEventListener("change", function () {
        console.log("Enzyme selected")
        select_enzyme(data);
    });
}

