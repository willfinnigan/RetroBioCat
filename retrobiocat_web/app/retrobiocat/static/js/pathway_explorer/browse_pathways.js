document.getElementById("left").addEventListener("click", function() {
    window.pathway_num --
    if (window.pathway_num <1) {
        window.pathway_num = 1
    }
    document.getElementById("pathway_number_select").value = window.pathway_num;
    window.pathway_varient = 1
    change_pathway(data)
});

document.getElementById("right").addEventListener("click", function() {
    window.pathway_num ++
    if (window.pathway_num > window.max_pathways) {
        window.pathway_num = window.max_pathways
    }
    document.getElementById("pathway_number_select").value = window.pathway_num;
    window.pathway_varient = 1
    change_pathway(data)
});


document.getElementById("left_varient").addEventListener("click", function() {
    window.pathway_varient --
    if (window.pathway_varient < 1) {
        window.pathway_varient = 1
    } else {
        change_pathway(data)
    }
});

document.getElementById("right_varient").addEventListener("click", function() {
    window.pathway_varient ++
    if (window.pathway_varient > window.max_varient) {
        window.pathway_varient = window.max_varient
    } else {
        change_pathway(data)
    }
});

function change_pathway() {
    window.pathway_num = document.getElementById("pathway_number_select").value
    document.getElementById("pathway_varient_select").innerHTML = "" + window.pathway_varient + ' of ' + window.max_varient;
    $.post($SCRIPT_ROOT + '/_next_pathway', {
        pathway_num: window.pathway_num,
        varient_num: window.pathway_varient,
        task_id: window.task_id
    }).done(function(response_data) {
        let new_nodes = response_data.result.nodes
        let new_edges = response_data.result.edges
        window.max_varient = response_data.result.max_varient
        document.getElementById("pathway_varient_select").innerHTML = "" + window.pathway_varient + ' of ' + window.max_varient;
        remove_all_nodes(data)
        addNodes(data, new_nodes, new_edges)
        network.stabilize()
    })
}

document.getElementById("pathway_number_select").onchange = function()
    {
        change_pathway(data)
    }