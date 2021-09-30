/*jshint esversion: 6 */
/*
__author__ = 'Thomas Duigou'
__license__ = 'MIT'
*/

class PathwayHandler {

    /**
     * 
     * @param {cytoscape.js object} cy
     * @param {json structure} pathways_info 
     */
    constructor(cy, pathways_info){
        // List of the class attributes
        this.cy = cy;
        this.all_path_ids = new Set()
        this.path_to_edges = new Object()
        this.path_to_nodes = new Object()
        this.path_to_scores = new Object()
        this.pinned_path_ids = new Set()
        
        for (let path_id in pathways_info){
            if (this.all_path_ids.has(path_id)){
                console.log('path_id already referenced: ' + path_id);
            } else {
                // Path ID itselft
                this.all_path_ids.add(path_id);
                // List involved edges and nodes
                let info = pathways_info[path_id];
                this.path_to_edges[path_id] = info['edge_ids'];
                this.path_to_nodes[path_id] = info['node_ids'];
                // Extract scores
                this.path_to_scores[path_id] = info['scores'];
                // Set specific data field in the cytoscape object
                cy.elements().data('pinned', 0);
            }
        }
    }

    /**
     * Get the list of pinned pathways
     */
    get_pinned_paths(){
        return [...this.pinned_path_ids];
    }

    /**
     * Add one or more pinned pathways
     * 
     * @param {Array} path_ids
     */
    add_pinned_paths(path_ids){
        // Update path ids
        path_ids.forEach((path_id) => {
            this.pinned_path_ids.add(path_id);
        }, this);
        // Update the "pinned" status of nodes and edges
        path_ids.forEach((path_id) => {
            let node_ids = this.path_to_nodes[path_id];
            let edge_ids = this.path_to_edges[path_id];
            let element_ids = new Set([...node_ids, ...edge_ids]);
            element_ids.forEach((ele_id) => {
                let element = this.cy.getElementById(ele_id);
                element.data('pinned', 1);
            }, this);
        }, this);
        return true
    }

    /**
     * Remove one or more pinned pathways
     * 
     * @param {Array} path_ids: list of pathway IDs
     */
    remove_pinned_paths(path_ids){
        // Update path ids
        path_ids.forEach((path_id) => {
            this.pinned_path_ids.delete(path_id);
        }, this);
        // Update the "pinned" status of nodes and edges
        path_ids.forEach((path_id) => {
            let node_ids = this.path_to_nodes[path_id];
            let edge_ids = this.path_to_edges[path_id];
            let element_ids = new Set([...node_ids, ...edge_ids]);
            element_ids.forEach((element_id) => {
                if (! this.involved_in_any_pinned_path(element_id)){  // knowing we already removed the path_id from the list pinned path 
                    this.cy.getElementById(element_id).data('pinned', 0);
                }
            }, this);
        }, this);
        return true
    }

    /**
     * To know if an element is involved in at least one pinned path
     * 
     * @param {string} element_id: cytoscape element ID
     */
    involved_in_any_pinned_path(element_id){
        let path_ids = this.cy.getElementById(element_id).data('path_ids')  // returns a list
        for (let i = 0; i < path_ids.length; i++){
            if (this.pinned_path_ids.has(path_ids[i])){
                return true
            }
        }
        return false
    }

    /** 
     * Update the visibility of pinned elements
     */
    update_pinned_elements(){
        // If no any element pinned, unfade everything
        if (this.cy.elements('[pinned = 1]').length == 0){
            this.cy.elements().forEach((element) => {
                element.removeClass('faded');
            })
            return true;
        }
        // Otherwise, separate the grain from the bran
        this.cy.elements('[pinned = 1]').forEach((element) => {
            element.removeClass('faded');
        });
        this.cy.elements('[pinned = 0]').forEach((element) => {
            element.addClass('faded');
        });
    }

    /**
     * Highlight "even more" a particular set of pathways, new algo
     * 
     * @param {Array} path_ids 
     */
    highlight_pathways(path_ids=[]){
        // No pathway to highlight
        if (path_ids.length == 0){
            this.cy.edges('[highlighted = 1]').forEach((edge) => {
                edge.data('highlighted', 0);
                edge.removeClass('highlighted');
            });
            this.cy.nodes('[highlighted = 1]').forEach((node) => {
                node.data('highlighted', 0);
            });
            this.update_pinned_elements();
            return true;
        }

        // No pinned pathways
        if (this.pinned_path_ids.size == 0){
            this.cy.elements().addClass('faded');
        }

        //
        path_ids.forEach((path_id) => {
            let edge_ids = this.path_to_edges[path_id];
            let node_ids = this.path_to_nodes[path_id];
            let element_ids = new Set([...node_ids, ...edge_ids]);
            element_ids.forEach((element_id) => {
                let element = this.cy.getElementById(element_id);
                element.data('highlighted', 1);
            }, this);
        }, this);

        this.cy.edges('[highlighted = 1]').forEach((edge) => {
            edge.addClass('highlighted');
            edge.removeClass('faded');
        });
        this.cy.nodes('[highlighted = 1]').forEach((node) => {
            node.removeClass('faded');
        });
    }

    /** Colourise one pathway
     * 
     * @param {String} path_id: pathway ID
     * @param {string} colour_hex: colour in HTML hexadecical notation
     */
    colourise_one_pathwat(path_id, colour_hex){
        return true;
    }

    /**
     * Colourise a list pathways
     * 
     * @param {Array} path_ids: dictionary provided as a JSON
     * @param {String} score_label: the score label to use within available scores
     */
    colourise_pathways(path_ids, score_label='global_score'){
        let score_values = Object();
        // Shortand for all pathways
        if (path_ids == '__ALL__'){
            path_ids = [...this.all_path_ids];  // all_path_ids is a Set
        }
        // Collect and refine scores
        for (let i = 0; i < path_ids.length; i++){
            let path_id = path_ids[i];
            let score = this.path_to_scores[path_id][score_label];
            if (! isNaN(score)){
                let score_value = parseFloat(score);
                score_values[path_id] = score_value;
            }
        }
        // Sort path IDs by their values
        let items = Object.keys(score_values).map(function(key) {
            return [key, score_values[key]];
        });
        items.sort(function(first, second) {  // by inceasing order
            return first[1] - second[1];
        });
        // Set up the scale
        let list_of_values = [];
        for (let i = 0; i < items.length; i++){
            list_of_values.push(items[i][1]);
        }
        let min_score = Math.min(...list_of_values);
        let max_score = Math.max(...list_of_values);
        let colour_maker = chroma.scale(['red', 'green', 'blue']).domain([max_score, min_score]);
        // Finally colourise
        for (let i = 0; i < items.length; i++){
            // Get values
            let path_id = items[i][0];
            let score = items[i][1];
            // Get colour according to scale 
            let score_hex = colour_maker(score).hex();
            // Apply colour on edges
            let edge_ids = this.path_to_edges[path_id];
            edge_ids.forEach((edge_id) => {
                this.cy.getElementById(edge_id).style({
                    'line-color': score_hex,
                    'target-arrow-color': score_hex
                })
            }, this);
            // Apply colour on colour picker
            let colour_input = $('td.path_colour[data-path_id=' + path_id + '] > input')
            colour_input.val(score_hex);
        }
        return true;
    }

    /**
     * Reset all pathway colours to a same colour
     * 
     * @param {String} colour_hex: the hexadecimal colour to apply
     */
    reset_pathway_colours(colour_hex='#A9A9A9'){
        // Reset edge colours
        this.cy.edges().style({
            'line-color': colour_hex,
            'target-arrow-color': colour_hex
        });
        // Reset colour pickers
        $('td.path_colour > input').val(colour_hex);
        return true;
    }

};

// Utils ///////////////////////////

/**
 * Build the pathway table
 *
 * Derived from: http://jsfiddle.net/manishmmulani/7MRx6
 */
function build_pathway_table(){
    console.assert(pathways_info);
    
    // Table skeleton
    let table_base = $('<table></table>');
    
    // Build the header
    let field_names = ['Pathway', 'Show', 'Info', 'Colour', 'Score'];
    let field_classes = ['path_id_head', 'path_checkbox_head', 'path_info_head', 'path_colour_head', 'path_value_head'];  // This is needed for tablesort
    let table_row = $('<tr></tr>');
    for (let i = 0; i < field_names.length; i++){
        let value = field_names[i];
        let class_ = field_classes[i];
        table_row.append($('<th class="' + class_ + '"></th>').html(value));
    }
    table_base.append($('<thead></thead>').append(table_row));
    
    // Build the body
    let table_body = $('<tbody></tbody>');
    for (let path_id in pathways_info){
        let info = pathways_info[path_id];
        let table_row = $('<tr></tr>');
        table_row.append($('<td class="path_id" data-path_id="' + path_id + '"></td>').html(path_id));
        table_row.append($('<td class="path_checkbox"></td>').append($('<input type="checkbox" name="path_checkbox" value=' + path_id + '>')));
        table_row.append($('<td class="path_info" data-path_id="' + path_id + '"></td>'));
        table_row.append($('<td class="path_colour" data-path_id="' + path_id + '"><input type="color" name="head" value="#A9A9A9"></td>'));
        table_row.append($('<td class="path_value" data-path_id="' + path_id + '"></td>'));
        table_body.append(table_row);
    }
    table_base.append(table_body);

    // Append the content to the HTML
    $("#table_choice").append(table_base);
}

/**
 * Collect checked pathways
 */
function get_checked_pathways(){
    let selected_paths=[];
    $('input[name=path_checkbox]:checked').each(function(){
        let path_id = $(this).val();
        selected_paths.push(path_id);
    });
    return selected_paths;
}

/**
 * Get the collection of edges involved in a given path_id
 *
 * @param path_id (str): pathway ID
 * @param cy (cytoscape object): Cytoscape object
 */
function get_edges_from_path_id(path_id, cy){
    edges_col = cy.collection();
    cy.edges().forEach(function(edge, index){
        let edge_path_ids = edge.data('path_ids');
        if (share_at_least_one(edge_path_ids, [path_id])){
            edges_col = edges_col.union(edge);
        }
    });
    return edges_col;
}

/**
 * Get pinned pathways IDs
 */
function get_pinned_pathway_IDs(){
    let pinned_paths = [];
    $('td.pinned').each(function(){
        let path_id = $(this).text();
        pinned_paths.push(path_id);
    });
    return pinned_paths
}

/**
 * Put chemical info into the information panel
 */
function panel_chemical_info(node, show=false){
    if (show){
        // Collect
        let node_id = node.data('id');
        let label = node.data('label');
        let svg = node.data('svg');
        let smiles = node.data('smiles');
        let inchi = node.data('inchi');
        let inchikey = node.data('inchikey');
        if (node.data('cofactor') == 1){
            var cofactor = 'Yes';
        } else {
            var cofactor = 'No';
        }
        if (node.data('sink_chemical')){
            var insink = 'Yes';
        } else {
            var insink = 'No';
        }
        let xlinks = node.data('xlinks');
        let path_ids = node.data('path_ids');
        // Inject
        $("span.chem_info_label").html(label);
        if (inchikey == "" || inchikey == null){
            $("span.chem_info_inchikey").html("NA");
            $("span.chem_info_inchikey_search").html("");
        } else {
            $("span.chem_info_inchikey").html(inchikey);
            $("span.chem_info_inchikey_search").html('<a target="_blank" href="http://www.google.com/search?q=' + encodeURI(inchikey) + '">Look for identical structure using Google</a>');
        }
        if (inchi == ""|| inchi == null){
            $("span.chem_info_inchi").html("NA");
            $("span.chem_info_inchi_search").html("");
        } else {
            $("span.chem_info_inchi").html(inchi);
            $("span.chem_info_inchi_search").html('<a target="_blank" href="https://pubchem.ncbi.nlm.nih.gov/search/#collection=compounds&query_type=structure&query_subtype=identity&query=' + encodeURI(inchi) + '">Look for identical structure using PubChem</a>');
        }
        if (smiles == ""|| smiles == null){
            $("span.chem_info_smiles").html("NA");
            $("span.chem_info_smiles_search").html("");
        } else {
            $("span.chem_info_smiles").html(smiles);
            $("span.chem_info_smiles_search").html('<a target="_blank" href="https://pubchem.ncbi.nlm.nih.gov/search/#collection=compounds&query_type=structure&query_subtype=identity&query=' + encodeURI(smiles) + '">Look for identical structure using PubChem</a>');
        }
        $("span.chem_info_iscofactor").html(cofactor);
        $("span.chem_info_isprecursor").html(insink);
        
        // Inject SVG depiction as a background image (if any)
        if (svg !== null && svg !== ""){
            $('div.img-box').show();
            $('div.chem_info_svg').css('background-image', "url('" + svg + "')");
        } else {
            $('div.img-box').hide();
        }
        // Inject crosslinks
        $("div.chem_info_xlinks").html('');  // Reset div content
        if (xlinks.length > 0){
            for (let i = 0; i < xlinks.length; i++){
                $("div.chem_info_xlinks").append('<a target="_blank" href="' + xlinks[i]['url'] + '">' + xlinks[i]['db_name'] + ':' + xlinks[i]['entity_id'] + '</a>');
                $("div.chem_info_xlinks").append('<br/>');
            }
        } else {
            $("div.chem_info_xlinks").append('None<br/>');
        }
        // Inject path IDs
        $("div.chem_info_pathids").html('');  // Reset div content
        if (path_ids.length > 0){
            for (let i = 0; i < path_ids.length; i++){
                $("div.chem_info_pathids").append(path_ids[i] + '<br/>');
            }
        } else {
            $("div.chem_info_pathids").append('None<br/>');
        }
        // Show
        $("#panel_chemical_info").show();
    } else {
        $("#panel_chemical_info").hide();
    }
}
    
/**
 * Put reaction info into the information panel
 */
function panel_reaction_info(node, show=true){
    if (show){
        // Collect
        let node_id = node.data('id');
        let label = node.data('label')
        let rsmiles = node.data('rsmiles');
        let rule_ids = node.data('rule_ids');  // TODO: handle list of IDs
        let rxn_template_id = node.data('rxn_template_id');  // TODO: handle list of IDs
        let path_ids = node.data('path_ids');
        let ec_numbers = node.data('ec_numbers');
        let xlinks = node.data('xlinks');
        let thermo_value = node.data('thermo_dg_m_gibbs');
        let rule_score = node.data('rule_score');
        let uniprot_ids = node.data('uniprot_ids');
        // Inject 
        $("span.reaction_info_rsmiles").html(rsmiles);
        // Reaction name
        $("span.reaction_info_name").html(label);
        // Rule ID
        $("div.reaction_info_ruleids").html('');  // Reset div content
        for (let i = 0; i < rule_ids.length; i++){
            $("div.reaction_info_ruleids").append(rule_ids[i] + '<br/>');
        }
        // Reaction template ID
        $("span.reaction_info_reaction_template_id").html('');  // Reset div content
        $("span.reaction_info_reaction_template_id").html(rxn_template_id);
        // EC numbers
        $("div.reaction_info_ecnumbers").html('');  // Reset div content
        if (ec_numbers == null || ec_numbers.length == 0){
            $("div.reaction_info_ecnumbers").append('None<br/>');
        } else {
            for (let i = 0; i < ec_numbers.length; i++){
                $("div.reaction_info_ecnumbers").append(ec_numbers[i] + '<br/>');
            }
        }
        // Inject crosslinks
        $("div.reaction_info_xlinks").html('');  // Reset div content
        if (xlinks != null && xlinks.length > 0){
            for (let i = 0; i < xlinks.length; i++){
                $("div.reaction_info_xlinks").append('<a target="_blank" href="' + xlinks[i]['url'] + '">' + xlinks[i]['db_name'] + ':' + xlinks[i]['entity_id'] + '</a>');
                $("div.reaction_info_xlinks").append('<br/>');
            }
        } else {
            $("div.reaction_info_xlinks").append('None<br/>');
        }
        // Inject path IDs
        $("div.reaction_info_pathids").html('');  // Reset div content
        if (path_ids.length > 0){
            for (let i = 0; i < path_ids.length; i++){
                $("div.reaction_info_pathids").append(path_ids[i] + '<br/>');
            }
        } else {
            $("div.reaction_info_pathids").append('None<br/>');
        }
        // Thermodynamic value
        if (isNaN(thermo_value)){
            thermo_value = "NaN";
        } else {
            thermo_value = parseFloat(thermo_value).toFixed(3);
        }
        $("span.reaction_info_thermo").html(thermo_value);
        // Rule score
        if (isNaN(rule_score)){
            rule_score = "NaN";
        } else {
            rule_score = parseFloat(rule_score).toFixed(3);
        }
        $("span.reaction_info_rule_score").html(rule_score);
        // Inject UniProt IDs
        $("div.reaction_info_uniprot_crosslinks").html('');  // Reset div content
        let nb_ids = 0;
        for (let uid in uniprot_ids) {
            ++nb_ids;
            let selenzy_score = parseFloat(uniprot_ids[uid]['score']).toFixed(1)
            $("div.reaction_info_uniprot_crosslinks").append(
                    '<a href="' + get_uniprot_xlink(uid) + '">' + uid + '</a> (' + selenzy_score + ')'
                );
            $("div.reaction_info_uniprot_crosslinks").append(
                '<br/>'
                );
        }
        if (nb_ids == 0){
            $("div.reaction_info_uniprot_crosslinks").append('None<br/>');
        }
        // Selenzyme crosslink
        $("span.reaction_info_selenzyme_crosslink").html('<a target="_blank" href="http://selenzyme.synbiochem.co.uk/results?smarts=' + encodeURIComponent( rsmiles ) + '">Crosslink to Selenzyme</a>');
        // Show
        $("#panel_reaction_info").show();
    } else {
        $("#panel_reaction_info").hide();
    }
}


/**
 * Generate links to UniProt web site
 * @param {uid} UniProt ID to be used 
 */
function get_uniprot_xlink(uid){
    return 'https://www.uniprot.org/uniprot/' + uid
}


/**
 * Write some default text message on the info panel
 */
function panel_startup_info(show=true){  // node
    if (show){
        $("#panel_startup_legend").show();
    } else {
        $("#panel_startup_legend").hide();
    }
    
}

/**
 * Put pathway info into the information panel
 *
 * @param path_id (str): pathway ID
 */
function panel_pathway_info(path_id, show=true){
    if (show){
        // Collect
        let global_score = pathways_info[path_id]['scores']['global_score'];
        let rule_score = pathways_info[path_id]['scores']['rule_score'];
        let fba_value = pathways_info[path_id]['scores']['fba_target_flux'];
        let thermo_value = pathways_info[path_id]['scores']['thermo_dg_m_gibbs'];
        let nb_steps = pathways_info[path_id]['nb_steps'];
        // Refine the global score value
        if (isNaN(global_score)){
            global_score = "NaN";
        } else {
            global_score = parseFloat(global_score).toFixed(3);
        }
        // Refines thermodynamic value
        if (isNaN(thermo_value)){
            thermo_value = "NaN";
        } else {
            thermo_value = parseFloat(thermo_value).toFixed(3);
        }
        // Refines rule score
        if (isNaN(rule_score)){
            rule_score = "NaN";
        } else {
            rule_score = parseFloat(rule_score).toFixed(3);
        }
        // Refines target's flux production
        if (isNaN(fba_value)){
            fba_value = "NaN";
        } else {
            fba_value = parseFloat(fba_value).toFixed(3);
        }
        // Inject
        $("span.pathway_info_path_id").html(path_id);
        $("span.pathway_info_global_score").html(global_score);
        $("span.pathway_info_thermo").html(thermo_value);
        $("span.pathway_info_rule_score").html(rule_score);
        $("span.pathway_info_target_flux").html(fba_value);
        $("span.pathway_info_nb_steps").html(nb_steps);
        // Show
        $("#panel_pathway_info").show();
    } else {
        $("#panel_pathway_info").hide();
    }
}

/**
 * Return true if the array have at least one common items
 *
 * @param array1 (array): items
 * @param array2 (array): items
 */
function share_at_least_one(array1, array2){
    for (let i = 0; i < array1.length; i++){
        for (let j = 0; j < array2.length; j++){
            // We have  a match
            if (array1[i] == array2[j]){
                return true;
            }
        }
    }
    return false;
}

/**
 * Make labels for chemicals
 *
 * @param {Integer} max_length: string size cutoff before label truncation
 */
function make_chemical_labels(max_length=6){
    let nodes = cy.nodes().filter('[type = "chemical"]');
    for (let i = 0; i < nodes.size(); i++){
        let node = nodes[i];
        let label = node.data('label');
        if ((typeof label != 'undefined') && (label != 'None') && (label != '')){
            if (label.length > max_length){
                short_label = label.substr(0, max_length-2)+'..';
            } else {
                short_label = label;
            }
        } else {
            short_label = '';
        }
        node.data('short_label', short_label);
    }
}

/**
 * Make labels for reactions
 *
 * @param {Integer} max_length: string size cutoff before label truncation
 */
function make_reaction_labels(max_length=10){
    let nodes = cy.nodes().filter('[type = "reaction"]');
    for (let i = 0; i < nodes.size(); i++){
        let node = nodes[i];
        let label = node.data('label');
        if ((typeof label != 'undefined') && (label != 'None') && (label != '')){
            if (label.length > max_length){
                short_label = label.substr(0, max_length-2)+'..';
            } else {
                short_label = label;
            }
        } else {
            short_label = '';
        }
        node.data('short_label', short_label);
    }
}

// Live ///////////////////////////


$(function(){

    // Cytoscape object to play with all along
    var cy = window.cy = cytoscape({
        container: document.getElementById('cy'),
        motionBlur: true
    });

    // Basic stuff to do only once
    build_pathway_table();
    panel_startup_info(true);
    panel_chemical_info(null, false);
    panel_reaction_info(null, false);
    panel_pathway_info(null, false);
    init_network(true);
    annotate_hiddable_cofactors();  // Need to be done after init_network so the network is already loaded
    refresh_layout();
    show_cofactors(false);
    put_pathway_values('global_score');
    make_pathway_table_sortable();  // Should be called only after the table has been populated with values

    // Pathway Handler stuff
    window.path_handler = new PathwayHandler(cy, pathways_info);
    path_handler.colourise_pathways('__ALL__', 'global_score');

    /**
     * Initialise the network, but hide everything
     *
     * @param show_graph (bool): show the loaded network
     */
    function init_network(show_graph=true){
        // Reset the graph
        cy.json({elements: {}});
        cy.minZoom(1e-50);
        
        // Load the full network
        cy.json({elements: network['elements']});
        
        // Create node labels
        make_chemical_labels(6);
        make_reaction_labels(9);
        
        // Hide them 'by default'
        if (! show_graph){
            show_pathways(selected_paths='__NONE__');
        } else {
            $('input[name=path_checkbox]').prop('checked', true);  // Check all
        }
        
        // Once the layout is done:
        // * set the min zoom level
        // * put default info
        cy.on('layoutstop', function(e){
            cy.minZoom(cy.zoom());
        });
        
        cy.style(
            cytoscape.stylesheet()
                .selector("node[type='reaction']")
                    .css({
                        'height': 60,
                        'width': 120,
                        'background-color': 'white',
                        'border-width': 5,
                        'border-color': 'darkgray',
                        'border-style': 'solid',
                        'content': 'data(short_label)',
                        'text-valign': 'center',
                        'text-halign': 'center',
                        'text-opacity': 1,
                        'color': '#575757',
                        'font-size': '20px',
                    })
                .selector("node[type='chemical']")
                    .css({
                        'background-color': '#52be80',
                        'background-fit':'contain',
                        'shape': 'roundrectangle',
                        'width': 80,
                        'height': 80,
                        'label': 'data(short_label)',
                        'font-size': '20px',
                        // 'font-weight': 'bold',
                        'text-valign': 'top',
                        'text-halign': 'center',
                        'text-margin-y': 8,
                        'text-opacity': 1,
                        'text-background-color': 'White',
                        'text-background-opacity': 0.85,
                        'text-background-shape': 'roundrectangle',
                    })
                .selector("node[type='chemical'][?target_chemical]")
                    .css({
                        'background-color': '#B22222',
                        'border-color': '#B22222',
                    })
                .selector("node[type='chemical'][?sink_chemical]")
                    .css({
                        'background-color': '#68956D',
                        'border-color': '#68956D'
                    })
                .selector("node[type='chemical'][!target_chemical][!sink_chemical]")  // ie: intermediates
                    .css({
                        'background-color': '#235789',
                        'border-color': '#235789',
                    })
                .selector("node[type='chemical'][?svg]")  // The beauty of it: "?" will match only non null values
                    .css({
                        'background-image': 'data(svg)',
                        'background-fit': 'contain',
                        'border-width': 8,
                    })
                .selector('edge')
                    .css({
                        'curve-style': 'bezier',
                        'line-color': 'darkgray',
                        'width': '5px',
                        'target-arrow-shape': 'triangle',
                        'target-arrow-color': 'darkgray',
                        'arrow-scale' : 2
                    })                    
                .selector('.faded')
                    .css({
                        'opacity': 0.15,
                        'text-opacity': 0.25
                    })
                .selector('.highlighted')
                    .css({
                        'width': '9px'
                    })
                .selector('node:selected')
                    .css({
                        'border-width': 5,
                        'border-color': 'black'
                    })
        );
        
        cy.on('tap', 'node', function(evt){
            let node = evt.target;
            // Dump into console
            console.log(node.data());
            // Print info
            if (node.is('[type = "chemical"]')){
                panel_startup_info(false);
                panel_reaction_info(null, false);
                panel_pathway_info(null, false);
                panel_chemical_info(node, true);
            } else if (node.is('[type = "reaction"]')){
                panel_startup_info(false);
                panel_chemical_info(null, false);
                panel_pathway_info(null, false);
                panel_reaction_info(node, true);
            }
        });

        cy.on('tap', 'edge', function(evt){
            let edge = evt.target;
            console.log(edge.data());
        });
        
    }
    
    /**
     * Trigger a layout rendering
     * 
     * @param {cytoscape collection} element_collection: a collection of elements.
     */
    function render_layout(element_collection){
        // Playing with zoom to get the best fit
        cy.minZoom(1e-50);
        cy.on('layoutstop', function(e){
            cy.minZoom(cy.zoom()*0.9);  // 0.9 to enable the user dezoom a little
        });
        // Layout
        let layout = element_collection.layout({
            name: 'breadthfirst',
            roots: cy.elements("node[?target_chemical]")
        });
        layout.run();
    }
        
    /** Load a metabolic network
     *
     * Only nodes and edges involved in 'selected_paths' will be displayed.
     *
     * @param selected_paths (array or str): path IDs or special flags
     */
    function show_pathways(selected_paths='__ALL__'){
      
        if (selected_paths == '__ALL__'){
            cy.nodes().css({visibility: 'visible'});
            cy.edges().css({visibility: 'visible'});
        } else if (selected_paths == '__NONE__'){
            cy.nodes().css({visibility: 'hidden'});
            cy.edges().css({visibility: 'hidden'});
        } else {
            // Nodes
            cy.nodes().forEach(function(node, index){
                let node_paths = node.data('path_ids');
                if (share_at_least_one(node_paths, selected_paths)){
                    node.css({visibility:'visible'});
                } else {
                    node.css({visibility:'hidden'});
                }
            });
            // Edges
            cy.edges().forEach(function(edge, index){
                let edge_paths = edge.data('path_ids');
                if (share_at_least_one(edge_paths, selected_paths)){
                    edge.css({visibility:'visible'});
                } else {
                    edge.css({visibility:'hidden'});
                }
            });
        }
    }

    /**
     * Tag cofactor weither there could be hidden or not 
     *
     * If hidding cofactors lead to lonely / unconnected reactions then 
     * cofactors related to sucbh reaction are marked as not hiddable.
     * Otherwise, cofactors are marked as hiddable.
     */
    function annotate_hiddable_cofactors(){
        cy.elements('node[type = "reaction"]').forEach((rxn_node, i) => {
            // Check
            let in_not_cof = rxn_node.incomers().filter('node[?cofactor]');
            let out_not_cof = rxn_node.outgoers().filter('node[?cofactor]');
            // Decide
            let hiddable;
            if (in_not_cof.length == 0 || out_not_cof.length == 0){
                hiddable = 0;
            } else {
                hiddable = 1;
            }
            // Tag
            let in_chems = rxn_node.incomers().filter('node');
            in_chems.forEach((chem_node, j) => {
                if (
                    chem_node.data('cofactor') == 1 &&
                    chem_node.data('hiddable_cofactor') != 0
                ){
                    chem_node.data('hiddable_cofactor', hiddable);
                }
            });
            let out_chems = rxn_node.outgoers().filter('node');
            out_chems.forEach((chem_node, j) => {
                if (
                    chem_node.data('cofactor' == 1) &&
                    chem_node.data('hiddable_cofactor') != 0
                ){
                    chem_node.data('hiddable_cofactor', hiddable);
                }
            });
        });
    }

    /** Handle cofactor display
     *
     * Hide of show all nodes annotated as cofactor
     *
     * @param show (bool): will show cofactors if true
     */
    function show_cofactors(show=true){
        if (show){
            cy.elements().style("display", "element");
        } else {
            cy.elements('node[?cofactor][?hiddable_cofactor]').style("display", "none");
        }
        refresh_layout();
    }

    /**
     * Make the pathway table sortable
     */
    function make_pathway_table_sortable(){
        $("#table_choice > table").tablesorter({
            theme : 'default',
            sortList: [[4,1],[0,0]],  // Sort on the fourth column (descending) and then on the first column (ascending order)
            headers : {  // Disable sorting for these columns
                '.path_checkbox_head, .path_info_head, .path_colour_head': {
                    sorter: false
                }
            }
        });
    }
    
    /**
     * Refresh layout according to visible nodes
     */
    function refresh_layout(){
        render_layout(cy.elements().not(':hidden'));
    }
    
    // When a pathway is checked
    $("input[name=path_checkbox]").change(function(){
        selected_paths = get_checked_pathways();
        show_pathways(selected_paths);
    });
    
    /** 
     * Pathway visibility is updated when a pathway label is hovered
     * 
     * Note: the hover CSS is handled in the CSS file.
     * Node: some vocabulary precisions, pinned stands for path ID locked "on",
     *      while highlighted stands for the path ID currently hovered
     */
    $("td.path_id").hover(function(){
        let path_id = $(this).data('path_id');
        path_handler.highlight_pathways([path_id]);

    }, function(){
        let path_id = $(this).data('path_id');
        path_handler.highlight_pathways([]);
    });
    
    /**
     * Pathway are pinned on click
     */
    $("td.path_id").click(function(){
        let path_id = $(this).data('path_id');
        // Removing
        if ($(this).hasClass('pinned')){
            $(this).removeClass('pinned');
            path_handler.remove_pinned_paths([path_id]);
            path_handler.update_pinned_elements();
        // Adding
        } else {
            path_handler.add_pinned_paths([path_id]);
            path_handler.update_pinned_elements();
            $(this).addClass('pinned');
        }
    });
    
    // When a pathway "info" is clicked
    $("td.path_info").click(function(){
        path_id = $(this).data('path_id');
        panel_startup_info(false);
        panel_chemical_info(null, false);
        panel_reaction_info(null, false);
        panel_pathway_info(path_id, true);
    });
        
    // Pathways selection
    $('#hide_all_pathways_button').on('click', function(event){
        show_pathways(selected_paths='__NONE__');  // Hide all
        $('input[name=path_checkbox]').prop('checked', false);  // Uncheck all
    });
    $('#view_all_pathways_button').on('click', function(event){
        show_pathways(selected_paths='__ALL__');  // Show all
        $('input[name=path_checkbox]').prop('checked', true);  // Check all
    });
    $('#redraw_pathways_button').on('click', function(event){
        refresh_layout();
    });
    
    // Cofactors handling
    $('#show_cofactors_button').on('click', function(event){
        show_cofactors(true);
        // Update visible pathways to update their cofactor nodes visibility
        selected_paths = get_checked_pathways();
        show_pathways(selected_paths);
        // Update hilighted pathways to update their cofactor nodes status
        path_handler.update_pinned_elements();
    });
    $('#remove_cofactors_button').on('click', function(event){
        show_cofactors(false);
    });
    
    // Manual colour handling
    colour_pickers = document.querySelectorAll(".path_colour");
    for (let i = 0; i < colour_pickers.length; i++){
        colour_pickers[i].addEventListener("input", live_update_colour, false);
    }
    
    /**
     * Set the colour of all edges involved in a pathway
     *
     * @param event: event related to a pathway
     */
    function live_update_colour(event) {
        let path_id = $(this).data('path_id');
        edges = get_edges_from_path_id(path_id, cy);
        edges.style({
            'line-color': event.target.value,
            'target-arrow-color': event.target.value
        });
    }

    /**
     * 
     * Fill table values
     * 
     * @param score_label (str): the score label to use within the path info
     */
    function put_pathway_values(score_label='global_score'){
        for (let path_id in pathways_info){
            // Collect the value
            let score = pathways_info[path_id]['scores'][score_label];
            if (! isNaN(score)){
                score = parseFloat(score).toFixed(3);
            } else {
                score = 'NaN';
            }
            // Push it into the pathway table
            let path_td = $('td.path_value[data-path_id=' + path_id + ']');
            path_td.html(score);
            
        }
    }

});
