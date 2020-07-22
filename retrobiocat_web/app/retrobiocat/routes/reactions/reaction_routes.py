from flask import render_template
from retrobiocat_web.retro.generation.load import load_rule_yamls
from retrobiocat_web.app.retrobiocat import bp
from retrobiocat_web.app.retrobiocat.functions.get_images import rxntosvg
from retrobiocat_web.retro.generation.load.rxn_class import RetroBioCat_Reactions

@bp.route('/reaction_rules')
def reaction_rules():
    rxns = RetroBioCat_Reactions()
    rxns.load_additional_info()

    images = {}
    strings = {}
    rules_by_type = rxns.rules_by_type
    for rule_type in rules_by_type:
        reaction_names = rxns.rules_by_type[rule_type]
        for name in reaction_names:
            list_rxns = rxns.rxns[name]
            images[name] = rxntosvg(list_rxns,rxnSize=(450,150))
            strings[name] = rxns.rxns_strings[name]
    return render_template('reactions/reaction_rules.html',
                            rule_types=rules_by_type,
                            rule_images=images,
                            reaction_enzymes=rxns.reaction_enzyme_map)
