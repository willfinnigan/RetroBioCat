from retrobiocat_web.retro.generation.load import load_rule_yamls, load_from_mongo
from pathlib import Path

yaml_path = str(Path(__file__).parents[3]) + '/data/rxn_yaml/reaction_rules.yaml'

class RetroBioCat_Reactions():
    def __init__(self, mode='mongo', yaml_path=yaml_path, include_experimental=False, include_two_step=False):
        self.rxns_strings = None
        self.rules_by_type = None
        self.mode = mode
        self.include_experimental=include_experimental

        if mode=='yaml':
            yaml_dict = load_rule_yamls.load_yamls(yaml_path)
            self.rxns = load_rule_yamls.load_rxns(yaml_dict)
            self.reactions, self.enzymes, self.reaction_enzyme_map = load_rule_yamls.load_reactions_and_enzymes(yaml_dict)
            self.reactionEnzymeCofactorDict = load_rule_yamls.load_cofactors(yaml_dict)

        elif mode=='mongo':
            query_result = load_from_mongo.get_reactions(include_experimental=self.include_experimental, include_two_step=include_two_step)
            self.rxns = load_from_mongo.load_rxns(query_result)
            self.reactions, self.enzymes, self.reaction_enzyme_map = load_from_mongo.load_reactions_and_enzymes(query_result)
            self.reactionEnzymeCofactorDict = load_from_mongo.load_cofactors(query_result)
        else:
            print(f'WARNING NO RULES LOADED - MODE {mode} NOT RECOGNISED')

    def load_additional_info(self):
        if self.mode=='yaml':
            yaml_dict = load_rule_yamls.load_yamls(yaml_path)
            self.rxns_strings = load_rule_yamls.load_rxn_strings(yaml_dict)
            self.rules_by_type = load_rule_yamls.load_rules_by_type(yaml_dict)
        elif self.mode=='mongo':
            query_result = load_from_mongo.get_reactions()
            self.rxns_strings = load_from_mongo.load_rxn_strings(query_result )
            self.rules_by_type = load_from_mongo.load_rules_by_type(query_result )

if __name__ == '__main__':
    import time
    for i in range(10):
        t0 = time.time()
        rxns_class = RetroBioCat_Reactions()
        t1 = time.time()
        print(t1-t0)