import requests
import xmltodict
import collections

def strip_uniref_name(name):
    name = name.replace('UniRef50_', '')
    name = name.replace('UniRef90_', '')
    name = name.replace('UniRef100_', '')
    return name

class UniProt_Parser(object):

    def __init__(self):
        self.uniprot_url = "https://www.uniprot.org/uniprot"
        self.xml_dict = {}

    def get_pfams(self):
        """ Return a list of the pfam id's in a uniprot or uniref representative sequence"""

        p_data = self.xml_dict['uniprot']['entry']['dbReference']
        pfams = dict()
        for list_item in p_data:
            if list_item.get("@type", '') == 'Pfam':
                pfam_id = ''
                pfam_name = ''
                if '@id' in list_item:
                    pfam_id = list_item['@id']
                    #print(list_item)
                if 'property' in list_item:
                    for property_item in list_item['property']:
                        if property_item.get('@type', '') == 'entry name':
                            if '@value' in property_item:
                                pfam_name = property_item['@value']

                if pfam_id != '':
                    pfams[pfam_id] = pfam_name

        return pfams

    def load_xml(self, name):
        url = f"{self.uniprot_url}/{name}.xml"
        req = requests.get(url)

        if req.status_code in [200]:
            xml = req.text
        else:
            xml = None

        self.xml_dict = xmltodict.parse(xml)

class UniRef_Parser(object):

    def __init__(self, log_level=0):
        self.uniref_url = "https://www.uniprot.org/uniref"
        self.xml_dict = {}
        self.log_level = log_level

    def get_cluster_name(self):
        return self.xml_dict['UniRef']['entry']['@id']

    def get_uniref_members(self):
        list_to_process = []
        uniref_90 = set()
        uniref_100 = set()
        uniprot_dict = dict()
        list_to_process.append(self.xml_dict['UniRef']['entry']['representativeMember']['dbReference']['property'])
        if 'member' in self.xml_dict['UniRef']['entry']:
            for member in self.xml_dict['UniRef']['entry']['member']:
                if type(member) == collections.OrderedDict:
                    list_to_process.append(member['dbReference']['property'])

        for member_list in list_to_process:
            uniprot_acc = ""
            org = ""
            protein_name = ""
            for rep in member_list:
                if rep.get('@type', '') == 'UniRef100 ID':
                    uniref_100.add(rep['@value'])
                elif rep.get('@type', '') == 'UniRef90 ID':
                    uniref_90.add(rep['@value'])
                elif rep.get('@type', '') == 'UniProtKB accession':
                    uniprot_acc = rep['@value']
                elif rep.get('@type', '') == 'protein name':
                    protein_name = rep['@value']
                elif rep.get('@type', '') == 'source organism':
                    org = rep['@value']

            if uniprot_acc != '':
                uniprot_dict[uniprot_acc] = [protein_name, org]

        return uniref_90, uniref_100, uniprot_dict

    def load_xml(self, name):
        uniprot_id = strip_uniref_name(name)
        url = f"{self.uniref_url}/?query=member:{uniprot_id}+AND+identity:0.5&format=xml"
        req = requests.get(url)

        if req.status_code in [200]:
            xml = req.text
        else:
            xml = None
            self.log(f"Failed to retrieve xml for {name}. Status code = {req.status_code}")
            return False

        try:
            self.xml_dict = xmltodict.parse(xml)
            self.log(f"XML dict parsed for {name}")
        except Exception as e:
            print(e)
            print(name)

    def check_id_match(self, original_name):
        retrieved_name = self.get_cluster_name()
        if retrieved_name != original_name:
            return False
        return True

    def log(self, msg, level=1):
        if level <= self.log_level:
            print(f"Uniref_Parser: {msg}")

if __name__ == "__main__":
    test_name = "UniRef50_Q7YR71"
    print('start')
    ref_parser = UniRef_Parser()
    ref_parser.load_xml(test_name)
    uni90, uni100, uniprot = ref_parser.get_uniref_members()
    ref_parser.get_cluster_name()
    print(uni100)

    rep_seq = list(uniprot.keys())[0]
    prot_parser = UniProt_Parser()
    prot_parser.load_xml(rep_seq)
    pfams = prot_parser.get_pfams()

    print(pfams)
