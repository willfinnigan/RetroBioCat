from retrobiocat_web.retro.generation import node_analysis
from retrobiocat_web.app.retrobiocat.functions import get_images

def df_to_json_tabulate(df, img_path='static/specificity_images/'):
    img_size = '50px'

    table_columns = [{'title': "Substrate 1", 'field':"Substrate 1"},
                     {'title': "Substrate 2", 'field': "Substrate 2"},
                     {'title': "Product 1", 'field': "Product 1"},

                     {'title': "Substrate 1 Structure", 'field': "Substrate 1 Structure" ,'formatter':"html", 'height':50},
                     {'title': "Substrate 2 Structure", 'field': "Substrate 2 Structure", 'formatter':"html", 'height':50, 'variableHeight':'true'},
                     {'title': "Product 1 Structure", 'field': "Product 1 Structure", 'formatter':"html", 'variableHeight':'true'},

                     {'title': "Similarity", 'field': "Similarity"},

                     {'title': "Active", 'field': "Active", 'formatter':"tickCross"},

                     ]

    table_data = []
    for index, row in df.iterrows():
        substrate_1 = node_analysis.rdkit_smile(str(row['Substrate 1 SMILES']))
        substrate_2 = node_analysis.rdkit_smile(str(row['Substrate 2 SMILES']))
        product_1 = node_analysis.rdkit_smile(str(row['Product 1 SMILES']))

        if substrate_1 == None:
            substrate_1 = ''
        if substrate_2 == None:
            substrate_2 = ''
        if product_1 == None:
            product_1 = ''

        for smiles in [substrate_1, substrate_2, product_1]:
            if smiles != '':
                try:
                    get_images.get_images_of_substrates([smiles], img_dir=img_path)
                except:
                    pass



        substrate_1_img = '<img src="static/specificity_images/' + str(
            get_images.apply_smiles_to_filename_check(substrate_1)) + '.png">'
        substrate_2_img = '<img src="static/specificity_images/' + str(
            get_images.apply_smiles_to_filename_check(substrate_2)) + '.png">'
        product_1_img = '<img src="static/specificity_images/' + str(
            get_images.apply_smiles_to_filename_check(product_1)) + '.png">'

        row_dict = {'id' : index,
                    'Substrate 1' : substrate_1,
                    'Substrate 1 Structure': substrate_1_img,
                    'Substrate 2': substrate_2,
                    'Substrate 2 Structure': substrate_2_img,
                    'Product 1' : product_1,
                    'Product 1 Structure': product_1_img,
         }


        table_data.append(row_dict)

    return table_columns, table_data

