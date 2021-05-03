import db_process.helpers as helper
import db_process.query as query
import pandas as pd


def get_subset_dict(dictionary, subkey):
    """Function to seperate inner key and its values - to extracting pharmacogenomics_therapeutics, 
    pharmacogenomics_combined_variants_therapeutics, adverse_effect_therapeutics, 
    adverse_effect_combined_variants_therapeutics."""
    sub = {}
    for key, value in dictionary.items():
        sub[key] = []
        for sub_dict in value:
            if not sub_dict[subkey]:
                continue
            sub[key].append({"variant": sub_dict["variant"], "variant_class": sub_dict["variant_class"], "variant_type": sub_dict["variant_type"], subkey: sub_dict[subkey],
                             "chromosome": sub_dict["chromosome"], "assembly_version": sub_dict["assembly_version"], "alteration_base": sub_dict["alteration_base"],
                             "reference_base": sub_dict["reference_base"], "start": sub_dict["start"], "stop": sub_dict["stop"]})
        if not sub[key]:
            del sub[key]
    return sub


def apply_add_tumor_string(dictionary, keyword):
    """Function to apply add_tumor_string to pharmaco-dynamics data
    Warning: it changes the input dictionary"""
    for value in dictionary.values():
        for variant in value:
            helper.add_tumor_string(keyword, variant)
    return dictionary


def filter_gene_pair(gene_list_snv, gene_list_cnv, dictionary, keyword):
    """Function to filter combined_variants_therapeutics. Both genes should be observed in the input VCF.
    This is for combined variant effect datasets, i.e. pharmacogenomics_combined_variants_therapeutics and 
    adverse_effect_combined_variants_therapeutics. Warning: this function is not for direct targets"""
    filtered_dict = {}
    for key, value in dictionary.items():
        filtered_dict[key] = []
        for element in value:
            filtered_items = []
            # filtered_items = [item for item in element[keyword] if str(item["paired_gene_hgnc"]) in gene_list]
            if ";MUT" in element["variant_class"]:
                filtered_items = [item for item in element[keyword] if str(
                    item["paired_gene_hgnc"]) in gene_list_snv]
            elif ";CNA" in element["variant_class"]:
                filtered_items = [item for item in element[keyword] if (
                    str(item["paired_gene_hgnc"]), item["paired_variant"], "cnv") in gene_list_cnv]
            elif ";FUS" in element["variant_class"]:
                # not implemented
                pass
            if filtered_items:
                filtered_dict[key].append({"variant": element["variant"],
                                           "variant_class": element["variant_class"], "variant_type": element["variant_type"],
                                           "chromosome": element["chromosome"], "assembly_version": element["assembly_version"],
                                           "alteration_base": element["alteration_base"], "reference_base": element["reference_base"],
                                           "start": element["start"], "stop": element["stop"], keyword: filtered_items})
        if not filtered_dict[key]:
            del filtered_dict[key]
    return filtered_dict


def dict_to_dataframe(pharm_dictionary, keyword):
    """Function for converting pharm dictipnaries into pandas dataframe.
    Warning: it changes the input dictionary"""
    pharm_list = []
    for key, value in pharm_dictionary.items():
        for element in value:
            for item in element[keyword]:
                item["hgnc_id"] = key
                item["variant"] = element["variant"]
                item["variant_class"] = element["variant_class"]
                item["variant_type"] = element["variant_type"]
                item["chromosome"] = element["chromosome"]
                item["assembly_version"] = element["assembly_version"]
                item["alteration_base"] = element["alteration_base"]
                item["reference_base"] = element["reference_base"]
                item["start"] = element["start"]
                item["stop"] = element["stop"]
                pharm_list.append(item)
    df_pharm = pd.DataFrame(pharm_list)
    df_pharm = df_pharm.drop_duplicates().reset_index(drop=True)
    return df_pharm


def direct_adverse_combined(variant_class="MUT"):
    """Database does not include any such cases for MUT and MUT;MUT category"""
    print("Not implemented. No such entries in the database.")
