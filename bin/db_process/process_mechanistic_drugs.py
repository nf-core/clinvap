import db_process.helpers as helper
import pandas as pd


def reshape(dictionary, temp):
    """Function to first remove fields that are not reported in the final product. Then, to reshape mechanistic 
    drug targets dictionary by creating nested dictionaries within genes (main keys), using drugbank ids as 
    keys for the nested part. This is created to merge information of same drugs coming from different 
    sources for a given gene."""
    for key, value in dictionary.items():
        temp[key] = {}
        for drug in value:
            helper.remove_key_value(drug, "drug_type", "interaction_type",
                                    "is_cancer_drug", "target_action", "target_known_action")  # Remove fields that will not be included in the final product.
            if drug["drugbank_id"] not in temp[key]:
                temp[key][drug["drugbank_id"]] = []
            temp_drug = temp[key][drug["drugbank_id"]]
            temp_drug.append(drug)
    dictionary = temp
    return dictionary


def assign_approval_status(role_set):
    """Function to collapse approval status coming from diffeent sources and solving the conradicting info 
    if there is any. """
    if role_set == {"Approved"}:
        approval_status = "approved"
    elif role_set == {"approved|vet_approved"}:
        approval_status = "approved"
    elif role_set == {"investigational|vet_approved"}:
        approval_status = "investigational"
    elif len(list(role_set)) == 1:
        approval_status = list(role_set)[0]
    elif role_set == {"approved", "Approved"}:
        approval_status = "approved"
    elif role_set == {"approved|vet_approved", "Approved"}:
        approval_status = "approved"
    elif role_set == {"approved|investigational", "Approved"}:
        approval_status = "approved|investigational"
    elif role_set == {"approved|investigational|withdrawn", "Approved"}:
        approval_status = "approved|investigational|withdrawn"
    elif role_set == {"approved|investigational|vet_approved", "Approved"}:
        approval_status = "approved|investigational"
    elif role_set == {"approved|experimental|investigational", "Approved"}:
        approval_status = "approved|experimental|investigational"
    elif role_set == {"investigational", "Phase 3"}:
        approval_status = "investigational"
    elif role_set == {"approved|investigational|withdrawn", "Withdrawn from market"}:
        approval_status = "Withdrawn from market"
    return approval_status


EXCLUDED_VALUES = ["withdrawn",
                   "Withdrawn from market", "Discontinued in Phase 3", "Discontinued in Phase 2", "Terminated", "experimental"]

SEPERATED_VALUES = ["investigational",
                    "Phase 3", "Phase 1", "Phase 2", "Phase 1/2", "Phase 2/3", "experimental|investigational"]


def filter_approval_info(dictionary):
    """Function to filter drugs which are not approved"""
    for key in list(dictionary):
        modified = [item for item in dictionary[key]
                    if item["approval_status"] not in EXCLUDED_VALUES]
        dictionary[key] = modified
        if not dictionary[key]:
            helper.remove_key_value(dictionary, key)
    return dictionary


def divide_approved_investigational(dictionary):
    """Function to divide data into two for approved drugs and investigational drugs."""
    investigational_dict = {}
    approved_dict = {}
    for key, value in dictionary.items():
        investigational = [
            item for item in value if item["approval_status"] in SEPERATED_VALUES]
        approved = [item for item in value if item["approval_status"]
                    not in SEPERATED_VALUES]
        if investigational:
            investigational_dict[key] = investigational
        if approved:
            approved_dict[key] = approved
    return approved_dict, investigational_dict


def dict_to_dataframe(mechanistic_dict):
    """Function to create dataframe from mechanistic drug targets. Warning it changes the input dictionary."""
    item_list = []
    for key, value in mechanistic_dict.items():
        for element in value:
            element["hgnc_id"] = key
            item_list.append(element)
    df_mechanistic = pd.DataFrame(item_list)
    df_mechanistic = df_mechanistic.drop_duplicates().reset_index(drop=True)

    return df_mechanistic
