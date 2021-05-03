import db_process.helpers as helper
import pandas as pd


def assign_driver_role(role_set):
    """Function to collapse driver info from different sources into one, take driver role set, 
    assign driver role to a string"""
    if len(list(role_set)) == 1:
        driver_role = list(role_set)[0]
    elif role_set == {"TSG", "Oncogene"}:
        driver_role = "Oncogene/TSG"
    elif role_set == {"Oncogene/TSG", "TSG", "Oncogene"}:
        driver_role = "Oncogene/TSG"
    elif role_set == {"Oncogene/TSG", "Unknown", "TSG", "Oncogene"}:
        driver_role = "Oncogene/TSG"
    elif role_set == {"Unknown", "TSG", "Oncogene"}:
        driver_role = "Oncogene/TSG"
    elif role_set == {"Unknown", "TSG"}:
        driver_role = "TSG"
    elif role_set == {"Unknown", "Oncogene"}:
        driver_role = "Oncogene"
    elif role_set == {"Oncogene/TSG", 'Unknown'}:
        driver_role = "Oncogene/TSG"
    elif role_set == {"Oncogene/TSG", "TSG"}:
        driver_role = "Oncogene/TSG"
    elif role_set == {"Oncogene/TSG", "Unknown", "TSG"}:
        driver_role = "Oncogene/TSG"
    elif role_set == {"Oncogene/TSG", "Oncogene"}:
        driver_role = "Oncogene/TSG"
    return driver_role


def dict_to_dataframe(driver_dict):
    """Function to convert driver information dictionary into dataframe."""
    df_driver = pd.DataFrame.from_dict(
        driver_dict, orient="index").rename_axis("hgnc_id").reset_index()
    df_driver = df_driver.drop_duplicates().reset_index()
    return df_driver


# Update driver dict
# collapse info from sources into one representation
def get_driver_df(driver_dict):
    """Function to process driver annotation results dictionary. 
    It collapses info from different sources into one representation. 
    Then creates dataframe from updated dictionary."""
    for key, annotation in driver_dict.items():
        role = set()
        var_type = set()
        source = []
        source_pmid = []
        reference = []
        reference_source = set()
        tumor = []
        db_tumor = []
        for driver_info in annotation:
            # fn updating global role_set
            helper.populate_set(driver_info, "driver_role", role)
            # fn updating global variant type
            helper.populate_set(driver_info, "var_type", var_type)
            # fn updating global source_list
            helper.populate_list(driver_info, "source_name", source)
            # fn updating global source_pmid_list
            helper.populate_list(driver_info, "source_pmid", source_pmid)
            # fn updating global reference_id list
            helper.populate_list(driver_info, "reference_id", reference)
            # fn updating global reference_source list
            helper.populate_set(
                driver_info, "reference_source", reference_source)
            # fn updating global tumor_list
            helper.populate_list(driver_info, "tumor_list", tumor)
            # fn updating global db disease keywords
            helper.populate_list(driver_info, "db_tumor_repr", db_tumor)

        driver_role = assign_driver_role(role)
        agg_var_type = "|".join(list(var_type))
        agg_source = "|".join(source)
        agg_pmid = "|".join(source_pmid)
        agg_ref = "|".join(reference)
        agg_ref_source = "|".join(list(reference_source))
        agg_tumor = "|".join([t for t in tumor if t != ""])
        agg_db_tumor = "|".join([t for t in db_tumor if t != ""])

        # replace the value of the dict. collapse list into one dict
        driver_dict[key] = {"driver_role": driver_role, "source_name": agg_source,
                            "source_pmid": agg_pmid, "reference_id": agg_ref, "reference_source": agg_ref_source, "tumor_list": agg_tumor, "db_tumor_repr": agg_db_tumor, "var_type": agg_var_type}

    #  driver information dictionary to dataframe
    driver_dataframe = dict_to_dataframe(driver_dict)
    return driver_dataframe
