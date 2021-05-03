import pandas as pd
import requests
import json
import time


KEY_LIST = ["reference_id", "reference_source"]


def get_references(dataframe, key_list=KEY_LIST):
    """Create reference content: Function to return dataframe of reference information."""
    df_ref = dataframe[key_list].drop_duplicates()
    return df_ref


def reshape(dataframe):
    """Create reference content: Function to reshape reference dataframe to seperate aggregated values into new columns."""
    df_lines = []
    for _, row in dataframe.iterrows():
        ref_sp = row["reference_id"].split('|')
        source_sp = row["reference_source"].split('|')
        for i in range(len(ref_sp)):
            if len(source_sp) == 1:
                line = [ref_sp[i], source_sp[0]]
            else:
                line = [ref_sp[i], source_sp[i]]
            df_lines.append(line)
    df = pd.DataFrame(df_lines, columns=[
                      "reference_id", "reference_source"]).drop_duplicates().reset_index(drop=True)
    return df


def get_pubmed_ids(dataframe):
    """Create reference content: Function to seperate pubmed ids to form eutils query string"""
    pubmed = dataframe[(dataframe["reference_source"] == "PubMed")
                       | (dataframe["reference_source"] == "NCCN")
                       | (dataframe["reference_source"] == "FDA")]
    return pubmed


def id_string(dataframe):
    """Create reference content: Function to create id string"""
    shaped = reshape(dataframe)
    pubmed = get_pubmed_ids(shaped)
    id_string = pubmed["reference_id"].str.cat(sep=',')
    return id_string, shaped


def get_ref_details(URL):
    """Create reference content: Function to get reference details from etuils API into a dataframe"""
    for i in range(3):
        try:
            response = requests.get(URL)
            value = response.json()
            references = value["result"]
            del references["uids"]
            df = pd.DataFrame.from_dict(references, orient='index', columns=[
                "sortfirstauthor", "title", "fulljournalname", "volume", "issue", "pages", "pubdate"])
            df["pubdate"] = df["pubdate"].str.replace("((?<=\d{4}).*)", "")
            df["sortfirstauthor"] = df["sortfirstauthor"].str.replace(
                "((?<=) .*)", " et al.")
            df["title"] = df["title"].str.replace("(\<.*?\>)", "")
            df["combined"] = df.apply(
                lambda row: ", ".join(row.values.astype(str)), axis=1)
            df = df.drop(columns=["sortfirstauthor", "title",
                                "fulljournalname", "volume", "issue", "pages", "pubdate"])
            return df
        except requests.exceptions.RequestException as e:
            if i != 2:
                time.sleep(10)
                continue
            else:
                raise SystemExit(e)



def replace_ref_nan(row):
    """Function to replace nan info values of some references coming from references without PubMed id
    with their reference id."""
    if isinstance(row["combined"], float):
        return row["reference_id"]
    else:
        return row["combined"]


def reference_table(reference_map_complete):
    """Function to create reference appendix table with full citation info and map values of references."""
    reference_map_complete.index += 1
    reference_map_complete = reference_map_complete.reset_index(level=0)
    reference_map_complete["combined"] = reference_map_complete.apply(
        replace_ref_nan, axis=1)
    table_content = reference_map_complete.drop(
        columns=["reference_source"])
    map_dict = reference_map_complete[["index", "reference_id"]].set_index(
        "reference_id").to_dict("dict")
    return table_content, map_dict


def get_drug_class(dataframe):
    """Function to replace drug name with drug class, if drug name is null."""
    df = dataframe.copy()
    if not df.empty:
        cond = (df["drug_name"] == "null") & (df["drug_class"] != "null")
        df.loc[cond, "drug_name"] = df.loc[cond, "drug_class"]
    return df


def get_drug_relation(dataframe):
    """Function to extend drug name drug drug relationship 
    is not null for variants with more than on drug."""
    df = dataframe.copy()
    if not df.empty:
        cond = df["drug_drug_relationship"] != "null"
        df.loc[cond, "drug_name"] = df.loc[cond, "drug_name"] + \
            " (" + df.loc[cond, "drug_drug_relationship"] + ")"
    return df


def assign_match_level_direct(dataframe, col_name):
    processed_info = dataframe.assign(
        col_name="Same gene, same variant, same consequence.").rename(columns={"col_name": col_name})
    return processed_info


def cnv_match_level_indirect(dataframe, col_name):
    processed = dataframe.assign(col_name="Same gene, different variant, same consequence.").rename(
        columns={"col_name": col_name})
    return processed


def pair_match_level(dataframe, mvld):
    cond_snv = dataframe["variant_class"].str.contains(";MUT")
    cond_cnv = dataframe["variant_class"].str.contains(";CNA")
    snv_processed = snv_match_level(
        dataframe[cond_snv], mvld, "info_pair", pair=1)
    cnv_processed = cnv_match_level_indirect(dataframe[cond_cnv], "info_pair")
    processed = pd.concat([snv_processed, cnv_processed]).drop(
        columns=["variant"]).rename(columns={"combined_variant": "variant"})
    return processed


def get_content_pharm_notcomb_dir(dataframe):
    processed = get_drug_class(dataframe)
    processed_1 = get_drug_relation(processed)
    processed_info = assign_match_level_direct(processed_1, "info")
    return processed_info


def cnv_get_content_pharm_notcomb_notdir(dataframe):
    processed = get_drug_class(dataframe)
    processed_1 = get_drug_relation(processed)
    processed_info = cnv_match_level_indirect(processed_1, "info")
    return processed_info


def cnv_get_content_pharm_comb_notdir(dataframe, mvld):
    processed = get_drug_class(dataframe)
    processed_1 = get_drug_relation(processed)
    processed_info = cnv_match_level_indirect(
        processed_1, "info")  # for main gene
    processed_info_comb = pair_match_level(processed_info, mvld)
    return processed_info_comb


def snv_get_content_pharm_comb_dir(dataframe, mvld):
    processed = get_drug_class(dataframe)
    processed_2 = get_drug_relation(processed)
    processed_info = assign_match_level_direct(
        processed_2, "info")  # for main gene
    processed_info_pair = pair_match_level(
        processed_info, mvld)  # for pair gene
    return processed_info_pair


def snv_get_content_pharm_comb_notdir(dataframe, mvld):
    processed = get_drug_class(dataframe)
    processed_1 = get_drug_relation(processed)
    processed_info = snv_match_level(
        processed_1, mvld, "info")  # for main gene
    processed_info_pair = pair_match_level(
        processed_info, mvld)   # for pair gene
    return processed_info_pair


def snv_get_content_pharm_notcomb_notdir(dataframe, mvld):
    processed = get_drug_class(dataframe)
    processed_1 = get_drug_relation(processed)
    processed_info = snv_match_level(processed_1, mvld, "info")
    return processed_info


def snv_get_content_adv_comb_dir(dataframe, mvld):
    processed = get_drug_class(dataframe)
    processed_info = assign_match_level_direct(
        processed, "info")  # for main gene
    processed_info_pair = pair_match_level(
        processed_info, mvld)    # for pair gene
    return processed_info_pair


def snv_get_content_adv_comb_notdir(dataframe, mvld):
    processed = get_drug_class(dataframe)
    processed_info = snv_match_level(processed, mvld, "info")  # for main gene
    processed_info_pair = pair_match_level(processed_info, mvld)
    return processed_info_pair


def snv_get_content_adv_notcomb_dir(dataframe):
    processed = get_drug_class(dataframe)
    processed_info = assign_match_level_direct(processed, "info")
    return processed_info


def snv_get_content_adv_notcomb_notdir(dataframe, mvld):
    processed = get_drug_class(dataframe)
    processed_info = snv_match_level(processed, mvld, "info")
    return processed_info


def snv_get_content(dataframe, mvld, cols, direct=None, combined=None, adverse=None, pharm=None, indirect=None, notcombined=None):
    """Function to create pharmacogenomics and adverse effect dataframes with required 
    table content for SNVs."""
    if direct and combined and pharm:
        processed = snv_get_content_pharm_comb_dir(dataframe, mvld)
    if indirect and combined and pharm:
        processed = snv_get_content_pharm_comb_notdir(dataframe, mvld)
    if direct and notcombined and pharm:
        processed = get_content_pharm_notcomb_dir(dataframe)
    if indirect and notcombined and pharm:
        processed = snv_get_content_pharm_notcomb_notdir(dataframe, mvld)
    if direct and combined and adverse:
        processed = snv_get_content_adv_comb_dir(dataframe, mvld)
    if indirect and combined and adverse:
        processed = snv_get_content_adv_comb_notdir(dataframe, mvld)
    if direct and notcombined and adverse:
        processed = snv_get_content_adv_notcomb_dir(dataframe)
    if indirect and notcombined and adverse:
        processed = snv_get_content_adv_notcomb_notdir(dataframe, mvld)

    content = pd.merge(processed, mvld, how="left",
                       left_on="hgnc_id", right_on="HGNC_ID")[cols]

    return content


def cnv_get_content(dataframe, cnv_mvld, cols, mvld=None, direct=None, combined=None, pharm=None, indirect=None, notcombined=None):
    """Function to create pharmacogenomics and adverse effect dataframes with required 
    table content for CNVs."""
    if pharm and direct and notcombined:
        processed = get_content_pharm_notcomb_dir(dataframe)
    if pharm and indirect and notcombined:
        processed = cnv_get_content_pharm_notcomb_notdir(dataframe)
    if pharm and indirect and combined:
        processed = cnv_get_content_pharm_comb_notdir(dataframe, mvld)

    content = pd.merge(processed, cnv_mvld, how="left",
                       left_on="hgnc_id", right_on="HGNC_ID")[cols]

    return content


def get_consequence_list(mvld):
    """Function to expand the consequences seperated by & and to create list of gene, 
    consequence for match_level() to use."""

    def consequence_reshape(dataframe, col1, col2, delimiter):
        """Function to reshape reference dataframe to seperate aggregated values into new columns."""
        df_lines = []
        for _, row in dataframe.iterrows():
            col_split = row[col2].split(delimiter)
            for i in range(len(col_split)):
                line = [row[col1], col_split[i]]
                df_lines.append(line)
        df = pd.DataFrame(df_lines, columns=[
            col1, col2]).drop_duplicates().reset_index(drop=True)
        return df

    df = mvld.copy()
    df = consequence_reshape(
        df[["HGNC_ID", "Consequence"]], "HGNC_ID", "Consequence", "&")
    consequence_list = df.values.tolist()
    return consequence_list


def snv_match_level(dataframe, mvld, col_name, pair=None):
    """Function the assign similarity levels of the genes, variants, consequences 
    between the database results and the input"""
    df = dataframe.copy()
    df[col_name] = ""
    variant = mvld[["HGNC_ID", "one_letter_repr"]].values.tolist()
    consequence = get_consequence_list(mvld)
    for index, row in df.iterrows():
        if pair:
            comp = [str(row["paired_gene_hgnc"]), row["paired_variant"],
                    row["variant_type"].split(',')[1]]
        elif "," in row["variant_type"]:
            comp = [row["hgnc_id"], row["variant"],
                    row["variant_type"].split(',')[0]]
        else:
            comp = [row["hgnc_id"], row["variant"], row["variant_type"]]

        if [comp[0], comp[1]] in variant:
            df.loc[index, col_name] = "Same gene, same variant, same consequence."
        elif [comp[0], comp[2]] in consequence:
            df.loc[index, col_name] = "Same gene, different variant, same consequence."
        else:
            df.loc[index, col_name] = "Same gene, different variant, different consequence."
    return df


def get_mechanistic_content(dataframe, mvld, cols):
    """Function to get mechanistic drugs content by getting gene symbols and selecting columns to be reported."""
    df = dataframe.copy()
    processed_1 = pd.merge(df, mvld, how="left",
                           right_on="HGNC_ID", left_on="hgnc_id")
    processed_2 = processed_1[cols].drop_duplicates()
    return processed_2


def get_mutation(row):
    """Helper function for appendix_variants()"""
    if row["HGVSp"] != "null":
        return row["HGVSp"]
    elif row["HGVSc"] != "null":
        return row["HGVSc"]
    else:
        return "null"


def appendix_variants(mvld, var_type):
    """Function the create the appendix table for every variant observed in the input file."""
    df = mvld.copy()
    if var_type == "SNV":
        df = df[["HGVSp", "HGVSc", "gene", "dbSNP",
                 "COSMIC", "Consequence", "vaf"]].fillna("null")
        df["HGVSp"] = df["HGVSp"].str.replace("EN.*:", "")
        df["HGVSc"] = df["HGVSp"].str.replace("EN.*:", "")
        df["mutation"] = df.apply(get_mutation, axis=1)
        df = df.drop(columns=["HGVSp", "HGVSc"]).drop_duplicates()
    elif var_type == "CNV":
        df = df[["type", "copy_number", "gene"]].fillna("null")
    return df


def reference_map_ids(dataframe, map_dict):
    """Function to match reference ids with the references in the appendix table through map ids. """
    dataframe["ref_map"] = ""
    for index, row in dataframe.iterrows():
        refs = row["reference_id"].split("|")
        refs = list(set(refs))
        maplist = [map_dict["index"][r]
                   for r in refs if r in map_dict["index"]]
        maplist.sort()
        dataframe.loc[index, "ref_map"] = ",".join([str(m) for m in maplist])
    return dataframe


def nonsynonymous_count(mvld):
    """Function to retrieve the number os non-synonymous mutations present in the vcf file."""
    missense = mvld["Consequence"].str.contains(
        ":?missense_variant|start_lost?|stop_lost?|stop_gained?")
    df_missense = mvld[missense][["chr", "start", "end",
                                  "ref", "alt"]].drop_duplicates().reset_index(drop=True)
    return len(df_missense.index)


def oncogenes_count(driver_dataframe):
    """Function to count number of oncogenes from final driver gene dataframe."""
    oncogene = driver_dataframe["driver_role"] == "Oncogene"
    df_oncogene = driver_dataframe[oncogene][[
        "gene", "driver_role"]].drop_duplicates().reset_index(drop=True)
    return len(df_oncogene.index)


def tsg_count(driver_dataframe):
    """Function to count number of tumor suppressor genes from final driver gene dataframe."""
    tsg = driver_dataframe["driver_role"] == "TSG"
    df_tsg = driver_dataframe[tsg][["gene", "driver_role"]
                                   ].drop_duplicates().reset_index(drop=True)
    return len(df_tsg.index)


def onco_tsg_count(driver_dataframe):
    """Function to count number of both onco and tumor suppressor genes from final driver gene dataframe."""
    onco_tsg = driver_dataframe["driver_role"] == "Oncogene/TSG"
    df_onco_tsg = driver_dataframe[onco_tsg][[
        "gene", "driver_role"]].drop_duplicates().reset_index(drop=True)
    return len(df_onco_tsg.index)


def unknown_count(driver_dataframe):
    """Function to count number of both onco and tumor suppressor genes from final driver gene dataframe."""
    unknown = driver_dataframe["driver_role"] == "Unknown"
    df_unknown = driver_dataframe[unknown][[
        "gene", "driver_role"]].drop_duplicates().reset_index(drop=True)
    return len(df_unknown.index)
