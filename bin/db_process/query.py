import json
import pandas
import sys
import copy


def read_knowledgebase(filename):
    """Helper fn, Read Knowledgebase"""
    with open(filename) as db:
        cancer_db = json.load(db)
    return cancer_db


def get_gene_list(dataframe):
    """Helper fn, function to get gene list from the dataframe. Function is extended to return 
    the dictionary of genes and their variant type. """
    if "type" in dataframe.columns:
        new_df = dataframe[["HGNC_ID", "type", "var_type"]]
    else:
        new_df = dataframe[["HGNC_ID", "var_type"]]
    new_df = new_df.drop_duplicates()
    gene_list = new_df["HGNC_ID"].tolist()
    gene_dict = new_df.to_dict('records')
    gene_tuple = [tuple(r) for r in new_df.to_numpy()]
    return gene_list, gene_dict, gene_tuple


def get_inner_dict(key1, key2, dictionary):
    """Helper fn, function to get second level inner dictionary from a nested one"""
    if key1 in dictionary:
        inner_dict = dictionary[key1].get(key2)
        inner_d = copy.deepcopy(inner_dict)
        if inner_dict:
            return inner_d
        else:
            return None


def get_variant_list(dataframe):
    """Get variant coordinates into a list. 
    Different fileds are selected for different variant types, i.e for CNA type is selected, ref and alt are not selected.
    Returns list of dictionaries of search subset"""
    if "type" in dataframe.columns:  # implies that it is etiher CNV or SV
        fields = ["HGNC_ID", "chr", "start", "end", "type"]
    else:  # implies that it is either SNP or INDEL
        fields = ["HGNC_ID", "chr", "start", "end", "ref", "alt"]
    rows = dataframe[fields].values.tolist()
    return rows


def get_driver_annotation(dataframe, knowledgebase):
    """Get driver gene annotation, query with hgnc is and variant class. Returns driver genes"""

    def tag_dict(list_dict, tag):
        for dictionary in list_dict:
            dictionary["var_type"] = tag

    driver_mutation_annotation = {}
    driver_indel_annotation = {}
    driver_copynumber_annotation = {}
    gene_dict = get_gene_list(dataframe)[1]
    for gene in gene_dict:
        driver = get_inner_dict(
            gene["HGNC_ID"], "driver_annotation", knowledgebase)
        if not isinstance(driver, dict):
            continue
        mutation = driver.get("mutation")
        indel = driver.get("indel")
        copynumber = driver.get("copy_number")
        if mutation:
            if gene["var_type"] == "snp":
                driver_mutation_annotation[gene["HGNC_ID"]] = mutation
                tag_dict(mutation, "snp")
        if indel:
            if gene["var_type"] == "indel":
                driver_indel_annotation[gene["HGNC_ID"]] = indel
                tag_dict(indel, "indel")
        if copynumber:
            if gene["var_type"] == "cnv":
                type_specific = [cnv_ann for cnv_ann in copynumber if (cnv_ann["alteration_type"]
                                                                       == gene["type"] or cnv_ann["alteration_type"] == "null")]
                if type_specific:
                    driver_copynumber_annotation[gene["HGNC_ID"]
                                                 ] = type_specific
                    tag_dict(type_specific, "cnv")
                else:
                    continue
        # TODO fusion
    return driver_mutation_annotation, driver_indel_annotation, driver_copynumber_annotation


def get_mechanistic_drugs(dataframe, knowledgebase):
    """query with hgnc id, select is_cancer_drug=true cases. Returns all the drugs mechanistically 
    targeting the affected gene."""
    gene_list = get_gene_list(dataframe)[0]
    mechanistic_drug_targets = {}
    for gene in gene_list:
        mechanistic_drug = get_inner_dict(
            gene, "mechanistic_therapeutics_annotation", knowledgebase)
        if not mechanistic_drug:
            continue
        mechanistic_cancer_drug = [  # first filtering condiitions
            drug for drug in mechanistic_drug if (drug["is_cancer_drug"] == True and drug["interaction_type"] == "target")]
        if not mechanistic_cancer_drug:
            continue
        mechanistic_drug_targets[gene] = mechanistic_cancer_drug
    return mechanistic_drug_targets


def ann_pharm_variant(dataframe, genome, knowledgebase, variant_type):
    """query with chr, start, stop, ref, alt, genome assembly version. Returns all the drugs 
    targeting the observed variant. """
    all_direct_targets = {}
    rows = get_variant_list(dataframe)
    for r in rows:
        direct_target_list = []
        gene = r[0]
        if variant_type == "MUT":
            search_subset = {"chromosome": r[1], "start": r[2], "stop": r[3],
                             "reference_base": r[4], "alteration_base": r[5], "assembly_version": genome}
        elif variant_type == "CNA":
            search_subset = {"chromosome": r[1], "start": r[2], "stop": r[3],
                             "global_type": r[4], "assembly_version": genome}
        superset = get_inner_dict(gene, "variant_annotation", knowledgebase)
        if not superset:
            continue
        for variant in superset:
            v_index = superset.index(variant)
            search_set = superset[v_index]
            coverage = all(item in search_set.items()
                           for item in search_subset.items())
            if coverage:
                direct_target_list.append(search_set)
        if not direct_target_list:
            continue
        if gene in all_direct_targets:
            all_direct_targets[gene].extend(direct_target_list)
        else:
            all_direct_targets[gene] = direct_target_list
    return all_direct_targets

# TODO generalize it for SVs


def ann_pharm_gene(dataframe, knowledgebase, genome):
    """query with hgnc id, returns all the drugs targeting the affected gene."""
    drug_targets_dict = {}
    gene_list = get_gene_list(dataframe)[1]
    for gene in gene_list:
        drug_targets = get_inner_dict(
            gene["HGNC_ID"], "variant_annotation", knowledgebase)
        if drug_targets and gene["var_type"] != "cnv":
            drug_targets_dict[gene["HGNC_ID"]] = drug_targets
        elif drug_targets and gene["var_type"] == "cnv":
            filtered = [annotation for annotation in drug_targets if annotation.get(
                "global_type") == gene["type"]]
            if filtered:
                drug_targets_dict[gene["HGNC_ID"]] = filtered
            else:
                continue
    assembly_filtered = filter_genome_assembly(drug_targets_dict, genome)
    return assembly_filtered


def filter_main_variant_class(dictionary, variant_class):
    """Function to filter pharmaco-dyamics results based on variant class. Variant class should be either MUT, 
    CNA or FUS."""
    filtered = {}
    for key, value in dictionary.items():
        filtered[key] = [
            variant for variant in value if variant["variant_class"] == "null" or variant["variant_class"] == variant_class or variant_class+";" in variant["variant_class"]]
    return filtered


def filter_genome_assembly(dictionary, genome):
    """Function to filter pharmacogenomics dictionary(indirect one) based on the genome assembly"""
    filtered = {}
    for key, value in dictionary.items():
        filtered[key] = [
            variant for variant in value if variant["assembly_version"] == genome]
    return filtered
