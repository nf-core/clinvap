import pandas as pd

def join_columns(row):
    """Join info level with evidence level"""
    row1 = row["evidence_level"]
    row2 = row["info_level"]
    if row1 != "null" and row2 != "null":
        match_level = "-".join([row1, row2])
    elif row1 == "null" and row2 != "null":
        match_level = row2
    elif row1 != "null" and row2 == "null":
        match_level = row1
    else:
        match_level = "null"
    return match_level



def combined_join_columns(row):
    """Function to combine match level and evidence level columns if dataframe has info_pair column."""
    row1 = row["info_level"]
    row2 = row["info_pair_level"]
    row3 = row["evidence_level"]

    if row1 != "null" and row2 != "null":
        level = "/".join([row["info_level"], row["info_pair_level"]])
    elif row1 != "null" and row2 == "null":
        level = row1
    elif row1 == "null" and row2 != "null":
        level = " /" + row2 
    else:
        level = "null"

    if row3 != "null":
        match_level = "-".join([row["evidence_level"], level])
    else:
        match_level = level

    return match_level


def map_level(dataframe): 
    """Function to assign a measure that describes the gene, variant and consequence similarity between 
    the database result and the observed variants."""

    mapper = {"Same gene, same variant, same consequence.": "1", "Same gene, different variant, same consequence.": "2",
              "Same gene, different variant, different consequence.": "3"}

    if "info" in dataframe.columns:
        dataframe["info_level"] = dataframe["info"].map(mapper)
    if "info_pair" in dataframe.columns:
        dataframe["info_pair_level"] = dataframe["info_pair"].map(mapper)
    dataframe.fillna("null", inplace=True)
    return dataframe


def result_compatibility(dataframe):
    """Function to join compability results and evidence level as  a match level. """
    dataframe = map_level(dataframe)
    if {"info_level", "info_pair_level"}.issubset(set(dataframe.columns)):
        dataframe["match_level"] = dataframe.apply(combined_join_columns, axis=1)
    elif {"info_level"}.issubset(set(dataframe.columns)):
        dataframe["match_level"] = dataframe.apply(join_columns, axis=1)
    cols = ["info_level", "info", "info_pair_level", "info_pair"]
    dataframe.drop(columns=cols, inplace=True, errors='ignore')
    return dataframe



def compatibility_sort(dataframe):
    """Function to sort the dataframe according to the match level column which is a combination of 
    evidence level and match level."""
    if not dataframe.empty:
        dataframe = result_compatibility(dataframe)
        rank_map = {"A-1": 1, "A-1/1": 1, "A-1/2": 2, "A-2/1": 2, "A-1/3": 3, "A-3/1": 3, "A-2": 4, "A-2/2": 4, "A-3/2": 5, "A-2/3": 5, "A-3":6 , "A-3/3": 6, "A or B-1": 7, "A or B-1/1": 7, "A or B-1/2": 8, "A or B-2/1": 8, "A or B-1/3": 9, "A or B-3/1": 9, "A or B-2": 10, "A or B-2/2": 10, "A or B-3/2": 11, "A or B-2/3": 11, "A or B-3": 12, "A or B-3/3": 12, "B-1": 13, "B-1/1": 13, "B-1/2": 14, "B-2/1": 14, "B-1/3": 15, "B-3/1": 15, "B-2": 16, "B-2/2": 16, "B-3/2": 17, "B-2/3": 17, "B-3": 18, "B-3/3": 18, "B or C-1": 19, "B or C-1/1": 19, "B or C-1/2 ": 20, "B or C-2/1": 20, "B or C-1/3": 21, "B or C-3/1": 21, "B or C-2": 22, "B or C-2/2": 22, "B or C-3/2": 23, "B or C-2/3": 23, "B or C-3": 24, "B or C-3/3": 24, "C-1": 25, "C-1/1": 25, "C-1/2": 26, "C-2/1": 26, "C-1/3": 27, "C-3/1": 27, "C-2": 28, "C-2/2": 28, "C-3/2": 29, "C-2/3": 29, "C-3": 30, "C-3/3": 30, "D-1": 31, "D-1/1": 31, "D-1/2": 32, "D-2/1": 32, "D-1/3": 33, "D-3/1": 33, "D-2": 34, "D-2/2": 34, "D-3/2 ": 35, "D-2/3": 35, "D-3": 36, "D-3/3": 36, "E-1": 37, "E-1/1": 37, "E-1/2": 38, "E-2/1": 38, "E-1/3": 39, "E-3/1": 39, "E-2": 40, "E-2/2": 40, "E-3/2": 41, "E-2/3": 41, "E-3": 42, "E-3/3": 42, "1": 43, "1/1": 43, "1/2 ": 44, "2/1": 44, "1/3": 45, "3/1": 45, "2": 46, "2/2": 46, "3/2": 47, "2/3": 47, "3": 48, "3/3": 48}
        dataframe["rank"] = dataframe["match_level"].map(rank_map)
        dataframe = dataframe.sort_values(by=["rank"]).drop(columns=["rank", "evidence_level"])
    else:
        dataframe = pd.DataFrame(columns=['gene', 'drug_name', 'variant_drug_association', 'tumor_list',
                                          'hgnc_id', 'variant', 'variant_type', 'match_level', 'reference_id'])
    return dataframe


def reference_sort(dataframe):
    """Function to sort mechanistic drug targets based on number of references."""
    if not dataframe.empty:
        dataframe["rank"] = dataframe["ref_map"].str.count(",")
        dataframe = dataframe.sort_values(by=["rank"], ascending=False).drop(columns=["rank"])
    return dataframe


def driver_type_sort(dataframe):
    """Function to sort driver genes based on driver type."""
    if not dataframe.empty:
        dataframe["driver_type"] = pd.Categorical(dataframe["driver_role"], ["Oncogene", "TSG", "Oncogene/TSG", "Unknown"])
        dataframe["rank"] = dataframe["ref_map"].str.count(",")
        dataframe = dataframe.sort_values(by=["rank", "driver_type"], ascending=[False, True]).drop(columns=["rank"])
    return dataframe


def unnest_ref(dataframe, group_cols):
    """Function to join references of same rows."""
    if not dataframe.empty:
        dataframe = dataframe.groupby(group_cols, sort=False)["reference_id"].apply("|".join).reset_index()
    return dataframe
