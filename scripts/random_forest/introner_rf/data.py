import pandas as pd


def load_dataset(matrix) -> pd.DataFrame:
    df: pd.DataFrame = pd.read_csv(matrix, sep="\t", header=0)
    return df


def prepare_dataset(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["group_id"] = df["gene"].fillna(df["sequence_id"])
    df["family"] = df["family"].astype("string")
    return df


def load_prepared_dataset(matrix) -> pd.DataFrame:
    return prepare_dataset(load_dataset(matrix))
