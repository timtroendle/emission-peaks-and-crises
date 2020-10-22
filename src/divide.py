import pandas as pd


def divide(path_to_dataset1, path_to_dataset2, path_to_output):
    df1 = pd.read_csv(path_to_dataset1, index_col=0)
    df2 = pd.read_csv(path_to_dataset2, index_col=0)
    result = df1 / df2
    result.to_csv(path_to_output, index=True, header=True)


if __name__ == "__main__":
    divide(
        path_to_dataset1=snakemake.input.df1,
        path_to_dataset2=snakemake.input.df2,
        path_to_output=snakemake.output[0]
    )
