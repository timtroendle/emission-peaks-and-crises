from importlib.resources import path
from PIL import Image


def preprocess_flag(path_to_flag: str, new_height: int, path_to_output: str):
    im = Image.open(path_to_flag)
    aspect_ratio = im.height / im.width
    new_width = int(new_height / aspect_ratio)
    im = im.resize((new_width, new_height), Image.ANTIALIAS)
    im.save(path_to_output)


if __name__ == "__main__":
    preprocess_flag(
        path_to_flag=snakemake.input.flag,
        new_height=int(snakemake.params.height),
        path_to_output=snakemake.output[0]
    )
