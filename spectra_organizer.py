# %% Imports
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os


# %% Functions
def get_filenames(par_pattern):
    """
    get_filenames(str) -> list[str
    Parameters:
    -----------
    An extension so that the program can look for files with such extension
    and make a list of those file
    Returns:
    --------
    A list with the names of the files containing the data that is going to be
    in the plo
    Example:
    --------
    get_filenames('.txt') -> ['file1.txt', 'file2.txt']
    """
    return list(filter(lambda x: par_pattern in x or "white" in x, os.listdir()))


def transform_name(par_str):
    """
    transform_name(str) -> st
    rename files so that white standard is distinguishable from
    sample spectra
    """
    if par_str.startswith("white"):
        return par_str
    elif par_str.startswith("m"):
        return par_str.replace("m", "-")
    elif par_str[0].isnumeric():
        return par_str.split("_")[0]
    else:
        names_ls = par_str.split("_")
        pn = names_ls[0][-1]
        tailnm = "".join(names_ls[3].split("-"))
        new_name = names_ls[1] + pn + tailnm
        new_name = new_name.replace("m", "-")
        return new_name


def read_data(par_list):
    """
    read_data(list[str]) -> pandas.DataFrame

    Parameters:
    -----------
    A list of strings corresponding to filenames to process

    Returns:
    --------
    A single pandas.DataFrame with the aggregation of the data from the file(s)
    in the list
    """
    holder_list = []
    white_df = None
    for i in par_list:
        curr_name = i.strip(".txt")
        curr_name = transform_name(curr_name)
        if i.startswith("white"):
            white_df = pd.read_csv(
                i,
                skiprows=50,
                sep="\t",
                names=["wavelength", "intensity"],
                engine="python",
                skipfooter=1,
            )
            white_df = white_df.assign(name=np.repeat(curr_name, white_df.shape[0]))
            # Apply corrections dividing by 1e07
            white_df["wavelength"] = white_df.wavelength.values / 1e07
        else:
            curr_data = pd.read_csv(
                i,
                skiprows=50,
                sep="\t",
                names=["wavelength", "intensity"],
                engine="python",
                skipfooter=1,
            )
            name_col = np.repeat(curr_name, curr_data.shape[0])
            curr_data = curr_data.assign(name=pd.to_numeric(name_col))

            # Apply corrections dividing by 1e07
            curr_data["wavelength"] = curr_data.wavelength.values / 1e07
            holder_list.append(curr_data)
    transformed_list = []
    for df in holder_list:
        new_intensity = (df.intensity / white_df.intensity) * 100
        transformed_list.append(df.assign(n_intensity=new_intensity))
    return pd.concat(transformed_list)


def plot_2d(par_df, savefig=False, showfig=False, fmt="pdf", legend=True):
    """
    plot_2d(pandas.Dataframe, bool, bool) -> None

    Parameters:
    -----------
    A pandas.DataFrame with 4 columns:
    <float> <float> <int> <float>

    Returns:
    --------
    None
    """
    unique_names = list(np.unique(par_df.name.values))
    n_colors = len(unique_names)
    list_colours = np.apply_along_axis(
        mpl.colors.rgb2hex, 1, mpl.cm.get_cmap(name="viridis", lut=n_colors).colors
    )
    title = os.getcwd().split("\\")[-1]
    for i, j in zip(unique_names, list_colours):
        my_subset = par_df.loc[par_df["name"] == i]
        plt.scatter(my_subset.wavelength, my_subset.n_intensity, label=i, s=8, color=j)
    plt.title(title)
    plt.xlabel("Wavelength (nm)")
    plt.ylabel("Reflectance (%)")
    if legend:
        plt.legend(markerscale=2)
    plt.tight_layout()
    if savefig:
        plt.savefig("{}.{}".format(title, fmt), format=fmt)
        print("{}.{} saved".format(title, fmt))
    if not showfig:
        plt.close()


def save_pavo(par_df):
    """
    Saves the data in a format suitable for pavo analysis
    """
    #
    angles = list(np.unique(par_df["name"]))
    for angle in angles:
        sub_df = par_df.loc[par_df["name"] == angle][["wavelength", "n_intensity"]]
        sub_df.to_csv(
            "{}_{}_PAVO.txt".format(str(angle), os.getcwd().split("\\")[-1]),
            sep="\t",
            index=False,
            float_format="%.3f",
        )


# %% Load data
# Alternative method commented out
# my_list = [str(i).strip() for i in open("list_of_files.txt",
#           "r").readlines()]
ind_name = os.getcwd().split("\\")[-1]
my_list = get_filenames(".txt")
all_data = read_data(my_list)
all_data.to_csv("{}.tsv".format(ind_name), sep="\t", index=False, float_format="%.3f")

all_data = all_data.loc[all_data["wavelength"] > 260]
all_data = all_data.loc[all_data["wavelength"] < 700]
# Plot
plot_2d(all_data, savefig=True, showfig=True, legend=True, fmt="png")

# Save to pavo friendly format
all_data = read_data(my_list)
save_pavo(all_data)
