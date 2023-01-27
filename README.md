# SIDM
Repo to start playing with SIDM analysis at coffea-casa. This may become the full-blown analysis repo, or I might wind up forking github.com/phylsix/Firefighter and/or github.com/phylsix/FireROOT

## Getting started
- Fork this repository ([here's a nice guide to follow](https://gist.github.com/Chaser324/ce0505fbed06b947d962))
- Log in to coffea.casa as described [here](https://coffea-casa.readthedocs.io/en/latest/cc_user.html#cms-authz-authentication-instance)
- Clone your fork of this repository using the [coffea-casa git interface](https://coffea-casa.readthedocs.io/en/latest/cc_user.html#using-git). Note the following:
  - Use the https link, not the ssh one (e.g. `https://github.com/btcardwell/SIDM.git`, not `git@github.com:btcardwell/SIDM.git`)
  - You will likely need to generate a [github personal access token](https://docs.github.com/en/enterprise-server@3.4/authentication/keeping-your-account-and-data-secure/creating-a-personal-access-token) to log in to your github account on coffea-casa
- Navigate to the newly created SIDM directory using the coffea-casa file browser
- You should be good to go! If you want to test that your environment is set up correctly, try running any of the existing notebooks in SIDM/analysis/studies and comparing your output with the output shown [here](https://github.com/btcardwell/SIDM/tree/main/analysis/studies)

## Code structure
All the interesting code in this repository is in `SIDM/analysis/`, which is organized into the following subdirectories:
- `tools` contains the classes and methods that form the backbone of SIDM. In particular, `ffschema.py` defines how [Firefighter ntuples](https://github.com/phylsix/Firefighter) are converted to a data format we can use, and `sidm_processor.py` defines how events produced by `ffschema.py` are analyzed.
- `definitions` contains files that define the current set of histograms one can make, the cuts one can apply, and some of the objects one can use when making histograms and applying cuts.
- `configs` contains yaml configuration files that define how histograms are grouped together into collections and cuts are grouped together into selections. When running `sidm_processor`, one provides names of selections and names of histogram collections to choose which cuts to apply and which histograms to make.
- `test_notebooks` contains notebooks to test new classes or functionalities as they are added. These notebooks can also serve as a form  of unit test: if you edit some code in a way that you think shouldn't affect the behavior, you can run these notebooks to confirm the output is unchanged.
- `studies` is where the physics happens. The notebooks in this directory are meant to serve effectively as pages in a lab notebook. The intention is to create a new notebook for each unique physics study and to include markdown comments to describe the intentions and observations of the person performing the study. These notebooks can also serve as unit tests in the same way as those in `test_notebooks`.

## Analysis how-tos

### General workflow
I suggest the following workflow for performing a physics study:
1. Create a new branch with a descriptive name (e.g. `leptonIsolation`)
2. Create a new notebook in `studies` with a descriptive name (e.g. `study_lepton_isolation.ipynb`)
3. Following the examples of the existing notebooks in `studies`, use `sidm_processor` to apply cuts and make histograms starting from an Firefighter ntuple of your choosing. Make sure to describe your reasoning and observations in text as you go.
4. If you find you need to define new selections, new cuts, or new histograms, follow the guides below.
5. Commit your changes as you go and submit a Pull Request once you have a reasonable standalone study or have added new selections, cuts, histograms, objects, classes, or features.

### How to define a new histogram and add it to a collection
1. Add an entry to the `hist_defs` dictionary inside [hists.py](https://github.com/btcardwell/SIDM/blob/440069c11e78814da88c86e67fe635d4b655ef6d/analysis/definitions/hists.py). One can potentially do this by mimicking the structure of the existing histograms, but here are some details for those who are interested:
    - Each histogram in `hist_defs` requires a name and Histogram object, which is created using the [Histogram constructor](https://github.com/btcardwell/SIDM/blob/440069c11e78814da88c86e67fe635d4b655ef6d/analysis/tools/histogram.py#L14-L18).
    - The only required argument when creating a Hist is a list of Axis objects, which are created using the [Axis constructor](https://github.com/btcardwell/SIDM/blob/440069c11e78814da88c86e67fe635d4b655ef6d/analysis/tools/histogram.py#L47-L50).
    - When creating an Axis, one must provide a [Hist axis](https://hist.readthedocs.io/en/latest/user-guide/axes.html) and a fill function, which is a lambda expression that defines the quantity used to fill the histogram axis. Note that the fill function will usually take the object from which the quantity is derived as an argument (this is the magic bit that allows us to define how histograms will be filled before running the analysis).
2. After defining your new histogram in `hists.py`, add the name of the histogram to one of the histogram collections in [hist_collections.yaml](https://github.com/btcardwell/SIDM/blob/440069c11e78814da88c86e67fe635d4b655ef6d/analysis/configs/hist_collections.yaml). Note that the new histogram will automatically be included in any collections that include the collection to which it was added (e.g. any histograms added to [electron_base](https://github.com/btcardwell/SIDM/blob/440069c11e78814da88c86e67fe635d4b655ef6d/analysis/configs/hist_collections.yaml#L12-L15) will automatically be included in [base](https://github.com/btcardwell/SIDM/blob/440069c11e78814da88c86e67fe635d4b655ef6d/analysis/configs/hist_collections.yaml#L77-L88))
3. That's it! Running `sidm_processor` while specifying the proper histogram collection will now produce the new histogram.

### How to define a new cut
Coming soon...

### How to define a new selection
Coming soon...

## Miscellaneous how-tos

### How to update requirements.txt
```
cd SIDM/
pip install pipreqs
pipreqs . --force # overwrites current requirements.txt
```
