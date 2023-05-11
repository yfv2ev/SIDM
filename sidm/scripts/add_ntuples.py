# Tool to add ntuple locations to sidm/configs/ntuple_locations.yaml. Developed and tested on FNAL
# LPC. Will be updated after issue https://github.com/CoffeaTeam/coffea-casa/issues/374 is
# resolved. Note that cmsenv or equivalent is needed to import XRootD.

# Usage: python add_ntuples.py -o OUTPUT_CONFIG -n NTUPLE_NAME -c NTUPLE_COMMENT -d NTUPLE_ROOT_DIR

from __future__ import print_function
import argparse
import yaml
from XRootD import client


def parse_name(name):
    """Parse sample directory name to produce simplified name
    
    Assumes structure like "SIDM_XXTo2ATo2Mu2E_mXX-100_mA-1p2_ctau-0p096_TuneCP..."
    """

    process_names = {
        "SIDM_XXTo2ATo2Mu2E_mXX": "2Mu2E_",
        # fixme: add 4mu and backgrounds as necessary
    }
    chunks = name.split("-")
    simplified_name = process_names[chunks[0]] # process name
    simplified_name += chunks[1].replace("_mA", "GeV_") # bound state mass
    simplified_name += chunks[2].replace("_ctau", "GeV_") # dark photon mass
    simplified_name += chunks[3].split("_TuneCP")[0] + "mm" # dark photon ctau

    return(simplified_name)

def descend(ntuple_path, sample_path):
    path = ntuple_path + "/" + sample_path
    dir_contents = xrd_client.dirlist(path)[1]
    num_found = dir_contents.size

    # Handle emtpy directories
    if num_found == 0:
        print("Found zero objects in {}. Skipping.".format(path))
        return None

    # Allow user to choose directory if more than one is found
    if num_found > 1:
        print("Unexpected directory structure. Found {} objects in {}".format(num_found, path))
        print("Please type the number of the directory you would like to use. Options are:")
        for i, x in enumerate(dir_contents):
            print(i, x.name)
        dir_ix = input() # fixme: check input
    else:
        dir_ix = 0

    return sample_path + "/" + dir_contents.dirlist[dir_ix].name


parser = argparse.ArgumentParser()
parser.add_argument("-o", "--output_cfg", dest="cfg", required=True,
                    help="Path to output config, e.g. '../configs/ntuple_locations.yaml'")
parser.add_argument("-n", "--name", dest="name", required=True,
                    help="Name of group of ntuples, e.g. 'ffntuple_v4'")
parser.add_argument("-c", "--comment", dest="comment", required=True,
                    help=("Comment to describe group of ntuples, e.g. "
                          "'Most recent ntuples from Weinan -- only includes 2mu2e'"))
parser.add_argument("-d", "--directory", dest="directory", required=True,
                    help=("Path to ntuple root directory, e.g. "
                    "'root://cmseos.fnal.gov//store/group/lpcmetx/SIDM/ffNtupleV4/2018/'"))
args = parser.parse_args()

# Set up xrd client
redirector = args.directory.split("//store")[0]
xrd_client = client.FileSystem(redirector)

output = {
    args.name: {
        "path": args.directory,
        "samples": {}
    }
}
ntuple_path = args.directory.split(redirector)[1]
samples = xrd_client.dirlist(ntuple_path)[1]

# Traverse ntuple directory and construct output dictionary
# Assumes same structure as root://cmseos.fnal.gov//store/group/lpcmetx/SIDM/ffNtupleV4/2018/
for sample in samples:
    if sample.name.endswith("tmp"): # fixme: would be better to check if dir
        continue
    simple_name = parse_name(sample.name)
    output[args.name]["samples"][simple_name] = {}
    sample_path = sample.name

    # Descend three layers, expecting to find a single directory at each
    try:
        for _ in range(3):
            sample_path = descend(ntuple_path, sample_path)
            if sample_path is None:
                raise StopIteration()
    except StopIteration:
        continue

    # If traversal was successful, add path and files to output dictionary
    output[args.name]["samples"][simple_name]["path"] = sample_path + "/"
    files = [f.name for f in xrd_client.dirlist(ntuple_path + sample_path)[1]]
    output[args.name]["samples"][simple_name]["files"] = files

with open('test.yaml', 'a') as out_file:
    out_file.write("\n\n#" + args.comment + "\n")
    yaml.dump(output, out_file, default_flow_style=False)
