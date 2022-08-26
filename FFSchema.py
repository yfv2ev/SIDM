from coffea.nanoevents.schemas.base import BaseSchema, zip_forms
from coffea.nanoevents import transforms
import importlib
import utilities
importlib.reload(utilities)


class FFSchema(BaseSchema):
    """
    FF schema builder

    The FF schema is built from all branches found in the supplied file, which
    is intended to be a FireFighter ntuple (https://github.com/phylsix/Firefighter).
    Analogous to the procedure in in TreeMakerSchema, array collections are 
    generated in three steps:

    1. The branches of vector-like objects are grouped into a single collection
       with the corresponding coordinate variables names mapped to the standard
       variable names for coffea.nanoevents.methods.vector behaviors. For example:
       - muons will behave as LorentzVectors
       - primary vertices will behave as ThreeVectors
       - missing energy will behave as a TwoVector
    
    2. The branches corresponding to attributes of a subobject are merged into 
       a single collection associated with that subobject. For example:
       - all pfjets_pfcands_<attribute> branches are merged into a single collection
    
    3. The branches corresponding to the attributes and subobjects associated with
       a given object are grouped into a single collection named <object>.
       For example:
       - all pfjet_<attribute> branches are merged into a single collection so
         muon attributes can be accessed like pfjet.<attribute>
       - all pfjet_<subobject> branches are merged into a single collection so
         <subobject> attributes can be accessed like pfjet.<subobject>.<attribute>


    All collections are then zipped into one `base.NanoEvents` record and returned.
    """

    def __init__(self, base_form):
        super().__init__(base_form)
        self._form["contents"] = self._build_collections(self._form["contents"])

    def _build_collections(self, branch_forms):
        # define special cases
        counts_branches = {
            "pfjet_pfcand" : "pfjet_pfcands_n",
            "pfjet_pfcands" : None,
        }
        
        # Turn any special classes (e.g. LorentzVectors) into the appropriate awkward form
        vector_objects = list(set(b.split("/")[0] for b in branch_forms if "/" in b))

        for obj in vector_objects:
            components = set(k.split('.')[-1] for k in branch_forms if k.startswith(obj + "/"))
            # optional fixme: add case for candidates (lorentz+charge)
            # optional fixme: add case for pfjet_pfcand, which could be PtEtaPhiELorentzVector
            # handle lorentz vectors
            if components == {"fX", "fY", "fZ", "fT"}:
                form = zip_forms(
                    {
                        "x": branch_forms.pop(f"{obj}/{obj}.fCoordinates.fX"),
                        "y": branch_forms.pop(f"{obj}/{obj}.fCoordinates.fY"),
                        "z": branch_forms.pop(f"{obj}/{obj}.fCoordinates.fZ"),
                        "t": branch_forms.pop(f"{obj}/{obj}.fCoordinates.fT"),
                    },
                    obj,
                    "LorentzVector",
                )
                branch_forms[obj] = form
            # handle three-vectors
            elif components == {"fX", "fY", "fZ"}:
                form = zip_forms(
                    {
                        "x": branch_forms.pop(f"{obj}/{obj}.fCoordinates.fX"),
                        "y": branch_forms.pop(f"{obj}/{obj}.fCoordinates.fY"),
                        "z": branch_forms.pop(f"{obj}/{obj}.fCoordinates.fZ"),
                    },
                    obj+"_p3",
                    "ThreeVector",
                )
                branch_forms[obj+"_p3"] = form
            # handle two-vectors
            elif components == {"fX", "fY"}:
                form = zip_forms(
                    {
                        "x": branch_forms.pop(f"{obj}/fCoordinates/fCoordinates.fX"),
                        "y": branch_forms.pop(f"{obj}/fCoordinates/fCoordinates.fY"),
                    },
                    obj+"_p2",
                    "TwoVector",
                )
                branch_forms[obj+"_p2"] = form
            else:
                raise ValueError(
                    f"Unrecognized class with split branches: {components}"
                )

        # identify branches with one value per event (e.g. "lumi")
        single_value_branches, remaining_branches = utilities.partition_list(
            branch_forms,
            lambda x: (
                "_" not in x
                or x.split("_")[1] == "n"
                or x.startswith(("HLT", "tomatchfilter", "akjet_ak4PFJetsCHS_n")) # special cases
            )
        )
        
        # identify base objects (e.g. "muon")
        objects = list(set(b.split("_")[0] for b in remaining_branches))
        objects = [o if o != "akjet" else "akjet_ak4PFJetsCHS" for o in objects] # special case
        
        for obj in objects:
            # identify all object attributes
            # e.g. "pt" from "muon_pt" or "pfcand_pt" from "pfjet_pfcand_pt"
            attributes = [b.split(obj+"_")[1] for b in remaining_branches if b.startswith(obj+"_")]
            
            # identify subobjects (e.g. "pfcand" from "pfjet_pfcand_pt")
            subobjects = [a.split("_")[0] for a in attributes if "_" in a]
            subobjects = list(set(s for s in subobjects if subobjects.count(s) > 1))
            
            # distinguish between attributes of objects and of subobjects
            attributes, subattributes = utilities.partition_list(
                attributes,
                lambda x: not x.startswith(tuple(subobjects))
            )
            
            # create subobject jagged arrays
            for subobj in subobjects:
                # get subobject offsets 
                counts_branch = counts_branches.get("_".join((obj, subobj)), "_".join((obj, subobj, "n")))
                if counts_branch in remaining_branches:
                    offsets = transforms.counts2offsets_form(branch_forms.pop(counts_branch))
                else:
                    offsets = None
            
                base_name = "_".join((obj, subobj))
                branch_forms[base_name] = zip_forms(
                    {
                        a.split(subobj+"_")[1] : branch_forms.pop("_".join((obj, a)))
                        for a in subattributes
                        if a.startswith(subobj+"_") and not a.endswith("_n")
                    },
                    base_name,
                    offsets,
                )
                attributes.append(subobj)
                
            # create object jagged arrays, including nesting of subobjects
            counts_name = "_".join((obj, "n"))
            if counts_name in remaining_branches:
                offsets = transforms.counts2offsets_form(branch_forms.pop("_".join((obj, "n"))))
            else:
                offsets = None
            branch_forms[obj] = zip_forms(
                {a : branch_forms.pop("_".join((obj, a))) for a in attributes},
                obj,
                offsets=offsets
            )
            
        return branch_forms

    @property
    def behavior(self):
        """Behaviors necessary to implement this schema"""
        from coffea.nanoevents.methods import base, vector

        behavior = {}
        behavior.update(base.behavior)
        behavior.update(vector.behavior)
        return behavior