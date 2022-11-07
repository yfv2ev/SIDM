"""Define all commonly used objects"""

obj_defs = {
    # 3=pfmuon, 8=dsamuon, 2=pfelectron, 4=pfphoton
    "mu_ljs" : lambda objs: objs["ljs"][(objs["ljs"]["type"] == 3) | (objs["ljs"]["type"] == 8)],
    "egm_ljs" : lambda objs: objs["ljs"][(objs["ljs"]["type"] == 2) | (objs["ljs"]["type"] == 4)],
}
