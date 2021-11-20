# 2019 © Remi Monthiller <remi.monthiller@gmail.com>
#-------------------#
# Measures
#-------------------#
json_save = False
measures = {
        "phi": "hydroEngine.phiPart[:]",
        "vx": "getVxPartProfile()",
        "vfx": "hydroEngine.vxFluid[0:-1]",
        "ReS": "hydroEngine.ReynoldStresses[:]",
        "lm": "hydroEngine.lm[:]",
        "shields": "getShields()",
        "dirs": "getOrientationHist(5.0, 0.0)",
        "mdirs": "getVectorMeanOrientation()",
        "ori": "getOrientationProfiles(pM.z_ground, pF.dz*20, int(pN.n_z/20))",
}
