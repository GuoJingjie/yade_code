# 2019 © Remi Monthiller <remi.monthiller@gmail.com>
#-------------------#
# Measures
#-------------------#
json_save = False
measures = {
        "z": "[i * pF.dz for i in range(pN.n_z)]",
        #"phi":"getDryProfile('phiPart')", # For dry cases
        "phi": "hydroEngine.phiPart[:]",  # For wet cases
        "vx": "getVxPartProfile()",
        "vfx": "hydroEngine.vxFluid[0:-1]",  # Only in wet cases
        "ReS": "hydroEngine.ReynoldStresses[:]",  # Only in wet cases
        "lm": "hydroEngine.lm[:]",  # Only in wet cases
        "shields": "getShields()",  # Only in wet cases
        "dirs": "getOrientationHist(5.0, 0.0)",
        "mdirs": "getVectorMeanOrientation()",
        "ori": "getOrientationProfiles(pM.z_ground, pF.dz*20, int(pN.n_z/20))",
}
