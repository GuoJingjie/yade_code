# 2019 © Remi Monthiller <remi.monthiller@gmail.com>
#-------------------#
# Measures
#-------------------#
json_save = False
measures = {
        "z": "[i * pF.dz for i in range(pN.n_z)]",
        "phi": "getDryProfile('phiPart')",
        "vx": "getVxPartProfile()",
        "dirs": "getOrientationHist(5.0, 0.0)",
        "mdirs": "getVectorMeanOrientation()",
        "ori": "getOrientationProfiles(pM.z_ground, pF.dz*20, int(pN.n_z/20))",
}
