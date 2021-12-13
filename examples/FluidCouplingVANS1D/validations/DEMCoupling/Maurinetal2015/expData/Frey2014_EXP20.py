from numpy import *
zData = array(
        [
                2.48366013, 2.54901961, 2.61437908, 2.67973856, 2.74509804, 2.81045752, 2.87581699, 2.94117647, 3.00653595, 3.07189542, 3.1372549, 3.20261438,
                3.26797386, 3.33333333, 3.39869281, 3.46405229, 3.52941176, 3.59477124, 3.66013072, 3.7254902, 3.79084967, 3.85620915, 3.92156863, 3.9869281,
                4.05228758, 4.11764706, 4.18300654, 4.24836601, 4.31372549, 4.37908497, 4.44444444, 4.50980392, 4.5751634, 4.64052288, 4.70588235, 4.77124183,
                4.83660131, 4.90196078, 4.96732026, 5.03267974, 5.09803922, 5.16339869, 5.22875817, 5.29411765, 5.35947712, 5.4248366, 5.49019608, 5.55555556,
                5.62091503, 5.68627451, 5.75163399, 5.81699346, 5.88235294, 5.94771242, 6.0130719, 6.07843137, 6.14379085, 6.20915033, 6.2745098, 6.33986928,
                6.40522876, 6.47058824, 6.53594771, 6.60130719, 6.66666667, 6.73202614, 6.79738562, 6.8627451, 6.92810458, 6.99346405, 7.05882353, 7.12418301,
                7.18954248, 7.25490196, 7.32026144, 7.38562092, 7.45098039, 7.51633987, 7.58169935, 7.64705882, 7.7124183, 7.77777778, 7.84313725, 7.90849673,
                7.97385621, 8.03921569, 8.10457516, 8.16993464, 8.23529412, 8.30065359, 8.36601307, 8.43137255, 8.49673203, 8.5620915, 8.62745098, 8.69281046,
                8.75816993, 8.82352941, 8.88888889, 8.95424837, 9.01960784, 9.08496732, 9.1503268, 9.21568627, 9.28104575, 9.34640523, 9.41176471, 9.47712418,
                9.54248366
        ]
)
alphasData = array(
        [
                0.5046172, 0.51226814, 0.52388062, 0.54253272, 0.56002043, 0.57348531, 0.58443779, 0.58344348, 0.57304945, 0.5557779, 0.53469547, 0.51128047,
                0.49038335, 0.47736993, 0.4772912, 0.48613011, 0.49834204, 0.51155588, 0.52308532, 0.5329939, 0.53868952, 0.54064241, 0.5375283, 0.52766685,
                0.5174465, 0.50544434, 0.49389126, 0.4864515, 0.48629442, 0.49366087, 0.50805599, 0.52599695, 0.54307099, 0.55499335, 0.55989816, 0.55776025,
                0.54845714, 0.53614654, 0.52124182, 0.50657195, 0.49456883, 0.48583609, 0.48265908, 0.48568332, 0.49443324, 0.50793791, 0.52486548, 0.54169568,
                0.5544221, 0.56002623, 0.5562887, 0.5425574, 0.52213226, 0.500436, 0.48028134, 0.4641239, 0.45301304, 0.44728913, 0.44628308, 0.44925339,
                0.45519523, 0.46239793, 0.46966561, 0.4740174, 0.47185341, 0.46119903, 0.44138215, 0.41377148, 0.38101508, 0.34678851, 0.3132361, 0.28215449,
                0.2543112, 0.22988094, 0.20907622, 0.19219254, 0.17924567, 0.1695847, 0.16199408, 0.15504848, 0.14756032, 0.13906461, 0.12934315, 0.11866579,
                0.10764874, 0.09700555, 0.08690017, 0.0773169, 0.06831156, 0.05993213, 0.05238007, 0.04579115, 0.04008645, 0.03508413, 0.03061402, 0.0265941,
                0.02293822, 0.01959042, 0.01653656, 0.01377699, 0.01132335, 0.00917354, 0.0073393, 0.00580755, 0.0045244, 0.00346685, 0.00261227, 0.00194493,
                0.00143133
        ]
)
usData = array(
        [
                1.23502690e-04, 1.31195680e-04, 1.57470210e-04, 1.71027280e-04, 1.81200370e-04, 1.83694970e-04, 1.81261390e-04, 1.79387760e-04, 1.81396240e-04,
                1.89638470e-04, 2.06923250e-04, 2.32725200e-04, 2.62505840e-04, 2.94988850e-04, 3.24570640e-04, 3.53877640e-04, 3.94443130e-04, 4.37602170e-04,
                4.75908550e-04, 5.10095490e-04, 5.34748080e-04, 5.45805420e-04, 5.55180720e-04, 5.77110430e-04, 6.13283150e-04, 6.71062740e-04, 7.58822550e-04,
                8.75724190e-04, 1.01541800e-03, 1.15989610e-03, 1.29460420e-03, 1.41483370e-03, 1.52768170e-03, 1.63763040e-03, 1.74461190e-03, 1.83765100e-03,
                1.92436960e-03, 2.01219400e-03, 2.13950340e-03, 2.32916870e-03, 2.59259080e-03, 2.96056600e-03, 3.42859720e-03, 3.97389320e-03, 4.56246580e-03,
                5.12976560e-03, 5.63158300e-03, 6.06266410e-03, 6.44196670e-03, 6.78798550e-03, 7.13755490e-03, 7.55058080e-03, 8.09798570e-03, 8.85392420e-03,
                9.91822630e-03, 1.13627590e-02, 1.32440560e-02, 1.55353650e-02, 1.81165590e-02, 2.08450760e-02, 2.35978580e-02, 2.62765190e-02, 2.87802150e-02,
                3.11168790e-02, 3.33835580e-02, 3.57248170e-02, 3.83646630e-02, 4.15697970e-02, 4.56183760e-02, 5.07671880e-02, 5.73698380e-02, 6.57669230e-02,
                7.62996540e-02, 8.92153260e-02, 1.04456210e-01, 1.21556070e-01, 1.39507840e-01, 1.57139260e-01, 1.73668550e-01, 1.89089370e-01, 2.03788340e-01,
                2.18056830e-01, 2.32505780e-01, 2.47390870e-01, 2.62377860e-01, 2.76594480e-01, 2.89949470e-01, 3.03023510e-01, 3.15903890e-01, 3.28855170e-01,
                3.41467940e-01, 3.53286460e-01, 3.64081700e-01, 3.73937150e-01, 3.82969100e-01, 3.91288760e-01, 3.99086560e-01, 4.06637650e-01, 4.14074800e-01,
                4.21336910e-01, 4.28544680e-01, 4.35403680e-01, 4.41416290e-01, 4.46874160e-01, 4.52121170e-01, 4.57144320e-01, 4.62375740e-01, 4.67121200e-01,
                4.71449160e-01
        ]
)