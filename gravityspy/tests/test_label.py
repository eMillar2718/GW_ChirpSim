"""Unit test for GravitySpy
"""

__author__ = 'Scott Coughlin <scott.coughlin@ligo.org>'

import os

import gravityspy.ml.read_image as read_image
import gravityspy.ml.labelling_test_glitches as label_glitches
import gravityspy.ml.train_classifier as train_classifier

import pandas as pd
import numpy

TEST_IMAGES_PATH = os.path.join(os.path.split(__file__)[0], 'data',
'images')
MODEL_NAME_CNN = os.path.join(os.path.split(__file__)[0], '..', '..', 'models',
                              'O3-multiview-classifer.h5')
MODEL_NAME_FEATURE_MULTIVIEW = os.path.join(os.path.split(__file__)[0], '..', '..', 'models',
                                            'similarity-model-O3.h5')
MULTIVIEW_FEATURES_FILE = os.path.join(os.path.split(__file__)[0], 'data',
                                       'MULTIVIEW_FEATURES.npy')

SCORE = 0.9987664222717285

FEATURES = numpy.array([[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
136.32681274414062, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
0.0, 0.0, 0.0, 0.0, 143.1021728515625, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
0.0, 0.0, 0.0, 0.0, 0.0, 74.27071380615234, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 110.8934326171875, 0.0, 0.0,
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
142.75534057617188, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
0.0, 120.54851531982422, 132.58575439453125, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 125.55241394042969, 0.0, 0.0, 0.0, 0.0, 0.0,
0.0, 0.0, 0.0, 154.26510620117188, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 34.673343658447266, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]])

MULTIVIEW_FEATURES = numpy.load(MULTIVIEW_FEATURES_FILE)

class TestGravitySpyML(object):
    """`TestCase` for the GravitySpy
    """
    def test_label(self):

        list_of_images = []
        for ifile in os.listdir(TEST_IMAGES_PATH):
            if 'spectrogram' in ifile:
                list_of_images.append(ifile)

        image_dataDF = pd.DataFrame()
        for idx, image in enumerate(list_of_images):
            image_data = read_image.read_grayscale(os.path.join(
                                                       TEST_IMAGES_PATH,
                                                       image),
                                          resolution=0.3)

            image_dataDF[image] = [image_data]

        # Now label the image
        scores, MLlabel, _, _, _, _, _, = label_glitches.label_glitches(
                                                        image_data=image_dataDF,
                                                        model_name='{0}'.format(MODEL_NAME_CNN),
                                                        order_of_channels="channels_last",
                                                        image_order=['0.5.png', '1.0.png', '2.0.png', '4.0.png'],
                                                        image_size=[140, 170],
                                                        verbose=False)

        confidence = float(scores[0][MLlabel])
        assert confidence == SCORE


    def test_multiview_rgb(self):

        list_of_images = []
        for ifile in os.listdir(TEST_IMAGES_PATH):
            if 'spectrogram' in ifile:
                list_of_images.append(ifile)

        image_dataDF = pd.DataFrame()
        for idx, image in enumerate(list_of_images):
            image_data_r, image_data_g, image_data_b = read_image.read_rgb(os.path.join(
                                                       TEST_IMAGES_PATH,
                                                       image),
                                          resolution=0.3)

            image_dataDF[image] = [[image_data_r, image_data_g, image_data_b]]

        # Now label the image
        features, _ = label_glitches.get_multiview_feature_space(
                                image_data=image_dataDF,
                                semantic_model_name='{0}'.format(MODEL_NAME_FEATURE_MULTIVIEW),
                                order_of_channels="channels_last",
                                image_size=[140, 170], verbose=False)

        numpy.testing.assert_array_almost_equal(features, MULTIVIEW_FEATURES,
                                                decimal=3)
