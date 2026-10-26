
import numpy as np

from pipnick.utils import masking


def test_polygons():
    vertices = np.array([[-0.5,33.5], [34.5, -0.5], [-0.5, -0.5], [-0.5,33.5]])
    shape = (900,1200)
    mask = masking.mask_polygons(shape, [vertices])

    assert mask[0,0], 'Should mask the first pixel'
    assert np.sum(mask) == 595, 'Should mask 595 pixels in total'


def test_create_mask():
    shape = (1024,1056)
    bad_columns = [255, 256, 783, 784, 1002]
    bad_slices = [np.s_[10:30,10:30], np.s_[1000:,1000:]]
    bad_regions = [
        np.array([[-0.5,33.5], [34.5, -0.5], [-0.5, -0.5], [-0.5,33.5]]),
        np.array([[-0.5, 960.5], [64.5, 1024.5], [-0.5, 1024.5]])
    ]

    # Mask only the bad columns
    mask = masking.create_mask(shape, bad_cols=bad_columns)
    assert np.sum(np.sum(mask, axis=0) > 0) == len(bad_columns), (
        'Incorrect number of columns masked'
    )
    assert np.sum(mask) == len(bad_columns) * shape[0], 'Incorrect number of pixels masked'

    # Mask only the bad slices
    mask = masking.create_mask(shape, bad_slices=bad_slices)
    assert mask[1001,1001], 'Masked slice incorrect'
    assert np.sum(mask) == 1744, 'Incorrect number of pixels masked'

    # Mask only the bad regions
    mask = masking.create_mask(shape, bad_regions=bad_regions)
    assert mask[-1,0] and mask[0,0], 'Incorrect region mask'
    assert np.sum(mask) == 2611, 'Incorrect number of pixels masked'

    # Mask everything
    mask = masking.create_mask(
        shape, bad_cols=bad_columns, bad_slices=bad_slices, bad_regions=bad_regions
    )
    assert np.all(mask[:,bad_columns[0]]), 'Bad column not fully masked'
    assert mask[1001,1001], 'Masked slice incorrect'
    assert mask[-1,0] and mask[0,0], 'Incorrect region mask'
    assert np.sum(mask) == 9346, 'Incorrect number of pixels masked'

