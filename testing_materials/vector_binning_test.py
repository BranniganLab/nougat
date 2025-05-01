import pytest
import sys
import os
sys.path.append(os.path.abspath('../python/'))
from vecbinning import VectorBinning, Multibinner, convert_cartesian_to_polar, convert_radians_and_degrees
import numpy as np

# Converting positive x and y values to polar coordinates
def test_positive_values_conversion():
    # Setup
    x_values = [1, 2, 3]
    y_values = [1, 2, 3]
    
    # Execute
    result = convert_cartesian_to_polar(x_values, y_values)
    
    # Assert
    expected_r = np.sqrt(np.array(x_values)**2 + np.array(y_values)**2)
    expected_theta = np.arctan(np.divide(np.array(y_values), np.array(x_values)))

    np.testing.assert_array_almost_equal(result[0], expected_r)
    np.testing.assert_array_almost_equal(result[1], expected_theta)


# Converting a list of radian values to degrees using default parameter
def test_radians_to_degrees_default_parameter():
    # Arrange
    radian_values = [0, np.pi/4, np.pi/2, np.pi, 2*np.pi]
    expected_degrees = [0, 45, 90, 180, 360]
    
    # Act
    result = convert_radians_and_degrees(radian_values)
    
    # Assert
    np.testing.assert_array_almost_equal(result, expected_degrees, decimal=10)

# Handling very large values that might cause precision issues
def test_large_values_precision():
    # Arrange
    large_radian_values = [1e10, 1e15, 1e20]
    expected_degrees = [np.array(large_radian_values) * (180/np.pi)]
    
    # Act
    result = convert_radians_and_degrees(large_radian_values)
    
    # Assert
    # For very large values, we use a relative tolerance
    # to account for floating point precision issues
    np.testing.assert_allclose(
        result, 
        expected_degrees[0], 
        rtol=1e-10
    )

# Initialization with regular data and default binning creates correct bins
def test_initialization_with_regular_data():
    # Arrange
    data = [1.2, 3.5, 7.8, 9.1, 4.6]
    n_bins = 5
    
    # Act
    binning = VectorBinning(data, n_bins)
    
    # Assert
    assert binning.prebinneddata.tolist() == data
    assert binning.N_bins == n_bins
    assert binning.bin_interval[0] == 1.0  # floor of min
    assert binning.bin_interval[1] == 10.0  # ceil of max
    assert len(binning.bins) == n_bins
    assert binning.bins[0] == 1.0
    assert binning.bins[-1] == 10.0
    assert binning.bin_indicies is not None
    assert len(binning.bin_indicies) == len(data)

# Correctly calculates bin intervals using floor of min and ceil of max values
def test_bin_intervals_calculation():
    # Create test data
    test_data = np.array([1.2, 3.7, 5.9, 2.4])
    n_bins = 5

    # Act
    binning = VectorBinning(test_data, n_bins)  # Use the correct class name
    binning.calculate_bins()

    # Assert
    assert binning.bin_interval[0] == np.floor(min(test_data))
    assert binning.bin_interval[1] == np.ceil(max(test_data))
    assert len(binning.bins) == n_bins
    assert binning.bins[0] == binning.bin_interval[0]
    assert binning.bins[-1] == binning.bin_interval[1]

# Handles empty prebinneddata array
def test_empty_prebinneddata_array():
    # Create empty test data
    test_data = np.array([])
    n_bins = 5

    # Act & Assert
    with pytest.raises(ValueError):  # Assuming min/max on empty array raises ValueError
        binning = VectorBinning(test_data, n_bins)
        binning.calculate_bins()

# Correctly assigns bin indices to data points within the range of bins
def test_assigns_correct_bin_indices():
    # Create test data
    test_data = np.array([1.0, 2.5, 3.7, 5.0])
    bins = np.array([0.0, 2.0, 4.0, 6.0])

    # Create instance with test data
    binning = VectorBinning(test_data, N_bins=3)

    # Override the bins for controlled testing
    binning.bins = bins
    binning.prebinneddata = test_data

    # Act
    binning.update_bin_indicies()

    # Assert
    expected_indices = np.array([1, 2, 2, 3])  # Based on np.digitize behavior
    np.testing.assert_array_equal(binning.bin_indicies, expected_indices)


# Handles empty prebinneddata array
def test_handles_empty_data():
    # Create non-empty test data to initialize the instance
    initial_data = np.array([1.0, 2.0, 3.0])
    binning = VectorBinning(initial_data, N_bins=3)

    # Override the data and bins for controlled testing
    empty_data = np.array([])
    bins = np.array([0.0, 2.0, 4.0, 6.0])
    binning.prebinneddata = empty_data
    binning.bins = bins

    # Act
    binning.update_bin_indicies()

    # Assert
    assert isinstance(binning.bin_indicies, np.ndarray)
    assert len(binning.bin_indicies) == 0

# Initialize Multibinner with valid multi-dimensional data and matching N_bins
def test_initialize_with_valid_data():
    # Arrange
    multi_bin_data = [[1, 2, 3, 4], [5, 6, 7, 8]]
    n_bins = [3, 4]
    
    # Act
    binner = Multibinner(multi_bin_data, n_bins)
    
    # Assert
    assert binner.multi_bin_data == multi_bin_data
    assert binner.N_bins == n_bins
    assert binner.theta == [None, None]
    assert binner.raw_data == []
    assert binner.raw_bin_indicies == []
    assert binner.raw_binning_data == []
    assert binner.multi_bin_dictionary == {}

# Initialize with empty data arrays
def test_initialize_with_empty_arrays():
    # Arrange
    multi_bin_data = [[], []]
    n_bins = [3, 4]
    
    # Act
    binner = Multibinner(multi_bin_data, n_bins)
    
    # Assert
    assert binner.multi_bin_data == [[], []]
    assert binner.N_bins == [3, 4]
    assert binner.theta == [None, None]
    assert binner.raw_data == []
    assert binner.raw_bin_indicies == []
    assert binner.raw_binning_data == []
    assert binner.multi_bin_dictionary == {}