#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 13 08:15:57 2023

@author: jje63
"""

import pytest
import sys
import os
import numpy as np
sys.path.insert(0,'../plotting/')
import utils



def test_multiple_lines_with_text():
        file = ["line 1", "line 2", "line 3"]
        result = list(utils.strip_blank_lines(file))
        assert result == ["line 1", "line 2", "line 3"]

def test_one_line_with_text():
        file = ["line 1"]
        result = list(utils.strip_blank_lines(file))
        assert result == ["line 1"]
        
def test_tabs_file():
        file = ["\t\t\t"]
        result = list(utils.strip_blank_lines(file))
        assert result == []

def test_multiple_key_value_pairs():
        # Create a temporary config file with multiple key/value pairs
        config_file = open("temp_config.txt", "w")
        config_file.write("key1=value1\n")
        config_file.write("key2=value2\n")
        config_file.write("key3=value3\n")
        config_file.close()

        # Call the read_config function with the temporary config file
        config_dict = utils.read_config("temp_config.txt")

        # Check if the config_dict contains the correct key/value pairs
        assert config_dict == {"key1": "value1", "key2": "value2", "key3": "value3"}

        # Remove the temporary config file
        os.remove("temp_config.txt")
        
def test_config_with_blank_lines():
        # Create a temporary config file with blank lines
        config_file = open("temp_config.txt", "w")
        config_file.write("\n")
        config_file.write("key1=value1\n")
        config_file.write("\n")
        config_file.write("key2=value2\n")
        config_file.close()

        # Call the read_config function with the temporary config file
        config_dict = utils.read_config("temp_config.txt")

        # Check if the config_dict contains the correct key/value pairs (without blank lines)
        assert config_dict == {"key1": "value1", "key2": "value2"}

        # Remove the temporary config file
        os.remove("temp_config.txt")

def test_empty_config_file():
        # Create an empty temporary config file
        config_file = open("temp_config.txt", "w")
        config_file.close()

        # Call the read_config function with the temporary config file
        config_dict = utils.read_config("temp_config.txt")

        # Check if the config_dict is empty
        assert config_dict == {}

        # Remove the temporary config file
        os.remove("temp_config.txt")
        
def test_config_with_only_comments():
        # Create a temporary config file with only comments
        config_file = open("temp_config.txt", "w")
        config_file.write("# This is a comment\n")
        config_file.write("# Another comment\n")
        config_file.close()

        # Call the read_config function with the temporary config file
        config_dict = utils.read_config("temp_config.txt")

        # Check if the config_dict is empty (no key/value pairs)
        assert config_dict == {}

        # Remove the temporary config file
        os.remove("temp_config.txt")
        
def test_multiple_non_nan_values():
        input_list = [1.0, 2.0, 3.0]
        expected_output = 1.0
        assert utils.find_first_val(input_list) == expected_output
        
        
def test_mix_of_nan_and_non_nan_values():
        input_list = [np.nan, 1.0, np.nan, 2.0]
        expected_output = 1.0
        assert utils.find_first_val(input_list) == expected_output
        
def test_multiple_non_nan_values_last():
        in_list = [1.0, 2.0, 3.0]
        assert utils.find_last_val(in_list) == 3.0
        
def test_nan_and_non_nan_values():
        in_list = [1.0, np.nan, 2.0, np.nan, 3.0]
        assert utils.find_last_val(in_list) == 3.0
        
def test_num_123_size_5():
        assert utils.gifformat(123, 5) == "00123"
        
def test_num_123456_size_3():
        assert utils.gifformat(123456, 3) == "123456"

def test_num_123_size_negative_1():
        assert utils.gifformat(123, -1) == "123"
        
def test_height_measurement_with_dat_type():
        assert utils.filename_generator("sys1", "lipid1", "field1", "bead1", "coordsys1", "height", "dat") == "sys1.field1.bead1.coordsys1.avgheight.dat"
        
def test_invalid_system_name():
        with pytest.raises(TypeError):
            utils.filename_generator(123, "lipid6", "field6", "bead6", "coordsys6", "height", "dat")
            
def test_nans():
        matrix_data = np.array([[[np.nan, np.nan, np.nan], [np.nan, np.nan, np.nan], [np.nan, np.nan, np.nan]],
                               [[np.nan, np.nan, np.nan], [np.nan, np.nan, np.nan], [np.nan, np.nan, np.nan]],
                               [[np.nan, np.nan, np.nan], [np.nan, np.nan, np.nan], [np.nan, np.nan, np.nan]]])
        avg = utils.calc_avg_over_time(matrix_data)
        expected_avg = np.full((3, 3), np.nan)
        assert np.array_equal(avg, expected_avg, equal_nan=True)
        
def test_positive_integers():
        matrix_data = np.array([[[1, 2, 3], [4, 5, 6], [7, 8, 9]],
                               [[10, 11, 12], [13, 14, 15], [16, 17, 18]],
                               [[19, 20, 21], [22, 23, 24], [25, 26, 27]]])
        avg = utils.calc_avg_over_time(matrix_data)
        expected_avg = np.array([[2.0, 5.0, 8.0],
                                 [11.0, 14.0, 17.0],
                                 [20.0, 23.0, 26.0]])
        assert np.array_equal(avg, expected_avg)
        
def test_negative_integers():
        matrix_data = np.array([[[-1, -2, -3], [-4, -5, -6], [-7, -8, -9]],
                               [[-10, -11, -12], [-13, -14, -15], [-16, -17, -18]],
                               [[-19, -20, -21], [-22, -23, -24], [-25, -26, -27]]])
        avg = utils.calc_avg_over_time(matrix_data)
        expected_avg = np.array([[-2.0, -5.0, -8.0],
                                 [-11.0, -14.0, -17.0],
                                 [-20.0, -23.0, -26.0]])
        assert np.array_equal(avg, expected_avg)
        
def test_zeros():
        matrix_data = np.array([[[0.0, -0.0, 0.0], [-0.0, 0.0, -0.0], [0.0, -0.0, 0.0]],
                               [[-0.0, 0.0, -0.0], [0.0, -0.0, 0.0], [-0.0, 0.0, -0.0]],
                               [[0.0, -0.0, 0.0], [-0.0, 0.0, -0.0], [0.0, -0.0, 0.0]]])
        avg = utils.calc_avg_over_time(matrix_data)
        expected_avg = np.zeros((3, 3))
        assert np.array_equal(avg, expected_avg)

def test_returns_two_ndarrays():
        bin_info = {'N1': 10, 'N2': 5, 'd1': 0.1, 'd2': 0.2}
        coordsys = "cart"
        dim1vals, dim2vals = utils.bin_prep(bin_info, coordsys)
        assert isinstance(dim1vals, np.ndarray)
        assert isinstance(dim2vals, np.ndarray)

    
def test_rectangular_positive_int_values():
        N1_bins = 5
        d1 = 2
        N2_bins = 3
        d2 = 1
        min_val = 0
        coordsys = "rectangular"
        sys_name = "test"
        
        os.mkdir('npy')
        utils.save_areas(N1_bins, d1, N2_bins, d2, min_val, coordsys, sys_name)
        expected_areas = np.ones([N1_bins, N2_bins]) * d1 * d2
        np.testing.assert_array_equal(np.load('npy/' + sys_name + "." + coordsys + ".areas.npy"), expected_areas)
        os.remove('npy/' + sys_name + "." + coordsys + ".areas.npy")
        os.rmdir('npy')
        
def test_rectangular_zero_bins():
        N1_bins = 0
        d1 = 1
        N2_bins = 0
        d2 = 1
        min_val = 0
        coordsys = "rectangular"
        sys_name = "test"
        
        os.mkdir('npy')
        utils.save_areas(N1_bins, d1, N2_bins, d2, min_val, coordsys, sys_name)
    
        expected_areas = np.ones([N1_bins, N2_bins]) * d1 * d2
        np.testing.assert_array_equal(np.load('npy/' + sys_name + "." + coordsys + ".areas.npy"), expected_areas)
        os.remove('npy/' + sys_name + "." + coordsys + ".areas.npy")
        os.rmdir('npy')
        
def test_polar_zero_bins():
        N1_bins = 0
        d1 = 1
        N2_bins = 0
        d2 = 1
        min_val = 0
        coordsys = "polar"
        sys_name = "test"
        
        os.mkdir('npy')
        utils.save_areas(N1_bins, d1, N2_bins, d2, min_val, coordsys, sys_name)
    
        expected_areas = np.ones([N1_bins, N2_bins]) * d1 * d2
        np.testing.assert_array_equal(np.load('npy/' + sys_name + "." + coordsys + ".areas.npy"), expected_areas)
        os.remove('npy/' + sys_name + "." + coordsys + ".areas.npy")
        os.rmdir('npy')
        
        
def test_only_nan_values():
        data_array = np.full((3, 3, 3), np.nan)
        result = utils.mostly_empty(data_array)
        assert np.array_equal(result, data_array, equal_nan=True)
        
        
def test_one_frame():
        data_array = np.ones((3, 3, 1))
        result = utils.mostly_empty(data_array)
        assert np.array_equal(result, data_array, equal_nan=True)

def test_convert_to_cart_positive():
        assert utils.convert_to_cart(1, 0) == (1, 0)
        assert ("{:.5f}".format(utils.convert_to_cart(2, np.pi/4)[0]),"{:.5f}".format(utils.convert_to_cart(2, np.pi/4)[1])) == ("{:.5f}".format(np.sqrt(2)), "{:.5f}".format(np.sqrt(2)))
        assert ("{:.0f}".format(utils.convert_to_cart(3, np.pi/2)[0]),"{:.0f}".format(utils.convert_to_cart(3, np.pi/2)[1])) == ('0', '3')

def test_convert_to_cart_type_error():
        with pytest.raises(TypeError):
            utils.convert_to_cart("1", np.pi/6)
        with pytest.raises(TypeError):
            utils.convert_to_cart(2, "0")
        with pytest.raises(TypeError):
            utils.convert_to_cart("3", "0")
        with pytest.raises(TypeError):
            utils.convert_to_cart("4", "1")

def test_positive_float_input():
        result = utils.coord_format(12.345)
        assert isinstance(result, str)
        assert len(result) == 8
        assert result[4] == '.'
        assert result[5:] == '345'

def test_negative_float_input():
        result = utils.coord_format(-12.345)
        assert isinstance(result, str)
        assert len(result) == 8
        assert result[4] == '.'
        assert result[5:] == '345'

def test_integer_value():
        assert utils.bin_format(10) == ' 10.00'

# Returns a string with no spaces and the value with .00 appended to the end when given a zero value.
def test_zero_value():
        assert utils.bin_format(0) == '  0.00'
    
# Returns a string with two spaces and the value with .00 appended to the end when given a single-digit integer value.
def test_single_digit_integer_value():
        assert utils.bin_format(5) == '  5.00'