#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 13 08:15:57 2023

@author: jje63
"""

import pytest
import sys
import os
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

def test_multiple_key_value_pairs(self):
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