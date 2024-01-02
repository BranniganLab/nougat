"""
Created on Thu Aug 17 10:44:25 2023.

@author: js2746
"""
import shutil


def pytest_sessionfinish(session, exitstatus):
    """Run this hook if all pytests pass or xfail."""
    if exitstatus == 0:
        print("\nRemoving test files upon successful test completion!")
        shutil.rmtree("E-protein_trajectory/test_cart_5_5_0_-1_1")
        shutil.rmtree("E-protein_trajectory/test_polar_3_12_0_-1_1")
        shutil.rmtree("flat_surface_test/test_cart_5_5_0_-1_1")
        shutil.rmtree("flat_surface_test/test_polar_3_12_0_-1_1")
