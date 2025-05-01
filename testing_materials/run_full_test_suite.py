#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  1 08:57:13 2025

@author: js2746
"""
import argparse
import subprocess
from pathlib import Path
import os


def run_pytests_and_log(target_path, output_file):
    """
    Runs pytest on a specified file or directory, prints output to screen,
    and saves output to a file.

    :param target_path: Path to the test file or directory
    :param output_file: Path to the output log file
    """
    if not os.path.exists(target_path):
        print(f"Error: The path '{target_path}' does not exist.")
        return

    try:
        # Run pytest and capture output
        process = subprocess.Popen(
            ['pytest', '-v', target_path],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True
        )

        with open(output_file, 'w', encoding='utf-8') as f:
            for line in process.stdout:
                print(line, end='')  # Print to console
                f.write(line)        # Write to file

        process.wait()
        print(f"\nPytest finished. Output saved to '{output_file}'.")

    except Exception as e:
        print(f"An error occurred while running pytest: {e}")


def run_vmd_and_log_output(vmd_alias, args, output_file, windows, command_desc):
    """
    Runs vmd with arguments and redirects output to a file.

    :param vmd_alias: The command to run VMD on your computer (str or Path)
    :param args: List of command-line arguments (list of str)
    :param output_file: Path to the file where output will be saved (str or Path)
    :param windows: If True, use WINDOWS commands (bool)
    :param command_desc: The description for this command (str)
    """
    output_file = Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)  # Ensure directory exists

    cmd = [str(vmd_alias)] + args
    if windows:
        cmd = ["start"] + cmd
    print(f"Starting {command_desc}")

    try:
        with output_file.open('w', encoding='utf-8') as f:
            process = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True
            )

            for line in process.stdout:
                # print(line, end='')  # Also print to console
                f.write(line)        # Write to file

        process.wait()
        print(f"\n{command_desc} completed. Output saved to '{output_file}'.")
    except Exception as e:
        print(f"An error occurred while running {command_desc}: {e}")


def search_for_keywords(file_path, keywords):
    try:
        with open(file_path, 'r', encoding='utf-8') as file:
            contents = file.read()

        # Normalize text for case-insensitive search
        contents_lower = contents.lower()

        found_keywords = {}
        for keyword in keywords:
            count = contents_lower.count(keyword.lower())
            found_keywords[keyword] = count

        return found_keywords

    except FileNotFoundError:
        print(f"Error: The file at '{file_path}' was not found.")
    except Exception as e:
        print(f"An error occurred: {e}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run unit tests on .tcl and .py portions of codebase and then run end-to-end tests.")
    parser.add_argument("VMD_alias", help="The alias for VMD on your computer")
    parser.add_argument("-w", "--windows_compatible", help="Use this flag if you are running this on a WINDOWS machine")
    args = parser.parse_args()

    # resolve current directory
    cwd = Path.cwd().resolve()

    # run tcl unit tests and divert output to file
    unit_test_log = cwd.joinpath("tcl_unit_test.log")
    run_vmd_and_log_output(args.VMD_alias, ["-dispdev", "none", "-eofexit", "-e", str(cwd.joinpath("tcltest", "unit_test.test").resolve()), "-args", str(cwd)], unit_test_log, args.windows_compatible, "TCL unit testing")

    # check to make sure VMD didn't error
    if not unit_test_log.exists():
        print("Error: Something went wrong. VMD may not have run unit tests.")
        exit()

    # check to make sure all tests passed
    kw_counts = search_for_keywords(unit_test_log, ["invalid command", "bad option", "nt load file", "unknown option", "FAILED"])
    for count in kw_counts.values():
        if count != 0:
            "Something went wrong during unit testing."
            exit()

    # run python unit tests and divert output to file
    pytest_log = cwd.joinpath("pyunittest.log")
    run_pytests_and_log(cwd.joinpath("tcltest", "Unit_Test.py"), pytest_log)

    # check to make sure tests didn't error out
    kw_counts = search_for_keywords(pytest_log, ["no module named", "FAILED"])
    for count in kw_counts.values():
        if count != 0:
            "Something went wrong during unit testing."
            exit()

    # run nougat.tcl and nougat.py on test systems
    acc_test_gen_log = cwd.joinpath("nougat_test_outputs.log")
    run_vmd_and_log_output(args.VMD_alias, ["-dispdev", "none", "-eofexit", "-e", str(cwd.joinpath("run_nougat_test.tcl").resolve())], acc_test_gen_log, args.windows_compatible, "generation of acceptance test data")

    # check to make sure nougat didn't error out
    kw_counts = search_for_keywords(acc_test_gen_log, ["can't", "couldn't"])
    for count in kw_counts.values():
        if count != 0:
            "Something went wrong during generation of acceptance test data."
            exit()

    # run acceptance tests on test system results
    run_pytests_and_log(cwd.joinpath("tests"), cwd.joinpath("pytest.log"))

    # sound the all clear
    print("Acceptance testing finished")
