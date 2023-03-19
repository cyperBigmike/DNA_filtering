@echo off

rem List of input values
set input_values=1 2 3 4 5 6 7 8

rem Loop through each input value and run the Python script with it
for %%i in (%input_values%) do (
    python implementations\filter.py %%i
)