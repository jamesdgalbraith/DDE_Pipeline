#!/usr/bin/env python
# adapted from https://wellsr.com/python/3-ways-to-calculate-python-execution-time/

# import the builtin time module
import time

# Grab Currrent Time Before Running the Code
start = time.time()

# Put commands here

# Grab Currrent Time After Running the Code
end = time.time()

#Subtract Start Time from The End Time
total_time = end - start
print("\n"+ str(total_time))