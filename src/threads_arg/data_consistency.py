# This to make a data-consistent extension

# need: 
#   - one pass through data
#   - age of every mutation
# (if not provided, estimate using Kimura’s age = -2N ln(1 - freq), and derived = minor)

# for every threading instruction, keep track of "new threading instructions" and "current constraints"
# to do this, for every site, update constraints and break segment (starting new segment) if constraints are broken

