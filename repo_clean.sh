# Assuming we're at the top of the repo

# Remove all cache directories from knitting
# This forces recalculation of any cached objects
find . -type d -name "*cache" | xargs rm -rf

# Remove all objects in the results/rdata folder
# This forces saving the geo_so, etc. objects at the end of lessons
find . -type d -name "rdata" | xargs rm -rf