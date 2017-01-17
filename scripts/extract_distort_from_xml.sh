# Take only distortion lines from list of input xml files, and remove the "
grep --no-filename "distortion model=" $* |tr -d \"
