# activate rclone conda env first
# Incorporate this into pipeline eventually
batch=$1

#conda activate rclone

rclone copy /well/bag/users/vbr851/projects/agnostic_diagnostic_v2/output/$batch neonatal:agnostic_diagnostic/$batch --ignore-existing

#conda deactivate