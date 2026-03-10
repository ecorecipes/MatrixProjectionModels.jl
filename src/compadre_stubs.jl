"""
Stub functions for COMPADRE/COMADRE database extension.
Implementations provided by MatrixProjectionModelsCompadreExt.
"""

# Types (defined here so they can be referenced without loading extension)
abstract type AbstractCompadreDB end

# IO
function cdb_fetch end
function cdb_load end
function cdb_save end

# Accessors
function cdb_matA end
function cdb_matU end
function cdb_matF end
function cdb_matC end
function cdb_metadata end
function cdb_id end

# Quality flags
function cdb_flag end

# Operations
function cdb_collapse end
function cdb_rbind end
function cdb_flatten end
function cdb_subset end
function cdb_build_cdb end

# Statistics
function mpm_mean end
function mpm_sd end
function mpm_median end

# Stage queries
function mpm_has_prop end
function mpm_has_active end
function mpm_has_dorm end
function mpm_first_active end
