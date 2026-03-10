"""
COMPADRE database I/O: fetch and load.
"""

const COMPADRE_URL = "https://compadre-db.org/Data/CompadreDownload"
const COMADRE_URL = "https://compadre-db.org/Data/ComadreDownload"

function MPM.cdb_fetch(db_type::Symbol=:compadre;
                       save::Bool=false, destination::Union{Nothing,String}=nothing)
    url = db_type == :compadre ? COMPADRE_URL : COMADRE_URL
    # TODO: Implement download and RData parsing
    error("COMPADRE extension not yet fully implemented. Use cdb_load with a local file.")
end

function MPM.cdb_load(path::String)
    # Load RData file
    rdata = RData.load(path)
    # TODO: Parse the R list structure into CompadreDB
    error("COMPADRE loading not yet fully implemented")
end

function MPM.cdb_save(cdb::CompadreDB, path::String)
    # TODO: Save as Julia-native format
    error("COMPADRE saving not yet fully implemented")
end
