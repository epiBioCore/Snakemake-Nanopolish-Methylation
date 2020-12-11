def get_bioc_TxDb_pkg(wildcards):
    """Get the package bioconductor package name for the the species in config.yaml"""
    species = config["txdb"]["species"].capitalize()
    Source = config["txdb"]["Source"]
    build = config["txdb"]["build"]
    version = config["txdb"]["version"]

    if Source == "UCSC":
        return "TxDb.{species}.UCSC.{build}.knownGene".format(species=species,build=build)
    elif Source == "Ensembl":
        return "EnsDb.{species}.{version}".format(species=species,version=version)   

def get_bioc_pkg_txdb_path(wildcards):
    return "resources/bioconductor/lib/R/library/{pkg}".format(pkg=get_bioc_txdb_pkg(wildcards))

def get_bioc_species_pkg(wildcards):
    """Get the Txdb bioconductor package name for the the species in config.yaml"""
    species_letters = config["species"][0:2].capitalize()
    return "org.{species}.eg.db".format(species=species_letters)

def get_bioc_pkg_species_path(wildcards):
    return "resources/bioconductor/lib/R/library/{pkg}".format(pkg=get_bioc_species_pkg(wildcards))    

