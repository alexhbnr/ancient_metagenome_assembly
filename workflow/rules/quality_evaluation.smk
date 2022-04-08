rule checkM_prepareDatabase:
    output:
        f"{config['resourcedir']}/checkM/.dmanifest"
    message: "Download the checkM databases and extract the tarballs"
    params:
        outdir = f"{config['resourcedir']}/checkM",
        url = "https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz"
    wrapper:
        "file:///home/alexander_huebner/github/snakemake-wrappers/bio/checkm/preparedatabase"

rule checkM_setRoot:
    input:
        f"{config['resourcedir']}/checkM/.dmanifest"
    output:
        touch(f"{config['resourcedir']}/checkM/setRoot.done")
    message: "Specify the database folder in checkM"
    params:
        dbdir = f"{config['resourcedir']}/checkM"
    wrapper:
        "file:///home/alexander_huebner/github/snakemake-wrappers/bio/checkm/setroot"



