module data_hifi:
    snakefile: "00-data/hifi/Snakefile"
    config: config
use rule * from data_hifi
