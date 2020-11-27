import yaml

"""
1. read in the default config file 
2. go through every option (with defaults)
3. get input, write new yaml to the path of new config file 

"""


def config_to_dict(default_config):
    config_dict = {}

    with open(default_config) as file:
        # The FullLoader parameter handles the conversion from YAML
        # scalar values to Python the dictionary format
        config_dict = yaml.load(file, Loader=yaml.FullLoader)
    return config_dict

def get_user_input(config_dict):
    project_dict = {}

    basic_mappings = {
        'projectName': "Project Name",
        'dataDir': "Path to raw data",
        'outDir': "Path to output directory",
        'sampleFile': "Path to file with sample names",
        'species': "Bacterial species",
        'fq_fwd': 'Fastq file forward suffix',
        'fq_rvr': 'Fastq file reverse suffix',
        'reference': 'Path to reference genome',



    }
    extra_mappings = {
        'qc': 'Run qc (True/False)',
        'mink': 'BBmap mink (?) ',
        'trimq': 'BBmap trimpq (?)',
        'mapq': 'BBmap mapq',
        'minlen': 'BBmap minlen',
        'merged': 'Merge reads for assembly? (True/False)',
        'fastqc': 'Run FastQC? (True/False)',
        'scaffold': 'Enter scaffold file (?)',
        'ram': 'RAM (for what?)',
        'phyloConfig': 'Phylophlan config',
        'fnaFolder': 'fnaFolder',
        'gbkFolder': 'gbkFolder', #todo what are these settings, come up with reasonable defaults,
        'panXpath': 'panX path',
        'cg': 'cg ?',
        'aribaDB': 'aribaDB',
        'adapters': 'adapters',

    }

    basic = []
    extra = []

    for config in config_dict.keys():
        if config in basic_mappings.keys():
            basic.append(config)
        else:
            extra.append(config)

    print("Configuring settings for the project")
    print("Basic Settings")
    print("Enter a new value, press enter to keep default, or type 'null' to leave blank")
    for setting in basic:
        new_val = get_new_value(basic_mappings[setting], config_dict[setting])
        project_dict[setting] = new_val
    process_extra = input("See extra settings? y/n\n\n")
    if process_extra == 'n':
        return project_dict
    else:
        for setting in extra:
            if setting in extra_mappings.keys():
                new_val = get_new_value(extra_mappings[setting], config_dict[setting])
            else:
                new_val = get_new_value(setting, config_dict[setting])
            project_dict[setting] = new_val
    print('Configuration Complete')
    return project_dict

def get_new_value(name, default):
    new_val = input(f'Enter {name}. Default= {default}.\n\n\t')
    if not new_val:
        new_val = default
    elif new_val == 'null':
        new_val = ''
    return new_val


def write_new_config(project_dict, new_config):
    with open(new_config, 'w') as fh:
        yaml.dump(project_dict, fh)

def configure(default_config, new_config):
    config_dict = config_to_dict(default_config)
    project_dict = get_user_input(config_dict)
    write_new_config(project_dict, new_config)
    return project_dict


if __name__ == "__main__":
    df = "configs/test_config.yaml"
    nc = "configs/test2_config.yaml"
    print(configure(df, nc))