import yaml
import os
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
        'projectName': ["Project Name", 'Test'],
        'dataDir': ["Path to raw data", os.path.abspath('data/raw')],
        'outDir': ["Path to output directory", os.path.abspath('data/processed')],
        'sampleFile': ["Path to file with sample names", os.path.abspath('code/configs/samples.txt')],
        'species': ["Bacterial species", "Salmonella_enterica"],
        'fq_fwd': ['Fastq file forward suffix', ".1.fq.gz"],
        'fq_rvr': ['Fastq file reverse suffix', ".2.fq.gz"],
        'reference': ['Path to reference genome', os.path.abspath("data/processed/ref/genome.fasta")]
    }
    extra_mappings = {
        'qc': ['Run qc (True/False)', 'True'],
        'mink': ['BBmap mink (?) ', '11'],
        'trimq': ['BBmap trimpq (?)', '14'],
        'mapq': ['BBmap mapq', '20'],
        'minlen': ['BBmap minlen', '45'],
        'merged': ['Merge reads for assembly? (True/False)', 'False'],
        'fastqc': ['Run FastQC? (True/False)', 'True'],
        'scaffold': ['Enter scaffold file (?)', ''],
        'ram': ['RAM (for what?)', '4'],
        'phyloConfig': ['Phylophlan config'],
        'fnaFolder': ['fnaFolder'],
        'gbkFolder': ['gbkFolder'], #todo what are these settings, come up with reasonable defaults,
        'panXpath': ['panX path'],
        'cg': ['cg ?'],
        'aribaDB': ['aribaDB'],
        'adapters': ['adapters'],

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
        new_val = get_new_value(basic_mappings[setting][0], basic_mappings[setting][1])
        project_dict[setting] = new_val
    process_extra = input("See extra settings? y/n\n\n")
    if process_extra == 'n':
        for setting in extra:
            project_dict[setting] = config_dict[setting]
        return project_dict
    else:
        for setting in extra:
            if setting in extra_mappings.keys():
                if len(extra_mappings[setting]) > 1:
                    default = extra_mappings[setting][1]
                else:
                    default = config_dict[setting]
                new_val = get_new_value(extra_mappings[setting][0], default)
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