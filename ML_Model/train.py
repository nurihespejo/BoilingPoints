import logging
import pathlib as pl
import sys

import pandas as pd

import mother.train as mother
import mother.utils as utils
from mother.config_definition import ModelFlag

# set current working directory
repo_dir: pl.Path = pl.Path(__file__).parent
logs_dir: pl.Path = repo_dir.joinpath("logs")
if not logs_dir.exists():
    logs_dir.mkdir()
# setting up logging
logger: logging.Logger = utils.setup_logging(logging.INFO, folder=logs_dir)

# loading data from disk
data_path: pl.Path = repo_dir.joinpath( "data", "compounds.csv")
compounds: pd.DataFrame = pd.read_csv(data_path)

# create a new "config" folder in the working directory with a default YAML configuration
# please comment this line or the config folder will be overwritten with default configuration
if not repo_dir.joinpath("config").exists():
    logger.info("Creating default config for training.")
    utils.create_config_folder(dest_dir=repo_dir)
    logger.warn("Please configure you project properly before doing the training.")
    sys.exit(1)

# define training configuration, will override defaults from config folder.
model_flags: ModelFlag = ModelFlag.PREDICT | ModelFlag.UNCERTAINTY | ModelFlag.EXPLAIN

logger.info("Starting Training")
pipeline, config = mother.train(
    compounds=compounds, model_config=model_flags, project_path=repo_dir
)

logger.info("Training done, exporting using export_model")
model_name = f"{repo_dir.name}_model.pkl.gz"
pipeline.export_model(config, model_name=model_name)
logger.info("Script done.")
sys.exit(0)
