import json
import pandas as pd
from typing import Iterable
from sqlalchemy import create_engine, text, Engine
import logging

logger = logging.getLogger(__name__)


def get_db_config(config_file: str):
    logger.info(f"Loading database configuration from {config_file}")
    with open(config_file) as f:
        config = json.load(f)

    # Check that all required fields are present
    required_fields = ["db_name", "db_user", "db_password", "db_host", "db_port"]
    for field in required_fields:
        if field not in config:
            raise ValueError(
                f"Error: {field} was not found. Please check the db_config.json file."
            )

    return config


def get_db_engine(db_config: dict) -> Engine:
    logger.info("Creating database engine for")
    engine = create_engine(
        f"postgresql://{db_config['db_user']}:{db_config['db_password']}@{db_config['db_host']}:{db_config['db_port']}/{db_config['db_name']}"
    )
    logger.info("Database engine created")
    return engine


def get_library_sequences(
    sp_names: Iterable[str], library: str, engine: Engine
) -> pd.DataFrame:
    logger.info("Getting SP sequences based on library from database")
    names_joined = ",".join(f"'{name}'" for name in sp_names)
    sp_library_df = pd.read_sql(
        text(
            "SELECT sp.name, amino_acid_sequence, sp_to_sp_library.dna_sequence_sp_only "
            "FROM sp_library "
            "JOIN sp_to_sp_library ON sp_library.id = sp_to_sp_library.sp_library_id "
            "JOIN sp ON sp_to_sp_library.sp_id = sp.id "
            f"WHERE sp_library.name = '{library}'"
            f"AND sp.name IN ({names_joined})"
        ),
        engine,
    )

    if sp_library_df.shape[0] != len(sp_names):
        raise ValueError(
            f"Expected {len(sp_names)} sequences, but got {sp_library_df.shape[0]} sequences from database"
        )

    return sp_library_df
