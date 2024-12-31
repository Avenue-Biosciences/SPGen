import json
from sqlalchemy import create_engine, text, Engine


def get_db_config(config_file: str):
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
    return create_engine(
        f"postgresql://{db_config['db_user']}:{db_config['db_password']}@{db_config['db_host']}:{db_config['db_port']}/{db_config['db_name']}"
    )
