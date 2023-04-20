from typing import Dict
from os.path import join, abspath, dirname
from credentials import get_credentials

FILENAME = join(abspath(dirname(__file__)), "ERS_credentials.txt")

def get_ERS_credentials(filename: str = FILENAME) -> Dict[str, str]:
    return get_credentials(
        filename=filename,
        displayed=["username", "header"],
        hidden=["password"],
        prompt="credentials for EROS Registration System https://ers.cr.usgs.gov/register"
    )
