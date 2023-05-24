from typing import Dict
from os.path import join, abspath, dirname, exists
from credentials import get_credentials

FILENAME = join(abspath(dirname(__file__)), "ERS_credentials.txt")

def get_ERS_credentials(filename: str = FILENAME) -> Dict[str, str]:
    if filename is None or not exists(filename):
        filename = FILENAME

    credentials = get_credentials(
        filename=filename,
        displayed=["username", "header"],
        hidden=["password"],
        prompt="credentials for EROS Registration System https://ers.cr.usgs.gov/register"
    )

    return credentials
