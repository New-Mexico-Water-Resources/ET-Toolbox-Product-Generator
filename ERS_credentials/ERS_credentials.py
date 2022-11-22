from typing import Dict

from credentials import get_credentials


def get_ERS_credentials(filename: str = "~/.ERS") -> Dict[str, str]:
    return get_credentials(
        filename=filename,
        displayed=["username", "header"],
        hidden=["password"],
        prompt="credentials for EROS Registration System https://ers.cr.usgs.gov/register"
    )
