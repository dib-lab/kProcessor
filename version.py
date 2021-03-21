import os
import urllib.request, json
import semantic_version
import sys

STAGE = str()
if len(sys.argv) > 1:
    if sys.argv[1] != "release":
        STAGE = "dev"
else:
    STAGE = "release"


def is_github_action():
    if "GITHUB_WORKFLOW" in dict(os.environ.items()):
        return True
    else:
        return False
    

def get_pypa_dev_latest():
    with urllib.request.urlopen("https://test.pypi.org/pypi/kProcessor/json") as url:
        data = json.loads(url.read().decode())
        return data["info"]["version"]


MAJOR = 1
MINOR = 1
PATCH = 0

dev_version = semantic_version.Version(major=MAJOR, minor=MINOR, patch=PATCH, prerelease=('dev', '0'))
release_version = semantic_version.Version(major=MAJOR, minor=MINOR, patch=PATCH)

def get_version():
    version_tag = str()

    if STAGE == "release":
        version_tag = release_version
    else:
        # If it's running on github action, increment the dev patch number
        if is_github_action():
            test_pypa_latest_version = get_pypa_dev_latest()
            version_tag = semantic_version.Version(test_pypa_latest_version).next_patch()
        
        # Running on local machine
        else:
            version_tag = dev_version

    return version_tag


__version__ = "1.1.1-dev0" #get_version()
