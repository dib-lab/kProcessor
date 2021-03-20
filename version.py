import os
import subprocess

def get_version():
    kProcessor_version = 1.1
    version_tag = str()

    if os.path.isdir(os.path.join(os.path.dirname(os.path.realpath(__file__)), ".git")):
        commit_hash_short_name = subprocess.getoutput("git rev-parse --short HEAD").split()[0]
        branch_name = subprocess.getoutput("git rev-parse --abbrev-ref HEAD").split()[0]
        if branch_name == "master":
            version_tag = kProcessor_version
        else:
            version_tag = f"{kProcessor_version}.dev0+{commit_hash_short_name}"

    else:
        version_tag = kProcessor_version
        pass

    return version_tag


__version__ = get_version()
print(__version__)