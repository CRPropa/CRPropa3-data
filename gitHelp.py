"""gitHelp

Functionality to access information on the git repository.
"""

import subprocess

def get_git_revision_hash() -> str:
    """ get the full hash of the repository """
    return subprocess.check_output(['git', 'rev-parse', 'HEAD']).decode('ascii').strip()
    
def get_git_revision_hash_short() -> str:
    """ get the short hash of the repository """
    return subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD']).decode('ascii').strip() 
