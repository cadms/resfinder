#!/usr/bin/env python3
import os.path

from git import Repo


# git_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
git_path = os.path.abspath(os.path.dirname(__file__))
print("Git path: {}".format(git_path))
repo = Repo(git_path)
commit = repo.commit()
print("commit: {}".format(commit.name_rev))
