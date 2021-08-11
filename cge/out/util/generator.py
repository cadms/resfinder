#!/usr/bin/env python3

from git import Repo
from git.exc import InvalidGitRepositoryError
from datetime import datetime, timezone

from ..result import Result


class Generator():
    """ """

    @staticmethod
    def init_software_result(name, gitdir):
        """"""
        version, commit = Generator.get_version_commit(gitdir)
        date = datetime.now(timezone.utc).date().isoformat()

        result_dict = {
            "type": "software_result",
            "software_name": name,
            "software_version": version,
            "software_commit": commit,
            "run_date": date,
            "key": "{}-{}".format(name, version)
        }
        return Result(**result_dict)

    @staticmethod
    def get_version_commit(gitdir):
        """
            Input: Path to git directory
            Return: (version, commmit_hash)

            commmit_hash: The 40 character hash describing the exact commit of
            the git directory.
            version: The tag of the current commit. If no tag exists the first
            7 characters of the commit hash will be returned instead.
        """
        try:
            repo = Repo(gitdir)
        except InvalidGitRepositoryError:
            return ("unknown", "unknown")

        com2tag = {}
        for tag in repo.tags:
            com2tag[tag.commit.hexsha] = str(tag)

        version = com2tag.get(repo.commit().hexsha, repo.commit().hexsha[:7])

        return (version, repo.commit().hexsha)
