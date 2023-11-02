# How to archive a local branch without pushing to a repo

Git has an archive ability that will not grab stuff that will otherwise 
not be in the repo.

```shell
git archive --format=tar.gz -o <name_you_want>.tar.gz --prefix=<folder_name_when_extracted>/ <branch_name>
```
Take special note of the forward slash after `<folder_name_when_extracted>`

## Example:
Say I want to archive a special branch with secret code. The branch name is
`do_not_push`. For this I would do, from within the `chi-tech` repo's main
folder.

```shell
git archive --format=tar.gz -o secret_stuff.tar.gz --prefix=yes_secret/ do_not_push
```

## Excluding folders
Maybe you don't want documentation in the archive. Well this is kind of hard.
Best move? Do the above. Extract it somewhere not controlled by git. Delete 
the folder you dont want, then re-tar.