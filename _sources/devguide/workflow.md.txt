# Developer Workflow

To participate in the OpenSn development, you will need to fork the upstream repo on GitHub.
Go to [OpenSn repo on GitHub](https://github.com/Open-Sn/opensn) and [create your personal fork](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/working-with-forks/fork-a-repo).

Then, clone the fork to the machine where you will be developing:
```shell
$ git clone git@github.com:<username>/opensn.git
```
where `<username>` is your GitHub username.
This, assumes, you [set up your SSH keys for accessing GitHub](https://docs.github.com/en/authentication/connecting-to-github-with-ssh).

Your fork will be referred to as `origin` (both in the git commands and in the text below).
Also set up an `upstream` remote for pulling changes from the upstream repo you forked from.
```shell
$ git remote add upstream git@github.com:Open-Sn/opensn.git
```

OpenSn development uses pull request workflow which is common on GitHub.
This means all development is done on so-called feature branches.
Do **not** develop on the `main` branch, you will most likely create problems for yourself down the road.

## Create a Branch

First, make sure you do not have any local changes:

```shell
$ git status
On branch main
Your branch is up to date with 'origin/main'.

nothing to commit, working tree clean
```

Branch off of `main`:

```shell
$ git checkout -b <branch name> main
Switched to a new branch '<branch name>'
```

Some tips on naming your branch:

- Use lower-case letters.
- Dashes (`-`) work better than underscores (`_`).
- If working on an issue, you can name your branch `issue/<issue number>`.


**Useful commands:**

| Action                    | git command     |
|:--------------------------|:----------------|
| To get a list of branches | `git branch -a` |
| To see your active branch | `git branch`    |

## Create a Commit

In git, you have to stage your changes first.
Then you check in those changes locally, and it will form a commit (also referred to as a patch).
Every commit has a SHA1 hash which uniquely identifies it in the repository.
Here is an example:

```shell
commit <SHA1>
Author: A. U. Thor <a.u.thor@somewhere.com>
Date:   <Date>
```

**Useful commands:**

| Action                                                             | git command                     |
|:-------------------------------------------------------------------|:--------------------------------|
| To add (stage) a new file                                          | `git add <file>`                |
| To add (stage) modification on a file                              | `git add <file>`                |
| To add all new files and local changes                             | `git add -A`                    |
| To add only modified files                                         | `git add -u`                    |
| To move a file                                                     | `git mv <source> <destination>` |
| To remove a file                                                   | `git rm <file>`                 |
| To see the status (staged, unstaged changes and non-tracked files) | `git status`                    |
| To commit locally                                                  | `git commit`                    |

The commit message should have this format:

```shell
Short description (up to 70 chars)
<empty line>
Detailed description. List as much useful and relevant information as possible.
This can have multiple lines and even paragraphs.
```

The detailed description will help a) the reviewer, b) in future when we forget what we did and will need to look it up.
Even though the detailed description is optional, it is *highly* recommended that you take the time to write it.

**Tips for creating commits:**

- You want your commits to be as small as possible, so they are easier (and therefore faster) to review.
- Stay on one topic per commit.
- Try **not** to do multiple things in one commit.
- Do a series of patches rather than one big patch.

**Useful commands:**

| Action                                        | git command  |
|:----------------------------------------------|:-------------|
| To see what files have changed                | `git status` |
| To look at list of commits                    | `git log`    |
| To look at code changes made in commits       | `git log -p` |
| To see local changes that were not staged yet | `git diff`   |

## Sending a Pull Request

When you think you are finished with your branch, send your changes for review.

Push the branch into your fork on GitHub

```shell
$ git push origin <branch name>
Enumerating objects: 38, done.
Counting objects: 100% (38/38), done.
Delta compression using up to 16 threads
Compressing objects: 100% (32/32), done.
Writing objects: 100% (32/32), 10.37 KiB | 5.19 MiB/s, done.
Total 32 (delta 13), reused 0 (delta 0)
remote: Resolving deltas: 100% (4/4), completed with 4 local objects.
To github.com:<username>/opensn.git
+ fffea64..f1b3ad0 <branch name> -> <branch name>
```

1. Go to `http://github.com/<username>/opensn.git`
2. [Create a pull request](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request-from-a-fork) to `main` branch

At this point, someone will review your branch and possibly leave comments on what to fix.
Otherwise, your branch will be merged.
You will get a notification (email) in both cases.
If your branch was merged, you should go update your `main` branch by pulling from the `upstream`

```shell
$ git checkout main
$ git pull upstream main
```

and you can delete your local branch with:

```shell
$ git branch -d <branch name>
```

## Fixing your branch

The easiest way how to fix your branch is by using fixup commits (see [fun with autosquash](https://technosorcery.net/blog/2010/02/fun-with-the-upcoming-1-7-release-of-git-rebase---interactive---autosquash/)).

1. Go patch by patch and fix the code according to the reviewer comments.
2. Stage your changes. Typically `git add -u` will do.
3. Do `git commit --fixup=<commit>` where `<commit>` is the SHA1 hash of the patch you are fixing.

   Alternatively, you can use the aliases mentioned in the article.

4. Continue until all is fixed.
5. Now, do `git rebase --interactive --autosquash main`. This will open an editor. Save the file you see. Quit the editor, and let git do the hard work.
6. If all went right, you will have your branch fixed. You can look at it by: `git log -p`.
7. Push your fixed up branch: `git push -f`

**Useful commands:**

| Action                                    | git command                    |
|:------------------------------------------|:-------------------------------|
| Create a fixup patch                      | `git commit --fixup=<sha1>`    |
| Add staged changes to the top-most commit | `git commit --amend --no-edit` |


## Updating Your Branch

If you wish to use the latest `main` branch (assuming it changed since you started your branch), you will have to rebase your branch on top of it, i.e.:

```shell
$ git checkout main
$ git pull upstream main
$ git checkout <my-branch-name>
$ git rebase main
```
Note: You will most likely run into conflicts that you will have to resolve.

Do **NOT** merge `main` into your feature branch.
